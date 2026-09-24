package cli

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"regexp"
	"slices"
	"strconv"
	"strings"

	"github.com/liminalpurple/vcfq/internal/ensembl"
	"github.com/liminalpurple/vcfq/internal/vcf"
)

const usageText = `vcfq — query a tabix-indexed VCF with Ensembl coordinate resolution.

Usage:
  vcfq [flags] <query>...
  vcfq cache clean
  vcfq version

Queries are auto-detected:
  rsNNNN              rsID (e.g. rs1801133)
  CHR:START[-END]     genomic region in GRCh38 (e.g. chr1:11796000-11796500)
  SYMBOL              gene symbol (e.g. HNF1A)

Multiple queries may be passed as arguments or piped on stdin (whitespace
separated). Stdin is consumed only when not attached to a terminal.

Flags:
`

// Run is the program entrypoint, separated from main() so it's testable.
func Run(args []string, stdin io.Reader, stdout, stderr io.Writer) int {
	if len(args) > 0 {
		switch args[0] {
		case "cache":
			return runCache(args[1:], stderr)
		case "version":
			if _, err := fmt.Fprintln(stdout, "vcfq", Version); err != nil {
				return 1
			}
			return 0
		case "help", "-h", "--help":
			printUsage(stderr)
			return 0
		}
	}

	var (
		vcfPath  string
		annotate bool
		format   string
		cacheDir string
	)

	fs := flag.NewFlagSet("vcfq", flag.ContinueOnError)
	fs.SetOutput(stderr)
	fs.Usage = func() { printUsageTo(stderr, fs) }
	fs.StringVar(&vcfPath, "vcf", os.Getenv("VCFQ_VCF"), "path to bgzipped tabix-indexed VCF (or set $VCFQ_VCF)")
	fs.BoolVar(&annotate, "a", false, "annotate variants with VEP (consequence, gene, AA change, AF)")
	fs.BoolVar(&annotate, "annotate", false, "annotate variants with VEP (long form of -a)")
	fs.StringVar(&format, "f", "tsv", "output format: tsv, table, json, vcf")
	fs.StringVar(&format, "format", "tsv", "output format (long form of -f)")
	fs.StringVar(&cacheDir, "cache", "", "override cache directory (default: $XDG_CACHE_HOME/vcfq)")
	if err := fs.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return 0
		}
		return 2
	}

	queries := append([]string(nil), fs.Args()...)
	queries = append(queries, readStdinTokens(stdin)...)
	if len(queries) == 0 {
		printUsageTo(stderr, fs)
		return 2
	}

	if vcfPath == "" {
		diag(stderr, "error: no VCF path; pass -vcf <path> or set VCFQ_VCF")
		return 2
	}

	cache, err := ensembl.NewCache(cacheDir)
	if err != nil {
		diag(stderr, "error:", err)
		return 1
	}

	reader, err := vcf.Open(vcfPath)
	if err != nil {
		diag(stderr, "error:", err)
		return 1
	}
	defer func() {
		if err := reader.Close(); err != nil {
			diag(stderr, "warning: close vcf:", err)
		}
	}()

	formatter, err := NewFormatter(format, stdout)
	if err != nil {
		diag(stderr, "error:", err)
		return 2
	}
	if err := formatter.Header(reader.Header, annotate); err != nil {
		diag(stderr, "error: write header:", err)
		return 1
	}

	client := ensembl.NewClient(Version)
	ctx := context.Background()

	var allRecords []Record
	for _, q := range queries {
		recs, err := resolveQuery(ctx, client, cache, reader, q)
		if err != nil {
			diag(stderr, "warning: "+q+":", err)
			continue
		}
		allRecords = append(allRecords, recs...)
	}

	if annotate && len(allRecords) > 0 {
		if err := annotateRecords(ctx, client, cache, allRecords); err != nil {
			diag(stderr, "warning: annotation failed:", err)
		}
	}

	for _, rec := range allRecords {
		if err := formatter.Write(rec); err != nil {
			diag(stderr, "error: write record:", err)
			return 1
		}
	}
	if err := formatter.Close(); err != nil {
		diag(stderr, "error: close formatter:", err)
		return 1
	}

	// After the output, so it reads as a footnote to the rows rather than
	// splitting a table header from its body on an interleaved terminal.
	if anyNoCall(allRecords) {
		diag(stderr, "note: no variant call is not the same as homozygous reference —")
		diag(stderr, "      verify coverage on the CRAM before interpreting an absence.")
	}
	return 0
}

// resolveQuery dispatches a single query token to the right resolver, returning
// the matching VCF records as cli.Records ready for formatting.
func resolveQuery(ctx context.Context, client *ensembl.Client, cache *ensembl.Cache, reader *vcf.Reader, q string) ([]Record, error) {
	switch Detect(q) {
	case KindRSID:
		return resolveRSID(ctx, client, cache, reader, q)
	case KindRegion:
		return resolveRegion(reader, q)
	default:
		return resolveSymbol(ctx, client, cache, reader, q)
	}
}

func resolveSymbol(ctx context.Context, client *ensembl.Client, cache *ensembl.Cache, reader *vcf.Reader, symbol string) ([]Record, error) {
	info, err := client.LookupSymbol(ctx, cache, symbol)
	if err != nil {
		return nil, fmt.Errorf("symbol lookup: %w", err)
	}
	lines, err := reader.ScanRegion(info.Chrom, info.Start, info.End)
	if err != nil {
		return nil, err
	}
	return splitToRecords(lines, symbol), nil
}

func resolveRSID(ctx context.Context, client *ensembl.Client, cache *ensembl.Cache, reader *vcf.Reader, rsid string) ([]Record, error) {
	loc, err := client.LookupRSID(ctx, cache, rsid)
	if err != nil {
		return nil, fmt.Errorf("rsID lookup: %w", err)
	}
	start, end := rsidWindow(loc)
	lines, err := reader.ScanRegion(loc.Chrom, start, end)
	if err != nil {
		return nil, err
	}
	// Prefer records where ID matches rsid; if none match the ID exactly, return
	// every variant at the location that could be this one (the VCF may not
	// annotate this rsID but the variant is still present). For an indel the
	// window includes the anchor base, so a SNV sitting on that base is a
	// neighbour, not this variant, and must not stand in for it.
	var matched, all []Record
	snv := isSNV(loc)
	for _, r := range splitToRecords(lines, rsid) {
		if strings.EqualFold(r.ID, rsid) {
			matched = append(matched, r)
		}
		if snv || len(r.Ref) != 1 || len(r.Alt) != 1 {
			all = append(all, r)
		}
	}
	if len(matched) > 0 {
		return matched, nil
	}
	if len(all) == 0 {
		return []Record{noCallRecord(loc.Chrom, loc.Start, loc.ID, loc.Ref, loc.Alts, rsid)}, nil
	}
	return all, nil
}

// noCallRecord builds the synthetic row emitted when a single-position query
// resolves cleanly but finds no line in the VCF. Reporting the resolved
// coordinates is the whole point: an empty result tells you nothing about where
// to look next, whereas "no call at chr6:32,638,107 (C/T)" is directly
// actionable against the CRAM.
func noCallRecord(chrom string, pos int, id, ref string, alts []string, query string) Record {
	return Record{
		Chrom:  chrom,
		Pos:    pos,
		ID:     id,
		Ref:    ref,
		Alt:    strings.Join(alts, ","),
		Query:  query,
		NoCall: true,
	}
}

// rsidWindow converts an Ensembl variant location into the VCF positions to
// scan. Ensembl reports indels without VCF's leading anchor base: a deletion
// spans only the deleted bases, and an insertion has start = end+1. VCF puts
// POS on the anchor base before them, so anything other than a plain SNV is
// widened one base to the left.
func rsidWindow(loc *ensembl.VariantLoc) (start, end int) {
	if loc.Start == loc.End && isSNV(loc) {
		return loc.Start, loc.End
	}
	return min(loc.Start-1, loc.End), max(loc.Start, loc.End)
}

func isSNV(loc *ensembl.VariantLoc) bool {
	if len(loc.Ref) != 1 || loc.Ref == "-" {
		return false
	}
	for _, a := range loc.Alts {
		if len(a) != 1 || a == "-" {
			return false
		}
	}
	return true
}

var regionParseRegex = regexp.MustCompile(`^((?:chr)?[\dXYMxym]+):(\d+)(?:-(\d+))?$`)

func resolveRegion(reader *vcf.Reader, region string) ([]Record, error) {
	m := regionParseRegex.FindStringSubmatch(region)
	if m == nil {
		return nil, fmt.Errorf("invalid region %q", region)
	}
	chrom := m[1]
	start, _ := strconv.Atoi(m[2])
	end := start
	if m[3] != "" {
		end, _ = strconv.Atoi(m[3])
	}
	if !strings.HasPrefix(chrom, "chr") {
		chrom = "chr" + chrom
	}
	lines, err := reader.ScanRegion(chrom, start, end)
	if err != nil {
		return nil, err
	}
	// Only a single-position query has an unambiguous "this exact site had no
	// call" reading. An empty range is just an empty range.
	if len(lines) == 0 && start == end {
		return []Record{noCallRecord(chrom, start, "", "", nil, region)}, nil
	}
	return splitToRecords(lines, region), nil
}

// splitToRecords expands multi-allelic lines into one Record per ALT and tags
// each with the originating query string. On multi-allelic lines only the
// called ALTs are emitted (both of them for a 1/2 het). Lines with no called
// ALT — hom-ref, missing, or single-allelic — emit every ALT, so users see
// ref/ref calls at positions they explicitly queried.
func splitToRecords(lines []vcf.DataLine, query string) []Record {
	out := make([]Record, 0, len(lines))
	for li := range lines {
		d := &lines[li]
		gt := vcf.ParseGT(d.Format, d.Sample)
		called := vcf.CalledAlts(gt)
		for i, alt := range d.Alts {
			if len(d.Alts) > 1 && len(called) > 0 && !slices.Contains(called, i+1) {
				continue
			}
			out = append(out, Record{
				Chrom:    d.Chrom,
				Pos:      d.Pos,
				ID:       d.ID,
				Ref:      d.Ref,
				Alt:      alt,
				Line:     d,
				AltIndex: i + 1,
				GT:       gt,
				Query:    query,
			})
		}
	}
	return out
}

func anyNoCall(recs []Record) bool {
	for _, r := range recs {
		if r.NoCall {
			return true
		}
	}
	return false
}

func annotateRecords(ctx context.Context, client *ensembl.Client, cache *ensembl.Cache, recs []Record) error {
	// No-call rows are deliberately excluded: annotating them would describe a
	// variant the sample does not carry, in the same columns used for variants it
	// does. targets keeps the mapping back from response index to record index.
	var (
		inputs  []ensembl.VEPInput
		targets []int
	)
	for i, r := range recs {
		if r.NoCall {
			continue
		}
		inputs = append(inputs, ensembl.VEPInput{Chrom: r.Chrom, Pos: r.Pos, Ref: r.Ref, Alt: r.Alt})
		targets = append(targets, i)
	}
	if len(inputs) == 0 {
		return nil
	}
	anns, err := client.AnnotateBatch(ctx, cache, inputs)
	if err != nil {
		return err
	}
	for i, a := range anns {
		if i >= len(targets) {
			break
		}
		rec := &recs[targets[i]]
		rec.HasAnnot = true
		rec.Consequence = a.Consequence
		rec.Gene = a.Gene
		rec.AAChange = a.AAChange
		rec.AF = a.AF
	}
	return nil
}

func runCache(args []string, stderr io.Writer) int {
	if len(args) == 0 {
		diag(stderr, "usage: vcfq cache clean")
		return 2
	}
	switch args[0] {
	case "clean":
		c, err := ensembl.NewCache("")
		if err != nil {
			diag(stderr, "error:", err)
			return 1
		}
		if err := c.Clean(); err != nil {
			diag(stderr, "error:", err)
			return 1
		}
		return 0
	default:
		diag(stderr, fmt.Sprintf("unknown cache subcommand %q", args[0]))
		return 2
	}
}

func printUsage(w io.Writer) {
	diag(w, strings.TrimSuffix(usageText, "\n"))
}

// diag writes one diagnostic line to w (normally stderr). Write errors are
// discarded deliberately: if stderr itself is broken there is nowhere left to
// report them, and the exit code still carries the outcome.
func diag(w io.Writer, a ...any) {
	_, _ = fmt.Fprintln(w, a...)
}

func printUsageTo(w io.Writer, fs *flag.FlagSet) {
	printUsage(w)
	fs.PrintDefaults()
}
