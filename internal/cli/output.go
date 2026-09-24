package cli

import (
	"encoding/json"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Formatter writes Records to an output stream. Header is called once with the
// source VCF header lines (only the vcf formatter uses them) and whether
// annotations are enabled. Write is called per record. Close flushes.
type Formatter interface {
	Header(vcfHeader []string, annotated bool) error
	Write(rec Record) error
	Close() error
}

// NewFormatter returns the formatter for the named format. Names: tsv, table,
// json, vcf.
func NewFormatter(name string, w io.Writer) (Formatter, error) {
	switch name {
	case "", "tsv":
		return &tsvFormatter{w: w}, nil
	case "table":
		return &tableFormatter{w: w}, nil
	case "json":
		return &jsonFormatter{w: w}, nil
	case "vcf":
		return &vcfFormatter{w: w}, nil
	default:
		return nil, fmt.Errorf("unknown format %q (want tsv, table, json, vcf)", name)
	}
}

// columns lists the output columns for tsv/table/json. Annotation columns are
// appended when annotated is true.
func columns(annotated bool) []string {
	cols := []string{"chrom", "pos", "id", "ref", "alt", "gt", "query"}
	if annotated {
		cols = append(cols, "consequence", "gene", "aa_change", "af")
	}
	return cols
}

func recordValues(rec Record, annotated bool) []string {
	vals := []string{
		rec.Chrom,
		strconv.Itoa(rec.Pos),
		emptyDot(rec.ID),
		emptyDot(rec.Ref),
		emptyDot(rec.Alt),
		gtDisplay(rec),
		rec.Query,
	}
	if annotated {
		af := "."
		if rec.HasAnnot && rec.AF > 0 {
			af = strconv.FormatFloat(rec.AF, 'g', 4, 64)
		}
		vals = append(vals,
			emptyDot(rec.Consequence),
			emptyDot(rec.Gene),
			emptyDot(rec.AAChange),
			af,
		)
	}
	return vals
}

// gtDisplay renders the genotype cell. No-call rows report "no-call" rather than
// a genotype, which is a weaker and more accurate claim than "0/0": vcfq knows
// only that the VCF has no line at the position, not why.
func gtDisplay(rec Record) string {
	if rec.NoCall {
		return "no-call"
	}
	return emptyDot(rec.GT)
}

func emptyDot(s string) string {
	if s == "" {
		return "."
	}
	return s
}

// --- tsv ---

type tsvFormatter struct {
	w           io.Writer
	annotated   bool
	wroteHeader bool
}

func (f *tsvFormatter) Header(_ []string, annotated bool) error {
	f.annotated = annotated
	_, err := fmt.Fprintln(f.w, strings.Join(columns(annotated), "\t"))
	f.wroteHeader = err == nil
	return err
}

func (f *tsvFormatter) Write(rec Record) error {
	if !f.wroteHeader {
		if err := f.Header(nil, f.annotated); err != nil {
			return err
		}
	}
	_, err := fmt.Fprintln(f.w, strings.Join(recordValues(rec, f.annotated), "\t"))
	return err
}

func (f *tsvFormatter) Close() error { return nil }

// --- table (markdown) ---

type tableFormatter struct {
	w           io.Writer
	annotated   bool
	wroteHeader bool
}

func (f *tableFormatter) Header(_ []string, annotated bool) error {
	f.annotated = annotated
	cols := columns(annotated)
	if _, err := fmt.Fprintf(f.w, "| %s |\n", strings.Join(cols, " | ")); err != nil {
		return err
	}
	seps := make([]string, len(cols))
	for i := range seps {
		seps[i] = "---"
	}
	_, err := fmt.Fprintf(f.w, "| %s |\n", strings.Join(seps, " | "))
	f.wroteHeader = err == nil
	return err
}

func (f *tableFormatter) Write(rec Record) error {
	if !f.wroteHeader {
		if err := f.Header(nil, f.annotated); err != nil {
			return err
		}
	}
	_, err := fmt.Fprintf(f.w, "| %s |\n", strings.Join(recordValues(rec, f.annotated), " | "))
	return err
}

func (f *tableFormatter) Close() error { return nil }

// --- json (newline-delimited) ---

type jsonFormatter struct {
	w         io.Writer
	annotated bool
}

func (f *jsonFormatter) Header(_ []string, annotated bool) error {
	f.annotated = annotated
	return nil
}

type jsonRow struct {
	Chrom       string   `json:"chrom"`
	Pos         int      `json:"pos"`
	ID          string   `json:"id,omitempty"`
	Ref         string   `json:"ref"`
	Alt         string   `json:"alt"`
	GT          string   `json:"gt,omitempty"`
	Query       string   `json:"query"`
	Source      string   `json:"source,omitempty"`
	Consequence string   `json:"consequence,omitempty"`
	Gene        string   `json:"gene,omitempty"`
	AAChange    string   `json:"aa_change,omitempty"`
	AF          *float64 `json:"af,omitempty"`
}

func (f *jsonFormatter) Write(rec Record) error {
	row := jsonRow{
		Chrom: rec.Chrom,
		Pos:   rec.Pos,
		ID:    rec.ID,
		Ref:   rec.Ref,
		Alt:   rec.Alt,
		GT:    rec.GT,
		Query: rec.Query,
	}
	if rec.NoCall {
		// ref/alt on these rows are Ensembl's, not the VCF's — say so, since a
		// consumer joining on this output has no other way to tell.
		row.GT = "no-call"
		row.Source = "ensembl"
	}
	if f.annotated && rec.HasAnnot {
		row.Consequence = rec.Consequence
		row.Gene = rec.Gene
		row.AAChange = rec.AAChange
		if rec.AF > 0 {
			af := rec.AF
			row.AF = &af
		}
	}
	b, err := json.Marshal(row)
	if err != nil {
		return err
	}
	_, err = fmt.Fprintln(f.w, string(b))
	return err
}

func (f *jsonFormatter) Close() error { return nil }

// --- vcf ---

// annotationInfoLines are appended to the source header (just before #CHROM) when
// annotations are enabled. Order is stable for reproducible output. Each field is
// Number=A: one value per ALT on the line, "." for ALTs that weren't annotated.
var annotationInfoLines = []string{
	`##INFO=<ID=VCFQ_CONSEQUENCE,Number=A,Type=String,Description="VEP consequence term">`,
	`##INFO=<ID=VCFQ_GENE,Number=A,Type=String,Description="VEP-resolved gene symbol">`,
	`##INFO=<ID=VCFQ_AA,Number=A,Type=String,Description="Amino acid change (ref/alt)">`,
	`##INFO=<ID=VCFQ_AF,Number=A,Type=Float,Description="Allele frequency from VEP (gnomAD/1000G)">`,
}

// vcfFormatter writes each source VCF line once, unsplit, with its original ALT
// list and sample column. Records arrive split per ALT, so consecutive records
// sharing a Line are buffered and emitted together.
type vcfFormatter struct {
	w         io.Writer
	annotated bool
	pending   []Record
}

func (f *vcfFormatter) Header(vcfHeader []string, annotated bool) error {
	f.annotated = annotated
	// Pass through header lines, injecting our INFO declarations just before the
	// final #CHROM... line. If the source has no header (empty slice), emit a
	// minimal one.
	if len(vcfHeader) == 0 {
		vcfHeader = []string{"##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}
	}
	for _, line := range vcfHeader {
		if strings.HasPrefix(line, "#CHROM") && annotated {
			for _, info := range annotationInfoLines {
				if _, err := fmt.Fprintln(f.w, info); err != nil {
					return err
				}
			}
		}
		if _, err := fmt.Fprintln(f.w, line); err != nil {
			return err
		}
	}
	return nil
}

func (f *vcfFormatter) Write(rec Record) error {
	// A no-call has no honest VCF data-line representation — every genotype this
	// format can express is a claim vcfq is not entitled to make. Skip the row and
	// let the stderr note carry it. (A ##vcfq_nocall header line would preserve it
	// for this format; see CLAUDE.md.)
	if rec.NoCall {
		return nil
	}
	if rec.Line == nil {
		return fmt.Errorf("vcf output: record %s:%d has no source line", rec.Chrom, rec.Pos)
	}
	if len(f.pending) > 0 && f.pending[0].Line != rec.Line {
		if err := f.flush(); err != nil {
			return err
		}
	}
	f.pending = append(f.pending, rec)
	return nil
}

func (f *vcfFormatter) Close() error { return f.flush() }

func (f *vcfFormatter) flush() error {
	if len(f.pending) == 0 {
		return nil
	}
	d := f.pending[0].Line
	info := d.Info
	if f.annotated {
		info = appendAnnotations(info, len(d.Alts), f.pending)
	}
	f.pending = f.pending[:0]
	cols := []string{
		d.Chrom,
		strconv.Itoa(d.Pos),
		emptyDot(d.ID),
		d.Ref,
		strings.Join(d.Alts, ","),
		emptyDot(d.Qual),
		emptyDot(d.Filter),
		emptyDot(info),
	}
	if d.Format != "" {
		cols = append(cols, d.Format)
		if d.Sample != "" {
			cols = append(cols, d.Sample)
		}
	}
	_, err := fmt.Fprintln(f.w, strings.Join(cols, "\t"))
	return err
}

// appendAnnotations injects Number=A VCFQ_* INFO fields into an existing INFO
// string, with one value per ALT (nAlts of them) taken from the annotated
// records for that line. A field with no value for any ALT is omitted. If the
// source INFO is empty or ".", the fields replace it; otherwise they're appended
// with ";" separators.
func appendAnnotations(info string, nAlts int, recs []Record) string {
	fields := []struct {
		key   string
		value func(Record) string
	}{
		{"VCFQ_CONSEQUENCE", func(r Record) string { return escapeInfoValue(r.Consequence) }},
		{"VCFQ_GENE", func(r Record) string { return escapeInfoValue(r.Gene) }},
		{"VCFQ_AA", func(r Record) string { return escapeInfoValue(r.AAChange) }},
		{"VCFQ_AF", func(r Record) string {
			if r.AF <= 0 {
				return ""
			}
			return strconv.FormatFloat(r.AF, 'g', 4, 64)
		}},
	}
	var parts []string
	for _, fld := range fields {
		vals := make([]string, nAlts)
		found := false
		for i := range vals {
			vals[i] = "."
		}
		for _, r := range recs {
			if !r.HasAnnot || r.AltIndex < 1 || r.AltIndex > nAlts {
				continue
			}
			if v := fld.value(r); v != "" {
				vals[r.AltIndex-1] = v
				found = true
			}
		}
		if found {
			parts = append(parts, fld.key+"="+strings.Join(vals, ","))
		}
	}
	added := strings.Join(parts, ";")
	if added == "" {
		return info
	}
	if info == "" || info == "." {
		return added
	}
	return info + ";" + added
}

// escapeInfoValue replaces characters that VCF reserves in INFO fields (space,
// semicolon, equals, comma) with safe alternatives.
func escapeInfoValue(s string) string {
	r := strings.NewReplacer(" ", "_", ";", "_", "=", "_", ",", "_")
	return r.Replace(s)
}
