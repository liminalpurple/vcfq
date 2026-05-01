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
		rec.Ref,
		rec.Alt,
		emptyDot(rec.GT),
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
// annotations are enabled. Order is stable for reproducible output.
var annotationInfoLines = []string{
	`##INFO=<ID=VCFQ_CONSEQUENCE,Number=1,Type=String,Description="VEP consequence term">`,
	`##INFO=<ID=VCFQ_GENE,Number=1,Type=String,Description="VEP-resolved gene symbol">`,
	`##INFO=<ID=VCFQ_AA,Number=1,Type=String,Description="Amino acid change (ref/alt)">`,
	`##INFO=<ID=VCFQ_AF,Number=1,Type=Float,Description="Allele frequency from VEP (gnomAD/1000G)">`,
}

type vcfFormatter struct {
	w         io.Writer
	annotated bool
}

func (f *vcfFormatter) Header(vcfHeader []string, annotated bool) error {
	f.annotated = annotated
	// Pass through header lines, injecting our INFO declarations just before the
	// final #CHROM... line. If the source has no header (empty slice), emit a
	// minimal one.
	if len(vcfHeader) == 0 {
		if _, err := fmt.Fprintln(f.w, "##fileformat=VCFv4.2"); err != nil {
			return err
		}
		if annotated {
			for _, line := range annotationInfoLines {
				if _, err := fmt.Fprintln(f.w, line); err != nil {
					return err
				}
			}
		}
		_, err := fmt.Fprintln(f.w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
		return err
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
	info := rec.Info
	if f.annotated && rec.HasAnnot {
		info = appendAnnotations(info, rec)
	}
	cols := []string{
		rec.Chrom,
		strconv.Itoa(rec.Pos),
		emptyDot(rec.ID),
		rec.Ref,
		rec.Alt,
		emptyDot(rec.Qual),
		emptyDot(rec.Filter),
		emptyDot(info),
	}
	if rec.Format != "" {
		cols = append(cols, rec.Format)
		if rec.Sample != "" {
			cols = append(cols, rec.Sample)
		}
	}
	_, err := fmt.Fprintln(f.w, strings.Join(cols, "\t"))
	return err
}

func (f *vcfFormatter) Close() error { return nil }

// appendAnnotations injects VCFQ_* INFO fields into an existing INFO string. If
// the source INFO is empty or ".", they replace it; otherwise they're appended
// with ";" separators.
func appendAnnotations(info string, rec Record) string {
	parts := []string{}
	if rec.Consequence != "" {
		parts = append(parts, "VCFQ_CONSEQUENCE="+escapeInfoValue(rec.Consequence))
	}
	if rec.Gene != "" {
		parts = append(parts, "VCFQ_GENE="+escapeInfoValue(rec.Gene))
	}
	if rec.AAChange != "" {
		parts = append(parts, "VCFQ_AA="+escapeInfoValue(rec.AAChange))
	}
	if rec.AF > 0 {
		parts = append(parts, "VCFQ_AF="+strconv.FormatFloat(rec.AF, 'g', 4, 64))
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
