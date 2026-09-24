// Package cli implements the vcfq command-line interface: query auto-detection,
// stdin pipe handling, and output formatters.
package cli

import "github.com/liminalpurple/vcfq/internal/vcf"

// Record is one called ALT allele produced by a query. A multi-allelic VCF line
// yields one Record per called ALT; all of them share the same Line, which the
// vcf formatter uses to reassemble the original line with annotations injected.
type Record struct {
	Chrom string
	Pos   int
	ID    string
	Ref   string
	Alt   string

	// Line is the source VCF line, and AltIndex is Alt's 1-based position in
	// Line.Alts.
	Line     *vcf.DataLine
	AltIndex int

	GT    string
	Query string

	// NoCall marks a synthetic row: the query resolved to a single position but
	// the VCF has no line there. Ref/Alt on these rows come from Ensembl rather
	// than the VCF, and the absence itself is ambiguous — it may be homozygous
	// reference, or the site may simply not have been callable. vcfq cannot tell
	// the two apart without read-level data, so it reports the absence and makes
	// no claim about the genotype.
	NoCall bool

	HasAnnot    bool
	Consequence string
	Gene        string
	AAChange    string
	AF          float64
}
