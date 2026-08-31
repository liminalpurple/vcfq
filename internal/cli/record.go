// Package cli implements the vcfq command-line interface: query auto-detection,
// stdin pipe handling, and output formatters.
package cli

// Record is one variant row produced by a query. VCF columns are kept separately
// so the vcf formatter can reassemble lines with annotations injected into INFO.
type Record struct {
	Chrom  string
	Pos    int
	ID     string
	Ref    string
	Alt    string
	Qual   string
	Filter string
	Info   string
	Format string
	Sample string

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
