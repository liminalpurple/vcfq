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

	HasAnnot    bool
	Consequence string
	Gene        string
	AAChange    string
	AF          float64
}
