package ensembl

import (
	"context"
	"fmt"
	"strings"
)

// VEPInput is one variant for batch annotation. Fields are 1-based inclusive
// GRCh38 coordinates matching VCF conventions.
type VEPInput struct {
	Chrom string
	Pos   int
	Ref   string
	Alt   string
}

// Key returns a stable cache key for this variant.
func (v VEPInput) Key() string {
	return fmt.Sprintf("%s:%d:%s>%s", v.Chrom, v.Pos, v.Ref, v.Alt)
}

// vcfInput encodes a VEP input as a VCF data line ("1 100 . A T . . ."). VCF
// format keeps REF/ALT exactly as the source VCF has them, including the anchor
// base on indels, so VEP normalises them itself. Ensembl names the mitochondrial
// sequence MT rather than M.
func (v VEPInput) vcfInput() string {
	chrom := strings.TrimPrefix(v.Chrom, "chr")
	if chrom == "M" {
		chrom = "MT"
	}
	return fmt.Sprintf("%s %d . %s %s . . .", chrom, v.Pos, v.Ref, v.Alt)
}

// vepAllele returns ALT as VEP reports it after trimming the anchor base that
// VCF indels share with REF: CTT>C becomes "-", A>AG becomes "G". SNVs are
// returned unchanged.
func vepAllele(ref, alt string) string {
	if ref == "" || alt == "" || ref[0] != alt[0] {
		return alt
	}
	if trimmed := alt[1:]; trimmed != "" {
		return trimmed
	}
	return "-"
}

// Annotation is the subset of VEP output vcfq surfaces. AF is gnomAD genomes
// where present, falling back to 1000G overall AF.
type Annotation struct {
	Consequence string  `json:"consequence"`
	Gene        string  `json:"gene"`
	AAChange    string  `json:"aa_change"`
	AF          float64 `json:"af"`
}

// rawVEPRecord matches the top level of an item in /vep/human/region's response.
type rawVEPRecord struct {
	Input               string                 `json:"input"`
	MostSevere          string                 `json:"most_severe_consequence"`
	TranscriptConsequen []rawTranscriptConsequ `json:"transcript_consequences"`
	ColocatedVariants   []rawColocated         `json:"colocated_variants"`
}

type rawTranscriptConsequ struct {
	GeneSymbol      string   `json:"gene_symbol"`
	Consequences    []string `json:"consequence_terms"`
	AminoAcids      string   `json:"amino_acids"`
	BiotypeCanonica int      `json:"canonical"`
}

// rawColocated.Frequencies is keyed by ALT allele then population:
//
//	{ "A": { "gnomadg": 0.275, "af": 0.245, ... }, ... }
type rawColocated struct {
	Frequencies map[string]map[string]float64 `json:"frequencies"`
}

// vepCacheKind names the cache directory for VEP results. Bumped from "vep"
// when inputs switched to VCF format: entries cached under the old encoding
// may hold wrong annotations for indels, so they're never read again.
const vepCacheKind = "vep2"

// AnnotateBatch sends up to 200 variants per request and merges results back in
// the input order. Variants found in cache are not re-fetched. Variants for
// which VEP returns no consequence get a zero-value Annotation.
func (c *Client) AnnotateBatch(ctx context.Context, cache *Cache, inputs []VEPInput) ([]Annotation, error) {
	out := make([]Annotation, len(inputs))
	missingIdx := []int{}
	missingInputs := []VEPInput{}

	for i, v := range inputs {
		if cache != nil {
			var hit Annotation
			ok, err := cache.Get(vepCacheKind, v.Key(), &hit)
			if err == nil && ok {
				out[i] = hit
				continue
			}
		}
		missingIdx = append(missingIdx, i)
		missingInputs = append(missingInputs, v)
	}

	const chunkSize = 200
	for start := 0; start < len(missingInputs); start += chunkSize {
		end := start + chunkSize
		if end > len(missingInputs) {
			end = len(missingInputs)
		}
		chunk := missingInputs[start:end]
		variants := make([]string, len(chunk))
		for j, v := range chunk {
			variants[j] = v.vcfInput()
		}
		var recs []rawVEPRecord
		body := map[string]any{"variants": variants}
		// af=1 enables 1000 Genomes frequencies; af_gnomadg=1 enables gnomAD
		// genomes frequencies. Without these, frequencies come back empty.
		if err := c.postJSON(ctx, "/vep/human/region?af=1&af_gnomadg=1", body, &recs); err != nil {
			return nil, err
		}
		// VEP response order isn't guaranteed to match input order, so match by
		// the "input" string it echoes back.
		byInput := make(map[string]rawVEPRecord, len(recs))
		for _, r := range recs {
			byInput[normaliseInput(r.Input)] = r
		}
		for j, v := range chunk {
			i := missingIdx[start+j]
			rec, ok := byInput[normaliseInput(v.vcfInput())]
			if !ok {
				continue // VEP returned nothing for this variant
			}
			ann := buildAnnotation(rec, vepAllele(v.Ref, v.Alt))
			out[i] = ann
			if cache != nil {
				_ = cache.Set(vepCacheKind, v.Key(), ann)
			}
		}
	}
	return out, nil
}

// normaliseInput collapses whitespace differences between what we send and what
// VEP echoes back ("1 100 100 A/T 1" vs tab-separated etc).
func normaliseInput(s string) string {
	return strings.Join(strings.Fields(s), " ")
}

// buildAnnotation picks the canonical-transcript consequence when available,
// falling back to the most-severe consequence for the variant overall. alt is
// the ALT allele in VEP's trimmed form (see vepAllele), used to pick the right
// per-allele frequency block.
func buildAnnotation(r rawVEPRecord, alt string) Annotation {
	a := Annotation{Consequence: r.MostSevere}
	for _, tc := range r.TranscriptConsequen {
		if tc.BiotypeCanonica == 1 {
			if len(tc.Consequences) > 0 {
				a.Consequence = tc.Consequences[0]
			}
			a.Gene = tc.GeneSymbol
			a.AAChange = tc.AminoAcids
			break
		}
	}
	// If no canonical transcript flagged, take the first transcript consequence
	// for gene/AA so we don't silently drop them.
	if a.Gene == "" && len(r.TranscriptConsequen) > 0 {
		tc := r.TranscriptConsequen[0]
		a.Gene = tc.GeneSymbol
		a.AAChange = tc.AminoAcids
	}
	a.AF = pickAF(r.ColocatedVariants, alt)
	return a
}

// pickAF chooses an allele frequency from VEP's colocated_variants. The
// frequencies map is keyed by ALT allele (e.g. "A"), and inside that by
// population key. Preference: gnomAD genomes (gnomadg), gnomAD exomes (gnomade),
// 1000G overall (af). If alt is empty we take the first allele found.
func pickAF(cvs []rawColocated, alt string) float64 {
	preferred := []string{"gnomadg", "gnomade", "af"}
	for _, cv := range cvs {
		freqsForAlt, ok := cv.Frequencies[alt]
		if !ok && alt != "" {
			continue
		}
		if !ok {
			for _, f := range cv.Frequencies {
				freqsForAlt = f
				break
			}
		}
		for _, k := range preferred {
			if v, ok := freqsForAlt[k]; ok {
				return v
			}
		}
	}
	return 0
}
