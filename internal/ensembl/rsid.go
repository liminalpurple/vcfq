package ensembl

import (
	"context"
	"fmt"
	"strings"
)

// VariantLoc is a GRCh38 mapping for an rsID. Allele strings are slash-separated
// in Ensembl's response (e.g. "C/T", "TC/T"); for vcfq's purposes Ref is the
// allele matching the reference and Alts contains every non-ref allele.
type VariantLoc struct {
	ID    string   `json:"id"`
	Chrom string   `json:"chrom"`
	Start int      `json:"start"`
	End   int      `json:"end"`
	Ref   string   `json:"ref"`
	Alts  []string `json:"alts"`
}

// rawRSID matches /variation/human/{rsID}.
type rawRSID struct {
	Name     string       `json:"name"`
	Mappings []rawMapping `json:"mappings"`
}

type rawMapping struct {
	AssemblyName  string `json:"assembly_name"`
	SeqRegionName string `json:"seq_region_name"`
	Start         int    `json:"start"`
	End           int    `json:"end"`
	Allele        string `json:"allele_string"` // "C/T" or "C/T/G"
	AncestralAll  string `json:"ancestral_allele"`
}

// LookupRSID resolves an rsID to its GRCh38 location, ref, and alt allele(s).
// Mappings to other assemblies are filtered out. Returns ErrNotFound for
// unknown rsIDs.
func (c *Client) LookupRSID(ctx context.Context, cache *Cache, rsid string) (*VariantLoc, error) {
	key := strings.ToLower(rsid)
	if cache != nil {
		var hit VariantLoc
		ok, err := cache.Get("rsid", key, &hit)
		if err == nil && ok {
			return &hit, nil
		}
	}
	var raw rawRSID
	if err := c.getJSON(ctx, "/variation/human/"+key, &raw); err != nil {
		return nil, err
	}
	for _, m := range raw.Mappings {
		if m.AssemblyName != "GRCh38" {
			continue
		}
		alleles := strings.Split(m.Allele, "/")
		if len(alleles) < 2 {
			continue
		}
		loc := &VariantLoc{
			ID:    raw.Name,
			Chrom: normaliseChrom(m.SeqRegionName),
			Start: m.Start,
			End:   m.End,
			Ref:   alleles[0],
			Alts:  alleles[1:],
		}
		if cache != nil {
			_ = cache.Set("rsid", key, loc)
		}
		return loc, nil
	}
	return nil, fmt.Errorf("%w: no GRCh38 mapping for %s", ErrNotFound, rsid)
}
