package ensembl

import (
	"context"
	"strings"
)

// SymbolInfo is the GRCh38 region of a gene. Strand is +1 or -1.
type SymbolInfo struct {
	Symbol string `json:"symbol"`
	Chrom  string `json:"chrom"`
	Start  int    `json:"start"`
	End    int    `json:"end"`
	Strand int    `json:"strand"`
}

// rawSymbol matches the fields vcfq cares about from /lookup/symbol.
type rawSymbol struct {
	SeqRegionName string `json:"seq_region_name"`
	Start         int    `json:"start"`
	End           int    `json:"end"`
	Strand        int    `json:"strand"`
	Display       string `json:"display_name"`
}

// LookupSymbol resolves a HGNC gene symbol to GRCh38 coordinates. Cached on disk
// after the first hit. Returns ErrNotFound for unknown symbols.
func (c *Client) LookupSymbol(ctx context.Context, cache *Cache, symbol string) (*SymbolInfo, error) {
	key := strings.ToUpper(symbol)
	if cache != nil {
		var hit SymbolInfo
		ok, err := cache.Get("symbol", key, &hit)
		if err == nil && ok {
			return &hit, nil
		}
	}
	var raw rawSymbol
	if err := c.getJSON(ctx, "/lookup/symbol/homo_sapiens/"+key, &raw); err != nil {
		return nil, err
	}
	info := &SymbolInfo{
		Symbol: key,
		Chrom:  normaliseChrom(raw.SeqRegionName),
		Start:  raw.Start,
		End:    raw.End,
		Strand: raw.Strand,
	}
	if cache != nil {
		_ = cache.Set("symbol", key, info)
	}
	return info, nil
}

// normaliseChrom returns the chr-prefixed form ("1" -> "chr1") to match the VCF.
// MT in Ensembl maps to chrM in GRCh38 VCFs from most aligners (including
// Nebula's mm2 pipeline).
func normaliseChrom(s string) string {
	if strings.HasPrefix(s, "chr") {
		return s
	}
	if s == "MT" {
		return "chrM"
	}
	return "chr" + s
}
