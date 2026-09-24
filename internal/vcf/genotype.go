package vcf

import (
	"slices"
	"strconv"
	"strings"
)

// ParseGT extracts the GT subfield from a VCF FORMAT/sample pair. Returns ""
// if the format doesn't include GT or the sample column is missing values for
// it. The raw genotype is returned verbatim ("0/1", "1|2", "1/1", etc).
func ParseGT(format, sample string) string {
	if format == "" || sample == "" {
		return ""
	}
	keys := strings.Split(format, ":")
	values := strings.Split(sample, ":")
	for i, k := range keys {
		if k == "GT" {
			if i >= len(values) {
				return ""
			}
			return values[i]
		}
	}
	return ""
}

// CalledAlts returns the distinct 1-based ALT indices a genotype calls, in
// ascending order: "1/2" gives [1 2], "0/1" gives [1]. It returns nil for
// homozygous reference, missing, or unparseable genotypes.
func CalledAlts(gt string) []int {
	var out []int
	for _, a := range strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' }) {
		if a == "." {
			continue
		}
		n, err := strconv.Atoi(a)
		if err != nil || n < 0 {
			return nil
		}
		if n > 0 && !slices.Contains(out, n) {
			out = append(out, n)
		}
	}
	slices.Sort(out)
	return out
}
