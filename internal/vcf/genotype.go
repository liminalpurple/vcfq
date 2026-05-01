package vcf

import "strings"

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

// AltIndex returns the 1-based ALT allele index that a genotype calls, or 0 if
// the genotype is homozygous reference, or -1 if the genotype is missing or
// can't be parsed unambiguously. For diploid genotypes, AltIndex returns the
// highest non-zero allele index — sufficient for vcfq's purpose of selecting
// which ALT to display when a record has multiple ALTs.
func AltIndex(gt string) int {
	if gt == "" || gt == "." || gt == "./." || gt == ".|." {
		return -1
	}
	sep := "/"
	if strings.Contains(gt, "|") {
		sep = "|"
	}
	max := 0
	for _, a := range strings.Split(gt, sep) {
		if a == "." {
			continue
		}
		n := 0
		for _, c := range a {
			if c < '0' || c > '9' {
				return -1
			}
			n = n*10 + int(c-'0')
		}
		if n > max {
			max = n
		}
	}
	return max
}
