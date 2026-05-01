package cli

import "testing"

func TestDetect(t *testing.T) {
	cases := []struct {
		in   string
		want QueryKind
	}{
		{"rs1801133", KindRSID},
		{"RS1801133", KindRSID},
		{"rs53576", KindRSID},
		{"chr1:11796321", KindRegion},
		{"chr1:11796000-11796500", KindRegion},
		{"1:11796000-11796500", KindRegion},
		{"chrX:1234-5678", KindRegion},
		{"chrM:100-200", KindRegion},
		{"HNF1A", KindSymbol},
		{"MTHFR", KindSymbol},
		{"BRCA1", KindSymbol},
		{"TP53", KindSymbol},
		// rsfoo is not a valid rsID (must be digits after rs)
		{"rsfoo", KindSymbol},
		// missing colon -> symbol
		{"chr1", KindSymbol},
	}
	for _, c := range cases {
		if got := Detect(c.in); got != c.want {
			t.Errorf("Detect(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}
