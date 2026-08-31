package cli

import (
	"bytes"
	"strings"
	"testing"
)

func TestNoCallRecord(t *testing.T) {
	cases := []struct {
		name     string
		rec      Record
		wantAlt  string
		wantPos  int
		wantNoCa bool
	}{
		{
			name:     "biallelic rsID",
			rec:      noCallRecord("chr6", 32638107, "rs2187668", "C", []string{"T"}, "rs2187668"),
			wantAlt:  "T",
			wantPos:  32638107,
			wantNoCa: true,
		},
		{
			name:     "multiallelic rsID joins alts",
			rec:      noCallRecord("chr6", 32690302, "rs7775228", "T", []string{"A", "C"}, "rs7775228"),
			wantAlt:  "A,C",
			wantPos:  32690302,
			wantNoCa: true,
		},
		{
			name:     "single-position region has no alleles to report",
			rec:      noCallRecord("chr6", 32638107, "", "", nil, "chr6:32638107"),
			wantAlt:  "",
			wantPos:  32638107,
			wantNoCa: true,
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if c.rec.Alt != c.wantAlt {
				t.Errorf("Alt = %q, want %q", c.rec.Alt, c.wantAlt)
			}
			if c.rec.Pos != c.wantPos {
				t.Errorf("Pos = %d, want %d", c.rec.Pos, c.wantPos)
			}
			if c.rec.NoCall != c.wantNoCa {
				t.Errorf("NoCall = %v, want %v", c.rec.NoCall, c.wantNoCa)
			}
			if c.rec.GT != "" {
				t.Errorf("GT = %q, want empty: a no-call record must not assert a genotype", c.rec.GT)
			}
		})
	}
}

func TestGTDisplay(t *testing.T) {
	cases := []struct {
		name string
		rec  Record
		want string
	}{
		{"called het", Record{GT: "0/1"}, "0/1"},
		{"called hom alt", Record{GT: "1/1"}, "1/1"},
		{"missing GT", Record{GT: ""}, "."},
		{"no-call", Record{NoCall: true}, "no-call"},
		// NoCall is the stronger claim: it must win even if a GT is somehow set,
		// so a stray "0/0" can never be reported as an observed genotype.
		{"no-call outranks a stray GT", Record{GT: "0/0", NoCall: true}, "no-call"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := gtDisplay(c.rec); got != c.want {
				t.Errorf("gtDisplay() = %q, want %q", got, c.want)
			}
		})
	}
}

func TestAnyNoCall(t *testing.T) {
	cases := []struct {
		name string
		recs []Record
		want bool
	}{
		{"nil", nil, false},
		{"all called", []Record{{GT: "0/1"}, {GT: "1/1"}}, false},
		{"one no-call among calls", []Record{{GT: "0/1"}, {NoCall: true}}, true},
		{"only no-call", []Record{{NoCall: true}}, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := anyNoCall(c.recs); got != c.want {
				t.Errorf("anyNoCall() = %v, want %v", got, c.want)
			}
		})
	}
}

// noCallFixture is the row today's HLA work would have produced: rs2187668
// resolves to a real GRCh38 position, and Morgan's VCF has no line there.
func noCallFixture() Record {
	return noCallRecord("chr6", 32638107, "rs2187668", "C", []string{"T"}, "rs2187668")
}

func TestFormatterNoCallTSV(t *testing.T) {
	var buf bytes.Buffer
	f := &tsvFormatter{w: &buf}
	if err := f.Header(nil, false); err != nil {
		t.Fatal(err)
	}
	if err := f.Write(noCallFixture()); err != nil {
		t.Fatal(err)
	}
	got := buf.String()
	want := "chr6\t32638107\trs2187668\tC\tT\tno-call\trs2187668\n"
	if !strings.HasSuffix(got, want) {
		t.Errorf("tsv row =\n%q\nwant suffix\n%q", got, want)
	}
	if strings.Contains(got, "0/0") {
		t.Error("tsv row reports a 0/0 genotype; absence must never be rendered as hom-ref")
	}
}

// A region no-call has no alleles to report; the tsv/table convention is "." for
// missing, not an empty cell.
func TestFormatterNoCallRegionTSVUsesDots(t *testing.T) {
	var buf bytes.Buffer
	f := &tsvFormatter{w: &buf}
	if err := f.Header(nil, false); err != nil {
		t.Fatal(err)
	}
	if err := f.Write(noCallRecord("chr6", 32638107, "", "", nil, "chr6:32638107")); err != nil {
		t.Fatal(err)
	}
	want := "chr6\t32638107\t.\t.\t.\tno-call\tchr6:32638107\n"
	if !strings.HasSuffix(buf.String(), want) {
		t.Errorf("tsv row =\n%q\nwant suffix\n%q", buf.String(), want)
	}
}

func TestFormatterNoCallJSON(t *testing.T) {
	var buf bytes.Buffer
	f := &jsonFormatter{w: &buf}
	if err := f.Header(nil, false); err != nil {
		t.Fatal(err)
	}
	if err := f.Write(noCallFixture()); err != nil {
		t.Fatal(err)
	}
	got := buf.String()
	for _, want := range []string{`"gt":"no-call"`, `"source":"ensembl"`, `"pos":32638107`} {
		if !strings.Contains(got, want) {
			t.Errorf("json row %s missing %s", got, want)
		}
	}
}

func TestFormatterJSONOmitsSourceWhenCalled(t *testing.T) {
	var buf bytes.Buffer
	f := &jsonFormatter{w: &buf}
	if err := f.Header(nil, false); err != nil {
		t.Fatal(err)
	}
	rec := Record{Chrom: "chr6", Pos: 32713706, ID: "rs7454108", Ref: "T", Alt: "C", GT: "0/1", Query: "rs7454108"}
	if err := f.Write(rec); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(buf.String(), "source") {
		t.Errorf("json row for a called variant carries a source field: %s", buf.String())
	}
}

func TestFormatterNoCallVCFSkipped(t *testing.T) {
	var buf bytes.Buffer
	f := &vcfFormatter{w: &buf}
	if err := f.Write(noCallFixture()); err != nil {
		t.Fatal(err)
	}
	if buf.Len() != 0 {
		t.Errorf("vcf formatter emitted a data line for a no-call: %q", buf.String())
	}
}
