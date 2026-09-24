package vcf

import (
	"os"
	"path/filepath"
	"slices"
	"testing"
)

// testdata/sample.vcf is the readable source of the fixture. The .gz and .tbi
// beside it were produced from it with htslib (bgzip + tabix -p vcf); rebuild
// them the same way after editing it.
const fixture = "testdata/sample.vcf.gz"

func openFixture(t *testing.T) *Reader {
	t.Helper()
	r, err := Open(fixture)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	t.Cleanup(func() {
		if err := r.Close(); err != nil {
			t.Errorf("Close: %v", err)
		}
	})
	return r
}

func TestOpenReadsHeader(t *testing.T) {
	r := openFixture(t)
	if len(r.Header) != 7 {
		t.Fatalf("got %d header lines, want 7", len(r.Header))
	}
	if r.Header[0] != "##fileformat=VCFv4.2" || r.Header[6][:6] != "#CHROM" {
		t.Errorf("unexpected header bounds: %q ... %q", r.Header[0], r.Header[6])
	}
}

func TestScanRegion(t *testing.T) {
	r := openFixture(t)
	cases := []struct {
		chrom      string
		start, end int
		wantPos    []int
	}{
		{"chr1", 11796321, 11796321, []int{11796321}},
		{"1", 11796321, 11796321, []int{11796321}}, // bare name resolves to chr1
		{"chr1", 20000100, 20000400, []int{20000100, 20000200, 20000300, 20000400}},
		{"chr1", 20000101, 20000299, []int{20000200}},
		{"chr1", 1, 1000, nil},
		{"chr2", 100, 100, []int{100}},
	}
	for _, c := range cases {
		lines, err := r.ScanRegion(c.chrom, c.start, c.end)
		if err != nil {
			t.Fatalf("ScanRegion(%s:%d-%d): %v", c.chrom, c.start, c.end, err)
		}
		var got []int
		for _, l := range lines {
			got = append(got, l.Pos)
		}
		if !slices.Equal(got, c.wantPos) {
			t.Errorf("ScanRegion(%s:%d-%d) positions = %v, want %v", c.chrom, c.start, c.end, got, c.wantPos)
		}
	}
}

func TestScanRegionParsesColumns(t *testing.T) {
	r := openFixture(t)
	lines, err := r.ScanRegion("chr1", 20000300, 20000300)
	if err != nil || len(lines) != 1 {
		t.Fatalf("got %d lines, err %v; want 1", len(lines), err)
	}
	d := lines[0]
	if d.Ref != "G" || !slices.Equal(d.Alts, []string{"A", "T"}) || d.Info != "AC=1,1" ||
		d.Format != "GT:DP" || d.Sample != "1/2:25" {
		t.Errorf("unexpected columns: %+v", d)
	}
}

func TestScanRegionUnknownChrom(t *testing.T) {
	r := openFixture(t)
	if _, err := r.ScanRegion("chrX", 1, 100); err == nil {
		t.Error("expected an error for a chrom absent from the index")
	}
}

func TestOpenMissingIndex(t *testing.T) {
	b, err := os.ReadFile(fixture)
	if err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "noindex.vcf.gz")
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := Open(path); err == nil {
		t.Error("expected an error when the .tbi is missing")
	}
}
