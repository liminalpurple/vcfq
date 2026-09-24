package cli

import (
	"bytes"
	"io"
	"slices"
	"strings"
	"testing"

	"github.com/liminalpurple/vcfq/internal/vcf"
)

// runVCFQ runs the CLI against the fixture and returns stdout. Region queries
// never touch Ensembl, so these tests are fully offline.
func runVCFQ(t *testing.T, args ...string) string {
	t.Helper()
	var stdout, stderr bytes.Buffer
	args = append([]string{"-vcf", fixtureVCF, "-cache", t.TempDir()}, args...)
	if code := Run(args, strings.NewReader(""), &stdout, &stderr); code != 0 {
		t.Fatalf("Run(%q) exited %d: %s", args, code, stderr.String())
	}
	return stdout.String()
}

func dataLines(out string) []string {
	var lines []string
	for _, l := range strings.Split(strings.TrimSuffix(out, "\n"), "\n") {
		if l != "" && !strings.HasPrefix(l, "#") {
			lines = append(lines, l)
		}
	}
	return lines
}

func TestTSVEmitsEveryCalledAlt(t *testing.T) {
	got := dataLines(runVCFQ(t, "chr1:20000300-20000400"))
	want := []string{
		"chrom\tpos\tid\tref\talt\tgt\tquery",
		"chr1\t20000300\t.\tG\tA\t1/2\tchr1:20000300-20000400", // 1/2 het: both ALTs
		"chr1\t20000300\t.\tG\tT\t1/2\tchr1:20000300-20000400",
		"chr1\t20000400\t.\tC\tA\t0/2\tchr1:20000300-20000400", // only ALT 2 called
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("tsv output:\n%s\nwant:\n%s", strings.Join(got, "\n"), strings.Join(want, "\n"))
	}
}

func TestVCFPreservesMultiAllelicLines(t *testing.T) {
	out := runVCFQ(t, "-f", "vcf", "chr1:20000300-20000400")
	if !strings.HasPrefix(out, "##fileformat=VCFv4.2\n") {
		t.Errorf("source header not passed through:\n%s", out)
	}
	got := dataLines(out)
	want := []string{
		"chr1\t20000300\t.\tG\tA,T\t50\tPASS\tAC=1,1\tGT:DP\t1/2:25",
		"chr1\t20000400\t.\tC\tG,A\t50\tPASS\tAC=0,1\tGT:DP\t0/2:27",
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("vcf output:\n%s\nwant:\n%s", strings.Join(got, "\n"), strings.Join(want, "\n"))
	}
}

func TestVCFAnnotationsAreNumberA(t *testing.T) {
	line := &vcf.DataLine{
		Chrom: "chr1", Pos: 100, ID: ".", Ref: "G", Alts: []string{"A", "T", "C"},
		Qual: "50", Filter: "PASS", Info: "AC=1,1,0", Format: "GT", Sample: "1/2",
	}
	recs := []Record{
		{Line: line, AltIndex: 1, HasAnnot: true, Gene: "GENE1", Consequence: "missense variant", AF: 0.25},
		{Line: line, AltIndex: 2, HasAnnot: true, Gene: "GENE1", Consequence: "synonymous_variant"},
	}
	var buf bytes.Buffer
	f := &vcfFormatter{w: &buf, annotated: true}
	for _, r := range recs {
		if err := f.Write(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	want := "chr1\t100\t.\tG\tA,T,C\t50\tPASS\t" +
		"AC=1,1,0;VCFQ_CONSEQUENCE=missense_variant,synonymous_variant,.;VCFQ_GENE=GENE1,GENE1,.;VCFQ_AF=0.25,.,.\tGT\t1/2\n"
	if buf.String() != want {
		t.Errorf("got  %q\nwant %q", buf.String(), want)
	}
}

func TestVCFHeaderDeclaresAnnotations(t *testing.T) {
	var buf bytes.Buffer
	f := &vcfFormatter{w: &buf}
	header := []string{"##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}
	if err := f.Header(header, true); err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSuffix(buf.String(), "\n"), "\n")
	if len(lines) != 6 || !strings.HasPrefix(lines[5], "#CHROM") {
		t.Fatalf("INFO declarations should sit between the header and #CHROM:\n%s", buf.String())
	}
	for _, l := range lines[1:5] {
		if !strings.Contains(l, "Number=A") {
			t.Errorf("annotation INFO line not Number=A: %s", l)
		}
	}
}

func TestRegionNoCallEndToEnd(t *testing.T) {
	got := dataLines(runVCFQ(t, "chr1:11796320", "chr1:1-1000"))
	want := []string{
		"chrom\tpos\tid\tref\talt\tgt\tquery",
		"chr1\t11796320\t.\t.\t.\tno-call\tchr1:11796320", // empty range stays empty
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("output:\n%s\nwant:\n%s", strings.Join(got, "\n"), strings.Join(want, "\n"))
	}
	if out := runVCFQ(t, "-f", "vcf", "chr1:11796320"); len(dataLines(out)) != 0 {
		t.Errorf("-f vcf should omit no-call rows, got:\n%s", out)
	}
}

// forbiddenStdin fails the test if read: it stands in for a stdin inherited
// from a script or cron job that is open but will never reach EOF.
type forbiddenStdin struct{ t *testing.T }

func (f forbiddenStdin) Read([]byte) (int, error) {
	f.t.Error("stdin was read although queries were given as arguments")
	return 0, io.EOF
}

func runWithStdin(t *testing.T, stdin io.Reader, args ...string) string {
	t.Helper()
	var stdout, stderr bytes.Buffer
	args = append([]string{"-vcf", fixtureVCF, "-cache", t.TempDir()}, args...)
	if code := Run(args, stdin, &stdout, &stderr); code != 0 {
		t.Fatalf("Run(%q) exited %d: %s", args, code, stderr.String())
	}
	return stdout.String()
}

func queryColumn(out string) []string {
	var qs []string
	for _, l := range dataLines(out)[1:] {
		f := strings.Split(l, "\t")
		qs = append(qs, f[len(f)-1])
	}
	return qs
}

func TestStdinOnlyReadWhenAsked(t *testing.T) {
	t.Run("arguments leave stdin alone", func(t *testing.T) {
		got := queryColumn(runWithStdin(t, forbiddenStdin{t}, "chr2:100"))
		if !slices.Equal(got, []string{"chr2:100"}) {
			t.Errorf("queries = %v", got)
		}
	})
	t.Run("no arguments reads stdin", func(t *testing.T) {
		got := queryColumn(runWithStdin(t, strings.NewReader("chr2:100\nchr1:11796321\n")))
		if !slices.Equal(got, []string{"chr2:100", "chr1:11796321"}) {
			t.Errorf("queries = %v", got)
		}
	})
	t.Run("dash splices stdin in place", func(t *testing.T) {
		got := queryColumn(runWithStdin(t, strings.NewReader("chr2:100"), "chr1:11796321", "-", "chr1:20000100"))
		if !slices.Equal(got, []string{"chr1:11796321", "chr2:100", "chr1:20000100"}) {
			t.Errorf("queries = %v", got)
		}
	})
}
