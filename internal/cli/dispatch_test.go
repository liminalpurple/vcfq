package cli

import (
	"context"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"sync"
	"testing"

	"github.com/liminalpurple/vcfq/internal/ensembl"
	"github.com/liminalpurple/vcfq/internal/vcf"
)

const fixtureVCF = "../vcf/testdata/sample.vcf.gz"

// fakeEnsembl serves canned /variation responses and echoes /vep requests back
// with a fixed annotation, recording every VEP input line it receives.
type fakeEnsembl struct {
	mu        sync.Mutex
	vepInputs []string
	mappings  map[string]map[string]any
}

func newFakeEnsembl(t *testing.T) (*fakeEnsembl, *ensembl.Client) {
	t.Helper()
	fe := &fakeEnsembl{
		// Ensembl reports indels without VCF's anchor base: the deletion covers
		// only the deleted bases, the insertion has start = end+1.
		mappings: map[string]map[string]any{
			"rs900001": {"start": 20000101, "end": 20000102, "allele_string": "TT/-"},
			"rs900002": {"start": 20000201, "end": 20000200, "allele_string": "-/G"},
			// Absent from the fixture, but its anchor base carries rs1801133.
			"rs900003": {"start": 11796322, "end": 11796323, "allele_string": "GC/-"},
		},
	}
	mux := http.NewServeMux()
	mux.HandleFunc("GET /variation/human/{id}", func(w http.ResponseWriter, r *http.Request) {
		id := r.PathValue("id")
		m, ok := fe.mappings[id]
		if !ok {
			http.NotFound(w, r)
			return
		}
		mapping := map[string]any{"assembly_name": "GRCh38", "seq_region_name": "1"}
		for k, v := range m {
			mapping[k] = v
		}
		_ = json.NewEncoder(w).Encode(map[string]any{"name": id, "mappings": []any{mapping}})
	})
	mux.HandleFunc("POST /vep/human/region", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Variants []string `json:"variants"`
		}
		if err := json.NewDecoder(r.Body).Decode(&body); err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		fe.mu.Lock()
		fe.vepInputs = append(fe.vepInputs, body.Variants...)
		fe.mu.Unlock()
		out := make([]any, 0, len(body.Variants))
		for _, v := range body.Variants {
			// VEP trims the shared anchor base from VCF input, so the deletion
			// CTT>C is keyed "-" in colocated frequencies.
			out = append(out, map[string]any{
				"input":                   v,
				"most_severe_consequence": "intron_variant",
				"transcript_consequences": []any{map[string]any{
					"gene_symbol": "FAKE", "consequence_terms": []string{"intron_variant"}, "canonical": 1,
				}},
				"colocated_variants": []any{map[string]any{
					"frequencies": map[string]any{"-": map[string]float64{"gnomadg": 0.125}},
				}},
			})
		}
		_ = json.NewEncoder(w).Encode(out)
	})
	srv := httptest.NewServer(mux)
	t.Cleanup(srv.Close)
	return fe, &ensembl.Client{BaseURL: srv.URL, HTTP: srv.Client(), UserAgent: "vcfq-test"}
}

func openFixture(t *testing.T) *vcf.Reader {
	t.Helper()
	r, err := vcf.Open(fixtureVCF)
	if err != nil {
		t.Fatalf("open fixture: %v", err)
	}
	t.Cleanup(func() {
		if err := r.Close(); err != nil {
			t.Errorf("close fixture: %v", err)
		}
	})
	return r
}

func TestResolveRSIDIndels(t *testing.T) {
	_, client := newFakeEnsembl(t)
	reader := openFixture(t)
	cases := []struct {
		rsid    string
		wantPos int
		wantAlt string
	}{
		{"rs900001", 20000100, "C"},  // deletion, matched by ID
		{"rs900002", 20000200, "AG"}, // insertion, VCF has no ID so matched by window
	}
	for _, c := range cases {
		cache, _ := ensembl.NewCache(t.TempDir())
		recs, err := resolveRSID(context.Background(), client, cache, reader, c.rsid)
		if err != nil {
			t.Fatalf("%s: %v", c.rsid, err)
		}
		if len(recs) != 1 || recs[0].Pos != c.wantPos || recs[0].Alt != c.wantAlt {
			t.Errorf("%s: got %+v, want one record at %d with ALT %s", c.rsid, recs, c.wantPos, c.wantAlt)
		}
	}
}

func TestAnnotateIndelUsesVCFInput(t *testing.T) {
	fe, client := newFakeEnsembl(t)
	cache, _ := ensembl.NewCache(t.TempDir())
	recs := []Record{{Chrom: "chr1", Pos: 20000100, Ref: "CTT", Alt: "C"}}
	if err := annotateRecords(context.Background(), client, cache, recs); err != nil {
		t.Fatalf("annotate: %v", err)
	}
	want := []string{"1 20000100 . CTT C . . ."}
	if strings.Join(fe.vepInputs, "|") != strings.Join(want, "|") {
		t.Errorf("VEP inputs = %q, want %q", fe.vepInputs, want)
	}
	if !recs[0].HasAnnot || recs[0].Gene != "FAKE" || recs[0].AF != 0.125 {
		t.Errorf("annotation = %+v, want gene FAKE and AF 0.125", recs[0])
	}
}

// An indel's scan window includes its anchor base. A SNV called on that base is
// a neighbour, so an indel absent from the VCF must still report no-call rather
// than letting the neighbour stand in for it.
func TestResolveRSIDAbsentIndelBesideSNV(t *testing.T) {
	_, client := newFakeEnsembl(t)
	reader := openFixture(t)
	cache, _ := ensembl.NewCache(t.TempDir())
	recs, err := resolveRSID(context.Background(), client, cache, reader, "rs900003")
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 1 || !recs[0].NoCall || recs[0].Pos != 11796322 {
		t.Errorf("got %+v, want a single no-call row at 11796322", recs)
	}
}
