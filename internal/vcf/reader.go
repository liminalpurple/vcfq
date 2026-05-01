// Package vcf reads a tabix-indexed bgzipped VCF for region-bounded scans. It's
// intentionally minimal: just what vcfq needs to extract Records given
// (chrom, start, end) inputs from Ensembl coordinate resolution.
package vcf

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bgzf"
	bgzfindex "github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/tabix"
)

// Reader holds an open bgzip file and the tabix index alongside it. Header
// holds every header line (## and the #CHROM line) for VCF passthrough.
type Reader struct {
	path   string
	f      *os.File
	bgz    *bgzf.Reader
	idx    *tabix.Index
	Header []string
	names  map[string]string // canonical chrom name as stored in the index, keyed by both "1" and "chr1"
}

// Open reads the VCF header and loads the tabix index from path+".tbi".
func Open(path string) (*Reader, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open vcf: %w", err)
	}
	bgz, err := bgzf.NewReader(f, 1)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("bgzf: %w", err)
	}

	idxF, err := os.Open(path + ".tbi")
	if err != nil {
		bgz.Close()
		f.Close()
		return nil, fmt.Errorf("open tabix index (%s.tbi): %w", path, err)
	}
	defer idxF.Close()
	// .tbi is bgzipped; gzip.Reader handles bgzf streams transparently.
	gz, err := gzip.NewReader(idxF)
	if err != nil {
		bgz.Close()
		f.Close()
		return nil, fmt.Errorf("decompress tabix index: %w", err)
	}
	idx, err := tabix.ReadFrom(gz)
	gz.Close()
	if err != nil {
		bgz.Close()
		f.Close()
		return nil, fmt.Errorf("read tabix index: %w", err)
	}

	r := &Reader{path: path, f: f, bgz: bgz, idx: idx}
	if err := r.readHeader(); err != nil {
		r.Close()
		return nil, err
	}
	r.buildNameMap()
	return r, nil
}

// Close releases the underlying file handles.
func (r *Reader) Close() error {
	if r.bgz != nil {
		_ = r.bgz.Close()
	}
	if r.f != nil {
		return r.f.Close()
	}
	return nil
}

// readHeader scans from the start of the bgzip stream and collects every line
// that starts with '#'. The bgzf reader is left at an arbitrary position
// afterwards; ScanRegion seeks to a tabix chunk before reading data.
func (r *Reader) readHeader() error {
	sc := bufio.NewScanner(r.bgz)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<22)
	for sc.Scan() {
		line := sc.Text()
		if !strings.HasPrefix(line, "#") {
			break
		}
		r.Header = append(r.Header, line)
	}
	if err := sc.Err(); err != nil && !errors.Is(err, io.EOF) {
		return fmt.Errorf("read header: %w", err)
	}
	if len(r.Header) == 0 {
		return errors.New("no VCF header found")
	}
	return nil
}

// buildNameMap records each chrom name from the index under both its bare and
// chr-prefixed forms so callers can pass either.
func (r *Reader) buildNameMap() {
	r.names = make(map[string]string, r.idx.NumRefs()*2)
	for _, n := range r.idx.Names() {
		r.names[n] = n
		if strings.HasPrefix(n, "chr") {
			r.names[strings.TrimPrefix(n, "chr")] = n
		} else {
			r.names["chr"+n] = n
		}
	}
}

// resolveChrom maps a user-supplied chrom name to the canonical form stored in
// the tabix index. Returns "" if the chrom isn't present.
func (r *Reader) resolveChrom(chrom string) string {
	if r.names == nil {
		return chrom
	}
	if c, ok := r.names[chrom]; ok {
		return c
	}
	return ""
}

// DataLine is a single VCF data record returned from ScanRegion. Fields are
// the raw VCF columns; conversion to internal/cli.Record happens in dispatch.
type DataLine struct {
	Chrom  string
	Pos    int
	ID     string
	Ref    string
	Alts   []string // multi-allelic split into separate alts
	Qual   string
	Filter string
	Info   string
	Format string
	Sample string
}

// ScanRegion returns every VCF data line whose POS falls within the closed
// 1-based interval [start, end] on chrom. Tabix narrows to a bin; we filter
// the precise interval ourselves. Multi-allelic ALT fields stay packed into
// DataLine.Alts and are split downstream.
func (r *Reader) ScanRegion(chrom string, start, end int) ([]DataLine, error) {
	canonical := r.resolveChrom(chrom)
	if canonical == "" {
		return nil, fmt.Errorf("chrom %q not present in tabix index for %s", chrom, r.path)
	}
	// Tabix's Chunks API uses 0-based half-open coords.
	chunks, err := r.idx.Chunks(canonical, start-1, end)
	if err != nil {
		if errors.Is(err, errNoRef()) {
			return nil, nil
		}
		return nil, fmt.Errorf("tabix chunks: %w", err)
	}

	if len(chunks) == 0 {
		return nil, nil
	}
	cr, err := bgzfindex.NewChunkReader(r.bgz, chunks)
	if err != nil {
		return nil, fmt.Errorf("bgzf chunk reader: %w", err)
	}
	defer cr.Close()
	sc := bufio.NewScanner(cr)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<22)

	var out []DataLine
	for sc.Scan() {
		line := sc.Text()
		if line == "" || line[0] == '#' {
			continue
		}
		d, ok := parseDataLine(line)
		if !ok {
			continue
		}
		if d.Chrom != canonical {
			continue
		}
		// VCFs are position-sorted within a chrom, so once past `end` we're done.
		if d.Pos > end {
			break
		}
		if d.Pos < start {
			continue
		}
		out = append(out, d)
	}
	if err := sc.Err(); err != nil && !errors.Is(err, io.EOF) {
		return nil, fmt.Errorf("scan chunk: %w", err)
	}
	return out, nil
}

// errNoRef returns the index.ErrNoReference sentinel as a generic error so we
// can compare without importing biogo's internal index package.
func errNoRef() error {
	return errors.New("no reference")
}

// parseDataLine splits a VCF data line into columns. Returns ok=false on lines
// that don't have at least the 8 mandatory VCF columns.
func parseDataLine(line string) (DataLine, bool) {
	parts := strings.Split(line, "\t")
	if len(parts) < 8 {
		return DataLine{}, false
	}
	pos, err := strconv.Atoi(parts[1])
	if err != nil {
		return DataLine{}, false
	}
	d := DataLine{
		Chrom:  parts[0],
		Pos:    pos,
		ID:     parts[2],
		Ref:    parts[3],
		Alts:   strings.Split(parts[4], ","),
		Qual:   parts[5],
		Filter: parts[6],
		Info:   parts[7],
	}
	if len(parts) >= 10 {
		d.Format = parts[8]
		d.Sample = parts[9]
	}
	return d, true
}
