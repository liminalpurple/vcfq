package ensembl

import (
	"path/filepath"
	"testing"
)

func TestCacheRoundTrip(t *testing.T) {
	dir := t.TempDir()
	c, err := NewCache(dir)
	if err != nil {
		t.Fatalf("NewCache: %v", err)
	}

	type payload struct {
		A string
		B int
	}
	in := payload{A: "hello", B: 42}

	// Miss before set
	var out payload
	ok, err := c.Get("symbol", "FOO", &out)
	if err != nil {
		t.Fatalf("Get miss: %v", err)
	}
	if ok {
		t.Fatal("expected cache miss before Set")
	}

	if err := c.Set("symbol", "FOO", in); err != nil {
		t.Fatalf("Set: %v", err)
	}
	ok, err = c.Get("symbol", "FOO", &out)
	if err != nil {
		t.Fatalf("Get hit: %v", err)
	}
	if !ok {
		t.Fatal("expected cache hit after Set")
	}
	if out != in {
		t.Errorf("round-trip mismatch: got %+v, want %+v", out, in)
	}

	// Subkey directory should exist on disk
	if _, err := filepath.Abs(filepath.Join(dir, "symbol", "FOO.json")); err != nil {
		t.Errorf("expected file path resolvable: %v", err)
	}

	// Clean wipes everything
	if err := c.Clean(); err != nil {
		t.Fatalf("Clean: %v", err)
	}
	ok, _ = c.Get("symbol", "FOO", &out)
	if ok {
		t.Fatal("expected miss after Clean")
	}
}

func TestCacheKeySanitisation(t *testing.T) {
	dir := t.TempDir()
	c, _ := NewCache(dir)
	// A key containing path separators must not escape the cache root.
	if err := c.Set("evil", "../../etc/passwd", "boom"); err != nil {
		t.Fatalf("Set: %v", err)
	}
	var s string
	ok, err := c.Get("evil", "../../etc/passwd", &s)
	if err != nil || !ok {
		t.Fatalf("round-trip with sanitised key failed: ok=%v err=%v", ok, err)
	}
	if s != "boom" {
		t.Errorf("expected 'boom', got %q", s)
	}
}
