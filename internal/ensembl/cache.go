// Package ensembl is a small client for the public Ensembl REST API, with a
// file-backed cache and polite rate limiting. It exposes only what vcfq needs:
// gene symbol lookup, rsID lookup, and VEP batch annotation. All coordinates are
// GRCh38 and 1-based inclusive.
package ensembl

import (
	"encoding/json"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
)

// Cache stores Ensembl responses on disk as JSON files keyed by query type and
// identifier. There is no TTL; entries are stable per Ensembl release in
// practice. Wipe with os.RemoveAll(c.Dir()) or `vcfq cache clean`.
type Cache struct {
	dir string
}

// NewCache returns a cache rooted at $XDG_CACHE_HOME/vcfq (or its OS equivalent).
// If dir is non-empty it overrides the default.
func NewCache(dir string) (*Cache, error) {
	if dir == "" {
		base, err := os.UserCacheDir()
		if err != nil {
			return nil, fmt.Errorf("resolve user cache dir: %w", err)
		}
		dir = filepath.Join(base, "vcfq")
	}
	return &Cache{dir: dir}, nil
}

// Dir returns the cache's root directory.
func (c *Cache) Dir() string { return c.dir }

// Get unmarshals the cached entry for (kind, key) into dest. Returns (false, nil)
// on a clean miss. Errors only on corrupted entries or filesystem failures.
func (c *Cache) Get(kind, key string, dest any) (bool, error) {
	path := c.path(kind, key)
	b, err := os.ReadFile(path)
	if errors.Is(err, fs.ErrNotExist) {
		return false, nil
	}
	if err != nil {
		return false, fmt.Errorf("read cache %s: %w", path, err)
	}
	if err := json.Unmarshal(b, dest); err != nil {
		return false, fmt.Errorf("decode cache %s: %w", path, err)
	}
	return true, nil
}

// Set writes val under (kind, key), creating directories as needed. Writes are
// atomic via tmp+rename so a partial entry never appears on a concurrent reader.
func (c *Cache) Set(kind, key string, val any) error {
	path := c.path(kind, key)
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	b, err := json.MarshalIndent(val, "", "  ")
	if err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(path), ".vcfq-*.tmp")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	if _, err := tmp.Write(b); err != nil {
		tmp.Close()
		os.Remove(tmpName)
		return err
	}
	if err := tmp.Close(); err != nil {
		os.Remove(tmpName)
		return err
	}
	return os.Rename(tmpName, path)
}

// Clean removes every cached entry.
func (c *Cache) Clean() error {
	if c.dir == "" {
		return errors.New("cache dir is empty")
	}
	return os.RemoveAll(c.dir)
}

// path returns the on-disk path for (kind, key). Keys are sanitised by
// replacing path separators with underscores so a malicious or weird identifier
// can't escape the cache root.
func (c *Cache) path(kind, key string) string {
	safeKey := strings.NewReplacer("/", "_", "\\", "_", "..", "_").Replace(key)
	return filepath.Join(c.dir, kind, safeKey+".json")
}
