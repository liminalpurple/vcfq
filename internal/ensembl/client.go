package ensembl

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"time"
)

// DefaultBaseURL is Ensembl's public REST endpoint.
const DefaultBaseURL = "https://rest.ensembl.org"

// DefaultRateLimit is requests per second. Ensembl publishes a 15 rps ceiling;
// 10 rps stays comfortably under it for shared infrastructure courtesy.
const DefaultRateLimit = 10

// Client is a rate-limited Ensembl REST client.
type Client struct {
	BaseURL    string
	HTTP       *http.Client
	UserAgent  string
	limiter    <-chan time.Time
	limiterDur time.Duration
}

// NewClient returns a client with the default base URL, a 30s HTTP timeout, and
// 10 rps rate limiting. The caller may override fields directly after
// construction.
func NewClient() *Client {
	dur := time.Second / time.Duration(DefaultRateLimit)
	return &Client{
		BaseURL:    DefaultBaseURL,
		HTTP:       &http.Client{Timeout: 30 * time.Second},
		UserAgent:  "vcfq/0.1 (+https://github.com/liminalpurple/vcfq)",
		limiter:    time.Tick(dur),
		limiterDur: dur,
	}
}

// wait blocks until the next rate-limit token is available, or the context is
// done.
func (c *Client) wait(ctx context.Context) error {
	if c.limiter == nil {
		return nil
	}
	select {
	case <-c.limiter:
		return nil
	case <-ctx.Done():
		return ctx.Err()
	}
}

// getJSON does a GET, decoding the JSON body into dest. Returns a typed error
// for HTTP 404s so callers can distinguish "no such symbol" from "network down".
func (c *Client) getJSON(ctx context.Context, path string, dest any) error {
	if err := c.wait(ctx); err != nil {
		return err
	}
	url := c.BaseURL + path
	if !contains(path, '?') {
		url += "?content-type=application/json"
	}
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, url, nil)
	if err != nil {
		return err
	}
	req.Header.Set("Accept", "application/json")
	req.Header.Set("User-Agent", c.UserAgent)
	resp, err := c.HTTP.Do(req)
	if err != nil {
		return fmt.Errorf("ensembl GET %s: %w", path, err)
	}
	defer resp.Body.Close()
	body, err := io.ReadAll(resp.Body)
	if err != nil {
		return fmt.Errorf("ensembl GET %s: read body: %w", path, err)
	}
	if resp.StatusCode == http.StatusNotFound {
		return ErrNotFound
	}
	if resp.StatusCode >= 400 {
		return fmt.Errorf("ensembl GET %s: status %d: %s", path, resp.StatusCode, truncate(string(body), 200))
	}
	if err := json.Unmarshal(body, dest); err != nil {
		return fmt.Errorf("ensembl GET %s: decode: %w", path, err)
	}
	return nil
}

// postJSON sends a JSON body and decodes the JSON response.
func (c *Client) postJSON(ctx context.Context, path string, body, dest any) error {
	if err := c.wait(ctx); err != nil {
		return err
	}
	bb, err := json.Marshal(body)
	if err != nil {
		return err
	}
	url := c.BaseURL + path
	if !contains(path, '?') {
		url += "?content-type=application/json"
	}
	req, err := http.NewRequestWithContext(ctx, http.MethodPost, url, bytes.NewReader(bb))
	if err != nil {
		return err
	}
	req.Header.Set("Accept", "application/json")
	req.Header.Set("Content-Type", "application/json")
	req.Header.Set("User-Agent", c.UserAgent)
	resp, err := c.HTTP.Do(req)
	if err != nil {
		return fmt.Errorf("ensembl POST %s: %w", path, err)
	}
	defer resp.Body.Close()
	rb, err := io.ReadAll(resp.Body)
	if err != nil {
		return fmt.Errorf("ensembl POST %s: read body: %w", path, err)
	}
	if resp.StatusCode >= 400 {
		return fmt.Errorf("ensembl POST %s: status %d: %s", path, resp.StatusCode, truncate(string(rb), 200))
	}
	if err := json.Unmarshal(rb, dest); err != nil {
		return fmt.Errorf("ensembl POST %s: decode: %w", path, err)
	}
	return nil
}

// ErrNotFound is returned when Ensembl responds with HTTP 404 — typically an
// unknown symbol or rsID.
var ErrNotFound = fmt.Errorf("not found")

func contains(s string, b byte) bool {
	for i := 0; i < len(s); i++ {
		if s[i] == b {
			return true
		}
	}
	return false
}

func truncate(s string, n int) string {
	if len(s) <= n {
		return s
	}
	return s[:n] + "..."
}
