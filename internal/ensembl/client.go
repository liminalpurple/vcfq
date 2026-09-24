package ensembl

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"strings"
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
// 10 rps rate limiting. version is reported in the User-Agent. The caller may
// override fields directly after construction.
func NewClient(version string) *Client {
	dur := time.Second / time.Duration(DefaultRateLimit)
	return &Client{
		BaseURL:    DefaultBaseURL,
		HTTP:       &http.Client{Timeout: 30 * time.Second},
		UserAgent:  "vcfq/" + version + " (+https://github.com/liminalpurple/vcfq)",
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

// getJSON does a GET, decoding the JSON body into dest.
func (c *Client) getJSON(ctx context.Context, path string, dest any) error {
	return c.do(ctx, http.MethodGet, path, nil, dest)
}

// postJSON sends a JSON body and decodes the JSON response.
func (c *Client) postJSON(ctx context.Context, path string, body, dest any) error {
	bb, err := json.Marshal(body)
	if err != nil {
		return err
	}
	return c.do(ctx, http.MethodPost, path, bb, dest)
}

// do performs one rate-limited request, decoding the JSON response into dest.
// Returns ErrNotFound for HTTP 404s so callers can distinguish "no such symbol"
// from "network down".
func (c *Client) do(ctx context.Context, method, path string, body []byte, dest any) (err error) {
	if err := c.wait(ctx); err != nil {
		return err
	}
	url := c.BaseURL + path
	if !strings.Contains(path, "?") {
		url += "?content-type=application/json"
	}
	var rb io.Reader
	if body != nil {
		rb = bytes.NewReader(body)
	}
	req, err := http.NewRequestWithContext(ctx, method, url, rb)
	if err != nil {
		return err
	}
	req.Header.Set("Accept", "application/json")
	if body != nil {
		req.Header.Set("Content-Type", "application/json")
	}
	req.Header.Set("User-Agent", c.UserAgent)
	resp, err := c.HTTP.Do(req)
	if err != nil {
		return fmt.Errorf("ensembl %s %s: %w", method, path, err)
	}
	defer func() {
		if cerr := resp.Body.Close(); cerr != nil {
			err = errors.Join(err, fmt.Errorf("ensembl %s %s: close body: %w", method, path, cerr))
		}
	}()
	respBody, err := io.ReadAll(resp.Body)
	if err != nil {
		return fmt.Errorf("ensembl %s %s: read body: %w", method, path, err)
	}
	if resp.StatusCode == http.StatusNotFound {
		return ErrNotFound
	}
	if resp.StatusCode >= 400 {
		return fmt.Errorf("ensembl %s %s: status %d: %s", method, path, resp.StatusCode, truncate(string(respBody), 200))
	}
	if err := json.Unmarshal(respBody, dest); err != nil {
		return fmt.Errorf("ensembl %s %s: decode: %w", method, path, err)
	}
	return nil
}

// ErrNotFound is returned when Ensembl responds with HTTP 404 — typically an
// unknown symbol or rsID.
var ErrNotFound = errors.New("not found")

func truncate(s string, n int) string {
	if len(s) <= n {
		return s
	}
	return s[:n] + "..."
}
