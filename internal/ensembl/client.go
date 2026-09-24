package ensembl

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"strconv"
	"strings"
	"time"
)

// DefaultBaseURL is Ensembl's public REST endpoint.
const DefaultBaseURL = "https://rest.ensembl.org"

// DefaultRateLimit is requests per second. Ensembl publishes a 15 rps ceiling;
// 10 rps stays comfortably under it for shared infrastructure courtesy.
const DefaultRateLimit = 10

// DefaultAttempts is how many times a request is tried before its last error is
// returned. Ensembl's load balancer intermittently fronts unhealthy backends, so
// a 5xx or timeout on one attempt says little about the next.
const DefaultAttempts = 3

// DefaultRetryBase is the delay before the first retry; each later retry
// doubles it.
const DefaultRetryBase = time.Second

// maxRetryAfter caps how long a 429's Retry-After header can make us wait.
const maxRetryAfter = time.Minute

// Client is a rate-limited Ensembl REST client.
type Client struct {
	BaseURL    string
	HTTP       *http.Client
	UserAgent  string
	Attempts   int
	RetryBase  time.Duration
	limiter    <-chan time.Time
	limiterDur time.Duration
}

// NewClient returns a client with the default base URL, a 15s timeout per
// attempt, 10 rps rate limiting, and up to 3 attempts per request. version is
// reported in the User-Agent. The caller may override fields directly after
// construction.
func NewClient(version string) *Client {
	dur := time.Second / time.Duration(DefaultRateLimit)
	return &Client{
		BaseURL:    DefaultBaseURL,
		HTTP:       &http.Client{Timeout: 15 * time.Second},
		UserAgent:  "vcfq/" + version + " (+https://github.com/liminalpurple/vcfq)",
		Attempts:   DefaultAttempts,
		RetryBase:  DefaultRetryBase,
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

// do performs a rate-limited request, decoding the JSON response into dest.
// Transport errors (including timeouts), 5xx and 429 responses are retried with
// exponential backoff, honouring Retry-After on 429; other 4xx responses fail
// immediately. After the last attempt the last error is returned unchanged.
// Returns ErrNotFound for HTTP 404s so callers can distinguish "no such symbol"
// from "network down".
func (c *Client) do(ctx context.Context, method, path string, body []byte, dest any) error {
	delay := c.RetryBase
	for attempt := 1; ; attempt++ {
		retryAfter, err := c.attempt(ctx, method, path, body, dest)
		if err == nil || retryAfter < 0 || attempt >= c.Attempts || ctx.Err() != nil {
			return err
		}
		wait := max(delay, retryAfter)
		delay *= 2
		t := time.NewTimer(wait)
		select {
		case <-t.C:
		case <-ctx.Done():
			t.Stop()
			return err
		}
	}
}

// attempt performs one rate-limited request. On failure, retryAfter reports
// whether it is worth retrying: negative means never, otherwise it is the
// minimum wait the server asked for (zero if it didn't say).
func (c *Client) attempt(ctx context.Context, method, path string, body []byte, dest any) (retryAfter time.Duration, err error) {
	if err := c.wait(ctx); err != nil {
		return -1, err
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
		return -1, err
	}
	req.Header.Set("Accept", "application/json")
	if body != nil {
		req.Header.Set("Content-Type", "application/json")
	}
	req.Header.Set("User-Agent", c.UserAgent)
	resp, err := c.HTTP.Do(req)
	if err != nil {
		return 0, fmt.Errorf("ensembl %s %s: %w", method, path, err)
	}
	defer func() {
		if cerr := resp.Body.Close(); cerr != nil {
			err = errors.Join(err, fmt.Errorf("ensembl %s %s: close body: %w", method, path, cerr))
		}
	}()
	respBody, err := io.ReadAll(resp.Body)
	if err != nil {
		return 0, fmt.Errorf("ensembl %s %s: read body: %w", method, path, err)
	}
	if resp.StatusCode == http.StatusNotFound {
		return -1, ErrNotFound
	}
	if resp.StatusCode >= 400 {
		err := fmt.Errorf("ensembl %s %s: status %d: %s", method, path, resp.StatusCode, truncate(string(respBody), 200))
		switch {
		case resp.StatusCode == http.StatusTooManyRequests:
			return parseRetryAfter(resp.Header.Get("Retry-After")), err
		case resp.StatusCode >= 500:
			return 0, err
		default:
			return -1, err
		}
	}
	if err := json.Unmarshal(respBody, dest); err != nil {
		return -1, fmt.Errorf("ensembl %s %s: decode: %w", method, path, err)
	}
	return 0, nil
}

// parseRetryAfter reads a Retry-After header in seconds. Ensembl sends
// fractional values, so it is parsed as a float. Missing or unparseable values
// give zero, leaving the backoff delay in charge.
func parseRetryAfter(v string) time.Duration {
	secs, err := strconv.ParseFloat(strings.TrimSpace(v), 64)
	if err != nil || secs <= 0 {
		return 0
	}
	return min(time.Duration(secs*float64(time.Second)), maxRetryAfter)
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
