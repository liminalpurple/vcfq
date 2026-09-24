package ensembl

import (
	"context"
	"errors"
	"net/http"
	"net/http/httptest"
	"strings"
	"sync/atomic"
	"testing"
	"time"
)

// scriptedServer answers each request with the next status in statuses,
// repeating the last one once the script runs out. A 200 carries {"ok":true}.
func scriptedServer(t *testing.T, statuses []int, header http.Header) (*Client, *atomic.Int32) {
	t.Helper()
	var hits atomic.Int32
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, _ *http.Request) {
		n := int(hits.Add(1))
		status := statuses[min(n, len(statuses))-1]
		for k, v := range header {
			w.Header()[k] = v
		}
		w.WriteHeader(status)
		if status == http.StatusOK {
			_, _ = w.Write([]byte(`{"ok":true}`))
		} else {
			_, _ = w.Write([]byte("<html>error</html>"))
		}
	}))
	t.Cleanup(srv.Close)
	return &Client{
		BaseURL:   srv.URL,
		HTTP:      srv.Client(),
		UserAgent: "vcfq-test",
		Attempts:  3,
		RetryBase: time.Millisecond,
	}, &hits
}

func TestDoRetries(t *testing.T) {
	cases := []struct {
		name     string
		statuses []int
		wantHits int32
		wantErr  string
	}{
		{"transient 5xx then success", []int{500, 502, 200}, 3, ""},
		{"5xx exhausts attempts", []int{500}, 3, "status 500"},
		{"4xx is not retried", []int{400}, 1, "status 400"},
		{"404 is not retried", []int{404}, 1, "not found"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			client, hits := scriptedServer(t, c.statuses, nil)
			var out struct{ OK bool }
			err := client.getJSON(context.Background(), "/x", &out)
			if got := hits.Load(); got != c.wantHits {
				t.Errorf("hits = %d, want %d", got, c.wantHits)
			}
			if c.wantErr == "" {
				if err != nil || !out.OK {
					t.Fatalf("err = %v, ok = %v; want success", err, out.OK)
				}
				return
			}
			if err == nil || !strings.Contains(err.Error(), c.wantErr) {
				t.Fatalf("err = %v, want containing %q", err, c.wantErr)
			}
		})
	}
}

func TestDoNotFoundIsSentinel(t *testing.T) {
	client, _ := scriptedServer(t, []int{404}, nil)
	var out struct{}
	if err := client.getJSON(context.Background(), "/x", &out); !errors.Is(err, ErrNotFound) {
		t.Fatalf("err = %v, want ErrNotFound", err)
	}
}

func TestDoHonoursRetryAfter(t *testing.T) {
	client, hits := scriptedServer(t, []int{429, 200}, http.Header{"Retry-After": {"0.05"}})
	start := time.Now()
	var out struct{ OK bool }
	if err := client.getJSON(context.Background(), "/x", &out); err != nil {
		t.Fatalf("getJSON: %v", err)
	}
	if elapsed := time.Since(start); elapsed < 50*time.Millisecond {
		t.Errorf("elapsed = %v, want at least the 50ms Retry-After", elapsed)
	}
	if got := hits.Load(); got != 2 {
		t.Errorf("hits = %d, want 2", got)
	}
}

func TestDoRetriesTimeouts(t *testing.T) {
	var hits atomic.Int32
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if hits.Add(1) == 1 {
			select {
			case <-time.After(time.Second):
			case <-r.Context().Done():
			}
			return
		}
		_, _ = w.Write([]byte(`{"ok":true}`))
	}))
	t.Cleanup(srv.Close)
	httpClient := srv.Client()
	httpClient.Timeout = 50 * time.Millisecond
	client := &Client{BaseURL: srv.URL, HTTP: httpClient, Attempts: 3, RetryBase: time.Millisecond}

	var out struct{ OK bool }
	if err := client.getJSON(context.Background(), "/x", &out); err != nil || !out.OK {
		t.Fatalf("err = %v, ok = %v; want success after a timed-out first attempt", err, out.OK)
	}
	if got := hits.Load(); got != 2 {
		t.Errorf("hits = %d, want 2", got)
	}
}

func TestDoStopsOnCancelledContext(t *testing.T) {
	client, hits := scriptedServer(t, []int{500}, nil)
	client.RetryBase = time.Hour
	ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
	defer cancel()
	var out struct{}
	if err := client.getJSON(ctx, "/x", &out); err == nil || !strings.Contains(err.Error(), "status 500") {
		t.Fatalf("err = %v, want the 500 that preceded cancellation", err)
	}
	if got := hits.Load(); got != 1 {
		t.Errorf("hits = %d, want 1", got)
	}
}

func TestParseRetryAfter(t *testing.T) {
	cases := map[string]time.Duration{
		"":      0,
		"junk":  0,
		"-1":    0,
		"2":     2 * time.Second,
		"0.25":  250 * time.Millisecond,
		" 1 ":   time.Second,
		"86400": maxRetryAfter,
	}
	for in, want := range cases {
		if got := parseRetryAfter(in); got != want {
			t.Errorf("parseRetryAfter(%q) = %v, want %v", in, got, want)
		}
	}
}
