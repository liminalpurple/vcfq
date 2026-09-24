package cli

import (
	"bufio"
	"io"
	"os"
)

// collectQueries returns the queries to run. Arguments are used as given, and
// stdin is read only when asked for: when there are no arguments at all (and
// stdin isn't a terminal, so an interactive shell never blocks), or where an
// argument is "-", whose place the stdin tokens take. A script or cron job that
// passes queries as arguments therefore never waits on an inherited stdin.
func collectQueries(args []string, stdin io.Reader) ([]string, error) {
	if len(args) == 0 {
		if isTTY(stdin) {
			return nil, nil
		}
		return readTokens(stdin)
	}
	var (
		out  []string
		read bool
	)
	for _, a := range args {
		if a != "-" {
			out = append(out, a)
			continue
		}
		if read {
			continue // stdin can only be consumed once
		}
		tokens, err := readTokens(stdin)
		if err != nil {
			return nil, err
		}
		out = append(out, tokens...)
		read = true
	}
	return out, nil
}

// readTokens returns every whitespace-separated token from r.
func readTokens(r io.Reader) ([]string, error) {
	var tokens []string
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 1<<20)
	scanner.Split(bufio.ScanWords)
	for scanner.Scan() {
		tokens = append(tokens, scanner.Text())
	}
	return tokens, scanner.Err()
}

// isTTY reports whether r is a terminal (interactive). Anything that isn't an
// *os.File, or whose mode lacks ModeCharDevice, counts as not-a-TTY.
func isTTY(r io.Reader) bool {
	f, ok := r.(*os.File)
	if !ok {
		return false
	}
	st, err := f.Stat()
	if err != nil {
		return false
	}
	return st.Mode()&os.ModeCharDevice != 0
}
