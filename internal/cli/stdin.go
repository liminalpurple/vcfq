package cli

import (
	"bufio"
	"io"
	"os"
)

// readStdinTokens returns whitespace-separated tokens from stdin if it's a pipe
// or redirected file. Returns nil if stdin is a TTY (interactive shell) so we
// never block waiting for keyboard input.
func readStdinTokens(stdin io.Reader) []string {
	if isTTY(stdin) {
		return nil
	}
	var tokens []string
	scanner := bufio.NewScanner(stdin)
	scanner.Buffer(make([]byte, 0, 64*1024), 1<<20)
	scanner.Split(bufio.ScanWords)
	for scanner.Scan() {
		if t := scanner.Text(); t != "" {
			tokens = append(tokens, t)
		}
	}
	return tokens
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
