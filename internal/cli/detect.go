package cli

import "regexp"

// QueryKind classifies a single token from CLI args or stdin.
type QueryKind int

const (
	KindSymbol QueryKind = iota
	KindRSID
	KindRegion
)

func (k QueryKind) String() string {
	switch k {
	case KindRSID:
		return "rsid"
	case KindRegion:
		return "region"
	default:
		return "symbol"
	}
}

var (
	rsIDRegex   = regexp.MustCompile(`^[rR][sS]\d+$`)
	regionRegex = regexp.MustCompile(`^(?:chr)?[\dXYMxym]+:\d+(?:-\d+)?$`)
)

// Detect classifies a query token. rsIDs and regions follow strict patterns;
// anything else is treated as a gene symbol and passed to Ensembl for resolution.
func Detect(q string) QueryKind {
	switch {
	case rsIDRegex.MatchString(q):
		return KindRSID
	case regionRegex.MatchString(q):
		return KindRegion
	default:
		return KindSymbol
	}
}
