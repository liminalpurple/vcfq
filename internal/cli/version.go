package cli

// Version is the build version, set via the linker at release time:
//
//	go build -ldflags "-X github.com/liminalpurple/vcfq/internal/cli.Version=v0.1.0"
//
// Defaults to "dev" for local builds.
var Version = "dev"
