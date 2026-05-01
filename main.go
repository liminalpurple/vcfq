// Command vcfq queries a tabix-indexed VCF with Ensembl coordinate resolution.
package main

import (
	"os"

	"github.com/liminalpurple/vcfq/internal/cli"
)

func main() {
	os.Exit(cli.Run(os.Args[1:], os.Stdin, os.Stdout, os.Stderr))
}
