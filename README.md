# vcfq

CLI toolkit for genome queries — gene-symbol, rsID, and region lookups against a tabix-indexed
VCF, with Ensembl coordinate resolution and optional VEP annotation. Single static binary, zero
runtime dependencies.

## Install

```sh
go install github.com/liminalpurple/vcfq@latest
```

Requires Go 1.22+. The binary needs network access for Ensembl REST calls; all results are cached
locally so repeat queries are offline-fast.

## Quick start

```sh
export VCFQ_VCF=/path/to/sample.vcf.gz   # tabix-indexed (.tbi alongside)

# Look up by gene symbol — Ensembl resolves the GRCh38 region, vcfq scans it
vcfq HNF1A

# By rsID — Ensembl resolves the location, vcfq grabs the row
vcfq rs1801133

# By genomic region (GRCh38, 1-based inclusive)
vcfq chr1:11796000-11796500

# Multiple queries in one go (CLI args, stdin, or both)
vcfq HNF1A MTHFR rs1801133
echo 'HNF1A rs1801133 chr12:120978000-120998000' | vcfq
```

## Annotation

Pass `-a` (or `--annotate`) to enrich each variant with VEP consequence, gene symbol, amino acid
change, and gnomAD/1000G allele frequency:

```sh
vcfq -a rs1801133
# chrom  pos       id          ref  alt  gt   query       consequence       gene   aa_change  af
# chr1   11796321  rs1801133   G    A    1/1  rs1801133   missense_variant  MTHFR  A/V        0.2752
```

VEP calls are batched (≤200 variants per request) and cached, so re-running an annotated query is
instant.

## Output formats

`-f tsv` (default), `-f table` (Markdown), `-f json` (newline-delimited), `-f vcf` (VCF
passthrough).

The `vcf` format is the headline: it produces a strictly richer VCF than the input. The source
header is passed through, `##INFO=<ID=VCFQ_*>` declarations are added when `-a` is set, and
annotations are injected into each line's INFO field. The output round-trips cleanly through
`bcftools view` and `bcftools query`:

```sh
vcfq -a -f vcf HNF1A > annotated.vcf
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/VCFQ_GENE\t%INFO/VCFQ_CONSEQUENCE\n' annotated.vcf
```

This makes vcfq useful as a preprocessing step for VCFs that ship without functional annotations
(e.g. raw consumer-genomics output).

## Cache

Ensembl symbol, rsID, and VEP results are cached as JSON files under `os.UserCacheDir()/vcfq/`
(typically `$XDG_CACHE_HOME/vcfq` or `~/.cache/vcfq`). There is no TTL — entries are stable per
Ensembl release in practice. Wipe with:

```sh
vcfq cache clean
```

## Scope

vcfq deliberately handles VCF-level work only. CRAM read-level queries (depth ratios, CIGAR
analysis, paralog disambiguation) belong with `samtools` / `bcftools`, which already handle them
ergonomically. vcfq complements those tools rather than reimplementing them.

All coordinates are GRCh38, 1-based inclusive.

## Licence

Apache-2.0.
