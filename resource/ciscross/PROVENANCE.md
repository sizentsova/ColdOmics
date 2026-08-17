# Vendored CisCross-local engine and indices

These files let ColdOmics run the cistromics step **offline**, without the
CisCross web service and without a checkout of the CisCross-local repository.

## Source

| | |
|---|---|
| Repository | https://github.com/VictoriaVMironova/CisCross (CisCross-local) |
| Commit | `641dd554f4ab05ef6fb3fbbbc9de36e8c218a506` (2026-07-29) |
| Licence | MIT — see the header of `enrichment.R` |
| Method | Lavrekha et al., *Front. Plant Sci.* 2022, [10.3389/fpls.2022.942710](https://doi.org/10.3389/fpls.2022.942710) |

## Contents

| Path | Origin |
|---|---|
| `enrichment.R` | verbatim copy of `R/enrichment.R` |
| `index/ciscross_curated608_up500.rds` | `Rscript R/build_index.R 500 curated608` |
| `index/ciscross_curated608_up1000.rds` | `Rscript R/build_index.R 1000 curated608` |
| `index/ciscross_curated608_up1500.rds` | `Rscript R/build_index.R 1500 curated608` |
| `index/ciscross_curated608_up2000.rds` | `Rscript R/build_index.R 2000 curated608` |
| `index/ciscross_cistrome_up1500.rds` | `Rscript R/build_index.R 1500 cistrome` |

Nothing else from CisCross-local is needed at run time: `enrichment.R` imports
only `data.table`, and the indices are self-contained. The DAP-seq peak BEDs
(~400 MB), the TAIR10 GFF3 (~100 MB), Quarto and `renv` are **build-time only**.

## What an index is

A list with `$universe` (all genes with a defined promoter), `$incidence`
(peak set → genes whose promoter overlaps a peak) and `$meta` (peak set → TF
name, AGI, family, DNA source, peak count).

The index does not depend on any input gene list, so one file serves any set of
DEGs, any gene count, and any custom background. Freezing it here also freezes
the DAP-seq collection, the promoter definition and the TAIR10 build used for
the published analysis.

## Collections

- **`curated608`** (default) — CisCross-MACS2. DAP-seq reads re-called with
  MACS2; per TF × DNA source profile the IDR-merged set is used when it has
  > 2000 peaks or ≥ 2× the best replica, otherwise the highest-FRIP replica;
  sets with < 200 peaks dropped. 628 sets / 413 TFs.
- **`cistrome`** — Plant Cistrome-GEM (O'Malley et al. 2016): the original
  GEM-called peaks from Salk. 568 sets / 387 TFs.

Both describe the same underlying DAP-seq experiments through different
peak-calling pipelines. Run them **separately** and intersect the results as a
robustness check — do not pool them into one test, which would double-count
each TF while inflating the multiple-testing correction.

Promoters are defined upstream of the TSS in a strand-aware manner. The default
1500-bp setting aligns most closely with the web service; curated608 indices for
500, 1000 and 2000 bp support promoter-length sensitivity analysis.

## Rebuilding

```bash
git clone https://github.com/VictoriaVMironova/CisCross && cd CisCross
# follow README "Data setup" to fetch the peak sets and TAIR10 GFF3
Rscript R/build_index.R 1500 curated608
Rscript R/build_index.R 1500 cistrome
```

The vendored sensitivity-analysis indices were built with:

```bash
Rscript R/build_index.R 500 curated608
Rscript R/build_index.R 1000 curated608
Rscript R/build_index.R 2000 curated608
```

ColdOmics derives report names from `config.yaml`; changing
`cistromics.upstream` or `cistromics.fdr` no longer requires R script edits.

## SHA-256 checksums

| Index | SHA-256 |
|---|---|
| `ciscross_curated608_up500.rds` | `d0a0af510a53abdcca75d8f806bf892e6e90242c1dabff5f1812f29bad7633f4` |
| `ciscross_curated608_up1000.rds` | `865c3ae8e84230c88fdbe6ca902282f5faa73bda61baa420606953e5752fa2eb` |
| `ciscross_curated608_up1500.rds` | `71a5f91640b30e7cd87ae257257034d5df2f3089d41d30a55c7c815234942e85` |
| `ciscross_curated608_up2000.rds` | `dbc2d6a9e70db019f6453a9fe57367bc6da121561e81f23b3188a566bf776a59` |
