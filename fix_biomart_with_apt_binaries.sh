#!/usr/bin/env bash
set -euo pipefail

# Use Ubuntu binary R packages for the fragile biomaRt/dbplyr/BiocFileCache stack.
# This avoids source-build failures for dbplyr and downstream Bioconductor packages.

sudo apt update
sudo apt install -y \
  r-cran-dbplyr \
  r-cran-dplyr \
  r-cran-rsqlite \
  r-cran-dbi \
  r-cran-filelock \
  r-cran-curl \
  r-cran-openssl \
  r-cran-httr \
  r-cran-xml2 \
  r-bioc-biocfilecache \
  r-bioc-annotationdbi \
  r-bioc-biomart

Rscript - <<'RSCRIPT'
cat("Testing core biomaRt stack...\n")
needed <- c("dbplyr", "BiocFileCache", "biomaRt")
failed <- character()
for (pkg in needed) {
  ok <- suppressWarnings(suppressPackageStartupMessages(require(pkg, character.only = TRUE)))
  if (!ok) failed <- c(failed, pkg)
}
if (length(failed) > 0) {
  stop("Still failed to load: ", paste(failed, collapse = ", "))
}
cat("biomaRt stack OK\n")
RSCRIPT

echo "Done. Now rerun: bash ./run_rnaseq_then_network_test.sh"
