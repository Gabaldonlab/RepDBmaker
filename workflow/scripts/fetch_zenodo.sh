#!/usr/bin/env bash
# Download one Zenodo-hosted file and verify it against a pinned sha256 before
# it's considered done, so a truncated/corrupt download never lands at the
# final path Snakemake tracks as complete. Resumable (wget -c) since these are
# routinely tens of GB.
#
# Usage: fetch_zenodo.sh <url> <sha256> <output_path>
set -euo pipefail

url="$1"
sum="$2"
out="$3"

tmp="${out}.part"
wget -c --tries=8 --retry-connrefused -O "$tmp" "$url"
echo "${sum}  ${tmp}" | sha256sum -c -
mv "$tmp" "$out"
