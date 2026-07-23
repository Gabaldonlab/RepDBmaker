# Release custom proteomes as a versioned bundle (the former standalone
# `new_genomes` pipeline, folded in here). See docs/split-pipeline.md.
#
# Two decoupled halves:
#   CURATION (author only) -- run `snakemake package_custom` to build the bundle
#     from locally-curated proteomes (data/custom_proteomes/CUS<id>.fa),
#     then upload it (eventually Zenodo) and pin it below.
#   CONSUMPTION (everyone) -- when dbs.build.repdb.custom_bundle points at a
#     bundle (a local tar.gz for now, a URL later), the sequences come from it,
#     so no local files are needed to reproduce a release. Unset -> the sequences
#     are taken from the files.custom_proteomes folder (CUS<id>.faa.gz | .fa).
#
# The bundle is a flat, self-contained tar.gz:
#   proteomes/CUS<id>.faa.gz   custom_metadata.tsv

import csv
import glob


def repdb_custom_bundle():
    repdb = (config.get("dbs", {}).get("build", {}) or {}).get("repdb") or {}
    return repdb.get("custom_bundle") if isinstance(repdb, dict) else None


CUSTOM_BUNDLE = repdb_custom_bundle()
# folder holding the local custom proteomes (CUS<id>.faa.gz | .fa | ...)
CUSTOM_DIR = config.get("files", {}).get("custom_proteomes", "resources/custom_proteomes")


def _custom_codes():
    """CUS ids from the custom table (drives curation/packaging)."""
    try:
        with open(config["files"]["new_genomes"]) as f:
            return [r["ID"].strip() for r in csv.DictReader(f, delimiter="\t")
                    if r.get("ID") and r["ID"].strip()]
    except (OSError, KeyError):
        return []


CUSTOM_CODES = _custom_codes()


def custom_proteome_path(code):
    """The local proteome file for a CUS id (any extension), from CUSTOM_DIR."""
    hits = sorted(glob.glob(f"{CUSTOM_DIR}/{code}.*"))
    return hits[0] if hits else f"{CUSTOM_DIR}/{code}.faa.gz"


rule unpack_custom_bundle:
    """Extract the release custom bundle (proteomes/ + custom_busco.tsv). For now
    `custom_bundle` is a local tar.gz; a URL+sha256 handler can be added later."""
    input:
        CUSTOM_BUNDLE if CUSTOM_BUNDLE else [],
    output:
        proteomes=directory("results/custom_bundle/proteomes")
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
odir=$(dirname {output.proteomes})
mkdir -p "$odir"
tar -xzf {input} -C "$odir"
"""


rule package_custom_bundle:
    """Assemble the versioned release-custom bundle (proteomes + metadata) into
    one tar.gz. Upload it (eventually Zenodo) and set dbs.build.repdb.custom_bundle
    to its path/URL."""
    input:
        table=config["files"]["new_genomes"],
        proteomes=lambda wc: [custom_proteome_path(c) for c in CUSTOM_CODES],
    output:
        "results/custom/repdb_custom_bundle.tar.gz",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
stage=$(mktemp -d)
mkdir -p "$stage/proteomes"
for f in {input.proteomes}; do
    code=$(basename "$f"); code=${{code%%.*}}
    case "$f" in
        *.gz) cp "$f" "$stage/proteomes/$code.faa.gz" ;;
        *)    gzip -c "$f" > "$stage/proteomes/$code.faa.gz" ;;
    esac
done
cp {input.table} "$stage/custom_metadata.tsv"
mkdir -p $(dirname {output})
tar -czf {output} -C "$stage" .
rm -rf "$stage"
"""
