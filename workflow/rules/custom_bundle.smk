# Release custom proteomes as a versioned bundle (the former standalone
# `new_genomes` pipeline, folded in here). See docs/split-pipeline.md.
#
# Two decoupled halves:
#   CURATION (author only) -- run `snakemake package_custom` to build the bundle
#     from locally-curated proteomes (data/custom_proteomes/CUS<id>.fa) + BUSCO,
#     then upload it (eventually Zenodo) and pin it below.
#   CONSUMPTION (everyone) -- when dbs.build.repdb.custom_bundle points at a
#     bundle (a local tar.gz for now, a URL later), the release custom sequences
#     and their BUSCO completeness are taken from it, so no local files are
#     needed to reproduce a release. Unset -> the old local-Fasta path is used.
#
# The bundle is a flat, self-contained tar.gz:
#   proteomes/CUS<id>.faa.gz   custom_busco.tsv   custom_metadata.tsv   methods.txt

import csv


def repdb_custom_bundle():
    repdb = (config.get("dbs", {}).get("build", {}) or {}).get("repdb") or {}
    return repdb.get("custom_bundle") if isinstance(repdb, dict) else None


CUSTOM_BUNDLE = repdb_custom_bundle()


def _custom_codes():
    """CUS ids from the custom table (drives curation/packaging)."""
    try:
        with open(config["files"]["new_genomes"]) as f:
            return [r["ID"].strip() for r in csv.DictReader(f, delimiter="\t")
                    if r.get("ID") and r["ID"].strip()]
    except (OSError, KeyError):
        return []


CUSTOM_CODES = _custom_codes()


# ---- consume the bundle ------------------------------------------------------
rule unpack_custom_bundle:
    """Extract the release custom bundle (proteomes/ + custom_busco.tsv). For now
    `custom_bundle` is a local tar.gz; a URL+sha256 handler can be added later."""
    input:
        CUSTOM_BUNDLE if CUSTOM_BUNDLE else [],
    output:
        proteomes=directory("results/custom_bundle/proteomes"),
        busco="results/custom_bundle/custom_busco.tsv",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
odir=$(dirname {output.busco})
mkdir -p "$odir"
tar -xzf {input} -C "$odir"
"""


# ---- produce the bundle (author-only curation) -------------------------------
# Raw curated proteomes live at data/custom_proteomes/CUS<id>.fa (author-placed,
# gitignored). BUSCO runs offline against resources/busco_db (as in new_genomes).
rule custom_busco:
    input:
        "data/custom_proteomes/{code}.fa",
    output:
        directory("results/custom_curation/busco/{code}"),
    log:
        "results/log/custom/{code}_busco.log",
    threads: 4
    conda:
        "../envs/busco.yaml"
    shell:
        """
busco -i {input} -l eukaryota -o {output} -m proteins \
--offline --download_path resources/busco_db -c {threads} > {log} 2>&1
"""


rule custom_busco_tsv:
    input:
        expand("results/custom_curation/busco/{code}", code=CUSTOM_CODES),
    output:
        "results/custom_curation/custom_busco.tsv",
    conda:
        "../envs/busco.yaml"
    shell:
        """
echo -e "file\\tcomplete\\tsingle\\tmulticopy\\tfragmented\\tmissing\\tn_markers" > {output}
cat {input}/*/*.json | jq -r \
'[.parameters.out, .results."Complete percentage", .results."Single copy percentage", \
.results."Multi copy percentage", .results."Fragmented percentage", \
.results."Missing percentage", .results.n_markers ] | @tsv' | \
sed 's|results/custom_curation/busco/||' >> {output}
"""


rule custom_stats:
    input:
        expand("data/custom_proteomes/{code}.fa", code=CUSTOM_CODES),
    output:
        "results/custom_curation/custom_stats.tsv",
    conda:
        "../envs/utils.yaml"
    shell:
        "seqkit stats -b -T {input} | sed 's/.fa//' > {output}"


rule package_custom_bundle:
    """Assemble the versioned release-custom bundle (proteomes + BUSCO + metadata
    + provenance) into one tar.gz. Upload it (eventually Zenodo) and set
    dbs.build.repdb.custom_bundle to its path/URL."""
    input:
        busco="results/custom_curation/custom_busco.tsv",
        stats="results/custom_curation/custom_stats.tsv",
        table=config["files"]["new_genomes"],
        proteomes=expand("data/custom_proteomes/{code}.fa", code=CUSTOM_CODES),
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
    gzip -c "$f" > "$stage/proteomes/$(basename "$f" .fa).faa.gz"
done
cp {input.busco} "$stage/custom_busco.tsv"
cp {input.stats} "$stage/custom_stats.tsv"
cp {input.table} "$stage/custom_metadata.tsv"
# provenance (source URLs per CUS id); move it under resources/ to version it
[ -f resources/custom_methods.txt ] && cp resources/custom_methods.txt "$stage/methods.txt" || true
mkdir -p $(dirname {output})
tar -czf {output} -C "$stage" .
rm -rf "$stage"
"""
