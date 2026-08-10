def repdb_universe():
    """Path to a pinned universe.tsv, or None.

    The universe is the enriched, versionable table of every available proteome
    (id, source, 7 ranks, completeness, annotated). When ``dbs.build.repdb.universe``
    is set, the universe is taken from this frozen file instead of being rebuilt
    from the live sources: the taxonomy harmonization is skipped, the taxdump is
    the full pinned universe, and the eukaryote selection re-runs deterministically
    on it.
    """
    repdb_conf = (config.get("dbs", {}).get("build", {}) or {}).get("repdb") or {}
    if isinstance(repdb_conf, dict):
        return repdb_conf.get("universe")
    return None


REPDB_UNIVERSE = repdb_universe()


def repdb_zenodo():
    """Pinned Zenodo record for repdb's heavy CONSTRUCTION artifacts, or None.

    When ``dbs.build.repdb.zenodo`` is set, `snakemake build` fetches the raw
    fasta, its maps, the decontamination report, the cluster membership table
    and the taxdump directly from Zenodo instead of assembling, decontaminating
    and clustering repdb from raw sources - a second, further-downstream seam
    than the universe pin above (that one still redoes all of that; this one
    skips it). See workflow/rules/zenodo_fetch.smk and docs/releasing.md.
    Independent of `dbs.build.repdb.universe`: nothing under `zenodo:` is
    derived from the universe, so Pipeline 1 (`snakemake sample`) isn't needed
    at all just to build repdb from a Zenodo-pinned release.

    Expected shape::

        dbs:
          build:
            repdb:
              zenodo:
                base_url: "https://zenodo.org/records/<record_id>/files"
                files:
                  fasta:            {name: repdb.fa.gz, sha256: "..."}
                  accession_map:    {name: repdb_accession_map.txt, sha256: "..."}
                  map:              {name: repdb.map, sha256: "..."}
                  nohead_map:       {name: repdb_nohead.map, sha256: "..."}
                  clusters:         {name: repdb_clusters.tsv, sha256: "..."}
                  contaminants_tsv: {name: repdb_contaminants.tsv, sha256: "..."}
                  contaminants_ids: {name: repdb_contaminants.txt, sha256: "..."}
                  taxdump:          {name: repdb_taxdump.tar.gz, sha256: "..."}
    """
    repdb_conf = (config.get("dbs", {}).get("build", {}) or {}).get("repdb") or {}
    if isinstance(repdb_conf, dict):
        z = repdb_conf.get("zenodo")
        return z if isinstance(z, dict) else None
    return None


REPDB_ZENODO = repdb_zenodo()


def _zenodo_file(key):
    """(url, sha256) for one pinned Zenodo asset under REPDB_ZENODO, or a
    clear config error. Shared by create_taxdump's zenodo branch below and by
    every fetch rule in zenodo_fetch.smk."""
    files = REPDB_ZENODO.get("files", {})
    entry = files.get(key)
    if not isinstance(entry, dict) or "name" not in entry or "sha256" not in entry:
        raise ValueError(
            f"dbs.build.repdb.zenodo.files.{key} is missing or incomplete "
            f"(need 'name' and 'sha256'; see repdb_zenodo() above)"
        )
    base_url = REPDB_ZENODO.get("base_url")
    if not base_url:
        raise ValueError("dbs.build.repdb.zenodo.base_url is not set")
    return f"{base_url}/{entry['name']}", entry["sha256"]


# ---- eukaryote downsampling parameters (dbs.build.repdb.downsample) ----------
# Defaults reproduce the original hard-coded RepDB behaviour so an absent
# `downsample:` block changes nothing.
DEFAULT_REDUCE_CLADES = [
    {"rank": "class", "taxon": "Opisthokonta", "n": 20},
    {"rank": "order", "taxon": "Ciliophora", "n": 20},
    {"rank": "family", "taxon": "Embryophyta", "n": 20},
]
_VALID_REDUCE_RANKS = {"phylum", "class", "order", "family"}


def _repdb_downsample():
    repdb = (config.get("dbs", {}).get("build", {}) or {}).get("repdb") or {}
    downsample = repdb.get("downsample") if isinstance(repdb, dict) else None
    return downsample if isinstance(downsample, dict) else {}


def _downsample_param(key, default):
    val = _repdb_downsample().get(key)
    return default if val is None else val


def _reduce_clades():
    """Validated list of {rank, taxon, n} clade-reduction entries."""
    entries = _downsample_param("reduce_abundant_clades", DEFAULT_REDUCE_CLADES) or []
    if not isinstance(entries, list):
        raise ValueError(
            "dbs.build.repdb.select.reduce_abundant_clades must be a list of "
            "{rank, taxon, n} entries."
        )
    for e in entries:
        if not isinstance(e, dict) or {"rank", "taxon", "n"} - e.keys():
            raise ValueError(
                f"Invalid reduce_abundant_clades entry {e!r}: needs rank, taxon, n."
            )
        if str(e["rank"]).lower() not in _VALID_REDUCE_RANKS:
            raise ValueError(
                f"reduce_abundant_clades: rank '{e['rank']}' must be one of "
                f"{sorted(_VALID_REDUCE_RANKS)}."
            )
    return entries


# validate the downsampling config once, at parse time
_reduce_clades()


rule virus_taxonomy:
    input:
        meta=rules.get_virus_genomes.output.meta,
        folder=rules.get_virus_genomes.output.folder,
        # declared as input (not params) so Snakemake enforces the dependency
        # on download_taxdump and avoids a race condition
        taxdump=rules.download_taxdump.output,
    output:
        ids="results/meta/refseq_virus_ids.txt",
        tax="results/taxonomies/virus_taxonomy.tsv",
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        """
unzip -l {input.folder} | awk '{{print $NF}}' | grep protein | cut -f3 -d'/' > {output.ids}
cut -f1,6 {input.meta} | grep -F -w -f {output.ids} | sed 's/_//' | sed 's/\\..*\\t/\\t/g' | \
taxonkit reformat -I 2 -P -F -p "unclassified_" -s "" --data-dir {input.taxdump} | cut -f1,3 | \
sed 's/k__unclassified_Viruses/d__Viruses/g' | \
awk -F'\\t' '$2 ~ /[^;]/' > {output.tax}
"""


rule p10k_taxonomy:
    input:
        lineage=rules.get_p10k.output.lineage,
        unieuk=rules.download_unieuk.output,
        ep=rules.get_eukprot.output.euk_included,
        notep=rules.get_eukprot.output.euk_excluded,
    output:
        "results/taxonomies/p10k_taxonomy.tsv",
    conda:
        "../envs/R.yaml"
    localrule: True
    script:
        "../scripts/get_tax.R"


rule uniprot_taxonomy:
    input:
        lineage=rules.get_uniprot_meta.output.lineage,
        unieuk=rules.download_unieuk.output,
        ep=rules.get_eukprot.output.euk_included,
        notep=rules.get_eukprot.output.euk_excluded,
    output:
        "results/taxonomies/uniprot_eukprot_taxonomy.tsv",
    conda:
        "../envs/R.yaml"
    localrule: True
    script:
        "../scripts/get_tax.R"


rule eukprot_taxonomy:
    input:
        rules.get_eukprot.output.euk_included,
    output:
        "results/taxonomies/eukprot_taxonomy.tsv",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
awk 'NR>1' {input} | sed 's/\\_/ /g' | \
awk -F'\\t' '{{split($11, a, ";"); print $1,"d__"a[1]";p__"a[2]";c__"$8";o__"$9";f__"$10";g__"$6";s__"$2}}' OFS='\\t' | \
sed  's/N\\/A//g' | tr -d \\'\\" | sed 's/\\ ;/;/g' | sed 's/other Gyrista/other_Gyrista/g'> {output}
"""


rule eukaryotes_taxonomy_ref:
    """Eukaryotic taxonomy from the public sources only (no custom proteomes).

    Serves as the reference against which custom lineages are checked for
    conflicts (a known taxon placed under a different parent).
    """
    input:
        up=rules.uniprot_taxonomy.output,
        ep=rules.eukprot_taxonomy.output,
        p10k=rules.p10k_taxonomy.output,
    output:
        "results/taxonomies/eukaryotes_taxonomy_ref.tsv",
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        "cat {input} | sort -k2,2 > {output}"


# Per-domain references for the custom-proteome conflict check. In pinned-universe
# mode these rules are simply not in the DAG (harmonization is skipped). The script
# only *parses* the reference for a domain that actually occurs in the custom
# table, so the huge GTDB taxonomy is not read unless a prokaryotic custom row
# exists.
def validate_reference_input(wildcards):
    return rules.eukaryotes_taxonomy_ref.output


def validate_reference_prok(wildcards):
    return rules.get_gtdb_tax.output.tax


def validate_reference_virus(wildcards):
    return rules.virus_taxonomy.output.tax


def _custom_fasta_inputs(wildcards):
    """Local per-id FASTA files to content-check, one per CUS id - only when NOT
    using a pinned release bundle (CUSTOM_BUNDLE/custom_proteome_path come from
    custom_bundle.smk, included later, but this resolves fine: Snakemake calls
    input functions only after the whole Snakefile, all includes, is loaded).
    A bundle's proteomes are validated once at packaging time via
    `package_custom` instead, since they are not local per-id files here."""
    if CUSTOM_BUNDLE:
        return []
    return [custom_proteome_path(c) for c in CUSTOM_CODES]


rule validate_custom_proteomes:
    """Validate the custom-proteome table; the run aborts on violations.

    Rejects (rather than silently dropping/mislabelling) rows with a non-CUS ID,
    a duplicate ID, or a lineage that is not exactly seven correctly prefixed
    non-empty ranks. Custom proteomes may be eukaryotic OR not; in normal mode
    each row is also checked against the reference taxonomy for its domain
    (Eukaryota, GTDB for Bacteria/Archaea, virus taxonomy for Viruses) and
    rejected if it places a known taxon under a conflicting parent. Also
    content-checks each row's FASTA file (non-empty, has records, no embedded
    whitespace or other non-residue characters in sequence lines) - the same
    class of defect that otherwise only surfaces much later as diamond's
    "Invalid character in sequence: ' '", deep inside the final RepDB FASTA.
    See workflow/scripts/check_custom_proteomes.R.
    """
    input:
        table=config["files"]["new_genomes"],
        reference=validate_reference_input,
        reference_prok=validate_reference_prok,
        reference_virus=validate_reference_virus,
        fasta=_custom_fasta_inputs,
    output:
        "results/meta/custom_proteomes_valid.txt",
    params:
        check_fasta=True,
    conda:
        "../envs/R.yaml"
    localrule: True
    script:
        "../scripts/check_custom_proteomes.R"


rule custom_taxonomy:
    input:
        table=config["files"]["new_genomes"],
        # gate: fails the run before any custom lineage is propagated
        valid=rules.validate_custom_proteomes.output,
    output:
        "results/taxonomies/custom_taxonomy.tsv",
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        """
csvtk cut -t -f ID,Lineage {input.table} | awk 'NR>1' > {output}
"""


rule eukaryotes_taxonomy:
    """Full eukaryotic taxonomy = public reference + validated custom proteomes."""
    input:
        reference=rules.eukaryotes_taxonomy_ref.output,
        custom=rules.custom_taxonomy.output,
    output:
        "results/taxonomies/eukaryotes_taxonomy.tsv",
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        "cat {input.reference} {input.custom} | sort -k2,2 > {output}"


rule select_repdb_eukaryotes:
    input:
        # the eukaryote downsampling reads the universe (ranks + completeness), so
        # it re-runs deterministically from a pinned universe. Literal path to
        # decouple from include order (build_universe lives in universe.smk).
        universe="results/universe/universe.tsv",
        exclude=config["files"]["genomes_to_exclude"],
        to_keep=config["files"]["clades_to_keep"],
    output:
        tax="results/taxonomies/selected_eukaryotes.tsv",
        plot="results/plots/euka_db.pdf",
    params:
        remove_duplicated_species=lambda wildcards: _downsample_param(
            "remove_duplicated_species", True
        ),
        top_n_genuses=lambda wildcards: _downsample_param("top_n_genuses", 3),
        # pass the clade-reduction table as three parallel lists so the nested
        # structure survives serialization into the R `snakemake` object.
        reduce_ranks=lambda wildcards: [e["rank"] for e in _reduce_clades()],
        reduce_taxa=lambda wildcards: [e["taxon"] for e in _reduce_clades()],
        reduce_ns=lambda wildcards: [int(e["n"]) for e in _reduce_clades()],
    conda:
        "../envs/R.yaml"
    localrule: True
    # log: "results/log/db/select_euka.log" # it does not automatically redirect to stderr.....
    script:
        "../scripts/filter_tax.R"


rule create_full_taxdump:
    """The full taxonomy of every available proteome (all sources), id + 7 ranks.
    Feeds build_universe; only built in normal mode (a pinned universe skips it)."""
    input:
        virus=rules.virus_taxonomy.output.tax,
        gtdb=rules.get_gtdb_tax.output.tax,
        all_euka=rules.eukaryotes_taxonomy.output,
    output:
        all_taxa="results/taxonomies/all_taxonomy.tsv",
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        """
cat {input.virus} {input.gtdb} {input.all_euka} | \
sed 's/d__//g' | sed 's/[a-z]__/\\t/g' | sed 's/;//g' > {output.all_taxa}
"""


if REPDB_ZENODO:

    rule create_taxdump:
        """Fetch the taxdump a Zenodo-pinned repdb release was built with
        (dbs.build.repdb.zenodo), instead of regenerating it from the
        universe. This has to be a branch on `create_taxdump` itself, not a
        competing ruleorder'd rule (like the other zenodo_fetch.smk fetches):
        its output has no wildcards, so make_diamonddb/make_mmseqsdb/get_clade
        - which reference it as `rules.create_taxdump.output.full_taxdump` -
        bind to this rule directly at parse time, bypassing ruleorder
        entirely. See zenodo_fetch.smk for the fasta/map/cluster equivalents,
        which ARE wildcarded and so can use ruleorder normally. The fetched
        taxdump is already taxid-compacted (see the else branch below and
        compact_taxdump_ids.py) - that happened once, at build time, before
        it was uploaded."""
        output:
            full_taxdump=directory("results/taxdump/repdb_taxdump/"),
        params:
            url=lambda wc: _zenodo_file("taxdump")[0],
            sha=lambda wc: _zenodo_file("taxdump")[1],
        log:
            "results/log/downloads/zenodo_repdb_taxdump.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
tmp=$(mktemp --suffix .tar.gz)
bash workflow/scripts/fetch_zenodo.sh {params.url} {params.sha} "$tmp" >{log} 2>&1
mkdir -p {output.full_taxdump}
tar -xzf "$tmp" -C {output.full_taxdump}
rm -f "$tmp"
"""

else:

    rule create_taxdump:
        """Build the taxdump from the UNIVERSE (id + 7 ranks), so it always covers
        the full set of available proteomes and is defined by the versioned universe
        artifact rather than by whichever subset a given database selects. taxids are
        hashed per-lineage by taxonkit, so they are identical to a fresh full run.
        """
        input:
            universe="results/universe/universe.tsv",
        output:
            full_taxdump=directory("results/taxdump/repdb_taxdump/"),
        conda:
            "../envs/utils.yaml"
        localrule: True
        shell:
            """
tmp=$(mktemp)
awk -F'\\t' 'NR>1' {input.universe} | cut -f1,3-9 > "$tmp"
taxonkit create-taxdump -A1 "$tmp" --out-dir {output.full_taxdump} --force \
--rank-names "superkingdom","phylum","class","order","family","genus","species"
rm -f "$tmp"
python3 workflow/scripts/compact_taxdump_ids.py {output.full_taxdump} {output.full_taxdump}/taxid.map
"""


# find $(dirname {input.gtdb_proteomes}) -type f -name "*faa.gz" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort > {output.ids}


checkpoint repdb_taxonomy:
    """RepDB composition = all prokaryotes + all viruses + the selected
    eukaryotes, taken as a subset of the universe (id + 7 ranks). Works
    identically for a freshly-built or a pinned universe."""
    input:
        universe="results/universe/universe.tsv",
        selected=rules.select_repdb_eukaryotes.output.tax,
    output:
        "results/taxonomies/repdb_taxonomy.tsv",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
keep=$(mktemp)
{{ cut -f1 {input.selected}; \
   awk -F'\\t' 'NR>1 && ($2=="gtdb" || $2=="virus"){{print $1}}' {input.universe}; }} \
   | sort -u > "$keep"
awk -F'\\t' 'NR==FNR{{k[$1];next}} FNR>1 && ($1 in k){{print $1"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9}}' \
   "$keep" {input.universe} > {output}
rm -f "$keep"
"""


# cat <(shuf -n 50 {input.gtdb}) <(shuf -n 50 {input.virus}) <(shuf -n 50 {input.filtered_euka}) | cut -f1 | \

# cat {input.gtdb} {input.virus} $(shuf -n 20 {input.filtered_euka}) | cut -f1 | \


def get_custom_codes(wildcards):
    """Return the path to the genome table for a custom database"""
    custom_dbs = config.get("dbs", {}).get("build", {}).get("custom", {})
    # print(f"Custom DBs in config: {custom_dbs}")
    if isinstance(custom_dbs, dict) and wildcards.db in custom_dbs:
        return custom_dbs[wildcards.db]["ids"]
    raise ValueError(f"Custom database '{wildcards.db}' not found in config")


# both checkpoints can produce results/taxonomies/repdb_taxonomy.tsv (the custom
# rule via its {db} wildcard); prefer the dedicated RepDB checkpoint for it.
ruleorder: repdb_taxonomy > custom_dbs_taxonomy


checkpoint custom_dbs_taxonomy:
    """Custom-database taxonomy: the requested ids joined to their ranks in the
    universe (id + 7 ranks). Works with a freshly-built or a pinned universe."""
    input:
        mnemonics=get_custom_codes,
        universe="results/universe/universe.tsv",
    output:
        "results/taxonomies/{db}_taxonomy.tsv",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
awk -F'\\t' 'NR>1' {input.universe} | cut -f1,3-9 | \
csvtk join -H -t {input.mnemonics} - > {output}
"""
