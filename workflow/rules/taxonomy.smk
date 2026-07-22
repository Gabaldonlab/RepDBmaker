def repdb_manifest():
    """Path to a frozen RepDB manifest, or None.

    A manifest is a previously produced ``repdb_taxonomy.tsv`` (ID + the seven
    tab-separated rank columns k,p,c,f,o,g,s). When set under
    ``dbs.build.repdb.manifest`` in the config, RepDB is *reproduced* from this
    exact composition instead of being re-selected from the live sources: the
    manifest replaces both the ``repdb_taxonomy`` checkpoint and the taxonomy
    that feeds ``create_full_taxdump``. This decouples the whole selection /
    taxonomy-harmonization phase from the drifting online metadata — the only
    remaining online access is fetching the sequences for the frozen IDs
    (composition-level, not bitwise, reproducibility).
    """
    repdb_conf = (config.get("dbs", {}).get("build", {}) or {}).get("repdb") or {}
    if isinstance(repdb_conf, dict):
        return repdb_conf.get("manifest")
    return None


REPDB_MANIFEST = repdb_manifest()


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


def validate_reference_input(wildcards):
    # In reproduce (manifest) mode there is no fresh eukaryotic taxonomy to
    # check against, so validation is schema-only; in normal mode the reference
    # (public eukaryotic taxonomy, no custom) enables the conflict check.
    if REPDB_MANIFEST:
        return []
    return rules.eukaryotes_taxonomy_ref.output


rule validate_custom_proteomes:
    """Validate the custom-proteome table; the run aborts on violations.

    Rejects (rather than silently dropping/mislabelling) rows with a non-CUS ID,
    a duplicate ID, or a lineage that is not exactly seven correctly prefixed
    non-empty ranks. When a reference eukaryotic taxonomy is available (normal
    mode) it also rejects lineages that place a known taxon under a parent that
    conflicts with the rest of the eukaryotic taxonomy.
    See workflow/scripts/check_custom_proteomes.R.
    """
    input:
        table=config["files"]["new_genomes"],
        reference=validate_reference_input,
    output:
        "results/meta/custom_proteomes_valid.txt",
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
        up=rules.uniprot_taxonomy.output,
        ep=rules.eukprot_taxonomy.output,
        p10k=rules.p10k_taxonomy.output,
        custom=rules.custom_taxonomy.output,
        up_stats=rules.get_uniprot_meta.output.stats,
        ep_stats=rules.get_eukprot.output.euk_busco,
        p10k_stats=rules.get_p10k.output.meta,
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


if REPDB_MANIFEST:

    rule create_full_taxdump:
        """Reproduce mode: build the taxdump directly from the frozen manifest.

        The manifest is already in ``all_taxonomy.tsv`` format (ID + 7 ranks),
        so no live taxonomy sources are needed.
        """
        input:
            manifest=REPDB_MANIFEST,
        output:
            all_taxa="results/taxonomies/all_taxonomy.tsv",
            full_taxdump=directory("results/taxdump/repdb_taxdump/"),
        conda:
            "../envs/utils.yaml"
        localrule: True
        shell:
            """
cp {input.manifest} {output.all_taxa}
taxonkit create-taxdump -A1 {output.all_taxa} --out-dir {output.full_taxdump} --force \
--rank-names "superkingdom","phylum","class","order","family","genus","species"
"""

else:

    rule create_full_taxdump:
        """Create taxdump from all genomes (used by all databases)"""
        input:
            # gtdb_proteomes=rules.decompress_gtdb_genomes.output,
            virus=rules.virus_taxonomy.output.tax,
            gtdb=rules.get_gtdb_tax.output.tax,
            all_euka=rules.eukaryotes_taxonomy.output,
        output:
            # ids="results/meta/all.ids",
            all_taxa="results/taxonomies/all_taxonomy.tsv",
            full_taxdump=directory("results/taxdump/repdb_taxdump/"),
        conda:
            "../envs/utils.yaml"
        localrule: True
        shell:
            """
cat {input.virus} {input.gtdb} {input.all_euka} | \
sed 's/d__//g' | sed 's/[a-z]__/\\t/g' | sed 's/;//g' > {output.all_taxa}
taxonkit create-taxdump -A1 {output.all_taxa} --out-dir {output.full_taxdump} --force \
--rank-names "superkingdom","phylum","class","order","family","genus","species"
"""


# find $(dirname {input.gtdb_proteomes}) -type f -name "*faa.gz" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort > {output.ids}


if REPDB_MANIFEST:

    checkpoint repdb_taxonomy:
        """Reproduce mode: use the frozen manifest as the RepDB composition."""
        input:
            manifest=REPDB_MANIFEST,
        output:
            "results/taxonomies/repdb_taxonomy.tsv",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            "cp {input.manifest} {output}"

else:

    checkpoint repdb_taxonomy:
        """Create RepDB-specific taxonomy with selected eukaryotes, all viruses and all GTDB genomes"""
        input:
            gtdb=rules.gtdb_species_clusters.output,
            virus=rules.virus_taxonomy.output.tax,
            filtered_euka=rules.select_repdb_eukaryotes.output.tax,
            full_table=rules.create_full_taxdump.output.all_taxa,
        output:
            "results/taxonomies/repdb_taxonomy.tsv",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
cat {input.gtdb} {input.virus} {input.filtered_euka} | cut -f1 | \
csvtk join -H -t - {input.full_table} > {output}
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
    """Create user-defined taxonomy for custom database"""
    input:
        mnemonics=get_custom_codes,
        full_table=rules.create_full_taxdump.output.all_taxa,
    output:
        "results/taxonomies/{db}_taxonomy.tsv",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
csvtk join -H -t {input.mnemonics} {input.full_table} > {output}
"""
