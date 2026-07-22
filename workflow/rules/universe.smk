# The UNIVERSE seam between the two pipelines (see docs/split-pipeline.md).
#
# Pipeline 1 (sampling/curation) produces:
#   results/universe/universe.tsv  - id, ranks, completeness, source, data_type
#   results/universe/repdb.ids     - the selected RepDB composition (id list)
# Pipeline 2 (construction) consumes them to build the taxdump and the databases.
#
# Freezing universe.tsv per release freezes taxonomy + BUSCO, so a rebuild
# reproduces the release deterministically without re-running harmonization.


if REPDB_UNIVERSE:

    rule build_universe:
        """Pinned universe (dbs.build.repdb.universe): use the frozen file
        directly. Taxonomy harmonization is skipped; the taxdump and the
        eukaryote selection are (re)built deterministically from this file."""
        input:
            REPDB_UNIVERSE,
        output:
            "results/universe/universe.tsv",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            "cp {input} {output}"

else:

    rule build_universe:
        """Assemble the enriched universe from the harmonized taxonomy + BUSCO
        stats: the full set of available proteomes (every source) with the ranks
        and completeness the selection and taxdump need."""
        input:
            tax=rules.create_full_taxdump.output.all_taxa,
            up_stats=rules.get_uniprot_meta.output.stats,
            ep_stats=rules.get_eukprot.output.euk_busco,
            ep=rules.get_eukprot.output.euk_included,
            p10k_stats=rules.get_p10k.output.meta,
            gtdb_meta=rules.get_gtdb_tax.output.meta,
            custom_table=config["files"]["new_genomes"],
            custom_busco=_custom_busco,
        output:
            "results/universe/universe.tsv",
        conda:
            "../envs/R.yaml"
        localrule: True
        script:
            "../scripts/build_universe.R"


rule emit_repdb_ids:
    """Emit the RepDB selection (id list) chosen from the universe.

    Composition = all prokaryotes + all viruses + the selected eukaryotes; this
    is exactly the id column of the repdb_taxonomy checkpoint.
    """
    input:
        "results/taxonomies/repdb_taxonomy.tsv",
    output:
        "results/universe/repdb.ids",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        "cut -f1 {input} | sort -u > {output}"


def _repdb_ids_reference(wildcards):
    # optional committed selection to diff the freshly-computed one against
    ref = config.get("files", {}).get("repdb_ids_reference")
    return [ref] if ref else []


rule repdb_ids_diff:
    """Soft, non-fatal report: how the freshly-selected RepDB composition differs
    from a committed reference id list (added / removed proteomes). Never gates
    the build - it is provenance, not validation."""
    input:
        current=rules.emit_repdb_ids.output,
        reference=_repdb_ids_reference,
    output:
        "results/qc/repdb_ids_diff.txt",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
if [ -n "{input.reference}" ] && [ -s "{input.reference}" ]; then
    sort -u "{input.reference}" > .ref.$$
    added=$(comm -13 .ref.$$ {input.current} | wc -l)
    removed=$(comm -23 .ref.$$ {input.current} | wc -l)
    {{
      echo "RepDB composition vs reference {input.reference}"
      echo "  added (in current, not reference): $added"
      echo "  removed (in reference, not current): $removed"
      echo "--- added ---";   comm -13 .ref.$$ {input.current}
      echo "--- removed ---"; comm -23 .ref.$$ {input.current}
    }} > {output}
    rm -f .ref.$$
else
    echo "No reference id list configured (files.repdb_ids_reference); nothing to diff." > {output}
fi
"""
