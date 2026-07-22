# Interactive Krona charts, one selectable dataset per database/taxonomy. For
# UniProt and P10K both the source schema and the UniEuk-harmonized version are
# charted, so the harmonization can be compared. Each entry is (file, format).
KRONA_TAX = {
    # name: (file, format, root)  -- root fills an empty domain (UniProt NCBI)
    "gtdb": (rules.get_gtdb_tax.output.tax, "prefixed", ""),
    "virus": (rules.virus_taxonomy.output.tax, "prefixed", ""),
    "eukprot": (rules.eukprot_taxonomy.output, "prefixed", ""),
    "uniprot_ncbi": (rules.get_uniprot_meta.output.ncbi_tax, "prefixed", "Eukaryota"),
    "uniprot_unieuk": (rules.uniprot_taxonomy.output, "prefixed", ""),
    "p10k_native": (rules.get_p10k.output.lineage, "p10k_native", ""),
    "p10k_unieuk": (rules.p10k_taxonomy.output, "prefixed", ""),
    "custom": (rules.custom_taxonomy.output, "prefixed", ""),
}


rule krona_text:
    """Convert one taxonomy table into KronaTools' ktImportText text format."""
    input:
        tax=lambda wc: KRONA_TAX[wc.name][0],
        functions="workflow/scripts/functions.R",
    output:
        "results/qc/krona/{name}.krona.txt",
    params:
        fmt=lambda wc: KRONA_TAX[wc.name][1],
        root=lambda wc: KRONA_TAX[wc.name][2],
    conda:
        "../envs/R.yaml"
    localrule: True
    script:
        "../scripts/taxonomy_to_krona.R"


rule krona_plot:
    """Interactive Krona chart with one selectable dataset per database/taxonomy."""
    input:
        texts=expand("results/qc/krona/{name}.krona.txt", name=list(KRONA_TAX)),
    output:
        "results/qc/krona/repdb_krona.html",
    params:
        specs=lambda wc, input: " ".join(
            f"{f},{name}" for name, f in zip(KRONA_TAX, input.texts)
        ),
    conda:
        "../envs/krona.yaml"
    localrule: True
    shell:
        """
# include only non-empty datasets (e.g. custom may be empty)
args=""
for spec in {params.specs}; do
    f="${{spec%%,*}}"
    [ -s "$f" ] && args="$args $spec"
done
ktImportText -o {output} $args
"""


rule taxonomy_harmonization_report:
    """QC report of what changes when the UniProt (NCBI) and
    P10K (native) taxonomies are translated into the UniEuk framework - coverage,
    phylum relabelling, genus changes and proteomes mapped to more than one
    UniEuk lineage (ambiguous). 
    """
    input:
        functions="workflow/scripts/functions.R",
        palettes="workflow/scripts/palettes.R",
        eukprot=rules.eukprot_taxonomy.output,
        uniprot_ncbi=rules.get_uniprot_meta.output.ncbi_tax,
        uniprot_unieuk=rules.uniprot_taxonomy.output,
        p10k_orig=rules.get_p10k.output.lineage,
        p10k_unieuk=rules.p10k_taxonomy.output,
        custom=rules.custom_taxonomy.output,
    output:
        # a .Rmd script rule may declare only the rendered document; the
        # per-source TSVs are written by the notebook into the same directory
        html="results/qc/harmonization_report.html",
    conda:
        "../envs/R.yaml"
    localrule: True
    script:
        "../notebooks/harmonization_report.Rmd"


rule db_stats:
    input:
        "results/dbs/{db}/genome_table.tsv",
    output:
        "results/stats/{db}_stats.tsv",
    threads: 48
    conda:
        "../envs/utils.yaml"
    # localrule: True
    # group: "create_db"
    shell:
        """
cut -f1 {input} | seqkit stats -j {threads} --infile-list - -T -b | \
cut -f1,4- | sed 's/.faa.gz//' > {output}
"""


rule clustdb_stats:
    input:
        repr_all_clades,
    output:
        "results/stats/{db}_clustered_stats.tsv",
    threads: 48
    conda:
        "../envs/utils.yaml"
    # group: "create_db"
    # localrule: True
    shell:
        """
seqkit stats -j {threads} -T -b {input} | \
cut -f1,4- | sed 's/_rep_seq.fasta//' > {output}
"""


def _db_meta_tax(wildcards):
    # each database's composition taxonomy: repdb_taxonomy.tsv for RepDB,
    # <db>_taxonomy.tsv for a custom database (both are ID + 7 ranks).
    return f"results/taxonomies/{wildcards.db}_taxonomy.tsv"


def _db_meta_contaminants(wildcards):
    # RepDBmaker's per-organism decontamination flagging, only when this database
    # is decontaminated (the script degrades gracefully without it).
    if _is_db_decon(wildcards.db):
        return f"results/dbs/{wildcards.db}/decontaminate/contaminants.txt"
    return []


def _custom_busco(wildcards):
    # optional user-provided BUSCO scores for custom proteomes
    path = config.get("files", {}).get("new_genomes_busco")
    return [path] if path else []


rule make_db_meta:
    """Per-organism provenance metadata for a database (RepDB or a custom db),
    written to results/meta/<db>_meta.tsv. Includes the taxdump taxid.
    """
    input:
        tax=_db_meta_tax,
        stats="results/stats/{db}_stats.tsv",
        taxdump=rules.create_full_taxdump.output.full_taxdump,
        up_stats=rules.get_uniprot_meta.output.stats,
        ep_stats=rules.get_eukprot.output.euk_busco,
        ep=rules.get_eukprot.output.euk_included,
        p10k_stats=rules.get_p10k.output.meta,
        gtdb_meta=rules.get_gtdb_tax.output.meta,
        custom_table=config["files"]["new_genomes"],
        custom_busco=_custom_busco,
        contaminants=_db_meta_contaminants,
    output:
        "results/meta/{db}_meta.tsv",
    conda:
        "../envs/R.yaml"
    localrule: True
    script:
        "../scripts/analyze_stats.R"
