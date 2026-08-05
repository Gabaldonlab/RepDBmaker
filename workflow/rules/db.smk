def _is_db_clustered(db_name):
    dbs_config = config.get("dbs", {}).get("build", {})
    if db_name == "repdb":
        db_conf = dbs_config.get("repdb")
    else:
        db_conf = dbs_config.get("custom", {}).get(db_name)
    if isinstance(db_conf, dict):
        return bool(db_conf.get("cluster"))
    return False


def _fa_raw(db):
    return f"results/dbs/{db}/{db}.fa.gz"


def _fa_decontaminated(db):
    return f"results/dbs/{db}/{db}_decontaminated.fa.gz"


def _fa_clustered(db):
    return f"results/dbs/{db}/{db}_clustered.fa.gz"


def get_cluster_input_fa(wildcards):
    """FASTA fed into clustering: the decontaminated one when the db is
    decontaminated (decontamination runs on the raw set first, so contaminant
    detection keeps full sensitivity), otherwise the raw assembled FASTA.
    Transform chain: raw -> [decontaminate] -> [cluster] -> final."""
    db = wildcards.db
    return _fa_decontaminated(db) if _is_db_decon(db) else _fa_raw(db)


def _final_db_fa(db):
    """FASTA the search indices are built from = the last enabled stage of the
    chain: clustered if the db is clustered, else decontaminated if it is
    decontaminated, else the raw assembled FASTA."""
    if _is_db_clustered(db):
        return _fa_clustered(db)
    if _is_db_decon(db):
        return _fa_decontaminated(db)
    return _fa_raw(db)


def get_final_db_fa(wildcards):
    """Input-function form of `_final_db_fa` (`_is_db_decon` lives in
    decontaminate.smk; both flags resolve at DAG time)."""
    return _final_db_fa(wildcards.db)


mmseqs_ext = [
    ".dbtype",
    "_h",
    "_h.dbtype",
    "_h.index",
    ".lookup",
    ".index",
    "_mapping",
    ".source",
    "_taxonomy",
]


rule make_db_fasta:
    input:
        table="results/dbs/{db}/genome_table.tsv",
        taxdump=rules.create_taxdump.output.full_taxdump,
        stats="results/stats/{db}_stats.tsv",  # if this failed it means some proteomes had problems while downloading!
    output:
        # fa=temp("results/dbs/{db}/{db}_raw.fa.gz"),
        fa="results/dbs/{db}/{db}.fa.gz",
        idmap="results/dbs/{db}/{db}_accession_map.txt",
    # log: "results/log/dbs/repdb/parse.log"
    threads: 112
    benchmark:
        "results/benchmarks/dbs/{db}/parse.txt"
    # group: "create_db"
    conda:
        "../envs/python.yaml"
    # group: "create_db"
    script:
        "../scripts/parse_gnm_v3.py"


rule make_db_map:
    input:
        idmap=rules.make_db_fasta.output.idmap,
    output:
        headermap="results/dbs/{db}/{db}.map",
        noheadermap="results/dbs/{db}/{db}_nohead.map",
    # localrule: True
    # group: "create_db"
    conda:
        "../envs/utils.yaml"
    shell:
        """
echo -e "accession.version\\ttaxid" > {output.headermap}
cut -f2 {input.idmap} | awk -F'\\t' '{{split($NF, a, "_"); print $0"\\t"a[1]}}' >> {output.headermap}
awk 'NR>1' {output.headermap} > {output.noheadermap}
"""


# using mmseqs was the fastest way I found to extract clade specific fastas.
rule make_mmseqsdb_clustering:
    input:
        # decontaminated fasta when the db is decontaminated, else the raw one
        # (see get_cluster_input_fa): clustering is the last transform, so it
        # consumes whatever decontamination produced.
        fa=get_cluster_input_fa,
        taxidmap=rules.make_db_map.output.noheadermap,
        taxdump=rules.create_taxdump.output.full_taxdump,
    output:
        db=temp("results/dbs/{db}/cluster/{db}_mmseqs"),
        db_extra=temp(
            expand("results/dbs/{{db}}/cluster/{{db}}_mmseqs{ext}", ext=mmseqs_ext)
        ),
    log:
        "results/log/dbs/{db}/make_mmseqs_clustering.log",
    benchmark:
        "results/benchmarks/dbs/{db}/make_mmseqs_clustering.txt"
    threads: 112
    # group: "create_db"
    conda:
        "../envs/homology.yaml"
    shell:
        """ 
mmseqs createdb {input.fa} {output.db} > {log}
mmseqs createtaxdb {output.db} $TMPDIR --ncbi-tax-dump {input.taxdump} \
--tax-mapping-file {input.taxidmap} --threads {threads} >> {log}
"""


# this is done in order to set the unclustered fasta as tmp
# rule unify_fasta:
#     input:
#         get_db_fa,
#     output:
#         "results/dbs/{db}/{db}.fa.gz",
#     localrule: True
#     conda:
#         "../envs/utils.yaml"
#     shell:
#         "cp {input} {output}"


rule make_blastdb:
    input:
        fa=get_final_db_fa,
        taxidmap=rules.make_db_map.output.noheadermap,
    output:
        "results/dbs/{db}/{db}_blastp",
    log:
        "results/log/dbs/{db}/make_blastp.log",
    benchmark:
        "results/benchmarks/dbs/{db}/make_blastp.txt"
    conda:
        "../envs/homology.yaml"
    # group: "create_db"
    shell:
        """
gunzip -c {input.fa} | makeblastdb -in - -parse_seqids -taxid_map {input.taxidmap} \
-dbtype prot -out {output} -title {wildcards.db} -logfile {log}
touch {output}
"""


rule make_diamonddb:
    input:
        fa=get_final_db_fa,
        taxidmap="results/dbs/{db}/{db}.map",
        taxdump=rules.create_taxdump.output.full_taxdump,
    output:
        "results/dbs/{db}/{db}_diamond",
    log:
        "results/log/dbs/{db}/make_diamond.log",
    benchmark:
        "results/benchmarks/dbs/{db}/make_diamond.txt"
    threads: 48
    conda:
        "../envs/homology.yaml"
    # group: "create_db"
    shell:
        """
diamond makedb --in {input.fa} -d {output} --threads {threads} \
--taxonnodes {input.taxdump}/nodes.dmp --taxonmap {input.taxidmap} 2> {log}
touch {output}
"""
# --taxonnodes {input.taxdump}/nodes.dmp --taxonnames {input.taxdump}/names.dmp


rule make_mmseqsdb:
    input:
        fa=get_final_db_fa,
        taxidmap="results/dbs/{db}/{db}_nohead.map",
        taxdump=rules.create_taxdump.output.full_taxdump,
    output:
        "results/dbs/{db}/{db}_mmseqs",
    log:
        "results/log/dbs/{db}/make_mmseqs.log",
    benchmark:
        "results/benchmarks/dbs/{db}/make_mmseqs.txt"
    threads: 112
    conda:
        "../envs/homology.yaml"
    # group: "create_db"
    shell:
        """
mmseqs createdb {input.fa} {output} > {log}
mmseqs createtaxdb {output} $TMPDIR --ncbi-tax-dump {input.taxdump} \
--tax-mapping-file {input.taxidmap} --threads {threads} >> {log}
"""
