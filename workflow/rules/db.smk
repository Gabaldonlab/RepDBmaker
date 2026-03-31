def _is_db_clustered(db_name):
    dbs_config = config.get("dbs", {}).get("build", {})
    if db_name == "repdb":
        db_conf = dbs_config.get("repdb")
    else:
        db_conf = dbs_config.get("custom", {}).get(db_name)
    if isinstance(db_conf, dict):
        return bool(db_conf.get("cluster"))
    return False


def get_db_fa(wildcards):
    if _is_db_clustered(wildcards.db):
        return f"results/dbs/{wildcards.db}/{wildcards.db}_clustered.fa.gz"
    return f"results/dbs/{wildcards.db}/{wildcards.db}.fa.gz"

mmseqs_ext = [".dbtype", "_h", "_h.dbtype", "_h.index", 
              ".lookup", ".index", "_mapping", ".source", "_taxonomy"]

rule make_db_fasta:
    input:
        table="results/dbs/{db}/genome_table.tsv",
        taxdump=rules.create_full_taxdump.output.full_taxdump,
        stats="results/stats/{db}_stats.tsv" # if this failed it means some proteomes had problems while downloading!
    output:
        fa="results/dbs/{db}/{db}.fa.gz",
        idmap="results/dbs/{db}/{db}_accession_map.txt"
    # log: "results/log/dbs/repdb/parse.log"
    threads: 112
    benchmark: "results/benchmarks/dbs/{db}/parse.txt"
    # group: "create_db"
    conda: "../envs/python.yaml"
    # group: "create_db"
    script: "../scripts/parse_gnm_v2.py"


rule make_db_map:
    input:
        idmap=rules.make_db_fasta.output.idmap
    output:
        headermap="results/dbs/{db}/{db}.map",
        noheadermap="results/dbs/{db}/{db}_nohead.map"
    # localrule: True
    # group: "create_db"
    shell:'''
echo -e "accession.version\\ttaxid" > {output.headermap}
cut -f2 {input.idmap} | awk -F'\\t' '{{split($NF, a, "_"); print $0"\\t"a[1]}}' >> {output.headermap}
awk 'NR>1' {output.headermap} > {output.noheadermap}
'''

# using mmseqs was the fastest way I found to extract clade specific fastas.
rule make_mmseqsdb_clustering:
    input:
        fa=rules.make_db_fasta.output.fa,
        taxidmap=rules.make_db_map.output.noheadermap,
        taxdump=rules.create_full_taxdump.output.full_taxdump
    output: 
        db=temp("results/dbs/{db}/cluster/{db}_mmseqs"),
        db_extra=temp(expand("results/dbs/{{db}}/cluster/{{db}}_mmseqs{ext}", ext=mmseqs_ext))
    log: "results/log/dbs/{db}/make_mmseqs_clustering.log"
    benchmark: "results/benchmarks/dbs/{db}/make_mmseqs_clustering.txt"
    threads: 112
    # group: "create_db"
    conda: "../envs/homology.yaml"
    shell:''' 
mmseqs createdb {input.fa} {output.db} > {log}
mmseqs createtaxdb {output.db} $TMPDIR --ncbi-tax-dump {input.taxdump} \
--tax-mapping-file {input.taxidmap} --threads {resources.cpus_per_task} >> {log}
'''

rule make_blastdb:
    input: 
        fa=get_db_fa,
        taxidmap=rules.make_db_map.output.noheadermap
    output: "results/dbs/{db}/{db}_blastp"
    log: "results/log/dbs/{db}/make_blastp.log"
    benchmark: "results/benchmarks/dbs/{db}/make_blastp.txt"
    conda: "../envs/homology.yaml"
    # group: "create_db"
    shell:'''
gunzip -c {input.fa} | makeblastdb -in - -parse_seqids -taxid_map {input.taxidmap} \
-dbtype prot -out {output} -title {wildcards.db} -logfile {log}
touch {output}
'''

rule make_diamonddb:
    input: 
        fa=get_db_fa,
        taxidmap="results/dbs/{db}/{db}.map",
        taxdump=rules.create_full_taxdump.output.full_taxdump
    output: "results/dbs/{db}/{db}_diamond"
    log: "results/log/dbs/{db}/make_diamond.log"
    benchmark: "results/benchmarks/dbs/{db}/make_diamond.txt"
    threads: 48
    conda: "../envs/homology.yaml"
    # group: "create_db"
    shell:'''
diamond makedb --in {input.fa} -d {output} --threads {resources.cpus_per_task} \
--taxonnodes {input.taxdump}/nodes.dmp --taxonmap {input.taxidmap} 2> {log}
touch {output}
'''
# --taxonnodes {input.taxdump}/nodes.dmp --taxonnames {input.taxdump}/names.dmp

rule make_mmseqsdb:
    input: 
        fa=get_db_fa,
        taxidmap="results/dbs/{db}/{db}_nohead.map",
        taxdump=rules.create_full_taxdump.output.full_taxdump
    output: "results/dbs/{db}/{db}_mmseqs"
    log: "results/log/dbs/{db}/make_mmseqs.log"
    benchmark: "results/benchmarks/dbs/{db}/make_mmseqs.txt"
    threads: 112
    conda: "../envs/homology.yaml"
    # group: "create_db"
    shell:'''
mmseqs createdb {input.fa} {output} > {log}
mmseqs createtaxdb {output} $TMPDIR --ncbi-tax-dump {input.taxdump} \
--tax-mapping-file {input.taxidmap} --threads {resources.cpus_per_task} >> {log}
'''

