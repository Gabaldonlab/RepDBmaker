configfile: "config/benchmark.yaml"


rule all:
    input:
        "results/comparison/report.html",


rule get_random_seqs:
    output:
        "results/comparison/input.fa",
    params:
        n_seqs=config["benchmark"]["n_seqs"],
        tax_annot=config["benchmark"]["tax_annot"],
        source_seqs=config["benchmark"]["source_seqs"],
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        """
shuf -n {params.n_seqs} {params.tax_annot} | cut -f2 | seqkit grep -f - {params.source_seqs} > {output}
"""


rule diamond_random_seqs:
    input:
        rules.get_random_seqs.output,
    output:
        "results/comparison/hits/{db}_matches.tsv",
    params:
        db=lambda wcs: str(config["benchmark"]["dbs"][wcs.db]["dmnd"]),
    conda:
        "../envs/benchmark.yaml"
    log:
        "results/log/comparison/{db}_dmnd.log",
    benchmark:
        "results/benchmarks/comparison/{db}_dmnd.txt"
    group:
        "dmnd_bench"
    threads: 28
    shell:
        """
diamond blastp -q {input} -d {params.db} --out {output} --threads {resources.cpus_per_task} \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen \
--evalue 0.0001 --max-target-seqs 100 --max-hsps 1 --tmpdir $TMPDIR 2> {log}
"""


rule mmseqstax_random_seqs:
    input:
        rules.get_random_seqs.output,
    output:
        "results/comparison/tax/{db}_report",
    params:
        db=lambda wcs: str(config["benchmark"]["dbs"][wcs.db]["mmseqs"]),
        lca=config["benchmark"]["lca_mode"],
    conda:
        "../envs/benchmark.yaml"
    log:
        "results/log/comparison/{db}_mmseqs.log",
    benchmark:
        "results/benchmarks/comparison/{db}_mmseqs.txt"
    # resources: slurm_extra="'--qos=gp_bscls' '--constraint=highmem'"
    # group: "mmseqs_bench"
    threads: 112
    shell:
        """
output=$(echo {output} | sed 's/_report$//g')
mmseqs easy-taxonomy {input} {params.db} $output $TMPDIR \
--lca-mode {params.lca} --threads {resources.cpus_per_task} &> {log}
"""


rule notebook_comparison:
    input:
        mmseqs=expand(
            "results/comparison/tax/{db}_report",
            db=list(config["benchmark"]["dbs"].keys()),
        ),
        dmnd=expand(
            "results/comparison/hits/{db}_matches.tsv",
            db=list(config["benchmark"]["dbs"].keys()),
        ),
    output:
        "results/comparison/report.html",
    # conda: "../envs/R.yaml"
    localrule: True
    script:
        "../notebooks/comparison.Rmd"
