import pandas as pd

checkpoint repdb_classes:
    input: "results/taxonomies/repdb.tsv"
    output: "results/taxonomies/repdb_classes.txt"
    localrule: True
    shell: "cut -f4 {input} | sort -u | grep . > {output}"

# Read taxonomic classes from GTDB taxonomy file
def repr_all_classes(wildcards):
  with checkpoints.repdb_classes.get(**wildcards).output[0].open() as f:
    classes = [line.strip() for line in f if line.strip()]
    return expand("results/clustered_db/tmp/{rank}/{rank}_rep_seq.fasta", rank=classes)

def cluster_all_classes(wildcards):
  with checkpoints.repdb_classes.get(**wildcards).output[0].open() as f:
    classes = [line.strip() for line in f if line.strip()]
    return expand("results/clustered_db/tmp/{rank}/{rank}_cluster.tsv", rank=classes)


rule get_class:
    input:
        mmseqs=rules.make_mmseqsdb.output,
        tax=rules.create_repdb_taxdump.output.all_taxa,
        taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
    output: "results/clustered_db/tmp/{class}/{class}.fasta"
    threads: 8
    localrule: True
    conda: "../envs/homology.yaml"
    # group: "cluster_db"
    shell: """
mkdir -p results/clustered_db/tmp/{wildcards.class}

if [[ {wildcards.class} == "unclassified_"*  ]]; then
    echo "{wildcards.class} is unclassified, proceeding with lookup mode"
    awk '$4=="{wildcards.class}"' {input.tax} | cut -f1 | grep -f - {input.mmseqs}.lookup | cut -f1 > {output}.lookup
    
    mmseqs createsubdb --subdb-mode 1 --id-mode 0 -v 3 {output}.lookup {input.mmseqs} {output}_db
    mmseqs convert2fasta {output}_db {output}
    rm {output}.lookup
else 
    echo "{wildcards.class} is ok, proceeding with filtertaxseqdb mode"
    taxid=$(echo {wildcards.class} | taxonkit name2taxid --data-dir {input.taxdump} -r | awk '$3=="class"' | cut -f2)
    
    mmseqs filtertaxseqdb {input.mmseqs} {output}_db --taxon-list $taxid --threads {threads}
    mmseqs convert2fasta {output}_db {output}
fi

rm {output}_db*
"""


rule cluster_class:
    input: rules.get_class.output
    output: 
        seqs="results/clustered_db/tmp/{class}/{class}_rep_seq.fasta",
        clusters="results/clustered_db/tmp/{class}/{class}_cluster.tsv"
    params:
        identity=config["clustering"]["identity"],
        coverage=config["clustering"]["coverage"]
    conda: "../envs/homology.yaml"
    localrule: True
    # group: "cluster_db"
    threads: 2
    shell: """
clusterdir=$(dirname {output.seqs})

mmseqs easy-linclust {input} $clusterdir/{wildcards.class} $TMPDIR \
--min-seq-id {params.identity} -c {params.coverage} --cluster-mode 2 -e 0.001 --threads {threads}
rm $clusterdir/{wildcards.class}_all_seqs.fasta
"""


rule merge_clustered:
    input: 
        seqs=repr_all_classes,
        clusters=cluster_all_classes
    output: 
        seqs="results/clustered_db/repdb.fa.gz",
        clusters="results/clustered_db/db_clusters.tsv"
    # localrule: True
    shell: """
cat {input.seqs} | gzip > {output.seqs}
cat {input.clusters} > {output.clusters}
"""
