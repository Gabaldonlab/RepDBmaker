import pandas as pd

dict_ranks = {
    "kingdom": 2,
    "phylum": 3,
    "class": 4,
    "order": 5,
    "family": 6,
    "genus": 7,
    "species": 8,
}


def _get_cluster_settings(db):
    build_conf = config.get("dbs", {}).get("build", {})
    if db == "repdb":
        db_conf = build_conf.get("repdb", {})
    else:
        db_conf = build_conf.get("custom", {}).get(db, {})
    if not isinstance(db_conf, dict):
        db_conf = {}
    clustered_conf = db_conf.get("cluster")
    if not isinstance(clustered_conf, dict):
        clustered_conf = {}
    return clustered_conf


def _validate_cluster_levels():
    if not isinstance(config.get("dbs", {}).get("build", {}), dict):
        return
    build_conf = config["dbs"]["build"]
    # check repdb cluster config
    repdb_conf = build_conf.get("repdb")
    if isinstance(repdb_conf, dict):
        cluster_conf = (
            repdb_conf.get("cluster")
            or repdb_conf.get("clustered")
            or repdb_conf.get("clustering")
        )
        if isinstance(cluster_conf, dict):
            level = cluster_conf.get("level")
            if level and level not in dict_ranks:
                raise ValueError(
                    f"Invalid cluster level for repdb: '{level}'. Allowed: {sorted(dict_ranks.keys())}"
                )
    # check custom dbs
    custom_conf = build_conf.get("custom", {})
    if isinstance(custom_conf, dict):
        for db_name, conf in custom_conf.items():
            if isinstance(conf, dict):
                cluster_conf = (
                    conf.get("cluster")
                    or conf.get("clustered")
                    or conf.get("clustering")
                )
                if isinstance(cluster_conf, dict):
                    level = cluster_conf.get("level")
                    if level and level not in dict_ranks:
                        raise ValueError(
                            f"Invalid cluster level for custom db '{db_name}': '{level}'. Allowed: {sorted(dict_ranks.keys())}"
                        )


_validate_cluster_levels()


rule write_cluster_params:
    output:
        params="results/dbs/{db}/cluster/cluster_params.yaml",
    localrule: True
    run:
        cluster_conf = _get_cluster_settings(wildcards.db)
        if not isinstance(cluster_conf, dict):
            cluster_conf = {}
        import os

        os.makedirs(os.path.dirname(output.params), exist_ok=True)
        with open(output.params, "w") as f:
            f.write("cluster:\n")
            for k, v in cluster_conf.items():
                f.write(f"  {k}: {v}\n")


checkpoint db_clades:
    input:
        "results/taxonomies/{db}_taxonomy.tsv",
    output:
        "results/dbs/{db}/cluster/{db}_clades.txt",
    params:
        rank=lambda wildcards: dict_ranks.get(
            _get_cluster_settings(wildcards.db).get("level", "class")
        ),
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        "cut -f{params.rank} {input} | sort -u | grep . > {output}"


# Read taxonomic clades from GTDB taxonomy file
def repr_all_clades(wildcards):
    with checkpoints.db_clades.get(**wildcards).output[0].open() as f:
        clades = [line.strip() for line in f if line.strip()]
        return expand(
            "results/dbs/{db}/cluster/tmp/{rank}/{rank}_rep_seq.fasta",
            db=wildcards.db,
            rank=clades,
        )


def cluster_all_clades(wildcards):
    with checkpoints.db_clades.get(**wildcards).output[0].open() as f:
        clades = [line.strip() for line in f if line.strip()]
        return expand(
            "results/dbs/{db}/cluster/tmp/{rank}/{rank}_cluster.tsv",
            db=wildcards.db,
            rank=clades,
        )


rule get_clade:
    input:
        mmseqs=rules.make_mmseqsdb_clustering.output.db,
        mmseqs_extra=rules.make_mmseqsdb_clustering.output.db_extra,
        # universe → id + 7 ranks (cut -f1,3-9), same column layout the awk below
        # expects (col2=superkingdom ... col8=species).
        tax="results/universe/universe.tsv",
        taxdump=rules.create_taxdump.output.full_taxdump,
    output:
        temp("results/dbs/{db}/cluster/tmp/{clade}/{clade}.fasta"),
    params:
        level=lambda wildcards: _get_cluster_settings(wildcards.db).get(
            "level", "class"
        ),
        rank=lambda wildcards: dict_ranks.get(
            _get_cluster_settings(wildcards.db).get("level", "class")
        ),
    threads: 8
    localrule: True
    conda:
        "../envs/homology.yaml"
    # group: "cluster_db"
    shell:
        """
mkdir -p $(dirname {output})

if [[ {wildcards.clade} == "unclassified_"*  ]]; then
    echo "{wildcards.clade} is unclassified, proceeding with lookup mode"
    awk -F'\\t' 'NR>1' {input.tax} | cut -f1,3-9 | awk '${params.rank}=="{wildcards.clade}"' | cut -f1 | grep -f - {input.mmseqs}.lookup | cut -f1 > {output}.lookup
    
    mmseqs createsubdb --subdb-mode 1 --id-mode 0 -v 3 {output}.lookup {input.mmseqs} {output}_db
    mmseqs convert2fasta {output}_db {output}
    rm {output}.lookup
else 
    echo "{wildcards.clade} is ok, proceeding with filtertaxseqdb mode"
    taxid=$(echo "{wildcards.clade}" | taxonkit name2taxid --data-dir {input.taxdump} -r | awk '$3=="{params.level}"' | cut -f2)
    
    mmseqs filtertaxseqdb {input.mmseqs} {output}_db --taxon-list $taxid --threads {threads}
    mmseqs convert2fasta {output}_db {output}
fi

rm {output}_db*
"""


rule cluster_clade:
    input:
        rules.get_clade.output,
    output:
        seqs=temp("results/dbs/{db}/cluster/tmp/{clade}/{clade}_rep_seq.fasta"),
        clusters=temp("results/dbs/{db}/cluster/tmp/{clade}/{clade}_cluster.tsv"),
    params:
        identity=lambda wildcards: _get_cluster_settings(wildcards.db).get(
            "identity", 0.9
        ),
        coverage=lambda wildcards: _get_cluster_settings(wildcards.db).get(
            "coverage", 0.9
        ),
    conda:
        "../envs/homology.yaml"
    localrule: True
    # group: "cluster_db"
    threads: 8
    shell:
        """
clusterdir=$(dirname {output.seqs})
tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/mmseqs.XXXXXX")

mmseqs easy-linclust {input} $clusterdir/{wildcards.clade} "$tmp" \
--min-seq-id {params.identity} -c {params.coverage} --cluster-mode 2 -e 0.001 --threads {threads}

rm -rf "$tmp"
rm $clusterdir/{wildcards.clade}_all_seqs.fasta
"""


rule merge_clustered:
    input:
        seqs=repr_all_clades,
        clusters=cluster_all_clades,
    output:
        seqs="results/dbs/{db}/{db}_clustered.fa.gz",
        clusters="results/dbs/{db}/{db}_clusters.tsv",
    # localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
cat {input.seqs} | gzip > {output.seqs}
cat {input.clusters} > {output.clusters}
"""
