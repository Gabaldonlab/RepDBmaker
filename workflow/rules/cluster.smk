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
    wildcard_constraints:
        db=_C_PARENTS,
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
    wildcard_constraints:
        db=_C_PARENTS,
    input:
        # id + 7 ranks for the FULL selection ...
        tax="results/taxonomies/{db}_taxonomy.tsv",
        # ... intersected with the genomes actually built into this db. In test
        # mode genome_table.tsv is subsampled (its inputs go through
        # _test_subset), so without this intersection db_clades would list every
        # clade in the full taxonomy and the empty ones would make cluster_clade
        # fail on empty input. In a full run the two sets coincide, so this is a
        # no-op. genome_table col2 (file_code) == taxonomy col1 (id).
        table="results/dbs/{db}/genome_table.tsv",
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
        """
cut -f2 {input.table} | sort -u > {output}.codes
awk -F'\\t' 'NR==FNR{{present[$1]; next}} ($1 in present)' {output}.codes {input.tax} | \
cut -f{params.rank} | sort -u | grep . > {output}
rm -f {output}.codes
"""


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
    wildcard_constraints:
        db=_C_PARENTS,
    input:
        mmseqs=rules.make_mmseqsdb_clustering.output.db,
        mmseqs_extra=rules.make_mmseqsdb_clustering.output.db_extra,
        # universe → id + 7 ranks (cut -f1,3-9), same column layout the awk below
        # expects (col2=superkingdom ... col8=species).
        tax="results/universe/universe.tsv",
        taxdump=rules.create_taxdump.output.full_taxdump,
    output:
        # NOT temp(): reclaimed by `cleanup` (which removes cluster/tmp), not by
        # Snakemake - see make_mmseqsdb_clustering. Keeps re-invoked builds no-op.
        "results/dbs/{db}/cluster/tmp/{clade}/{clade}.fasta",
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
    wildcard_constraints:
        db=_C_PARENTS,
    input:
        rules.get_clade.output,
    output:
        # NOT temp(): reclaimed by `cleanup` (cluster/tmp), not by Snakemake.
        # clustdb_stats + merge_clustered read these directly, so temp() here is
        # what dragged the whole clustering chain on every re-invocation.
        seqs="results/dbs/{db}/cluster/tmp/{clade}/{clade}_rep_seq.fasta",
        clusters="results/dbs/{db}/cluster/tmp/{clade}/{clade}_cluster.tsv",
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
    # keyed on the PARENT db; produces that db's clustered fasta as an
    # intermediate, which publish_clustered_variant exposes as a sibling db.
    wildcard_constraints:
        db=_C_PARENTS,
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


rule publish_clustered_variant:
    wildcard_constraints:
        parent=_C_PARENTS,
    input:
        fa="results/dbs/{parent}/{parent}_clustered.fa.gz",
        headmap="results/dbs/{parent}/{parent}.map",
        noheadmap="results/dbs/{parent}/{parent}_nohead.map",
    output:
        fa="results/dbs/{parent}_clustered/{parent}_clustered.fa.gz",
        headmap="results/dbs/{parent}_clustered/{parent}_clustered.map",
        noheadmap="results/dbs/{parent}_clustered/{parent}_clustered_nohead.map",
    localrule: True
    shell:
        """
mkdir -p $(dirname {output.fa})
# relative symlinks (ln -rs), not hard links: own inode/mtime, no shared-inode
# mtime poisoning, ~0 bytes. cp fallback just in case.
ln -rsf {input.fa} {output.fa} 2>/dev/null || cp {input.fa} {output.fa}
ln -rsf {input.headmap} {output.headmap} 2>/dev/null || cp {input.headmap} {output.headmap}
ln -rsf {input.noheadmap} {output.noheadmap} 2>/dev/null || cp {input.noheadmap} {output.noheadmap}
"""
