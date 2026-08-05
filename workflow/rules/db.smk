import re

CLUSTER_SUFFIX = "_clustered"


def _db_conf(db_name):
    """The build config dict for a REAL db (repdb or a custom entry), or None."""
    dbs_config = config.get("dbs", {}).get("build", {})
    if db_name == "repdb":
        return dbs_config.get("repdb")
    return dbs_config.get("custom", {}).get(db_name)


def _has_cluster_variant(db_name):
    """True for a REAL db that declares a `cluster` block. Clustering is now
    ADDITIVE: such a db spawns a parallel, first-class sibling db
    `<db>_clustered` (built by clustering this db's final fasta), instead of the
    clustering replacing this db's own final in place."""
    conf = _db_conf(db_name)
    return isinstance(conf, dict) and bool(conf.get("cluster"))


def _cluster_parent(db_name):
    """If `db_name` is a synthetic clustered variant (`<parent>_clustered` whose
    parent declares a `cluster` block), return the parent; else None."""
    if db_name.endswith(CLUSTER_SUFFIX):
        parent = db_name[: -len(CLUSTER_SUFFIX)]
        if _has_cluster_variant(parent):
            return parent
    return None


def _cluster_keep_full(db_name):
    """Whether the full (unclustered) db is kept by `cleanup` when it also has a
    clustered variant. `cluster.keep_full: false` -> the full db is treated as a
    disposable intermediate (its indices are not even built)."""
    conf = _db_conf(db_name)
    cl = conf.get("cluster") if isinstance(conf, dict) else None
    return bool(cl.get("keep_full", True)) if isinstance(cl, dict) else True


def _real_db_names():
    dbs = config.get("dbs", {}).get("build", {})
    names = ["repdb"] if "repdb" in dbs else []
    custom = dbs.get("custom", {})
    if isinstance(custom, dict):
        names.extend(custom.keys())
    return names


REAL_DBS = _real_db_names()
# a real db may not collide with the reserved `<parent>_clustered` namespace
for _d in REAL_DBS:
    if _d.endswith(CLUSTER_SUFFIX):
        raise ValueError(
            f"db name '{_d}' ends in '{CLUSTER_SUFFIX}', which is reserved for the "
            f"clustered variant of another db. Please rename it."
        )
PARENT_DBS = [d for d in REAL_DBS if _has_cluster_variant(d)]
CLUSTERED_DBS = [f"{d}{CLUSTER_SUFFIX}" for d in PARENT_DBS]
ALL_DBS = REAL_DBS + CLUSTERED_DBS

# enumerated alternations (unambiguous) for wildcard fencing. `$.^` never matches.
_C_REAL = "|".join(re.escape(d) for d in REAL_DBS) or r"$.^"
_C_PARENTS = "|".join(re.escape(d) for d in PARENT_DBS) or r"$.^"
_C_ALL = "|".join(re.escape(d) for d in ALL_DBS) or r"$.^"


wildcard_constraints:
    db=_C_ALL,
    dbtype="|".join(sorted(valid_db_types)),


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
    """FASTA the search indices are built from.

    Real db: decontaminated if it has a `decontaminate` block, else the raw
    assembled FASTA. Clustering no longer replaces the final in place - it
    produces a separate `<db>_clustered` sibling db.

    Synthetic clustered variant (`<parent>_clustered`): its own published fasta
    at results/dbs/<db>/<db>.fa.gz (written by publish_clustered_variant); since
    such a db has no decontaminate block, `_fa_raw` already points there."""
    if _is_db_decon(db):
        return _fa_decontaminated(db)
    return _fa_raw(db)


def get_final_db_fa(wildcards):
    """Input-function form of `_final_db_fa` (`_is_db_decon` lives in
    decontaminate.smk; resolved at DAG time)."""
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
    # real dbs only; a `<db>_clustered` variant gets its fasta from clustering
    # (publish_clustered_variant), not from a genome table.
    wildcard_constraints:
        db=_C_REAL,
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
    # real dbs only; the clustered variant reuses its parent's maps (the
    # representative IDs are a subset), hard-linked by publish_clustered_variant.
    wildcard_constraints:
        db=_C_REAL,
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
    # keyed on the PARENT db (the one with a cluster block); consumes the
    # parent's final fasta (decontaminated or raw) as the clustering substrate.
    wildcard_constraints:
        db=_C_PARENTS,
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
