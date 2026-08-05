def _get_decon_settings(db):
    build_conf = config.get("dbs", {}).get("build", {})
    if db == "repdb":
        db_conf = build_conf.get("repdb", {})
    else:
        db_conf = build_conf.get("custom", {}).get(db, {})
    if not isinstance(db_conf, dict):
        db_conf = {}
    decontaminate_conf = db_conf.get("decontaminate")
    if not isinstance(decontaminate_conf, dict):
        decontaminate_conf = {}
    return decontaminate_conf


def _is_db_decon(db_name):
    dbs_config = config.get("dbs", {}).get("build", {})
    if db_name == "repdb":
        db_conf = dbs_config.get("repdb")
    else:
        db_conf = dbs_config.get("custom", {}).get(db_name)
    if isinstance(db_conf, dict):
        return bool(db_conf.get("decontaminate"))
    return False


# def get_decon_fa(wildcards):
#     if _is_db_decon(wildcards.db):
#         if wildcards.db == "repdb":
#             return "results/repdb/repdb.fa.gz"
#         return f"results/{wildcards.db}/decontaminate/repdb_{wildcards.db}_decon.fa.gz"


rule write_decontamination_params:
    output:
        params="results/dbs/{db}/decontaminate/decontaminate_params.yaml",
    localrule: True
    run:
        decon_conf = _get_decon_settings(wildcards.db)
        if not isinstance(decon_conf, dict):
            decon_conf = {}
        import os

        os.makedirs(os.path.dirname(output.params), exist_ok=True)
        with open(output.params, "w") as f:
            f.write("decontaminate:\n")
            for k, v in decon_conf.items():
                f.write(f"  {k}: {v}\n")


rule concat_fasta_decont:
    input:
        repdb="results/dbs/repdb/repdb.fa.gz",
        other="results/dbs/{db}/{db}.fa.gz",
    output:
        # NOT temp(): reclaimed by `cleanup` (which removes _decon.fa.gz), not by
        # Snakemake. contaminants.tsv/pair_counts.tsv sit downstream of this, so
        # temp() here dragged the whole decontamination chain on re-invocation.
        "results/dbs/{db}/decontaminate/{db}_decon.fa.gz",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
if [ "{wildcards.db}" == "repdb" ]; then
    # hard link, NOT a symlink: Snakemake treats symlink outputs as perpetually
    # out of date (they resolve to the same file/mtime as the input), which made
    # concat_fasta_decont rerun on every invocation and cascade downstream.
    # A hard link is a regular file with a stable mtime (as in decontaminate_db).
    ln -f {input.repdb} {output} 2>/dev/null || cp {input.repdb} {output}
else
    cat {input.repdb} {input.other} > {output}
fi
"""


rule cluster_decontaminate:
    input:
        rules.concat_fasta_decont.output,
    output:
        "results/dbs/{db}/decontaminate/{db}_cluster.tsv",
    params:
        identity=lambda wildcards: _get_decon_settings(wildcards.db).get(
            "identity", 0.9
        ),
        coverage=lambda wildcards: _get_decon_settings(wildcards.db).get(
            "coverage", 0.5
        ),
        cov_mode=lambda wildcards: _get_decon_settings(wildcards.db).get("cov_mode", 3),
    threads: 112
    conda:
        "../envs/homology.yaml"
    benchmark:
        "results/benchmarks/decontamination/{db}_cluster.txt"
    group:
        "decontaminate_clust"
    shell:
        """
clusterdir=$(dirname {output})
mkdir -p $clusterdir

mmseqs easy-linclust {input} $clusterdir/{wildcards.db} $TMPDIR/{wildcards.db}test \
--min-seq-id {params.identity} -c {params.coverage} --cov-mode {params.cov_mode} \
-e 0.001 --threads {threads}

rm $clusterdir/{wildcards.db}_all_seqs.fasta $clusterdir/{wildcards.db}_rep_seq.fasta
"""


rule remove_singletons:
    input:
        clusters=rules.cluster_decontaminate.output,
    output:
        dups=temp("results/dbs/{db}/decontaminate/dups.ids"),
        clusters="results/dbs/{db}/decontaminate/non_singletons_clusters.tsv",
    # threads: 24
    conda:
        "../envs/utils.yaml"
    group:
        "decontaminate"
    shell:
        """
cont_dir=$(dirname {output.dups})
echo "getting the non singletons representative"
cut -f1 {input.clusters} | uniq -d > {output.dups}

echo "splitting the file"
split {input.clusters} -n l/24 ${{cont_dir}}/chunk_

mkdir -p ${{cont_dir}}/results/

> {output.clusters}
for file in ${{cont_dir}}/chunk_*; do
    echo "processing chunk $file"
    csvtk join -H -t -f 1 "$file" {output.dups} >> {output.clusters}
    echo "done!"
done

rm ${{cont_dir}}/chunk_*
"""


# extract_lines_with_duplicates() {{
#     local chunk="$1"
#     csvtk join -H -t -f 1 "$chunk" {output.dups} > "test/mmseqs/results/$(basename $chunk)_result.txt"
# }}

# export -f extract_lines_with_duplicates

# echo "grepping the duplicates ids"
# # Use parallel to process each chunk
# parallel -j 24 extract_lines_with_duplicates ::: test/mmseqs/chunk_*

# cat test/mmseqs/results/*_result.txt > {output.clusters}

# rm -r test/mmseqs/results/


rule get_mixed_clusters:
    input:
        clusters=rules.remove_singletons.output.clusters,
        taxdump=rules.create_taxdump.output.full_taxdump,
    output:
        interesting="results/dbs/{db}/decontaminate/mixed.ids",
        mixed="results/dbs/{db}/decontaminate/mixed_cluster.tsv",
    conda:
        "../envs/utils.yaml"
    group:
        "decontaminate"
    shell:
        """
awk '{{print $0"\\t"substr($2, 1, index($2, "_")-1)}}' {input.clusters} | \
taxonkit reformat -I 3 --data-dir {input.taxdump} -f "{{k}}"  | csvtk uniq -H -t -f 1,4 | \
sort -k4,4 | csvtk fold -t -H -f 1 -v 4 -s"|" | grep "|" | grep Eukaryota | cut -f1 > {output.interesting}

awk 'NR==FNR {{ dup[$1]; next }} $1 in dup' {output.interesting} {input.clusters} | \
awk '{{print $0"\\t"substr($2, 1, index($2, "_")-1)}}' | \
taxonkit reformat -I 3 --data-dir {input.taxdump} > {output.mixed}
"""


# rule get_viral_clusters:
#     input:
#         clusters=rules.remove_singletons.output.clusters,
#         taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
#     output:
#         interesting="results/contamination/viral.ids",
#         mixed="results/contamination/viral_cluster.tsv"
#     conda: "../envs/utils.yaml"
#     group: "decontaminate"
#     shell: '''
# awk '{{print $0"\\t"substr($2, 1, index($2, "_")-1)}}' {input.clusters} | \
# taxonkit reformat -I 3 --data-dir {input.taxdump} -f "{{k}}"  | csvtk uniq -H -t -f 1,4 | \
# sort -k4,4 | csvtk fold -t -H -f 1 -v 4 -s"|" | grep Viruses | cut -f1 > {output.interesting}

# awk 'NR==FNR {{ dup[$1]; next }} $1 in dup' {output.interesting} {input.clusters} | \
# awk '{{print $0"\\t"substr($2, 1, index($2, "_")-1)}}' | \
# taxonkit reformat -I 3 --data-dir {input.taxdump} > {output.mixed}
# '''


rule get_pairwise_combination:
    input:
        rules.remove_singletons.output.clusters,
    output:
        "results/dbs/{db}/decontaminate/pair_counts.tsv",
    localrule: True
    # group: "decontaminate"
    conda:
        "../envs/utils.yaml"
    shell:
        """
cut -f2,4 -d'_' {input} | sed 's/_/\\t/g' | sort | uniq -c | sed -E 's/^[[:space:]]*([0-9]+)[[:space:]]+(.*)/\\1\t\\2/' > {output}
"""


# should you try different mmseqs parameters?
# get eukaryotic sequences in these clusters and flag them as contaminatnts
# check how many of these are in your trees.


rule get_contaminants:
    input:
        rules.get_mixed_clusters.output.mixed,
    output:
        # a data frame of the flagged proteins with their cluster properties
        # (produced regardless of the hard/soft filter mode) and the plain ID
        # list used by the filter step
        df="results/dbs/{db}/decontaminate/contaminants.tsv",
        ids="results/dbs/{db}/decontaminate/contaminants.txt",
    params:
        prop_euka=lambda wildcards: _get_decon_settings(wildcards.db).get(
            "prop_euka", 0.5
        ),
    conda:
        "../envs/R.yaml"
    group:
        "decontaminate"
    script:
        "../scripts/get_contaminants.R"


rule decontaminate_db:
    """Apply the contamination filter to the assembled DB FASTA.

    filter: hard  -> remove the flagged (contaminant) sequences from the DB
    filter: soft  -> keep them in the DB (they remain flagged in contaminants.tsv)
    In both modes the contaminants data frame (rule get_contaminants) is produced.
    """
    input:
        # decontamination is a filter on the RAW assembled fasta (detection ran
        # on the raw set too); clustering, if enabled, consumes this output.
        fa=lambda w: f"results/dbs/{w.db}/{w.db}.fa.gz",
        contaminants=rules.get_contaminants.output.ids,
    output:
        "results/dbs/{db}/{db}_decontaminated.fa.gz",
    params:
        mode=lambda wildcards: _get_decon_settings(wildcards.db).get("filter", "soft"),
    conda:
        "../envs/utils.yaml"
    group:
        "decontaminate"
    shell:
        """
if [ "{params.mode}" = "hard" ]; then
    echo "hard filter: removing $(wc -l < {input.contaminants}) contaminant sequences"
    seqkit grep -v -f {input.contaminants} {input.fa} -o {output}
elif [ "{params.mode}" = "soft" ]; then
    echo "soft filter: keeping all sequences (contaminants flagged in contaminants.tsv)"
    # no sequences removed -> hard-link instead of copying (saves a full-size
    # duplicate of the raw fasta). Falls back to cp across filesystems. A hard
    # link survives later deletion of the raw file (unlike a symlink).
    ln {input.fa} {output} 2>/dev/null || cp {input.fa} {output}
else
    echo "ERROR: invalid decontaminate.filter '{params.mode}' (expected 'soft' or 'hard')" >&2
    exit 1
fi
"""


# rule diamond_noneuk:
#     input:
#         query="results/proteomes/cus/CUS00001.faa.gz",
#         db=rules.make_diamonddb.output,
#         taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
#     output: "test/decon/test.out"
#     threads: 12
#     shell: '''
# diamond blastp -q {input.query} -d {input.db} --out {output} --outfmt 6 \
# --threads {threads} --fast \
# --evalue 1e-10 --max-target-seqs 5 --id 80 --query-cover 70
# '''
# euka_taxid=$(echo Eukaryota | taxonkit name2taxid --data-dir {input.taxdump} | cut -f2)
# --taxon-exclude $euka_taxid
# the idea is to blast non-uniprot proteomes to the first iteration of repdb. If a sequence has a very \
# close homologs to non euka it could be excluded
