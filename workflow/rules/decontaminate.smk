rule cluster_repdb:
    input: rules.make_repdb_fasta.output.fa
    output: "results/contamination/repdb_cluster.tsv"
    params:
        identity=config["mm_identity"],
        coverage=config["mm_coverage"]
    threads: 112
    shell:'''
clusterdir=$(dirname {output})
mkdir -p $clusterdir

mmseqs easy-linclust {input} $clusterdir/repdb $TMPDIR/repdbtest \
--min-seq-id {params.identity} -c {params.coverage} --cluster-mode 2 -e 0.001 --threads {threads}
rm $clusterdir/repdb_all_seqs.fasta $clusterdir/repdb_rep_seq.fasta
'''

rule remove_singletons:
    input: 
        clusters=rules.cluster_repdb.output
    output:
        dups=temp("results/contamination/dups.ids"),
        clusters=temp("results/contamination/non_singletons_clusters.tsv")
    # threads: 24
    shell: '''
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
'''
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
        taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
    output:
        interesting="results/contamination/mixed.ids",
        mixed="results/contamination/mixed_cluster.tsv"
    shell: '''
awk '{{print $0"\\t"substr($2, 1, index($2, "_")-1)}}' {input.clusters} | \
taxonkit reformat -I 3 --data-dir {input.taxdump} -f "{{k}}"  | csvtk uniq -H -t -f 1,4 | \
sort -k4,4 | csvtk fold -t -H -f 1 -v 4 -s"|" | grep "|" | grep Eukaryota | cut -f1 > {output.interesting}

awk 'NR==FNR {{ dup[$1]; next }} $1 in dup' {output.interesting} {input.clusters} | \
awk '{{print $0"\\t"substr($2, 1, index($2, "_")-1)}}' | \
taxonkit reformat -I 3 --data-dir {input.taxdump} > {output.mixed}
'''

# should you try different mmseqs parameters?
# get eukaryotic sequences in these clusters and flag them as contaminatnts
# check how many of these are in your trees.


rule get_contaminants:
    input: rules.get_mixed_clusters.output.mixed
    output: "results/contamination/contaminants.txt"
    params:
        prop_euka=config["prop_euka"],
        size_cluster=config["size_cluster"]
    script: '../scripts/get_contaminants.R'



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
