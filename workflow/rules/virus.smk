checkpoint check_virus_genomes:
    input:
        meta=rules.get_virus_genomes.output.meta,
        folder=rules.get_virus_genomes.output.folder
    output:
        ids="results/meta/refseq_virus_ids.txt",
        tax="results/taxonomies/virus_taxonomy.tsv"
    params:
        taxdump=rules.download_taxdump.output
    shell: '''
unzip -l {input.folder} | awk '{{print $NF}}' | grep protein | cut -f3 -d'/' > {output.ids}
cut -f1,6 {input.meta} | grep -F -w -f {output.ids} | sed 's/_//' | sed 's/\\..*\\t/\\t/g' | \
taxonkit reformat -I 2 -P --data-dir {params.taxdump} | cut -f1,3 | sed 's/k__/d__/g'> {output.tax}
'''

def viruses_to_extract(wildcards):
    with open(str(checkpoints.check_virus_genomes.get(**wildcards).output.ids)) as virus:
        genomes = [gn.strip().replace('_','').split('.')[0] for gn in virus]
        return expand("results/proteomes/virus/{i}.faa.gz", i=genomes)

rule extract_virus:
    input:
        ids=rules.check_virus_genomes.output.ids,
        folder=rules.get_virus_genomes.output.folder
    output: "results/proteomes/virus/{virus}.faa.gz",
    shell:'''
virus_id=$(grep $(echo {wildcards.virus} | sed 's/^GCF//') {input.ids})
unzip -p {input.folder} ncbi_dataset/data/$virus_id/protein.faa | gzip > {output}
'''


rule create_virus_table:
    input: viruses_to_extract
    output: "results/meta/virus_genome_table.tsv"
    shell:'''
> {output}
for genome in {input}; do
    bn=$(basename $genome ".faa.gz")
    path=$(realpath $genome)
    echo -e "$path\\t$bn\\tViruses"
done >> {output}
'''
