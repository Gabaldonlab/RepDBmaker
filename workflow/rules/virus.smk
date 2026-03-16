checkpoint check_virus_genomes:
    input:
        meta=rules.get_virus_genomes.output.meta,
        folder=rules.get_virus_genomes.output.folder
    output:
        ids="results/meta/refseq_virus_ids.txt",
        tax="results/taxonomies/virus_taxonomy.tsv"
    params:
        taxdump=rules.download_taxdump.output
    conda: "../envs/utils.yaml"
    localrule: True
    shell: '''
unzip -l {input.folder} | awk '{{print $NF}}' | grep protein | cut -f3 -d'/' > {output.ids}
cut -f1,6 {input.meta} | grep -F -w -f {output.ids} | sed 's/_//' | sed 's/\\..*\\t/\\t/g' | \
taxonkit reformat -I 2 -P -F -p "unclassified_" -s "" --data-dir {params.taxdump} | cut -f1,3 | \
sed 's/k__unclassified_Viruses/d__Viruses/g'> {output.tax}
'''

# because of this https://ncbiinsights.ncbi.nlm.nih.gov/2025/02/27/new-ranks-ncbi-taxonomy/
# now I have to add viruses tag to domain
# probably would be a good idea to have the new realm as domain but it seems a bit out of scope for this

def viruses_to_extract(wildcards):
    with open(str(checkpoints.check_virus_genomes.get(**wildcards).output.ids)) as virus:
        genomes = [gn.strip().replace('_','').split('.')[0] for gn in virus]
        return expand("results/proteomes/virus/{i}.faa.gz", i=genomes)

rule extract_virus:
    input:
        ids=rules.check_virus_genomes.output.ids,
        folder=rules.get_virus_genomes.output.folder
    output: "results/proteomes/virus/{virus}.faa.gz",
    localrule: True
    shell:'''
virus_id=$(grep $(echo {wildcards.virus} | sed 's/^GCF//') {input.ids})
unzip -p {input.folder} ncbi_dataset/data/$virus_id/protein.faa | gzip > {output}
'''


rule create_virus_table:
    input: viruses_to_extract
    output: "results/meta/virus_genome_table.tsv"
    localrule: True
    shell:'''
> {output}
for genome in {input}; do
    bn=$(basename $genome ".faa.gz")
    path=$(realpath $genome)
    echo -e "$path\\t$bn\\tViruses"
done >> {output}
'''
