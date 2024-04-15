## Custom new genomes

rule get_custom_tax:
    input: config["new_genomes"]
    output: "results/taxonomies/custom_taxonomy.tsv"
    shell:'''
cut -f1,4 {input} | awk 'NR>1' > {output}
'''

def custom_to_extract(wildcards):
    with open(str(checkpoints.get_euka_tax.get(**wildcards).output.tax)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("CUS")]
    return expand("results/proteomes/cus/{i}.faa.gz", i=genomes)


rule get_custom_genomes:
    input:
        custom_table=config["new_genomes"]
    output: "results/proteomes/cus/{genome}.faa.gz"
    shell:'''
file=$(grep {wildcards.genome} {input.custom_table} | cut -f3)
cat $file | gzip > {output} 
'''

rule create_custom_table:
    input: custom_to_extract
    output: "results/meta/cus_genome_table.tsv"
    shell:'''
> {output}
for genome in {input}; do
    bn=$(basename $genome ".faa.gz")
    path=$(realpath $genome)
    echo -e "$path\\t$bn\\tEukaryota"
done >> {output}
'''
