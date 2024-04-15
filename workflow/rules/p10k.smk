rule get_p10k_tax:
    input: 
        lineage=rules.get_p10k.output.lineage,
        unieuk=rules.download_unieuk.output,
        ep=rules.get_eukprot.output.euk_included,
        notep=rules.get_eukprot.output.euk_excluded
    output: "results/taxonomies/p10k_taxonomy.tsv"
    script: "../scripts/get_tax.R"
#     shell: '''
# cat {input.p10k} | Rscript src/get_p10k_tax.R -u {input.unieuk} --ep {input.ep} --notep {input.notep} > {output}
# '''

def p10k_to_extract(wildcards):
    with open(str(checkpoints.get_euka_tax.get(**wildcards).output.tax)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("P10K")]
    return expand("results/proteomes/p10k/{i}.faa.gz", i=genomes)


rule get_p10k_genomes:
    output: "results/proteomes/p10k/{genome}.faa.gz"
    shell:'''
wget -O - https://ngdc.cncb.ac.cn/p10k/static/Protein/{wildcards.genome}_protein.fa | gzip > {output}
'''

rule create_p10k_table:
    input: p10k_to_extract
    output: "results/meta/p10k_genome_table.tsv"
    shell:'''
> {output}
for genome in {input}; do
    bn=$(basename $genome ".faa.gz")
    path=$(realpath $genome)
    echo -e "$path\\t$bn\\tEukaryota"
done >> {output}
'''
