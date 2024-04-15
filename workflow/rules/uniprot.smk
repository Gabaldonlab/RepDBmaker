rule get_uniprot_tax:
    input: 
        lineage=rules.get_uniprot_meta.output.lineage,
        unieuk=rules.download_unieuk.output,
        ep=rules.get_eukprot.output.euk_included,
        notep=rules.get_eukprot.output.euk_excluded
    output: "results/taxonomies/uniprot_eukprot_taxonomy.tsv"
    script: "../scripts/get_tax.R"

def uniprot_to_extract(wildcards):
    with open(str(checkpoints.get_euka_tax.get(**wildcards).output.tax)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("UP")]
    return expand("results/proteomes/up/{i}.faa.gz", i=genomes)

rule get_uniprot_genomes:
    input:
        up=rules.get_uniprot_meta.output.meta,
    output: "results/proteomes/up/{genome}.faa.gz"
    shell:'''
taxid=$(grep {wildcards.genome} {input.up} | cut -f2)
wget -nc "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{wildcards.genome}/{wildcards.genome}_$taxid.fasta.gz" \
-O {output}
'''

rule create_uniprot_table:
    input: uniprot_to_extract
    output: "results/meta/up_genome_table.tsv"
    shell:'''
> {output}
for genome in {input}; do
    bn=$(basename $genome ".faa.gz")
    path=$(realpath $genome)
    echo -e "$path\\t$bn\\tEukaryota"
done >> {output}
'''
