rule online_resources:
    input:
        rules.download_taxdump.output,
        rules.download_unieuk.output,
        rules.get_virus_genomes.output,
        rules.get_gtdb_tax.output,
        rules.get_gtdb_genomes.output,
        rules.get_uniprot_meta.output,
        rules.get_eukprot.output,
        rules.get_p10k.output,
        # rules.create_p10k_table.output,
        # rules.create_uniprot_table.output,
        rules.get_virus_class.output
        # rules.get_nr_accessionmap.output
    output: "results/meta/check_resources.txt"
    localrule: True
    shell:'''
for file in {input}; do
    if [ -s "$file" ]; then
        echo "$file ok"
    else
        echo "$file NOT FOUND"
    fi
done > {output}
'''

rule available_proteomes:
    input:
        gtdb=rules.get_gtdb_tax.output.meta,
        gtdb_tax=rules.get_gtdb_tax.output.tax,
        viruses_tax=rules.virus_taxonomy.output.tax,
        viruses=rules.get_virus_genomes.output.meta,
        euka_tax=rules.eukaryotes_taxonomy.output,
        ep_included=rules.get_eukprot.output.euk_included,
        ep_busco=rules.get_eukprot.output.euk_busco,
        p10k=rules.get_p10k.output.meta,
        up=rules.get_uniprot_meta.output.stats,
    output: "results/meta/available_proteomes.tsv"
    conda: "../envs/R.yaml"
    localrule: True
    script: "../scripts/get_available_proteomes.R"

def custom_repdb(wildcards):
    with open(str(checkpoints.repdb_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("CUS")]
    return expand("results/proteomes/cus/{i}.faa.gz", i=genomes)

# eheh so good at programming yeas
def custom_custom(wildcards):
    with open(str(checkpoints.custom_dbs_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("CUS")]
    return expand("results/proteomes/cus/{i}.faa.gz", i=genomes)

def eukprot_repdb(wildcards):
    with open(str(checkpoints.repdb_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("EP")]
    return expand("results/proteomes/ep/{i}.faa.gz", i=genomes)

def eukprot_custom(wildcards):
    with open(str(checkpoints.custom_dbs_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("EP")]
    return expand("results/proteomes/ep/{i}.faa.gz", i=genomes)

def p10k_repdb(wildcards):
    with open(str(checkpoints.repdb_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("P10K")]
    return expand("results/proteomes/p10k/{i}.faa.gz", i=genomes)

def p10k_custom(wildcards):
    with open(str(checkpoints.custom_dbs_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("P10K")]
    return expand("results/proteomes/p10k/{i}.faa.gz", i=genomes)

def uniprot_repdb(wildcards):
    with open(str(checkpoints.repdb_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("UP")]
    return expand("results/proteomes/up/{i}.faa.gz", i=genomes)

def uniprot_custom(wildcards):
    with open(str(checkpoints.custom_dbs_taxonomy.get(**wildcards).output)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("UP")]
    return expand("results/proteomes/up/{i}.faa.gz", i=genomes)

def viruses_repdb(wildcards):
    with open(str(checkpoints.repdb_taxonomy.get(**wildcards).output)) as virus:
        viruses = [gn.strip().split("\t")[0] for gn in virus if gn.strip().split("\t")[1] == "Viruses"]
    return expand("results/proteomes/virus/{i}.faa.gz", i=viruses)

def viruses_custom(wildcards):
    with open(str(checkpoints.custom_dbs_taxonomy.get(**wildcards).output)) as virus:
        viruses = [gn.strip().split("\t")[0] for gn in virus if gn.strip().split("\t")[1] == "Viruses"]
    return expand("results/proteomes/virus/{i}.faa.gz", i=viruses)

def gtdb_repdb(wildcards):
    with open(str(checkpoints.repdb_taxonomy.get(**wildcards).output)) as gtdb:
        prokas = [gn.strip().split("\t")[0] for gn in gtdb if gn.strip().split("\t")[1] in ["Bacteria", "Archaea"]]
    return expand("results/proteomes/gtdb/{i}.faa.gz", i=prokas)

def gtdb_custom(wildcards):
    with open(str(checkpoints.custom_dbs_taxonomy.get(**wildcards).output)) as gtdb:
        prokas = [gn.strip().split("\t")[0] for gn in gtdb if gn.strip().split("\t")[1] in ["Bacteria", "Archaea"]]
    return expand("results/proteomes/gtdb/{i}.faa.gz", i=prokas)

# GTDB is not treated like this because it would mean a huuuuuge snakemake workflow,
# far from ideal but ok for now

rule get_custom_genomes:
    input:
        custom_table=config["files"]["new_genomes"]
    output: "results/proteomes/cus/{genome}.faa.gz"
    localrule: True
    shell:'''
file=$(grep {wildcards.genome} {input.custom_table} | cut -f4)
cat $file | gzip > {output} 
'''

rule get_eukprot_genomes:
    input:
        euk_fa=rules.get_eukprot.output.euk_fa
    output: "results/proteomes/ep/{genome}.faa.gz"
    localrule: True
    shell:'''
odir=$(dirname {output})
tar -vxzf {input.euk_fa} -C $odir --strip-components 1 --wildcards "proteins/{wildcards.genome}*" --occurrence=1 --anchored
cat $odir/{wildcards.genome}*.fasta | gzip > {output}
rm $odir/{wildcards.genome}*.fasta
'''

rule get_p10k_genomes:
    output: "results/proteomes/p10k/{genome}.faa.gz"
    localrule: True
    shell:'''
wget -O - https://ngdc.cncb.ac.cn/p10k/static/Protein/{wildcards.genome}_protein.fa | gzip > {output}
'''

rule get_uniprot_genomes:
    input:
        up=rules.get_uniprot_meta.output.meta,
    output: "results/proteomes/up/{genome}.faa.gz"
    localrule: True
    shell:'''
taxid=$(grep {wildcards.genome} {input.up} | cut -f2)
wget -nc "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{wildcards.genome}/{wildcards.genome}_$taxid.fasta.gz" \
-O {output}
'''

rule extract_virus:
    input:
        ids=rules.virus_taxonomy.output.ids,
        folder=rules.get_virus_genomes.output.folder
    output: "results/proteomes/virus/{virus}.faa.gz",
    localrule: True
    shell:'''
virus_id=$(grep $(echo {wildcards.virus} | sed 's/^GCF//') {input.ids})
unzip -p {input.folder} ncbi_dataset/data/$virus_id/protein.faa | gzip > {output}
'''

rule repdb_online_genomes:
    """touch files for those dbs that need internet access"""
    input:
        uniprot_repdb,
        p10k_repdb,
    output: "results/dbs/repdb/online_gnms.txt"
    localrule: True
    shell:'''
touch {output}
'''

rule repdb_genome_table:
    """Create comprehensive genome table with selected RepDB genomes (unfiltered)"""
    input:
        custom_repdb,
        uniprot_repdb,
        p10k_repdb,
        eukprot_repdb,
        gtdb_repdb,
        viruses_repdb
    output: "results/dbs/repdb/genome_table.tsv"
    localrule: True
    conda: "../envs/python.yaml"
    script: "../scripts/make_genome_table.py"


rule custom_genome_table:
    """Create comprehensive genome table for custom database with requested genomes (unfiltered)"""
    input:
        uniprot_custom,
        p10k_custom,
        eukprot_custom,
        custom_custom,
        gtdb_custom,
        viruses_custom,
    output: "results/dbs/{db}/genome_table.tsv"
    localrule: True
    conda: "../envs/python.yaml"
    script: "../scripts/make_genome_table.py"

