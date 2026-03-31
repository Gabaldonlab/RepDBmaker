rule virus_taxonomy:
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

rule p10k_taxonomy:
    input: 
        lineage=rules.get_p10k.output.lineage,
        unieuk=rules.download_unieuk.output,
        ep=rules.get_eukprot.output.euk_included,
        notep=rules.get_eukprot.output.euk_excluded
    output: "results/taxonomies/p10k_taxonomy.tsv"
    conda: "../envs/R.yaml"
    localrule: True
    script: "../scripts/get_tax.R"


rule custom_taxonomy:
    input: config["files"]["new_genomes"]
    output: "results/taxonomies/custom_taxonomy.tsv"
    conda: "../envs/utils.yaml"
    localrule: True
    shell:'''
csvtk cut -t -f ID,Lineage {input} | awk 'NR>1' > {output}
'''


rule uniprot_taxonomy:
    input: 
        lineage=rules.get_uniprot_meta.output.lineage,
        unieuk=rules.download_unieuk.output,
        ep=rules.get_eukprot.output.euk_included,
        notep=rules.get_eukprot.output.euk_excluded
    output: "results/taxonomies/uniprot_eukprot_taxonomy.tsv"
    conda: "../envs/R.yaml"
    localrule: True
    script: "../scripts/get_tax.R"


rule eukprot_taxonomy:
    input: rules.get_eukprot.output.euk_included
    output: "results/taxonomies/eukprot_taxonomy.tsv"
    localrule: True
    shell: '''
awk 'NR>1' {input} | sed 's/\\_/ /g' | \
awk -F'\\t' '{{split($11, a, ";"); print $1,"d__"a[1]";p__"a[2]";c__"$8";o__"$9";f__"$10";g__"$6";s__"$2}}' OFS='\\t' | \
sed  's/N\\/A//g' | tr -d \\'\\" | sed 's/\\ ;/;/g' | sed 's/other Gyrista/other_Gyrista/g'> {output}
'''


rule eukaryotes_taxonomy:
    input:
        up=rules.uniprot_taxonomy.output,
        ep=rules.eukprot_taxonomy.output,
        p10k=rules.p10k_taxonomy.output,
        custom=rules.custom_taxonomy.output
    output: "results/taxonomies/eukaryotes_taxonomy.tsv"
    localrule: True
    shell: "cat {input} | sort -k2,2 > {output}"


rule select_repdb_eukaryotes:
    input:
        up=rules.uniprot_taxonomy.output,
        ep=rules.eukprot_taxonomy.output,
        p10k=rules.p10k_taxonomy.output,
        custom=rules.custom_taxonomy.output,
        up_stats=rules.get_uniprot_meta.output.stats,
        ep_stats=rules.get_eukprot.output.euk_busco,
        p10k_stats=rules.get_p10k.output.meta,
        exclude=config["files"]["genomes_to_exclude"],
        to_keep=config["files"]["clades_to_keep"]
    output:
        tax="results/taxonomies/selected_eukaryotes.tsv",
        plot="results/plots/euka_db.pdf"
    conda: "../envs/R.yaml"
    localrule: True 
    script:"../scripts/filter_tax.R"


rule create_full_taxdump:
    """Create taxdump from all genomes (used by all databases)"""
    input:
        # gtdb_proteomes=rules.decompress_gtdb_genomes.output,
        virus=rules.virus_taxonomy.output.tax,
        gtdb=rules.get_gtdb_tax.output.tax,
        all_euka=rules.eukaryotes_taxonomy.output
    output:
        # ids="results/meta/all.ids",
        all_taxa="results/taxonomies/all_taxonomy.tsv",
        full_taxdump=directory("results/taxdump/repdb_taxdump/")
    conda: "../envs/utils.yaml"
    localrule: True
    shell:'''
cat {input.virus} {input.gtdb} {input.all_euka} | \
sed 's/d__//g' | sed 's/[a-z]__/\\t/g' | sed 's/;//g' > {output.all_taxa}
taxonkit create-taxdump -A1 {output.all_taxa} --out-dir {output.full_taxdump} --force \
--rank-names "superkingdom","phylum","class","order","family","genus","species"
'''
# find $(dirname {input.gtdb_proteomes}) -type f -name "*faa.gz" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort > {output.ids}


checkpoint repdb_taxonomy:
    """Create RepDB-specific taxonomy with selected eukaryotes, all viruses and all GTDB genomes"""
    input:
        gtdb=rules.gtdb_species_clusters.output,
        virus=rules.virus_taxonomy.output.tax,
        filtered_euka=rules.select_repdb_eukaryotes.output.tax,
        full_table=rules.create_full_taxdump.output.all_taxa
    output: "results/taxonomies/repdb_taxonomy.tsv"
    localrule: True
    shell:'''
cat {input.gtdb} {input.virus} {input.filtered_euka} | cut -f1 | \
csvtk join -H -t - {input.full_table} > {output}
'''
# cat <(shuf -n 50 {input.gtdb}) <(shuf -n 50 {input.virus}) <(shuf -n 50 {input.filtered_euka}) | cut -f1 | \

# cat {input.gtdb} {input.virus} $(shuf -n 20 {input.filtered_euka}) | cut -f1 | \

def get_custom_codes(wildcards):
    """Return the path to the genome table for a custom database"""
    custom_dbs = config.get("dbs", {}).get("build", {}).get("custom", {})
    # print(f"Custom DBs in config: {custom_dbs}")
    if isinstance(custom_dbs, dict) and wildcards.db in custom_dbs:
        return custom_dbs[wildcards.db]["ids"]
    raise ValueError(f"Custom database '{wildcards.db}' not found in config")


checkpoint custom_dbs_taxonomy:
    """Create user-defined taxonomy for custom database"""
    input: 
        mnemonics=get_custom_codes,
        full_table=rules.create_full_taxdump.output.all_taxa
    output: "results/taxonomies/{db}_taxonomy.tsv"
    localrule: True
    shell:'''
csvtk join -H -t {input.mnemonics} {input.full_table} > {output}
'''

