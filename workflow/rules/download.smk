# Download unieuk taxonomy

rule download_unieuk:
    output:
        "results/meta/unieuk_taxonomy.tsv"
    localrule: True
    shell:'''
wget https://eukmap.unieuk.net/exports/unieuk/1.0.0-first_release/unieuk-1.0.0-first_release.tsv -O {output}
'''

rule download_taxdump:
    output:
        td=directory("results/taxdump/ncbi_taxdump")
        # acc2taxid="results/tmp/prot.accession2taxid.FULL.gz"
    localrule: True
    shell:'''
mkdir -p {output.td}
wget -O {output.td}/ncbi_taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xf {output.td}/ncbi_taxdump.tar.gz -C {output.td}
'''
# wget -O {output.acc2taxid} https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz


rule get_gtdb_tax:
    output:
        bac="results/tmp/bac.tmp",
        ar="results/tmp/ar.tmp",
        tax="results/taxonomies/gtdb_taxonomy.tsv",
        bac_meta="results/tmp/bac_meta.tmp",
        ar_meta="results/tmp/ar_meta.tmp",
        meta="results/meta/gtdb_meta.tsv"
    conda: "../envs/utils.yaml"
    params: gtdb=config["dbs"]["gtdb_version"]
    localrule: True
    shell:'''
wget -O {output.bac} https://data.ace.uq.edu.au/public/gtdb/data/releases/{params.gtdb}/bac120_taxonomy.tsv
wget -O {output.ar}  https://data.ace.uq.edu.au/public/gtdb/data/releases/{params.gtdb}/ar53_taxonomy.tsv
wget -O {output.bac_meta} https://data.ace.uq.edu.au/public/gtdb/data/releases/{params.gtdb}/bac120_metadata.tsv.gz
wget -O {output.ar_meta} https://data.ace.uq.edu.au/public/gtdb/data/releases/{params.gtdb}/ar53_metadata.tsv.gz

zcat {output.bac_meta} {output.ar_meta} | csvtk filter2 -t -f'$gtdb_representative=="t"' > {output.meta}

cat {output.bac} {output.ar} | awk 'BEGIN{{OFS=FS="\\t"}} {{ $1 = substr($1, 4, 13) }} 1' | \
sed 's/_//' > {output.tax}
'''
# | src/filter_gtdb.R

rule get_gtdb_genomes:
    output: "results/tmp/gtdb_proteins_aa_reps.tar.gz"
    params: gtdb=config["dbs"]["gtdb_version"]
    shell:'''
wget -O {output} https://data.ace.uq.edu.au/public/gtdb/data/releases/{params.gtdb}/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz
'''

# Get refseq viruses
rule get_virus_genomes:
    output:
        meta="results/meta/refseq_virus_meta.txt",
        folder="results/tmp/refseq.zip",
    conda: "../envs/utils.yaml"
    localrule: True
    shell:'''
wget -O - https://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt | awk 'NR>2' > {output.meta}
cut -f1 {output.meta} | datasets download genome accession --inputfile - --filename {output.folder} --include protein
'''

rule get_virus_class:
    input: rules.get_virus_genomes.output.meta
    output: "results/meta/refseq_virus_groups.tsv"
    localrule: True
    conda: "../envs/utils.yaml"
    shell: """
cut -f6 {input} | taxonkit reformat -I 1 -f "{{r}}\\t{{K}}\\t{{p}}" | cut -f2- | sort -u | grep . > {output}
"""

# Get EukProt metadata
rule get_eukprot:
    output:
        euk_included="results/meta/EukProt_included.tsv",
        euk_excluded="results/meta/EukProt_not_included.tsv",
        euk_busco="results/meta/EukProt_busco.tsv",
        euk_fa="results/tmp/eukprot.tgz",
    localrule: True
    shell:'''
wget -O - https://figshare.com/ndownloader/files/34436249 | sed 's/\\"//g' > {output.euk_excluded}
wget -O - https://figshare.com/ndownloader/files/34436246 | sed 's/\\"//g' > {output.euk_included}
wget --no-check-certificate -O {output.euk_busco} https://evocellbio.com/SAGdb/images/EukProtv3.busco.output.txt
wget -O {output.euk_fa} https://figshare.com/ndownloader/files/34434377
'''

# Get UniProt taxonomy and annotate it with unieuk/EukProt
rule get_uniprot_meta:
    output:
        meta="results/meta/uniprot_proteomes.tsv",
        stats="results/meta/uniprot_busco.tsv",
        ncbi_tax="results/taxonomies/uniprot_ncbi_taxonomy.tsv",
        lineage="results/tmp/uniprot_lineage.tsv"
    params:
        taxdump=rules.download_taxdump.output
    conda: "../envs/utils.yaml"
    localrule: True
    shell:'''
wget -O - https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README | \
grep ^UP | awk '$4=="eukaryota"' > {output.meta}

wget "https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd&format=tsv&query=%28*%29+AND+%28proteome_type%3A1%29+AND+%28superkingdom%3AEukaryota%29" \
-O {output.stats}
cut -f1,2 {output.meta} | taxonkit reformat -I 2 -r "" -P --data-dir {params.taxdump} | \
cut -f 1,3- | sed 's/k__/d__/' > {output.ncbi_tax}

cut -f1,2 {output.meta} | taxonkit lineage -i 2 --data-dir {params.taxdump} > {output.lineage}
'''
# cat {output.meta} | cut -f1,2 | taxonkit lineage -i 2 --data-dir {input.taxdump} | \
# Rscript src/get_tax.R -u {input.unieuk} --ep {input.ep} --notep {input.notep} | tr ';' '\\t' | python src/fill_na.py | \
# tr '\\t' ';' | sed 's/;/\\t/' > {output.eukprot_tax}


# Download p10k metadata
rule get_p10k:
    output:
        sample="results/tmp/p10k_sample.tsv",
        tax="results/tmp/p10k_tax.tsv",
        assembly="results/tmp/p10k_assembly.tsv",
        genome="results/tmp/p10k_genome.tsv",
        annotation="results/tmp/p10k_annotation.tsv",
        meta="results/meta/p10k_meta.tsv",
        lineage="results/meta/p10k_lineage.tsv"
    conda: "../envs/utils.yaml"
    localrule: True
    shell:'''
echo -e "p10k_id\\tbiosample\\tsource\\tspecies" > {output.sample}
wget https://ngdc.cncb.ac.cn/p10k/api/sample/list -O - | \
jq -r  '.data[] | [.sampleId, .biosampleId, .source, .species] | @tsv' >> {output.sample}

echo -e "p10k_id\\tssu_identity\\tlineage" > {output.tax}
wget https://ngdc.cncb.ac.cn/p10k/api/taxonomy/list -O - | \
jq -r  '.data[] | [.sampleId, .ssuIdentity, .referableLineage] | @tsv' >> {output.tax}

echo -e "p10k_id\\tassembly_id\\tsize\\tn_contigs\\tN50\\tcompleteness\\tn_genes\\tCDS_completeness\\tannotation_level" > {output.assembly}
wget https://ngdc.cncb.ac.cn/p10k/api/assembly/list -O - | \
jq -r  '.[] | @tsv' | cut -f2-10 >> {output.assembly}

echo -e "p10k_id\\tplatform\\tstrategy" > {output.genome}
wget https://ngdc.cncb.ac.cn/p10k/api/sequencing/list -O - | \
jq -r '.[] | [.sampleId,.sequencingPlatform,.sequencingStrategy] | @tsv' >> {output.genome}

echo -e "p10k_id\\taverage_gene_length\\taverage_cds_length\\taverage_exon_per_gene\\taverage_exon_length\\tcodon_table" > {output.annotation}
wget https://ngdc.cncb.ac.cn/p10k/api/annotation/list -O - | \
jq -r '.[] | [.sampleId,.averageGeneLength,.averageCdsLength,.averageExonPerGene,.averageExonLength,.codonTable] | @tsv' >> {output.annotation}

csvtk join -t -f1 {output.sample} {output.tax} {output.genome} {output.assembly} {output.annotation} | \
csvtk uniq -t -f1 > {output.meta}

awk 'NR>1' {output.meta} | cut -f1,4,6 > {output.lineage}
'''
# echo -e "id\\tp10k_id\\tbiosample\\tsource\\tspecies\\tlineage\
# \\tassembly_id\\tsize\\tn_contigs\\tN50\\tcompleteness\\t\
# n_genes\\tCDS_completeness\\tannotation_level" > {output.meta}

# rule get_nr_accessionmap:
#     output: "results/tmp/prot.accession2taxid.gz"
#     shell:'''
# wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
# '''
