GTDB_BASE = "https://data.ace.uq.edu.au/public/gtdb/data/releases"
NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"

def gtdb_url(filename, version=config["versions"]["gtdb"]):
    """Build a GTDB download URL, handling the two layouts GTDB uses.

    `filename` is the base name as it appears under the `latest/` directory,
    e.g. "bac120_taxonomy.tsv"

        latest        -> releases/latest/<name>.<ext>
        release<NNN>  -> releases/release<NNN>/<NNN>.<minor>/<name>_r<NNN>.<ext>
    """
    if version == "latest":
        return f"{GTDB_BASE}/latest/{filename}"
    num = "".join(c for c in version if c.isdigit())
    head, _, tail = filename.rpartition("/")
    stem, dot, ext = tail.partition(".")
    tail_r = f"{stem}_r{num}{dot}{ext}"
    subpath = f"{head}/{tail_r}" if head else tail_r
    return f"{GTDB_BASE}/release{num}/{num}.0/{subpath}"

def taxdump_url(version=config["versions"]["taxdump"]):
    """Build a TaxDump download URL.

    `filename` is the base name as it appears under the `latest/` directory,
    e.g. "bac120_taxonomy.tsv"

        latest        -> https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        release<NNN>  -> https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_{version}.zip
    """
    if version == "latest":
        return f"{NCBI_BASE}/taxdump.tar.gz"
    else:
        return f"{NCBI_BASE}/taxdump_archive/taxdmp_{version}.zip"

# Download unieuk taxonomy


rule download_unieuk:
    params:
        unieuk_version=config["versions"]["unieuk"]
    output:
        "results/meta/unieuk_taxonomy.tsv",
    localrule: True
    log: 
        "results/log/downloads/unieuk.log"
    conda: 
        "../envs/utils.yaml"
    shell:
        """
wget https://eukmap.unieuk.net/exports/unieuk/{params.unieuk_version}/unieuk-{params.unieuk_version}.tsv -O {output} 2> {log}
"""


rule download_taxdump:
    output:
        td=directory("results/taxdump/ncbi_taxdump"),
        # acc2taxid="results/tmp/prot.accession2taxid.FULL.gz"
    localrule: True
    params: 
        taxdump_url=taxdump_url(version=config["versions"]["taxdump"]),
        taxdump_version=config["versions"]["taxdump"]
    log: 
        "results/log/downloads/taxdump.log"
    conda: 
        "../envs/utils.yaml"
    shell:
        """
mkdir -p {output.td}
if [[ {params.taxdump_version} == "latest" ]];
then
    wget -O {output.td}/ncbi_taxdump.tar.gz {params.taxdump_url} 2> {log}
    tar xf {output.td}/ncbi_taxdump.tar.gz -C {output.td}
else
    wget -O {output.td}/ncbi_taxdump.zip {params.taxdump_url} 2> {log}
    unzip -o {output.td}/ncbi_taxdump.zip -d {output.td}
fi
"""
# wget -O {output.acc2taxid} https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz


rule get_gtdb_tax:
    output:
        bac="results/tmp/bac.tmp",
        ar="results/tmp/ar.tmp",
        tax="results/taxonomies/gtdb_taxonomy.tsv",
        bac_meta="results/tmp/bac_meta.tmp",
        ar_meta="results/tmp/ar_meta.tmp",
        meta="results/meta/gtdb_meta.tsv",
    conda:
        "../envs/utils.yaml"
    log: 
        "results/log/downloads/gtdb.log"
    params:
        bac_url=gtdb_url("bac120_taxonomy.tsv"),
        ar_url=gtdb_url("ar53_taxonomy.tsv"),
        bac_meta_url=gtdb_url("bac120_metadata.tsv.gz"),
        ar_meta_url=gtdb_url("ar53_metadata.tsv.gz"),
    localrule: True
    shell:
        """
wget -O {output.bac} {params.bac_url} 2> {log}
wget -O {output.ar}  {params.ar_url} 2>> {log}
wget -O {output.bac_meta} {params.bac_meta_url} 2>> {log}
wget -O {output.ar_meta} {params.ar_meta_url} 2>> {log}

zcat {output.bac_meta} {output.ar_meta} | csvtk filter2 -t -f'$gtdb_representative=="t"' > {output.meta}

cat {output.bac} {output.ar} | awk 'BEGIN{{OFS=FS="\\t"}} {{ $1 = substr($1, 4, 13) }} 1' | \
sed 's/_//' > {output.tax}
"""


# | src/filter_gtdb.R
rule get_gtdb_genomes:
    output:
        "results/tmp/gtdb_proteins_aa_reps.tar.gz",
    params:
        url=gtdb_url("genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"),
    log:
        "results/log/downloads/gtdb_genomes.log"
    conda:
        "../envs/utils.yaml"
    shell:
        """
wget -O {output} {params.url} 2> {log}
"""


# Get refseq viruses
rule get_virus_genomes:
    output:
        meta="results/meta/refseq_virus_meta.txt",
        folder="results/tmp/refseq.zip",
    log:
        "results/log/downloads/virus_summary.log"
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        """
wget -O - https://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt | awk 'NR>2' > {output.meta} 2> {log}
cut -f1 {output.meta} | datasets download genome accession --inputfile - --filename {output.folder} --include protein
"""


rule get_virus_class:
    input:
        rules.get_virus_genomes.output.meta,
    output:
        "results/meta/refseq_virus_groups.tsv",
    localrule: True
    conda:
        "../envs/utils.yaml"
    shell:
        """
cut -f6 {input} | taxonkit reformat -I 1 -f "{{r}}\\t{{K}}\\t{{p}}" | cut -f2- | sort -u | grep . > {output}
"""


if config["versions"]["EukProt"] == "3":
    eukprot_url = "https://ndownloader.figshare.com/files/34434377"
    eukprot_meta_url = "https://ndownloader.figshare.com/files/34436246"
    eukprot_busco_url = "https://evocellbio.com/SAGdb/images/EukProtv3.busco.output.txt"
else:
    exit(f"EukProt version {config['versions']['EukProt']} not supported")

# Get EukProt metadata
rule get_eukprot:
    output:
        euk_included="results/meta/EukProt_included.tsv",
        euk_excluded="results/meta/EukProt_not_included.tsv",
        euk_busco="results/meta/EukProt_busco.tsv",
        euk_fa="results/tmp/eukprot.tgz",
    localrule: True
    log: 
        "results/log/downloads/eukprot.log"
    conda: 
        "../envs/utils.yaml"
    shell:
        """
wget -O - {eukprot_url} | sed 's/\\"//g' > {output.euk_excluded} 2> {log}
wget -O - {eukprot_meta_url} | sed 's/\\"//g' > {output.euk_included} 2>> {log}
wget --no-check-certificate -O {output.euk_busco} {eukprot_busco_url} 2>> {log}
wget -O {output.euk_fa} {eukprot_url} 2>> {log}
"""


# Get UniProt taxonomy and annotate it with unieuk/EukProt
rule get_uniprot_meta:
    output:
        meta="results/meta/uniprot_proteomes.tsv",
        stats="results/meta/uniprot_busco.tsv",
        ncbi_tax="results/taxonomies/uniprot_ncbi_taxonomy.tsv",
        lineage="results/tmp/uniprot_lineage.tsv",
    input:
        # declared as input (not params) so Snakemake enforces the dependency
        # on download_taxdump and avoids a race condition
        taxdump=rules.download_taxdump.output,
    log:
        "results/log/downloads/uniprot.log"
    conda:
        "../envs/utils.yaml"
    localrule: True
    shell:
        """
wget -O - https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README | \
grep ^UP | awk '$4=="eukaryota"' > {output.meta} 2> {log}

wget "https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd&format=tsv&query=%28*%29+AND+%28proteome_type%3A1%29+AND+%28superkingdom%3AEukaryota%29" \
-O {output.stats} 2>> {log}
cut -f1,2 {output.meta} | taxonkit reformat -I 2 -r "" -P --data-dir {input.taxdump} | \
cut -f 1,3- | sed 's/k__/d__/' > {output.ncbi_tax}

cut -f1,2 {output.meta} | taxonkit lineage -i 2 --data-dir {input.taxdump} > {output.lineage}
"""


rule get_p10k:
    output:
        sample="results/tmp/p10k_sample.tsv",
        tax="results/tmp/p10k_tax.tsv",
        assembly="results/tmp/p10k_assembly.tsv",
        genome="results/tmp/p10k_genome.tsv",
        annotation="results/tmp/p10k_annotation.tsv",
        meta="results/meta/p10k_meta.tsv",
        lineage="results/meta/p10k_lineage.tsv",
    conda:
        "../envs/utils.yaml"
    log: 
        "results/log/downloads/unieuk.log"
    localrule: True
    shell:
        """
echo -e "p10k_id\\tbiosample\\tsource\\tspecies" > {output.sample}
wget https://ngdc.cncb.ac.cn/p10k/api/sample/list -O - | \
jq -r  '.data[] | [.sampleId, .biosampleId, .source, .species] | @tsv' >> {output.sample} 2> {log}

echo -e "p10k_id\\tssu_identity\\tlineage" > {output.tax}
wget https://ngdc.cncb.ac.cn/p10k/api/taxonomy/list -O - | \
jq -r  '.data[] | [.sampleId, .ssuIdentity, .referableLineage] | @tsv' >> {output.tax} 2>> {log}

echo -e "p10k_id\\tassembly_id\\tsize\\tn_contigs\\tN50\\tcompleteness\\tn_genes\\tCDS_completeness\\tannotation_level" > {output.assembly}
wget https://ngdc.cncb.ac.cn/p10k/api/assembly/list -O - | \
jq -r  '.[] | @tsv' | cut -f2-10 >> {output.assembly} 2>> {log}

echo -e "p10k_id\\tplatform\\tstrategy" > {output.genome}
wget https://ngdc.cncb.ac.cn/p10k/api/sequencing/list -O - | \
jq -r '.[] | [.sampleId,.sequencingPlatform,.sequencingStrategy] | @tsv' >> {output.genome} 2>> {log}

echo -e "p10k_id\\taverage_gene_length\\taverage_cds_length\\taverage_exon_per_gene\\taverage_exon_length\\tcodon_table" > {output.annotation}
wget https://ngdc.cncb.ac.cn/p10k/api/annotation/list -O - | \
jq -r '.[] | [.sampleId,.averageGeneLength,.averageCdsLength,.averageExonPerGene,.averageExonLength,.codonTable] | @tsv' >> {output.annotation} 2>> {log}

csvtk join -t -f1 {output.sample} {output.tax} {output.genome} {output.assembly} {output.annotation} | \
csvtk uniq -t -f1 > {output.meta}

awk 'NR>1' {output.meta} | cut -f1,4,6 > {output.lineage}
"""


# echo -e "id\\tp10k_id\\tbiosample\\tsource\\tspecies\\tlineage\
# \\tassembly_id\\tsize\\tn_contigs\\tN50\\tcompleteness\\t\
# n_genes\\tCDS_completeness\\tannotation_level" > {output.meta}
# rule get_nr_accessionmap:
#     output: "results/tmp/prot.accession2taxid.gz"
#     shell:'''
# wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
# '''
