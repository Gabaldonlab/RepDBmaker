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
        rules.create_p10k_table.output,
        rules.create_uniprot_table.output
    output: "results/meta/check_resources.txt"
    shell:'''
for file in {input}; do
    if [ -s "$file" ]; then
        echo "$file ok"
    else
        echo "$file NOT FOUND"
    fi
done > {output}
'''


checkpoint get_euka_tax:
    input:
        up=rules.get_uniprot_tax.output,
        ep=rules.get_eukprot_tax.output,
        p10k=rules.get_p10k_tax.output,
        custom=rules.get_custom_tax.output,
        up_stats=rules.get_uniprot_meta.output.stats,
        ep_stats=rules.get_eukprot.output.euk_busco,
        p10k_stats=rules.get_p10k.output.meta,
        exclude=config["genomes_to_exclude"],
        to_keep=config["clades_to_keep"]
    output:
        tax="results/taxonomies/eukaryotes.tsv",
        plot="results/plots/euka_db.pdf"
    script:"../scripts/filter_tax.R"

rule get_all_euka_tax:
    input:
        up=rules.get_uniprot_tax.output,
        ep=rules.get_eukprot_tax.output,
        p10k=rules.get_p10k_tax.output,
        custom=rules.get_custom_tax.output
    output: "results/taxonomies/all_eukaryotes.tsv"
    shell: "cat {input} | sort -k2,2 > {output}"

# final merging of all the dbs

rule create_repdb_taxdump:
    input:
        gtdb_proteomes=rules.decompress_gtdb_genomes.output,
        virus=rules.check_virus_genomes.output.tax,
        gtdb=rules.get_gtdb_tax.output.tax,
        euka=rules.get_euka_tax.output.tax
    output:
        ids="results/meta/repdb.ids",
        all_taxa="results/taxonomies/repdb.tsv",
        repdb_taxdump=directory("results/taxdump/repdb_taxdump/")
    shell:'''
find $(dirname {input.gtdb_proteomes}) -type f -name "*faa.gz" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort > {output.ids}
cat {input.virus} {input.gtdb} {input.euka} | csvtk join -H -t {output.ids} - | \
sed 's/d__//g' | sed 's/[a-z]__/\\t/g' | sed 's/;//g' > {output.all_taxa}
taxonkit create-taxdump -A1 {output.all_taxa} --out-dir {output.repdb_taxdump} --force \
--rank-names "superkingdom","phylum","class","order","family","genus","species"
'''

# create genome table for all 3 databases then user inputs custom_genome table, concatenate and then parse!
rule create_repdb_table:
    input:
        rules.create_custom_table.output,
        rules.create_eukprot_table.output,
        rules.create_uniprot_table.output,
        rules.create_p10k_table.output,
        rules.create_virus_table.output,
        rules.create_gtdb_table.output
    output: "results/meta/repdb_genome_table.tsv"
    shell:'''
cat {input} > {output}
'''

rule make_repdb_fasta:
    input:
        table=rules.create_repdb_table.output,
        taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
    output:
        fa="results/db/repdb.fa",
        idmap="results/db/repdb_accession_map.txt"
    script: "../scripts/parse_gnm.py"

rule make_repdb_map:
    input:
        idmap=rules.make_repdb_fasta.output.idmap
    output:
        headermap="results/db/repdb.map",
        noheadermap="results/db/repdb_nohead.map"
    shell:'''
echo -e "accession.version\\ttaxid" > {output.headermap}
cut -f2 {input.idmap} | awk -F'\\t' '{{split($NF, a, "_"); print $0"\\t"a[1]}}' >> {output.headermap}
awk 'NR>1' {output.headermap} > {output.noheadermap}
'''

rule make_blastdb:
    input: 
        fa=rules.make_repdb_fasta.output.fa,
        taxidmap=rules.make_repdb_map.output.noheadermap
    output: "results/db/repdb_blast"
    shell:'''
makeblastdb -in {input.fa} -parse_seqids -taxid_map {input.taxidmap} -dbtype prot -out {output}
touch {output}
'''

rule make_diamonddb:
    input: 
        fa=rules.make_repdb_fasta.output.fa,
        taxidmap=rules.make_repdb_map.output.headermap,
        taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
    output: "results/db/repdb_diamond"
    threads: 48
    shell:'''
diamond makedb --in {input.fa} -d {output} --taxonmap {input.taxidmap} --taxonnodes {input.taxdump}/nodes.dmp --taxonnames {input.taxdump}/names.dmp --threads {threads}
touch {output}
'''

rule make_mmseqsdb:
    input: 
        fa=rules.make_repdb_fasta.output.fa,
        taxidmap=rules.make_repdb_map.output.noheadermap,
        taxdump=rules.create_repdb_taxdump.output.repdb_taxdump
    output: "results/db/repdb_mmseqs"
    threads: 48
    shell:'''
mmseqs createdb {input.fa} {output}
mmseqs createtaxdb {output} $TMPDIR --ncbi-tax-dump {input.taxdump} --tax-mapping-file {input.taxidmap} --threads {threads}
'''
