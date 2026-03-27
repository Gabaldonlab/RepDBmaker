rule paths_gtdb:
    input: rules.get_gtdb_tax.output.meta
    output: "results/meta/gtdb_paths.txt",
    conda: "../envs/utils.yaml"
    # localrule: True
    group: "gtdb"
    shell:'''
csvtk cut -t -f accession,gtdb_taxonomy {input} | awk 'NR>1' | sed 's/;.*\\|d__//g' | \
awk -F"\\t" '{{l = $0; sub($2, "", l); print "protein_faa_reps/"tolower($2)"/"$1"_protein.faa.gz"}}' > {output}
'''

# def prokaryotes_to_extract(wildcards):
#     with open(str(checkpoints.paths_gtdb.get(**wildcards).output)) as proka:
#         genomes = [''.join(gn.strip().split('/')[-1].split('.')[0].split('_')[1:3]) for gn in proka][1:10]
#         return expand("results/proteomes/gtdb/{i}.faa.gz", i=genomes)

rule decompress_gtdb_genomes:
    input: 
        db=rules.get_gtdb_genomes.output,
        files=rules.paths_gtdb.output
    output: 
        folder=directory("results/proteomes/gtdb"),
        # is_done="proteomes/gtdb/.is_extracted"
    group: "gtdb"
    shell:'''
mkdir -p {output.folder}
tar -xzvf {input.db} -C {output.folder} --files-from={input.files} \
--strip-components=2 --transform 's/RS_\\|GB_\\|_\\|.[0-9]_protein//g'
'''

rule gtdb_species_clusters:
    input: rules.decompress_gtdb_genomes.output
    output: "results/meta/gtdb.ids"
    conda: "../envs/utils.yaml"
    localrule: True
    shell:'''
find {input} -type f -name "*faa.gz" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort > {output}
'''

# create genome table for all 3 databases then user inputs custom_genome table, concatenate and then parse!
# rule create_gtdb_table:
#     input: 
#         gtdb_tax=rules.get_gtdb_tax.output.tax,
#         gtdb=rules.decompress_gtdb_genomes.output.folder # prokaryotes_to_extract in theory
#     output: "results/meta/gtdb_genome_table.tsv"
#     conda: "../envs/utils.yaml"
#     group: "gtdb"
#     shell:'''
# > {output}
# find $(realpath {input.gtdb}) -type f | \
# awk '{{split($1,a,"/"); print $1"\\t"a[length(a)]}}' | \
# sed 's/.faa.gz$//' | csvtk join -f"2;1" -H -t - {input.gtdb_tax} | \
# sed 's/d__\\|;.*//g' >> {output}
# '''
# for genome in {input.gtdb}; do
#     bn=$(basename $genome ".faa.gz")
#     path=$(realpath $genome)
#     echo -e "$path\\t$bn\\tProkaryota"
# done >> {output}
