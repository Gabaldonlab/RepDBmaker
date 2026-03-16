### Get EukProt

rule get_eukprot_tax:
    input: rules.get_eukprot.output.euk_included
    output: "results/taxonomies/eukprot_taxonomy.tsv"
    localrule: True
    shell: '''
awk 'NR>1' {input} | sed 's/\\_/ /g' | \
awk -F'\\t' '{{split($11, a, ";"); print $1,"d__"a[1]";p__"a[2]";c__"$8";o__"$9";f__"$10";g__"$6";s__"$2}}' OFS='\\t' | \
sed  's/N\\/A//g' | tr -d \\'\\" | sed 's/\\ ;/;/g' | sed 's/other Gyrista/other_Gyrista/g'> {output}
'''

def eukprot_to_extract(wildcards):
    with open(str(checkpoints.get_euka_tax.get(**wildcards).output.tax)) as euka:
        genomes = [gn.strip().split('\t')[0] for gn in euka if gn.startswith("EP")]
    return expand("results/proteomes/ep/{i}.faa.gz", i=genomes)


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

rule create_eukprot_table:
    input: eukprot_to_extract
    output: "results/meta/ep_genome_table.tsv"
    localrule: True
    shell:'''
> {output}
for genome in {input}; do
    bn=$(basename $genome ".faa.gz")
    path=$(realpath $genome)
    echo -e "$path\\t$bn\\tEukaryota"
done >> {output}
'''
