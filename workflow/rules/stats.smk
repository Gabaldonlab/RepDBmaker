rule make_repdb_stats:
    input: rules.create_repdb_table.output
    output: "results/meta/repdb_stats.tsv"
    threads: 48
    shell:'''
cut -f1 {input} | seqkit stats -j {threads} --infile-list - -T -b | \
cut -f1,4- | sed 's/.faa.gz//' > {output}
'''


rule make_repdb_meta:
    input:
        tax=rules.create_repdb_taxdump.output.all_taxa,
        stats=rules.make_repdb_stats.output,
        up_stats=rules.get_uniprot_meta.output.stats,
        ep_stats=rules.get_eukprot.output.euk_busco,
        ep=rules.get_eukprot.output.euk_included,
        p10k_stats=rules.get_p10k.output.meta,
        custom_table=config["new_genomes"],
        custom_busco=config["new_genomes_busco"],
        contaminants=rules.get_contaminants.output
    output: "results/meta/repdb_meta.tsv"
    # log: "results/log/db/parse.log"
    # benchmark: "results/benchmarks/db/parse.txt"
    script: "../scripts/analyze_stats.R"

