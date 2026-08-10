# Fetch repdb's heavy CONSTRUCTION artifacts from a pinned Zenodo record
# instead of rebuilding them from raw sources (dbs.build.repdb.zenodo).
#
# This is a second, further-downstream seam than the universe pin in
# universe.smk. A pinned universe still re-downloads and re-parses every
# source proteome and re-runs decontamination-detection + taxon-aware
# clustering on the deterministic selection; a pinned Zenodo record skips all
# of that, retrieving its outputs directly. Every fetch rule below has the
# SAME (non-wildcarded) output path as the construction rule it replaces, and
# wins over it via `ruleorder` - so nothing downstream (make_diamonddb,
# make_mmseqsdb, make_blastdb, publish_clustered_variant, or a custom db's own
# decontamination, which reads results/dbs/repdb/repdb.fa.gz as its reference)
# needs to know or care whether repdb's files were built or fetched.
#
# The clustered fasta is deliberately NOT hosted separately - it is derived
# locally from the raw fasta + the cluster membership table
# (repdb_clusters.tsv): a cluster's representative is, by mmseqs' own
# convention, whichever ID appears in column 1 of its row (the same IDs
# `{clade}_rep_seq.fasta` would contain), so `cut -f1 | sort -u` on that table
# reproduces the exact representative-ID set `merge_clustered` would otherwise
# get by actually re-running `mmseqs easy-linclust` per clade (~2h at repdb's
# scale - see docs/releasing.md). One hosted fasta thus reproduces BOTH repdb
# and repdb_clustered, and this is also why the paper's own point that the
# clustered set loses little (mirroring NCBI's nr -> ClusteredNR move) matters
# here beyond the manuscript: it's what makes a single raw-fasta download
# sufficient to reconstruct both released databases.
#
# Every fetch is checksum-verified against the sha256 pinned in
# dbs.build.repdb.zenodo.files (not fetched live from Zenodo's API), so a
# release config is self-verifying without depending on Zenodo's API still
# behaving the same way whenever this is rebuilt.
#
# The taxdump fetch is NOT here: create_taxdump's output has no wildcards, so
# rules referencing it as `rules.create_taxdump.output.full_taxdump`
# (make_diamonddb, make_mmseqsdb, get_clade) bind to that rule directly at
# parse time - ruleorder can't redirect them. Its zenodo branch lives at the
# rule's own definition in taxonomy.smk instead (same `if REPDB_ZENODO:`
# pattern universe.smk already uses for build_universe). `_zenodo_file` is
# defined there too, right next to REPDB_ZENODO, and reused here.


if REPDB_ZENODO:

    rule fetch_repdb_fasta:
        """Replaces make_db_fasta for repdb: the raw assembled fasta (also the
        decontamination reference custom dbs cluster against, via
        concat_fasta_decont) and its accession map, fetched instead of parsed
        from source proteomes."""
        output:
            fa="results/dbs/repdb/repdb.fa.gz",
            idmap="results/dbs/repdb/repdb_accession_map.txt",
        params:
            fa_url=_zenodo_file("fasta")[0],
            fa_sha=_zenodo_file("fasta")[1],
            map_url=_zenodo_file("accession_map")[0],
            map_sha=_zenodo_file("accession_map")[1],
        log:
            "results/log/downloads/zenodo_repdb_fasta.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
bash workflow/scripts/fetch_zenodo.sh {params.fa_url} {params.fa_sha} {output.fa} >{log} 2>&1
bash workflow/scripts/fetch_zenodo.sh {params.map_url} {params.map_sha} {output.idmap} >>{log} 2>&1
"""

    ruleorder: fetch_repdb_fasta > make_db_fasta

    rule fetch_repdb_map:
        """Replaces make_db_map for repdb."""
        output:
            headermap="results/dbs/repdb/repdb.map",
            noheadermap="results/dbs/repdb/repdb_nohead.map",
        params:
            h_url=_zenodo_file("map")[0],
            h_sha=_zenodo_file("map")[1],
            n_url=_zenodo_file("nohead_map")[0],
            n_sha=_zenodo_file("nohead_map")[1],
        log:
            "results/log/downloads/zenodo_repdb_map.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
bash workflow/scripts/fetch_zenodo.sh {params.h_url} {params.h_sha} {output.headermap} >{log} 2>&1
bash workflow/scripts/fetch_zenodo.sh {params.n_url} {params.n_sha} {output.noheadermap} >>{log} 2>&1
"""

    ruleorder: fetch_repdb_map > make_db_map

    rule fetch_repdb_contaminants:
        """Replaces get_contaminants for repdb: skips re-running the whole
        decontamination-DETECTION clustering chain (concat_fasta_decont ->
        cluster_decontaminate -> remove_singletons -> get_mixed_clusters ->
        get_pairwise_combination), which is independent of, and comparable in
        cost to, the taxon-aware clustering step below."""
        output:
            df="results/dbs/repdb/decontaminate/contaminants.tsv",
            ids="results/dbs/repdb/decontaminate/contaminants.txt",
        params:
            df_url=_zenodo_file("contaminants_tsv")[0],
            df_sha=_zenodo_file("contaminants_tsv")[1],
            ids_url=_zenodo_file("contaminants_ids")[0],
            ids_sha=_zenodo_file("contaminants_ids")[1],
        log:
            "results/log/downloads/zenodo_repdb_contaminants.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
bash workflow/scripts/fetch_zenodo.sh {params.df_url} {params.df_sha} {output.df} >{log} 2>&1
bash workflow/scripts/fetch_zenodo.sh {params.ids_url} {params.ids_sha} {output.ids} >>{log} 2>&1
"""

    ruleorder: fetch_repdb_contaminants > get_contaminants

    rule derive_repdb_decontaminated_fasta:
        """Replaces decontaminate_db for repdb. RepDB v1.0 uses the default
        *soft* filter (flag, don't remove - see docs/releasing.md), so the
        decontaminated fasta is byte-identical to the raw one; matching
        decontaminate_db's own soft-mode behaviour, this is a relative
        symlink, not a copy. A release built with `filter: hard` would need
        its already-filtered fasta hosted under its own Zenodo key instead of
        reusing "fasta" here - not needed for v1.0, so not implemented."""
        input:
            fa=rules.fetch_repdb_fasta.output.fa,
            # contaminants.txt isn't actually read by the soft-mode shell logic
            # below (decontaminate_db doesn't read it in soft mode either - it's
            # a dependency edge only), but is declared so Snakemake still
            # schedules fetch_repdb_contaminants when this file is requested.
            contaminants=rules.fetch_repdb_contaminants.output.ids,
        output:
            "results/dbs/repdb/repdb_decontaminated.fa.gz",
        localrule: True
        shell:
            "ln -rsf {input.fa} {output} 2>/dev/null || cp {input.fa} {output}"

    ruleorder: derive_repdb_decontaminated_fasta > decontaminate_db

    rule fetch_repdb_clusters:
        """Replaces merge_clustered's `clusters` output for repdb."""
        output:
            clusters="results/dbs/repdb/repdb_clusters.tsv",
        params:
            url=_zenodo_file("clusters")[0],
            sha=_zenodo_file("clusters")[1],
        log:
            "results/log/downloads/zenodo_repdb_clusters.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            "bash workflow/scripts/fetch_zenodo.sh {params.url} {params.sha} {output.clusters} >{log} 2>&1"

    ruleorder: fetch_repdb_clusters > merge_clustered

    rule derive_repdb_clustered_fasta:
        """Replaces merge_clustered's `seqs` output for repdb: reconstructs the
        clustered fasta from the raw one instead of re-running
        `mmseqs easy-linclust` per clade - see module docstring above for why
        column 1 of the cluster tsv is exactly the representative-ID set."""
        input:
            fa=rules.fetch_repdb_fasta.output.fa,
            clusters=rules.fetch_repdb_clusters.output.clusters,
        output:
            seqs="results/dbs/repdb/repdb_clustered.fa.gz",
        threads: 8
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
repids=$(mktemp)
cut -f1 {input.clusters} | sort -u > "$repids"
seqkit grep -f "$repids" {input.fa} -o {output.seqs} --threads {threads}
rm -f "$repids"
"""

    ruleorder: derive_repdb_clustered_fasta > merge_clustered
