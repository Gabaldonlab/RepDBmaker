if REPDB_ZENODO:

    rule fetch_repdb_fasta:
        """Replaces make_db_fasta's `fa` output for repdb: the raw assembled
        fasta, also the decontamination reference custom dbs cluster against
        (via concat_fasta_decont). Its `idmap` output (repdb_accession_map.txt,
        original-source-ID -> repdb-ID) is deliberately NOT fetched - nothing
        downstream of it needs it (it only ever fed make_db_map, which is now
        derived straight from the fasta instead, see derive_repdb_map below),
        so it's pure traceability with no functional use in this DAG. It's
        still produced by a real (non-Zenodo) build; a Zenodo-only rebuild
        just doesn't have it."""
        output:
            fa="results/dbs/repdb/repdb.fa.gz",
        params:
            fa_url=_zenodo_file("fasta")[0],
            fa_sha=_zenodo_file("fasta")[1],
        log:
            "results/log/downloads/zenodo_repdb_fasta.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            "bash workflow/scripts/fetch_zenodo.sh {params.fa_url} {params.fa_sha} {output.fa} >{log} 2>&1"

    ruleorder: fetch_repdb_fasta > make_db_fasta

    rule derive_repdb_map:
        """Replaces make_db_map for repdb - not fetched at all. Both `.map`
        and `_nohead.map` are just `{new_id}\\t{taxid}` per record, and the
        fasta headers already ARE `new_id` (={taxid}_{code}_AA###) - the taxid
        is simply the header's own first `_`-split token, exactly what
        make_db_map's `awk '{{split($NF,a,"_"); print $0"\\t"a[1]}}'` extracts
        from repdb_accession_map.txt today. No Zenodo asset needed for this."""
        input:
            fa=rules.fetch_repdb_fasta.output.fa,
        output:
            headermap="results/dbs/repdb/repdb.map",
            noheadermap="results/dbs/repdb/repdb_nohead.map",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
echo -e "accession.version\\ttaxid" > {output.headermap}
zcat {input.fa} | grep '^>' | sed 's/^>//' | awk -F'_' '{{print $0"\\t"$1}}' >> {output.headermap}
awk 'NR>1' {output.headermap} > {output.noheadermap}
"""

    ruleorder: derive_repdb_map > make_db_map

    rule fetch_repdb_contaminants:
        """Replaces get_contaminants for repdb: skips re-running the whole
        decontamination-DETECTION clustering chain (concat_fasta_decont ->
        cluster_decontaminate -> remove_singletons -> get_mixed_clusters ->
        get_pairwise_combination), which is independent of, and comparable in
        cost to, the taxon-aware clustering step below. Only the dataframe
        (`df`) is fetched; `contaminants.txt` is just its `seq` column (col 1),
        the same as get_contaminants.R's own `writeLines(contaminants$seq, ...)`
        - deriving it locally avoids hosting a second, redundant file."""
        output:
            df="results/dbs/repdb/decontaminate/contaminants.tsv",
            ids="results/dbs/repdb/decontaminate/contaminants.txt",
        params:
            url=_zenodo_file("contaminants_tsv")[0],
            sha=_zenodo_file("contaminants_tsv")[1],
        log:
            "results/log/downloads/zenodo_repdb_contaminants.log",
        localrule: True
        conda:
            "../envs/utils.yaml"
        shell:
            """
bash workflow/scripts/fetch_zenodo.sh {params.url} {params.sha} {output.df} >{log} 2>&1
awk -F'\\t' 'NR>1{{print $1}}' {output.df} > {output.ids}
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
        """Replaces merge_clustered's `clusters` output for repdb. Hosted
        gzipped: a plain 2-column TSV that repeats each cluster's
        representative ID once per member compresses enormously (~10x on
        RepDB v1.0), so it's fetched compressed, checksum-verified against
        that compressed form, then decompressed to the plain TSV everything
        downstream (derive_repdb_clustered_fasta) expects."""
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
            """
bash workflow/scripts/fetch_zenodo.sh {params.url} {params.sha} {output.clusters}.gz >{log} 2>&1
gunzip -f {output.clusters}.gz
"""

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
