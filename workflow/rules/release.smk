# Freeze a completed run into a versioned, reproducible RepDB release.
#
#   snakemake resources/releases/v1/release.yaml   # rule make_release
#   snakemake resources/releases/v1/zenodo_stage/zenodo_config.yaml   # rule stage_zenodo_assets
#
# Three hosting tiers, split by what each artifact IS rather than by size -
# see docs/releasing.md for the full reasoning:
#   git (this repo)   - config.yaml, repdb.ids, release.yaml: small, versioned
#                        text, committed directly at the release tag.
#   GitHub release     - universe.tsv, custom_bundle.tar.gz (rule make_release):
#                        Pipeline 1 *inputs* - needed only to re-run curation
#                        from scratch, tied to a specific software version,
#                        too heavy for git history.
#   Zenodo              - repdb.fa.gz + friends (rule stage_zenodo_assets): the
#                        actual data and results, DOI'd and citable.
#
# make_release writes, under resources/releases/<version>/:
#   config.yaml   - a self-contained config pinning the universe + custom bundle
#   repdb.ids     - the selected composition (tracked; human-diffable changelog)
#   release.yaml  - the manifest (software version, git commit, asset checksums,
#                   and the GitHub-release asset URLs)
#   universe.tsv, custom_bundle.tar.gz  - the heavy GitHub-release assets (gitignored)
#
# Software version = VERSION (currently 0.0.1, pre-release); database version =
# the <version> wildcard (e.g. v1).

REPO_SLUG = "Gabaldonlab/RepDBmaker"


rule make_release:
    input:
        universe="results/universe/universe.tsv",
        ids="results/universe/repdb.ids",
        base_config="config/repdb.yaml",
    output:
        config="resources/releases/{version}/config.yaml",
        ids="resources/releases/{version}/repdb.ids",
        manifest="resources/releases/{version}/release.yaml",
    localrule: True
    run:
        import os, shutil, hashlib, subprocess, yaml

        version = wildcards.version
        reldir = f"resources/releases/{version}"
        os.makedirs(reldir, exist_ok=True)

        def sha256(path):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            return h.hexdigest()

        # --- heavy assets (gitignored; uploaded to the GitHub release) ---------
        shutil.copyfile(input.universe, f"{reldir}/universe.tsv")
        shutil.copyfile(input.ids, output.ids)
        assets = {"universe": {"file": "universe.tsv",
                               "sha256": sha256(f"{reldir}/universe.tsv")}}

        bundle_pinned = None
        if CUSTOM_BUNDLE and os.path.exists(CUSTOM_BUNDLE):
            shutil.copyfile(CUSTOM_BUNDLE, f"{reldir}/custom_bundle.tar.gz")
            assets["custom_bundle"] = {"file": "custom_bundle.tar.gz",
                                       "sha256": sha256(f"{reldir}/custom_bundle.tar.gz")}
            bundle_pinned = f"{reldir}/custom_bundle.tar.gz"

        # --- self-contained release config: pin universe (+ bundle) ------------
        with open(input.base_config) as fh:
            cfg = yaml.safe_load(fh)
        repdb = cfg.setdefault("dbs", {}).setdefault("build", {}).setdefault("repdb", {})
        repdb["universe"] = f"{reldir}/universe.tsv"
        if bundle_pinned:
            repdb["custom_bundle"] = bundle_pinned
        cfg["test"] = False  # never freeze a release in test mode
        with open(output.config, "w") as fh:
            yaml.safe_dump(cfg, fh, sort_keys=False)

        # --- manifest ----------------------------------------------------------
        software = open("VERSION").read().strip()
        try:
            commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
        except Exception:
            commit = "unknown"
        base_url = f"https://github.com/{REPO_SLUG}/releases/download/{version}"
        manifest = {
            "repdb_version": version,
            "software": {"version": software, "git_commit": commit},
            "config": f"{reldir}/config.yaml",
            "selection": {"file": "repdb.ids",
                          "n": sum(1 for _ in open(output.ids))},
            "assets": {k: {**v, "url": f"{base_url}/{v['file']}"}
                       for k, v in assets.items()},
        }
        with open(output.manifest, "w") as fh:
            yaml.safe_dump(manifest, fh, sort_keys=False)

        # --- publish: print copy-paste commands --------------------------------
        upload = " ".join(f"{reldir}/{a['file']}" for a in assets.values())
        light = f"{output.config} {output.ids} {output.manifest}"
        print(f"\n=== RepDB {version} staged in {reldir}/ (software {software}) ===\n")
        print("1) upload the heavy assets to the GitHub release (not git):")
        print(f"   gh release create {version} {upload} "
              f"--title 'RepDB {version}' --generate-notes\n")
        print("2) commit the light, versionable files (git):")
        print(f"   git add {light} && git commit -m 'RepDB {version} release'\n")


rule stage_zenodo_assets:
    """After `snakemake build` has produced repdb's outputs, stage everything
    for a Zenodo upload: symlink/copy the files into one directory, checksum
    them, and write SHA256SUMS.txt + a ready-to-paste zenodo_config.yaml.

    Staged: repdb.fa.gz, repdb_clusters.tsv.gz, repdb_contaminants.tsv,
    repdb_taxdump.tar.gz, plus repdb_meta.tsv and the two
    stats tables

    repdb_clusters.tsv itself isn't compressed at its source path, but it's
    a highly repetitive 2-column TSV (each cluster's representative ID
    repeated once per member) that gzips down enormously, so it's gzipped
    during staging rather than symlinked as-is. fetch_repdb_clusters (in
    zenodo_fetch.smk) decompresses it back on the way in.
    """
    input:
        fa="results/dbs/repdb/repdb.fa.gz",
        clusters="results/dbs/repdb/repdb_clusters.tsv",
        contaminants_tsv="results/dbs/repdb/decontaminate/contaminants.tsv",
        taxdump="results/taxdump/repdb_taxdump/",
        repdb_meta="results/meta/repdb_meta.tsv",
        stats="results/stats/repdb_stats.tsv",
        clustered_stats="results/stats/repdb_clustered_stats.tsv",
    output:
        sha256sums="resources/releases/{version}/zenodo_stage/SHA256SUMS.txt",
        config_block="resources/releases/{version}/zenodo_stage/zenodo_config.yaml",
    localrule: True
    run:
        import gzip
        import hashlib
        import os
        import shutil
        import tarfile

        import yaml

        stagedir = f"resources/releases/{wildcards.version}/zenodo_stage"
        os.makedirs(stagedir, exist_ok=True)

        def sha256(path):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            return h.hexdigest()

        def stage(dest_name, src, symlink):
            dest = f"{stagedir}/{dest_name}"
            if os.path.lexists(dest):
                os.remove(dest)
            if symlink:
                os.symlink(os.path.abspath(src), dest)
            else:
                shutil.copyfile(src, dest)
            return sha256(src)

        def stage_gzip(dest_name, src):
            """Gzip `src` into the stage dir, return the compressed file's
            own sha256 (not src's - that's what a downloader actually gets)."""
            dest = f"{stagedir}/{dest_name}"
            if os.path.lexists(dest):
                os.remove(dest)
            with open(src, "rb") as fin, gzip.open(dest, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            return sha256(dest)

        # key -> (Zenodo filename, source path) - what reproduce mode fetches
        fetched = {
            "fasta": ("repdb.fa.gz", input.fa),
            "clusters": ("repdb_clusters.tsv.gz", input.clusters),
            "contaminants_tsv": ("repdb_contaminants.tsv", input.contaminants_tsv),
        }
        # key -> (Zenodo filename, source path) - staged for reference only
        reference = {
            "repdb_meta": ("repdb_meta.tsv", input.repdb_meta),
            "repdb_stats": ("repdb_stats.tsv", input.stats),
            "repdb_clustered_stats": ("repdb_clustered_stats.tsv", input.clustered_stats),
        }

        checksums = {}
        for key, (name, src) in fetched.items():
            if key == "clusters":
                checksums[key] = stage_gzip(name, src)
            else:
                checksums[key] = stage(name, src, symlink=True)
        for key, (name, src) in reference.items():
            checksums[key] = stage(name, src, symlink=False)

        taxdump_name = "repdb_taxdump.tar.gz"
        taxdump_tar = f"{stagedir}/{taxdump_name}"
        with tarfile.open(taxdump_tar, "w:gz") as tf:
            tf.add(input.taxdump, arcname=".")
        checksums["taxdump"] = sha256(taxdump_tar)
        fetched["taxdump"] = (taxdump_name, taxdump_tar)

        # standard `sha256sum -c`-checkable file, uploaded to Zenodo alongside
        # the data - covers everything staged, so anyone downloading straight
        # from the record (not through our fetch mechanism at all) can verify
        # their own download rather than trusting it blindly.
        with open(output.sha256sums, "w") as fh:
            for key, (name, _) in {**fetched, **reference}.items():
                fh.write(f"{checksums[key]}  {name}\n")

        # ready to paste into resources/releases/<version>/config.yaml under
        # dbs.build.repdb, once the record exists and <record_id> is filled in.
        block = {
            "zenodo": {
                "base_url": "https://zenodo.org/records/<record_id>/files",
                "files": {key: {"name": name, "sha256": checksums[key]}
                          for key, (name, _) in fetched.items()},
            }
        }
        with open(output.config_block, "w") as fh:
            yaml.safe_dump(block, fh, sort_keys=False)

        print(f"\n=== Zenodo assets staged in {stagedir}/ ===")
        print("Upload every file in that directory to a new Zenodo deposit "
              "(web UI or the REST API - no official CLI): "
              "https://zenodo.org/deposit/new")
        print(f"Then merge {output.config_block} (filling in <record_id>) into "
              f"resources/releases/{wildcards.version}/config.yaml under dbs.build.repdb.\n")
