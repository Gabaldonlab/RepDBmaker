# Freeze a completed run into a versioned, reproducible RepDB release.
#
#   snakemake resources/releases/v1/release.yaml
#
# writes, under resources/releases/<version>/:
#   config.yaml   - a self-contained config pinning the universe + custom bundle
#   repdb.ids     - the selected composition (tracked; human-diffable changelog)
#   release.yaml  - the manifest (software version, git commit, asset checksums,
#                   and the GitHub-release asset URLs)
#   universe.tsv, custom_bundle.tar.gz  - the heavy assets (gitignored; upload as
#                   GitHub-release assets now, Zenodo later)
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
    """After `snakemake build` has produced repdb's construction artifacts,
    stage + checksum them for a Zenodo upload and print the `dbs.build.repdb.zenodo`
    config block to paste into resources/releases/<version>/config.yaml once
    the record exists (see docs/releasing.md and zenodo_fetch.smk - the
    reverse of this rule is what fetch_repdb_* consume at reproduce time).

    Symlinks rather than copies into the staging dir (these are tens of GB;
    staging must not double disk usage), except the taxdump, which is a
    directory and has to be tarred to be a single uploadable file.
    """
    input:
        fa="results/dbs/repdb/repdb.fa.gz",
        idmap="results/dbs/repdb/repdb_accession_map.txt",
        headermap="results/dbs/repdb/repdb.map",
        noheadermap="results/dbs/repdb/repdb_nohead.map",
        clusters="results/dbs/repdb/repdb_clusters.tsv",
        contaminants_tsv="results/dbs/repdb/decontaminate/contaminants.tsv",
        contaminants_ids="results/dbs/repdb/decontaminate/contaminants.txt",
        taxdump="results/taxdump/repdb_taxdump/",
    output:
        manifest="resources/releases/{version}/zenodo_stage/MANIFEST.txt",
    localrule: True
    run:
        import hashlib
        import os
        import tarfile

        stagedir = f"resources/releases/{wildcards.version}/zenodo_stage"
        os.makedirs(stagedir, exist_ok=True)

        def sha256(path):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            return h.hexdigest()

        # key -> (Zenodo filename, source path)
        assets = {
            "fasta": ("repdb.fa.gz", input.fa),
            "accession_map": ("repdb_accession_map.txt", input.idmap),
            "map": ("repdb.map", input.headermap),
            "nohead_map": ("repdb_nohead.map", input.noheadermap),
            "clusters": ("repdb_clusters.tsv", input.clusters),
            "contaminants_tsv": ("repdb_contaminants.tsv", input.contaminants_tsv),
            "contaminants_ids": ("repdb_contaminants.txt", input.contaminants_ids),
        }

        checksums = {}
        for key, (name, src) in assets.items():
            dest = f"{stagedir}/{name}"
            if os.path.lexists(dest):
                os.remove(dest)
            os.symlink(os.path.abspath(src), dest)
            checksums[key] = sha256(src)

        taxdump_name = "repdb_taxdump.tar.gz"
        taxdump_tar = f"{stagedir}/{taxdump_name}"
        with tarfile.open(taxdump_tar, "w:gz") as tf:
            tf.add(input.taxdump, arcname=".")
        checksums["taxdump"] = sha256(taxdump_tar)
        assets["taxdump"] = (taxdump_name, taxdump_tar)

        with open(output.manifest, "w") as mf:
            for key, (name, _) in assets.items():
                mf.write(f"{key}\t{name}\t{checksums[key]}\n")

        print(f"\n=== Zenodo assets staged in {stagedir}/ ===\n")
        print("1) create a new Zenodo deposit and upload every file in that "
              "directory (web UI, or the Zenodo REST API - there's no official "
              "CLI): https://zenodo.org/deposit/new\n")
        print("2) once published, paste this into "
              f"resources/releases/{wildcards.version}/config.yaml under "
              "dbs.build.repdb, filling in <record_id>:\n")
        print("  zenodo:")
        print('    base_url: "https://zenodo.org/records/<record_id>/files"')
        print("    files:")
        for key, (name, _) in assets.items():
            print(f'      {key}: {{name: {name}, sha256: "{checksums[key]}"}}')
        print()
