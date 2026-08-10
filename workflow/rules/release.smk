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
