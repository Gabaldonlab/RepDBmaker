# RepDB releases

Each `resources/releases/<db_version>/` freezes one reproducible RepDB release.
The **software** version is in the top-level `VERSION` file (currently `0.0.1`,
pre-release); the **database** version is the directory name (`v1`, `v2`, …). They
are tracked separately — a release records which software version built it.

## Layout

Tracked (light, in git):

- `config.yaml` — a self-contained config that pins the universe + custom bundle;
  reproduce with `snakemake build --configfile resources/releases/<v>/config.yaml`.
- `repdb.ids` — the selected composition (human-diffable changelog across releases).
- `release.yaml` — the manifest: software version, git commit, and each heavy
  asset's `sha256` + download URL.

Not tracked (heavy, gitignored → GitHub Release assets, later Zenodo):

- `universe.tsv` — the frozen taxonomy + completeness for every proteome.
- `custom_bundle.tar.gz` — the release custom proteomes + BUSCO.

## Cutting a release

After a full run (so `results/universe/` exists):

```bash
snakemake resources/releases/v1/release.yaml     # stages the release
gh release create v1 resources/releases/v1/universe.tsv \
   resources/releases/v1/custom_bundle.tar.gz \
   --title 'RepDB v1' --generate-notes           # command is printed by the rule
git add resources/releases/v1/{config.yaml,repdb.ids,release.yaml} && git commit
```

Tag the **software** too (semver), and — once enabled — the Zenodo–GitHub
integration mints a DOI per tag:

```bash
git tag -a 0.0.1 -m "RepDBmaker 0.0.1"; git push origin 0.0.1
```

## Reproducing a release

```bash
git checkout <software git_commit from release.yaml>      # or: docker run <digest>
# download the assets named in release.yaml into resources/releases/<v>/ (verify sha256)
snakemake build --configfile resources/releases/<v>/config.yaml
```

The pinned universe skips harmonization and keeps the taxdump full; the selection
re-runs deterministically; release custom is fetched from the bundle — no local
files, no drifting sources.
