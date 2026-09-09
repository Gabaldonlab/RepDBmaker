# Reproducing the paper figures

Every figure in the manuscript is produced by three R Markdown notebooks under
`workflow/notebooks/`. 

All figures land in `results/plots/` (override with the `figdir` param).

## Step 1

```bash
Rscript -e 'root <- getwd()
  dir.create(file.path(root, "results/qc"), recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render("workflow/notebooks/eda_prepare.Rmd",
                    knit_root_dir = root,
                    output_dir    = file.path(root, "results/qc"),
                    output_file   = "eda_prepare.md",
                    run_pandoc    = FALSE)'
```

## Step 2


```bash
Rscript -e 'root <- getwd()
  dir.create(file.path(root, "results/qc"), recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render("workflow/notebooks/comparison.Rmd",
                    knit_root_dir = root,
                    output_dir    = file.path(root, "results/qc"),
                    output_file   = "comparison.html")'

Rscript -e 'root <- getwd()
  rmarkdown::render("workflow/notebooks/eda_report.Rmd",
                    knit_root_dir = root,
                    output_dir    = file.path(root, "results/qc"),
                    output_file   = "eda.html")'
```


## Known gaps

- **`FigS1.pdf` requires `results/plots/recall.RDS`, which nothing in this
  repository produces.** It comes from a sibling decontamination project. Copy
  it into `results/plots/` before rendering step 2, or Fig. S1 is skipped.
- The NR count tables are hard-coded to BSC paths in the `eda_prepare.Rmd`
  defaults. Override `nr_counts` / `clustnr_counts` / `ecoli_ids` elsewhere.
- The reduced GTDB backbone trees are built by the `obtain-reduced-prokatree`
  chunk and expected in `test/tree/`; without them the tree panels of Fig. 2
  and Fig. S2 are skipped.
