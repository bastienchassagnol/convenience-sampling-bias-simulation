# Convenience sampling bias simulation

Simulate the impact of convenience / range-restricted sampling on common effect-size estimators.

This repository is an **R analysis project** (not an R package): no `DESCRIPTION` / `NAMESPACE` / `devtools::check`.

## Quarto book

HTML and PDF editions are defined in `_quarto.yml` (technical chapter: `scripts/Simulation_study.qmd`).

```bash
Rscript -e 'renv::restore()'
quarto render
```

## GitHub Pages

On push to `main`, `.github/workflows/publish.yml`:

1. restores packages from `renv.lock`
2. runs `quarto render` (HTML + PDF into `_book/`)
3. deploys `_book/` to the **`gh-pages`** branch (with `.nojekyll`)

Keep **`main`** as the default branch. In GitHub:

**Settings → Pages → Build and deployment**

* Source: `Deploy from a branch`
* Branch: `gh-pages` / `(root)`

Also allow Actions write access: **Settings → Actions → General → Workflow permissions → Read and write**.

## Materials

* `scripts/Helper.R` — Monte Carlo, bias, and power helpers
* `scripts/Simulation_study.qmd` — confirmatory + exploratory DGPs
* `dependencies.R` — packages recorded by `renv`
* `biblio/` — instructions and manuscript notes
* `images/` — DAG assets
