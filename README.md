# kat-interactive

Interactive companion to the ELAIA-2 trial paper ([Nature Communications 2023; 14:2826](https://doi.org/10.1038/s41467-023-38532-3)). Readers can restrict cohort-level tables and figures to pre-specified subgroups while staying anchored to the published results.

## Structure

- `index.html` — single-page site, served as-is by GitHub Pages
- `css/styles.css`, `js/app.js` — presentation + interactivity (vanilla, no build step)
- `data/` — precomputed JSON for each table/figure (committed)
- `scripts/precompute.py` — regenerates `data/*.json` from the deidentified dataset

## Data

The deidentified trial dataset is **not** in this repository. To regenerate `data/*.json`, download it from [Dryad kh189327p](https://doi.org/10.5061/dryad.kh189327p) (ORCID login required) and place the CSV under `files/`. The column names expected by `scripts/precompute.py` match the Dryad deposit.

## Local development

```
python scripts/precompute.py     # regenerate data/*.json
python -m http.server 8000       # fetch() needs http://, not file://
```

Then open http://localhost:8000.

## Deploy

Push to `main`; enable GitHub Pages → "Deploy from a branch" → `main` / root.
