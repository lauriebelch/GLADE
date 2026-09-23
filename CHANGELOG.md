# Changelog

## v1.0.0 (2026)

First tagged release. Changes since GitHub commit 7f558f7:

- Uses ete4 (an OrthoFinder v3 dependency) instead of ete3, so no extra installs are needed in an OrthoFinder environment
- Supports OrthoFinder runs made with `-X`
- Species names with dots (e.g. `Rozella.allomycis`) now work (thanks to Alan Beavan for reporting this)
- Gene names with `(` `)` (which OrthoFinder changes to `_`) and `+` now work
- Gene trees are now read as trees rather than with a regular expression, so every leaf is converted. GLADE stops with a clear message if a leaf can't be matched, instead of failing later
- Species names starting with "n" are now converted correctly
- Gene-tree polytomies: all child clades are now used when finding duplications and losses after duplication (previously only the first two)
- Species-tree polytomies: GLADE now stops with a clear message (a rooted, fully bifurcating species tree is required, as in OrthoFinder) instead of silently giving too few losses
- New `--seed` option: ancestral gene sets are now reproducible (same input + seed = same output, whatever the number of threads)
- `extant_OG_counts.tsv` and `OrthogroupBranchChange.tsv` are now written to `GainsLossDuplication/` with species names (previously left in `WorkingDirectory/GladeWD/`); `*_bybranch.tsv` files now use species names
- Empty results tables (e.g. no losses after duplication) no longer crash GLADE
- If `Log.txt` is missing GLADE assumes a standard (non `-X`) OrthoFinder run instead of stopping
- New `--version` option, version shown when GLADE runs, and a `GLADE_run_info.txt` file for each run
- README: installation, input requirements, and a full description of every output file
- Added tests (`tests/`)
