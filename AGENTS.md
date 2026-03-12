# Repository Guidelines

## Project Structure & Module Organization
`scICER` is an R package with standard package layout:
- `R/`: core source code (clustering pipeline, ECS/IC utilities, visualization, installation helpers).
- `tests/testthat/`: unit tests (`test-*.R`) and snapshots in `_snaps/`.
- `man/`: generated `.Rd` documentation.
- `vignettes/`: long-form usage docs (`scICER-quickstart.Rmd`).
- `examples/`: runnable examples.
- Root metadata: `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `README.md`.

Prefer adding new user-facing functions in `R/*.R` with roxygen comments, then regenerate docs.

## Build, Test, and Development Commands
Run commands from repository root:

```bash
Rscript install.R
```
Installs package locally from source and runs dependency setup helper.

```bash
Rscript install_and_test.R
```
Performs dependency checks, `devtools::document()`, `devtools::check()`, install, and a quick functional smoke test.

```bash
R -q -e "devtools::test()"
R -q -e "devtools::check()"
```
Use during development for unit tests and full package checks.

## Coding Style & Naming Conventions
- Use idiomatic R style: 2-space indentation, clear argument names, and explicit input validation (`stop(...)`).
- Prefer `snake_case` for internal helpers and test names; preserve existing exported API names (for example `scICE_clustering`).
- Keep files focused by concern (for example clustering logic in `R/clustering_core.R`, plotting in `R/visualization.R`).
- Document exported functions with roxygen2; do not hand-edit `man/*.Rd`.

## Testing Guidelines
- Framework: `testthat` (see `tests/testthat.R`).
- Add tests under `tests/testthat/test-<feature>.R`.
- Cover normal behavior, edge cases, and error paths (`expect_error`, bounds checks, reproducibility-sensitive paths).
- For algorithm/performance changes, add a small deterministic fixture (set `seed`) to keep tests stable.

## Commit & Pull Request Guidelines
Current history favors short, imperative commit messages such as `Update README.md`, `fix namespace problem`, or `Optimize nested parallel worker allocation in clustering`.

- Keep commit subjects concise and action-oriented; one logical change per commit.
- PRs should include: purpose, key code paths changed, test evidence (`devtools::test()`/`devtools::check()`), and any user-visible behavior changes.
- Link related issues and include updated examples/plots when visualization or output interpretation changes.
