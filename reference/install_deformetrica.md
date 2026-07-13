# Install Deformetrica into a managed Python environment

`deformetricar` shells out to the `deformetrica` (\>= 4.3) command-line
tool, so you need that tool available once. This is the smooth path: it
uses [reticulate](https://rstudio.github.io/reticulate/) to create (or
reuse) a dedicated Python environment, `pip install`s Deformetrica into
it, verifies the `deformetrica` console script, and caches its path in
`options(deformetricar.exe=)` so every subsequent
[`find_deformetrica()`](https://natverse.github.io/deformetricar/reference/find_deformetrica.md)
call resolves it automatically (including in new sessions if you persist
the option).

## Usage

``` r
install_deformetrica(
  envname = "deformetrica",
  method = c("conda", "virtualenv", "auto"),
  python_version = "3.8",
  version = "deformetrica",
  extra_packages = "numpy",
  x86 = NA,
  force = FALSE,
  ...
)
```

## Arguments

- envname:

  Name of the conda env / virtualenv to install into (default
  `"deformetrica"`). A dedicated env is recommended because Deformetrica
  pins specific `torch`/`pykeops` versions.

- method:

  One of `"conda"` (default, recommended - Deformetrica's own install
  instructions use a conda env), `"virtualenv"`, or `"auto"` (conda if a
  conda binary is found, else virtualenv).

- python_version:

  Python for the new environment (default `"3.8"`, which matches the
  `torch==1.6` / `pykeops==1.4.1` that Deformetrica 4.3 pins).

- version:

  Package spec passed to pip (default `"deformetrica"`; pin with e.g.
  `"deformetrica==4.3.0"`).

- extra_packages:

  Extra packages to co-install (default `"numpy"`, per Deformetrica's
  install notes).

- x86:

  Build an x86-64 (`osx-64`) conda env that runs under Rosetta. `NA`
  (default) auto-enables this on Apple Silicon and disables it
  elsewhere. Needed on arm64 macOS because Deformetrica's `torch==1.6`
  pin has no arm64 wheels; the x86-64 wheels install and run fine under
  Rosetta.

- force:

  Reinstall even if a `deformetrica` executable is already present in
  the environment.

- ...:

  Passed to
  [`reticulate::conda_install()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  /
  [`reticulate::virtualenv_install()`](https://rstudio.github.io/reticulate/reference/virtualenv-tools.html).

## Value

The path to the installed `deformetrica` executable, invisibly.

## Apple Silicon (arm64 macOS)

Deformetrica 4.3 pins `torch==1.6`, which has no native arm64 wheels, so
a plain arm64 `pip install` fails. This function handles it for you: on
an M-series Mac it defaults to `x86 = TRUE`, building an `osx-64` conda
env (Python 3.8) whose interpreter runs under Rosetta and pulls the
x86-64 wheels - verified end-to-end (estimate + compute) on this
hardware. You need Rosetta (`softwareupdate --install-rosetta`) and a
conda binary; the first `deformetrica` call per process pays a ~45 s
Rosetta cold-start to import torch/vtk, then runs normally.

## See also

[`find_deformetrica()`](https://natverse.github.io/deformetricar/reference/find_deformetrica.md)
to locate an existing install.

## Examples

``` r
if (FALSE) { # \dontrun{
# One-time setup (needs conda/Miniconda; reticulate::install_miniconda() gets one):
install_deformetrica()
# ... then the rest of the package just works:
find_deformetrica()
} # }
```
