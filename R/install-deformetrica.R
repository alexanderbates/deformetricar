# Smooth, R-native setup of the Deformetrica (>= 4.3) command-line tool that this
# package wraps. Modelled on the `influencer` package's `install_python_*()` path:
# reticulate creates/uses a managed Python environment, pip-installs Deformetrica
# into it, and we then locate the `deformetrica` console script it drops in that
# environment's bin/ so `find_deformetrica()` resolves it with zero further config.

#' Install Deformetrica into a managed Python environment
#'
#' `deformetricar` shells out to the `deformetrica` (>= 4.3) command-line tool, so
#' you need that tool available once. This is the smooth path: it uses
#' [reticulate](https://rstudio.github.io/reticulate/) to create (or reuse) a
#' dedicated Python environment, `pip install`s Deformetrica into it, verifies the
#' `deformetrica` console script, and caches its path in
#' `options(deformetricar.exe=)` so every subsequent [find_deformetrica()] call
#' resolves it automatically (including in new sessions if you persist the option).
#'
#' @param envname Name of the conda env / virtualenv to install into
#'   (default `"deformetrica"`). A dedicated env is recommended because Deformetrica
#'   pins specific `torch`/`pykeops` versions.
#' @param method One of `"conda"` (default, recommended - Deformetrica's own install
#'   instructions use a conda env), `"virtualenv"`, or `"auto"` (conda if a conda
#'   binary is found, else virtualenv).
#' @param python_version Python for the new environment (default `"3.8"`, which
#'   matches the `torch==1.6` / `pykeops==1.4.1` that Deformetrica 4.3 pins).
#' @param version Package spec passed to pip (default `"deformetrica"`; pin with e.g.
#'   `"deformetrica==4.3.0"`).
#' @param extra_packages Extra packages to co-install (default `"numpy"`, per
#'   Deformetrica's install notes).
#' @param x86 Build an x86-64 (`osx-64`) conda env that runs under Rosetta. `NA`
#'   (default) auto-enables this on Apple Silicon and disables it elsewhere. Needed on
#'   arm64 macOS because Deformetrica's `torch==1.6` pin has no arm64 wheels; the
#'   x86-64 wheels install and run fine under Rosetta.
#' @param force Reinstall even if a `deformetrica` executable is already present in
#'   the environment.
#' @param ... Passed to [reticulate::conda_install()] / [reticulate::virtualenv_install()].
#' @return The path to the installed `deformetrica` executable, invisibly.
#' @section Apple Silicon (arm64 macOS): Deformetrica 4.3 pins `torch==1.6`, which has
#'   no native arm64 wheels, so a plain arm64 `pip install` fails. This function
#'   handles it for you: on an M-series Mac it defaults to `x86 = TRUE`, building an
#'   `osx-64` conda env (Python 3.8) whose interpreter runs under Rosetta and pulls the
#'   x86-64 wheels - verified end-to-end (estimate + compute) on this hardware. You
#'   need Rosetta (`softwareupdate --install-rosetta`) and a conda binary; the first
#'   `deformetrica` call per process pays a ~45 s Rosetta cold-start to import
#'   torch/vtk, then runs normally.
#' @seealso [find_deformetrica()] to locate an existing install.
#' @examples
#' \dontrun{
#' # One-time setup (needs conda/Miniconda; reticulate::install_miniconda() gets one):
#' install_deformetrica()
#' # ... then the rest of the package just works:
#' find_deformetrica()
#' }
#' @export
install_deformetrica <- function(envname = "deformetrica",
                                 method = c("conda", "virtualenv", "auto"),
                                 python_version = "3.8",
                                 version = "deformetrica",
                                 extra_packages = "numpy",
                                 x86 = NA,
                                 force = FALSE, ...) {
  method <- match.arg(method)
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("install_deformetrica() needs the 'reticulate' package. ",
         "install.packages('reticulate') first.", call. = FALSE)

  # On Apple Silicon the Deformetrica torch==1.6 pin has no arm64 wheels, so we build
  # an x86-64 (osx-64) env that runs under Rosetta and DOES get the x86 wheels. This
  # is verified working. Default to it automatically on arm64 macOS.
  is_arm_mac <- Sys.info()[["sysname"]] == "Darwin" && Sys.info()[["machine"]] == "arm64"
  if (is.na(x86)) x86 <- is_arm_mac
  if (isTRUE(x86)) {
    old_subdir <- Sys.getenv("CONDA_SUBDIR", unset = NA)
    Sys.setenv(CONDA_SUBDIR = "osx-64")
    on.exit(if (is.na(old_subdir)) Sys.unsetenv("CONDA_SUBDIR")
            else Sys.setenv(CONDA_SUBDIR = old_subdir), add = TRUE)
  }

  if (!force) {
    exe <- .deformetrica_env_exe(envname)
    if (!is.null(exe)) {
      message("Deformetrica is already installed in environment '", envname, "':\n  ", exe)
      options(deformetricar.exe = exe)
      return(invisible(exe))
    }
  }

  have_conda <- !inherits(try(reticulate::conda_binary(), silent = TRUE), "try-error")
  if (method == "auto") method <- if (have_conda) "conda" else "virtualenv"
  if (method == "conda" && !have_conda)
    stop("No conda binary found. Install Miniconda with ",
         "reticulate::install_miniconda(), or call install_deformetrica(method = 'virtualenv').",
         call. = FALSE)
  if (isTRUE(x86) && method != "conda")
    stop("x86 = TRUE (needed on Apple Silicon) requires method = 'conda'.", call. = FALSE)

  pkgs <- unique(c(extra_packages, version))
  message("Installing Deformetrica into the '", envname, "' ", method, " environment",
          if (isTRUE(x86)) " (x86-64 / Rosetta)" else "",
          " (this pulls torch/vtk and can take a few minutes)...")
  if (method == "conda") {
    if (!envname %in% reticulate::conda_list()$name)
      reticulate::conda_create(envname, python_version = python_version)
    # Pin the env itself to osx-64 so later conda ops in it stay x86 too.
    if (isTRUE(x86)) {
      cl <- reticulate::conda_list()
      prefix <- dirname(dirname(cl$python[match(envname, cl$name)]))
      if (!is.na(prefix) && dir.exists(prefix))
        try(writeLines("subdir: osx-64", file.path(prefix, ".condarc")), silent = TRUE)
    }
    # pip = TRUE: Deformetrica is only distributed on PyPI, not conda-forge.
    reticulate::conda_install(envname, packages = pkgs, pip = TRUE, ...)
  } else {
    if (!envname %in% reticulate::virtualenv_list())
      reticulate::virtualenv_create(envname, version = python_version)
    reticulate::virtualenv_install(envname, packages = pkgs, ...)
  }

  exe <- .deformetrica_env_exe(envname)
  if (is.null(exe))
    stop("Install completed but the 'deformetrica' console script was not found in ",
         "environment '", envname, "'. Check the pip output above for errors ",
         "(on arm64 macOS the torch==1.6 pin commonly fails - see ?install_deformetrica).",
         call. = FALSE)

  ok <- identical(as.integer(suppressWarnings(
    system2(exe, "--help", stdout = FALSE, stderr = FALSE))), 0L)
  options(deformetricar.exe = exe)
  if (ok)
    message("Deformetrica installed and located:\n  ", exe,
            "\nAdd  options(deformetricar.exe = '", exe, "')  to your .Rprofile to persist it.")
  else
    message("Deformetrica installed at:\n  ", exe,
            "\nbut `deformetrica --help` returned a non-zero status, so the Python ",
            "environment may be broken (check the pip output above).")
  invisible(exe)
}

# Path to the `deformetrica` console script inside a reticulate-managed conda env or
# virtualenv, or NULL if the env / executable is absent. reticulate is optional, so
# resolve silently to NULL when it (or the env) is missing.
.deformetrica_env_exe <- function(envname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  py <- suppressWarnings(tryCatch(reticulate::conda_python(envname), error = function(e) NA_character_))
  if (length(py) != 1L || is.na(py) || !file.exists(py))
    py <- suppressWarnings(tryCatch(reticulate::virtualenv_python(envname), error = function(e) NA_character_))
  if (length(py) != 1L || is.na(py) || !file.exists(py)) return(NULL)
  exe <- file.path(dirname(py),
                   if (.Platform$OS.type == "windows") "deformetrica.exe" else "deformetrica")
  if (file.exists(exe)) exe else NULL
}
