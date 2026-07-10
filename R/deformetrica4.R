# Deformetrica >= 4 client. The historical shootflow() drives the Deformetrica 2.1
# `ShootAndFlow3` C++ binary (paramDiffeos.xml / CP_final.txt / Mom_final.txt),
# which no longer exists in Deformetrica 4.x. Deformetrica 4 instead exposes a
# single `deformetrica` CLI whose `compute` subcommand runs a <model-type>Shooting
# model that flows a template shape through a fitted diffeomorphism defined by
# initial control points + momenta. These functions wrap that path.

#' Locate the Deformetrica (>= 4.3) command-line executable
#'
#' Searches, in order: `options(deformetricar.exe=)`, the `deformetrica` on the
#' `PATH`, then the common conda location `~/.conda/envs/deformetrica/bin/deformetrica`.
#'
#' @param deformetrica Optional explicit path to (or name of) the executable.
#' @return The resolved executable path.
#' @export
find_deformetrica <- function(deformetrica = getOption("deformetricar.exe")) {
  if (!is.null(deformetrica)) {
    if (file.exists(deformetrica)) return(deformetrica)
    w <- Sys.which(deformetrica)
    if (nzchar(w)) return(unname(w))
  }
  w <- Sys.which("deformetrica")
  if (nzchar(w)) return(unname(w))
  cand <- path.expand("~/.conda/envs/deformetrica/bin/deformetrica")
  if (file.exists(cand)) return(cand)
  stop("Cannot find the 'deformetrica' (>= 4.3) executable. Install Deformetrica, ",
       "set options(deformetricar.exe = '/path/to/deformetrica'), or activate its ",
       "conda environment before starting R.", call. = FALSE)
}

#' Apply a fitted Deformetrica diffeomorphism to 3D points (geodesic shooting)
#'
#' The Deformetrica 4 replacement for the 2.1 `ShootAndFlow3` path in
#' [shootflow()]. Flows arbitrary points through the diffeomorphism defined by a
#' fitted registration's control points and initial momenta, by writing a
#' `<model-type>Shooting` model and calling `deformetrica compute`.
#'
#' @param x An N x 3 matrix of coordinates, or anything [nat::xyzmatrix()] accepts.
#' @param control_points,momenta Path to the estimated Deformetrica
#'   `...ControlPoints.txt` / `...Momenta.txt`, or a numeric matrix. Passing the
#'   files Deformetrica itself wrote is strongly preferred: it writes momenta with
#'   a leading `n_subjects n_cp dim` header that its own reader expects, so a bare
#'   matrix is best-effort only.
#' @param kernel_width Deformation kernel width, in the coordinate units. Larger =
#'   stiffer / more global; must match the width the registration was fit with.
#' @param timepoints Number of geodesic integration steps; must match the fit.
#' @param object_type Deformetrica `deformable-object-type` (default "Landmark").
#' @param device `kernel-device`: "auto", "cpu" or "cuda".
#' @param deformetrica Optional path to the executable (auto-detected otherwise).
#' @param workdir Working directory for the run (a fresh tempdir by default).
#' @param verbose Show `deformetrica` output.
#' @return An N x 3 matrix of deformed coordinates, in the input units and row order.
#' @seealso [shootflow()] for the legacy Deformetrica 2.1 path.
#' @export
deformetrica_shoot <- function(x, control_points, momenta, kernel_width,
                               timepoints = 10L, object_type = "Landmark",
                               device = c("auto", "cpu", "cuda"),
                               deformetrica = NULL,
                               workdir = tempfile("dfca_shoot_"),
                               verbose = FALSE) {
  device <- match.arg(device)
  pts <- nat::xyzmatrix(x)
  exe <- find_deformetrica(deformetrica)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  cp  <- .dfca_as_txt(control_points, file.path(workdir, "control_points.txt"))
  mom <- .dfca_as_txt(momenta,        file.path(workdir, "momenta.txt"))
  write.vtk(pts, file.path(workdir, "points.vtk"))
  .dfca_write_shooting_model(file.path(workdir, "model.xml"),
    points_file = "points.vtk",
    control_points = normalizePath(cp), momenta = normalizePath(mom),
    kernel_width = kernel_width, timepoints = timepoints,
    object_type = object_type, device = device)
  .dfca_write_min_optimization(file.path(workdir, "optimization_parameters.xml"))

  owd <- setwd(workdir); on.exit(setwd(owd), add = TRUE)
  status <- system2(exe, c("compute", "model.xml", "-p", "optimization_parameters.xml",
                           "--output=output/"),
                    stdout = if (verbose) "" else FALSE,
                    stderr = if (verbose) "" else FALSE)
  if (!identical(as.integer(status), 0L))
    stop("`deformetrica compute` failed (exit ", status, "). Re-run with verbose=TRUE.",
         call. = FALSE)

  out <- .dfca_read_final_flow(file.path(workdir, "output"))
  if (nrow(out) != nrow(pts))
    stop("Deformed point count (", nrow(out), ") != input (", nrow(pts), ").", call. = FALSE)
  unname(out)
}

# ---- internal helpers -------------------------------------------------------

.dfca_as_txt <- function(x, path) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) return(x)
  utils::write.table(as.matrix(x), path, row.names = FALSE, col.names = FALSE)
  path
}

.dfca_write_shooting_model <- function(file, points_file, control_points, momenta,
                                       kernel_width, timepoints, object_type, device) {
  writeLines(sprintf(
'<?xml version="1.0" encoding="UTF-8"?>
<model>
  <model-type>Shooting</model-type>
  <dimension>3</dimension>
  <template>
    <object id="shape">
      <deformable-object-type>%s</deformable-object-type>
      <kernel-width>%g</kernel-width>
      <kernel-type>torch</kernel-type>
      <kernel-device>%s</kernel-device>
      <filename>%s</filename>
    </object>
  </template>
  <initial-control-points>%s</initial-control-points>
  <initial-momenta>%s</initial-momenta>
  <deformation-parameters>
    <kernel-width>%g</kernel-width>
    <kernel-type>torch</kernel-type>
    <number-of-timepoints>%d</number-of-timepoints>
  </deformation-parameters>
</model>',
    object_type, kernel_width, device, points_file,
    control_points, momenta, kernel_width, as.integer(timepoints)), file)
}

.dfca_write_min_optimization <- function(file) {
  writeLines(
'<?xml version="1.0" encoding="UTF-8"?>
<optimization-parameters>
  <optimization-method-type>GradientAscent</optimization-method-type>
  <max-iterations>1</max-iterations>
</optimization-parameters>', file)
}

.dfca_read_final_flow <- function(outdir) {
  fs <- list.files(outdir, pattern = "GeodesicFlow.*\\.vtk$", full.names = TRUE)
  if (!length(fs)) fs <- list.files(outdir, pattern = "\\.vtk$", full.names = TRUE)
  if (!length(fs)) stop("No flow VTK produced in ", outdir, call. = FALSE)
  tp <- suppressWarnings(as.integer(sub(".*tp_([0-9]+).*", "\\1", basename(fs))))
  final <- fs[which.max(ifelse(is.na(tp), -1L, tp))]
  read.vtk(final, item = "points")
}
