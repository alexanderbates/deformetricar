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
  # Return the same class as the input: a warped mesh3d (faces preserved), a
  # neuron/dotprops with coordinates replaced, or a bare matrix.
  if (inherits(x, "mesh3d")) { x$vb <- rbind(t(out), 1); x$normals <- NULL; return(x) }
  if (inherits(x, c("neuron", "neuronlist", "dotprops"))) { nat::xyzmatrix(x) <- out; return(x) }
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

#' Fit a Deformetrica diffeomorphism between two point sets (registration)
#'
#' The Deformetrica 4 estimate path: fit a diffeomorphism that deforms `source`
#' (the moving template) onto `target` (the fixed subject). With ordered `source`
#' <-> `target` row correspondences a point-to-point Landmark attachment is used.
#' Returns the estimated control points + momenta, which [deformetrica_shoot()]
#' then applies to arbitrary new points.
#'
#' @param source,target N x 3 matrices of corresponding points (row *i* of
#'   `source` matches row *i* of `target`). Must have the same number of rows.
#' @param kernel_width Deformation kernel width (larger = stiffer / more global).
#' @param timepoints Geodesic integration steps.
#' @param object_type Deformetrica `deformable-object-type` (default "Landmark").
#' @param attachment_type "Landmark" (ordered L2, default), "Varifold" or "Current".
#' @param noise_std Data-attachment noise standard deviation.
#' @param max_iterations Optimiser iterations.
#' @param device `kernel-device`: "auto", "cpu" or "cuda".
#' @param deformetrica Optional path to the executable (auto-detected otherwise).
#' @param workdir Working directory for the run (a fresh tempdir by default).
#' @param verbose Show `deformetrica` output.
#' @return A list with `control_points` and `momenta` (file paths), `kernel_width`
#'   and `output_dir`; pass the paths straight to [deformetrica_shoot()].
#' @seealso [deformetrica_shoot()] to apply the fitted transform.
#' @export
deformetrica_register <- function(source, target, kernel_width,
                                  timepoints = 10L, object_type = "Landmark",
                                  attachment_type = c("Landmark", "Varifold", "Current"),
                                  noise_std = 0.01, max_iterations = 150L,
                                  device = c("auto", "cpu", "cuda"),
                                  deformetrica = NULL,
                                  workdir = tempfile("dfca_reg_"), verbose = FALSE) {
  attachment_type <- match.arg(attachment_type); device <- match.arg(device)
  src_is_mesh <- inherits(source, "mesh3d"); tgt_is_mesh <- inherits(target, "mesh3d")
  src <- nat::xyzmatrix(source); tgt <- nat::xyzmatrix(target)
  is_mesh <- src_is_mesh && tgt_is_mesh
  # Surfaces register unlabelled (Current/Varifold): no row correspondence needed
  # and the mesh FACES must go into the VTK. Point sets use ordered Landmark L2.
  if (is_mesh && object_type == "Landmark") object_type <- "SurfaceMesh"
  if (is_mesh && attachment_type == "Landmark") attachment_type <- "Current"
  if (!is_mesh && attachment_type == "Landmark" && nrow(src) != nrow(tgt))
    stop("For a Landmark fit source and target must have the same number of rows; ",
         "for unlabelled surfaces pass mesh3d objects (attachment_type='Current').", call. = FALSE)
  exe <- find_deformetrica(deformetrica)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  write.vtk(src, file.path(workdir, "source.vtk"), polygons = if (src_is_mesh) t(source$it) - 1L else NULL)
  write.vtk(tgt, file.path(workdir, "target.vtk"), polygons = if (tgt_is_mesh) t(target$it) - 1L else NULL)
  # Explicit MATCHING object id in both XMLs ("shape") — never derive it by
  # stripping a filename suffix (that silently desyncs when the suffix recurs).
  writeLines(sprintf(
'<?xml version="1.0" encoding="UTF-8"?>
<model>
  <model-type>Registration</model-type>
  <dimension>3</dimension>
  <template>
    <object id="shape">
      <deformable-object-type>%s</deformable-object-type>
      <attachment-type>%s</attachment-type>
      <noise-std>%g</noise-std>
      <kernel-width>%g</kernel-width>
      <kernel-type>torch</kernel-type>
      <kernel-device>%s</kernel-device>
      <filename>source.vtk</filename>
    </object>
  </template>
  <deformation-parameters>
    <kernel-width>%g</kernel-width>
    <kernel-type>torch</kernel-type>
    <number-of-timepoints>%d</number-of-timepoints>
  </deformation-parameters>
</model>',
    object_type, attachment_type, noise_std, kernel_width, device,
    kernel_width, as.integer(timepoints)), file.path(workdir, "model.xml"))
  writeLines(
'<?xml version="1.0" encoding="UTF-8"?>
<data-set>
  <subject id="subject">
    <visit id="visit">
      <filename object_id="shape">target.vtk</filename>
    </visit>
  </subject>
</data-set>', file.path(workdir, "data_set.xml"))
  writeLines(sprintf(
'<?xml version="1.0" encoding="UTF-8"?>
<optimization-parameters>
  <optimization-method-type>GradientAscent</optimization-method-type>
  <max-iterations>%d</max-iterations>
  <convergence-tolerance>1e-5</convergence-tolerance>
  <freeze-template>On</freeze-template>
</optimization-parameters>', as.integer(max_iterations)),
    file.path(workdir, "optimization_parameters.xml"))

  owd <- setwd(workdir); on.exit(setwd(owd), add = TRUE)
  status <- system2(exe, c("estimate", "model.xml", "data_set.xml", "-p",
                           "optimization_parameters.xml", "--output=output/"),
                    stdout = if (verbose) "" else FALSE,
                    stderr = if (verbose) "" else FALSE)
  if (!identical(as.integer(status), 0L))
    stop("`deformetrica estimate` failed (exit ", status, "). Re-run with verbose=TRUE.", call. = FALSE)
  od <- file.path(workdir, "output")
  cp  <- list.files(od, pattern = "ControlPoints\\.txt$", full.names = TRUE)
  mom <- list.files(od, pattern = "Momenta\\.txt$", full.names = TRUE)
  if (!length(cp) || !length(mom))
    stop("estimate produced no ControlPoints/Momenta in ", od, call. = FALSE)
  list(control_points = cp[[1]], momenta = mom[[1]],
       kernel_width = kernel_width, output_dir = od)
}

#' Write a neuron as a VTK NonOrientedPolyLine
#'
#' Neuron backbones are the natural registration object for connectome work
#' (cf. the FAFB left-right bridging registration in flyconnectome/deformetricaLR).
#'
#' @param x A `nat` neuron (or anything with `xyzmatrix` + a `SegList`).
#' @param file Output path.
#' @return "complete" (invisibly via the file written).
#' @export
write_neuron_vtk <- function(x, file) {
  pts <- nat::xyzmatrix(x)
  sl <- x$SegList
  if (is.null(sl)) sl <- nat::as.seglist(x)
  # Deformetrica 4 reads VTK LINES as 2-point segments ([2 i j], reshaped to (-1,3)),
  # so split each backbone path into consecutive edges rather than one long polyline.
  edges <- do.call(rbind, lapply(sl, function(s)
    if (length(s) >= 2L) cbind(s[-length(s)], s[-1L]) else NULL))
  ne <- nrow(edges)
  con <- file(file, "w"); on.exit(close(con))
  writeLines(c("# vtk DataFile Version 2.0", "deformetricar polyline", "ASCII",
               "DATASET POLYDATA", paste("POINTS", nrow(pts), "float")), con)
  utils::write.table(pts, con, row.names = FALSE, col.names = FALSE)
  writeLines(paste("LINES", ne, ne * 3L), con)
  utils::write.table(cbind(2L, edges[, 1] - 1L, edges[, 2] - 1L), con,
                     row.names = FALSE, col.names = FALSE)
  "complete"
}

# Classify an object and write it to VTK; returns its Deformetrica type + attachment.
.dfca_write_object <- function(x, file) {
  if (inherits(x, "mesh3d")) {
    write.vtk(nat::xyzmatrix(x), file, polygons = t(x$it) - 1L)   # 0-indexed faces for Deformetrica
    list(dtype = "SurfaceMesh", attach = "Current")
  } else if (inherits(x, c("neuron", "neuronlist"))) {
    write_neuron_vtk(x, file)
    list(dtype = "PolyLine", attach = "Varifold")   # Deformetrica 4.x name (was NonOrientedPolyLine in 2.1)
  } else {
    write.vtk(nat::xyzmatrix(x), file)
    list(dtype = "Landmark", attach = "Landmark")
  }
}

#' Fit ONE diffeomorphism to many matched objects (multi-object registration)
#'
#' The flyconnectome/deformetricaLR recipe, ported to Deformetrica 4: register a
#' whole SET of matched cognate objects (e.g. left/right neuron tracts) with a
#' single diffeomorphism, optionally anchored by a shared `landmarks` point cloud.
#' Objects may mix meshes, neurons (polylines) and point sets; each `sources[[i]]`
#' is matched to `targets[[i]]` by name.
#'
#' @param sources,targets Named, equal-length lists of matched objects. `sources`
#'   are templates (moving), `targets` the subject (fixed). Names become object ids.
#' @param kernel_width Deformation kernel width.
#' @param object_kernel_width Per-object data kernel width (defaults to `kernel_width`).
#' @param landmarks Optional `list(source=, target=)` of anchoring point matrices,
#'   added as a shared Landmark object (the `inc_tracts_lmarks` regulariser).
#' @param data_sigma Per-object data-attachment sigma.
#' @param timepoints,max_iterations,device,deformetrica,workdir,verbose As
#'   [deformetrica_register()].
#' @return A list with `control_points`, `momenta`, `kernel_width`, `output_dir`.
#' @seealso [deformetrica_register()] for the single-object case.
#' @export
deformetrica_register_multi <- function(sources, targets, kernel_width,
                                        object_kernel_width = kernel_width,
                                        landmarks = NULL, data_sigma = 0.5,
                                        timepoints = 10L, max_iterations = 150L,
                                        device = c("auto", "cpu", "cuda"),
                                        deformetrica = NULL,
                                        workdir = tempfile("dfca_multi_"), verbose = FALSE) {
  device <- match.arg(device)
  if (is.null(names(sources)) || length(sources) != length(targets))
    stop("sources and targets must be equal-length, named lists of matched objects.", call. = FALSE)
  ids <- make.names(names(sources), unique = TRUE)
  exe <- find_deformetrica(deformetrica)
  dir.create(file.path(workdir, "data"), recursive = TRUE, showWarnings = FALSE)

  tmpl <- character(); obj_xml <- character(); subj_xml <- character()
  for (i in seq_along(ids)) {
    sf <- file.path("data", paste0(ids[i], "_src.vtk"))
    tf <- file.path("data", paste0(ids[i], "_tgt.vtk"))
    spec <- .dfca_write_object(sources[[i]], file.path(workdir, sf))
    .dfca_write_object(targets[[i]], file.path(workdir, tf))
    obj_xml <- c(obj_xml, sprintf(
'    <object id="%s">
      <deformable-object-type>%s</deformable-object-type>
      <attachment-type>%s</attachment-type>
      <data-sigma>%g</data-sigma>
      <kernel-width>%g</kernel-width>
      <kernel-type>torch</kernel-type>
      <kernel-device>%s</kernel-device>
      <filename>%s</filename>
    </object>', ids[i], spec$dtype, spec$attach, data_sigma, object_kernel_width, device, sf))
    subj_xml <- c(subj_xml, sprintf('      <filename object_id="%s">%s</filename>', ids[i], tf))
  }
  if (!is.null(landmarks)) {
    write.vtk(nat::xyzmatrix(landmarks$source), file.path(workdir, "data/landmarks_src.vtk"))
    write.vtk(nat::xyzmatrix(landmarks$target), file.path(workdir, "data/landmarks_tgt.vtk"))
    obj_xml <- c(obj_xml, sprintf(
'    <object id="landmarks">
      <deformable-object-type>Landmark</deformable-object-type>
      <attachment-type>Landmark</attachment-type>
      <data-sigma>%g</data-sigma>
      <kernel-width>%g</kernel-width>
      <kernel-type>torch</kernel-type>
      <kernel-device>%s</kernel-device>
      <filename>data/landmarks_src.vtk</filename>
    </object>', data_sigma, object_kernel_width, device))
    subj_xml <- c(subj_xml, '      <filename object_id="landmarks">data/landmarks_tgt.vtk</filename>')
  }
  writeLines(c('<?xml version="1.0" encoding="UTF-8"?>', '<model>',
    '  <model-type>Registration</model-type>', '  <dimension>3</dimension>',
    '  <template>', obj_xml, '  </template>',
    '  <deformation-parameters>',
    sprintf('    <kernel-width>%g</kernel-width>', kernel_width),
    '    <kernel-type>torch</kernel-type>',
    sprintf('    <number-of-timepoints>%d</number-of-timepoints>', as.integer(timepoints)),
    '  </deformation-parameters>', '</model>'), file.path(workdir, "model.xml"))
  writeLines(c('<?xml version="1.0" encoding="UTF-8"?>', '<data-set>',
    '  <subject id="subject">', '    <visit id="visit">', subj_xml,
    '    </visit>', '  </subject>', '</data-set>'), file.path(workdir, "data_set.xml"))
  writeLines(sprintf(
'<?xml version="1.0" encoding="UTF-8"?>
<optimization-parameters>
  <optimization-method-type>GradientAscent</optimization-method-type>
  <max-iterations>%d</max-iterations>
  <convergence-tolerance>1e-5</convergence-tolerance>
  <freeze-template>On</freeze-template>
</optimization-parameters>', as.integer(max_iterations)),
    file.path(workdir, "optimization_parameters.xml"))

  owd <- setwd(workdir); on.exit(setwd(owd), add = TRUE)
  status <- system2(exe, c("estimate", "model.xml", "data_set.xml", "-p",
                           "optimization_parameters.xml", "--output=output/"),
                    stdout = if (verbose) "" else FALSE, stderr = if (verbose) "" else FALSE)
  if (!identical(as.integer(status), 0L))
    stop("`deformetrica estimate` (multi-object) failed (exit ", status, ").", call. = FALSE)
  od <- file.path(workdir, "output")
  cp  <- list.files(od, pattern = "ControlPoints\\.txt$", full.names = TRUE)
  mom <- list.files(od, pattern = "Momenta\\.txt$", full.names = TRUE)
  if (!length(cp) || !length(mom)) stop("multi-object estimate produced no CP/Momenta.", call. = FALSE)
  list(control_points = cp[[1]], momenta = mom[[1]], kernel_width = kernel_width, output_dir = od)
}
