# Deformetrica >= 4 client (estimate + compute). Deformetrica 4 replaced the 2.1
# `ShootAndFlow3` C++ binary (paramDiffeos.xml / CP_final.txt / Mom_final.txt),
# which no longer exists in Deformetrica 4.x. Deformetrica 4 instead exposes a
# single `deformetrica` CLI whose `compute` subcommand runs a <model-type>Shooting
# model that flows a template shape through a fitted diffeomorphism defined by
# initial control points + momenta. These functions wrap that path.

#' Locate the Deformetrica (>= 4.3) command-line executable
#'
#' Searches, in order: `options(deformetricar.exe=)`, the `deformetrica` on the
#' `PATH`, the reticulate-managed environments [install_deformetrica()] creates
#' (`"deformetrica"` then `"r-reticulate"`), then the common conda location
#' `~/.conda/envs/deformetrica/bin/deformetrica`.
#'
#' @param deformetrica Optional explicit path to (or name of) the executable.
#' @return The resolved executable path.
#' @seealso [install_deformetrica()] to set one up from R.
#' @export
find_deformetrica <- function(deformetrica = getOption("deformetricar.exe")) {
  if (!is.null(deformetrica)) {
    if (file.exists(deformetrica)) return(deformetrica)
    w <- Sys.which(deformetrica)
    if (nzchar(w)) return(unname(w))
  }
  w <- Sys.which("deformetrica")
  if (nzchar(w)) return(unname(w))
  # Environments created by install_deformetrica() (reticulate optional -> NULL).
  for (env in c("deformetrica", "r-reticulate")) {
    exe <- .deformetrica_env_exe(env)
    if (!is.null(exe)) return(exe)
  }
  cand <- path.expand("~/.conda/envs/deformetrica/bin/deformetrica")
  if (file.exists(cand)) return(cand)
  stop("Cannot find the 'deformetrica' (>= 4.3) executable. Run install_deformetrica() ",
       "to set it up in a managed Python environment, set ",
       "options(deformetricar.exe = '/path/to/deformetrica'), or add it to your PATH.",
       call. = FALSE)
}

#' Apply a fitted Deformetrica diffeomorphism to 3D points (geodesic shooting)
#'
#' Flows arbitrary points through the diffeomorphism defined by a fitted
#' registration's control points and initial momenta, by writing a
#' `<model-type>Shooting` model and calling `deformetrica compute`.
#'
#' @param x An N x 3 matrix of coordinates, or anything [nat::xyzmatrix()] accepts.
#' @param control_points,momenta Path to the estimated Deformetrica
#'   `...ControlPoints.txt` / `...Momenta.txt`, or a numeric matrix. Passing the
#'   files Deformetrica itself wrote is strongly preferred: it writes momenta with
#'   a leading `n_subjects n_cp dim` header that its own reader expects, so a bare
#'   matrix is best-effort only. Alternatively pass a `deformetricareg` object (from
#'   [deformetrica_register()]) as `control_points` and leave `momenta`/`kernel_width`
#'   unset — they are taken from the registration.
#' @param kernel_width Deformation kernel width, in the coordinate units. Larger =
#'   stiffer / more global; must match the width the registration was fit with.
#' @param timepoints Number of geodesic integration steps; must match the fit.
#' @param object_type Deformetrica `deformable-object-type` (default "Landmark").
#' @param device `kernel-device`: "auto", "cpu" or "cuda".
#' @param deformetrica Optional path to the executable (auto-detected otherwise).
#' @param workdir Working directory for the run (a fresh tempdir by default).
#' @param flow If `TRUE`, return the *whole* geodesic flow (every timepoint) as a
#'   list of deformed objects (one per timepoint, named `tp_00`..`tp_NN`) rather
#'   than just the final shape. Useful for animating a warp.
#' @param verbose Show `deformetrica` output.
#' @return By default an N x 3 matrix of deformed coordinates (or a warped object
#'   of the input class), in the input units and row order. With `flow = TRUE`, a
#'   list of such objects, one per geodesic timepoint.
#' @export
deformetrica_shoot <- function(x, control_points, momenta = NULL, kernel_width = NULL,
                               timepoints = 10L, object_type = "Landmark",
                               device = c("auto", "cpu", "cuda"),
                               deformetrica = NULL,
                               workdir = tempfile("dfca_shoot_"),
                               flow = FALSE,
                               verbose = FALSE) {
  device <- match.arg(device)
  # A fitted registration can be passed directly as `control_points`: unpack its
  # momenta / kernel width / timepoints (rewriting its files if they were cleaned).
  if (inherits(control_points, "deformetricareg")) {
    reg <- control_points; ff <- .dfca_reg_files(reg)
    control_points <- ff$cp
    if (is.null(momenta)) momenta <- ff$mom
    if (is.null(kernel_width)) kernel_width <- reg$kernel_width
    if (missing(timepoints) && !is.null(reg$timepoints)) timepoints <- reg$timepoints
  }
  if (is.null(momenta) || is.null(kernel_width))
    stop("Supply `momenta` and `kernel_width`, or pass a deformetricareg as `control_points`.",
         call. = FALSE)
  pts <- nat::xyzmatrix(x)
  exe <- find_deformetrica(deformetrica)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  cp  <- .dfca_as_txt(control_points, file.path(workdir, "control_points.txt"))
  mom <- .dfca_as_txt(momenta,        file.path(workdir, "momenta.txt"))
  # A single control point makes Deformetrica's torch kernel convolve a 1-D array and
  # die with a cryptic "Dimension out of range" IndexError. Catch it here with an
  # actionable message: it means the fit's kernel_width exceeded the object's extent
  # (control points sit on a grid spaced by kernel_width), so too few were laid down.
  if (nrow(utils::read.table(cp)) < 2L)
    stop("The fitted diffeomorphism has fewer than 2 control points, which ",
         "Deformetrica cannot shoot. Refit with a smaller `kernel_width` (it must be ",
         "smaller than the object's spatial extent so more than one control point is ",
         "placed).", call. = FALSE)
  write_vtk(pts, file.path(workdir, "points.vtk"))
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

  # Re-cast a deformed coordinate matrix back to the input's class.
  as_input <- function(out) {
    if (nrow(out) != nrow(pts))
      stop("Deformed point count (", nrow(out), ") != input (", nrow(pts), ").", call. = FALSE)
    if (inherits(x, "mesh3d")) { y <- x; y$vb <- rbind(t(out), 1); y$normals <- NULL; return(y) }
    if (inherits(x, c("neuron", "neuronlist", "dotprops"))) { y <- x; nat::xyzmatrix(y) <- out; return(y) }
    unname(out)
  }
  if (flow) {
    outs <- .dfca_read_all_flow(file.path(workdir, "output"))
    res <- lapply(outs, as_input)
    names(res) <- sprintf("tp_%02d", seq_along(res) - 1L)
    return(res)
  }
  as_input(.dfca_read_final_flow(file.path(workdir, "output")))
}

# Build a portable `deformetricareg` from a fit's control-point / momenta files. It keeps
# BOTH the file paths (fast, this session) and their verbatim contents, so a saveRDS'd
# registration survives the tempdir being cleaned or being moved to another machine.
.dfca_reg <- function(control_points, momenta, kernel_width, timepoints, output_dir = NULL) {
  structure(list(control_points = control_points, momenta = momenta,
                 kernel_width = kernel_width, timepoints = as.integer(timepoints),
                 output_dir = output_dir,
                 cp_text = readLines(control_points), mom_text = readLines(momenta)),
            class = "deformetricareg")
}

# Resolve a reg's control-point / momenta to readable files, rewriting from the stored
# contents when the original tempfiles are gone (e.g. after saveRDS + reload).
.dfca_reg_files <- function(reg) {
  cp <- reg$control_points; mom <- reg$momenta
  if (is.null(cp) || !file.exists(cp)) { cp <- tempfile("dfca_cp_", fileext = ".txt")
    writeLines(reg$cp_text, cp) }
  if (is.null(mom) || !file.exists(mom)) { mom <- tempfile("dfca_mom_", fileext = ".txt")
    writeLines(reg$mom_text, mom) }
  list(cp = cp, mom = mom)
}

#' Apply a Deformetrica registration as a natverse transform
#'
#' [deformetrica_register()] and [deformetrica_register_multi()] return a
#' `deformetricareg` — a portable handle on the fitted diffeomorphism, and a
#' [nat::xformpoints()] method. That makes a fitted registration a first-class natverse
#' transform: `nat::xform(x, reg)` warps any points / neuron / neuronlist / mesh through
#' it, and it composes with other registrations (e.g. in a `nat.templatebrains` reglist).
#' The control points and momenta are stored *inline*, so a `deformetricareg` survives
#' `saveRDS()` and being moved to another machine.
#'
#' @param reg A `deformetricareg`.
#' @param points An N x 3 matrix of coordinates.
#' @param ... Passed to [deformetrica_shoot()] (e.g. `object_type`, `device`).
#' @return An N x 3 matrix of deformed coordinates.
#' @seealso [deformetrica_shoot()], which does the actual geodesic shooting.
#' @exportS3Method nat::xformpoints
xformpoints.deformetricareg <- function(reg, points, ...) {
  # nat::xform passes control args (FallBackToAffine, na.action, ...) that shoot does
  # not take; keep only those deformetrica_shoot() understands.
  dots <- list(...)
  dots <- dots[names(dots) %in% names(formals(deformetrica_shoot))]
  do.call(deformetrica_shoot, c(list(points, reg), dots))
}

#' @export
print.deformetricareg <- function(x, ...) {
  cat("<deformetricareg> a fitted Deformetrica diffeomorphism\n")
  cat(sprintf("  ~%d control points | kernel width %g | %d timepoints\n",
              length(x$cp_text), x$kernel_width, x$timepoints))
  cat("  apply with nat::xform(object, reg) or deformetrica_shoot(object, reg)\n")
  invisible(x)
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

.dfca_flow_files <- function(outdir) {
  fs <- list.files(outdir, pattern = "GeodesicFlow.*\\.vtk$", full.names = TRUE)
  if (!length(fs)) fs <- list.files(outdir, pattern = "\\.vtk$", full.names = TRUE)
  if (!length(fs)) stop("No flow VTK produced in ", outdir, call. = FALSE)
  tp <- suppressWarnings(as.integer(sub(".*tp_([0-9]+).*", "\\1", basename(fs))))
  fs[order(ifelse(is.na(tp), 0L, tp))]
}

.dfca_read_final_flow <- function(outdir) {
  fs <- .dfca_flow_files(outdir)
  read_vtk(fs[[length(fs)]], item = "points")
}

.dfca_read_all_flow <- function(outdir) {
  lapply(.dfca_flow_files(outdir), read_vtk, item = "points")
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
#' @return A `deformetricareg` object (a portable handle on the fitted diffeomorphism,
#'   storing the control points + momenta inline). Apply it with `nat::xform(x, reg)` or
#'   [deformetrica_shoot()]; it survives `saveRDS()`. Its elements `control_points`,
#'   `momenta`, `kernel_width`, `timepoints` and `output_dir` are still accessible.
#' @seealso [xformpoints.deformetricareg()] / [deformetrica_shoot()] to apply the fit.
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
  write_vtk(src, file.path(workdir, "source.vtk"), polygons = if (src_is_mesh) t(source$it) - 1L else NULL)
  write_vtk(tgt, file.path(workdir, "target.vtk"), polygons = if (tgt_is_mesh) t(target$it) - 1L else NULL)
  # Explicit MATCHING object id in both XMLs ("shape") - never derive it by
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
  .dfca_reg(cp[[1]], mom[[1]], kernel_width, timepoints, od)
}

#' Write a neuron as a VTK NonOrientedPolyLine
#'
#' Neuron backbones are the natural registration object for connectome work
#' (cf. the FAFB left-right bridging registration in flyconnectome/deformetricaLR).
#'
#' @param x A `nat` neuron, or a `neuronlist` (all backbones are written into one
#'   VTK, so a whole set of neurons can be registered as a single PolyLine object).
#' @param file Output path.
#' @return "complete" (invisibly via the file written).
#' @export
write_neuron_vtk <- function(x, file) {
  # Points + consecutive-vertex edges for one neuron.
  one <- function(n) {
    sl <- n$SegList; if (is.null(sl)) sl <- nat::as.seglist(n)
    list(pts = nat::xyzmatrix(n),
         edges = do.call(rbind, lapply(sl, function(s)
           if (length(s) >= 2L) cbind(s[-length(s)], s[-1L]) else NULL)))
  }
  if (inherits(x, "neuronlist")) {
    # Concatenate every neuron, offsetting each one's vertex indices.
    pts <- list(); edges <- list(); off <- 0L
    for (n in x) {
      p <- one(n); pts[[length(pts) + 1L]] <- p$pts
      if (!is.null(p$edges)) edges[[length(edges) + 1L]] <- p$edges + off
      off <- off + nrow(p$pts)
    }
    pts <- do.call(rbind, pts); edges <- do.call(rbind, edges)
  } else {
    p <- one(x); pts <- p$pts; edges <- p$edges
  }
  # Deformetrica 4 reads VTK LINES as 2-point segments ([2 i j], reshaped to (-1,3)),
  # so split each backbone path into consecutive edges rather than one long polyline.
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
    write_vtk(nat::xyzmatrix(x), file, polygons = t(x$it) - 1L)   # 0-indexed faces for Deformetrica
    list(dtype = "SurfaceMesh", attach = "Current")
  } else if (inherits(x, c("neuron", "neuronlist"))) {
    write_neuron_vtk(x, file)
    list(dtype = "PolyLine", attach = "Varifold")   # Deformetrica 4.x name (was NonOrientedPolyLine in 2.1)
  } else {
    write_vtk(nat::xyzmatrix(x), file)
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
#' @param object_kernel_width Data-attachment kernel width - the spatial scale at which
#'   each object's surface/curve mismatch is measured (Current/Varifold). *Smaller =
#'   finer, more local matching*, which pulls small structures (galls, noduli) onto
#'   their target more tightly. A single value for every object, or a per-object vector
#'   (recycled if length 1; matched by name to `sources` if named, else in order).
#'   Defaults to `kernel_width`.
#' @param landmarks Optional `list(source=, target=)` of anchoring point matrices,
#'   added as a shared Landmark object (the `inc_tracts_lmarks` regulariser).
#' @param data_sigma Data-attachment noise sigma. A single value applied to every
#'   object, or a per-object vector (recycled if length 1; matched by name to
#'   `sources` if named, else taken in order). *Smaller sigma weights that object
#'   more strongly* in the fit - e.g. give homologous neuropils a small sigma and
#'   the outer hull a larger one so the neuropils drive the internal alignment.
#' @param landmark_sigma Data-attachment sigma for the optional shared `landmarks`
#'   object (defaults to the mean of `data_sigma`).
#' @param timepoints,max_iterations,device,deformetrica,workdir,verbose As
#'   [deformetrica_register()].
#' @return A `deformetricareg` object; see [deformetrica_register()]. Apply it with
#'   `nat::xform(x, reg)` or [deformetrica_shoot()].
#' @seealso [deformetrica_register()] for the single-object case.
#' @export
deformetrica_register_multi <- function(sources, targets, kernel_width,
                                        object_kernel_width = kernel_width,
                                        landmarks = NULL, data_sigma = 0.5,
                                        landmark_sigma = mean(data_sigma),
                                        timepoints = 10L, max_iterations = 150L,
                                        device = c("auto", "cpu", "cuda"),
                                        deformetrica = NULL,
                                        workdir = tempfile("dfca_multi_"), verbose = FALSE) {
  device <- match.arg(device)
  if (is.null(names(sources)) || length(sources) != length(targets))
    stop("sources and targets must be equal-length, named lists of matched objects.", call. = FALSE)
  ids <- make.names(names(sources), unique = TRUE)
  # Resolve data_sigma to one value per object (recycle scalar, match names, or in order).
  sigma <- if (length(data_sigma) == 1L) rep(data_sigma, length(sources))
    else if (!is.null(names(data_sigma))) unname(data_sigma[names(sources)])
    else if (length(data_sigma) == length(sources)) data_sigma
    else stop("data_sigma must be length 1, length(sources), or named by sources.", call. = FALSE)
  if (anyNA(sigma)) stop("named data_sigma is missing an entry for some source object.", call. = FALSE)
  # Resolve object_kernel_width per object the same way (scalar, named, or in order).
  okw <- if (length(object_kernel_width) == 1L) rep(object_kernel_width, length(sources))
    else if (!is.null(names(object_kernel_width))) unname(object_kernel_width[names(sources)])
    else if (length(object_kernel_width) == length(sources)) object_kernel_width
    else stop("object_kernel_width must be length 1, length(sources), or named by sources.", call. = FALSE)
  if (anyNA(okw)) stop("named object_kernel_width is missing an entry for some source object.", call. = FALSE)
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
      <noise-std>%g</noise-std>
      <kernel-width>%g</kernel-width>
      <kernel-type>torch</kernel-type>
      <kernel-device>%s</kernel-device>
      <filename>%s</filename>
    </object>', ids[i], spec$dtype, spec$attach, sigma[i], okw[i], device, sf))
    subj_xml <- c(subj_xml, sprintf('      <filename object_id="%s">%s</filename>', ids[i], tf))
  }
  if (!is.null(landmarks)) {
    write_vtk(nat::xyzmatrix(landmarks$source), file.path(workdir, "data/landmarks_src.vtk"))
    write_vtk(nat::xyzmatrix(landmarks$target), file.path(workdir, "data/landmarks_tgt.vtk"))
    obj_xml <- c(obj_xml, sprintf(
'    <object id="landmarks">
      <deformable-object-type>Landmark</deformable-object-type>
      <attachment-type>Landmark</attachment-type>
      <noise-std>%g</noise-std>
      <kernel-width>%g</kernel-width>
      <kernel-type>torch</kernel-type>
      <kernel-device>%s</kernel-device>
      <filename>data/landmarks_src.vtk</filename>
    </object>', landmark_sigma, object_kernel_width, device))
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
  .dfca_reg(cp[[1]], mom[[1]], kernel_width, timepoints, od)
}
