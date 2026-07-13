# Helpers that make a whole cross-brain warp workflow expressible in plain R:
# an affine pre-alignment, left/right splitting of neuropil surfaces (by a plane
# or, more carefully, by a self-mirror registration), and a nat.ggplot animation
# of a geodesic flow. Used by the package vignettes.

#' Affine pre-alignment of one object onto another
#'
#' A diffeomorphism should not have to absorb a gross change of pose and scale;
#' that is this step's job. Both objects are centred, `source` is isotropically
#' scaled to `target`'s bounding-box diagonal, and a rigid (or affine) ICP refines
#' the orientation. The returned `apply()` closure replays the *same* transform on
#' any other object, so a whole set of sub-objects (e.g. neuropils) can be carried
#' into the aligned frame consistently.
#'
#' @param source,target Objects accepted by [nat::xyzmatrix()] (`mesh3d`, neuron,
#'   dotprops or an N x 3 matrix). `source` moves onto `target`.
#' @param type ICP type passed to [Morpho::icpmat()]: "rigid", "similarity" or "affine".
#' @param iterations ICP iterations.
#' @return A list with `aligned` (the transformed `source`), `apply` (a function
#'   mapping any object through the identical transform) and the components
#'   `scale`, `centre_source`, `centre_target` and `icp` (the 4x4 ICP matrix).
#' @seealso [deformetrica_register()] for the non-rigid step that follows.
#' @export
affine_prealign <- function(source, target, type = c("rigid", "similarity", "affine"),
                            iterations = 30L) {
  type <- match.arg(type)
  if (!requireNamespace("Morpho", quietly = TRUE))
    stop("affine_prealign() needs the 'Morpho' package.", call. = FALSE)
  Vs <- nat::xyzmatrix(source); Vt <- nat::xyzmatrix(target)
  diag_len <- function(V) sqrt(sum((apply(V, 2, max) - apply(V, 2, min))^2))
  s  <- diag_len(Vt) / diag_len(Vs)
  cs <- colMeans(Vs); ct <- colMeans(Vt)
  prescale <- function(P) sweep(sweep(P, 2, cs) * s, 2, -ct)   # centre, scale, recentre on target
  Ps <- prescale(Vs)
  # ICP the pre-scaled source onto the target, then recover the 4x4 it applied so
  # the same transform can be replayed on other objects (Morpho::icp is not exported;
  # icpmat + computeTransform is the supported path).
  aligned <- Morpho::icpmat(Ps, Vt, iterations = iterations, type = type)
  Tic <- Morpho::computeTransform(aligned, Ps, type = type)
  apply_fun <- function(x) {
    P <- Morpho::applyTransform(prescale(nat::xyzmatrix(x)), Tic)
    nat::xyzmatrix(x) <- P
    x
  }
  list(aligned = apply_fun(source), apply = apply_fun,
       scale = s, centre_source = cs, centre_target = ct, icp = Tic)
}

# Build a sub-mesh from a subset of face indices, dropping now-unused vertices.
.submesh_faces <- function(mesh, keep) {
  it <- mesh$it[, keep, drop = FALSE]
  used <- sort(unique(as.vector(it)))
  remap <- integer(max(used)); remap[used] <- seq_along(used)
  out <- mesh
  out$vb <- mesh$vb[, used, drop = FALSE]
  out$it <- matrix(remap[it], nrow = 3)
  out$normals <- NULL
  out
}

#' Split a surface mesh into left and right halves by a midline plane
#'
#' Faces are assigned to a side by their centroid position along `axis` relative
#' to `mid`. For a cleaner split of a shape whose symmetry surface is curved
#' rather than a flat plane, see [mirror_lr_split()].
#'
#' @param x A `mesh3d`.
#' @param mid The midline coordinate on `axis`; defaults to the midpoint of the
#'   mesh's extent on that axis.
#' @param axis Split axis, "X" (default), "Y" or "Z".
#' @param left A face on the side with the *smaller* (`"low"`, default) or *larger*
#'   (`"high"`) `axis` coordinate is called the left half.
#' @return A list with `L` and `R` sub-meshes and the `mid` used.
#' @export
split_mesh_lr <- function(x, mid = NULL, axis = c("X", "Y", "Z"),
                          left = c("low", "high")) {
  axis <- match.arg(axis); left <- match.arg(left)
  if (!inherits(x, "mesh3d")) stop("split_mesh_lr() expects a mesh3d.", call. = FALSE)
  ax <- match(axis, c("X", "Y", "Z"))
  V <- t(x$vb[1:3, , drop = FALSE])
  if (is.null(mid)) mid <- mean(range(V[, ax]))
  fc <- colMeans(matrix(V[x$it, ax], nrow = 3))          # per-face centroid coord
  low <- fc < mid
  Lkeep <- if (left == "low") low else !low
  list(L = .submesh_faces(x, which(Lkeep)),
       R = .submesh_faces(x, which(!Lkeep)), mid = mid)
}

#' Assign left/right to an (even unpaired, midline) object by self-mirror warping
#'
#' A mini use-case of the package. An object that straddles the midline (a fly
#' ellipsoid body, protocerebral bridge, ...) has no left/right label, yet every
#' point still belongs to a hemisphere. We recover that the way a mirroring
#' registration does (cf. `bancr::banc_lr`): reflect the object across a nominal
#' midline plane and, optionally, `refine` the reflection with a Deformetrica warp
#' of the object onto its own reflection so the true (possibly curved) symmetry
#' surface is found. Each point is then labelled by *how far and in which direction
#' it must travel to reach its mirror partner* - points on the symmetry surface
#' barely move; points to one side travel toward the other.
#'
#' @param x A `mesh3d`.
#' @param axis Mirror axis, "X" (default), "Y" or "Z".
#' @param mid Midline coordinate on `axis` (defaults to the extent midpoint).
#' @param refine If `TRUE`, refine the flat reflection with a Deformetrica warp of
#'   `x` onto its reflection (finds a curved symmetry surface). Needs a working
#'   Deformetrica install; passes `...` (e.g. `kernel_width`, `device`) through to
#'   [deformetrica_register()]/[deformetrica_shoot()].
#' @param band Fraction of the axis extent within which |travel| is treated as
#'   midline (`"M"`) rather than L/R. Default 0 (every point is L or R).
#' @param positive_travel Which hemisphere a point that travels in the `+axis`
#'   direction belongs to ("L" default, or "R").
#' @param ... Passed to the Deformetrica calls when `refine = TRUE`.
#' @return A list with `L`/`R` sub-meshes, a per-vertex `side` factor and the
#'   signed `travel` along `axis`.
#' @seealso [split_mesh_lr()] for the plain-plane split.
#' @export
mirror_lr_split <- function(x, axis = c("X", "Y", "Z"), mid = NULL, refine = FALSE,
                            band = 0, positive_travel = c("L", "R"), ...) {
  axis <- match.arg(axis); positive_travel <- match.arg(positive_travel)
  if (!inherits(x, "mesh3d")) stop("mirror_lr_split() expects a mesh3d.", call. = FALSE)
  ax <- match(axis, c("X", "Y", "Z"))
  V <- t(x$vb[1:3, , drop = FALSE])
  if (is.null(mid)) mid <- mean(range(V[, ax]))
  # Reflection of x across the midline plane.
  refl <- x; Vr <- V; Vr[, ax] <- 2 * mid - Vr[, ax]
  refl$vb[1:3, ] <- t(Vr); if (!is.null(refl$it)) refl$it <- refl$it[c(1, 3, 2), , drop = FALSE]
  # W(p): where each original point lands when carried toward its mirror partner.
  if (refine) {
    fit <- deformetrica_register(x, refl, ...)
    W <- nat::xyzmatrix(deformetrica_shoot(V, fit$control_points, fit$momenta,
                                           kernel_width = fit$kernel_width, ...))
  } else {
    W <- Vr
  }
  travel <- W[, ax] - V[, ax]
  Lsign <- if (positive_travel == "L") 1 else -1
  tol <- band * diff(range(V[, ax]))
  side <- ifelse(abs(travel) <= tol, "M", ifelse(sign(travel) == Lsign, "L", "R"))
  vside <- factor(side, levels = c("L", "R", "M"))
  fside <- apply(matrix(as.integer(vside)[x$it], nrow = 3), 2,
                 function(v) which.max(tabulate(v, 3)))    # face = majority vertex side
  list(L = .submesh_faces(x, which(fside == 1L)),
       R = .submesh_faces(x, which(fside == 2L)),
       side = vside, travel = travel, mid = mid)
}

#' Animate a Deformetrica geodesic flow as a nat.ggplot GIF
#'
#' A thin convenience wrapper around [nat.ggplot::ggneuron_gif()] (where the general
#' "animate a sequence of nat objects to a GIF" machinery now lives). It renders a set
#' of warping objects across the timepoints of a geodesic flow (as returned by
#' [deformetrica_shoot()] with `flow = TRUE`) to a ping-pong GIF, each object in its
#' own colour over an optional reference volume and static greyscale targets. Kept for
#' convenience/back-compatibility; new code can call [nat.ggplot::ggneuron_gif()]
#' directly.
#'
#' @param flows A named list; each element is the per-timepoint list of objects for
#'   one structure (all lists the same length). Colours are taken per structure.
#' @param cols Named vector of colours (one per `flows` entry); recycled from a
#'   default palette if `NULL`.
#' @param volume Optional fixed reference object drawn translucent under every frame
#'   (e.g. the target brain hull).
#' @param volume_col Colour of the reference volume (default `"grey80"`; pass e.g. a
#'   light pink to render the brain hull as a soft envelope).
#' @param targets Optional named list of fixed reference objects drawn translucent
#'   *above* the volume but *below* the flow (e.g. the matched target neuropils each
#'   warping object should land on). Static across frames.
#' @param target_cols Named vector of colours for `targets` (one per entry); defaults
#'   to a transparent greyscale ramp so the coloured flow reads clearly on top.
#' @param volume_alpha,target_alpha,alpha Alphas for the reference volume, the target
#'   objects and the flow objects.
#' @param points Optional named list of point sets drawn as circles on top of the
#'   neurons - e.g. neuron root points / somata. Each entry is either a static N x 3
#'   matrix or a per-timepoint list of matrices (circles that move with the warp),
#'   forwarded to [nat.ggplot::ggneuron_gif()].
#' @param point_cols Named fill colours for `points` (one per entry).
#' @param point_size,point_alpha,point_stroke Circle size, opacity and outline width
#'   for `points`.
#' @param rotation_matrix Optional 4x4 view matrix (as `nat.ggplot`/`rgl` use).
#' @param file Output GIF path. Frames are written next to it.
#' @param width,height,delay,dpi GIF frame size (px), per-frame delay (s) and dpi.
#' @return The GIF path if written (needs `gifski`, or falls back to `magick`), else
#'   the vector of frame PNGs.
#' @seealso [nat.ggplot::ggneuron_gif()] for the underlying, general animation function.
#' @export
ggplot_flow_gif <- function(flows, cols = NULL, volume = NULL, volume_alpha = 0.12,
                            volume_col = "grey80", targets = NULL, target_cols = NULL,
                            target_alpha = 0.18, alpha = 0.6, rotation_matrix = NULL,
                            points = NULL, point_cols = NULL, point_size = 3,
                            point_alpha = 1, point_stroke = 0.8,
                            file = NULL, width = 900, height = 800, delay = 0.14,
                            dpi = 96) {
  if (!requireNamespace("nat.ggplot", quietly = TRUE))
    stop("ggplot_flow_gif() needs 'nat.ggplot' (>= 1.1.6).", call. = FALSE)
  nat.ggplot::ggneuron_gif(
    flows, cols = cols, volume = volume, volume_col = volume_col,
    volume_alpha = volume_alpha, targets = targets, target_cols = target_cols,
    target_alpha = target_alpha, alpha = alpha, rotation_matrix = rotation_matrix,
    points = points, point_cols = point_cols, point_size = point_size,
    point_alpha = point_alpha, point_stroke = point_stroke,
    file = file, width = width, height = height, delay = delay, dpi = dpi)
}
