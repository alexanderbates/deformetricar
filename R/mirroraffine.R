# Perform the initial mirroring, rigid affine and non-rigid affine transformations

#To generate the full transformation matrix, in preparation for the deformetrica deformation, I have generated a series
# of point clouds that represent different neuroanatomical features of the larval brain of interest, incluidng the longitudinal
#tracts down the ventrla nerve cord (VNC) that would be marked by a FasII stain, the three major olfactory tracts and
# the most prominent neuropils observable at a light level and within the EM data, namely the mushroom bodies (with their
#calyxes generated separately) the antennal lobes, and the entire neuropil as would be seen in an Nc82 (bruchpilot) stain.

#' @export
#' @rdname apply.mirror.affine
apply.mirror.affine<-function(x, ...) UseMethod("apply.mirror.affine")

#' @export
#' @rdname apply.mirror.affine
apply.mirror.affine.neuron<-function(x, ...) {
  moved.points<-apply.mirror.affine.default(x, ...)
  x$d$X <- moved.points[,1]
  x$d$Y <- moved.points[,2]
  x$d$Z <- moved.points[,3]
  x
}

#' @export
#' @rdname apply.mirror.affine
apply.mirror.affine.neuronlist<-function(x, ...){
  nat::nlapply(x, apply.mirror.affine, ...)
}

#' Apply a mirroring and affine transformation to an object in 3D space
#'
#' @param x a set of 3D coordinates (Nx3 matrix), a neuron, neuronlist, dotprops, hxsurf, igraph or mesh3d object.
#' currently point replacement is only supported for neuron and neuronlist objects.
#' @param calculatetransform whether or not to re-calculate the mirorring transformation, and the affine transformation for
#' the mirrored object onto the original ones. Null defaults to a list of transformation matrices generated from flipping a set
#' of neuroanatomical structures in the L1 larval central nervous system, then performing an affine transformation. Otherwise,
#' @param ... additional arguments eventually passed to: \code{\link{calculate.full.transformation}} \code{\link{mirrormat}} \code{\link{transform.points}}
#' @return a set of 3D cordinates (or a neuron/neuronlist object if that was given as input), that has undergone sequential mirroring
#' and affine transformations. The default is to flip points in the l1 larval central nervous system from one side to the
#' other.
#' @export
#' @seealso \code{\link{calculate.full.transformation}} \code{\link{mirrormat}} \code{\link{transform.points}}
#' \code{\link{trafoicpmat}}
#' @rdname apply.mirror.affine
apply.mirror.affine.default <- function (x, calculatetransform = NULL, ...){
  xyz = nat::xyzmatrix(x)
  if (is.null(calculatetransform))
    fullmirror = readRDS("inst/extdata/fullmirror.rds")
  else
    fullmirror = calculate.full.transformation()
  transform.points(xyz, fullmirror)
}


#' Calculate a mirroring, then a pair of affine transformations
#' @description Calculates a mirroring transformation, then a rigid, then a non-rigid affine transformation on a given 3D R object, or
#' one(s) saved as RDS files By default, the mirroring occurs across the X axis of the provided data.
#' @param objs An R object () or the location of (an) R object(s) for the basis of calculating the transformation matrices.
#' Defaults to a set of files describing key neuroanatomical structures of the L1 larval central nervous system.
#' @param pattern search term for RDS files to read in the given location.
#' @param ... additional arguments eventually passed to: \code{\link{trafoicpmat}} \code{\link{mirrormat}}
#'
#' @return a list of transformation matrices
#' @export

calculate.full.transformation <- function (objs = "inst/extdata/point_objects/", pattern = ".rds$", ...){
  if (is.character(objs)){
    filelist = list.files(objs, pattern = pattern)
    all.structures = matrix(ncol=3)
    for(i in 1:length(filelist))
    {
      oname = paste(gsub(pattern,".d",filelist[i]) )
      all.structures = rbind(all.structures, assign(oname, readRDS(paste(objs,filelist[i], sep = ""))))
    }
    all.structures = all.structures[-1,]
  } else
      all.structures = xyzmatrix(objs)
  # Generate transformation matrices
  flipmatrix = mirrormat(all.structures)
  all.structures.flipped = nat::mirror(all.structures, boundingbox(apply(all.structures, 2, range, na.rm = T)))
  t.rigid =trafoicpmat(all.structures.flipped, all.structures, 100, subsample = NULL, type='rigid')
  t.affine =trafoicpmat(t.rigid$xt, all.structures, 100, subsample = NULL, type='affine')
  # This is the matrix for the full, global transformation
  list(flipmatrix, t.rigid$trafo, t.affine$trafo)
}


#' Generate a transformation matrix for a mirroring registration
#' @description Generates
#' @param x a set of 3D coordinates (Nx3 matrix), a neuron, neuronlist, dotprops, hxsurf, igraph or mesh3d object.
#' @param mirrorAxis the axis along which to mirror the data. Defaults to the X axis.
#' @param mirrorAxisSize this defaults to the bounding box for x.
#' @param ... additional arguments eventually passed to methods
#'
#' @return A transformation matrix
#' @export
mirrormat <- function (x, mirrorAxis = c("X", "Y", "Z"), mirrorAxisSize = NULL, ...)
{
  if (is.null(mirrorAxisSize))
    bb = apply(x, 2, range, na.rm = T)
    mirrorAxisSize = nat::boundingbox(bb)
  if (is.character(mirrorAxis)) {
    mirrorAxis = match.arg(mirrorAxis)
    mirrorAxis = match(mirrorAxis, c("X", "Y", "Z"))
  }
  if (length(mirrorAxis) != 1 || is.na(mirrorAxis) || mirrorAxis <
      0 || mirrorAxis > 3)
    stop("Invalid mirror axis")
  lma = length(mirrorAxisSize > 1)
  if (lma > 1) {
    if (lma == 6)
      mirrorAxisSize = mirrorAxisSize[, mirrorAxis]
    else if (lma != 2)
      stop("Unrecognised mirrorAxisSize specification!")
    mirrorAxisSize = sum(mirrorAxisSize)
  }
  mirrormat = diag(4)
  mirrormat[mirrorAxis, 4] = mirrorAxisSize
  mirrormat[mirrorAxis, mirrorAxis] = -1
  return(mirrormat)
}


#' Get affine transformation matrices
#' @description A modification of the icpmat function from the package Morpho. This function matches
#' two landmark configurations using iteratively closest point search.
#' @param x moving landmarks
#' @param y target landmarks
#' @param iterations integer number of iterations
#' @param mindist restrict valid points to be within this distance
#' @param subsample use a subsample determined by kmean clusters to speed up computation
#' @param type character: select the transform to be applied, can be "rigid","similarity" or "affine"
#' @param ... additional arguments eventually passed to methods
#'
#' @return a data frame containing the transformed object, and the transformation matrices describing the performed
#' transformation
#' @export
trafoicpmat <- function (x, y, iterations, mindist = 1e+15, subsample = NULL,
                      type = c("rigid", "similarity", "affine"), ...) {
   m <- ncol(x)
  if (m == 2) {
    x <- cbind(x, 0)
    y <- cbind(y, 0)
  }
  type <- match.arg(type, c("rigid", "similarity", "affine"))
  if (!is.null(subsample)) {
    subsample <- min(nrow(x) - 1, subsample)
    subs <- duplicated(kmeans(x, centers = subsample, iter.max = 100)$cluster)
    xtmp <- x[!subs, ]
  }
  else {
    xtmp <- x
  }
  for (i in 1:iterations) {
    clost <- Rvcg::vcgKDtree(y, xtmp, 1)
    good <- which(clost$distance < mindist)
    trafo <- Morpho::computeTransform(y[clost$index[good], ], xtmp[good,
                                                                   ], type = type)
    xtmp <- Morpho::applyTransform(xtmp[, ], trafo)
  }
  if (!is.null(subsample)) {
    fintrafo <- Morpho::computeTransform(xtmp[, ], x[!subs, ], type = type)
    xtmp <- Morpho::applyTransform(x, fintrafo)
  }
  if (m == 2)
    xtmp <- xtmp[, 1:2]
  # compute transform of original floating points onto final floating points
  trafo2= Morpho::computeTransform(xtmp, x, type=type)
  return(list(xt=xtmp, trafo=trafo2))
}


#' Transform 3D coordinates using a list of transformation matrices
#'
#' @param positions 3D coodinates
#' @param transformations list of transformation matrices
#' @param ... additional arguments eventually passed to methods
#'
#' @return the transformed 3D coordinates
#' @export
transform.points = function (positions, transformations, ...){
  if (is.list(transformations) == F){
    positions = nat::xform(positions, transformations)
  }
  if (is.list(transformations) == T){
    for (transformation in transformations){
      positions <- xform(positions, transformation)
    }
  }
  return (positions)
}



