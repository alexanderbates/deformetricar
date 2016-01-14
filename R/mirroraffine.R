# Perform the initial mirroring, rigid affine and non-rigid affine transformations

#To generate the full transformation matrix, in preparation for the deformetrica deformation, I have generated a series
# of point clouds that represent different neuroanatomical features of the larval brain of interest, incluidng the longitudinal
#tracts down the ventrla nerve cord (VNC) that would be marked by a FasII stain, the three major olfactory tracts and
# the most prominent neuropils observable at a light level and within the EM data, namely the mushroom bodies (with their
#calyxes generated separately) the antennal lobes, and the entire neuropil as would be seen in an Nc82 (bruchpilot) stain.


mirror.affine <- function (){

}


calculate.full.transformation <- function (dir = "inst/extdata/point_objects/", pattern = ".rds$", ...){
  filelist = list.files(dir, pattern = pattern)
  all.structures = matrix(ncol=3)
  for(i in 1:length(filelist))
  {
    oname = paste(gsub(".rds",".d",filelist[i]) )
    all.structures = rbind(all.structures, assign(oname, readRDS(paste(dir,filelist[i], sep = ""))))
  }
  all.structures = all.structures[-1,]
  # Generate transformation matrices
  flipmatrix = mirrormat(all.structures)
  all.structures.flipped = nat::mirror(all.structures, boundingbox(apply(all.structures, 2, range, na.rm = T)))
  t.rigid =trafoicpmat(all.structures.flipped, all.structures, 100, subsample = NULL, type='rigid')
  t.affine =trafoicpmat(t.rigid$xt, all.structures, 100, subsample = NULL, type='affine')
  # This is the matrix for the full, global transformation
  list(flipmatrix, t.rigid$trafo2, t.affine$trafo2)
}


mirrormat <- function (x, mirrorAxis = c("X", "Y", "Z"), mirrorAxisSize = NULL, ...)
{
  if (is.null(mirrorAxisSize))
    bb = apply(x, 2, range, na.rm = T)
    mirrorAxisSize = boundingbox(bb)
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


trafoicpmat <- function (x, y, iterations, mindist = 1e+15, subsample = NULL,
                      type = c("rigid", "similarity", "affine"), ...) {
  require(Morpho)
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
  return(list(xt=xtmp, trafo=trafo, trafo2=trafo2))
}


transform.points = function (positions, transformations, ...){
  require(nat)
  if (is.list(transformations) == F){
    positions = xform(positions, transformations)
  }
  if (is.list(transformations) == T){
    for (transformation in transformations){
      positions <- xform(positions, transformation)
    }
  }
  return (positions)
}



