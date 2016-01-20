#' Apply a list of mirroring, affine and warping registrations to an object in 3D space
#'
#' @param x a set of 3D coordinates (Nx3 matrix), a neuron, neuronlist, dotprops, hxsurf, igraph or mesh3d object.
#' currently point replacement is only supported for neuron and neuronlist objects.
#' @param ... additional arguments passed to: \code{\link{shootflow}} \code{\link{apply.mirror.affine}}
#'
#' @return a set of 3D cordinates (or a neuron/neuronlist object if that was given as input), that has undergone sequential mirroring
#' affine and warping (via Deformetrica, see http://www.deformetrica.org/?page_id=232) transformations. The default is to flip points in
#' the l1 larval central nervous system from one side to the other, and apply a previousy generated Deformetrica deformation based on the
#' registering of a set of flipped, large neuroanatomical regions unto their un-flipped counterparts.
#' @export
#' @rdname otherside
#' @seealso \code{\link{shootflow}} \code{\link{apply.mirror.affine}}
otherside.default <- function (x, ...){
  x = apply.mirror.affine(x, ...)
  x = shootflow(x, ...)
  x
}

#' @export
#' @rdname otherside
otherside.neuron<-function(x, ...) {
  points=nat::xyzmatrix(x)
  morph.points<-otherside.default(x=points, ...)
  x = apply.mirror.affine(x, ...)
  x = shootflow(x, ...)
  x$d$X <- morph.points[,1]
  x$d$Y <- morph.points[,2]
  x$d$Z <- morph.points[,3]
  x
}

#' @export
#' @rdname otherside
otherside.neuronlist<-function(x, ...){
  points=nat::xyzmatrix(x)
  morph.points<-otherside.default(x=points, ...)
  count = 1
  for (neuron in 1:length(x)){
    n = count + nrow(xyzmatrix(x[[neuron]])) -1
    nat::xyzmatrix(x[[neuron]]) <- morph.points[count:n,]
    count = count + nrow(nat::xyzmatrix(x[[neuron]]))
  }
  return(x)
}
