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
  otherside<-function(x, method = c("tps3d", "deformetrica"),...) UseMethod("otherside")

#' @export
#' @rdname otherside
otherside.default <- function (x, method = c("tps3d", "deformetrica"), ...){
  x = apply.mirror.affine(x, ...)
  method = match.arg(method)
  if (method == "deformetrica")
    x = shootflow(x, ...)
  else if (method == "tps3d")
    finals = read.vtk(system.file("extdata/reg_output/finals.vtk", package = 'deformetricar'))
    cps = read.points(system.file("extdata/reg_output/finals.vtk", package = 'deformetricar'))
    x = tps3d(x, cps, finals, lambda = 0)
  x
}

#' @export
#' @rdname otherside
otherside.neuron<-function(x, object.type = "NonOrientedPolyLine", ...) {
  points=xyzmatrix(x)
  morph.points<-otherside.default(x=points, ...)
  xyzmatrix(x)<-morph.points
  x
}

#' @export
#' @rdname otherside
otherside.neuronlist<-function(x, method = c("tps3d", "deformetrica"), object.type = "NonOrientedPolyLine", ...){
  points=xyzmatrix(x)
  morph.points<-otherside.default(x=points, ...)
  count = 1
  for (neuron in 1:length(x)){
    n = count + nrow(xyzmatrix(x[[neuron]])) -1
    xyzmatrix(x[[neuron]]) <- morph.points[count:n,]
    count = count + nrow(xyzmatrix(x[[neuron]]))
  }
  return(x)
}
