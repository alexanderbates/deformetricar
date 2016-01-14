


otherside <- function (x, ...){
  x = apply.mirror.affine(x, ...)
  x = shootflow(x, ...)
  x
}
