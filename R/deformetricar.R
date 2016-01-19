#' R wrapper for elements of the Deformetrica registration tool
#'
#' \bold{deformetricar} wraps the Deformetrica tool for statistical analysis of
#' 2D and 3D shape data.
#'
#' @section Package Options: The following options can be set to specify default
#'   behaviour.
#'
#'   \itemize{
#'
#'   \item{\code{deformeticar.bindir}}{ Location of Deformetrica binaries}. You
#'   may wish to set this in your \code{\link[base]{Rprofile}} if you do not
#'   have deformeticra installed in a standard location.
#'
#'   }
#'
#' @name deformetricar-package
#' @references See \url{http://www.deformetrica.org/}
#' @aliases deformetricar, deformetrica
#' @seealso \code{\link{shootflow}}
#' @examples
#' # Show state of deformetricar package options
#' options()[grep('^deformetricar', names(options()))]
#' @docType package
#' @keywords package
NULL
