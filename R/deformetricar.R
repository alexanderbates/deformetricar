#' deformetricar: an R client for the Deformetrica shape-registration toolkit
#'
#' \bold{deformetricar} wraps the Deformetrica (>= 4.3) command-line tools for
#' statistical shape analysis: it fits diffeomorphisms between point clouds, neuron
#' backbones and surface meshes, and applies them by geodesic shooting. It also
#' reads and writes VTK object formats.
#'
#' @section Package options:
#' \describe{
#'   \item{\code{deformetricar.exe}}{Path to (or name of) the \code{deformetrica}
#'     executable; consulted by \code{\link{find_deformetrica}}. Set it in your
#'     \code{\link[base]{Rprofile}} if Deformetrica is not on the \code{PATH} or in
#'     a conda environment.}
#' }
#'
#' @references \url{https://www.deformetrica.org}
#' @seealso \code{\link{find_deformetrica}}, \code{\link{deformetrica_register}},
#'   \code{\link{deformetrica_shoot}}
#' @keywords internal
#' @importFrom nat xyzmatrix xyzmatrix<- xform
#' @importFrom stats kmeans
#' @importFrom utils write.table
"_PACKAGE"
