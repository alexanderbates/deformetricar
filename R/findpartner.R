#' Find neuron cognates on other side of a structure
#'
#' @param x a neuron or neuronlist from the Drosphila L1 connectome.
#' @param db a batabse of neurons against which to search for partners. Preferably a dotprops object.
#' @param ... additional arguments passed to: \code{\link{shootflow}} \code{\link{apply.mirror.affine}}
#'
#' @return "Scores" -  a matrix tabulating similarity scores (NBlast) for neurons on the otherside of the brain,
#' neurons being identified by skid, or "neuronlist" - a neuronlist containing the neurons with the highest similarity scores.
#' @export
#' @rdname findpartner
#' @seealso \code{\link{shootflow}} \code{\link{apply.mirror.affine}} \code{\link{otherside}}
findpartner<-function(x, ...) UseMethod("findpartner")

#' @export
#' @rdname findpartner
findpartner.default <- function (x, db, output = c("scores","neuronlist"), ...){
  x = otherside(x)
  x.dps = nat::dotprops(x, resample=0.1, .progress = 'text')
  if (!nat::is.dotprops(db))
    cat("Converting database neurons into dotprops. This can slow things down if the neuronlist is large, better to provide a dotprops neuronlist.")
    db = nlapply(db, dotprops, resample=0.1, .progress = 'text')
  l1_smat = readRDS(system.file("extdata/l1_smat.rds", package = 'deformetricar'))
  results = nat::nblast(x.dps, target = db, smat = l1_smat,  normalised = T, .progress = 'text', UseAlpha = T)
  if (output == "neuronlist")
    partners = c()
  for (col in 1:length(x))
    c = results[col]
  partners = c(partners, names(sort(c, decreasing = T))[1])
  partners = unique(partners)
  require(catmaid)
  results = catmaid::read.neurons.catmaid(partners)
  results
}

#' @export
#' @rdname findpartner
findpartner.neuron <- function (x, db, output = c("scores","neuronlist"), ...){
  x = otherside(x)
  x.dps = nat::dotprops(x, resample=0.1, .progress = 'text')
  if (!nat::is.dotprops(db))
    cat("Converting database neurons into dotprops. This can slow things down if the neuronlist is large, better to provide a dotprops neuronlist.")
    db = nlapply(db, dotprops, resample=0.1, .progress = 'text')
  l1_smat = readRDS(system.file("extdata/l1_smat.rds", package = 'deformetricar'))
  results = nat::nblast(x.dps, target = db, smat = l1_smat,  normalised = T, .progress = 'text', UseAlpha = T)
  if (output == "neuronlist")
    partners = c()
    for (col in 1:length(x))
      c = results[col]
      partners = c(partners, names(sort(c, decreasing = T))[1])
    partners = unique(partners)
    require(catmaid)
    results = catmaid::read.neurons.catmaid(partners)
  results
}

#' @export
#' @rdname findpartner
findpartner.neuronlist <- function (x, db ...){
  results = list()
  for (someneuron in 1:length(x))
    result = findpartner.neuron(x = x[[someneuron]], db = db, output = "scores")
    results = list(results, list(result))
    results[[length(results)+1]] <- result
  result
}
