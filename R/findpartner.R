#' Find neuron cognates on other side of a nervous system.
#'
#' @param x a neuron or neuronlist from the Drosphila L1 connectome.
#' @param db a database of neurons against which to search for partners. Preferably a dotprops object.
#' @param entries number of possible ranked cognates to display
#' @param possiblecognates object returned from calling findcognate
#' @param ... additional arguments passed to: \code{\link{otherside}} Necessary if not using defauly nervous system fro this package, the larval L1 reconstruction by Cardona and colleagues.
#'
#' @return "Scores" -  a matrix tabulating similarity scores (NBlast) for neurons on the otherside of the brain,
#' neurons being identified by skid, or "neuronlist" - a neuronlist containing the neurons with the highest similarity scores.
#' @export
#' @rdname findcognate
#' @seealso \code{\link{shootflow}} \code{\link{apply.mirror.affine}} \code{\link{otherside}}

#' @export
#' @rdname findcognate
findcognate <- function (x, db, entries = 10, ...){
  message("Remember, give all neuron in microns")
  cat("Transforming neuron into the other side of the brain")
  if (nat::is.neuronlist(x)||nat::is.neuron(x))
    x = otherside(x)
  if (is.character(x))
    x = otherside(db[x])
  cat("Making vector cloud objects from neurons")
  if (!nat::is.dotprops(x))
    x.dps = nat::nlapply(x, nat::dotprops, resample=0.1, .progress = 'text', OmitFailures = T)
  else{ x.dps = x }
  if (!nat::is.dotprops(db[[1]])){
    warning("Converting database neurons into dotprops. This can slow things down if the neuronlist is large, better to provide a dotprops neuronlist.")
    db = nat::nlapply(db, nat::dotprops, resample=0.1, .progress = 'text', OmitFailures = T)
  }
  cat("Running NBlast")
  l1_smat = readRDS(system.file("extdata/l1_smat.rds", package = 'deformetricar'))
  results = nat.nblast::nblast(x.dps, target = db, smat = l1_smat,  normalised = T, .progress = 'text', UseAlpha = T)
  output = list()
  for (n in 1:ncol(results)){
    scores = sort(results[,n], decreasing = T)[1:entries]
    df= data.frame(pid = as.data.frame(db[names(scores)])$pid, skid = names(scores), name = as.data.frame(db[names(scores)])$name, score = scores)
    data = list(as.data.frame(db[colnames(results)[n]]),df)
    names(data) <- c("Query","Matching")
    output = append(output, list(data))
  }
  names(output) <- colnames(results)
  output
}

#' @export
#' @rdname findcognate
checkL1cognates <- function(possiblecognates, db, entries = 10, ...){
  cat("GREY query neuron \nBLACK transformed neuron used for search \nCOLOURS possible cognates, the redder the better")
  jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))
  colours = jet.colors(entries)[entries:1]
  rgl::open3d(userMatrix = structure(c(0.999989449977875, -0.00419542612507939, -0.00183705519884825,
                                0, -0.00426193978637457, -0.999274253845215, -0.0378473773598671,
                                0, -0.00167691614478827, 0.037854727357626, -0.999281764030457,
                                0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.5050682, windowRect = c(1580L, 136L, 2768L, 1051L))
  models=catmaid::catmaid_fetch("1/stack/5/models")
  l1surf=matrix(as.numeric(models$cns$vertices), ncol=3, byrow = TRUE)
  l1mesh = nat::xyzmatrix(l1surf)/1000 #l1 mesh points in microns
  for (n in 1:length(possiblecognates)){
    cat("Current query neuron:", possiblecognates[[n]]$Query$name, "(", n, "/", length(possiblecognates),
        ")\n")
    message("Overview of matches")
    print(data.frame(possiblecognates[[n]]$Matching, colour = colours))
    rgl::points3d(l1mesh, col = 'lightgrey')
    query = db[as.character(possiblecognates[[n]]$Query$skid)]
    query.transformed = otherside(db[as.character(possiblecognates[[n]]$Query$skid)])
    rgl::plot3d(query, col = 'darkgrey')
    rgl::plot3d(query.transformed, col = 'black')
    matches = db[as.character(possiblecognates[[n]]$Matching$skid[1:entries])]
    rgl::plot3d(matches, col=colours, scalevecs=0.1)
    # Continue
    print ("Press [enter] to continue and 'b' to go back")
    number <- scan(n=1)
    message("Looking closer")
    e = 0
    progress = 'forward'
    while(e < entries){
      if (progress == 'b'){e = e - 1; progress = 'forward'
      }else{e = e + 1}
      if (e<1){e = 1}
      nat::npop3d()
      rgl::plot3d(matches[e], col = colours[e])
      print(data.frame(possiblecognates[[n]]$Matching[e,]))
      progress = readline()
    }
    number <- scan(n=1)
    rgl::clear3d()
  }
    rgl::clear3d()
}
