# shoot and flow functions
# idea is that we want to be able to apply an existing registratino to
# an abitrary set of coordinates

# see http://www.deformetrica.org/?page_id=232


#' Apply deformetrica registration to a set of 3D points (ShootAndFlow3)
#' @description Apply a Deformetrica deformation using calculated momentum vectors and control points (in files MOM_final.txt
#' and CP_final.txt respectively), by calling the executable ShootAndFlow3. This requires installation of Deformetrica version 2.1 or later (preferably in Applications).
#' @export
#' @rdname shootflow
shootflow<-function(x, ...) UseMethod("shootflow")

#' @export
#' @rdname shootflow
shootflow.neuron<-function(x, ...) {
  points=nat::xyzmatrix(x)
  deform.points<-shootflow.default(x=points, ...)
  x$d$X <- deform.points[,1]
  x$d$Y <- deform.points[,2]
  x$d$Z <- deform.points[,3]
  x
}

#' @export
#' @rdname shootflow
shootflow.neuronlist<-function(x, ...){
  if (mean(xyzmatrix(x)[,3]) > 1000)
    warning("Object appears to be in nanometers")
  nat::nlapply(x, shootflow, ...)
}

#' @param x object to be transformed (for \code{shootflow.default} method an Nx3
#'   matrix of 3D coordinates)
#' @param regdir Path to directory containing deformetrica registration
#' @param kernel.width the width of the deformation kernel. The larger the more “rigid”
#' the deformation. The smaller, the more local variations of the space is allowed.
#' @param object.type type of object to be deformed. PointCloud, OrientedPolyLine, NonOrientedPolyLine,
#' OrientedSurfaceMesh and NonOrientedSurfaceMesh
#' @param ... additional arguments eventually passed by methods
#'   \code{shootflow.default}
#'
#' @return Matrix of the same dimensions as \code{x}
#' @export
#' @rdname shootflow
#' @references See \url{http://www.deformetrica.org/?page_id=232} for details of
#'   the \bold{ShootAndFlow3} command line tool.
shootflow.default<-function(x, regdir = "inst/extdata/reg_output", kernel.width=5,
                    object.type = "NonOrientedPolyLine", ...){
  # we need to make a command line like this
  # ShootAndFlow3 paramsDiffeos.xml Direction CP.txt Mom.txt paramsObject1.xml object1 paramsObject2.xml object2 …
  # make a temp dir
  td=tempfile()
  dir.create(td)
  # clean up when done
  on.exit(unlink(td, recursive = TRUE))
  reg_files_we_want=c("paramDiffeos.xml", "CP_final.txt", "Mom_final.txt")
  owd=setwd(regdir)
  on.exit(setwd(owd), add = T)
  if(!all(file.exists(reg_files_we_want))){
    stop("cannot find some of the registration input files:", reg_files_we_want)
  }
  file.copy(reg_files_we_want, td)

  # now change to our temporary directory
  setwd(td)
  params_file=make_params_file(kernel.width = kernel.width, object.type = object.type)
  steps=read.paramdiffeos("paramDiffeos.xml")$steps
  write.vtk(x, "points.vtk")
  ShootAndFlow3("paramDiffeos.xml", 1, "CP_final.txt", "Mom_final.txt", params_file, "points.vtk")
  output.file = paste("points_flow__t_", steps, ".vtk", sep = "")
  read.vtk(output.file)
}


#' @rdname shootflow
read.paramdiffeos<-function(infile){
  ll=readLines(infile)
  rl=list()
  rl$steps = as.integer(unlist(strsplit(gsub("[^0-9]", "", ll[25]), "")))-1
  rl
}

#' @rdname shootflow
make_params_file<-function(outfile="paramsObject1.xml", kernel.width=5, object.type = c("PointCloud, OrientedPolyLine, NonOrientedPolyLine, OrientedSurfaceMesh, NonOrientedSurfaceMesh")){
  data.sigma=1
  lines=c("<?xml version=\"1.0\"?>", "<deformable-object-parameters>",
    "", "    <!-- Type of the deformable object (See DeformableObject::DeformableObjectType for details) -->",
    sprintf("\t<deformable-object-type>%s</deformable-object-type>", object.type),
    "    <!-- weight of the object in the fidelity-to-data term -->",
    sprintf("\t<data-sigma>%d</data-sigma>",data.sigma), "    <!-- kernel width -->",
    sprintf("\t<kernel-width>%d</kernel-width>",kernel.width),"", "</deformable-object-parameters>"
  )
  writeLines(lines, outfile)
  return(outfile)
}


#' @rdname shootflow
# Call the ShootAndFlow3 with optional arguments
ShootAndFlow3<-function(...){
  # we need to find the deformetrica executable
  # is it in the path?
  app=Sys.which('ShootAndFlow3')
  if(nchar(app)==0){
    # is it in a standard location?
    app="/Applications/deformetrica-2.1/deformetrica/bin/ShootAndFlow3"
    if(!file.exists(app))
      stop("Sorry I can't find the deformetrica ShootAndFlow3 executable!",
           "Please add it to your PATH or put it in ", dirname(app))
  }
  args=paste(..., collapse=" ")
  system(paste(app, args), ignore.stdout=TRUE)
}


