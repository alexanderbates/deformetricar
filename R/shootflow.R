# shoot and flow functions
# idea is that we want to be able to apply an existing registratino to
# an abitrary set of coordinates

# see http://www.deformetrica.org/?page_id=232

# first thing to achieve support for bare vertices (3d coordinates)



#' Apply deformetrica registration to a set of 3D points (ShootAndFlow3)
#'
#' @inheritParams write.vtk
#' @param regdir Path to directory containing deformetrica registration
#' @param data.sigma See deformetrica docs
#' @param kernel.width See deformetrica docs
#' @param object.type Type of object to be deformed
#'
#' @return Matrix of the same dimensions as \code{points}
#' @export
#' @references See \url{http://www.deformetrica.org/?page_id=232} for
#' details of the \bold{ShootAndFlow3} command line tool.
shootflow<-function(points, regdir = getwd(), data.sigma = 1, kernel.width=5,
                    object.type = "NonOrientedPolyLine"){
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
  params_file=make_params_file(data.sigma = data.sigma, kernel.width = kernel.width, object.type = object.type)
  steps=read.paramdiffeos("paramDiffeos.xml")$steps
  write.vtk(points, "points.vtk")
  ShootAndFlow3("paramDiffeos.xml", 1, "CP_final.txt", "Mom_final.txt", params_file, "points.vtk")
  output.file = paste("points_flow__t_", steps, ".vtk", sep = "")
  read.vtk(output.file)
}

read.paramdiffeos<-function(infile){
  ll=readLines(infile)
  rl=list()
  rl$steps = as.integer(unlist(strsplit(gsub("[^0-9]", "", ll[25]), "")))-1
  rl
}

make_params_file<-function(outfile="paramsObject1.xml", data.sigma=1, kernel.width=5, object.type = "NonOrientedPolyLine"){
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
  system(paste(app, args))
}
