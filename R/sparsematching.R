#' Run a Deformatrica sparseMatching3 deformation on 3D objects
#'
#' @param regdir Path to directory for performing the deformetrica registration. The .vtk files produced for 3D objects and the registrations generated will be saved here
#' @param moving.list list of mesh3D objects to be deformed
#' @param target.list list of mesh3D objects the moving objects are goign to try to match
#' @param choose.params whether or not to manually select values for kernel.width, data.sigma and object.type for eaxch object in the registration. Else, default values will be used
#' @param plotting whether to plot the 3D objects in question ot help you decide on the above
#' @param verbose whether or not to see the progress of the sparseMatching3 system command in your R console
#' @param kernel.width the width of the deformation kernel for individual objects
#' @param data.sigma matching factor. The lower ir is, the tighter the matching
#' @param object.type type of object to be deformed. PointCloud,
#'   OrientedPolyLine, NonOrientedPolyLine, OrientedSurfaceMesh and
#'   NonOrientedSurfaceMesh
#' @param k.width the width of the deformation kernel for the whoel deformation. The larger the more
#'   rigid the deformation. The smaller, the more local variations of the space
#'   is allowed
#' @param kernel.type the kernel type used. See the deformetrica website for details
#' @param no.threads the number of thread to use on your machine during the calculation
#' @param no.timepoints the number of time stesp at which ti generate a .vtk files for the moving objects. The higher this value, the slower the calculation
#' @param no.iterations how many iterations to go through before terminating the registration attempt
#' @param smooth.k.width.ratio the degree of smoothening to kernel width ratio
#' @param sparsity prior see deformetrica website for details
#' @param The path for the deformetrica app on your machine
#' @param ... additional arguments eventually passed by methods
#'
#' @return A list of deofrmed mesh3D objects, with the deformation stored in its attributes
#' @export
#' @rdname sparseMatching3
#' @references See \url{http://www.deformetrica.org/?page_id=189} for details of
#'   the \bold{sparseMatching3} command line tool.
sparsematching3<-function(regdir, choose.params = T, verbose = T, moving.list, target.list, k.width=10, kernel.type = c("Exact","cuda","fgt","p3m"), no.threads = 4, no.timepoints = 5, sparsity.prior = 0, max.iterations = 100, smooth.k.width.ratio = 1.0, data.sigma = 1, kernel.width = 10, data.type = "NoneOrientatedPolyline",plotting = F, app=sparseMatching3("sparseMatching3"), ...){
  owd = getwd()
  on.exit(setwd(owd), add = T)
  if (length(moving.list!=length(target.list))){
    stop("The moving and target lists contain a different number of objects")
  }
  # go to reg dir
  setwd(regdir)
  # Make param diffeos file
  make_diffeo_file(outfile="paramDiffeos.xml", k.width=k.width, kernel.type = kernel.type, no.threads = no.threads, no.timepoints = no.timepoints, sparsity.prior = sparsity.prior, max.iterations = max.iterations, smooth.k.width.ratio = smooth.k.width.ratio)
  if (choose.params){
    for (object in 1:length(moving.list)){
      if (plotting){
        rgl::clear3d()
        rgl::plot3d(moving.list[object], col = "green")
        rgl::plot3d(target.list[object], col = "red", add =T)
      }
      kernel.width <- as.numeric (readline(prompt="Select a value for the kernel width  "))
      type <- as.chracter(readline(prompt="Select an object type (n for neuronlist, m for mesh, p for pointcloud)  "))
      if (type ==n){
       object.type = "NonOrientedPolyLine"
       write.vtk.neurons(moving.list[object], filename = paste(object,"moving.vtk"))
       write.vtk.neurons(target.list[object], filename = paste(object,"target.vtk"))
      }else if (type == p){
       object.type = "NonOrientedSurfaceMesh"
        write.vtk.mesh(moving.list[object], filename = paste(object,"moving.vtk"), trafo = NULL)
        write.vtk.mesh(target.list[object], filename = paste(object,"target.vtk"), trafo = NULL)
      }else{ object.type = "PointCloud"
        write.vtk.mesh(moving.list[object], filename = paste(object,"moving.vtk"), trafo = NULL)
       write.vtk.mesh(target.list[object], filename = paste(object,"target.vtk"), trafo = NULL)
      }
      data.sigma <- as.numeric (readline(prompt="Select a value for data sigma (the lower, the higher the matching)  "))
      outfile = paste(object,"paramObject.xml")
      make_params_file(outfile =outfile, kernel.width = kernel.width, object.type = object.type, data.sigma = data.sigma)
    }
  }else{ # Auto make files using parameters given to the function
    for (object in 1:length(moving.list)){
      outfile = paste(object,"paramObject.xml")
      make_params_file(outfile=outfile, kernel.width = kernel.width, object.type = object.type, data.sigma = data.sigma)
      if (data.type == "NoneOrientatedPolyline"){
        write.vtk.neurons(moving.list[object], filename = paste(object,"moving.vtk"))
        write.vtk.neurons(target.list[object], filename = paste(object,"target.vtk"))
      }else{
        write.vtk.mesh(moving.list[object], filename = paste(object,"moving",".vtk"), trafo = NULL)
        write.vtk.mesh(target.list[object], filename = paste(object,"target",".vtk"), trafo = NULL)
      }
    }
  }
  make_registrationsh_file(path = app, outfile = "registration.sh", movingobjects = sort(list.files(pattern = "moving.vtk")), targetobjects=sort(list.files(pattern = "target.vtk")), paramfiles = sort(list.files(pattern = "paramObject.xml")))
  system("registration.sh", ignore.stdout=!verbose)
  print("files .vtk, control point, parameter and momenta saved in ", regdir)
  cpsi = read.points("CPs_initial.txt")
  cpsf = read.points("CPs_final.txt")
  deformation = Morpho::computeTransform(cpsf, cpsi, type = "tps")
  saveRDS(deformation,"deformation.rds")
  filelist = listfiles(pattern = paste("__t_", no.timepoints, ".vtk", sep = ""))
  deformed = lapply(filelist,read.vtk.to.mesh(x))
  attr(deformed,"trafo") = deformation
  return(deformed)
}

# Hidden functions

make_diffeo_file<-function(outfile="paramDiffeos.xml", k.width=10,
                           kernel.type = c("Exact","cuda","fgt","p3m"), no.threads = 4, no.timepoints = 5, sparsity.prior = 0, max.iterations = 100, smooth.k.width.ratio = 1.0){
  lines=c("<?xml version=\"1.0\"?>", "<sparse-diffeo-parameters>",
          "",
          "<!--=======================-->",
          "<!-- Compulsory parameters -->",
          "<!--=======================-->",
          "<!-- Size of the kernel (default : 0.0) -->",
          sprintf("\t<kernel-width>%d</kernel-width>", k.width),
          "",
          "<!--===================================-->",
          "<!-- Optional parameters (recommended) -->",
          "<!--===================================-->",
          "<!-- Choice of the evaluation method of the kernel : exact, fgt, p3m (default : p3m) -->",
          sprintf("\t<kernel-type>%s</kernel-type>",kernel.type),
          "",
          "<!--===================================-->",
          "<!-- Optional parameters (advanced) -->",
          "<!--===================================-->",
          "	<!-- weight of the sparsity prior (default: 0.0) -->",
          sprintf("\t<sparsity-prior>%d</sparsity-prior>",sparsity.prior),
          "",
          "<!-- Number of threads (default : 1) -->",
          sprintf("\t<number-of-threads>%d</number-of-threads>",no.threads),
          "<!-- Choice of the number of time points between 0 and 10 (default : 5) -->",
          sprintf("\t<number-of-timepoints>%d</number-of-timepoints>",no.timepoints),
          "",
          "<!-- ratio between the smoothing kernel width used in the gradient of template gradient and kernel width (default: 1.0)  -->",
          sprintf("\t<smoothing-kernel-width-ratio>%d</smoothing-kernel-width-ratio>",smooth.k.width.ratio),
          "",
          "<!-- Step of the regular lattice of control points (default : kernel-width) -->",
          "<!-- <initial-cp-spacing>40</initial-cp-spacing> -->",
          "<!-- Enables to freeze the control points or not (default : Off) -->",
          "<!-- <freeze-cp>On</freeze-cp> -->",
          "<!-- File name containg control points position (default: void) -->",
          sprintf("\t<initial-cp-position>CP_initial.txt</initial-cp-position>"),
          "<!-- File name containg initial momenta (default: void) -->",
          sprintf("\t<initial-momenta>MOM_initial.txt</initial-momenta>"),
          "",
          "<!-- Optimization method (default : fista) -->",
          "<!-- <optimization-method-type>GradientDescent</optimization-method-type> -->",
          "",
          "<!--===================================-->",
          "<!-- Optional parameters (avoid changing) -->",
          "<!--===================================-->",
          "<!-- Maximum of descent iterations (default : 100) -->",
          sprintf("\t<max-iterations>%d</max-iterations>",max.iterations),
          " <!-- Adaptive tolerance (default : 1e-04) -->",
          "<!-- <adaptive-tolerance>1e-4</adaptive-tolerance> -->",
          "<!-- Maximum number for the search of the optimal step on an iteration (default : 10) -->",
          "<!-- <max-line-search-iterations>20</max-line-search-iterations> -->",
          "<!-- Step expand for the gradient descent (default : 2.0) -->",
          "<!-- <step-expand>2.0</step-expand> -->",
          "<!-- Step shrink for the gradient descent  (default : 0.5) -->",
          "<!-- <step-shrink>0.5</step-shrink> -->",
          "<!-- Initial step multiplier for gradient descent (default : 10.0) -->",
          "<!-- <initial-step-multiplier>0.001</initial-step-multiplier> -->",
          "",
          "<!-- Enlarges the grid by P3MPaddingFactor x KernelWidth to avoid side effects (default : 3.0) -->",
          "<!-- <p3m-padding-factor>3.0</p3m-padding-factor> -->",
          "<!-- 0.2 gives a relative approximation error of about 5%. Use 0.3 for increased speed (default : 0.2) -->",
          "<!-- <p3m-working-spacing-ratio>0.2</p3m-working-spacing-ratio> -->",
          "",
          "</sparse-diffeo-parameters>"
  )
  writeLines(lines, outfile)
  return(outfile)
}

make_registrationsh_file <- function(outfile = "registration.sh",path,movingobjects, targetobjects, paramfiles){
  file.create(outfile)
  cat(path, "paramDiffeos.xml",file = outfile, sep = " ")
  for (object in 1:length(movingobjects)){
    cat(paramfiles[object],movingobjects[object],targetobjects[object], file = outfile, append = T, sep = " ")
  }
}

# Call the sparseMatching3 with optional arguments
sparseMatching3<-function(..., verbose=FALSE){
  app=sparseMatching3("sparseMatching3")
  args=paste(..., collapse=" ")
  system(paste(app, args), ignore.stdout=!verbose)
}

# Find sparsematching3
find_sparseMatching3<-function(app){
  deformetricar.bindir=getOption("deformetricar.bindir")
  if(!is.null(deformetricar.bindir)) {
    app=file.path(path.expand(deformetricar.bindir), app)
  } else {
    # is it in the path?
    app=Sys.which('sparseMatching3')
    if(nchar(app)==0){
      # is it in a standard location?
      app="/Applications/deformetrica-2.1/deformetrica/bin/sparseMatching3"
      if(!file.exists(app))
        stop("Sorry I can't find the deformetrica ",basename(app)," executable! ",
             'Please add it to your PATH or set ',
             'options(deformetricar.bindir="/path/to/deformetrica/apps")')
    }
    options(deformetricar.bindir=dirname(app))
  }
  app
}

