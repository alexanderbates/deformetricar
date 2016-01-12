# shoot and flow functions
# idea is that we want to be able to apply an existing registratino to
# an abitrary set of coordinates

# see http://www.deformetrica.org/?page_id=232

# first thing to achieve support for bare vertices (3d coordinates)

shootflow<-function(someneuronlist, regdir = getwd(), data.sigma = 1, kernel.width=5, object.type = "NonOrientedPolyLine"){
  # we need to make a command line like this
  # ShootAndFlow3 paramsDiffeos.xml Direction CP.txt Mom.txt paramsObject1.xml object1 paramsObject2.xml object2 …
  start = getwd()
  setwd(regdir)
  steps = as.integer(unlist(strsplit(gsub("[^0-9]", "", unlist(readLines("paramDiffeos.xml")[25])), "")))-1
  params_file=make_params_file(data.sigma = data.sigma, kernel.width = kernel.width, object.type = )
  WriteVTKNeurons("object.vtk", someneuronlist, "float")
  ShootAndFlow3("paramDiffeos.xml", 1, "CP_final.txt", "Mom_final.txt", params_file, "object.vtk")
  output.file = paste("object_flow__t_", steps, ".vtk", sep = "")
  deform.neurons = ReadVTKNeurons(output.file, someneuronlist)
  file.remove("object.vtk", params_file)
  for (i in 0:steps){
    file.remove(paste("object_flow__t_", i, ".vtk", sep = ""))
  }
  setwd(start)
  return(deform.neurons)
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
