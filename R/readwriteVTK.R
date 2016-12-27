#' Read a VTK format file
#'
#' @param filename The path to the file on disk
#' @param item The element(s) within the file to read (defaults to points)
#'
#' @return A matrix of points, indices (polygons) or normals
#' @export
read.vtk<-function(filename, item = c("points","triangles", "normals")){
  item=match.arg(item)

  if(!file.exists(filename)) stop("Cannot read: ",filename)
  con=file(filename,open='rb',encoding='ASCII')
  on.exit(close(con))
  magic=readLines(con,n=1)
  if(regexpr("# vtk DataFile Version [234]",magic,ignore.case =T)<0)
    stop("Bad header line in file: ",filename)

  title=readLines(con,1)
  encoding=readLines(con,1)
  if(regexpr("ASCII",encoding,ignore.case=TRUE)<0)
    stop("Can only read ASCII encoded VTK pointsets")

  datasetLine=toupper(readLines(con,1))
  if(regexpr("^DATASET",datasetLine)<0)
    stop("Missing DATASET line")

  datasetType=sub("DATASET\\s+(\\w+)","\\1",datasetLine)

  validDatasetTypes<-c("STRUCTURED_POINTS", "STRUCTURED_GRID",
                       "UNSTRUCTURED_GRID", "POLYDATA", "RECTILINEAR_GRID", "FIELD")

  if(!datasetType%in%validDatasetTypes)
    stop(datasetType," is not a valid VTK dataset type")
  if(datasetType!="POLYDATA")
    stop("ReadVTKLandmarks can currently only read POLYDATA.",
         " See http://www.vtk.org/VTK/img/file-formats.pdf for details.")

  pointsLine=toupper(readLines(con,1))
  if(regexpr("POINTS",pointsLine)<0)
    stop("Missing POINTS definition line")
  ptinfo=unlist(strsplit(pointsLine,"\\s+",perl=TRUE))
  if(length(ptinfo)!=3)
    stop("Unable to extract points information from POINTS line",pointsLine)
  nummarkers=as.integer(ptinfo[2])
  if(is.na(nummarkers))
    stop("Unable to extract number of points from POINTS line:",pointsLine)
  datatype=ptinfo[3]
  if(!datatype%in%toupper(c("unsigned_char", "char", "unsigned_short", "short", "unsigned_int", "int",
                            "unsigned_long", "long", "float", "double")))
    stop("Unrecognised VTK datatype: ",datatype)

  points=scan(con,what=1.0,n=3*nummarkers,quiet=TRUE)

  # VTK seems to be hardcoded for 3D
  if (item == "points"){
    m=matrix(points,ncol=3,byrow=T)
    colnames(m)=c("X","Y","Z")
    attr(m,"file")=filename
    attr(m,"title")=title
    attr(m,"vtk_datatype")=datatype
  }
  if (item != "points"){
    triangLine=toupper(readLines(con,1))
    if(length(triangLine)==0){
      warning("No data on polygons found")
      return(NULL)
    }
    if(regexpr("POLYGONS",triangLine)<0)
      stop("Missing POLYGONS definition line")
    lninfo=unlist(strsplit(triangLine,"\\s+",perl=TRUE))
    if(length(lninfo)!=3)
      stop("Unable to extract connection information from POLYGONS line:",triangLine)
    nummconns=as.integer(lninfo[2])
    if(is.na(nummconns))
      stop("Unable to extract number of connections from POLYGONS line:",triangLine)
    datatype=lninfo[3]
    triang=scan(con,what=1.0,n=4*nummconns,quiet=TRUE)
  }
  if (item == "triangles"){
    m=matrix(triang,ncol=4,byrow=T)[,-1]
    attr(m,"file")=filename
    attr(m,"title")=title
    attr(m,"vtk_datatype")= "int"
  }
  if (item == "normals"){
    normalsLine=toupper(readLines(con,1))
    if(length(normalsLine)==0){
      warning("No data on Normals found")
      return(NULL)
    }
    if(regexpr("NORMALS",normalsLine)<0)
      stop("Missing NORMALS definition line")
    ninfo=unlist(strsplit(normalsLine,"\\s+",perl=TRUE))
    if(length(ninfo)!=3)
      stop("Unable to extract connection information from POLYGONS line",triangLine)
    datatype=ninfo[3]
    if(!datatype%in%toupper(c("unsigned_char", "char", "unsigned_short", "short", "unsigned_int", "int",
                              "unsigned_long", "long", "float", "double")))
      stop("Unrecognised VTK datatype: ",datatype)
    normals=scan(con,what=1.0,n=3*nummarkers,quiet=TRUE)
    m=matrix(normals,ncol=3,byrow=T)
    attr(m,"file")=filename
    attr(m,"title")=title
    attr(m,"vtk_datatype")=datatype
  }
  m
}


#' Write VTK file
#'
#' @param points 3D coordinates (Nx3 matrix)
#' @param filename  Path to output file
#' @param polygons Triangle indices (Mx3 matrix)
#' @param normals Normals (Nx3 matrix)
#' @param datatype .vtk datatype (defaults to float)
#' @param title Title of the .vtk file (defaults to filename)
#' @param mesh mesh3d or surf object
#' @param trafo a transformation that can be given to Morpho::applyTransform
#' @export
#'
write.vtk <-function(points, filename, polygons = NULL, normals = NULL,
                     datatype=c("float","double"), title = filename){
  file.create(filename)
  if(ncol(points)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(points)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())

  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)

  write.table(points,col.names=F,row.names=F,file=filename,append=TRUE)
  if(!is.null(polygons)){
    if(ncol(polygons)!=3) stop("Expect N rows x 3 cols for polygons")
    # if(any(0%in%polygons)  == FALSE) { polygons = polygons -1 }# VTK files are 0 indexed
    numpoints = rep(3, nrow(polygons))
    mx = cbind(numpoints, polygons)
    cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
    write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
  }

  if(!is.null(normals)){
    cat(paste("NORMALS","clean",datatype),sep="\n",file=filename, append = TRUE)
    if(ncol(normals)!=3) stop("Expect N rows x 3 cols for normals")
    write.table(normals,col.names=F,row.names=F,file=filename,append=TRUE)
  }
  return("complete")
}

#' @export
write.vtk.mesh3d <- function(mesh, filename){
  write.vtk(nat::xyzmatrix(mesh), filename = filename, polygons = t(mesh$it)-1)
}

#' @export
write.vtk.surf<- function(mesh, filename, trafo = NULL){
  if (is.null(trafo)){
    write.vtk(nat::xyzmatrix(mesh$Vertices), filename = filename, polygons = mesh$Regions[[1]]-1)
  }else{
    write.vtk(Morpho::applyTransform(nat::xyzmatrix(mesh$Vertices), trafo), filename = filename, polygons = mesh$Regions[[1]]-1)
  }
}

#' Read txt files
#'
#' @param filename read CPS or MOM .txt file
#'
#'
read.points<-function(filename){
  if(!file.exists(filename)) stop("Cannot read: ",filename)
  con=file(filename,open='rb',encoding='txt')
  on.exit(close(con))
  magic=readLines(con,n=0)
  points=scan(con,what=1.0,quiet=TRUE)
  m=matrix(points,ncol=3,byrow=T)
  colnames(m)=c("X","Y","Z")
  attr(m,"file")=filename
  attr(m,"title")=filename
  m
}

#'  Transform a .vtk object
#'
#' @param vtk path to vtk object
#' @param transformations transformation(s)
#'
transform.vtk = function (vtk, transformations){
  positions = xyzmatrix(read.vtk(vtk, item = "points"))
  if (is.list(transformations) == F){
    cat("Single transformation")
    positions <- xform(positions, transformations)
  }
  if (is.list(transformations) == T){
    for (transformation in transformations){
      positions <- nat::xform(positions, transformation)
    }
  }
  return (positions)
}


#'  Read a mesh3d object from a .vtk file
#'
#' @param filename a list of filenames to covnert into a single mesh3d object
#'
read.vtk.to.mesh <- function (filenames){
  vertices = c()
  indices = c()
  for (file in filenames){
    vtk = read.vtk(file, "points")
    triang = read.vtk(file, "triangles") + 1
    indices = rbind(indices, triang)
    vertices = rbind(vertices, vtk)
  }
  vertices = cbind(vertices, 1)
  rgl::tmesh3d(t(vertices), t(indices), homogeneous = TRUE)
}

#' Write a mesh3d object, removing internal points
#'
#' @param Whereto save the object
#' @param Value of alpha to be used with alphashape3d::ashape3d() to re-define the surface
#' @param mesh3d a mesh3d object
write.vtk.light <- function(mesh3d, alpha = 2, filename){
  mesh3d.d = t(mesh3d$vb)[,-4][mesh3d$it,]
  mesh3d.alpha = alphashape3d::ashape3d(unique(mesh3d.d), alpha)
  triangles = mesh3d.alpha$triang[apply(mesh3d.alpha$triang, 1, function(x) {( any(as.numeric(x[9]) > 0))} ),][,1:3]
  newmesh = rgl::tmesh3d(t(mesh3d.alpha$x), t(triangles), homogeneous = F)
  write.vtk(nat::xyzmatrix(newmesh), filename = filename, polygons = t(newmesh$it)-1)
  return(newmesh)
}


#' @export
write.vtk.mesh <- function(mesh, filename, trafo = NULL){
  if (!is.null(trafo)){
    if(nat::is.amiramesh(mesh)){write.vtk(Morpho::applyTransform(nat::xyzmatrix(mesh$Vertices),trafo), filename = filename, polygons = mesh$Regions[[1]]-1)
    }else write.vtk(Morpho::applyTransform(nat::xyzmatrix(mesh),trafo), filename = filename, polygons = t(mesh$it)-1)
  }else if(nat::is.amiramesh(mesh)){write.vtk(nat::xyzmatrix(mesh$Vertices), filename = filename, polygons = mesh$Regions[[1]]-1)
  }else write.vtk(nat::xyzmatrix(mesh), filename = filename, polygons = t(mesh$it)-1)
}


write.vtk.neurons<-function(someneuronlist,filename,datatype=c("float","double"), connect = T){
  title = filename
  d = xyzmatrix(someneuronlist)
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())

  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)

  for (neuron in someneuronlist){
    s = xyzmatrix(neuron)
    write.table(s,col.names=F,row.names=F,file=filename,append=TRUE)
    # Then add a space to identify neurons later?
  }

  if (connect == T){ # Include the connections between points in each neuron in the VTK file
    count = 0
    nopolypoints = nrow(d) - length(someneuronlist)
    cat(paste("LINES",nopolypoints,nopolypoints*3),sep="\n",file=filename, append = TRUE)
    for (neuron in someneuronlist){
      r = 0
      rows = nrow(xyzmatrix(neuron))
      mx = matrix(, nrow = rows -1, ncol = 3)
      segments = neuron$SegList
      for(segment in segments){
        for (n in 1:(length(segment)-1)){
          r = r + 1
          mx[r,] <- c(2, (as.integer(segment[n]) -1 + count), (as.integer(segment[n+1]) -1 + count))
        }
      }
      count = count + rows
      write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
    }
  }
}
