ReadVTKLandmarks<-function(filename, item = "points"){
  if(!file.exists(filename)) stop("Cannot read: ",filename)
  con=file(filename,open='rb',encoding='ASCII')
  on.exit(close(con))
  magic=readLines(con,n=1)
  if(regexpr("# vtk DataFile Version [23]",magic,ignore=T)<0)
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
      print("yo")
      warning("No data on polygons found")
      return(NULL)
    }  
    if(regexpr("POLYGONS",triangLine)<0)
      stop("Missing POLYGONS definition line")
    lninfo=unlist(strsplit(triangLine,"\\s+",perl=TRUE))
    if(length(lninfo)!=3)
      stop("Unable to extract connection information from POLYGONS line",linesLine)
    nummconns=as.integer(lninfo[2])
    if(is.na(nummconns))
      stop("Unable to extract number of connections from POLYGONS line:",linesLine)
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
      stop("Unable to extract connection information from POLYGONS line",linesLine)
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


WriteVTKLandmarks<-function(filename,d,title,datatype=c("float","double")){
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
}

WriteVTKMesh<-function(filename,d,title = filename,datatype=c("float","double"), surf = NULL, ashape = NULL){
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
  
  if (is.null(surf) == F){
    mx = data.matrix(surf$Regions$Inside[, 1:3] - 1) # VTK files are 0 indexed
    numpoints = rep(3, nrow(mx))
    mx = cbind(numpoints, mx)
    cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
    write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
  }
  
  if (is.null(ashape) == F){
    data = ashape$triang
    keeps <- apply(data, 1, function(x) {( any(as.numeric(x[9]) > 1))} ) # Removes rows for triangles not included in the alphashape, for the chosen alpha. Includes other simplexes: interior, regular and singular.
    mx = data.matrix(data[keeps,][,1:3]-1) # VTK files are 0 indexed
    numpoints = rep(3, nrow(mx))
    mx = cbind(numpoints, mx)
    cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
    write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
  }
}

WriteVTKalphashapecombined <-function(filename, title = filename, datatype=c("float","double"), ashapelist){
  if(ncol(ashapelist[[1]]$x)!=3) stop("Expect N rows x 3 cols of 3d points")
  
  initial = ashapelist[[1]]
  positions = initial$x
  triangles = initial$triang
  for (a in 2:length(ashapelist))
  {
    count = nrow(positions)
    ashape = ashapelist[[a]]
    positions = rbind(positions, ashape$x)
    ashape$triang[,1:3] <- ashape$triang[,1:3] + count
    triangles = rbind(triangles, ashape$triang)
  }
  d = positions
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)

  data = triangles
  keeps <- apply(data, 1, function(x) {( any(as.numeric(x[9]) > 1))} ) # Removes rows for triangles not included in the alphashape, for the chosen alpha. Includes other simplexes: interior, regular and singular.
  mx = data.matrix(data[keeps,][,1:3]-1) # VTK files are 0 indexed
  numpoints = rep(3, nrow(mx))
  mx = cbind(numpoints, mx)
  cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
  write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
}

WriteVTKalphashape <-function(filename, title = filename, datatype=c("float","double"), ashape){
  d = ashape$x
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
  
  data = ashape$triang
  keeps <- apply(data, 1, function(x) {( any(as.numeric(x[9]) > 1))} ) # Removes rows for triangles not included in the alphashape, for the chosen alpha. Includes other simplexes: interior, regular and singular.
  mx = data.matrix(data[keeps,][,1:3]-1) # VTK files are 0 indexed
  numpoints = rep(3, nrow(mx))
  mx = cbind(numpoints, mx)
  cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
  write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
}




WriteVTKSpines<-function(filename,someneuronlist,title,datatype=c("float","double"), connect = F, is.spine = F, trim = F){
  if (is.spine == F) {spines = nlapply(someneuronlist, nat::spine, UseStartPoint = T, .progress = 'text')}
  else{spines = someneuronlist}
  d = xyzmatrix(spines)
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  for (spine in spines){
    s = xyzmatrix(spine)
    write.table(s,col.names=F,row.names=F,file=filename,append=TRUE)
    # Then add a space to identify neurons later?
  }
  
  if (connect == T){ # Include the connections between points in each neuron in the VTK file
    count = 0
    nopolypoints = nrow(d) - length(spines)
    cat(paste("LINES",nopolypoints,nopolypoints*3),sep="\n",file=filename, append = TRUE)
    for (spine in spines){
      r = 0
      rows = nrow(xyzmatrix(spine))
      mx = matrix(, nrow = rows -1, ncol = 3)
      segments = spine$SegList
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

WriteVTKSpines<-function(filename,someneuronlist,title,datatype=c("float","double"), connect = F, is.spine = F, trim = F){
  if (is.spine == F) {spines = nlapply(someneuronlist, nat::spine, UseStartPoint = T, .progress = 'text')}
  else{spines = someneuronlist}
  d = xyzmatrix(spines)
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  for (spine in spines){
    s = xyzmatrix(spine)
    write.table(s,col.names=F,row.names=F,file=filename,append=TRUE)
    # Then add a space to identify neurons later?
  }
  
  if (connect == T){ # Include the connections between points in each neuron in the VTK file
    count = 0
    nopolypoints = nrow(d) - length(spines)
    cat(paste("LINES",nopolypoints,nopolypoints*3),sep="\n",file=filename, append = TRUE)
    for (spine in spines){
      r = 0
      rows = nrow(xyzmatrix(spine))
      mx = matrix(, nrow = rows -1, ncol = 3)
      segments = spine$SegList
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


ReadVTKNeurons <- function(filename, someneuronlist = NULL){
  if(!file.exists(filename)) stop("Cannot read: ",filename)
  con=file(filename,open='rb',encoding='ASCII')
  on.exit(close(con))
  magic=readLines(con,n=1)
  if(regexpr("# vtk DataFile Version [234]",magic,ignore=T)<0)
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
  
  # VTK seems to be hardcoded for 3D
  points=scan(con,what=1.0,n=3*nummarkers,quiet=TRUE) 
  m=matrix(points,ncol=3,byrow=T)
  colnames(m)=c("X","Y","Z")
  if (is.null(someneuronlist) == F){  
    count = 1
    for (neuron in 1:length(someneuronlist)){
      n = count + nrow(xyzmatrix(someneuronlist[[neuron]])) -1
      someneuronlist[[neuron]]$d$X <- m[count:n,1]
      someneuronlist[[neuron]]$d$Y <- m[count:n,2]
      someneuronlist[[neuron]]$d$Z <- m[count:n,3]
      count = count + nrow(xyzmatrix(someneuronlist[[neuron]]))
    }
  }
  
  linesLine=toupper(readLines(con,1))
  if(regexpr("LINES",linesLine)<0)
    stop("Missing LINES definition line")
  lninfo=unlist(strsplit(linesLine,"\\s+",perl=TRUE))
  if(length(lninfo)!=3)
    stop("Unable to extract connection information from LINES line",linesLine)
  nummarkers=as.integer(lninfo[2])
  if(is.na(nummarkers))
    stop("Unable to extract number of connections from LINES line:",linesLine)
  
  return(someneuronlist)
}

WriteVTKNeurons<-function(filename,someneuronlist,datatype=c("float","double"), connect = T, is.spine = F, trim = F){
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


WriteVTKmesh3d <-function(filename, title = filename, datatype=c("float","double"), mesh3d){
  d = t(mesh3d$vb)[,-4]
  if(ncol(d)!=3) stop("Expect N rows x 3 cols of 3d points")
  nummarkers=nrow(d)
  datatype=match.arg(datatype)
  if(missing(title)) title=paste("Data written from R by WriteVTKLandmarks at",Sys.time())
  
  cat("# vtk DataFile Version 2.0",
      title,
      "ASCII",
      "DATASET POLYDATA",
      paste("POINTS",nummarkers,datatype),sep="\n",file=filename)
  
  write.table(d,col.names=F,row.names=F,file=filename,append=TRUE)
  mx = t(mesh3d$it)-1 # VTK files are 0 indexed
  numpoints = rep(3, nrow(mx))
  mx = cbind(numpoints, mx)
  cat(paste("POLYGONS",nrow(mx),nrow(mx)*4),sep="\n",file=filename, append = TRUE)
  write.table(mx,col.names=F,row.names=F,file=filename,append=TRUE)
  if (!is.null(mesh3d$normals)){ 
    cat(paste("NORMALS","clean",datatype),sep="\n",file=filename, append = TRUE)
    normals = t(mesh3d$normals)[,-4]
    write.table(normals,col.names=F,row.names=F,file=filename,append=TRUE)
  }
}

WriteVTK <-function(filename, title = filename, points, polygons = NULL, normals = NULL, datatype=c("float","double")){
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

ReadPoints<-function(filename){
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
