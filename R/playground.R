# Playground

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

regdir = "/Users/abates/Desktop/testreg3/"

# Set default neuronlist as the full data from catmaid
fulln = convert(readRDS('fulln-2015-12-14.rds'))
options('nat.default.neuronlist'='fulln')

l1mesh$vb <- l1mesh$vb/1000

moms = ReadPoints("/Users/abates/Desktop/testreg3/MOM_final.txt")
cps = ReadPoints("/Users/abates/Desktop/testreg3/CP_final.txt")
finals = read.vtk("/Users/abates/Desktop/testreg3/finals.vtk")
finalslines = read.vtk("/Users/abates/Desktop/testreg3/finalsline.vtk")
target = cps+(moms)

target1 = cps+(moms)
target2 = cps+(2*moms)
target3 = cps+(5*moms)
target4 = cps+(10*moms)
target5 = cps+(20*moms)

c = apply.mirror.affine(l1, pathtomatrix = "/Users/abates/Desktop/testreg3/fullmirror.rds")
c = tps3d(c, cps, finals)



# Visualise deformation field.



vector_field <- function(
  f,  # Function describing the vector field
  xmin=0, xmax=1, ymin=0, ymax=1,
  width=600, height=600,
  iterations=50,
  epsilon=.01,
  trace=TRUE
) {
  z <- matrix(runif(width*height),nr=height)
  i_to_x <- function(i) xmin + i / width  * (xmax - xmin)
  j_to_y <- function(j) ymin + j / height * (ymax - ymin)
  x_to_i <- function(x) pmin( width,  pmax( 1, floor( (x-xmin)/(xmax-xmin) * width  ) ) )
  y_to_j <- function(y) pmin( height, pmax( 1, floor( (y-ymin)/(ymax-ymin) * height ) ) )
  i <- col(z)
  j <- row(z)
  x <- i_to_x(i)
  y <- j_to_y(j)
  res <- z
  for(k in 1:iterations) {
    v <- matrix( f(x, y), nc=2 )
    x <- x+.01*v[,1]
    y <- y+.01*v[,2]
    i <- x_to_i(x)
    j <- y_to_j(y)
    res <- res + z[cbind(i,j)]
    if(trace) {
      cat(k, "/", iterations, "\n", sep="")
      dev.hold()
      image(res)
      dev.flush()
    }
  }
  if(trace) {
    dev.hold()
    image(res>quantile(res,.6), col=0:1)
    dev.flush()
  }
  res
}

# Sample data
van_der_Pol <- function(x,y, mu=1) c(y, mu * ( 1 - x^2 ) * y - x )
res <- vector_field(
  van_der_Pol,
  xmin=-3, xmax=3, ymin=-3, ymax=3,
  width=800, height=800,
  iterations=50,
  epsilon=.01
)
image(-res)


# Write control points to a .vtk file
write.vtk(cps, "latticecps.vtk")





x = l1
points3d(x)
x = apply.mirror.affine(x)
x = tps3d(x, cps, finals, lambda = 0)
points3d(x, col = 'red')
x = xyzmatrix(fulln[1])/1000
x = apply.mirror.affine(x)
x = shootflow(x)
points3d(x, col = 'orange')



