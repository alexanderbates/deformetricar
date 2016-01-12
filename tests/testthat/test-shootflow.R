context("shoot and flow")
test_that("we can run shoot and flow on a matrix",{
  vtkfile=system.file("extdata/mushroom_body_right.vtk", package = 'deformetricar')
  m<-read.vtk(vtkfile)
  regdir=system.file("extdata/reg_output", package = 'deformetricar')
  expect_is(deform.m<-shootflow(m, regdir), 'matrix')
})

test_that("we can run shoot and flow on a neuonlist",{
  n3=nat::Cell07PNs[1:3]
  regdir=system.file("extdata/reg_output", package = 'deformetricar')
  expect_is(deform.n3<-shootflow(n3, regdir), 'neuronlist')
})
