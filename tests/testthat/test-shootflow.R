context("shoot and flow")
test_that("we can run shoot and flow on a matrix",{
  vtkfile=system.file("extdata/mushroom_body_right.vtk", package = 'deformetricar')
  m<-read.vtk(vtkfile)
  regdir=system.file("extdata/reg_output", package = 'deformetricar')
  expect_is(deform.m<-shootflow(m, regdir), 'matrix')
})
