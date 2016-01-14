context("shoot and flow")
test_that("we can run shoot and flow on a matrix",{
  vtkfile=system.file("extdata/mushroom_body_right.vtk", package = 'deformetricar')
  m<-read.vtk(vtkfile)
  regdir=system.file("extdata/reg_output", package = 'deformetricar')
  expect_is(deform.m<-shootflow(m, regdir), 'matrix')

  first3_baseline=matrix(c(26.0168, 28.79, 26.0277, 28.3175, 24.1387, 33.2848,
                           41.1493, 45.1458, 35.4634), nrow=3,
                         dimnames = list(NULL, c("X", "Y", "Z")))
  expect_equal(deform.m[1:3,], first3_baseline)
})

test_that("we can run shoot and flow on a neuronlist",{
  n3=nat::Cell07PNs[1:3]
  regdir=system.file("extdata/reg_output", package = 'deformetricar')
  expect_is(deform.n3<-shootflow(n3, regdir), 'neuronlist')
})
