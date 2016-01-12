context("VTK IO")
test_that("we can read vtk files",{
  vtkfile=system.file("extdata/mushroom_body_right.vtk", package = 'deformetricar')
  expect_is(m<-read.vtk(vtkfile), 'matrix')
  expect_equal(ncol(m), 3)
  expect_equal(colnames(m), c("X", "Y", "Z"))
  expect_equal(nrow(m), 3088)
})
