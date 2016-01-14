context("VTK IO")
test_that("we can read vtk files",{
  vtkfile=system.file("extdata/vtk_files/larval_central_nervous_system.vtk", package = 'deformetricar')
  expect_is(m<-read.vtk(vtkfile), 'matrix')
  expect_equal(ncol(m), 3)
  expect_equal(colnames(m), c("X", "Y", "Z"))
  expect_equal(nrow(m), 3088)
})

test_that("we can write vtk files",{
  vtkfile=system.file("extdata/vtk_files/larval_central_nervous_system.vtk", package = 'deformetricar')
  m<-read.vtk(vtkfile)

  tf=tempfile(fileext = '.vtk')
  on.exit(unlink(tf))
  write.vtk(tf, points=m)
  n = read.vtk(tf)
  # nb don't check attributes since they include e.g. input file path
  expect_equivalent(m, n)
})
