context("mirroraffine")
test_that("we can perform a mirroring and affine transformation",{
  n3=nat::Cell07PNs[1:3]
  file=system.file("extdata/fullmirror.rds", package = 'deformetricar')
  expect_is(t.n3<-apply.mirror.affine(n3, pathtomatrix = file), 'neuronlist')
})
