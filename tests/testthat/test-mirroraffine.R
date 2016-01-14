context("mirroraffine")
test_that("we can perform a mirroring and affine transformation",{
  n3=nat::Cell07PNs[1:3]
  expect_is(t.n3<-apply.mirror.affine(n3), 'neuronlist')
})
