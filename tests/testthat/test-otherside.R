context("otherside")
test_that("Hello from the other side...",{
  n3=nat::Cell07PNs[1:3]
  expect_is(t.n3<-otherside(n3), 'neuronlist')
})
