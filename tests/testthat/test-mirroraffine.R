context("mirroraffine")
test_that("we can perform a mirroring and affine transformation",{
  app=try(find_app("ShootAndFlow3"))
  if(inherits(app, 'try-error')) skip("Can't find ShootAndFlow3!")

  n3=nat::Cell07PNs[1:3]
  file=system.file("extdata/fullmirror.rds", package = 'deformetricar')
  expect_is(t.n3<-apply.mirror.affine(n3, pathtomatrix = file), 'neuronlist')
})

test_that("transform3d can apply multiple transforms",{
  m=matrix(101:115, nrow=5, ncol=3)
  fullmirror=readRDS(system.file("extdata/fullmirror.rds", package = 'deformetricar'))
  baseline=structure(c(11.5662779431759, 10.702386783413, 9.83849562365012,
                       8.97460446388725, 8.11071330412435, 98.6329065200601, 99.4694416625517,
                       100.305976805043, 101.142511947535, 101.979047090026, 103.344312853856,
                       104.185382590343, 105.026452326829, 105.867522063315, 106.708591799802
  ), .Dim = c(5L, 3L))
  expect_equal(transform3dpoints(m, fullmirror), baseline)
})
