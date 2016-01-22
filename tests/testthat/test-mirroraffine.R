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
  baseline=matrix(c(12.5957462419454, 11.7470158497127, 10.89828545748,
                    10.0495550652473, 9.20082467301456, 97.8315721943438, 98.6296324297916,
                    99.4276926652394, 100.225752900687, 101.023813136135, 106.200027839867,
                    107.070026731711, 107.940025623554, 108.810024515398, 109.680023407242),
                  ncol=3)
  expect_equal(transform3dpoints(m, fullmirror), baseline)
})
