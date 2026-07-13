skip_if_no_deformetrica <- function() {
  ok <- tryCatch({ find_deformetrica(); TRUE }, error = function(e) FALSE)
  testthat::skip_if_not(ok, "deformetrica (>= 4.3) executable not found")
}

# Running deformetrica needs its (conda) environment activated for the torch/CUDA
# backend, which a bare `R CMD check` process does not provide — so live-compute
# tests are opt-in via DEFORMETRICAR_TEST_LIVE=1.
skip_if_no_live_deformetrica <- function() {
  skip_if_no_deformetrica()
  testthat::skip_if_not(nzchar(Sys.getenv("DEFORMETRICAR_TEST_LIVE")),
                        "set DEFORMETRICAR_TEST_LIVE=1 to run live Deformetrica compute tests")
}

test_that("write_vtk / read_vtk round-trips a point cloud (nat neuron data)", {
  skip_if_not_installed("nat")
  pts <- nat::xyzmatrix(nat::Cell07PNs[[1]])
  f <- tempfile(fileext = ".vtk"); on.exit(unlink(f))
  expect_equal(write_vtk(pts, f), "complete")
  back <- read_vtk(f, item = "points")
  expect_equal(nrow(back), nrow(pts))
  expect_equal(back, pts, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("write_vtk / read_vtk round-trips a triangular surface mesh", {
  skip_if_not_installed("Rvcg")
  m <- Rvcg::vcgSphere()                       # mesh3d: 3 x V verts, 3 x F faces
  v <- t(m$vb[1:3, ]); faces <- t(m$it)
  f <- tempfile(fileext = ".vtk"); on.exit(unlink(f))
  expect_equal(write_vtk(v, f, polygons = faces), "complete")
  expect_equal(nrow(read_vtk(f, item = "points")), nrow(v))
  tri <- read_vtk(f, item = "triangles")
  expect_equal(nrow(tri), nrow(faces))
})

test_that("find_deformetrica errors informatively when nothing is available", {
  withr::local_options(deformetricar.exe = "/no/such/deformetrica/binary")
  # Only assert the error path when there is genuinely no deformetrica resolvable on
  # this host (PATH, a managed reticulate/conda env, ...) -- find_deformetrica() knows
  # all its own search locations, so ask it rather than re-enumerate them here.
  if (!is.null(tryCatch(find_deformetrica(), error = function(e) NULL)))
    skip("a real deformetrica is resolvable on this host")
  expect_error(find_deformetrica(), "Cannot find")
})

test_that("deformetrica_register + deformetrica_shoot recover a known translation", {
  skip_if_no_live_deformetrica()
  skip_if_not_installed("nat")
  set.seed(1)
  src <- nat::xyzmatrix(nat::kcs20[[1]])
  src <- src[seq_len(min(60L, nrow(src))), , drop = FALSE]
  shift <- matrix(rep(c(3, -2, 1), each = nrow(src)), ncol = 3)  # a small rigid-ish move
  tgt <- src + shift
  # kernel_width must be smaller than this 60-point subset's ~16um extent, or the fit
  # lays down a single control point that Deformetrica cannot shoot (see the guard in
  # deformetrica_shoot()).
  fit <- deformetrica_register(src, tgt, kernel_width = 8, max_iterations = 40, device = "auto")
  expect_true(file.exists(fit$control_points) && file.exists(fit$momenta))
  warped <- deformetrica_shoot(src, fit$control_points, fit$momenta,
                               kernel_width = fit$kernel_width, device = "auto")
  expect_equal(dim(warped), dim(src))
  # the warp should move src markedly TOWARD tgt (not necessarily exactly)
  d_before <- mean(sqrt(rowSums((src    - tgt)^2)))
  d_after  <- mean(sqrt(rowSums((warped - tgt)^2)))
  expect_lt(d_after, d_before)
})

test_that("deformetrica_register returns a deformetricareg usable as a nat transform", {
  skip_if_no_live_deformetrica()
  skip_if_not_installed("nat")
  set.seed(1)
  src <- matrix(stats::rnorm(60), 20, 3) * 10
  tgt <- sweep(src, 2, c(6, 0, 0))                       # a pure translation
  reg <- deformetrica_register(src, tgt, kernel_width = 9, max_iterations = 40, device = "auto")
  expect_s3_class(reg, "deformetricareg")

  p <- matrix(stats::rnorm(15), 5, 3) * 8
  # nat::xform(x, reg) dispatches to xformpoints.deformetricareg == deformetrica_shoot()
  expect_equal(nat::xform(p, reg), deformetrica_shoot(p, reg), tolerance = 1e-5)
  expect_s3_class(nat::xform(nat::kcs20[1:2], reg), "neuronlist")

  # a saved reg survives its tempfiles being cleaned (contents are stored inline)
  a  <- nat::xform(p, reg)
  tf <- tempfile(fileext = ".rds"); saveRDS(reg, tf)
  unlink(reg$output_dir, recursive = TRUE)
  unlink(c(reg$control_points, reg$momenta))
  expect_equal(nat::xform(p, readRDS(tf)), a, tolerance = 1e-5)
})
