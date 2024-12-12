test_that("Coords has correct dimension", {
  n <- 100
  x <- runif(n)
  y <- runif(n)

  types <- sample(c("Type1","Type2"),size = n, replace = TRUE)
  image_ids <- sample(c("Image1","Image2"),size = n, replace = TRUE)
  nx <- ny <- 20
  out <- pixellate_grid(x,y,types,image_ids,nx,ny)

  expect_equal(nrow(out$coords),nx*ny)
})
