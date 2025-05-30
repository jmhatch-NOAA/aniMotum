context("test emf")

## emf builds a suitable data.frame for use in fit_ssm call
data(ellie)

test_that("f has s3 classes `ssm_df`, `sf`, `data.frame`", {
  f <- fit_ssm(ellie, vmax=5, model = "rw", time.step = 48, 
               control = ssm_control(optim = "nlminb", 
                                     verbose = 0), 
               emf = 
                 emf(emf.x = c(1, 1, 1, 1, 1, 10), 
                     emf.y  = c(1, 1, 1, 1, 1, 10)
                     )
               )
  expect_s3_class(f$ssm[[1]], "ssm")
  expect_equal(length(f$ssm[[1]]), 15)
  expect_s3_class(f, "ssm_df")
})