context("sensitize")

sample_df <- data.frame(
  id = 1,
  treat = c(0, 1),
  p_trt = c(.5, .5, .1, .1, .9, .9),
  resp_ctl = c(.2, .2),
  resp_trt = c(.3, .2)
)

test_that("sensitize() returns correct values for zero-effect params", {
  # Zero prevalence of u = 1
  generated <- sensitize(sample_df, q = 0, dp = 5, d0 = 5, d1 = 5)
  expect_equal(generated$resp_ctl, sample_df$resp_ctl)
  expect_equal(generated$resp_trt, sample_df$resp_trt)

  # Zero treatment effect
  generated <- sensitize(sample_df, q = .5, dp = 0, d0 = 5, d1 = 5)
  expect_equal(generated$resp_ctl, sample_df$resp_ctl)
  expect_equal(generated$resp_trt, sample_df$resp_trt)

  # Zero effect on outcome given control
  generated <- sensitize(sample_df, q = .5, dp = 5, d0 = 0, d1 = 5)
  expect_equal(generated$resp_ctl, sample_df$resp_ctl)

  # Zero effect on outcome given treatment
  generated <- sensitize(sample_df, q = .5, dp = 5, d0 = 5, d1 = 0)
  expect_equal(generated$resp_trt, sample_df$resp_trt)
})

test_that("sensitize() returns non-negative estimates for range of q", {
  qs <- seq(0, 1, 0.01)

  # Positive effects of u
  for (q in qs) {
    generated <- sensitize(sample_df, q = q, dp = 5, d0 = 5, d1 = 5)
    expect_gte(min(generated$resp_ctl), 0)
    expect_gte(min(generated$resp_trt), 0)
  }

  # Negative effects of u
  for (q in qs) {
    generated <- sensitize(sample_df, q = q, dp = -5, d0 = -5, d1 = -5)
    expect_gte(min(generated$resp_ctl), 0)
    expect_gte(min(generated$resp_trt), 0)
  }
})
