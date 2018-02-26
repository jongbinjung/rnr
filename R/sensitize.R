#' Compute the sensitivity-adjusted estimates of predicted outcome given
#' treatment/control
#'
#' @param df data frame to analyze; must include columns $treat: Observed
#'   (binary) treatment, e.g., bail_set $resp_ctl: Predicted probability of
#'   positive resp given control, $resp_trt: Predicted probability of positive
#'   resp given treatment, $p_trt: predicted probability of treatment
#' @param q p(u = 1 | x)
#' @param dp change in log-odds of treat = 1 if u = 1
#' @param d0 change in log-odds of response = 1 if treat = 0 and u = 1
#' @param d1 change in log-odds of response = 1 if treat = 1 and u = 1
#' @param debug logical, whether or not to return columns of intermediate
#'   variables for debugging purposes
#'
#' @return A data frame with the columns resp_ctl and resp_trt updated according
#'   to the sensitivity parameters. If debug = TRUE, returned data frame will
#'   also contain columns of intermediate variables computed for sensitivity,
#'   appended with "__" (e.g., gamma__), with the original response estimates
#'   renamed as resp_trt_trt__ = resp_trt resp_ctl_ctl__ = resp_ctl
#' @export
#'
#' @examples
#' d <- data.frame(treat = 0, resp_ctl = .2, resp_trt = .3, p_trt = .5)
#' sensitize(d, q = .5, dp = log(2), d0 = log(2), d1 = log(2))
sensitize <- function(df, q, dp, d0, d1, debug = FALSE) {
  dep_cols <- c("treat", "resp_ctl", "resp_trt", "p_trt")

  MSG <- paste0("Missing column: %s in argument df\n",
                "RnR sensitivity requires a data frame with columns: [",
                paste0(dep_cols, collapse = ","), "]\n",
                "See function definition of sensitize() for details")

  purrr::walk(dep_cols, function(cn) {
    assertthat::assert_that(cn %in% colnames(df), msg = sprintf(MSG, cn))
  })

  # Rename columns to reflect that, under unobserved confounding,
  # E[resp=1 | treat=1] is actually resp(treat=1) | treat=1, and
  # E[resp=1 | treat=0] is actually resp(treat=0) | treat=0
  resp_ctl_ctl__ <- df$resp_ctl
  resp_trt_trt__ <- df$resp_trt

  # trt = 1 means bail = 1
  gamma__ = solve(prob = 1 - q,
                  delta = dp,
                  lhs = 1 - df$p_trt)

  # Compute p(trt|trt,u)
  ptrt_u1__ <- inv_logit(gamma__ + dp)
  ptrt_u0__ <- inv_logit(gamma__)

  # Compute p(u=0|a)
  # u = 1 is more likely to have positive responce (resp=1)
  pu0_ctl__ <- (1 - ptrt_u0__) * (1 - q) /
    ((1 - ptrt_u0__) * (1 - q) + (1 - ptrt_u1__) * q)
  pu0_trt__ <- ptrt_u0__ * (1 - q) / (ptrt_u0__ * (1 - q) + ptrt_u1__ * q)
  pu0__ <- ifelse(df$treat, pu0_trt__, pu0_ctl__)

  # lhs r_t = 1 means resp = 1
  beta_ctl__ <-
    solve(prob = pu0_ctl__,
          delta = d0,
          lhs = 1 - resp_ctl_ctl__)
  beta_trt__ <-
    solve(prob = pu0_trt__,
          delta = d1,
          lhs = 1 - resp_trt_trt__)

  # resp(trt=0) | trt=1
  resp_ctl_trt__ = pu0__ * inv_logit(beta_ctl__) +
    (1 - pu0__) * inv_logit(beta_ctl__ + d0)
  # resp(trt=1) | trt=0
  resp_trt_ctl__ = pu0__ * inv_logit(beta_trt__) +
    (1 - pu0__) * inv_logit(beta_trt__ + d1)

  # Sensitivity adjusted measure of resp(treat) | observed treatment
  df$resp_trt <- ifelse(df$treat, resp_trt_trt__, resp_trt_ctl__)
  df$resp_ctl <- ifelse(df$treat, resp_ctl_trt__, resp_ctl_ctl__)

  # Record params
  df$q = q
  df$dp=dp
  df$d0=d0
  df$d1=d1

  if (debug) {
    debug_df <- data.frame(
      resp_ctl_ctl__ = resp_ctl_ctl__,
      resp_ctl_trt__ = resp_ctl_trt__,
      resp_trt_trt__ = resp_trt_trt__,
      resp_trt_ctl__ = resp_trt_ctl__,
      gamma__ = gamma__,
      ptrt_u1__ = ptrt_u1__,
      ptrt_u0__ = ptrt_u0__,
      pu0_trt__ = pu0_trt__,
      pu0_ctl__ = pu0_ctl__,
      pu0__ = pu0__,
      beta_trt__ = beta_trt__,
      beta_ctl__ = beta_ctl__
    )

    df <- cbind(df, debug_df)
  }

  return(df)
}

inv_logit <- stats::binomial()$linkinv
