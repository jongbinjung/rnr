#' Compute the sensitivity-adjusted estimates of predicted outcome given
#' treatment/control
#'
#' @param df data frame to analyze; must include columns $id: column used to
#'   identify each row (for future debugging) $resp: Observed (binary) response,
#'   e.g., FTA $treat: Observed (binary) treatment, e.g., bail_set $resp_ctl:
#'   Predicted probability of positive resp given control, e.g., E[FTA(bail=0) |
#'   bail = 0] $resp_trt: Predicted probability of positive resp given
#'   treatment, e.g., E[FTA(bail=1) | bail = 1] $p_trt: predicted probability of
#'   setting bail for observed x
#' @param q p(u = 1 | x)
#' @param dp change in log-odds of treat = 1 if u = 1
#' @param d0 change in log-odds of resp = 1 if treat = 0 and u = 1
#' @param d1 change in log-odds of resp = 1 if treat = 1 and u = 1
#' @param debug logical, whether or not to return columns of intermediate
#'   variables for debuging purposes
#'
#' @return A data frame with the columns resp_ctl and resp_trt updated according
#'   to the sensitivity parameters.
#'   If debug = TRUE, returned data frame will also contain columns of
#'   intermediate variables computed for sensitivity, starting with "v_rnr_"
#'   (e.g., v_rnr_gamma), with the original response estimates renamed as
#'      v_rnr_resp_trt_trt = resp_trt
#'      v_rnr_resp_ctl_ctl = resp_ctl
#' @export
#'
#' @note If debug = FALSE (default), all columns that start with "v_rnr_" will
#'   be removed from the final result
#' @examples
#' d <- data.frame(id = 1, resp = 1, treat = 0,
#'                 resp_ctl = .2, resp_trt = .3, p_trt = .5)
#' sensitize(d, q = .5, dp = log(2), d0 = log(2), d1 = log(2))
sensitize <- function(df, q, dp, d0, d1, debug = FALSE) {
  dep_cols <- c("id", "resp", "treat", "resp_ctl", "resp_trt", "p_trt")

  MSG <- paste0("Missing column: %s in argument df\n",
                "RnR sensitivity requires a data frame with columns: [",
                paste0(dep_cols, collapse = ","), "]\n",
                "See function definition of sensitize() for details")

  purrr::walk(dep_cols, function(cn) {
         assertthat::assert_that(cn %in% colnames(df), msg = sprintf(MSG, cn))
  })

  ret <- df %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      # Rename columns to reflect that, under unobserved confounding,
      # E[resp=1 | treat=1] is actually resp(treat=1) | treat=1, and
      # E[resp=1 | treat=0] is actually resp(treat=0) | treat=0
      v_rnr_resp_ctl_ctl = resp_ctl,
      v_rnr_resp_trt_trt = resp_trt) %>%
    dplyr::mutate(# Compute gamma
      # trt = 1 means bail = 1
      v_rnr_gamma = solve(prob = 1-q, delta = dp, lhs = 1-p_trt),

      # Compute p(trt|trt,u)
      v_rnr_ptrt_u1 = inv_logit(v_rnr_gamma + dp),
      v_rnr_ptrt_u0 = inv_logit(v_rnr_gamma),

      # Compute p(u=0|a)
      # u = 1 is more likely to have positive responce (resp=1)
      v_rnr_pu0_ctl = (1 - v_rnr_ptrt_u0)*q /
        ( (1 - v_rnr_ptrt_u0)*(1 - q) + (1 - v_rnr_ptrt_u1)*q ),
      v_rnr_pu0_trt = v_rnr_ptrt_u0*q /
        (v_rnr_ptrt_u0*(1 - q) + v_rnr_ptrt_u1*q ),
      v_rnr_pu0 = ifelse(treat, v_rnr_pu0_trt, v_rnr_pu0_ctl),

      # lhs r_t = 1 means resp = 1
      v_rnr_beta_ctl = solve(prob=v_rnr_pu0_ctl, delta=d0,
                             lhs=1-v_rnr_resp_ctl_ctl),
      v_rnr_beta_trt = solve(prob=v_rnr_pu0_trt, delta=d1,
                             lhs=1-v_rnr_resp_trt_trt),

      # resp(trt=0) | trt=1
      v_rnr_resp_ctl_trt = v_rnr_pu0 * inv_logit(v_rnr_beta_ctl) +(1-v_rnr_pu0)*inv_logit(v_rnr_beta_ctl+d0),
      # resp(trt=1) | trt=0
      v_rnr_resp_trt_ctl = v_rnr_pu0*inv_logit(v_rnr_beta_trt) +
                               (1-v_rnr_pu0)*(inv_logit(v_rnr_beta_trt+d1)),

      # Sensitivity adjusted measure of resp(treat) | observed treatment
      resp_trt = ifelse(treat, v_rnr_resp_trt_trt, v_rnr_resp_trt_ctl),
      resp_ctl = ifelse(treat, v_rnr_resp_ctl_trt, v_rnr_resp_ctl_ctl),

      # Record params
      q=q, dp=dp, d0=d0, d1=d1
    )

  if (!debug) {
    ret <- ret %>%
      dplyr::select(-starts_with("v_rnr_"))
  }

  return(ret)
}

inv_logit <- binomial()$linkinv
