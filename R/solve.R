solve_closed <- function(prob, delta, lhs) {
  # Closed-form solution for theta, given
  # 1 - lhs = prob * inv_logit(theta) + (1 - prob) * inv_logit(theta + delta)
  # as presented in Rosenbaum & Rubin
  a <- lhs * exp(delta)
  b <- (lhs - prob) * exp(delta) + lhs - 1 + prob
  c <- lhs - 1
  w <- (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
  theta <- log(w)

  return(theta)
}
solve <- Vectorize(solve_closed)
