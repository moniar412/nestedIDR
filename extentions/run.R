
source("nestedIDR_extension.R")

L <- 3
R <- 2
total_dim <- L * R

# Provide your X: n × (L×R)
# Provide your par0_vec with the right length
# Provide bounds matching your vector length

result <- opt_iterative(
  par0_vec = par0_vec, X = X, L = L, R = R,
  nlm_Bound = your_bounds,
  nlm_control = list(iter.max = 100),
  out_control = list(verbose = list(basic = TRUE, par = TRUE),
                     iterMax = 100,
                     eps_loglik = 1,
                     eps_parVec = rep(0.01, length(par0_vec)))
)

