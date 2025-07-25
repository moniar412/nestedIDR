

library(gumbel)
library(mvtnorm)


# -----------------------------------------------
# Generalized marginal quantile mapping
# -----------------------------------------------
qgmm.marginal2 <- function (u, theta, res = 1000, spread = 5) {
  d <- dim(u)[2]     # R dimensions for this block
  m <- theta$m       # 3 components per block
  n.samples <- round(res * theta$pie)
  n.samples[n.samples == 0] <- 2
  
  # Create grid of evaluation
  s <- NULL
  for (i in 1:d) {
    for (j in 1:m) {
      m.ij <- theta$mu[[j]][i]
      sd.ij <- theta$sigma[[j]][i]
      s <- c(s, seq(m.ij - spread * sd.ij, m.ij + spread * sd.ij, length.out = n.samples[j]))
    }
  }
  dim(s) <- c(sum(n.samples), d)
  
  # Evaluate CDF on grid
  eval <- array(0, dim(s))
  for (i in 1:d) {
    for (j in 1:m) {
      eval[, i] <- eval[, i] + theta$pie[j] * pnorm(s[, i], mean = theta$mu[[j]][i], sd = theta$sigma[[j]][i])
    }
  }
  
  # Invert
  z.out <- NULL
  for (j in 1:d) {
    z.out <- c(z.out, approxfun(eval[, j], s[, j], rule = 2)(u[, j]))
  }
  z.out.is.na <- is.na(z.out)
  if (any(z.out.is.na)) {
    z.out[z.out.is.na & u >= 1] <- Inf
    z.out[z.out.is.na & u <= 0] <- -Inf
  }
  dim(z.out) <- c(nrow(u), d)
  return(z.out)
}

# -----------------------------------------------
# Generalized pseudo data generator
# -----------------------------------------------
getPseduMix <- function(u, para, L, R) {
  Z <- matrix(0, nr = nrow(u), nc = ncol(u))
  theta <- vector("list", L)
  
  for (i in 1:L) {
    theta[[i]]$m <- 3
    theta[[i]]$R <- R
    theta[[i]]$pie <- c((1 - para[1]), para[1] * c(1 - para[i + 1], para[i + 1]))
    theta[[i]]$mu[[1]] <- rep(0,R)
    theta[[i]]$mu[[2]] <- rep(para[i+1+L*2],R)
    theta[[i]]$mu[[3]] <- rep(para[i+1+L],R)
    theta[[i]]$sigma[[1]] <- rep(1,R)
    theta[[i]]$sigma[[2]] <- rep(1,R)
    theta[[i]]$sigma[[3]] <- rep(para[i+1+L*3],R)
    
    idx_start <- (i - 1) * R + 1
    idx_end <- i * R
    Z[, idx_start:idx_end] <- qgmm.marginal2(u[, idx_start:idx_end], theta = theta[[i]])
  }
  
  return(Z)
}

# -----------------------------------------------
# Parameter pack/unpack (unchanged)
# -----------------------------------------------
par_list2vec <- function(par_l) {
  par_v <- c(par_l$pi_g[2], par_l$pi_k[, 2], par_l$mu, par_l$mu_k0, par_l$sgm, par_l$rho)
  return(par_v)
}

par_vec2list <- function(par_v, L) {
  par_l <- list()
  par_l$pi_g <- c(1 - par_v[1], par_v[1])
  par_l$pi_k <- matrix(NA, nrow = L, ncol = 2)
  par_l$mu <- rep(NA, L)
  par_l$mu_k0 <- rep(NA, L)
  par_l$sgm <- rep(NA, L)
  par_l$rho <- rep(NA, L)
  
  for (i in 1:L) {
    par_l$pi_k[i, ] <- c(1 - par_v[1 + i], par_v[1 + i])
    par_l$mu[i] <- par_v[L + 1 + i]
    par_l$mu_k0[i] <- par_v[L + L + 1 + i]
    par_l$sgm[i] <- par_v[L + L + L + 1 + i]
    par_l$rho[i] <- par_v[L + L + L + L + 1 + i]
  }
  return(par_l)
}

# -----------------------------------------------
# Generalized negative log-likelihood
# -----------------------------------------------
getNegLoglik <- function(par_v, Z, L, R) {
  par <- par_vec2list(par_v, L)
  total_dim <- L * R
  
  temp_lik_g1 <- matrix(0, nc = L, nr = nrow(Z))
  for (li in 1:L) {
    idx <- ((li - 1) * R + 1):(li * R)
    temp_lik_g1[, li] <- par$pi_k[li, 1] * dmvnorm(Z[, idx], mean = rep(par$mu_k0[li], R), sigma = diag(R)) +
                         par$pi_k[li, 2] * dmvnorm(Z[, idx], mean = rep(par$mu[li], R),
                                                    sigma = sgmToCovm(rep(par$sgm[li], R), par$rho[li]))
  }
  lik_g1 <- apply(temp_lik_g1, 1, prod)
  lik_g0 <- dmvnorm(Z, mean = rep(0, total_dim), sigma = diag(total_dim))
  
  neg_log_lik <- -sum(log(par$pi_g[1] * lik_g0 + par$pi_g[2] * lik_g1))
  return(neg_log_lik)
}

# -----------------------------------------------
# Covariance matrix builder for general R
# -----------------------------------------------
sgmToCovm <- function(sgm, rho) {
  R <- length(sgm)
  Sigma <- diag(sgm^2)
  if (R > 1) {
    for (i in 1:(R - 1)) {
      for (j in (i + 1):R) {
        Sigma[i, j] <- Sigma[j, i] <- sgm[i] * sgm[j] * rho
      }
    }
  }
  return(Sigma)
}

# -----------------------------------------------
# Rank transform
# -----------------------------------------------
Uhat <- function(x) {
  if (is.vector(x)) {
    x <- matrix(x, length(x), 1)
  }
  apply(x, 2, rank, ties.method = "max") / (nrow(x) + 1)
}

# -----------------------------------------------
# Iterative optimizer (needs L,R as input)
# -----------------------------------------------
opt_iterative <- function(par0_vec, X, L, R,
                          nlm_Bound, nlm_control,
                          out_control) {
  u <- Uhat(X)
  epar <- par0_vec
  epar0 <- rep(0, length(epar))
  
  flag_exit <- F
  itern <- 0
  neg_loglik_trace <- rep(0, out_control$iterMax)
  neg_loglik_trace_nlm <- rep(0, out_control$iterMax)
  
  while (TRUE) {
    itern <- itern + 1
    Z <- getPseduMix(u, epar, L, R)
    neg_loglik_trace[itern] <- getNegLoglik(epar, Z, L, R)
    
    if (out_control$verbose$basic) print(paste(itern, neg_loglik_trace[itern]))
    if (out_control$verbose$par) print(epar)
    
    if (itern >= out_control$iterMax) flag_exit <- T
    if (sum(abs(epar - epar0) > out_control$eps_parVec) == 0) flag_exit <- T
    if ((itern > 1) && (abs(neg_loglik_trace[itern - 1] - neg_loglik_trace[itern]) < out_control$eps_loglik)) flag_exit <- T
    if (flag_exit) break
    
    output <- nlminb(epar, getNegLoglik, Z = Z, L = L, R = R,
                     lower = nlm_Bound$low, upper = nlm_Bound$up, control = nlm_control)
    
    neg_loglik_trace_nlm[itern] <- output$objective
    epar0 <- epar
    epar <- output$par
  }

  idr_lab <- getidrLab(epar,Z,L,R)
  idr_all <- getidrAll(epar,Z,L,R)
  
  return(list(para = epar, n_iter = itern,
              neg_loglik_trace = neg_loglik_trace[neg_loglik_trace != 0],
              neg_loglik_trace_nlm = neg_loglik_trace_nlm[neg_loglik_trace_nlm != 0],
              idr_lab = idr_lab,
              idr_all = idr_all))
}



getNegLoglikCp <- function(par_v, Z, L, R) {
  par <- par_vec2list(par_v, L)
  total_dim <- L * R
  
  temp_lik_g1 <- matrix(0, nc = L, nr = nrow(Z))
  for (li in 1:L) {
    idx <- ((li - 1) * R + 1):(li * R)
    temp_lik_g1[, li] <- par$pi_k[li, 1] * dmvnorm(Z[, idx], mean = rep(par$mu_k0[li], R), sigma = diag(R)) +
                         par$pi_k[li, 2] * dmvnorm(Z[, idx], mean = rep(par$mu[li], R),
                                                    sigma = sgmToCovm(rep(par$sgm[li], R), par$rho[li]))
  }
  lik_g1 <- apply(temp_lik_g1, 1, prod)
  lik_g0 <- dmvnorm(Z, mean = rep(0, total_dim), sigma = diag(total_dim))
  log_lik <- sum(log(par$pi_g[1] * lik_g0 + par$pi_g[2] * lik_g1))
  
  # Marginal likelihoods (independence)
  temp_lik_g1_mg <- matrix(0, nc = L * R, nr = nrow(Z))
  for (li in 1:L) {
    for (rj in 1:R) {
      idx <- (li - 1) * R + rj
      temp_lik_g1_mg[, idx] <- par$pi_k[li, 1] * dnorm(Z[, idx], mean = par$mu_k0[li], sd = 1) +
                               par$pi_k[li, 2] * dnorm(Z[, idx], mean = par$mu[li], sd = par$sgm[li])
    }
  }
  temp_lik_g0_mg <- matrix(0, nc = L * R, nr = nrow(Z))
  for (ri in 1:(L * R)) {
    temp_lik_g0_mg[, ri] <- dnorm(Z[, ri], mean = 0, sd = 1)
  }
  temp_lik0 <- par$pi_g[1] * temp_lik_g0_mg + par$pi_g[2] * temp_lik_g1_mg
  lik_mg <- apply(temp_lik0, 1, prod)
  
  neg_log_lik_cp <- -(log_lik - sum(log(lik_mg)))
  return(neg_log_lik_cp)
}



getNegLoglik_warpper <- function(par0_vec, u, L, R) {
  Z <- getPseduMix(u, par0_vec, L, R)
  negLogLik <- getNegLoglik(par0_vec, Z, L, R)
  return(negLogLik)
}

getNegLoglikCp_warpper <- function(par0_vec, u, L, R) {
  Z <- getPseduMix(u, par0_vec, L, R)
  negLogLikCp <- getNegLoglikCp(par0_vec, Z, L, R)
  return(negLogLikCp)
}



getidrLab <- function(par_v, Z, L, R) {
  par <- par_vec2list(par_v, L)
  idrm <- matrix(0, nr = nrow(Z), nc = L)
  
  for (li in 1:L) {
    idx <- ((li - 1) * R + 1):(li * R)
    temp0 <- par$pi_k[li, 1] * dmvnorm(Z[, idx], mean = rep(par$mu_k0[li], R), sigma = diag(R))
    temp1 <- par$pi_k[li, 2] * dmvnorm(Z[, idx], mean = rep(par$mu[li], R),
                                       sigma = sgmToCovm(rep(par$sgm[li], R), par$rho[li]))
    idrm[, li] <- temp0 / (temp0 + temp1)
  }
  return(idrm)
}



getidrAll <- function(par_v, Z, L, R) {
  par <- par_vec2list(par_v, L)
  total_dim <- L * R
  
  temp_g1 <- matrix(0, nc = L, nr = nrow(Z))
  for (li in 1:L) {
    idx <- ((li - 1) * R + 1):(li * R)
    temp_g1[, li] <- par$pi_k[li, 1] * dmvnorm(Z[, idx], mean = rep(par$mu_k0[li], R), sigma = diag(R)) +
                     par$pi_k[li, 2] * dmvnorm(Z[, idx], mean = rep(par$mu[li], R),
                                                sigma = sgmToCovm(rep(par$sgm[li], R), par$rho[li]))
  }
  lik_g1 <- apply(temp_g1, 1, prod)
  lik_g0 <- dmvnorm(Z, mean = rep(0, total_dim), sigma = diag(total_dim))
  idro <- lik_g0 / (lik_g1 + lik_g0)
  return(idro)
}



idr2IDR <- function(idrv) {
  o <- order(idrv)
  idrv.o <- idrv[o]
  idrv.rank <- rank(idrv.o, ties.method = "max")
  
  top.mean <- function(index, x) mean(x[1:index])
  IDRv.o <- sapply(idrv.rank, top.mean, idrv.o)
  IDRv <- rep(NA, length(IDRv.o))
  IDRv[o] <- IDRv.o
  return(IDRv)
}


dataSimulating <- function(para, sampsize, L, R) {
  total_dim <- L * R
  
  ### Group g = 1 proportion
  ng.temp <- sum(rbinom(sampsize, 1, para$pi_g[2]))
  ng <- c(sampsize - ng.temp, ng.temp)
  
  ### Generate g = 0
  X0 <- rmvnorm(ng[1], mean = rep(0, total_dim), sigma = diag(total_dim))
  
  ### Generate g = 1
  lbk <- matrix(0, nr = ng[2], nc = L)  # k labels for each block
  X1 <- matrix(0, nr = ng[2], nc = total_dim)
  
  for (i in 1:L) {
    # Pick block membership k in {0,1}
    lbk[, i] <- sample(c(0, 1), ng[2], prob = para$pi_k[i, ], replace = TRUE)
    nk <- table(factor(lbk[, i], levels = c(0,1)))  # Always has 2 entries
    
    # Indices for this block
    idx <- ((i - 1) * R + 1):(i * R)
    
    # k = 0: standard
    covm0 <- diag(R)  # identity covariance for k = 0
    if (nk[1] > 0) {
      X1[lbk[, i] == 0, idx] <- rmvnorm(n = nk[1], mean = rep(para$mu_k0[i], R), sigma = covm0)
    }
    
    # k = 1: shifted, with correlation
    covm1 <- sgmToCovm(rep(para$sgm[i], R), para$rho[i])
    if (nk[2] > 0) {
      X1[lbk[, i] == 1, idx] <- rmvnorm(n = nk[2], mean = rep(para$mu[i], R), sigma = covm1)
    }
  }
  
  ### Combine X
  X <- as.data.frame(rbind(X0, X1))
  
  temp <- NULL
  for (i in 1:L) {
    for (j in 1:R) {
      temp <- c(temp, paste("M", i, "_R", j, sep = ""))
    }
  }
  colnames(X) <- temp
  
  ### Labels
  lb <- data.frame(matrix(0, nr = sampsize, nc = 1 + L))
  lb[, 1] <- c(rep(0, ng[1]), rep(1, ng[2]))
  temp <- "g"
  for (i in 1:L) {
    lb[lb[, 1] == 1, i + 1] <- lbk[, i]
    temp <- c(temp, paste("k_M", i, sep = ""))
  }
  colnames(lb) <- temp
  
  return(list(X = X, lb = lb))
}


sgmToCovm <- function(sgm, rho) {
  R <- length(sgm)
  Sigma <- diag(sgm^2)
  if (R > 1) {
    for (i in 1:(R - 1)) {
      for (j in (i + 1):R) {
        Sigma[i, j] <- Sigma[j, i] <- sgm[i] * sgm[j] * rho
      }
    }
  }
  return(Sigma)
}


