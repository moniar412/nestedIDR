

## idr multilab
anly_idrm<- function(X,par0,nlm_Bound,nlm_control,out_control){
  res <- opt_iterative(X = X, par0_vec = par0, nlm_Bound = nlm_Bound, nlm_control = nlm_control, out_control = out_control )
  return(list(all=res$idr_all,lab=res$idr_lab,epar=res$para))
}

## idr for all lab
anly_idra <- function(X,par0,ntry=10){
  res <- getIDR(X,ntry = ntry, 
                 muRange = par0$muRange, sigRange = par0$sgmRange,
                 rhoRange = par0$rhoRange, proRange = par0$p0Range,
                 eps = c(param = 1.e-3, loglik = 1))
  return(res)
}

## idr single lab
library(idr)

anly_idrs <- function(X,par0){
  idr_lab <- matrix(0,nc=2,nr=nrow(X))
  epar <- NULL
  
  for (l in 1:2){
    temp <- est.IDR(x = X[,c(l*2-1,l*2)],mu=par0$mu[l],sigma = par0$sgm[l],rho = par0$rho[l],p = par0$p[l])
    idr_lab[,l] <- temp$idr
    epar <- c(epar,unlist(temp$para))
  }
  names(epar)<-c("p1","rho1","mu1","sgm1","p2","rho2","mu2","sgm2")
  
  return(list(idr=idr_lab,epar=epar))
}

## rankprod
getRankprod <- function(x){
  rank.x <- 1
  ncol.x <- ncol(x)
  for(i in 1:ncol.x){
    rank.x <- rank(x[,i])*rank.x	
  }
  
  rank.prod <- (rank.x)^{1/ncol.x}
  return(rank.prod)
}

anly_rankprod <- function(X0){
  X <- -X0 
  temp <- getRankprod(X)
  res <- temp/max(temp)
  return(res)
}





