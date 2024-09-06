


dataSimulating <- function(para,sampsize){ 
  
  ### generating data g=0
  ng.temp <- sum(rbinom(sampsize,1,para$pi_g[2]))
  ng <- c(sampsize-ng.temp,ng.temp)
  X0 <- rmvnorm(ng[1],rep(0,2+2),sigma=diag(2+2));
  
  ### generating data g=1
  lbk <- data.frame(matrix(0,nr=ng[2], nc=2))
  X1 <- matrix(0, nr=ng[2], nc=2+2)
  for (i in 1:2){ 
    lbk[i] <- sample(c(0,1),ng[2],prob = para$pi_k[i,],replace = T)
    nk <- table(lbk[i])
    covm0 <- sgmToCovm(c(1,1),0)
    X1[which(lbk[,i]==0),c(i*2-1,i*2)] <- rmvnorm(n = nk[1],rep(para$mu_k0[i],2),sigma=covm0)  #k=0
    covm1 <- sgmToCovm(c(para$sgm[i],para$sgm[i]),para$rho[i])
    X1[which(lbk[,i]==1),c(i*2-1,i*2)] <- rmvnorm(n = nk[2],rep(para$mu[i],2),sigma=covm1)  #k=1
  }
  
  ### combine X
  X <- as.data.frame(rbind(X0,X1))
  temp <- NULL
  for (i in 1:2){
    for(j in 1:2){
      temp <- c(temp, paste("M",i,"_R",j,sep=""))
    }
  }
  colnames(X) <- temp
  
  ### combine lable
  lb <- data.frame(matrix(0,nr=sampsize,nc=1+2))
  temp <- "g"
  lb[,1] <- c(rep(0,ng[1]),rep(1,ng[2]))
  for (i in 2:ncol(lb)){
    lb[which(lb[,1]==1),i] <- lbk[,i-1]
    temp <- c(temp,paste("k_M",i-1,sep=""))
  }
  colnames(lb) <- temp
  ##
  
  return(list(X=X,lb=lb))
}


sgmToCovm <- function(sgm, rho){
  return(matrix(c(sgm[1]^2,sgm[1]*sgm[2]*rho,sgm[1]*sgm[2]*rho,sgm[2]^2),nc=2,nr=2))
}

getROC <- function(elb,tlb,cutoffn=100){
  cutoffs <- seq(0,1,length.out=cutoffn)
  tpn <- rep(0,cutoffn)
  fpn <- rep(0,cutoffn)
  
  for (i in 1:cutoffn){
    p <- elb<=cutoffs[i]
    tpn[i] <- sum(p&tlb)
    fpn[i] <- sum(p&(!tlb))
  }
  tpr <- tpn/sum(tlb)
  fpr <- fpn/sum(!tlb)
  
  return(list(tpn=tpn,fpn=fpn,tpr=tpr,fpr=fpr))
}


getNSP <- function(idrv,NSPcutn){
  
  IDRv <- idr2IDR(idrv)
  IDRr <- rank(IDRv,ties.method = "random")
  temp <- rep(0,length(NSPcutn))
  
  for (i in 1:length(NSPcutn)){
    temp[i] <- IDRv[which(IDRr==NSPcutn[i])]
  }
  
  res <- cbind(NSP=NSPcutn,IDR=temp)
  return(res)
}


getARI <-function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}


get_IDR <- function(idr){
  rk <- rank(idr,ties.method = "random")
  rkdidr <- rep(0, length(idr))
  rkdidr[rk] <- idr
  tmpIDR <- rep(0, length(rkdidr))
  for (i in 1:length(rkdidr)){
    tmpIDR[i] <- mean(rkdidr[1:i])
  }
  IDR <- tmpIDR[rk]
  return(IDR)
}

# Calculates Davies-Bouldin's cluster separation measure
index.DB <- function (x, cl, d = NULL, centrotypes = "centroids", p = 2, q = 2) {
  if (sum(c("centroids", "medoids") == centrotypes) == 0) 
    stop("Wrong centrotypes argument")
  if ("medoids" == centrotypes && is.null(d)) 
    stop("For argument centrotypes = 'medoids' d cannot be null")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  if (is.null(dim(x))) {
    dim(x) <- c(length(x), 1)
  }
  x <- as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  dAm <- d
  centers <- matrix(nrow = k, ncol = ncol(x))
  if (centrotypes == "centroids") {
    for (i in 1:k) {
      for (j in 1:ncol(x)) {
        centers[i, j] <- mean(x[cl == i, j])
      }
    }
  }
  else if (centrotypes == "medoids") {
    for (i in 1:k) {
      clAi <- dAm[cl == i, cl == i]
      if (is.null(clAi)) {
        centers[i, ] <- NULL
      }
      else {
        centers[i, ] <- .medoid(x[cl == i, ], dAm[cl == 
                                                    i, cl == i])
      }
    }
  }
  else {
    stop("wrong centrotypes argument")
  }
  S <- rep(0, k)
  for (i in 1:k) {
    ind <- (cl == i)
    if (sum(ind) > 1) {
      centerI <- centers[i, ]
      centerI <- rep(centerI, sum(ind))
      centerI <- matrix(centerI, nrow = sum(ind), ncol = ncol(x), 
                        byrow = TRUE)
      S[i] <- mean(sqrt(apply((x[ind, ] - centerI)^2, 1, 
                              sum))^q)^(1/q)
    }
    else S[i] <- 0
  }
  M <- as.matrix(dist(centers, p = p))
  R <- array(Inf, c(k, k))
  r = rep(0, k)
  for (i in 1:k) {
    for (j in 1:k) {
      R[i, j] = (S[i] + S[j])/M[i, j]
    }
    r[i] = max(R[i, ][is.finite(R[i, ])])
  }
  DB = mean(r[is.finite(r)])
  resul <- list(DB = DB, r = r, R = R, d = M, S = S, centers = centers)
  resul
}

