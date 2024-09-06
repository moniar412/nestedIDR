
run_simu <- function(simu_setting,nlm_Bound,nlm_control,out_control,par0.idrm,par0.idra,par0.idrs){
  rocs <- list()
  nsps <- list()
  aris <- as.data.frame(matrix(0,nc=7,nr=simu_setting$simuTimes))
  colnames(aris) <- c("idrm","idrml1","idrml2","idra","idrl1","idrl2","rankprod")
  epar_idrm <- as.data.frame(matrix(0,nc=11,nr=simu_setting$simuTimes))
  colnames(epar_idrm) <- c("pi_g1","pi_g1k1m1","pi_g1kim2","mu_k1m1","mu_k1m2","mu_k0m1","mu_k0m2","sgm_k1m1","sgm_k1m2","rho_k1m1","rho_k2m2")
  
  epar_idra <- as.data.frame(matrix(0,nc=4,nr=simu_setting$simuTimes))
  colnames(epar_idra) <- c("p0","mu","sigma","rho")
  
  epar_idrs <- as.data.frame(matrix(0,nc=8,nr=simu_setting$simuTimes))
  colnames(epar_idrs) <- c("p1","rho1","mu1","sgm1","p2","rho2","mu2","sgm2")
  
  
  for (simui in 1:simu_setting$simuTimes){
    
    ## data generating
    dat <- dataSimulating(tpar,simu_setting$sampleSize)
    
    #idr nested
    res.idrm <- anly_idrm(dat$X,par0.idrm,nlm_Bound,nlm_control,out_control)
    epar_idrm[simui,] <- res.idrm$epar
    
    #idr for all lab
    res.idra <- anly_idra(dat$X,par0.idra)
    epar_idra[simui,] <- unlist(res.idra$para)
    
    #idr for single lab
    res.idrs <- anly_idrs(dat$X,par0.idrs)
    epar_idrs[simui,] <- res.idrs$epar
    
    #rankprod
    res.rp <- anly_rankprod(dat$X)
    
    #### ROCs
    temp <- getROC(res.idrm$all,dat$lb[,1])
    d1 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="idrm")
    temp <- getROC(res.idrm$lab[,1],dat$lb[,1])
    d11 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="idrml1")
    temp <- getROC(res.idrm$lab[,2],dat$lb[,1])
    d12 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="idrml2")
    temp <- getROC(res.idra$idr,dat$lb[,1])
    d2 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="idra")
    temp <- getROC(res.idrs$idr[,1],dat$lb[,1])
    d3 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="idrl1")
    temp <- getROC(res.idrs$idr[,2],dat$lb[,1])
    d4 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="idrl2")
    temp <- getROC(res.rp,dat$lb[,1])
    d5 <- data.frame(tpr=temp$tpr,fpr=temp$fpr,method="rankprod")
    
    rocs[[simui]] <- rbind(d1,d11,d12,d2,d3,d4,d5)
    
    #### NSPs
    temp <- getNSP(res.idrm$all,simu_setting$NSPcutn)
    d1 <- data.frame(NSP=temp[,1],IDR=temp[,2],method="idrm")
    temp <- getNSP(res.idrm$lab[,1],simu_setting$NSPcutn)
    d11 <- data.frame(NSP=temp[,1],IDR=temp[,2],method="idrml1")
    temp <- getNSP(res.idrm$lab[,2],simu_setting$NSPcutn)
    d12 <- data.frame(NSP=temp[,1],IDR=temp[,2],method="idrml2")
    temp <- getNSP(res.idra$idr,simu_setting$NSPcutn)
    d2 <- data.frame(NSP=temp[,1],IDR=temp[,2],method="idra")
    temp <- getNSP(res.idrs$idr[,1],simu_setting$NSPcutn)
    d3 <- data.frame(NSP=temp[,1],IDR=temp[,2],method="idrl1")
    temp <- getNSP(res.idrs$idr[,2],simu_setting$NSPcutn)
    d4 <- data.frame(NSP=temp[,1],IDR=temp[,2],method="idrl2")
    
    nsps[[simui]] <- rbind(d1,d11,d12,d2,d3,d4)
    
    #### ARIs
    
    aris[simui,]$idrm <- getARI(res.idrm$all<simu_setting$clcut,dat$lb[,1])
    aris[simui,]$idrml1 <- getARI(res.idrm$lab[,1]<simu_setting$clcut,dat$lb[,1])
    aris[simui,]$idrml2 <- getARI(res.idrm$lab[,2]<simu_setting$clcut,dat$lb[,1])
    aris[simui,]$idra <- getARI(res.idra$idr<simu_setting$clcut,dat$lb[,1])
    aris[simui,]$idrl1 <- getARI(res.idrs$idr[,1]<simu_setting$clcut,dat$lb[,1])
    aris[simui,]$idrl2 <- getARI(res.idrs$idr[,2]<simu_setting$clcut,dat$lb[,1])
    aris[simui,]$rankprod <- getARI(res.rp<simu_setting$clcut,dat$lb[,1])
    
    if(0){
      sum((res.idrm$all<simu_setting$clcut)==dat$lb[,1])
      sum((res.idra$idr<simu_setting$clcut)==dat$lb[,1])
      sum((res.idrs$idr[,1]<simu_setting$clcut)==dat$lb[,1])
      sum((res.idrs$idr[,2]<simu_setting$clcut)==dat$lb[,1])
      sum((res.rp<simu_setting$clcut)==dat$lb[,1])
    }
    
  }
  
  return(list(epar_idrm=epar_idrm,epar_idra=epar_idra,epar_idrs=epar_idrs,rocs=rocs,nsps=nsps,aris=aris))
}

