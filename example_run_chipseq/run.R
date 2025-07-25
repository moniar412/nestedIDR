
load("all_proc.Rdata")
source("funs_methods_warpper.R")


X <- apply(dat,MARGIN = 2,FUN = rank,ties.method="random")
par(mfrow=c(2,2))
plot(dat[,c(1,2)],cex=0.5,ylim=c(0,400),xlim=c(0,400))
plot(dat[,c(3,4)],cex=0.5,ylim=c(0,400),xlim=c(0,400))


plot(X[,c(1,2)],cex=0.5)
plot(X[,c(3,4)],cex=0.5)

cor(dat$l1r1,dat$l1r2)
cor(dat$l2r1,dat$l2r2)



par0 <- list()
# g:1-signal 0-noise
# k:1-repducible 0-irre

par0$pi_g <- c(0.2,0.8)                      #pi g=0,1               
par0$pi_k <- rbind(c(0.2,0.8),c(0.3,0.7))    #pi k=0,1|g=1 for each lab
par0$mu    <- c(3, 3)                         #mu k=1|g=1 for each lab
par0$mu_k0 <- c(1, 1)                         #mu k=0|g=1 for each lab
par0$sgm <- c(1, 1)                        #sigma2 k=1|g=1 for each lab
par0$rho  <- c(0.9, 0.8)                    #rho k=1|g=1 for each lab

## nlminb bounds
nlm_Bound <- list()
################ pi_g1,pi_g1k1m1,pi_g1kim2,mu_k1m1,mu_k1m2,mu_k0m1,mu_k0m2,sgm_k1m1,sgm_k1m2,rho_k1m1,rho_k2m2
nlm_Bound$low <- c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.05, 0.05, -0.99, -0.99)
nlm_Bound$up  <- c(0.95, 0.95, 0.95,  20,  20,  20,  20,   10,   10,  0.99,  0.99)
nlm_control <- list(eval.max=200,iter.max=100,trace=0,rel.tol=0.01)

out_control <- list()
out_control$verbose <- list(basic=F,par=F) # basic, par
out_control$iterMax <- 100
out_control$eps_loglik <- 0.01
out_control$eps_parVec <- c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)

###
## par0 for different method

par0

par0.idrm <- par_list2vec(par0)

par0.idra <- list()
par0.idra$muRange <- c(sum(par0$mu)/2+0.2,sum(par0$mu)/2-0.2)
par0.idra$sgmRange <- c(sum(par0$sgm)/2+0.1,sum(par0$sgm)/2-0.1)
par0.idra$rhoRange <- c(sum(par0$rho)/2-0.05,sum(par0$rho)/2+0.05)
par0.idra$p0Range <-  c(min(par0$pi_k[,2]*par0$pi_g[2]-0.1),max(par0$pi_k[,2]*par0$pi_g[2]+0.1))

par0.idrs <- list()
par0.idrs$mu <- c(3,3)
par0.idrs$sgm <- par0$sgm
par0.idrs$rho <- par0$rho
par0.idrs$p <- par0$pi_k[,2]*par0$pi_g[2]

#########################

res.idrm <- anly_idrm(X,par0.idrm,nlm_Bound,nlm_control,out_control)

#idr for all lab
res.idra <- anly_idra(X,par0.idra)

#idr for single lab
res.idrs <- anly_idrs(X,par0.idrs)

#rankprod
res.rp <- anly_rankprod(X)


save(file = "rd_chipseq/result/res_uta_uw.Rdata",dat,X,res.idra,res.idrm,res.idrs,res.rp,par0,id)











