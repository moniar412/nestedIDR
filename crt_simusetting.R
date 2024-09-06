source("funs_general.R")
source("funs_nested_idr.R")

## setting

simu_setting <- list()
simu_setting$id <- 4
simu_setting$tag <- NULL
simu_setting$randSeed <- 1
simu_setting$simuTimes <- 200
simu_setting$sampleSize <- 5000

simu_setting$NSPcutn <- seq(1,5000,50)

simu_setting$clcut <- 0.5


## true parameters
#M  <- 2                                 #number of labs
#R  <- c(2,2)                            #number of reps in each lab

tpar <- list()
# g:1-signal 0-noise
# k:1-repducible 0-irre

tpar$pi_g <- c(0.5,0.5)                      #pi g=0,1               
tpar$pi_k <- rbind(c(0.2,0.8),c(0.5,0.5))    #pi k=0,1|g=1 for each lab
tpar$mu    <- c(2, 2)                         #mu k=1|g=1 for each lab
tpar$mu_k0 <- c(1, 1)                         #mu k=0|g=1 for each lab
tpar$sgm <- c(1, 1.4)                        #sigma2 k=1|g=1 for each lab
tpar$rho  <- c(0.9, 0.7)                    #rho k=1|g=1 for each lab       # Monia: replace the parameter rho with theta. Also check whether other paramether should be changed?

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

#
###


## par0 for different method

tpar

par0.idrm <- par_list2vec(tpar)                #Monia: Here we specify the initial parameter for nestedIDR, check "par_list2vec" function in "funs_nested_idr.R" and create a vector using your initial parameter.

par0.idra <- list()
par0.idra$muRange <- c(sum(tpar$mu)/2+0.2,sum(tpar$mu)/2-0.2)
par0.idra$sgmRange <- c(sum(tpar$sgm)/2+0.1,sum(tpar$sgm)/2-0.1)
par0.idra$rhoRange <- c(sum(tpar$rho)/2-0.05,sum(tpar$rho)/2+0.05)    # Monia:change the rho range for IDR method, use a rho corresponding to your theta
par0.idra$p0Range <-  c(min(tpar$pi_k[,2]*tpar$pi_g[2]-0.1),max(tpar$pi_k[,2]*tpar$pi_g[2]+0.1))

par0.idrs <- list()
par0.idrs$mu <- c(3,3)
par0.idrs$sgm <- tpar$sgm
par0.idrs$rho <- tpar$rho      # Monia:use a rho corresponding to your theta
par0.idrs$p <- tpar$pi_k[,2]*tpar$pi_g[2]


##

save(file = paste("simulation/data/smset_",simu_setting$id,".Rdata",sep=""),
     simu_setting,tpar,nlm_Bound,nlm_control,out_control,par0.idra,par0.idrm,par0.idrs)




























