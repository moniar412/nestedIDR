library(mixtools)

source("funs/funs_general.R")
source("funs/funs_nested_idr.R")
source("funs/funs_methods_warpper.R")
source("package/newIDRcode/estIDRopt.R")
source("package/newIDRcode/getIDR.R")
dyn.load("package/newIDRcode/f_funs/idr.dll")

source("simulation/fun_simu.R")

load("simulation/data/smset_4.Rdata")

simu_setting$simuTimes=200
set.seed(simu_setting$randSeed)
simu_res <- run_simu(simu_setting,nlm_Bound,nlm_control,out_control,par0.idrm,par0.idra,par0.idrs)


save(file =paste("simulation/simu_res/smres_",simu_setting$id,".Rdata",sep=""),
     simu_res,simu_setting,nlm_Bound,nlm_control,out_control,par0.idra,par0.idrm,par0.idrs,tpar)



