rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(parallel)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))
fixed_w <- analysis_options$fixed_w


## make data in long format for survival analysis
fledge_dat<-THBW_noRep
colnames(fledge_dat)[which(colnames(fledge_dat)=="fledge")]<-"survival"
colnames(fledge_dat)[which(colnames(fledge_dat)=="recruit")]<-"age"
fledge_dat$age<-15

recruit_dat<-THBW_noRep[which(THBW_noRep$fledge),]
colnames(recruit_dat)[which(colnames(recruit_dat)=="recruit")]<-"survival"
colnames(recruit_dat)[which(colnames(recruit_dat)=="fledge")]<-"age"
recruit_dat$age<-25

surv_dat<-rbind(fledge_dat, recruit_dat)
surv_dat$age<-as.factor(surv_dat$age)


prior_bv <- list( R = list(V = 1, fix = 1), 
			   G = list(
			   	G1 = list(V = diag(2), nu = 2, alpha.mu=c(0,0), alpha.V=diag(2)*100))
			 )
it_scale <- 100

mclapply(c("tarsus_mm","headbill_mm","wing_mm","weight_g"), function(trait){

	traitC <- paste0(trait,"C")
	traitC2 <- paste0(trait,"C2")

	fixed <- as.formula(paste0("survival ~ -1+ at.level(age,1):(age + ", traitC," + ", traitC2, " + ",fixed_w,")+ at.level(age,2):(age + ", traitC," + ", traitC2, " + ",fixed_w,")"))


	model_w <- MCMCglmm(fixed=fixed, random=~us(age):nest, data=surv_dat, prior=prior_bv, family="threshold", nitt=13000*it_scale, thin=10*it_scale, burnin=3000*it_scale)

	save(model_w, file = paste0(wd,"Data/Intermediate/day15_survival_model_",trait,".Rdata"))

}, mc.cores = 4)

