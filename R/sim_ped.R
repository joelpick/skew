rm(list=ls())

# source("~/Work/Skew/R/selection_gradients_skew_jarrod.R")

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(parallel)
library(asreml)
library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)
library(rstan)
library(cubature)

library(doppelgangR)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

trait<-"weight_g"
#posterior_mean<-TRUE
cond_term<-c("year", "sex") # terms to condition on 
cont_term<-c("timeC")   	# terms to control for


############################
## load in data
############################

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))
# survival models (including 2010 when eggs weren't weight) and THBW is the data

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))
# stan_data_weight$X is the design matrix and THBW_egg_noRep the data
# excluding years in which eggs weren't measured and repeat measures

if(trait=="weight_g"){
 model_w<-mod_weight_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModWeightPedN20191220_0905.Rdata")))
}
if(trait=="tarsus_mm"){
 model_w<-mod_tarsus_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20191217_1716.Rdata")))
}
rm("mod_weight_bv", "mod_tarsus_bv")


load(paste0(wd,"Data/Intermediate/stan_summary_data.Rdata"))

############################
## simulation parameters
############################

mu_pred<-stan_data_weight$X %*% pars(model_z,"beta")[,1]

e_st<-list(
	n_st= pars_ST(model_z,"nest")[,1],
	e_st=pars_ST(model_z,if(trait=="weight_g"){"E"}else{"ind"} )[,1], 
	fixed_st=doppelgangR::st.mle(y = mu_pred)$dp
	)

g_st<-c(0, pars(model_z,"sigma_A")[1], 0, 1e+16)

beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model_w$Sol))
# linear coefficients for trait values from the survival model

gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model_w$Sol))
# quadratic coefficients for trait values from the survival model

eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model_w$Sol))
eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))
# positions of non-trait terms in the survival model that pertain to Day15 to fledging. 

eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model_w$Sol))
eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))
# positions of non-trait terms in the survival model that pertain to fledging to recruitment. 

fixedEffects <- formula(~ sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC)
# single trait version of the fixed effect formula used in the survival model

X_eta1<-model.matrix(fixedEffects, THBW)
X_eta2<-model.matrix(fixedEffects, THBW)
# design matrix for non-trait effects in survival model

eta1 <- X_eta1%*%colMeans(model_w$Sol[,eta1_pos])
eta2 <- X_eta2%*%colMeans(model_w$Sol[,eta2_pos])
# non-trait fixed effects predictions from survival model

beta<-colMeans(model_w$Sol[,beta_pos])
gamma<-colMeans(model_w$Sol[,gamma_pos])
# trait coefficients from survival model

V_nest<-matrix(colMeans(model_w$VCV[,grep("nest", colnames(model_w$VCV))]), 2,2)

mu_etaz<-c(mean(eta1), mean(eta2), mean(THBW[,paste0(trait, "C")])) 
V_etaz<-cov(cbind(eta1, eta2, THBW[,paste0(trait, "C")]))
# mean and covariance matrices of linear predictor and traits
mu<-attr(THBW[,paste0(trait, "C")], "scaled:center")



############################
## simulation
############################


nf<-1000
no<-10
ng<-3


n_sims=100

# sims <- mclapply(1:n_sims,function(j,nf=1000,no=10,ng=3){

sims<-matrix(NA,5,n_sims)
for(j in 1:n_sims){


	ni<-nf*no	
	
	ped<-matrix(NA, nf*no*ng, 5)

	## generate breeding values, nest and residual effects
	g<-rz(ni, list(g_st))
	n<-rz(ni, list(e_st$n_st))
	e<-rz(ni, e_st[-which(names(e_st)=="n_st")])

	ped[,1][1:ni]<-1:ni
	ped[,5][1:ni]<-1

	for(i in 2:ng){

		W<-sapply(g+n+e-mu, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
		    # fitness evaluated for each (mean-centred) trait value

		## first half of previous gen pedigree are females/second half males
		## sample who survives - n1/2 of each sex, with their survival probability based in their
		dam<-sample(1:(ni/2), nf, prob=W[1:(ni/2)], replace=FALSE)
		sire<-sample((ni/2)+1:(ni/2), nf, prob=W[(ni/2)+1:(ni/2)], replace=FALSE)

		## put surviving individual into pedigree as parents with no offspring each (random mating)
		ped[(i-1)*ni+1:ni,1]<-(i-1)*ni+1:ni
		ped[(i-1)*ni+1:ni,2]<-rep(dam+ni*(i-2), each=no)
		ped[(i-1)*ni+1:ni,3]<-rep(sire+ni*(i-2), each=no)

		## generate new breeding values, nest and residual effects
		g<-(g[rep(dam, each=no)]+g[rep(sire, each=no)])/2+rz(ni, list(g_st))/sqrt(2)
		n<-rep(rz(nf, list(e_st$n_st)), each=no)
		e<-rz(ni, e_st[-which(names(e_st)=="n_st")])

		ped[(i-1)*ni+1:ni,4]<-g+n+e
		ped[(i-1)*ni+1:ni,5]<-i
	}
	    # fitness evaluated for each (mean-centred) trait value

	colnames(ped)<-c("id", "dam", "sire", "z", "generation")
	ped_df<-as.data.frame(ped)
	ped_df$id<-as.factor(ped_df$id)
	ped_df$dam<-as.factor(ped_df$dam)
	ped_df$sire<-as.factor(ped_df$sire)
	ped_df$generation<-as.factor(ped_df$generation)

	## animal model
	# ped.ainv1 <- asreml.Ainverse(ped[,1:3])
	# assign("ped.ainv", ped.ainv1$ginv, envir = .GlobalEnv) 

	ped.ainv <- asreml.Ainverse(ped[,1:3])$ginv

	m1<-asreml(z~1, random=~dam+giv(id), data=subset(ped_df, generation!=1), ginverse=list(id=ped.ainv),rcov = ~idv(units),trace=FALSE)
	## parent-offspring regression
	ped_df$pz <- (ped[ped[,"dam"],"z"] + ped[ped[,"sire"],"z"])/2
	po_df <- aggregate(cbind(pz,z)~dam,ped_df,mean)
	

	out <- c(summary(m1)$var[c(1,2,4),2],h2_animal=pin(m1, h2_animal~id/(id+dam+R))[1,1], h2_po =coef(lm(z~pz,po_df))[2])
	names(out)[1:3] <- c("Vnest","Va","Ve")
	print(j)
	sims[,j]<-out
	rm("ped.ainv")

	# return(out)
}
# , mc.cores = 7)

rowMeans(sims)


skewModPed_M$var[,1]
skewModPed_M$var["A",1]/sum(skewModPed_M$var[,1])

hist(sims[2,],breaks=20);abline(v=skewModPed_M$var["A",1])
hist(sims[4,],breaks=20);abline(v=skewModPed_M$var["A",1]/sum(skewModPed_M$var[,1]))
hist(sims[5,],breaks=20);abline(v=skewModPed_M$var["A",1]/sum(skewModPed_M$var[,1]))
