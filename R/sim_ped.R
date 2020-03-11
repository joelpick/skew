rm(list=ls())

# source("~/Work/Skew/R/selection_gradients_skew_jarrod.R")

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)
library(rstan)
library(cubature)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

trait<-"weight_g"
posterior_mean<-TRUE
cond_term<-c("year", "sex") # terms to condition on 
cont_term<-c("timeC")   	# terms to control for

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

t_terms<-c("xi", "omega", "alpha", "nu")

model_zpost<-extract(model_z)

mu_pred<-stan_data_weight$X%*%colMeans(model_zpost$beta)

e_st<-list(n_st=colMeans(as.data.frame(model_zpost)[,paste0(t_terms, "_nest")]),
		 e_st=colMeans(as.data.frame(model_zpost)[,paste0(t_terms, if(trait=="weight_g"){"_E"}else{"_ind"})]), 
		 fixed_st=doppelgangR::st.mle(y = mu_pred)$dp)

g_st<-c(0, mean(model_zpost$sigma_A), 0, 1e+16)

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

fixedEffects <- formula(~ sex+male_present + year+hatch_day +clutch_sizeC + nest_hatch_dateC)
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

nf<-1000
no<-10
ng<-2

ni<-nf*no

ped<-matrix(NA, nf*no*ng, 5)

g<-rz(ni, list(g_st))
n<-rz(ni, list(e_st$n_st))
e<-rz(ni, e_st[-which(names(e_st)=="n_st")])


ped[,1][1:ni]<-1:ni
ped[,5][1:ni]<-1

for(i in 2:ng){

W<-sapply(g+n+e-mu, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
    # fitness evaluated for each (mean-centred) trait value

dam<-sample(1:(ni/2), nf, prob=W[1:(ni/2)], replace=FALSE)
sire<-sample((ni/2)+1:(ni/2), nf, prob=W[(ni/2)+1:(ni/2)], replace=FALSE)


ped[(i-1)*ni+1:ni,1]<-(i-1)*ni+1:ni
ped[(i-1)*ni+1:ni,2]<-rep(dam+ni*(i-2), each=no)
ped[(i-1)*ni+1:ni,3]<-rep(sire+ni*(i-2), each=no)

g<-(g[rep(dam, each=no)]+g[rep(dam, each=no)])/2+rz(ni, list(g_st))/sqrt(2)
n<-rep(rz(nf, list(e_st$n_st)), each=no)
e<-rz(ni, e_st[-which(names(e_st)=="n_st")])



ped[(i-1)*ni+1:ni,4]<-g+n+e
ped[(i-1)*ni+1:ni,5]<-i
}
    # fitness evaluated for each (mean-centred) trait value

ped.ainv <- asreml.Ainverse(ped[,1:3])$ginv

colnames(ped)<-c("id", "dam", "sire", "z", "generation")
ped<-as.data.frame(ped)
ped$id<-as.factor(ped$id)
ped$dam<-as.factor(ped$dam)
ped$sire<-as.factor(ped$sire)
ped$generation<-as.factor(ped$generation)

m1<-asreml(z~1, random=~dam+giv(id), data=subset(ped, generation!=1), ginverse=list(id=ped.ainv))

