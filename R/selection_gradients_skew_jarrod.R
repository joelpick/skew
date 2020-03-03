rm(list=ls())

# source("~/Work/Skew/R/selection_gradients_skew_jarrod.R")

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

trait<-"weight_g"
save<-TRUE
posterior_mean<-FALSE
by<-c("year", "sex")
pred_cond<-c("year2012", "year2013", "year2014", "year2015","year2016","year2017","year2018","sexM", "timeC") # terms to omit from the variance

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))

if(trait=="weight_g"){
 model_w<-mod_weight_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModWeightPedN20191220_0905.Rdata")))
}
if(trait=="tarsus_mm"){
 model_w<-mod_tarsus_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20191217_1716.Rdata"))) 
}
rm("mod_weight_bv")
rm("mod_tarsus_bv")

pred_cond_pos<-match(pred_cond, colnames(stan_data_weight$X))
pred_pos<-match(setdiff(colnames(stan_data_weight$X), pred_cond), colnames(stan_data_weight$X))

mu_pred_cond<-stan_data_weight$X[,pred_cond_pos]%*%pars(model_z, "beta")[pred_cond_pos,1]
mu_pred<-stan_data_weight$X[,pred_pos]%*%pars(model_z, "beta")[pred_pos,1]

nestD <- pars_ST(model_z, "nest")
residD <- pars_ST(model_z, if(trait=="weight_g"){"E"}else{"ind"})
geneD <- pars(model_z,"_A")[1]
fixedD<-doppelgangR::st.mle(y = mu_pred)$dp
names(fixedD)<-paste0(c("xi", "omega", "alpha", "nu"), "_fixed")

e_st<-list(
  n_st=nestD[,1],
  e_st=residD[,1],
  fixed_st=fixedD
)
# list of environmental distribution parameters: xi, omega, alpha, nu

g_st=c(0, geneD, 0, 1e+16)
# list of genetic distribution parameters: xi, omega, alpha, nu


if(is.null(by)){
	THBW$category <- as.factor(rep(1, nrow(THBW)))
}else{
	THBW$category <- as.factor(apply(THBW[,by], 1, paste, collapse=""))
}

beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model_w$Sol))
gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model_w$Sol))

eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model_w$Sol))
eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))

eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model_w$Sol))
eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))

fixedEffects <- formula(~ sex+male_present + year+hatch_day +clutch_sizeC + nest_hatch_dateC)

X_eta1<-model.matrix(fixedEffects, THBW)
X_eta2<-model.matrix(fixedEffects, THBW)

if(any(!mapply(colnames(X_eta1), colnames(model_w$Sol)[eta1_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta1")}

if(any(!mapply(colnames(X_eta2), colnames(model_w$Sol)[eta2_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta2")}

z<-THBW[,paste0(trait, "C")]

if(posterior_mean){
	n_it<-1
}else{
   	n_it<-nrow(model_w$Sol)
}

n_comb<-nlevels(THBW$category)

beta1 <- matrix(NA, n_it, n_comb)  # linear term in linear best fit
beta2 <- matrix(NA, n_it, n_comb)  # linear term in quadratic best fit
beta3 <- matrix(NA, n_it, n_comb)  # average gradient

int_opt<- matrix(NA, n_it, n_comb)

for(i in 1:n_it){

    if(posterior_mean){
    	eta1 <- X_eta1%*%colMeans(model_w$Sol[,eta1_pos])
		eta2 <- X_eta2%*%colMeans(model_w$Sol[,eta2_pos])

	    beta<-colMeans(model_w$Sol[,beta_pos])
	    gamma<-colMeans(model_w$Sol[,gamma_pos])

	    V_nest<-matrix(colMeans(model_w$VCV[,grep("nest", colnames(model_w$VCV))]), 2,2)
    }else{
		eta1 <- X_eta1%*%model_w$Sol[i,eta1_pos]
		eta2 <- X_eta2%*%model_w$Sol[i,eta2_pos]

	    beta<-model_w$Sol[i,beta_pos]
	    gamma<-model_w$Sol[i,gamma_pos]

	    V_nest<-matrix(model_w$VCV[i,grep("nest", colnames(model_w$VCV))], 2,2)
    }
   

    for(j in 1:n_comb){

    	cat<-levels(THBW$category)[j]

    	eta1_sub<-eta1[THBW$category==cat]
    	eta2_sub<-eta2[THBW$category==cat]
    	z_sub<-z[THBW$category==cat]

	    mu_etaz<-c(mean(eta1_sub), mean(eta2_sub), mean(z_sub)) 
	    V_etaz<-cov(cbind(eta1_sub, eta2_sub, z_sub))

	    W<-sapply(z_sub, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
	    WD<-sapply(z_sub, wD_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)

        WDmin<-wD_func_norm(min(z_sub), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
        WDmax<-wD_func_norm(max(z_sub), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)

        S<-mean(W*z_sub/mean(W))-mean(z_sub)
        C<-mean(W*(z_sub-mean(z_sub)^2)/mean(W))-var(z_sub)

	    beta1[i,j]<-S/var(z_sub)  
	    beta2[i,j]<-betaLA_2(c(list(g_st), e_st), S, C, family="ST")
	    beta3[i,j]<-mean(WD)/mean(W)

	    int_opt[i,j]<-(WDmax<0 & WDmin>0)
	}
	print(i)
}

stop()


if(save){
	save(beta1,beta2, beta3, int_opt, file=paste0(wd,"Data/Intermediate/selection_gradient_",if(is.null(by)){""}else{"by_"}, trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
}

beta_skew<-rowMeans(beta_skew)
beta_normal<-rowMeans(beta_normal)
s_skew<-rowMeans(s_skew)

THBW$traitCat <- bin(z,n=20)+attr(z, "scaled:center")

par(mar=c(5,5,1,1),mfrow=c(3,1), bty="l")
hist(z_sub+attr(z, "scaled:center"), col="grey",main="", xlab=paste(trait), breaks=30)

binPlot(recruit~traitCat, subset(THBW, category==cat), xlab=paste(trait), ylab="Fitness", ylim=c(0,0.1), xlim=range(z)+attr(z, "scaled:center"),text.cex=0.7);
abline(v=mean(z_sub)+attr(z, "scaled:center"), col="grey")
lines(W[order(z_sub)]~I(z_sub[order(z_sub)]+attr(z, "scaled:center")), col="red", lwd=2)

plot(WD[order(z_sub)]~I(z_sub[order(z_sub)]+attr(z, "scaled:center")), type="l", ylab="Fitness Derivative", xlab=paste(trait), lwd=2, col="red")
abline(v=mean(z_sub)+attr(z, "scaled:center"), col="grey")

par(mfrow=c(3,1))
hist(beta_skew, xlim=range(c(beta_normal,beta_skew)), breaks=50)
hist(beta_normal, xlim=range(c(beta_normal,beta_skew)), breaks=50)
hist(beta_normal-beta_skew, breaks=50)

sum(beta_skew<beta_normal)/n_it
sum(beta_skew<0)/n_it
sum(beta_normal<0)/n_it

S<-mean(tapply(THBW$tarsus_mmC[which(THBW$recruit)], THBW$category[which(THBW$recruit)], mean)-tapply(THBW$tarsus_mmC, THBW$category, mean), na.rm=TRUE)


par(mfrow=c(3,1))
hist(beta_skew, xlim=range(c(beta_normal,beta_skew, s_skew)), breaks=50)
hist(s_skew, xlim=range(c(beta_normal,beta_skew, s_skew)), breaks=50)
hist(s_skew-beta_skew, breaks=50)

sum(beta_skew<s_skew)/n_it
sum(beta_skew<0)/n_it
sum(s_skew<0)/n_it

