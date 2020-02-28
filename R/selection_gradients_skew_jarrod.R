rm(list=ls())

# source("~/Work/Skew/R/selection_gradients_skew_jarrod.R")

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)
library(mvtnorm)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

trait<-"weight_g"
save<-FALSE
posterior_mean<-TRUE
by<-c("year", "sex")

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))

if(trait=="weight_g"){
 model<-mod_weight_bv
}
if(trait=="tarsus_mm"){
 model<-mod_tarsus_bv 
}
rm("mod_weight_bv")
rm("mod_tarsus_bv")

if(is.null(by)){
	THBW$category <- as.factor(rep(1, nrow(THBW)))
}else{
	THBW$category <- as.factor(apply(THBW[,by], 1, paste, collapse=""))
}

beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model$Sol))
gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model$Sol))

eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model$Sol))
eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))

eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model$Sol))
eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))

fixedEffects <- formula(~ sex+male_present + year+hatch_day +clutch_sizeC + nest_hatch_dateC)

X_eta1<-model.matrix(fixedEffects, THBW)
X_eta2<-model.matrix(fixedEffects, THBW)

if(any(!mapply(colnames(X_eta1), colnames(model$Sol)[eta1_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta1")}

if(any(!mapply(colnames(X_eta2), colnames(model$Sol)[eta2_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta2")}

z<-THBW[,paste0(trait, "C")]

if(posterior_mean){
	n_it<-1
}else{
   	n_it<-nrow(model$Sol)
}

n_comb<-nlevels(THBW$category)

beta_skew <- matrix(NA, n_it, n_comb)
beta_normal <- matrix(NA, n_it, n_comb)
s_skew <- matrix(NA, n_it, n_comb)
int_opt<- matrix(NA, n_it, n_comb)

for(i in 1:n_it){

    if(posterior_mean){
    	eta1 <- X_eta1%*%colMeans(model$Sol[,eta1_pos])
		eta2 <- X_eta2%*%colMeans(model$Sol[,eta2_pos])

	    beta<-colMeans(model$Sol[,beta_pos])
	    gamma<-colMeans(model$Sol[,gamma_pos])

	    V_nest<-matrix(colMeans(model$VCV[,grep("nest", colnames(model$VCV))]), 2,2)
    }else{
		eta1 <- X_eta1%*%model$Sol[i,eta1_pos]
		eta2 <- X_eta2%*%model$Sol[i,eta2_pos]

	    beta<-model$Sol[i,beta_pos]
	    gamma<-model$Sol[i,gamma_pos]

	    V_nest<-matrix(model$VCV[i,grep("nest", colnames(model$VCV))], 2,2)
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

	    Wnorm<-integrate(w_func_norm, lower=-Inf, upper=Inf, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)$value
	    WDnorm<-integrate(wD_func_norm, lower=-Inf, upper=Inf, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)$value

        WDmin<-wD_func_norm(min(z_sub), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
        WDmax<-wD_func_norm(max(z_sub), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)

	    beta_skew[i,j]<-mean(WD)/mean(W)
	    beta_normal[i,j]<-WDnorm/Wnorm
	    s_skew[i,j]<-(mean(W*z_sub/mean(W))-mean(z_sub))/var(z_sub)
	    int_opt[i,j]<-(WDmax<0 & WDmin>0)
	}
	print(i)
}

stop()


if(save){
	save(beta_skew,beta_normal, s_skew, int_opt, file=paste0(wd,"Data/Intermediate/selection_gradient_",if(is.null(by)){""}else{"by_"}, trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
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

