rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)
library(mvtnorm)

# wd <- "~/Work"
if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

trait<-"weight_g"

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))

if(trait=="weight_g"){
 model<-mod_weight_bv
}
if(trait=="tarsus_mm"){
 model<-mod_tarsus_bv 
}
rm("mod_weight_bv")
rm("mod_tarsus_bv")

THBW$year_sex <- as.factor(paste(THBW$year,THBW$sex))
THBW$year_sex <- as.factor(rep(1, nrow(THBW)))

## THBW is the dataset used for the survival analysis
## mod_tarsus_fledge is model of survival from day 15 to fledging for tarsus
## mod_tarsus_recruit is model of survival from fledging to recruitment for tarsus

cmvnorm<-function(mean=NULL, sigma=NULL, cond=NULL, df=NULL, keep_var, cond_var=NULL){

	if(!is.null(df)){
		mean <- colMeans(df)
		sigma <- cov(df)
		cond <- df
	}else{
		if(is.null(mean)){stop("mean needs to be specified if data.frame not passed in 'df'")}
		if(is.null(sigma)){stop("sigma needs to be specified if data.frame not passed in 'df'")}
		if(is.null(cond)){stop("cond needs to be specified if data.frame not passed in 'df'")}
	}	
	cV <- sigma[keep_var,keep_var]-sigma[keep_var, -keep_var]%*%solve(sigma[-keep_var, -keep_var])%*%sigma[-keep_var, keep_var]
	if(length(dim(cond))==2){
 		cM <- sapply(1:nrow(cond), function(x) mean[keep_var]+sigma[keep_var, -keep_var]%*%solve(sigma[-keep_var, -keep_var])%*%(cond[x,]-mean)[-keep_var])
 	}else{
 		cM <- mean[keep_var]+sigma[keep_var, -keep_var]%*%solve(sigma[-keep_var, -keep_var])%*%(cond-mean)[-keep_var]
 	}	
 	return(list(cM=cM, cV=cV))
}

w_func<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	g_s<-cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)$cM+beta*z+gamma*z^2
    V_s<-cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)$cV+V_nest+diag(2)

	return(pmvnorm(lower=c(0,0), mean=c(g_s), sigma=V_s))
}

w_func_norm<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	w_func(z, mu_etaz, V_etaz, beta, gamma,  V_nest)*dnorm(z, mu_etaz[3], sqrt(V_etaz[3,3]))
}

w_func_norm<-Vectorize(w_func_norm, "z")

wD_func<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	g_s<-cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)$cM+beta*z+gamma*z^2
    V_s<-cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)$cV+V_nest+diag(2)

    g_sc<-g_s[1]-cmvnorm(mean=c(0,0), sigma=V_s, cond=c(g_s), keep_var=1, cond_var=2)$cM

	V_sc<-cmvnorm(mean=c(0,0), sigma=V_s, cond=c(g_s), keep_var=1, cond_var=2)$cV

    ch1<-dnorm(g_sc[1], 0, sqrt(V_sc))*pnorm(g_s[2], 0, sqrt(V_s[2,2]))*(V_etaz[1,3]/V_etaz[3,3]+beta[1]+2*gamma[1]*z-(V_s[1,2]/V_s[2,2])*(V_etaz[2,3]/V_etaz[3,3]+beta[2]+2*gamma[2]*z))
    ch2<-pnorm(g_sc[1], 0, sqrt(V_sc))*dnorm(g_s[2], 0, sqrt(V_s[2,2]))*(V_etaz[2,3]/V_etaz[3,3]+beta[2]+2*gamma[2]*z)

    return(ch1+ch2)
}

wD_func_norm<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	wD_func(z, mu_etaz, V_etaz, beta, gamma,  V_nest)*dnorm(z, mu_etaz[3], sqrt(V_etaz[3,3]))
}
wD_func_norm<-Vectorize(wD_func_norm, "z")


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

n_it<-nrow(model$Sol)
n_it<-1

n_comb<-nlevels(THBW$year_sex)

beta_skew <- matrix(NA, n_it, n_comb)
beta_normal <- matrix(NA, n_it, n_comb)
s_skew <- matrix(NA, n_it, n_comb)

for(i in 1:n_it){

	eta1 <- X_eta1%*%model$Sol[i,eta1_pos]
	eta2 <- X_eta2%*%model$Sol[i,eta2_pos]

    beta<-model$Sol[i,beta_pos]
    gamma<-model$Sol[i,gamma_pos]
    
    V_nest<-matrix(model$VCV[i,grep("nest", colnames(model$VCV))], 2,2)

    for(j in 1:n_comb){

    	eta1_sub<-eta1[THBW$year_sex==levels(THBW$year_sex)[j]]
    	eta2_sub<-eta2[THBW$year_sex==levels(THBW$year_sex)[j]]
    	z_sub<-z[THBW$year_sex==levels(THBW$year_sex)[j]]

	    mu_etaz<-c(mean(eta1_sub), mean(eta2_sub), mean(z_sub)) 
	    V_etaz<-cov(cbind(eta1_sub, eta2_sub, z_sub))

	    W<-sapply(z_sub, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
	    WD<-sapply(z_sub, wD_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)

	    Wnorm<-integrate(w_func_norm, lower=-Inf, upper=Inf, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)$value
	    WDnorm<-integrate(wD_func_norm, lower=-Inf, upper=Inf, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)$value

	    beta_skew[i,j]<-mean(WD)/mean(W)
	    beta_normal[i,j]<-WDnorm/Wnorm
	    s_skew[i,j]<-mean(W*z_sub/mean(W))-mean(z_sub)
	}
	print(i)
}

beta_skew<-rowMeans(beta_skew)
beta_normal<-rowMeans(beta_normal)
s_skew<-rowMeans(s_skew)

q_pos<-seq(min(z_sub), max(z_sub), length=15)
rec_prob<-tapply(THBW$recruit, cut(z_sub, q_pos), mean)
mid_pos<-q_pos[-length(q_pos)]+diff(q_pos)/2

plot(W[order(z_sub)]~z_sub[order(z_sub)], type="l")
points(rec_prob~mid_pos, col="red")
for(j in seq(1,length(z_sub), length=20)){
arrows(z_sub[order(z_sub)][j], W[order(z_sub)][j], z_sub[order(z_sub)][j]+0.1, W[order(z_sub)][j]+WD[order(z_sub)][j]*0.1, length=0.1)
}
# Check everything looks look OK

save(beta_skew,beta_normal, s_skew,file=paste0(wd,"Data/Intermediate/selection_gradient_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

par(mfrow=c(3,1))
hist(beta_skew, xlim=range(c(beta_normal,beta_skew)), breaks=50)
hist(beta_normal, xlim=range(c(beta_normal,beta_skew)), breaks=50)
hist(beta_normal-beta_skew, breaks=50)

sum(beta_skew<beta_normal)/n_it
sum(beta_skew<0)/n_it
sum(beta_normal<0)/n_it


par(mfrow=c(3,1))
hist(beta_skew, xlim=range(c(beta_normal,beta_skew)), breaks=50)
hist(s_skew, xlim=range(c(beta_normal,beta_skew)), breaks=50)
hist(s_skew-beta_skew, breaks=50)

sum(beta_skew<s_skew)/n_it
sum(beta_skew<0)/n_it
sum(s_skew<0)/n_it

