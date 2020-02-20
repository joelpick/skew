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

THBW$year_sex <- paste(THBW$year,THBW$sex)

## THBW is the dataset used for the survival analysis
## mod_tarsus_fledge is model of survival from day 15 to fledging for tarsus
## mod_tarsus_recruit is model of survival from fledging to recruitment for tarsus

## function to apply a function to a subset in a dataframe and return a dataframe
groupFunc_dataFrame <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, c(y, make.row.names=FALSE))
}

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

	return(pmvnorm(lower=c(0,0), mean=g_s, sigma=V_s))
}

wD_func<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	g_s<-cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)$cM+beta*z+gamma*z^2
    V_s<-cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)$cV+V_nest+diag(2)

    g_sc<-cmvnorm(mean=c(0,0), sigma=V_s, cond=g_s, keep_var=1, cond_var=2)$cM

	V_sc<-cmvnorm(mean=c(0,0), sigma=V_s, cond=g_s, keep_var=1, cond_var=2)$cV

    ch1<-dnorm(g_sc, 0, sqrt(V_sc))*pnorm(g_s[2], 0, sqrt(V_s[2,2]))*(beta[1]+2*gamma[1]*z+(V_s[1,2]/V_s[2,2])*(beta[2]+2*gamma[2]*z))
    ch2<-pnorm(g_sc, 0, sqrt(V_sc))*dnorm(g_s[2], 0, sqrt(V_s[2,2]))*(beta[2]+2*gamma[2]*z)

    return(ch1+ch2)
}


beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model$Sol))
gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model$Sol))

eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model$Sol))
eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))

eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model$Sol))
eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))

X_eta1<-model$X[which(model$X[,"at.level(age, 1):age15"]==1),eta1_pos]
X_eta2<-model$X[which(model$X[,"at.level(age, 1):age15"]==0),eta2_pos]

z<-THBW[,paste0(trait, "C")]

for(i in 1:nrow(model$Sol)){

	eta1 <- X_eta1%*%model$Sol[i,eta1_pos]
	eta2 <- X_eta2%*%model$Sol[i,eta2_pos]
    beta<-model$Sol[i,beta_pos]
    gamma<-model$Sol[i,gamma_pos]
    
    mu_etaz<-c(mean(eta1), mean(eta2), mean(z)) 
    V_etaz<-cov(cbind(eta1, eta2, z))

	V_nest<-matrix(model$VCV[i,grep("nest", colnames(model$VCV))], 2,2)

    w_func(z[1], mu_etaz, V_etaz, beta, gamma,  V_nest)
    wD_func(z[1], mu_etaz, V_etaz, beta, gamma,  V_nest)

}

## empty vectors
S_skew <- vector(length=n_it)
S_normal <- vector(length=n_it)

##  partial derivative of fitness function



## fixed effects shared between morph and survival analyses
fixedEffects <- formula(~ male_present + hatch_day + clutch_sizeC + nest_hatch_dateC + sex + year)

V<-MCMCglmm::rIW(diag(2), 3)
mu<-rnorm(2)
P<-mvtnorm::pmvnorm(lower=as.vector(mu), mean=as.vector(c(0,0)), sigma=V)


for(j in 1:n_it){

	## for every year-sex combination, sample with replacement from fixed effects and trait, to keep dependencies between them
	DM_permute <- groupFunc_dataFrame(THBW, "year_sex", function(x){
		samp <- sample(x=1:nrow(x), size=N, replace=TRUE) 		# permuted rows
		trait_skew <- x[samp,trait] 							# permute trait
		DM_year_permute <- model.matrix(fixedEffects,x)[samp,]	# permute design matrix
		i1 <- DM_year_permute %*% betaS1_full[j,colnames(DM_year_permute)]	# intercept for fledging model
		i2 <- DM_year_permute %*% betaS2_full[j,colnames(DM_year_permute)]	# intercepts for recruitment model
		u<-mvrnorm(N, c(0,0), matrix(Vn[j,],2,2))
		i1<-i1+u[,1]
		i2<-i2+u[,2]
		trait_normal <- rcmvnorm(df= cbind(trait_skew,i1,i2), keep = 1, cond = cbind(NA,i1,i2))	# simulate 'normally distributed' trait values with same dependencies among covariates
		return(data.frame(trait_skew=trait_skew,trait_normal=trait_normal,i1=i1,i2=i2,year_sex=unique(x$year_sex)))
	})

	## data for skewed trait
	pars_skew <- data.frame(
		x = DM_permute[,"trait_skew"],
		i1 = DM_permute[,"i1"],
		s1 = betaS1_full[j,trait],
		q1 = betaS1_full[j,paste0(trait,2)],
		i2 = DM_permute[,"i2"],
		s2 = betaS2_full[j,trait],
		q2 = betaS2_full[j,paste0(trait,2)])
	
	## data for normal trait - just replace x
	pars_normal <- pars_skew 
	pars_normal$x <- DM_permute[,"trait_normal"]

	## partial derivatives over skewed trait distribution
	wD_skew <- eval(wD_func, pars_skew)

	## expected fitness over skewed trait distribution
	mW_skew <- eval(w_func, pars_skew)

	## selection gradient for skewed trait
	S_skew[j] <- mean(tapply(wD_skew,DM_permute$year_sex,mean) / tapply(mW_skew,DM_permute$year_sex,mean))

	## partial derivatives over normal trait distribution
	wD_normal <- eval(wD_func, pars_normal)
	
	## expected fitness over normal trait distribution
	mW_normal <- eval(w_func, pars_normal)
	
	## selection gradient for normal trait
	S_normal[j] <- mean(tapply(wD_normal,DM_permute$year_sex,mean) / tapply(mW_normal,DM_permute$year_sex,mean))

	## progress bar
	if((j/n_it) %in% ((1:10)/10)) cat(paste0((j/n_it) * 100,"% "))
}

save(S_normal,S_skew,file=paste0(wd,"Data/Intermediate/selection_gradient_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

par(mfrow=c(3,1))
hist(S_normal, xlim=range(c(S_normal,S_skew)), breaks=50)
hist(S_skew, xlim=range(c(S_normal,S_skew)), breaks=50)
hist((S_skew-S_normal), breaks=50)

sum(S_skew>S_normal)/n_it
sum(S_normal<0)/n_it
sum(S_skew<0)/n_it