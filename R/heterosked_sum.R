rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(sn)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/Dropbox/0_blue_tits/skew"
source(paste0(wd,"/R/functions.R"))
load(paste0(wd,"/Data/Intermediate/chick_data.Rdata"))

load(paste0(wd,"/Data/Intermediate/stan_summary_data.Rdata"))
ls()

rskt <- function(n,mu,sigma,alpha,nu){
	delta <- alpha/sqrt(1+alpha^2)
	b_nu <- (sqrt(nu)*gamma(0.5*(nu-1)))/(sqrt(pi)*gamma(0.5*nu))
	sigma_z <- sqrt(nu/(nu-2) - (b_nu*delta)^2)
	omega <- sigma/sigma_z
	xi <- mu - omega*b_nu*delta 
	rst(n,xi,omega,alpha,nu)
}



## residual level - normal

nest_means <- aggregate(wing_mm~male_present+nest,THBW_egg,mean)

sds <- with(nest_means,tapply(wing_mm,male_present,sd))
means <- with(nest_means,tapply(wing_mm,male_present,mean))
ns <- table(nest_means$male_present)
x <- as.factor(rep(1:3,ns))
X <- model.matrix(~x-1)

total_sd <- sd(nest_means[,"wing_mm"] - means[nest_means[,"male_present"]+1])

stanModel <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_reparam.stan"))


out_R_N <- replicate(100,{
	y <- rnorm(sum(ns),rep(means,ns), rep(sds,ns))

	# hist(y, breaks=30)

	stan_dat <- list(N=sum(ns),y=y, K=3, X=X)

	mod_stan <- sampling(stanModel, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","alpha_E","nu_E"))
	stan_out <- summary(mod_stan)$summary[,c(1,4,6,8,9)]

	# tapply(y - X%*%stan_out[1:3,1],x,mean)
	# tapply(y - X%*%stan_out[1:3,1],x,sd)
	# tapply(y,x,sd)

	return(rbind(
		c(tapply(y,x,mean),sigma=NA,alpha=NA,nu=NA,skew=stand_skew(y)),
		c(coef(lm(y~x-1)),summary(lm(y~x-1))$sigma,NA,NA,NA),
		c(stan_out[1:6,1],NA)
	))
})
#save(out_R_N,file= paste0(wd,"/Data/Intermediate/MP_sims_N.Rdata"))



## residual level - skew

out_R_S <- replicate(100,{
	
	y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],rep(sds,ns)[i],-5,10) )

	# hist(y, breaks=30)

	stan_dat <- list(N=sum(ns),y=y, K=3, X=X)

	mod_stan <- sampling(stanModel, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","alpha_E","nu_E"))
	stan_out <- summary(mod_stan)$summary

ee<-extract(mod_stan)
hist(ee$beta[,3])
hist(extract(mod_stan)[,,"beta[3]"])
	return(rbind(
			obs=c(tapply(y,x,mean),sigma=NA,alpha=NA,nu=NA,skew=stand_skew(y)),
			lm=c(coef(lm(y~x-1)),summary(lm(y~x-1))$sigma,NA,NA,NA),
			stan=c(stan_out[1:6,1],NA),
			q2.5=c(stan_out[1:6,4],NA),
			q25=c(stan_out[1:6,5],NA),
			q75=c(stan_out[1:6,7],NA),
			q97.5=c(stan_out[1:6,8],NA)
	))
})
save(out_R_S,file= paste0(wd,"/Data/Intermediate/MP_sims_S.Rdata"))

## residual level - skew - homogenous var

out_R_S_hom <- replicate(100,{
	
	y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],total_sd,-5,10) )

	# hist(y, breaks=30)

	stan_dat <- list(N=sum(ns),y=y, K=3, X=X)

	mod_stan <- sampling(stanModel, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","alpha_E","nu_E"))
	stan_out <- summary(mod_stan)$summary

	return(rbind(
			obs=c(tapply(y,x,mean),sigma=NA,alpha=NA,nu=NA,skew=stand_skew(y)),
			lm=c(coef(lm(y~x-1)),summary(lm(y~x-1))$sigma,NA,NA,NA),
			stan=c(stan_out[1:6,1],NA),
			q2.5=c(stan_out[1:6,4],NA),
			q25=c(stan_out[1:6,5],NA),
			q75=c(stan_out[1:6,7],NA),
			q97.5=c(stan_out[1:6,8],NA)
	))
})
#out_R_S_hom <- out_R_S
save(out_R_S_hom,file= paste0(wd,"/Data/Intermediate/MP_sims_S_hom.Rdata"))


## nest level - normal residuals

nest_effects <- rnorm(sum(ns),rep(means,ns), rep(sds,ns))
nest_id <- rep(nest_effects,each=3)
 (rep(nest_effects,each=3))
rskt(1000,0,sqrt(2),-5,10)


stanModel_RE <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_RE_reparam.stan"))
#stanModel_ME <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_ME_reparam.stan"))







load(file= paste0(wd,"/Data/Intermediate/MP_sims_N.Rdata"))
load(file= paste0(wd,"/Data/Intermediate/MP_sims_S.Rdata"))
load(file= paste0(wd,"/Data/Intermediate/MP_sims_S_hom.Rdata"))

out <- out_R_S

apply(out,c(1,2),mean)
apply(out,c(1,2),sd)
hist(out[2,1,] - out[3,1,])
hist(out[2,2,] - out[3,2,])
hist(out[2,3,] - out[3,3,])
plot(out[2,1,] - out[3,1,],out[2,2,] - out[3,2,])


get_gamma<-function(alpha, nu){
    delta = alpha / sqrt(1 + alpha^2); 
    b_nu = sqrt(nu/pi) * gamma((nu-1)/2)/gamma(nu/2);
    sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
    gamma = (b_nu*delta)/sigma_z^3 * ( (nu*(3-delta^2))/(nu-3) - (3*nu)/(nu-2)  + 2*(b_nu*delta)^2 );
    return(gamma);
  }
hist(out[3,"alpha",])
hist(out[3,"nu",])

get_gamma(-5,10)
resid_skew<-get_gamma(out[3,"alpha",],out[3,"nu",])

par(mfrow=c(3,1))
plot(out[2,1,] - out[3,1,],out[3,"nu",])
plot(out[2,1,] - out[3,1,],out[3,"alpha",])
plot(out[2,1,] - out[3,1,],resid_skew)

