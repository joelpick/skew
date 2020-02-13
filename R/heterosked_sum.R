rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

wd <- "~/Dropbox/0_blue_tits/skew"
source(paste0(wd,"/R/functions.R"))
load(paste0(wd,"/Data/Intermediate/chick_data.Rdata"))


nest_means <- aggregate(wing_mm~male_present+nest,THBW_egg,mean)

sds <- with(nest_means,tapply(wing_mm,male_present,sd))
means <- with(nest_means,tapply(wing_mm,male_present,mean))
ns <- table(nest_means$male_present)
x <- as.factor(rep(1:3,ns))
X <- model.matrix(~x-1)

stanModel <- stan_model(file = paste0(wd,"/stan/skew_t/skew_t_PSY_reparam.stan"))


out<-replicate(100,{
	y<-rnorm(sum(ns),rep(means,ns), rep(sds,ns))

	# hist(y, breaks=30)

	stan_dat <- list(N=sum(ns),y=y, K=3, X=X)

	mod_stan <- sampling(stanModel, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","alpha_E","nu_E"))
	stan_out <- summary(mod_stan)$summary[,c(1,4,6,8,9)]

	return(rbind(
		c(tapply(y,x,mean),sigma=NA,alpha=NA,nu=NA,skew=stand_skew(y)),
		c(coef(lm(y~x-1)),summary(lm(y~x-1))$sigma,NA,NA,NA),
		c(stan_out[1:6,1],NA)
	))
})

#save(out,out2,file= paste0(wd,"/Data/Intermediate/MP_sims.Rdata"))
library(abind)
load(file= paste0(wd,"/Data/Intermediate/MP_sims.Rdata"))
out<-abind(out,out2)
dim(out)


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

resid_skew<-get_gamma(out[3,"alpha",],out[3,"nu",])

par(mfrow=c(3,1))
plot(out[2,1,] - out[3,1,],out[3,"nu",])
plot(out[2,1,] - out[3,1,],out[3,"alpha",])
plot(out[2,1,] - out[3,1,],resid_skew)

