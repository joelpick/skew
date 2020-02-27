rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(sn)
library(parallel)
library(rstan)
rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
parallel::detectCores(logical = FALSE)
if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

load(paste0(wd,"Data/Intermediate/stan_summary_data.Rdata"))

# rskt <- function(n,mu,sigma,alpha,nu){
# 	delta <- alpha/sqrt(1+alpha^2)
# 	b_nu <- (sqrt(nu)*gamma(0.5*(nu-1)))/(sqrt(pi)*gamma(0.5*nu))
# 	sigma_z <- sqrt(nu/(nu-2) - (b_nu*delta)^2)
# 	omega <- sigma/sigma_z
# 	xi <- mu - omega*b_nu*delta 
# 	rst(n,xi,omega,alpha,nu)
# }

se <- function(x){sd(x)/sqrt(length(x))}

hetVarSim <- function(y,x,stanModel_skew_t,stanModel_skew,stanModel_t){

	stan_dat <- list(N=length(y),y=y, K=3, X=model.matrix(~x-1))

	mod_stan_skew_t <- sampling(stanModel_skew_t, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","gamma_E","nu_E"), refresh=0)
	mod_stan_skew <- sampling(stanModel_skew, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","gamma_E"), refresh=0)
	mod_stan_t <- sampling(stanModel_t, data = stan_dat, chains=1, iter = 4000, warmup = 2000, pars=c("beta","sigma_E","nu_E"), refresh=0)
	stanST_out <- summary(mod_stan_skew_t)$summary
	stanSN_out <- summary(mod_stan_skew)$summary
	stanT_out <- summary(mod_stan_t)$summary

	lm_out <- summary(lm(y~x-1))

	return(rbind(
			# means and ses of observed data
			obs=c(tapply(y,x,mean),sigma=NA,gamma=stand_skew(y),nu=NA),
			obs_se=c(tapply(y,x,se),NA,NA,NA),
			
			# means and ses of lm estimates
			lm=c(lm_out$coef[,1],lm_out$sigma,NA,NA),
			lm_se=c(lm_out$coef[,2],NA,NA,NA),
			
			# means and ses of estimates from 3 stan models
			stanST=c(stanST_out[1:6,1]),
			stanST_se=c(stanST_out[1:6,3]),
			stanSN=c(stanSN_out[1:5,1],NA),
			stanSN_se=c(stanSN_out[1:5,3],NA),
			stanT=c(stanT_out[1:4,1],NA,stanT_out[5,1]),
			stanT_se=c(stanT_out[1:4,3],NA,stanT_out[5,3])
			# stan_q2.5=c(stan_out[1:6,4]),
			# stan_q25=c(stan_out[1:6,5]),
			# stan_q75=c(stan_out[1:6,7]),
			# stan_q97.5=c(stan_out[1:6,8])
	))
}



nest_means <- aggregate(wing_mm~male_present+nest,THBW_egg,mean)

sds <- with(nest_means,tapply(wing_mm,male_present,sd))
means <- with(nest_means,tapply(wing_mm,male_present,mean))
ns <- table(nest_means$male_present)
x <- as.factor(rep(1:3,ns))

total_sd <- sd(nest_means[,"wing_mm"] - means[nest_means[,"male_present"]+1])

stanModel_skew_t <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_reparam.stan"))
stanModel_t <- stan_model(file = paste0(wd,"/stan/t1.stan"))
stanModel_skew <- stan_model(file = paste0(wd,"/stan/skew_normal.stan"))




## residual level - normal

out_R_N <- mclapply(1:100,function(i){
	y <- rnorm(sum(ns),rep(means,ns), rep(sds,ns))
	out <- hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
	print(i)
	return(out)
}, mc.cores = 7)
save(out_R_N,file= paste0(wd,"Data/Intermediate/MP_sims_N.Rdata"))

## residual level - skew

out_R_S <- mclapply(1:100,function(i){
	#y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],rep(sds,ns)[i],-5,10) )
	y <- sapply(1:sum(ns), function(i) rskn(1,mean=rep(means,ns)[i],sd=rep(sds,ns)[i],skew=0.95) )
	out<-hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
	print(i)
	return(out)
}, mc.cores = 7)
save(out_R_S,file= paste0(wd,"Data/Intermediate/MP_sims_S.Rdata"))

## residual level - skew - homogenous var

out_R_S_hom <- mclapply(1:100,function(i){
	# y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],total_sd,-5,10) )
	y <- sapply(1:sum(ns), function(i) rskn(1,mean=rep(means,ns)[i],sd=total_sd,skew=0.95) )
	out<-hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
	print(i)
	return(out)
}, mc.cores = 7)
save(out_R_S_hom,file= paste0(wd,"Data/Intermediate/MP_sims_S_hom.Rdata"))

## residual level - skew - homogenous MEANS

out_R_S_homM <-  mclapply(1:100,function(i){
	# y <- sapply(1:sum(ns), function(i) rskt(1,mean(nest_means$wing_mm),rep(sds,ns)[i],-5,10) )
	y <- sapply(1:sum(ns), function(i) rskn(1,mean=mean(nest_means$wing_mm),sd=rep(sds,ns)[i],skew=0.95) )
	out <- hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
	print(i)
	return(out)
}, mc.cores = 7)	
save(out_R_S_homM,file= paste0(wd,"Data/Intermediate/MP_sims_S_homM.Rdata"))





### sims to do 
## - skew t - same means, different variances
## - skew normal - 'residual level - skew' - are fixed effects also biased
## - t - different means, different variances - same means, different variances  - different means, same variances

### output - are significant differences in either direction predicted more than expected

## nest level - normal residuals

nest_effects <- rnorm(sum(ns),rep(means,ns), rep(sds,ns))
nest_id <- rep(nest_effects,each=3)
 (rep(nest_effects,each=3))
rskt(1000,0,sqrt(2),-5,10)


stanModel_RE <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_RE_reparam.stan"))
#stanModel_ME <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_ME_reparam.stan"))







load(file= paste0(wd,"Data/Intermediate/MP_sims_N.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S_hom.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S_homM.Rdata"))

out1 <- array(do.call(c,out_R_N), dim = c(10, 6, 100))
out2 <- array(do.call(c,out_R_S), dim = c(10, 6, 100))
out3 <- array(do.call(c,out_R_S_hom), dim = c(10, 6, 100))
out4 <- array(do.call(c,out_R_S_homM), dim = c(10, 6, 100))

library(abind)
outMeans<- abind(apply(out1,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out2,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out3,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out4,c(1,2),mean)[c(1,3,5,7,9),1:3],
	along=3)
outSD<- abind(apply(out1,c(1,2),sd)[c(1,3,5,7,9),1:3],
	apply(out2,c(1,2),sd)[c(1,3,5,7,9),1:3],
	apply(out3,c(1,2),sd)[c(1,3,5,7,9),1:3],
	apply(out4,c(1,2),sd)[c(1,3,5,7,9),1:3],
	along=3)

outMeansSE<- abind(apply(out1,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out2,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out3,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out4,c(1,2),mean)[c(2,4,6,8,10),1:3],
	along=3)
outSDSE<- abind(apply(out1,c(1,2),sd)[c(2,4,6,8,10),1:3],
	apply(out2,c(1,2),sd)[c(2,4,6,8,10),1:3],
	apply(out3,c(1,2),sd)[c(2,4,6,8,10),1:3],
	apply(out4,c(1,2),sd)[c(2,4,6,8,10),1:3],
	along=3)

x<-array(1:60,dim=c(5,3,4)) + array(rep(0:3,each=15),dim=c(5,3,4))
{
par(mfrow=c(2,1),mar=c(3,5,2,1))

plot(x,outMeans,xaxt="n",main="Estimate",xlab="", ylim=c(32,45), pch=21:25, bg=1:5)
arrows(x,outMeans+outSD,x,outMeans-outSD,angle=90,code=3, length=0.01)
abline(v=c(15,30,45)+1:3)
axis(1,x,rep(c("obs","lm","ST","SN","T"),12),las=2)

for(i in 1:3){
	arrows(x[1,i,1:3],means[i],x[5,i,1:3],means[i],code=0) 
	arrows(x[1,i,4],mean(nest_means$wing_mm),x[5,i,4],mean(nest_means$wing_mm),code=0) 
}

plot(x,outMeansSE,xaxt="n",main="SE",xlab="", ylim=c(0,2.5), pch=21:25, bg=1:5)
arrows(x,outMeansSE+outSDSE,x,outMeansSE-outSDSE,angle=90,code=3, length=0.01)
abline(v=c(15,30,45)+1:3)
axis(1,x,rep(c("obs","lm","ST","SN","T"),12),las=2)
}


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

