rm(list=ls())

options(stringsAsFactors=FALSE)

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

run <- FALSE


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




	nest_means <- aggregate(wing_mm~male_present+nest,THBW,mean)

	sds <- with(nest_means,tapply(wing_mm,male_present,sd))
	means <- with(nest_means,tapply(wing_mm,male_present,mean))
	ns <- table(nest_means$male_present)
	x <- as.factor(rep(1:3,ns))

	total_sd <- sd(nest_means[,"wing_mm"] - means[nest_means[,"male_present"]+1])

if(run){

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

	## residual level - negative skew

	out_R_S <- mclapply(1:100,function(i){
		#y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],rep(sds,ns)[i],-5,10) )
		y <- sapply(1:sum(ns), function(i) rskn(1,mean=rep(means,ns)[i],sd=rep(sds,ns)[i],skew=-0.95) )
		out<-hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
		print(i)
		return(out)
	}, mc.cores = 7)
	save(out_R_S,file= paste0(wd,"Data/Intermediate/MP_sims_S.Rdata"))


	## residual level - positive skew

	out_R_S_P <- mclapply(1:100,function(i){
		#y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],rep(sds,ns)[i],-5,10) )
		y <- sapply(1:sum(ns), function(i) rskn(1,mean=rep(means,ns)[i],sd=rep(sds,ns)[i],skew=0.95) )
		out<-hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
		print(i)
		return(out)
	}, mc.cores = 7)
	save(out_R_S_P,file= paste0(wd,"Data/Intermediate/MP_sims_S_P.Rdata"))


	## residual level - skew - homogenous var

	out_R_S_hom <- mclapply(1:100,function(i){
		# y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],total_sd,-5,10) )
		y <- sapply(1:sum(ns), function(i) rskn(1,mean=rep(means,ns)[i],sd=total_sd,skew=-0.95) )
		out<-hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
		print(i)
		return(out)
	}, mc.cores = 7)
	save(out_R_S_hom,file= paste0(wd,"Data/Intermediate/MP_sims_S_hom.Rdata"))

	## residual level - skew - homogenous MEANS

	out_R_S_homM <-  mclapply(1:100,function(i){
		# y <- sapply(1:sum(ns), function(i) rskt(1,mean(nest_means$wing_mm),rep(sds,ns)[i],-5,10) )
		y <- sapply(1:sum(ns), function(i) rskn(1,mean=mean(nest_means$wing_mm),sd=rep(sds,ns)[i],skew=-0.95) )
		out <- hetVarSim(y,x,stanModel_skew_t,stanModel_skew,stanModel_t)
		print(i)
		return(out)
	}, mc.cores = 7)	
	save(out_R_S_homM,file= paste0(wd,"Data/Intermediate/MP_sims_S_homM.Rdata"))


	## residual level - skew same N
	ns2 <- rep(round(sum(ns)/3),3) + c(0,0,1)
	x2 <- as.factor(rep(1:3,ns2))

	out_R_S_homN <- mclapply(1:100,function(i){
		#y <- sapply(1:sum(ns), function(i) rskt(1,rep(means,ns)[i],rep(sds,ns)[i],-5,10) )
		y <- sapply(1:sum(ns), function(i) 
			rskn(1,mean=rep(means,ns2)[i],sd=rep(sds,ns2)[i],skew=-0.95) )
		out<-hetVarSim(y,x2,stanModel_skew_t,stanModel_skew,stanModel_t)
		print(i)
		return(out)
	}, mc.cores = 7)
	save(out_R_S_homN,file= paste0(wd,"Data/Intermediate/MP_sims_S_homN.Rdata"))



	### sims to do 
	## - skew t - same means, different variances
	## - skew normal - 'residual level - skew' - are fixed effects also biased
	## - t - different means, different variances - same means, different variances  - different means, same variances

	### output - are significant differences in either direction predicted more than expected

	## nest level - normal residuals

	# nest_effects <- rnorm(sum(ns),rep(means,ns), rep(sds,ns))
	# nest_id <- rep(nest_effects,each=3)
	#  (rep(nest_effects,each=3))
	# rskt(1000,0,sqrt(2),-5,10)


	# stanModel_RE <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_RE_reparam.stan"))
	#stanModel_ME <- stan_model(file = paste0(wd,"/stan/skew_t_PSY_ME_reparam.stan"))
}else{
load(file= paste0(wd,"Data/Intermediate/MP_sims_N.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S_hom.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S_homM.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S_homN.Rdata"))
load(file= paste0(wd,"Data/Intermediate/MP_sims_S_P.Rdata"))
}

out1 <- array(do.call(c,out_R_S), dim = c(10, 6, 100))
out2 <- array(do.call(c,out_R_N), dim = c(10, 6, 100))
out3 <- array(do.call(c,out_R_S_hom), dim = c(10, 6, 100))
out4 <- array(do.call(c,out_R_S_homM), dim = c(10, 6, 100))
out5 <- array(do.call(c,out_R_S_homN), dim = c(10, 6, 100))
out6 <- array(do.call(c,out_R_S_P), dim = c(10, 6, 100))

library(abind)

orderMeans <-c(1,2,5,4,3)

outMeans<- abind(apply(out1,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out2,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out3,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out4,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out5,c(1,2),mean)[c(1,3,5,7,9),1:3],
	apply(out6,c(1,2),mean)[c(1,3,5,7,9),1:3],
	along=3)[orderMeans,,]
outSE<- abind(apply(out1,c(1,2),se)[c(1,3,5,7,9),1:3],
	apply(out2,c(1,2),se)[c(1,3,5,7,9),1:3],
	apply(out3,c(1,2),se)[c(1,3,5,7,9),1:3],
	apply(out4,c(1,2),se)[c(1,3,5,7,9),1:3],
	apply(out5,c(1,2),se)[c(1,3,5,7,9),1:3],
	apply(out6,c(1,2),se)[c(1,3,5,7,9),1:3],
	along=3)[orderMeans,,]


# outSkew <- abind(
# 			apply(out1,c(1,2),mean)[c(1,3,5,7,9),4:5],
# 			apply(out2,c(1,2),mean)[c(1,3,5,7,9),4:5],
# 			apply(out3,c(1,2),mean)[c(1,3,5,7,9),4:5],
# 			apply(out4,c(1,2),mean)[c(1,3,5,7,9),4:5],
# 			apply(out5,c(1,2),mean)[c(1,3,5,7,9),4:5],
# 			apply(out6,c(1,2),mean)[c(1,3,5,7,9),4:5],
# 		along=3)

outMeansSE<- abind(apply(out1,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out2,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out3,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out4,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out5,c(1,2),mean)[c(2,4,6,8,10),1:3],
	apply(out6,c(1,2),mean)[c(2,4,6,8,10),1:3],
	along=3)[orderMeans,,]
outSESE<- abind(apply(out1,c(1,2),se)[c(2,4,6,8,10),1:3],
	apply(out2,c(1,2),se)[c(2,4,6,8,10),1:3],
	apply(out3,c(1,2),se)[c(2,4,6,8,10),1:3],
	apply(out4,c(1,2),se)[c(2,4,6,8,10),1:3],
	apply(out5,c(1,2),se)[c(2,4,6,8,10),1:3],
	apply(out6,c(1,2),se)[c(2,4,6,8,10),1:3],
	along=3)[orderMeans,,]

x<-array(1:90,dim=c(5,3,6)) + array(rep(0:5,each=15),dim=c(5,3,6))




bias_fun <- function(output, FUN){
	t(apply(abind(
		output[3,1:3,] - output[1,1:3,],
		output[5,1:3,] - output[1,1:3,],
		output[7,1:3,] - output[1,1:3,],
		output[9,1:3,] - output[1,1:3,]
	,along=3),c(1,3),FUN))
}

bias_fun2 <- function(output, FUN){
	output2 <- abind(
		output[c(1,3,5,7,9),2,] - output[c(1,3,5,7,9),1,],
		output[c(1,3,5,7,9),3,] - output[c(1,3,5,7,9),1,]
	,along=3)
	t(apply(abind(output2[2,,] - output2[1,,],
		output2[3,,] - output2[1,,],
		output2[4,,] - output2[1,,],
		output2[5,,] - output2[1,,]
	,along=3),c(2,3),FUN))
}


orderBias <-c(1,4,3,2)
outBias <- abs(abind(bias_fun2(out1,mean),bias_fun2(out2,mean),bias_fun2(out3,mean),bias_fun2(out4,mean),bias_fun2(out5,mean),bias_fun2(out6,mean),along=3))[orderBias,,]
outBiasSE <- abind(bias_fun2(out1,se),bias_fun2(out2,se),bias_fun2(out3,se),bias_fun2(out4,se),bias_fun2(out5,se),bias_fun2(out6,se),along=3)[orderBias,,]
xBias<-array(1:72,dim=c(4,2,6)) + array(rep(0:5,each=8),dim=c(4,2,6))

# outBias <- abs(abind(bias_fun(out1,mean),bias_fun(out2,mean),bias_fun(out3,mean),bias_fun(out4,mean),bias_fun(out5,mean),bias_fun(out6,mean),along=3))[orderBias,,]
# outBiasSD <- abind(bias_fun(out1,sd),bias_fun(out2,sd),bias_fun(out3,sd),bias_fun(out4,sd),bias_fun(out5,sd),bias_fun(out6,sd),along=3)[orderBias,,]

# xBias<-array(1:72,dim=c(4,3,6)) + array(rep(0:5,each=12),dim=c(4,3,6))




setEPS()
pdf(paste0(wd,"R/plots/figure_SM_fixed_sim.pdf"), height=8, width=15)

{
layout(matrix(c(1:4),ncol=1), height=c(1,5,5,5))

par(mar=c(0,5,0,0))
plot(NULL,ylim=c(-1,1),xlim=range(xBias), xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
text(c(4.5+c(0:5)*9),0,c(1:6),cex=2)

par(mar=c(5,5,1,1),cex.lab=1.5)
plot(x,outMeans,xaxt="n",ylab="Estimate",xlab="", ylim=c(32,50), pch=21:25, cex=1.5, bg=1:5)
arrows(x,outMeans+outSE,x,outMeans-outSE,angle=90,code=3, length=0.01)
abline(v=15*(1:5)+(1:5))
abline(v=c(5.5+c(0:5)*16,10.5+c(0:5)*16),col="grey")
axis(1,x,rep(c("obs","lm","ST","SN","T")[orderMeans],18),las=2)

arrows(x[1,,],27.5,x[5,,],27.5,code=0, xpd=TRUE) 
text(x[3,,],26.4,c(0,1,2), xpd=TRUE)
text(-1,48,"A)",cex=2)

for(i in 1:3){
	arrows(x[1,i,c(1:3,5,6)],means[i],x[5,i,c(1:3,5,6)],means[i],code=0) 
	arrows(x[1,i,4],mean(nest_means$wing_mm),x[5,i,4],mean(nest_means$wing_mm),code=0) 
}



plot(x,outMeansSE,xaxt="n",ylab="SE",xlab="", ylim=c(0,2.5), pch=21:25, cex=1.5, bg=1:5)
arrows(x,outMeansSE+outSESE,x,outMeansSE-outSESE,angle=90,code=3, length=0.01)
abline(v=15*(1:5)+(1:5))
abline(v=c(5.5+c(0:5)*16,10.5+c(0:5)*16),col="grey")
axis(1,x,rep(c("obs","lm","ST","SN","T")[orderMeans],18),las=2)

arrows(x[1,,],-0.6,x[5,,],-0.6,code=0, xpd=TRUE) 
text(x[3,,],-0.75,c(0,1,2), xpd=TRUE)
text(-1,2.25,"B)",cex=2)



plot(xBias,outBias,xaxt="n",xlab="", ylim=c(-1,6), pch=22:25, bg=2:5, cex=1.5, ylab="% Absolute Bias in Contrasts")
abline(h=0)
arrows(xBias,outBias+outBiasSE,xBias,outBias-outBiasSE,angle=90,code=3, length=0.01)
abline(v=8*(1:5)+(1:5))
abline(v=c(4.5+c(0:5)*9),col="grey")
axis(1,xBias,rep(c("lm","ST","SN","T")[orderBias],12),las=2)
#abline(v=c(4.5+c(0:5)*13,8.5+c(0:5)*13),col="grey")

arrows(xBias[1,,],-2.7,xBias[4,,],-2.7,code=0, xpd=TRUE) 
text(c(2.9+c(0:5)*9,6.5+c(0:5)*9),-3.1,rep(c("0-1","0-2"),each=6), xpd=TRUE)
text(-0,5,"C)",cex=2)

}
dev.off()




#### bias












x2<-array(1:90,dim=c(5,2,6)) + array(rep(0:5,each=10),dim=c(5,2,6))

plot(x2,outSkew,xaxt="n",main="Estimate",xlab="", pch=21:25, bg=1:5)


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

