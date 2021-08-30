rm(list=ls())

options( stringsAsFactors=FALSE)

library(sn)
library(scales)
library(MASS)
library(mvtnorm)
library(viridis)


## load in data
if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

data_wd <- paste0(wd,"Data/Intermediate/")


reduced <- TRUE

load(paste0(data_wd,"chick_data.Rdata"))
load(paste0(data_wd,"stan_data",if(reduced)"_reduced",".Rdata"))
source(paste0(wd,"R/00_functions.R"))


 {

trait <- "wing_mm"

	X<-stan_data_ped_noRep$X
	

	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_\\d"),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS_delta <- model_z

	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_alpha10_"),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS_alpha10 <- model_z
	
	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_alpha1_"),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS_alpha1 <- model_z

	model_zpostDS_delta <-do.call(rbind,model_zDS_delta)[,-(1:7)]
	betasDS_delta <- model_zpostDS_delta[,grep("beta",colnames(model_zpostDS_delta))]
	mu_allDS_delta <- apply(betasDS_delta, 1, function(x) X%*%x)
	sigmasDS_delta <- cbind(sigma_fix=apply(mu_allDS_delta,2,sd),model_zpostDS_delta[,grep("sigma",colnames(model_zpostDS_delta))])
	gammasDS_delta <- cbind(gamma_fix=apply(mu_allDS_delta,2,stand_skew),model_zpostDS_delta[,grep("gamma",colnames(model_zpostDS_delta))])

	gammasDS_delta <- cbind(gammasDS_delta,gamma_ME=0)
	colnames(sigmasDS_delta)[colnames(sigmasDS_delta)=="sigma_E"] <- "sigma_ME"
	colnames(sigmasDS_delta)[colnames(sigmasDS_delta)=="sigma_ind"] <- "sigma_E"
	colnames(gammasDS_delta)[colnames(gammasDS_delta)=="gamma_ind"] <- "gamma_E"
		

	model_zpostDS_alpha10 <-do.call(rbind,model_zDS_alpha10)[,-(1:7)]
	betasDS_alpha10 <- model_zpostDS_alpha10[,grep("beta",colnames(model_zpostDS_alpha10))]
	mu_allDS_alpha10 <- apply(betasDS_alpha10, 1, function(x) X%*%x)
	sigmasDS_alpha10 <- cbind(sigma_fix=apply(mu_allDS_alpha10,2,sd),model_zpostDS_alpha10[,grep("sigma",colnames(model_zpostDS_alpha10))])
	gammasDS_alpha10 <- cbind(gamma_fix=apply(mu_allDS_alpha10,2,stand_skew),model_zpostDS_alpha10[,grep("gamma",colnames(model_zpostDS_alpha10))])

	gammasDS_alpha10 <- cbind(gammasDS_alpha10,gamma_ME=0)
	colnames(sigmasDS_alpha10)[colnames(sigmasDS_alpha10)=="sigma_E"] <- "sigma_ME"
	colnames(sigmasDS_alpha10)[colnames(sigmasDS_alpha10)=="sigma_ind"] <- "sigma_E"
	colnames(gammasDS_alpha10)[colnames(gammasDS_alpha10)=="gamma_ind"] <- "gamma_E"
	
	
	model_zpostDS_alpha1 <-do.call(rbind,model_zDS_alpha1)[,-(1:7)]
	betasDS_alpha1 <- model_zpostDS_alpha1[,grep("beta",colnames(model_zpostDS_alpha1))]
	mu_allDS_alpha1 <- apply(betasDS_alpha1, 1, function(x) X%*%x)
	sigmasDS_alpha1 <- cbind(sigma_fix=apply(mu_allDS_alpha1,2,sd),model_zpostDS_alpha1[,grep("sigma",colnames(model_zpostDS_alpha1))])
	gammasDS_alpha1 <- cbind(gamma_fix=apply(mu_allDS_alpha1,2,stand_skew),model_zpostDS_alpha1[,grep("gamma",colnames(model_zpostDS_alpha1))])

	gammasDS_alpha1 <- cbind(gammasDS_alpha1,gamma_ME=0)
	colnames(sigmasDS_alpha1)[colnames(sigmasDS_alpha1)=="sigma_E"] <- "sigma_ME"
	colnames(sigmasDS_alpha1)[colnames(sigmasDS_alpha1)=="sigma_ind"] <- "sigma_E"
	colnames(gammasDS_alpha1)[colnames(gammasDS_alpha1)=="gamma_ind"] <- "gamma_E"
	


	col_orderDS <- c("_fix", "_dam_sire","_nest", "_E", "_ME")

	sigmasDS_delta <- sigmasDS_delta[,paste0("sigma",col_orderDS)]^2
	sigmasDS_alpha10 <- sigmasDS_alpha10[,paste0("sigma",col_orderDS)]^2
	sigmasDS_alpha1 <- sigmasDS_alpha1[,paste0("sigma",col_orderDS)]^2

	gammasDS_delta <- gammasDS_delta[,paste0("gamma",col_orderDS)]
	gammasDS_alpha10 <- gammasDS_alpha10[,paste0("gamma",col_orderDS)]
	gammasDS_alpha1 <- gammasDS_alpha1[,paste0("gamma",col_orderDS)]


	stan_betasDS_delta <- t(apply(betasDS_delta,2,mean_CI))
	stan_varDS_delta <- t(apply(sigmasDS_delta,2,mean_CI))
	stan_gammaDS_delta <- t(apply(gammasDS_delta,2,mean_CI))
	
	stan_betasDS_alpha10 <- t(apply(betasDS_alpha10,2,mean_CI))
	stan_varDS_alpha10 <- t(apply(sigmasDS_alpha10,2,mean_CI))
	stan_gammaDS_alpha10 <- t(apply(gammasDS_alpha10,2,mean_CI))
	
	stan_betasDS_alpha1 <- t(apply(betasDS_alpha1,2,mean_CI))
	stan_varDS_alpha1 <- t(apply(sigmasDS_alpha1,2,mean_CI))
	stan_gammaDS_alpha1 <- t(apply(gammasDS_alpha1,2,mean_CI))


}	


setEPS()
pdf(paste0(wd,"R/plots/figure_SM_prior_sensitivity1.pdf"), , height=10, width=10)

	
{
	y_names <- c("Fixed", "Genetic","Nest", "Residual", "ME")
	beta_names <- c( "Intercept", "Year (12)","Year (13)","Year (14)","Year (15)","Year (16)","Year (17)","Year (18)","Time","Sex (M)", "Egg Mass")
layout(matrix(c(1,1,2,3),ncol=2))

cex.axis=1
ylims1=c(0,1.2)
cex.heading = 1.75

	par(mar=c(5,7.5,1,1))
	effectPlot(stan_betasDS_delta[-1,],pch=21,  names=beta_names[-1], xlab="beta", offset=0.15)
	effectPlot(stan_betasDS_alpha10[-1,], add=TRUE, pch=21, bg="grey")
	effectPlot(stan_betasDS_alpha1[-1,], add=TRUE, offset=-0.15)
	legend("topright", c("delta - uniform(-1,1) ","alpha - normal(0,10)","alpha - normal(0,1)"), pch=c(21,19,19), col=c(1,"grey",1))

	effectPlot(stan_varDS_delta[,],pch=21, names=y_names, ylim=c(0.5,5.5), xlab="variance", offset=0.15)
	effectPlot(stan_varDS_alpha10[,], add=TRUE, pch=21, bg="grey")
	effectPlot(stan_varDS_alpha1[,], add=TRUE, offset=-0.15)

	effectPlot(stan_gammaDS_delta[-5,], xlim=c(-4,1), pch=21,names=y_names[1:4], xlab="skew", offset=0.15)
	effectPlot(stan_gammaDS_alpha10[-5,], add=TRUE, pch=21, bg="grey")
	effectPlot(stan_gammaDS_alpha1[-5,], add=TRUE, offset=-0.15)
}

 
 dev.off()




get_gamma<-function(alpha10=NULL, delta=NULL, nu){
	if(is.null(delta)) delta = alpha10 / sqrt(1 + alpha10^2); 
	b_nu = sqrt(nu/pi) * gamma((nu-1)/2)/gamma(nu/2);
	sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
	skew = (b_nu*delta)/sigma_z^3 * ( (nu*(3-delta^2))/(nu-3) - (3*nu)/(nu-2)  + 2*(b_nu*delta)^2 );
	return(skew);
}
delta2alpha <- function(delta) delta/(sqrt(1-delta^2))
alpha2delta <- function(alpha) alpha/sqrt(1+alpha^2)


pp_plot<-function(x,prior,posterior,main="",xlab="",ylim=NULL,xlim=NULL){
	scale_d<-(length(posterior))/(sum(prior))
	prior_y <- prior*scale_d
	hist_out<-hist(posterior,breaks=x, plot=FALSE)
	if(is.null(ylim)) ylim <- c(0,max(prior_y,hist_out$counts))
	if(is.null(xlim)) xlim <- range(x)
	hist_out<-hist(posterior, xlim=xlim,breaks=x, col="grey",main=main,xlab="", ylim=ylim)
	lines(prior_y~x, col="red", lwd=2)
	mtext(xlab, 1, line=3, cex=1.5)
}

	gamma_prior1 <- density(get_gamma(delta=runif(10000,-1,1),nu=runif(10000,4,40)), from=-4,to=4,n=201)
	gamma_prior2 <- density(get_gamma(alpha=rnorm(10000,0,10),nu=runif(10000,4,40)), from=-4,to=4,n=201)
	gamma_prior3 <- density(get_gamma(alpha=rnorm(10000,0,1),nu=runif(10000,4,40)), from=-4,to=4,n=201)


	# alpha_prior1 <- density(delta2alpha(delta=runif(10000,-1,1)), from=-100,to=100,n=201)
	# delta_prior2 <- density(alpha2delta(alpha=rnorm(10000,0,10)), from=-1,to=1,n=201)
	# delta_prior3 <- density(alpha2delta(alpha=rnorm(10000,0,1)), from=-1,to=1,n=201)


	setEPS()
	pdf(paste0(wd,"R/plots/figure_SM_prior_sensitivity2.pdf"), , height=8, width=12)

{
	par(mfrow=c(3,3), mar=c(5,4,1,1))

	# pp_plot(x = alpha_prior1$x, 
	# 	prior = alpha_prior1$y,
	# 	posterior=model_zpostDS_delta[,"alpha_dam_sire"],main="alpha_dam_sire"
	# )
	pp_plot(x = seq(-1,1,length.out=201), 
		prior = dunif(seq(-1,1,length.out=201), -1,1),
		posterior=model_zpostDS_delta[,"delta_dam_sire"],xlab=expression(delta)
	)
	pp_plot(x = seq(4,40,length.out=201), 
		prior = dunif(seq(4,40,length.out=201), 4,40),
		posterior=model_zpostDS_delta[,"nu_dam_sire"],xlab=expression(nu)
	)
	pp_plot(x = gamma_prior1$x, 
		prior = gamma_prior1$y,
		posterior=model_zpostDS_delta[,"gamma_dam_sire"], xlab=expression(gamma)
	)


	pp_plot(x = seq(-40,40,length.out=201), 
		prior = dnorm(seq(-40,40,length.out=201), 0,10),
		posterior=model_zpostDS_alpha10[,"alpha_dam_sire"],xlab=expression(alpha)
	)
	# pp_plot(x = delta_prior2$x, 
	# 	prior = delta_prior2$y,
	# 	posterior=model_zpostDS_alpha10[,"delta_dam_sire"], xlab=expression(delta)
	# )
	pp_plot(x = seq(4,40,length.out=201), 
		prior = dunif(seq(4,40,length.out=201), 4,40),
		posterior=model_zpostDS_alpha10[,"nu_dam_sire"],xlab=expression(nu)
	)
	pp_plot(x = gamma_prior2$x, 
		prior = gamma_prior2$y,
		posterior=model_zpostDS_alpha10[,"gamma_dam_sire"], xlab=expression(gamma)
	)

	pp_plot(x = seq(-5,5,length.out=201), 
		prior = dnorm(seq(-5,5,length.out=201), 0,1),
		posterior=model_zpostDS_alpha1[,"alpha_dam_sire"],xlab=expression(alpha)
	)
	# pp_plot(x = delta_prior3$x, 
	# 	prior = delta_prior3$y,
	# 	posterior=model_zpostDS_alpha1[,"delta_dam_sire"], xlab=expression(delta)
	# )
	pp_plot(x = seq(4,40,length.out=201), 
		prior = dunif(seq(4,40,length.out=201), 4,40),
		posterior=model_zpostDS_alpha1[,"nu_dam_sire"],xlab=expression(nu)
	)
	pp_plot(x = gamma_prior3$x, 
		prior = gamma_prior3$y,
		posterior=model_zpostDS_alpha1[,"gamma_dam_sire"], xlab=expression(gamma)
	)
}
 dev.off()
