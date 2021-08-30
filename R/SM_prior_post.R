
rm(list=ls())

options( stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

library(sn)
library(scales)
library(coda)
library(MCMCvis)
# library(rstan)
library(MCMCglmm)
# library(parallel)

source(paste0(wd,"R/00_functions.R"))


get_gamma<-function(alpha=NULL, delta=NULL, nu){
	if(is.null(delta)) delta = alpha / sqrt(1 + alpha^2); 
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

reduced <- TRUE




for(trait in c("tarsus_mm","headbill_mm","wing_mm","weight_g")){
#trait <- "wing_mm"

	model_files <- list.files(paste0(wd,"Data/Intermediate"))[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_\\d"),list.files(paste0(wd,"Data/Intermediate")))]
	load(paste0(wd,"Data/Intermediate/",model_files[length(model_files)]))


	model_zpost<-do.call(rbind,model_z)[,-(1:7)]


	setEPS()
	pdf(paste0(wd,"R/plots/figure_SM_prior_post_",trait,".pdf"), , height=8, width=15)

	{

	gamma_prior1<-density(get_gamma(delta=runif(10000,-1,1),nu=runif(10000,4,40)), from=-5,to=5,n=101)

	layout(matrix(1:12,byrow=TRUE,nrow=3),width=c(1,2,2,2))
	
	par(mar=c(0,0,0,0))
	blankPlot(); text(0,0,"Nest", cex=2)
	par(mar=c(5,4,2,1))
	pp_plot(x = seq(-1,1,length.out=101), 
		prior = dunif(seq(-1,1,length.out=101), -1,1),
		posterior=model_zpost[,"delta_nest"],xlab=expression(delta)
	)
	pp_plot(x = seq(4,40,length.out=101), 
		prior = dunif(seq(4,40,length.out=101), 4,40),
		posterior=model_zpost[,"nu_nest"],xlab=expression(nu)
	)
	pp_plot(x = gamma_prior1$x, 
		prior = gamma_prior1$y,
		posterior=model_zpost[,"gamma_nest"], xlab=expression(gamma)
	)

	blankPlot(); text(0,0,"Dam/Sire", cex=2)
	pp_plot(x = seq(-1,1,length.out=101), 
		prior = dunif(seq(-1,1,length.out=101), -1,1),
		posterior=model_zpost[,"delta_dam_sire"],xlab=expression(delta)
	)
	pp_plot(x = seq(4,40,length.out=101), 
		prior = dunif(seq(4,40,length.out=101), 4,40),
		posterior=model_zpost[,"nu_dam_sire"],xlab=expression(nu)
	)
	pp_plot(x = gamma_prior1$x, 
		prior = gamma_prior1$y,
		posterior=model_zpost[,"gamma_dam_sire"], xlab=expression(gamma)
	)

	blankPlot(); text(0,0,"Residual", cex=2)
	pp_plot(x = seq(-1,1,length.out=101), 
		prior = dunif(seq(-1,1,length.out=101), -1,1),
		posterior=model_zpost[,if(trait=="weight_g"){"delta_E"}else{"delta_ind"}],xlab=expression(delta)
	)
	pp_plot(x = seq(4,40,length.out=101), 
		prior = dunif(seq(4,40,length.out=101), 4,40),
		posterior=model_zpost[,if(trait=="weight_g"){"nu_E"}else{"nu_ind"}],xlab=expression(nu)
	)
	pp_plot(x = gamma_prior1$x, 
		prior = gamma_prior1$y,
		posterior=model_zpost[,if(trait=="wing_mm"){"gamma_ind"}else{"gamma_E"}], xlab=expression(gamma)
	)
	}
 dev.off()
}

