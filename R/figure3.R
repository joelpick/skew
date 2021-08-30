rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

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


setEPS()
pdf(paste0(wd,"R/plots/figure3.pdf"), height=7.5, width=15)

{
layout(matrix(1:20,ncol=5,byrow=FALSE), width=c(1,rep(3,4)))#, height=c(1.5,1,1)

traits <-  c("tarsus_mm","headbill_mm","weight_g","wing_mm")
traits_lab <- paste0(sub("_"," (",traits),")")
traits_lab <- sub("weight","mass",traits_lab)
substr(traits_lab,1,1) <- LETTERS[match(substr(traits_lab,1,1),letters)]
cex.axis=1
ylims1=c(0,1.2)
cex.heading = 1.75

par(mar=c(0,0,0,0))
blankPlot(); #text(0,0,"Phenotype",cex=cex.heading)
blankPlot(); text(0,0,"Variance",cex=cex.heading)
blankPlot(); text(0,0,"Skew",cex=cex.heading)
blankPlot(); text(0,0,"Parent-\nOffpsring\nRegression",cex=cex.heading)


for(trait in traits){
#trait <- "tarsus_mm"

	X<-stan_data_ped_noRep$X
	
	# if(trait!="wing_mm") 
		load(paste0(wd,"Data/Intermediate/nonLinearPO_",trait,".Rdata"))

	model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	## posterior distributions for stan model for the trait (model_z), as an mcmc.list object
	model_zA <- model_z

	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_\\d"),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS <- model_z

	load(paste0(data_wd,"nonLinearPO_",trait,".Rdata"))

	model_zpostA<-do.call(rbind,model_zA)[,-(1:7)]
	betasA <- model_zpostA[,grep("beta",colnames(model_zpostA))]
	mu_allA <- apply(betasA, 1, function(x) X%*%x)
	sigmasA <- cbind(sigma_fix=apply(mu_allA,2,sd),model_zpostA[,grep("sigma",colnames(model_zpostA))])
	gammasA <- cbind(gamma_fix=apply(mu_allA,2,stand_skew),model_zpostA[,grep("gamma",colnames(model_zpostA))],gamma_A=0)


	if(trait!="weight_g") {
		gammasA <- cbind(gammasA,gamma_ME=0)
		colnames(sigmasA)[colnames(sigmasA)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasA)[colnames(sigmasA)=="sigma_ind"] <- "sigma_E"
		colnames(gammasA)[colnames(gammasA)=="gamma_ind"] <- "gamma_E"

		col_orderA <- c("_fix", "_A","_nest", "_E", "_ME")
		col_orderDS <- c("_fix", "_dam_sire","_nest", "_E", "_ME")
	}else{
		col_orderA <- c("_fix", "_A","_nest", "_E")
		col_orderDS <- c("_fix", "_dam_sire","_nest", "_E")			
	}


	sigmasA <- sigmasA[,paste0("sigma",col_orderA)]
	gammasA <- gammasA[,paste0("gamma",col_orderA)]

	stan_varA <- t(apply(sigmasA,2,mean_CI))^2
	stan_gammaA <- t(apply(gammasA,2,mean_CI))
	stan_gammaA["gamma_A",] <- as.numeric(pars(model_zDS, "gamma_dam_sire"))

	z <- THBW[,trait]

	x<- seq(-4,4 ,length.out=500)
	trait_col <- inferno(5)[which(traits%in%trait)]
    trait_lab <- traits_lab[which(traits%in%trait)]
	text_x <- min(z) + (max(z)-min(z))*0.2
	y_names <- c("Fixed", "Genetic","Nest", "Residual", "ME")
	if(trait=="weight_g")y_names <- y_names[1:4]
	y <- if(trait=="weight_g"){ 5:2 }else{ 5:1 }


	
	par(mar=c(3,1,1,1),cex.axis=cex.axis)
	# trait_d <- if(trait!="wing_mm") {
		trait_d <- density(z_corrected, from=POreg_out[1,"z"],to=max(POreg_out[,"z"]),n=100,adjust=2)
	# }else{
		# density(z , adjust=2 )
	# }
	# trait_d <- density(z , adjust=2.5 )
	# trait_d$y <- trait_d$y /max(trait_d$y)
	ylims1=c(0,max(trait_d$y)*1.2)

	plot(trait_d$x,trait_d$y, type="l", bty="n", xlim=range(z), ylim=ylims1,  yaxt="n", xlab="", ylab="",col=trait_col)
	polygon(trait_d$x,trait_d$y, col=alpha(trait_col,0.5), border=NA)
	abline(h=0, col="grey")
	text(text_x,max(trait_d$y)*0.8,trait_lab, cex=cex.heading)
	# text(text_x,0.8,trait_lab, cex=cex.heading)
	text(text_x,max(trait_d$y)*0.65,paste0("skew = ",round(stand_skew(THBW[,trait]),3)), cex=cex.heading*0.8)
	abline(v=mean(z),col="grey")
	# if(trait!="wing_mm")
	 lines(dz_p~z,POreg_out,col="red",lwd=2)


	par(mar=c(2,6,1,1))
	effectPlot(stan_varA[,],y, col=trait_col,names=y_names, ylim=c(0.5,5.5))

	effectPlot(stan_gammaA[-5,], xlim=c(-4,1), col=trait_col,names=y_names[1:4])

	par(mar=c(1,3,1,1))

	plot(I(POreg_out[,"Eg_p"] + mean(z))~POreg_out[,"z"], ylim=c(mean(z)-sd(z),mean(z)+1.5*sd(z)),xlim=range(z),  xlab="", ylab= "", lwd=2, col=trait_col, type="l", xaxt="n", yaxt="n", bty="l")


}

}
dev.off()






##### all models for supplements:
#sum_stat <- function(x)c(mean(x),sd(x),stand_skew(x))

setEPS()
pdf(paste0(wd,"R/plots/figure_SM_skewt_Gaussian.pdf"), height=7.5, width=15)

{
	# layout(matrix(1:20,ncol=5,byrow=FALSE), width=c(1,rep(3,4)), height=c(1,10,5,4))
layout(matrix(1:15,ncol=5,byrow=FALSE), width=c(1,rep(3,4)), height=c(1,5,4))
traits <-  c("tarsus_mm","headbill_mm","weight_g","wing_mm")
traits_lab <- paste0(sub("_"," (",traits),")")
traits_lab <- sub("weight","mass",traits_lab)
substr(traits_lab,1,1) <- LETTERS[match(substr(traits_lab,1,1),letters)]
cex.axis=1
ylims1=c(0,1.2)
cex.heading = 1.75

par(mar=c(0,0,0,0))
blankPlot()
# blankPlot(); text(0,0,"Fixed\nEffects",cex=cex.heading)

blankPlot(); text(0,0,"Variance",cex=cex.heading)
blankPlot(); text(0,0,"Skew",cex=cex.heading)


for(trait in traits){
#trait <- "wing_mm"

	X<-stan_data_ped_noRep$X
	
	# if(trait!="wing_mm") 
		load(paste0(wd,"Data/Intermediate/nonLinearPO_",trait,".Rdata"))

	model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	## posterior distributions for stan model for the trait (model_z), as an mcmc.list object
	model_zA <- model_z

	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_\\d"),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS <- model_z

	model_files <- list.files(data_wd)[grep(paste0("stanModNormal2_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zN2 <- model_z
	

	model_zpostA<-do.call(rbind,model_zA)[,-(1:7)]
	betasA <- model_zpostA[,grep("beta",colnames(model_zpostA))]
	mu_allA <- apply(betasA, 1, function(x) X%*%x)
	sigmasA <- cbind(sigma_fix=apply(mu_allA,2,sd),model_zpostA[,grep("sigma",colnames(model_zpostA))])
	gammasA <- cbind(gamma_fix=apply(mu_allA,2,stand_skew),model_zpostA[,grep("gamma",colnames(model_zpostA))],gamma_A=0)
	
	model_zpostDS<-do.call(rbind,model_zDS)[,-(1:7)]
	betasDS <- model_zpostDS[,grep("beta",colnames(model_zpostDS))]
	mu_allDS <- apply(betasDS, 1, function(x) X%*%x)
	sigmasDS <- cbind(sigma_fix=apply(mu_allDS,2,sd),model_zpostDS[,grep("sigma",colnames(model_zpostDS))])
	gammasDS <- cbind(gamma_fix=apply(mu_allDS,2,stand_skew),model_zpostDS[,grep("gamma",colnames(model_zpostDS))])

	model_zpostN2<-do.call(rbind,model_zN2)[,-(1:7)]
	betasN2 <- model_zpostN2[,grep("beta",colnames(model_zpostN2))]
	mu_allN2 <- apply(betasN2, 1, function(x) X%*%x)
	sigmasN2 <- cbind(sigma_fix=apply(mu_allN2,2,sd),model_zpostN2[,grep("sigma",colnames(model_zpostN2))])
	gammasN2 <- cbind(gamma_fix=apply(mu_allN2,2,stand_skew),gamma_nest=0,gamma_E=0,gamma_A=0)

	if(trait!="weight_g") {
		gammasA <- cbind(gammasA,gamma_ME=0)
		colnames(sigmasA)[colnames(sigmasA)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasA)[colnames(sigmasA)=="sigma_ind"] <- "sigma_E"
		colnames(gammasA)[colnames(gammasA)=="gamma_ind"] <- "gamma_E"

		gammasDS <- cbind(gammasDS,gamma_ME=0)
		colnames(sigmasDS)[colnames(sigmasDS)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasDS)[colnames(sigmasDS)=="sigma_ind"] <- "sigma_E"
		colnames(gammasDS)[colnames(gammasDS)=="gamma_ind"] <- "gamma_E"
		
		gammasN2 <- cbind(gammasN2,gamma_ME=0)
		colnames(sigmasN2)[colnames(sigmasN2)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasN2)[colnames(sigmasN2)=="sigma_ind"] <- "sigma_E"

		col_orderA <- c("_fix", "_A","_nest", "_E", "_ME")
		col_orderDS <- c("_fix", "_dam_sire","_nest", "_E", "_ME")
	}else{
		col_orderA <- c("_fix", "_A","_nest", "_E")
		col_orderDS <- c("_fix", "_dam_sire","_nest", "_E")			
	}


	sigmasA <- sigmasA[,paste0("sigma",col_orderA)]
	gammasA <- gammasA[,paste0("gamma",col_orderA)]
	stan_betasA <- t(apply(betasA,2,mean_CI))
	stan_varA <- t(apply(sigmasA,2,mean_CI))^2
	stan_gammaA <- t(apply(gammasA,2,mean_CI))


	mix_sd<-  sigmasDS[,"sigma_E"]
	sigmasDS <- sigmasDS[,paste0("sigma",col_orderDS)]^2
	#old
	sigmasDS[,"sigma_E"] <-sigmasDS[,"sigma_E"] - sigmasDS[,"sigma_dam_sire"]*2
	sigmasDS[,"sigma_dam_sire"] <-sigmasDS[,"sigma_dam_sire"]*4
	gammasDS <- gammasDS[,paste0("gamma",col_orderDS)]
	gammasDS[,"gamma_E"] <- (gammasDS[,"gamma_E"]*mix_sd^3) / sqrt(sigmasDS[,"sigma_E"])^3

	stan_betasDS <- t(apply(betasDS,2,mean_CI))
	stan_varDS <- t(apply(sigmasDS,2,mean_CI))
	stan_gammaDS <- t(apply(gammasDS,2,mean_CI))


	sigmasN2 <- sigmasN2[,paste0("sigma",col_orderA)]
	gammasN2 <- gammasN2[,paste0("gamma",col_orderA)]
	stan_betasN2 <- t(apply(betasN2,2,mean_CI))
	stan_varN2 <- t(apply(sigmasN2,2,mean_CI))^2
	stan_gammaN2 <- t(apply(gammasN2,2,mean_CI))

	z <- THBW[,trait]

	x<- seq(-4,4 ,length.out=500)
	trait_col <- inferno(5)[which(traits%in%trait)]
    trait_lab <- traits_lab[which(traits%in%trait)]
	text_x <- min(z) + (max(z)-min(z))*0.2
	y_names <- c("Fixed", "Genetic","Nest", "Residual", "ME")
	beta_names <- c( "Intercept", "Year (12)","Year (13)","Year (14)","Year (15)","Year (16)","Year (17)","Year (18)","Time","Sex (M)", "Egg Weight")

	if(trait=="weight_g")y_names <- y_names[1:4]
	y <- if(trait=="weight_g"){ 5:2 }else{ 5:1 }


	par(mar=c(0,0,0,0),cex.axis=cex.axis)
	blankPlot(); text(0,0,trait_lab, cex=cex.heading)

	par(mar=c(2,6,1,1))
	effectPlot(stan_varDS[,],y, col=trait_col,pch=21, bg=alpha(trait_col,0.5),names=y_names, ylim=c(0.5,5.5))
	effectPlot(stan_varA[,],y, add=TRUE, offset=-0.15, col=trait_col)
	effectPlot(stan_varN2[,],y, add=TRUE, offset=0.15, col=trait_col,pch=21, bg="white")

	effectPlot(stan_gammaDS[-5,], xlim=c(-4,1), pch=21, bg=alpha(trait_col,0.5), col=trait_col,names=y_names[1:4])
	effectPlot(stan_gammaA[-5,], add=TRUE, offset=-0.15, col=trait_col)
		effectPlot(stan_gammaN2[-5,], add=TRUE, offset=0.15, col=trait_col,pch=21, bg="white")

}

}
dev.off()
