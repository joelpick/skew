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
		# row_orderA <- c("nest", "animal", "bird_id", "units")
	}else{
		col_orderA <- c("_fix", "_A","_nest", "_E")
		col_orderDS <- c("_fix", "_dam_sire","_nest", "_E")			
		# row_orderA <- c("nest", "animal", "units")		
	}


	sigmasA <- sigmasA[,paste0("sigma",col_orderA)]
	gammasA <- gammasA[,paste0("gamma",col_orderA)]
	# post_predA <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_allA)[x],rep(0,if(trait=="weight_g"){3}else{4})), sigmasA[x,], gammasA[x,]))	
	# post_gA <- t(apply(post_predA,1,mean_CI))
	stan_betasA <- t(apply(betasA,2,mean_CI))
	stan_varA <- t(apply(sigmasA,2,mean_CI))^2
	stan_gammaA <- t(apply(gammasA,2,mean_CI))
# stan_gammaA <- rbind(stan_gammaA,ME=NA)


	mix_sd<-  sigmasDS[,"sigma_E"]
	sigmasDS <- sigmasDS[,paste0("sigma",col_orderDS)]^2
	#old
	sigmasDS[,"sigma_E"] <-sigmasDS[,"sigma_E"] - sigmasDS[,"sigma_dam_sire"]*2
	sigmasDS[,"sigma_dam_sire"] <-sigmasDS[,"sigma_dam_sire"]*4
	gammasDS <- gammasDS[,paste0("gamma",col_orderDS)]
	gammasDS[,"gamma_E"] <- (gammasDS[,"gamma_E"]*mix_sd^3) / sqrt(sigmasDS[,"sigma_E"])^3

	# post_predDS <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_allDS)[x],rep(0,if(trait=="weight_g"){3}else{4})), sigmasDS[x,], gammasDS[x,]))
	# post_gDS <- t(apply(post_predDS,1,mean_CI))
	stan_betasDS <- t(apply(betasDS,2,mean_CI))
	stan_varDS <- t(apply(sigmasDS,2,mean_CI))
	stan_gammaDS <- t(apply(gammasDS,2,mean_CI))
	
# rv <- sigmasDS[,"sigma_E"]^2-sigmasDS[,"sigma_dam_sire"]^2*2
# # from eq.8 in SM
# mean_CI(rv)
# rs <- (gammasDS[,"gamma_E"]*sigmasDS[,"sigma_E"]^3) / sqrt(rv)^3
# mean_CI(rs)
	## code below converts residual skew in DS model to residual skew after accounting for Mendelian sampling variance 
	## the problem is that sometimes Ve<0 (if damsire variance is estimates as particularly high in an iteration). Then sqrt(rv)^3 doesnt work, and so cannot get a gamma for that iteration	
	# if(DS){
	# 	resid <- apply(mod_out, c(1,2), function(x){ 
	# 		rv <- x["sigma_E"]^2-x["sigma_dam_sire"]^2*2
	# 		if(rv<0) rv <- 0
	# 		# from eq.8 in SM
	# 		rs <- (x[ paste0("gamma_",resid_name)]*x[paste0("sigma_",resid_name)]^3) / sqrt(rv)^3
	# 		return(c(rv,rs))
	# 	})
	# 	resid_sum <- apply(resid,1,mean_CI)
	# 	resid_skew <- resid_sum[,2]
 # 		resid_var <- resid_sum[,1]
	# }else{
	# 	resid_var <- pars(stan_mod, paste0("sigma_",resid_name))^2
	# 	resid_skew <- pars(stan_mod, paste0("gamma_",resid_name))
	# }

	sigmasN2 <- sigmasN2[,paste0("sigma",col_orderA)]
	gammasN2 <- gammasN2[,paste0("gamma",col_orderA)]
	stan_betasN2 <- t(apply(betasN2,2,mean_CI))
	stan_varN2 <- t(apply(sigmasN2,2,mean_CI))^2
	stan_gammaN2 <- t(apply(gammasN2,2,mean_CI))


	 # post_g <- matrix(rep(sum_stat(THBW[,trait]),3),3,3)
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

	# par(mar=c(2,6,1,1))
	# effectPlot(stan_betasDS[-1,], col=trait_col,pch=21, bg=alpha(trait_col,0.5), names=beta_names[-1])
	# effectPlot(stan_betasA[-1,], add=TRUE, offset=-0.15, col=trait_col)
	# effectPlot(stan_betasN2[-1,], add=TRUE, offset=0.15, col=trait_col,pch=21, bg="white")

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













# se2ci <- function(x) c(x[1],x[1]-x[2]*1.96,x[1]+x[2]*1.96 )

# stPlot <- function(x, mu, dp, col, xlim, r2=1,...){
# 	plot(I(r2*dst(x, dp=dp))~I(x+mu), type="l", xaxt="n", yaxt="n", bty="n", xlim=xlim, col=col,...)
# 	 abline(h=0, col="grey")
# 	 polygon(I(x+mu),I(r2*dst(x,dp=dp)), col=alpha(col,0.5), border=NA)
# }
# stPlot_fixedy <- function(x, mu, dp, col, xlim, ylim=c(0,1),...){
# 	y <- dst(x, dp=dp)
# 	plot(I(y/max(y))~I(x+mu), type="l", xaxt="n", yaxt="n", bty="n", xlim=xlim, col=col,ylim=ylim, xlab="", ylab="")
# 	 abline(h=0, col="grey")
# 	 polygon(I(x+mu),I(y/max(y)), col=alpha(col,0.5), border=NA)
# }

# ## plotting variables

# x<- seq(-4,4 ,length.out=500)
# cex.axis=1
# ylims1=c(0,1.2)
# cex.heading = 2
# ##plots

# cols <- c(1,2)

# setEPS()
# pdf(paste0(wd,"R/plots/figure3.pdf"), height=6, width=12)
# {
# layout(matrix(1:18,nrow=3,byrow=FALSE), height=c(2,rep(4,2)))
# # layout(matrix(1:18,ncol=3,byrow=TRUE), width=c(2,rep(4,2)))

# #layout(matrix(1:18,nrow=6,byrow=TRUE), width=c(2,rep(4,2)))
# #layout(matrix(1:18,ncol=6,byrow=FALSE), height=c(2,rep(4,2)))
# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Phenotype",cex=cex.heading)
# par(mar=c(3,1,0,1),cex.axis=cex.axis)
# tarsus_d <- density(THBW$tarsus_mm , adjust=2.5 )
# tarsus_d$y <- tarsus_d$y /max(tarsus_d$y)
# plot(tarsus_d$x,tarsus_d$y, type="l", bty="n", xlim=range(THBW$tarsus_mm), ylim=ylims1,  yaxt="n", xlab="", ylab="",col=cols[1]); polygon(tarsus_d$x,tarsus_d$y, col=alpha(cols[1],0.5), border=NA);abline(h=0, col="grey");text(14,0.8,"Tarsus", cex=cex.heading)
# abline(v=mean(THBW$tarsus_mm),col="grey")

# # HB_d <- density(THBW$headbill_mm , adjust=2.5 )
# # HB_d$y <- HB_d$y /max(HB_d$y)
# # plot(HB_d$x,HB_d$y, type="l", bty="n", xlim=range(THBW$headbill_mm), ylim=ylims1,  yaxt="n", xlab="", ylab="",col=3); polygon(HB_d$x,HB_d$y, col=alpha(3,0.5), border=NA);abline(h=0, col="grey");text(19,0.8,"Headbill", cex=cex.heading)
# # abline(v=mean(THBW$headbill_mm),col="grey")

# # wing_d <- density(THBW$wing_mm , adjust=2.5 )
# # wing_d$y <- wing_d$y /max(wing_d$y)
# # plot(wing_d$x,wing_d$y, type="l", bty="n", xlim=range(THBW$wing_mm), ylim=ylims1,  yaxt="n", xlab="", ylab="",col=4); polygon(wing_d$x,wing_d$y, col=alpha(4,0.5), border=NA);abline(h=0, col="grey");text(21,0.8,"Wing", cex=cex.heading)
# # abline(v=mean(THBW$wing_mm),col="grey")

# weight_d <- density(THBW$weight_g , adjust=2.5 )
# weight_d$y <- weight_d$y /max(weight_d$y)
# plot(weight_d$x,weight_d$y, type="l", bty="n", xlim=range(THBW$weight_g), ylim=ylims1, yaxt="n", xlab="", ylab="",col=cols[2]); polygon(weight_d$x,weight_d$y, col=alpha(cols[2],0.5), border=NA);text(6.5,0.8,"Weight", cex=cex.heading)
# abline(v=mean(THBW$weight_g),col="grey")

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Additive\nGenetic",cex=cex.heading)
# par(mar=c(3,1,0,1))
# stPlot_fixedy(x, mu=muT, dp=c(0, varT[3,1], 0, 1000), xlim=range(THBW$tarsus_mm), col=cols[1], ylim=ylims1)
# # stPlot_fixedy(x, mu=muHB, dp=c(0, varHB[3,1], 0, 1000), xlim=range(THBW$headbill_mm), col=3, ylim=ylims1)
# # stPlot_fixedy(x*4, mu=muW, dp=c(0, varW[3,1], 0, 1000), xlim=range(THBW$wing_mm), col=4, ylim=ylims1)
# stPlot_fixedy(x, mu=muM, dp=c(0, varM[3,1], 0, 1000), xlim=range(THBW$weight_g), col=cols[2], ylim=ylims1)

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Nest",cex=cex.heading)
# par(mar=c(3,1,0,1))
# stPlot_fixedy(x, mu=muT, dp=nestT[,1], xlim=range(THBW$tarsus_mm), col=cols[1], cex.axis=cex.axis, ylim=ylims1)
# # stPlot_fixedy(x, mu=muHB, dp=nestHB[,1], xlim=range(THBW$headbill_mm), col=3, cex.axis=cex.axis, ylim=ylims1)
# # stPlot_fixedy(x*4, mu=muW, dp=nestW[,1], xlim=range(THBW$wing_mm), col=4, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x, mu=muM, dp=nestM[,1], xlim=range(THBW$weight_g), col=cols[2], cex.axis=cex.axis, ylim=ylims1)

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Residual",cex=cex.heading)
# par(mar=c(3,1,0,1))
# stPlot_fixedy(x, mu=muT, dp=residT[,1], xlim=range(THBW$tarsus_mm), col=cols[1], cex.axis=cex.axis, ylim=ylims1)
# # stPlot_fixedy(x, mu=muHB, dp=residHB[,1], xlim=range(THBW$headbill_mm), col=3, cex.axis=cex.axis, ylim=ylims1)
# # stPlot_fixedy(x*4, mu=muW, dp=residW[,1], xlim=range(THBW$wing_mm), col=4, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x, mu=muM, dp=residM[,1], xlim=range(THBW$weight_g), col=cols[2], cex.axis=cex.axis, ylim=ylims1)

# ## make it proportion
# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Variance",cex=cex.heading)
# par(mar=c(3,1,1,3))
# effectPlot(varT, col=cols[1], xlim=c(0,0.2), cex.axis=cex.axis,yaxis=4)
# # effectPlot(varHB, col=3, xlim=c(0,0.25), cex.axis=cex.axis,yaxis=4)
# # effectPlot(varW, col=4, xlim=c(0,9), cex.axis=cex.axis,yaxis=4)
# effectPlot(varM, col=cols[2], xlim=c(0,0.9), cex.axis=cex.axis,yaxis=4)

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Standardised\nSkew",cex=cex.heading)
# par(mar=c(3,3,1,1))
# effectPlot(skewT, col=cols[1], xlim=c(-6,1.5), cex.axis=cex.axis)
# # effectPlot(skewHB, col=3, xlim=c(-6,1.5), cex.axis=cex.axis)
# # effectPlot(skewW, col=4, xlim=c(-6,1.5), cex.axis=cex.axis)
# effectPlot(skewM, col=cols[2], xlim=c(-6,1.5), cex.axis=cex.axis)

# }
# dev.off()


# # setEPS()
# # pdf(paste0(wd,"R/plots/figure3.pdf"), height=7.5, width=12)
# {
# layout(matrix(1:30,nrow=5,byrow=FALSE), height=c(2,rep(4,5)))
# # layout(matrix(1:30,ncol=5,byrow=TRUE), width=c(2,rep(4,5)))

# #layout(matrix(1:18,nrow=6,byrow=TRUE), width=c(2,rep(4,2)))
# #layout(matrix(1:18,ncol=6,byrow=FALSE), height=c(2,rep(4,2)))
# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Phenotype",cex=cex.heading)
# par(mar=c(3,1,0,1),cex.axis=cex.axis)
# tarsus_d <- density(THBW$tarsus_mm , adjust=2.5 )
# tarsus_d$y <- tarsus_d$y /max(tarsus_d$y)
# plot(tarsus_d$x,tarsus_d$y, type="l", bty="n", xlim=range(THBW$tarsus_mm), ylim=ylims1,  yaxt="n", xlab="", ylab=""); polygon(tarsus_d$x,tarsus_d$y, col=alpha(1,0.5));abline(h=0, col="grey");text(14,0.8,"Tarsus", cex=cex.heading)
# abline(v=mean(THBW$tarsus_mm),col="grey")

# HB_d <- density(THBW$headbill_mm , adjust=2.5 )
# HB_d$y <- HB_d$y /max(HB_d$y)
# plot(HB_d$x,HB_d$y, type="l", bty="n", xlim=range(THBW$headbill_mm), ylim=ylims1,  yaxt="n", xlab="", ylab="",col=3); polygon(HB_d$x,HB_d$y, col=alpha(3,0.5));abline(h=0, col="grey");text(19,0.8,"Headbill", cex=cex.heading)
# abline(v=mean(THBW$headbill_mm),col="grey")

# wing_d <- density(THBW$wing_mm , adjust=2.5 )
# wing_d$y <- wing_d$y /max(wing_d$y)
# plot(wing_d$x,wing_d$y, type="l", bty="n", xlim=range(THBW$wing_mm), ylim=ylims1,  yaxt="n", xlab="", ylab="",col=4); polygon(wing_d$x,wing_d$y, col=alpha(4,0.5));abline(h=0, col="grey");text(21,0.8,"Wing", cex=cex.heading)
# abline(v=mean(THBW$wing_mm),col="grey")

# weight_d <- density(THBW$weight_g , adjust=2.5 )
# weight_d$y <- weight_d$y /max(weight_d$y)
# plot(weight_d$x,weight_d$y, type="l", bty="n", xlim=range(THBW$weight_g), ylim=ylims1, yaxt="n", xlab="", ylab="",col=2); polygon(weight_d$x,weight_d$y, col=alpha(2,0.5), border=NA);text(6.5,0.8,"Weight", cex=cex.heading)
# abline(v=mean(THBW$weight_g),col="grey")

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Additive\nGenetic",cex=cex.heading)
# par(mar=c(3,1,0,1))
# stPlot_fixedy(x, mu=muT, dp=c(0, varT[3,1], 0, 1000), xlim=range(THBW$tarsus_mm), col=1, ylim=ylims1)
# stPlot_fixedy(x, mu=muHB, dp=c(0, varHB[3,1], 0, 1000), xlim=range(THBW$headbill_mm), col=3, ylim=ylims1)
# stPlot_fixedy(x*4, mu=muW, dp=c(0, varW[3,1], 0, 1000), xlim=range(THBW$wing_mm), col=4, ylim=ylims1)
# stPlot_fixedy(x, mu=muM, dp=c(0, varM[3,1], 0, 1000), xlim=range(THBW$weight_g), col=2, ylim=ylims1)

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Nest",cex=cex.heading)
# par(mar=c(3,1,0,1))
# stPlot_fixedy(x, mu=muT, dp=nestT[,1], xlim=range(THBW$tarsus_mm), col=1, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x, mu=muHB, dp=nestHB[,1], xlim=range(THBW$headbill_mm), col=3, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x*4, mu=muW, dp=nestW[,1], xlim=range(THBW$wing_mm), col=4, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x, mu=muM, dp=nestM[,1], xlim=range(THBW$weight_g), col=2, cex.axis=cex.axis, ylim=ylims1)

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Residual",cex=cex.heading)
# par(mar=c(3,1,0,1))
# stPlot_fixedy(x, mu=muT, dp=residT[,1], xlim=range(THBW$tarsus_mm), col=1, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x, mu=muHB, dp=residHB[,1], xlim=range(THBW$headbill_mm), col=3, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x*4, mu=muW, dp=residW[,1], xlim=range(THBW$wing_mm), col=4, cex.axis=cex.axis, ylim=ylims1)
# stPlot_fixedy(x, mu=muM, dp=residM[,1], xlim=range(THBW$weight_g), col=2, cex.axis=cex.axis, ylim=ylims1)

# ## make it proportion
# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Variance",cex=cex.heading)
# par(mar=c(3,1,1,3))
# effectPlot(varT, col=1, xlim=c(0,0.2), cex.axis=cex.axis,yaxis=4)
# effectPlot(varHB, col=3, xlim=c(0,0.25), cex.axis=cex.axis,yaxis=4)
# effectPlot(varW, col=4, xlim=c(0,9), cex.axis=cex.axis,yaxis=4)
# effectPlot(varM, col=2, xlim=c(0,0.9), cex.axis=cex.axis,yaxis=4)

# par(mar=c(0,0,0,0)); blankPlot(); text(0,0,"Standardised\nSkew",cex=cex.heading)
# par(mar=c(3,3,1,1))
# effectPlot(skewT, col=1, xlim=c(-6,1.5), cex.axis=cex.axis)
# effectPlot(skewHB, col=3, xlim=c(-6,1.5), cex.axis=cex.axis)
# effectPlot(skewW, col=4, xlim=c(-6,1.5), cex.axis=cex.axis)
# effectPlot(skewM, col=2, xlim=c(-6,1.5), cex.axis=cex.axis)

# }
# # dev.off()
