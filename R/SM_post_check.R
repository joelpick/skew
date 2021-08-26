rm(list=ls())

options( stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")

library(sn)
library(scales)
library(coda)
library(MCMCvis)
# library(rstan)
library(cubature)
library(MCMCglmm)
# library(parallel)

source(paste0(wd,"R/00_functions.R"))

sum_stat <- function(x)c(mean(x),sd(x),stand_skew(x))

se2ci <- function(x) c(x[1],x[1]-x[2]*1.96,x[1]+x[2]*1.96 )

reduced <- TRUE


traits <-  c("tarsus_mm","headbill_mm","weight_g","wing_mm")
traits_lab <- sub("_.+","",traits)
substr(traits_lab,1,1) <- LETTERS[match(substr(traits_lab,1,1),letters)]

# model_files <- list.files(paste0(wd,"Data/Intermediate"))[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(paste0(wd,"Data/Intermediate")))]
# load(paste0(wd,"Data/Intermediate/",model_files[length(model_files)]))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))

load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))


# par(mfrow=c(4,3))

{
setEPS()
pdf(paste0(wd,"R/plots/figure_SM_post_check.pdf"), , height=8, width=15)


layout(matrix(c(1:20),4,5), height=c(1,4,4,4))
par(mar=c(0,0,0,0))
blankPlot();
blankPlot();text(0,0,"Predicted \n Phenotypic Mean", cex=1.5)
blankPlot();text(0,0,"Predicted \n Phenotypic Variance", cex=1.5)
blankPlot();text(0,0,"Predicted \n Phenotypic Skew", cex=1.5)

for(trait in traits){
	#trait <- "weight_g"

	model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zA <- model_z

	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS <- model_z
	
	model_files <- list.files(data_wd)[grep(paste0("stanModNormal2_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zN2 <- model_z



	X<-stan_data_ped_noRep$X
	pred_pos <- which(!grepl((paste(c("year", "sex","timeC"),collapse="|")),colnames(X)))


	model_zpostA<-do.call(rbind,model_zA)[,-(1:7)]
	betasA <- model_zpostA[,grep("beta",colnames(model_zpostA))]
	mu_allA <- apply(betasA, 1, function(x) X%*%x)
	sigmasA <- cbind(sigma_fix=apply(mu_allA,2,sd),model_zpostA[,grep("sigma",colnames(model_zpostA))])
	gammasA <- cbind(gamma_fix=apply(mu_allA,2,stand_skew),model_zpostA[,grep("gamma",colnames(model_zpostA))],gamma_A=0)
	muA_C <- apply(betasA[,pred_pos,drop=FALSE], 1, function(x) X[,pred_pos,drop=FALSE]%*%x)
	sigmasA_C <- cbind(sigma_fix=apply(muA_C,2,sd),model_zpostA[,grep("sigma",colnames(model_zpostA))])
	h2bA <- apply(sigmasA_C, 1, function(x) x["sigma_A"]^2/sum(x^2))

	
	model_zpostDS<-do.call(rbind,model_zDS)[,-(1:7)]
	betasDS <- model_zpostDS[,grep("beta",colnames(model_zpostDS))]
	mu_allDS <- apply(betasDS, 1, function(x) X%*%x)
	sigmasDS <- cbind(sigma_fix=apply(mu_allDS,2,sd),model_zpostDS[,grep("sigma",colnames(model_zpostDS))])
	gammasDS <- cbind(gamma_fix=apply(mu_allDS,2,stand_skew),model_zpostDS[,grep("gamma",colnames(model_zpostDS))])
	muDS_C <- apply(betasDS[,pred_pos,drop=FALSE], 1, function(x) X[,pred_pos,drop=FALSE]%*%x)
	sigmasDS_C <- cbind(sigma_fix=apply(muDS_C,2,sd),model_zpostDS[,grep("sigma",colnames(model_zpostDS))])
	h2bDS <- apply(sigmasDS_C, 1, function(x) (x["sigma_dam_sire"]^2 *4)/(sum(x^2)+x["sigma_dam_sire"]^2*2))



	model_zpostN2<-do.call(rbind,model_zN2)[,-(1:7)]
	betasN2 <- model_zpostN2[,grep("beta",colnames(model_zpostN2))]
	mu_allN2 <- apply(betasN2, 1, function(x) X%*%x)
	sigmasN2 <- cbind(sigma_fix=apply(mu_allN2,2,sd),model_zpostN2[,grep("sigma",colnames(model_zpostN2))])
	gammasN2 <- cbind(gamma_fix=apply(mu_allN2,2,stand_skew),gamma_nest=0,gamma_E=0,gamma_A=0)
	muN2_C <- apply(betasN2[,pred_pos,drop=FALSE], 1, function(x) X[,pred_pos,drop=FALSE]%*%x)
	sigmasN2_C <- cbind(sigma_fix=apply(muN2_C,2,sd),model_zpostN2[,grep("sigma",colnames(model_zpostN2))])


	if(trait!="weight_g") {
		gammasA <- cbind(gammasA,gamma_ME=0)
		colnames(sigmasA)[colnames(sigmasA)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasA)[colnames(sigmasA)=="sigma_ind"] <- "sigma_E"
		colnames(gammasA)[colnames(gammasA)=="gamma_ind"] <- "gamma_E"

		gammasDS <- cbind(gammasDS,gamma_ME=0)
		colnames(sigmasDS)[colnames(sigmasDS)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasDS)[colnames(sigmasDS)=="sigma_ind"] <- "sigma_E"
		colnames(gammasDS)[colnames(gammasDS)=="gamma_ind"] <- "gamma_E"

		# gammasN1 <- cbind(gammasN1,gamma_ME=0)
		# colnames(sigmasN1)[colnames(sigmasN1)=="sigma_E"] <- "sigma_ME"
		# colnames(sigmasN1)[colnames(sigmasN1)=="sigma_ind"] <- "sigma_E"
		
		gammasN2 <- cbind(gammasN2,gamma_ME=0)
		colnames(sigmasN2)[colnames(sigmasN2)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasN2)[colnames(sigmasN2)=="sigma_ind"] <- "sigma_E"
		colnames(sigmasN2_C)[colnames(sigmasN2_C)=="sigma_E"] <- "sigma_ME"
		colnames(sigmasN2_C)[colnames(sigmasN2_C)=="sigma_ind"] <- "sigma_E"


		col_orderA <- c("_fix","_nest", "_A", "_E", "_ME")
		col_orderDS <- c("_fix","_nest", "_dam_sire", "_E", "_ME")
		row_orderA <- c("nest", "animal", "bird_id", "units")
	}else{
		col_orderA <- c("_fix","_nest", "_A", "_E")
		col_orderDS <- c("_fix","_nest", "_dam_sire", "_E")			
		row_orderA <- c("nest", "animal", "units")		
	}

	

	sigmasA <- sigmasA[,paste0("sigma",col_orderA)]
	

	sigmasDS <- sigmasDS[,paste0("sigma",col_orderDS)]
	# sigmasN1 <- sigmasN1[,paste0("sigma",col_orderA)]
	sigmasN2 <- sigmasN2[,paste0("sigma",col_orderA)]
	sigmasN2_C <- sigmasN2_C[,paste0("sigma",col_orderA)]
	h2bN2 <- apply(sigmasN2_C, 1, function(x) x["sigma_A"]^2/sum(x^2))

	gammasA <- gammasA[,paste0("gamma",col_orderA)]
	gammasDS <- gammasDS[,paste0("gamma",col_orderDS)]
	# gammasN1 <- gammasN1[,paste0("gamma",col_orderA)]
	gammasN2 <- gammasN2[,paste0("gamma",col_orderA)]

	post_predA <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_allA)[x],rep(0,if(trait=="weight_g"){3}else{4})), sigmasA[x,], gammasA[x,]))	
	post_predDS <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_allDS)[x],rep(0,if(trait=="weight_g"){3}else{4})), sigmasDS[x,], gammasDS[x,]))
	# post_predN1 <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_allN1)[x],rep(0,if(trait=="weight_g"){3}else{4})), sigmasN1[x,], gammasN1[x,]))
	post_predN2 <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_allN2)[x],rep(0,if(trait=="weight_g"){3}else{4})), sigmasN2[x,], gammasN2[x,]))


	# post_pred2 <- sapply(1:1000,function(x) add_skew(c(colMeans(mu_all2)[x],0,0,0,0), sigmas2[x,], gammas2[x,]))

	post_g <- matrix(rep(sum_stat(THBW[,trait]),3),3,3)
	post_gA <- t(apply(post_predA,1,mean_CI))
	post_gDS <- t(apply(post_predDS,1,mean_CI))
	# post_gN1 <- t(apply(post_predN1,1,mean_CI))
	post_gN2 <- t(apply(post_predN2,1,mean_CI))


	{

	mod_names <-c("G","DS","A")
	mod_col <- c("black","blue","red")
	trait_lab <- traits_lab[which(traits%in%trait)]


	par(mar=c(0,0,0,0))

	blankPlot();text(0,0,trait_lab,cex=1.5)

	par(mar=c(2,4,1,1))

	means <- rbind(post_gN2[1,],post_gDS[1,],post_gA[1,])
	vars <- rbind(post_gN2[2,],post_gDS[2,],post_gA[2,])
	skews <- rbind(post_gN2[3,],post_gDS[3,],post_gA[3,])

	effectPlot(means, xlim=range(c(means,post_g[1,]))+mean(means)*c(-0.01,0.01), names=mod_names, col=mod_col, pch=19)
	abline(v=post_g[1,],lwd=2,lty=2)

	effectPlot(vars, xlim=range(c(vars,post_g[2,]))+mean(vars)*c(-0.05,0.05), names=mod_names, col=mod_col, pch=19)
	abline(v=post_g[2,],lwd=2,lty=2)

	effectPlot(skews, xlim=c(-2.5,0.0), names=mod_names, col=mod_col, pch=19)
	abline(v=post_g[3,],lwd=2,lty=2)
	}
}

 dev.off()
}
