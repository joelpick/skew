rm(list=ls())

options(stringsAsFactors=FALSE)

library(rstan)
library(MCMCglmm)
library(sn)
library(scales)
library(pbapply)

wd <- "~/github/skew/"
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"chick_data.Rdata"))

load(paste0(data_wd,"stanModWing_DS20191210_1850.Rdata"))
load(paste0(data_wd,"stanModWing_pedN20191220_0127.Rdata"))

load(file= paste0(data_wd,"dam_sire_egg.Rdata"))

X <- model.matrix(~ malePresent + clutchSizeC + nestHatchDateC + hatchDay + year + timeC + sex + eggWeightC,THBW)


sumDatOut <- function(stan_mod,X,skew=TRUE,ME=TRUE,DS=TRUE){
#stan_mod=mod_stan_wing_dam
	nTime <- which(colnames(X)=="timeC")

	mod_out <- extract(stan_mod,permute=FALSE)
	betas_all <- mod_out[,,grepl("beta",dimnames(mod_out)[[3]])]
	fixed_sum_all <- apply(betas_all, c(1,2), function(x){
		y <- X[,-nTime]%*%x[-nTime]
		c(var = var(y), skew = stand_skew(y))
	})
	fixed_sum<-apply(fixed_sum_all, c(1), mean_CI)

	out <- list()
	out$fixed <- pars(stan_mod,"beta")
	rownames(out$fixed) <- colnames(X)

	resid_name <- if(ME){"ind"}else{"E"}
	## code below converts residual skew in DS model to residual skew after accounting for Mendelian sampling variance 
	## the problem is that sometimes Ve<0 (if damsire variance is estimates as particularly high in an iteration). Then sqrt(rv)^3 doesnt work, and so cannot get a gamma for that iteration	
	if(DS){
		resid <- apply(mod_out, c(1,2), function(x){ 
			rv <- x[paste0("sigma_",resid_name)]^2-x["sigma_dam_sire"]^2*2
			if(rv<0) rv <- 0
			# from eq.8 in SM
			rs <- (x[ paste0("gamma_",resid_name)]*x[paste0("sigma_",resid_name)]^3) / sqrt(rv)^3
			return(c(rv,rs))
		})
		resid_sum <- apply(resid,1,mean_CI)
		resid_skew <- resid_sum[,2]
 		resid_var <- resid_sum[,1]
	}else{
		resid_var <- pars(stan_mod, paste0("sigma_",resid_name))^2
		resid_skew <- pars(stan_mod, paste0("gamma_",resid_name))
	}

	A_var <- if(DS){
		pars(stan_mod, "sigma_dam_sire")^2*4
	}else{
		pars(stan_mod,"sigma_A")^2}

	out$var <- rbind(
		Fixed = fixed_sum[,1],
		A = as.vector(A_var),
		Nest = as.vector(pars(stan_mod, "sigma_nest")^2),
		Resid = as.vector(resid_var)
		)

	if(ME) out$var <- rbind(out$var,ME = as.vector(pars(stan_mod, "sigma_E")^2))

	if(skew){
		A_skew <- if(DS){pars(stan_mod, "gamma_dam_sire")}else{0}

		out$skew <- rbind(
			Fixed = fixed_sum[,2],
			A = as.vector(A_skew),
			Nest = as.vector(pars(stan_mod, "gamma_nest")),
			Resid = as.vector(resid_skew))
	}

	return(out)
}

SE2CI <- function(mean,se) cbind(mean,mean-se*1.96,mean+se*1.96)

asreml_sumDatOut <- function(asreml_mod,DS=TRUE,ME=TRUE){
	mod <- if(DS){"mod_egg_DS"}else{"mod_egg_A"}
	out <- list()
	out$fixed <- SE2CI(asreml_mod$fixed_mean[mod,],asreml_mod$fixed_SE[mod,])
	out$var <- rbind(fixed=asreml_mod$fixed_var[mod], SE2CI(asreml_mod$vars_mean[mod,],asreml_mod$vars_SE[mod,])[-3,])
	out$skew <- if(ME){
		matrix(0,nrow(out$var)-1,ncol(out$var))
	}else{
		matrix(0,nrow(out$var),ncol(out$var))
	}
	out$skew[1,] <- rep(stand_skew(X %*% asreml_mod$fixed_mean[mod,]),3)
	return(out)
}

skewMod <- sumDatOut(mod_stan_wing_dam,X)
skewModPed <- sumDatOut(mod_stan_wing_pedN,X,DS=FALSE)
asremlModDS <- asreml_sumDatOut(modW)
asremlModPed <- asreml_sumDatOut(modW,DS=FALSE)

setEPS()
pdf(paste0(wd,"Plots/figure_SM_wing_skew.pdf"), height=10, width=10)
{
	layout(matrix(c(1,1,2,3),2), width=c(8,6))
	par(mar=c(3,12,3,1))

	effectPlot(asremlModDS$fixed[-1,],offset=-0.1,col="darkgrey",main="Fixed Effects")
	effectPlot(asremlModPed$fixed[-1,],offset=-0.3, add=TRUE,col="darkgrey",pch=17)
	effectPlot(skewMod$fixed[-1,], add=TRUE,offset=0.3)#,add=TRUE)
	effectPlot(skewModPed$fixed[-1,],offset=0.1,add=TRUE,pch=17)

par(mar=c(3,6,3,1))
	effectPlot(skewMod$var,main="Variance",offset=0.3)#,add=TRUE)
	effectPlot(skewModPed$var,offset=0.1,add=TRUE,pch=17)
	effectPlot(asremlModDS$var,offset=-0.1,add=TRUE,col="darkgrey")
	effectPlot(asremlModPed$var,offset=-0.3,add=TRUE,col="darkgrey",pch=17)

	effectPlot(skewMod$skew, main="Skew",offset=0.3)
	effectPlot(skewModPed$skew,offset=0.1,add=TRUE,pch=17,fixed=2)
	effectPlot(asremlModDS$skew,col="darkgrey", offset=-0.1,add=TRUE,fixed=c(2:4))
	effectPlot(asremlModPed$skew,col="darkgrey", offset=-0.3,add=TRUE,fixed=c(2:4),pch=17)

}
dev.off()


