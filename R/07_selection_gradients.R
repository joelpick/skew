rm(list=ls())

options(stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)
#library(rstan)
library(cubature)
library(parallel)
library(viridis)
library(scales)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))
reduced <- analysis_options$reduced
cond_term<-analysis_options$cond_term  # terms to condition on 
fixed_w <- analysis_options$fixed_w
ncores <- 8
n_it <- 1000
nplot_points<- 100
re_run <- FALSE
save_plot <- TRUE

load(paste0(data_wd,"chick_data.Rdata"))
# THBW_noRep is data excluding repeat measures, used for survival analysis


if(is.null(cond_term)){
	THBW_noRep$category <- as.factor(rep(1, nrow(THBW_noRep)))
}else{
	THBW_noRep$category <- as.factor(apply(THBW_noRep[,cond_term], 1, paste, collapse=""))
}

n_comb<-nlevels(THBW_noRep$category)
# add a category column to the data.frames used in the fitness and trait models (THBW_noRep)
# that designate data into the (n_comb) categories based on the conditioning variables

if(re_run){
	for (trait in c("tarsus_mm","headbill_mm","wing_mm","weight_g")){	

#trait <- "weight_g"
		# load(paste0(data_wd,"day15_survival_model_ME_",trait,".Rdata"))
		# read in fitness models (model_w)
	model_files <- list.files(data_wd)[grep(paste0("day15_survival_ME_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_wpost<-do.call(rbind,model_w)[,-(1:7)]

		z<-THBW_noRep[,trait]
		# trait values


		if(trait=="weight_g"){
			zC<-THBW_noRep[,paste0(trait, "C")]
			# trait values in the survival model (note that it was globally mean centred)

			zmean <- mean(THBW_noRep[,trait])
			# value trait is centered at
		}else{
			
			zmean <- mean(model_wpost[,"beta0_x1"])
			# value trait is centered at

			zC<-THBW_noRep[,trait] - zmean
			# trait values in the survival model (note that it was globally mean centred)
		}


		Wplot_points <- seq(min(z),max(z),length.out=nplot_points)

		## extract relevant parameters from survival model 
		model_w_out<- extract_w2(model_wpost, trait, fixedEffects=formula(paste("~", fixed_w)), data=THBW_noRep)
		attach(model_w_out)

		beta1 <- matrix(NA, n_it, n_comb)  # linear term in linear best fit
		beta2 <- matrix(NA, n_it, n_comb)  # linear term in quadratic best fit
		beta3 <- matrix(NA, n_it, n_comb)  # average gradient
		    S <- matrix(NA, n_it, n_comb)  # selection differential

		dmax<- matrix(NA, n_it, n_comb)    # derivative of the fitness function at the smallest trait value
		dmin<- matrix(NA, n_it, n_comb)    # derivative of the fitness function at the largest trait value

		Wplot<-matrix(NA, n_it, nplot_points)   # fitness function averaged over categories
		WDplot<-matrix(NA, n_it, nplot_points)   # fitness function averaged over categories

		for(i in 1:n_it){
		#i=1
		    beta <- beta_all[i,]
		    gamma <- gamma_all[i,]
		    # trait coefficients from survival model

		    V_nest<-matrix(V_nest_all[i,], 2,2)
		    # nest covariance matrix from survival model

		    ## calculate selection gradients for each combination of conditioning variables
			out_sub <- mclapply(1:n_comb, function(j){
				
		       	cat <- levels(THBW_noRep$category)[j]
		       	# name of the category defined by the conditioning variables 
		   
		       	eta1_sub <- eta1_all[THBW_noRep$category==cat,i]
		       	eta2_sub <- eta2_all[THBW_noRep$category==cat,i]
		       	z_sub <- zC[THBW_noRep$category==cat]
		       	# get (mean-centred) trait and linear predictors for the right category
		   
		   	    mu_etaz <- c(mean(eta1_sub), mean(eta2_sub), mean(z_sub)) 
		   	    V_etaz <- cov(cbind(eta1_sub, eta2_sub, z_sub))
		           # mean and covariance matrices of linear predictor and traits
		   
		   	    W <- sapply(z_sub, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
		   	    # fitness evaluated for each (mean-centred) trait value
		   
		   	    WD <- sapply(z_sub, wD_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
		           # derivative of fitness with respect to  (mean-centred) trait value evaluated for each (mean-centred) trait value
		   
				dmin_sub <- wD_func(min(zC), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
				dmax_sub <- wD_func(max(zC), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
				# slope of the fitness function at extreme (mean-centered) trait values to tell whether there is there an internal stationary point

				Wplot_sub<-sapply(Wplot_points-zmean, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
				# fitness evaluated for (mean-centred) trait values used for visualisation of the fitness function
		   		WDplot_sub<-sapply(Wplot_points-zmean, wD_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)

				mu <- scm(z_sub)

				S_sub<-mean(W*z_sub/mean(W))-mu[1]
				# selection differential calculated conditional on observed (mean-centred) trait value

				C<-mean(W*((z_sub-mu[1])^2)/mean(W))-mu[2]
				# quadratic selection differential calculated conditional on observed (mean-centred) trait value
		   
		   	    beta1_sub<-S_sub/mu[2]  
		   	    beta2_sub<-((mu[4]-mu[2]^2)*S_sub-mu[3]*C)/(mu[2]*(mu[4]-mu[2]^2)-mu[3]^2)
		   	    ## don't know what family does in this function?
		   
		   	    beta3_sub<-mean(WD)/mean(W)

		   	    return(list(S=S_sub,beta1=beta1_sub, beta2=beta2_sub, beta3=beta3_sub, dmin=dmin_sub, dmax=dmax_sub, Wplot=Wplot_sub, WDplot_sub=WDplot_sub))
		   	}, mc.cores = ncores)

			S[i,] <- sapply(out_sub,function(x)x$S)
			beta1[i,] <- sapply(out_sub,function(x)x$beta1)
			beta2[i,] <- sapply(out_sub,function(x)x$beta2)
			beta3[i,] <- sapply(out_sub,function(x)x$beta3)
			dmin[i,] <- sapply(out_sub,function(x)x$dmin)
			dmax[i,] <- sapply(out_sub,function(x)x$dmax)
			Wplot[i,] <- rowMeans(sapply(out_sub,function(x)x$Wplot))
			WDplot[i,] <- rowMeans(sapply(out_sub,function(x)x$WDplot))

			if(((i/n_it)*100) %in% (c(1:10)*10)) cat((i/n_it)*100,"% ",sep="")

		}

		save(S,beta1, beta2, beta3, dmin, dmax, Wplot, WDplot, Wplot_points, file=paste0(data_wd,"selection_gradients_ME_",trait,".Rdata"))
		detach(model_w_out)
	}
}


p_star <- function(x) ifelse(x>0.05,"NS", ifelse(x>0.005,"*","**"))


beta_plot<-function(beta, beta1, beta2, col=1,names=c(expression(beta),expression(beta[1]),expression(beta[2]))){

	beta_all <- rbind(beta=mean_CI(beta),beta_1=mean_CI(beta1),beta_2=mean_CI(beta2))
	beta_p<-c(pMCMC(beta-beta1), pMCMC(beta1-beta2), pMCMC(beta-beta2))

	max_x <- max(beta_all[,1])*2.5
	x <- max_x * c(0.7,0.7,0.9)
	effectPlot(beta_all, xlim=c(0,max_x), col=col,names=names)
	arrows(x,c(3,1.9,3),x,c(2.1,1,1),length=0)
	text(x + max_x*0.075, c(2.5,1.5,2), p_star(beta_p))
}




if(save_plot){
setEPS()
pdf(paste0(wd,"R/plots/figure4_ME.pdf"), height=10, width=8)
}
{
layout(matrix(1:8,ncol=2,byrow=TRUE),width=c(4,2))
	par(mar=c(4,5,2,1),cex.axis=1, cex.lab=1.25)

scales <- c(0.25,0.25,0.4,1.2)
trait <-  c("tarsus_mm","headbill_mm","weight_g","wing_mm")
trait_lab <- paste0(sub("_"," (",trait),")")
trait_lab <- sub("weight","mass",trait_lab)
substr(trait_lab,1,1) <- LETTERS[match(substr(trait_lab,1,1),letters)]
cols <- inferno(5)

internal_optimum <- list()

for(i in 1:length(trait)){
	# i=1
	model_files <- list.files(data_wd)[grep(paste0("day15_survival_ME_",trait[i]),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_wpost<-do.call(rbind,model_w)[,-(1:7)]
	load(paste0(data_wd,"selection_gradients_ME_",trait[i],".Rdata"))

	coef_rows <-grep("beta_x1_hat|beta\\.31|beta\\.32|beta\\.33|beta\\.34",colnames(model_wpost))
	coefs <- t(apply(model_wpost[,coef_rows],2,mean_CI,pMCMC=TRUE))[c(1,3,2,4),]
	mod_pred <- cbind(Wplot_points,t(apply(Wplot,2,mean_CI)))
	colnames(mod_pred) <- c("z","y","y_l","y_u")

	z <- THBW_noRep[,trait[i]]
	z_d <- density(z , adjust=2.5)
	z_d$y <- z_d$y * scales[i]
	zCat <- bin(z,n=10)
	bindedMeans<- aggregate(recruit~zCat,THBW_noRep,mean)
	bindedCounts<- aggregate(recruit~zCat,THBW_noRep,length)

	plot(recruit~zCat,bindedMeans, pch=19, cex=0.05*sqrt(bindedCounts[,2]), ylab="Recruitment probability", ylim=c(0,0.2), xlab=trait_lab[i], xlim=range(z), col=cols[i])

	abline(h=0, v=mean(z), col="grey")
	polygon(z_d$x,z_d$y, col=alpha(cols[i],0.3),border=NA)
	lines(y ~ I(z),mod_pred, col=cols[i])
	lines(y_u ~ I(z),mod_pred, col=cols[i], lty=3)
	lines(y_l ~ I(z),mod_pred, col=cols[i], lty=3)

	text_x <- min(z) + diff(range(z))*c(0.05,0.15,0.25,0.35)
	text(text_x, 0.175, c("F-L","F-Q","R-L","R-Q"))
	text(text_x, 0.15, p_star(coefs[,4]))

	mtext(paste0(letters[i],")"),side=3,adj=0, line=0.5)	

	beta_all <- rbind(beta=mean_CI(rowMeans(beta3)),beta_1=mean_CI(rowMeans(beta1)),beta_2=mean_CI(rowMeans(beta2)))
	beta_p<-c(pMCMC(rowMeans(beta3)-rowMeans(beta1)), pMCMC(rowMeans(beta1)-rowMeans(beta2)), pMCMC(rowMeans(beta3)-rowMeans(beta2)))

	max_x <- max(beta_all[,1])*2.5
	x <- max_x * c(0.7,0.7,0.9)
	effectPlot(beta_all, xlim=c(0,max_x), col=cols[i],names=c(expression(beta),expression(beta[1]),expression(beta[2])))
	arrows(x,c(3,1.9,3),x,c(2.1,1,1),length=0)
	text(x + max_x*0.075, c(2.5,1.5,2), p_star(beta_p))

	mtext(paste0(letters[i+4],")"),side=3,adj=0, line=0.5)		
	internal_optimum[[trait[i]]] <-  sum(rowMeans(dmin) >0 & rowMeans(dmax) <0)/nrow(dmax)

}
}

if(save_plot) dev.off()

save(internal_optimum, file=paste0(data_wd,"internal_optimum_ME.Rdata"), version=2)
)
