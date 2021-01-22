rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

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


# trait<-"headbill_mm"
#cont_term<-c("timeC")   	# terms to control for
#model_moments<-FALSE 		# when calculating selection gradients should model-based or sample moments be used




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


		load(paste0(data_wd,"day15_survival_model_",trait,".Rdata"))
		# read in fitness models (model_w)

		# load(paste0(wd,"Data/Intermediate/Wplot_points.Rdata"))
		# trait values at which to evaluate the fitness function for visualisation purposes
		# Wplot_points<-get(paste0("Wplot_points_",trait))
		# nplot_points<-length(Wplot_points)


		z<-THBW_noRep[,trait]
		# trait values

		zC<-THBW_noRep[,paste0(trait, "C")]
		# trait values in the survival model (note that it was globally mean centred)

		zmean <- mean(THBW_noRep[,trait])
		# value trait is centered at

		Wplot_points <- seq(min(z),max(z),length.out=nplot_points)

		## extract relevant parameters from survival model 
		attach(extract_w(model_w, trait, fixedEffects=formula(paste("~", fixed_w)), data=THBW_noRep))

		beta1 <- matrix(NA, n_it, n_comb)  # linear term in linear best fit
		beta2 <- matrix(NA, n_it, n_comb)  # linear term in quadratic best fit
		beta3 <- matrix(NA, n_it, n_comb)  # average gradient
		    S <- matrix(NA, n_it, n_comb)  # selection differential

		dmax<- matrix(NA, n_it, n_comb)    # derivative of the fitness function at the smallest trait value
		dmin<- matrix(NA, n_it, n_comb)    # derivative of the fitness function at the largest trait value

		Wplot<-matrix(NA, n_it, nplot_points)   # fitness function averaged over categories

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
		   
				mu <- scm(z_sub)

				S_sub<-mean(W*z_sub/mean(W))-mu[1]
				# selection differential calculated conditional on observed (mean-centred) trait value

				C<-mean(W*((z_sub-mu[1])^2)/mean(W))-mu[2]
				# quadratic selection differential calculated conditional on observed (mean-centred) trait value
		   
		   	    beta1_sub<-S_sub/mu[2]  
		   	    beta2_sub<-betaLA_2(mu, S=S_sub, C=C, family="ST")
		   	    ## don't know what family does in this function?
		   
		   	    beta3_sub<-mean(WD)/mean(W)

		   	    return(list(S=S_sub,beta1=beta1_sub, beta2=beta2_sub, beta3=beta3_sub, dmin=dmin_sub, dmax=dmax_sub, Wplot=Wplot_sub))
		   	}, mc.cores = ncores)

			S[i,] <- sapply(out_sub,function(x)x$S)
			beta1[i,] <- sapply(out_sub,function(x)x$beta1)
			beta2[i,] <- sapply(out_sub,function(x)x$beta2)
			beta3[i,] <- sapply(out_sub,function(x)x$beta3)
			dmin[i,] <- sapply(out_sub,function(x)x$dmin)
			dmax[i,] <- sapply(out_sub,function(x)x$dmax)
			Wplot[i,] <- rowMeans(sapply(out_sub,function(x)x$Wplot))

			if(((i/n_it)*100) %in% (c(1:10)*10)) cat((i/n_it)*100,"% ",sep="")

		}

		save(S,beta1,beta2, beta3, dmin, dmax, Wplot, Wplot_points, file=paste0(data_wd,"selection_gradients_",trait,".Rdata"))
	}
}




pMCMC <- function(x) 2*pmax(0.5/length(x), pmin(sum(x > 0)/length(x), 1 - sum(x > 0)/length(x)))

p_star <- function(x) ifelse(x>0.05,"NS", ifelse(x>0.005,"*","**"))


surv_plot <- function(trait, surivial, mod_pred, coefs, xlab="", ylab="Recruitment probability", ylim=c(0,0.2), text=FALSE, col=1, scale=0.25){
	trait_d <- density(trait , adjust=2.5)
	trait_d$y <- trait_d$y * scale
	traitCat <- bin(trait,n=10)
	binPlot(surivial~traitCat, xlab=xlab, ylab=ylab, ylim=ylim, xlim=range(trait),text.cex=0.7, text=text, col=col)
	abline(h=0, v=mean(trait), col="grey")
	polygon(trait_d$x,trait_d$y, col=alpha(col,0.3),border=NA)
	lines(y ~ I(z),mod_pred, col=col)
	lines(y_u ~ I(z),mod_pred, col=col, lty=3)
	lines(y_l ~ I(z),mod_pred, col=col, lty=3)

	text_x <- min(trait) + diff(range(trait))*c(0.05,0.15,0.25,0.35)
	text(text_x, 0.175, c("F-L","F-Q","R-L","R-Q"))
	text(text_x, 0.15, p_star(coefs[,4]))

}


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
pdf(paste0(wd,"R/plots/figure4.pdf"), height=12, width=7)
}
{
layout(matrix(1:8,ncol=2,byrow=TRUE),width=c(4,2))
	par(mar=c(5,5,1,1),cex.axis=1, cex.lab=1.25)

scales <- c(0.25,0.25,1.2,0.4)
trait <-  c("tarsus_mm","headbill_mm","wing_mm","weight_g")
trait_lab <- paste0(sub("_"," (",trait),")")
substr(trait_lab,1,1) <- LETTERS[match(substr(trait_lab,1,1),letters)]

for(i in 1:length(trait)){
	load(paste0(data_wd,"day15_survival_model_",trait[i],".Rdata"))
	load(paste0(data_wd,"selection_gradients_",trait[i],".Rdata"))
	T_rows <-grep(trait[i],rownames(summary(model_w)$sol))
	T_sol <- summary(model_w)$sol[T_rows,c(1,2,3,5)]
	t_pred <- cbind(Wplot_points,t(apply(Wplot,2,mean_CI)))
	colnames(t_pred) <- c("z","y","y_l","y_u")

	surv_plot(THBW_noRep[,trait[i]],THBW_noRep$recruit,t_pred,T_sol, xlab=trait_lab[i], scale=scales[i], col=inferno(5)[i])
	beta_plot(rowMeans(beta3),rowMeans(beta1),rowMeans(beta2), col=inferno(5)[i])
}
}

if(save_plot) dev.off()






# load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))
# # data used in stan models, stan_data_ped_noRep$X is the design matrix (without repeat measures)

# # read in trait models (model_z):
# model_files <- list.files(paste0(wd,"Data/Intermediate"))[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(paste0(wd,"Data/Intermediate")))]
# load(paste0(wd,"Data/Intermediate/",model_files[length(model_files)]))
# ## this gives the posterior distributions for stan model for the trait, as an mcmc.list object

# ## design matrix for trait model
# X <- stan_data_ped_noRep$X


# pred_cond_pos <- grep((paste(c(cond_term),collapse="|")),colnames(X))
# pred_cont_pos <- grep((paste(c(cont_term),collapse="|")),colnames(X))
# # gets the position of the effects associated with conditioning and controlling variables in the design matrix for the trait models 

# # pred_cond<-colnames(X)[pred_cond_pos] # terms to condition on 
# # pred_cont<-colnames(X)[pred_cont_pos] # terms to control for 
# # gets the names of those

# pred_pos <- which(!grepl((paste(c(cond_term,cont_term),collapse="|")),colnames(X)))
# # gets the position of the effects that are neither controlling or conditioning (marginal variables)

# mu_pred_condX<-X[,pred_cond_pos,drop=FALSE]
# mu_pred_contX<-X[,pred_cont_pos,drop=FALSE]
# mu_predX<-X[,pred_pos,drop=FALSE]
# # gets the design matrices for the trait models for conditioning, controlling and marginal variables. 


# model_zpost<-matrix(posterior::as_draws_array(model_z), prod(dim(posterior::as_draws_array(model_z))[1:2]),dim(posterior::as_draws_array(model_z))[3], dimnames=list(NULL, dimnames(posterior::as_draws_array(model_z))[[3]]))[,-(1:7)]
# head(model_zpost)

# model_zbeta <- model_zpost[,grep("beta",colnames(model_zpost))]
# head(model_zbeta)
# # posterior for location effects from trait model

# model_znbeta<-model_zpost[,-grep("beta",colnames(model_zpost))]
# # posterior for distribution effects from trait model

# t_terms<-c("xi", "omega", "alpha", "nu")




# 	if(model_moments){
# 		mu_pred_cond<-mu_pred_condX%*%model_zbeta[i,pred_cond_pos]
# 		mu_pred_cont<-mu_pred_contX%*%model_zbeta[i,pred_cont_pos]
# 		mu_pred<-mu_predX%*%model_zbeta[i,pred_pos]
# 		# gets the predictions for the trait models for conditioning, controlling and marginal variables. 

# 		e_st<-list(n_st=model_znbeta[i,paste0(t_terms, "_nest")],
# 			       e_st=model_znbeta[i,paste0(t_terms, if(trait=="weight_g"){"_E"}else{"_ind"})], 
# 			       fixed_st=st.mple(y = mu_pred)$dp)

# 		if(trait!="weight_g") e_st$me_st<-c(0, model_znbeta[i,"sigma_E"], 0, 1e+16)

# 	    # list of environmental distribution parameters: xi, omega, alpha, nu
# 	    # fixed_st is a skew-t approximation for the marginal variable predictions 

# 		g_st<-c(0, model_znbeta[i,"sigma_A"], 0, 1e+16)
# 	# list of genetic distribution parameters
# 	}

#    		if(model_moments){
# 			mu<-dp2cm(c(list(g_st), e_st), family="ST")
# 			# central moments of the deviation of the (non-mean-centred) trait values from the conditional/controlling predictions as inferred from the trait model

# 			mu[1] <- mu[1] + mean(mu_pred_cond[THBW_noRep$category==cat])
# 			# mean adjusted for conditional prediction for category j 
#    		}else{
# 			mu <- scm(z_sub)
# 			# If model_moments=FALSE use the observed moments
#    		}