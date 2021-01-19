rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

# library(MCMCglmm)
# library(MASS)
# library(mvtnorm)
library(sn)
library(rstan)
library(cmdstanr)
library(coda)
library(cubature)
library(parallel)


if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/functions.R"))

ncores <- 8
reduced <- TRUE
cond_term<-c("year", "sex") # terms to condition on 
cont_term<-c("timeC")   	# terms to control for
re_run_POreg<-TRUE
re_run_h2 <- TRUE
n_it <- 1000
h2a_it<-10000              # number of simulated data to approximate h2a        





load(paste0(data_wd,"chick_data.Rdata"))
# THBW is the full data
# THBW_noRep excluding repeat measures

load(paste0(data_wd,"stan_data",if(reduced)"_reduced",".Rdata"))
# data used in stan models, stan_data_ped_noRep$X is the design matrix (without repeat measures)

# load(paste0(wd,"Data/Intermediate/Wplot_points.Rdata"))
# # trait values at which to evaluate the fitness function for visualisation purposes
# Wplot_points<-get(paste0("Wplot_points_",trait))
# nplot_points<-length(Wplot_points)
if(re_run_h2){
	for (trait in c("tarsus_mm","headbill_mm","wing_mm","weight_g")){
		model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
		load(paste0(data_wd,model_files[length(model_files)]))
		## posterior distributions for stan model for the trait (model_z), as an mcmc.list object

		model_zpost<-matrix(posterior::as_draws_array(model_z), prod(dim(posterior::as_draws_array(model_z))[1:2]),dim(posterior::as_draws_array(model_z))[3], dimnames=list(NULL, dimnames(posterior::as_draws_array(model_z))[[3]]))[,-(1:7)]
		##combine into one matrix

		model_zbeta <- model_zpost[,grep("beta",colnames(model_zpost))]
		# posterior for location effects from trait model

		model_znbeta<-model_zpost[,-grep("beta",colnames(model_zpost))]
		# posterior for distribution effects from trait model

		load(paste0(data_wd,"day15_survival_model_",trait,".Rdata"))
		# read in fitness models (model_w)

		X<-stan_data_ped_noRep$X
		#design matrix for trait model

		pred_pos <- which(!grepl((paste(c(cond_term,cont_term),collapse="|")),colnames(X)))
		#gets the position of the effects that are neither controlling or conditioning (marginal variables)

		mu_pred_all <- apply(model_zbeta[,pred_pos,drop=FALSE], 1, function(x) X[,pred_pos,drop=FALSE]%*%x)
		# gets the posterior of the predictions for the trait models for marginal variables. 


		t_terms<-c("xi", "omega", "alpha", "nu")


		# pred_cond_pos <- grep((paste(c(cond_term),collapse="|")),colnames(X))
		# pred_cont_pos <- grep((paste(c(cont_term),collapse="|")),colnames(X))
		# pred_pos <- setdiff(1:ncol(X),c(pred_cond_pos,pred_cont_pos))

		# mu_pred_condX <- X[,pred_cond_pos,drop=FALSE]
		# mu_pred_contX <- X[,pred_cont_pos,drop=FALSE]
		# mu_predX <- X[,pred_pos,drop=FALSE]

		# mu_pred_cond_all <- apply(model_zbeta[,pred_cond_pos,drop=FALSE], 1, function(x) mu_pred_condX%*%x)
		# mu_pred_cont_all <- apply(model_zbeta[,pred_cont_pos,drop=FALSE], 1, function(x) mu_pred_contX%*%x)
		# mu_pred_all <- apply(model_zbeta[,pred_pos,drop=FALSE], 1, function(x) mu_predX%*%x)





		z<-THBW_noRep[[trait]]
		zC<-THBW_noRep[,paste0(trait, "C")]
		# trait values in the survival model (note that it was globally mean centred)
		zmean <- mean(THBW_noRep[,trait])
		# value trait is centered at


		if(is.null(cond_term)){
			THBW_noRep$category <- as.factor(rep(1, nrow(THBW_noRep)))
		}else{
			THBW_noRep$category <- as.factor(apply(THBW_noRep[,cond_term], 1, paste, collapse=""))
		}

		n_comb<-nlevels(THBW_noRep$category)
		# add a category column to the data.frames used in the fitness and trait models (THBW_noRep)
		# that designate data into the (n_comb) categories based on the conditioning variables

		attach(extract_w(model_w, trait, fixedEffects=formula(~ sex + male_present + year+hatch_day + clutch_sizeC + nest_hatch_dateC), data=THBW_noRep))



		###################################
		##
		## Calculate heritability before and after selection
		##
		###################################


		h2a <- matrix(NA, n_it, n_comb)    # heritability after selection
		h2b <- matrix(NA, n_it, 1)         # heritability before selection
		Va <-  matrix(NA, n_it, 1) 

		for(i in 1:n_it){
		#i=1
			if(i==1) time_start <- Sys.time() 

			e_st<-list(n_st=model_znbeta[i,paste0(t_terms, "_nest")],
				       e_st=model_znbeta[i,paste0(t_terms, if(trait=="weight_g"){"_E"}else{"_ind"})], 
				       fixed_st=st.mple(y = mu_pred_all[,i])$dp)

			if(trait!="weight_g") e_st$me_st<-c(0, model_znbeta[i,"sigma_E"], 0, 1e+16)
		    # list of environmental distribution parameters: xi, omega, alpha, nu
		    # fixed_st is a skew-t approximation for the marginal variable predictions 

			g_st<-c(0, model_znbeta[i,"sigma_A"], 0, 1e+16)
			# list of genetic distribution parameters

			mu <- dp2cm(c(list(g_st), e_st), family="ST")
		    Va[i] <- model_znbeta[i,"sigma_A"]^2
		    h2b[i] <- Va[i]/mu[2]

			beta <- beta_all[i,]
		    gamma <- gamma_all[i,]
		    # trait coefficients from survival model

		    V_nest<-matrix(V_nest_all[i,], 2,2)
		   
			h2a_sub <- mclapply(1:n_comb, function(j){
		       	cat <- levels(THBW_noRep$category)[j]
		       	# name of the category defined by the conditioning variables 
		       	
		       	eta1_sub <- eta1_all[THBW_noRep$category==cat,i]
		       	eta2_sub <- eta2_all[THBW_noRep$category==cat,i]
		       	z_sub <- zC[THBW_noRep$category==cat]
		       	# get (mean-centred) trait and linear predictors for the right category
		   
		   	    mu_etaz <- c(mean(eta1_sub), mean(eta2_sub), mean(z_sub)) 
		   	    V_etaz <- cov(cbind(eta1_sub, eta2_sub, z_sub))
		   
				# if(model_moments){
				# 	mu<-dp2cm(c(list(g_st), e_st), family="ST")
				# 	# central moments of the deviation of the (non-mean-centred) trait values from the conditional/controlling predictions as inferred from the trait model

				# 	mu[1] <- mu[1] + mean(mu_pred_cond[THBW_noRep$category==cat])
				# 	# mean adjusted for conditional prediction for category j 
		  #  		}else{
				# 	mu <- scm(z_sub)
				# 	# If model_moments=FALSE use the observed moments
		  #  		}
				mu <- scm(z_sub)
		   		return(h2(h2a_it, g_st=g_st, e_st=e_st, adj_mean=mu[1], mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest))
		   	}, mc.cores = ncores)
			
			h2a[i,]<-do.call(c,h2a_sub)
			
			if(((i/n_it)*100) %in% (c(1:10)*10)) {
				cat((i/n_it)*100,"% ", difftime(Sys.time(),time_start, units="mins")," mins \n",sep="")
				time_start <- Sys.time()		
			}


		}

		save(h2a, h2b, Va, file=paste0(data_wd,"h2_",trait,".Rdata"))
	}
}

