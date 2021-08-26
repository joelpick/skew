rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

# library(MCMCglmm)
# library(MASS)
library(mvtnorm)
library(sn)
#library(rstan)
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

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced
cond_term<-analysis_options$cond_term  # terms to condition on 
cont_term<-analysis_options$cont_term  # terms to control for 
fixed_w <- analysis_options$fixed_w

ncores <- 8
n_it <- 1000
h2a_it<-10000              # number of simulated data to approximate h2a        
re_run_h2 <- TRUE
include_ME <- FALSE 		## should Vp include measurement error

load(paste0(data_wd,"chick_data.Rdata"))
# THBW is the full data
# THBW_noRep excluding repeat measures

load(paste0(data_wd,"stan_data",if(reduced)"_reduced",".Rdata"))
# data used in stan models, stan_data_ped_noRep$X is the design matrix (without repeat measures)

if(re_run_h2){
	for (trait in c("tarsus_mm","headbill_mm","wing_mm","weight_g")){
#trait="tarsus_mm"

		cat("Starting",trait,"... \n")

		model_files <- list.files(data_wd)[grep(paste0("stanModNormal2_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
		load(paste0(data_wd,model_files[length(model_files)]))
		model_zN <- model_z
		model_zpostN<-do.call(rbind,model_zN)[,-(1:7)]

		model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
		load(paste0(data_wd,model_files[length(model_files)]))
		## posterior distributions for stan model for the trait (model_z), as an mcmc.list object
		model_zpost<-do.call(rbind,model_z)[,-(1:7)]
		##combine into one matrix

		model_zbeta <- model_zpost[,grep("beta",colnames(model_zpost))]
		# posterior for location effects from trait model

		model_zbetaN <- model_zpost[,grep("beta",colnames(model_zpostN))]


		model_files <- list.files(data_wd)[grep(paste0("day15_survival_ME_",trait),list.files(data_wd))]
		load(paste0(data_wd,model_files[length(model_files)]))
		model_wpost<-do.call(rbind,model_w)[,-(1:7)]
		# read in fitness models (model_w)

		X<-stan_data_ped_noRep$X
		#design matrix for trait model

		pred_pos <- which(!grepl((paste(c(cond_term,cont_term),collapse="|")),colnames(X)))
		#gets the position of the effects that are neither controlling or conditioning (marginal variables)

		mu_pred_all <- apply(model_zbeta[,pred_pos,drop=FALSE], 1, function(x) X[,pred_pos,drop=FALSE]%*%x)
		mu_pred_allN <- apply(model_zbetaN[,pred_pos,drop=FALSE], 1, function(x) X[,pred_pos,drop=FALSE]%*%x)
		# gets the posterior of the predictions for the trait models for marginal variables. 


		pred_cond_pos <- grep((paste(c(cond_term),collapse="|")),colnames(X))
		mu_pred_cond_all <- apply(model_zbeta[,pred_cond_pos,drop=FALSE], 1, function(x) X[,pred_cond_pos,drop=FALSE]%*%x)
		## same for conditioning variables

		sigmas <- cbind(sigma_fix=apply(mu_pred_all,2,sd),model_zpost[,grep("sigma",colnames(model_zpost))])
		sigmasN <- cbind(sigma_fix=apply(mu_pred_allN,2,sd),model_zpostN[,grep("sigma",colnames(model_zpostN))])

# colMeans(sigmas)
# colMeans(sigmasN)

		## sds of random effects

		# heritability before selection
		h2b <- apply(sigmas[,paste0("sigma",c("_fix","_nest","_A",if(trait=="weight_g"){"_E"}else{"_ind"}))], 1, function(x) x["sigma_A"]^2/sum(x^2))
		Va <- sigmas[,"sigma_A"]^2

		h2bN <- apply(sigmasN[,paste0("sigma",c("_fix","_nest","_A",if(trait=="weight_g"){"_E"}else{"_ind"}))], 1, function(x) x["sigma_A"]^2/sum(x^2))
		h2bN_ME <- apply(sigmasN, 1, function(x) x["sigma_A"]^2/sum(x^2))


		z<-THBW_noRep[[trait]]

		if(trait=="weight_g"){
			zC<-THBW_noRep[,paste0(trait, "C")]
			# trait values in the survival model (note that it was globally mean centred)

			zmean <- mean(THBW_noRep[,trait])
			# value trait is centered at
		}else{
			zmean <- mean(model_wpost[,"beta0_x1"])
			zC<-THBW_noRep[,trait] - zmean
		}

		if(is.null(cond_term)){
			THBW_noRep$category <- as.factor(rep(1, nrow(THBW_noRep)))
		}else{
			THBW_noRep$category <- as.factor(apply(THBW_noRep[,cond_term], 1, paste, collapse=""))
		}

		n_comb<-nlevels(THBW_noRep$category)
		# add a category column to the data.frames used in the fitness and trait models (THBW_noRep)
		# that designate data into the (n_comb) categories based on the conditioning variables

		model_w_out<- extract_w2(model_wpost, trait, fixedEffects=formula(paste("~", fixed_w)), data=THBW_noRep)
		attach(model_w_out)
		## posteriors of parameters from the fitness model


		###################################
		##
		## Calculate heritability after selection
		##
		###################################


		h2a <- matrix(NA, n_it, n_comb)    # heritability after selection

		t_terms<-c("xi", "omega", "alpha", "nu")

		for(i in 1:n_it){
		#i=1
			if(i==1) time_start <- Sys.time() 

			e_st<-list(n_st=model_zpost[i,paste0(t_terms, "_nest")],
				       e_st=model_zpost[i,paste0(t_terms, if(trait=="weight_g"){"_E"}else{"_ind"})], 
				       fixed_st=st.mple(y = mu_pred_all[,i])$dp)
		    # list of environmental distribution parameters: xi, omega, alpha, nu
		    # fixed_st is a skew-t approximation for the marginal variable predictions 

			if(include_ME & trait!="weight_g"){
				e_st$me_st<-c(0, model_zpost[i,"sigma_E"], 0, 1e+16)
			}
		    # add on measurement error for traits that aren't weight 

			g_st<-c(0, model_zpost[i,"sigma_A"], 0, 1e+16)
			# list of genetic distribution parameters
		   
		    ##for each sex/year combination, work out the heritability after selection 
			h2a_sub <- mclapply(1:n_comb, function(j){
		      #j=2 	
		       	cat <- levels(THBW_noRep$category)[j]
		       	# name of the category defined by the conditioning variables 
		       	
		       	eta1_sub <- eta1_all[THBW_noRep$category==cat,i]
		       	eta2_sub <- eta2_all[THBW_noRep$category==cat,i]
				# z_sub <- z[THBW_noRep$category==cat]
		       	zC_sub <- zC[THBW_noRep$category==cat]
		       	# get (mean-centred) trait and linear predictors for the right category
		   
		   	    mu_etaz <- c(mean(eta1_sub), mean(eta2_sub), mean(zC_sub)) 
		   	    V_etaz <- cov(cbind(eta1_sub, eta2_sub, zC_sub))
 				# mean and covariance matrices of linear predictor and traits

				adj_mean <- mean(mu_pred_cond_all[THBW_noRep$category==cat])
 				# mean adjusted for conditional prediction for category j if that category was in the trait model 
				
				g<-rst(h2a_it, dp=g_st)
	  			zp<-g + rz(h2a_it, e_st) + adj_mean - zmean
	  			zo<-0.5*g+rst(h2a_it, dp=g_st)*sqrt(0.75) +rz(h2a_it, e_st) + adj_mean - zmean
	  			## phenotypes of parents and offspring

	  	## why are offspring given the year-sex mean?
		## should I centre on the global phenotypic mean, or on the predicted phenotypic mean from the model? (they are different, and will vary over iterations)

	  			Wz<-sapply(zp, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta_all[i,], gamma=gamma_all[i,],  V_nest=matrix(V_nest_all[i,], 2,2))
	  			## fitness of parents

		## this function then takes mu_etaz and V_etaz which are based on the observed data, does this matter? Because it includes a phenotypic mean for each year sex combo, which will differ from the model prediction

	  			wz<-Wz/mean(Wz)
	  			## relative fitness of parents

      			h2<-2*(mean(wz*zp*zo)-mean(wz*zp)*mean(wz*zo))/mean((wz*zp^2)-mean(wz*zp)^2)

		   		return(h2)
		   	}, mc.cores = ncores)
			
			h2a[i,]<-do.call(c,h2a_sub)
			
			if(((i/n_it)*100) %in% (c(1:10)*10)) {
				cat((i/n_it)*100,"% ", difftime(Sys.time(),time_start, units="mins")," mins \n",sep="")
				time_start <- Sys.time()		
			}
		}
		save(h2a, h2b, Va, h2bN, h2bN_ME, file=paste0(data_wd,"h2_",trait,".Rdata"))
		detach(model_w_out)

	}
}


# load(paste0(data_wd,"h2_tarsus_mm.Rdata"))
# h2a
	