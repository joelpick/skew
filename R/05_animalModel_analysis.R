rm(list=ls())

options(stringsAsFactors=FALSE)

# library(rstan)
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(cmdstanr)
library(coda)
library(mvtnorm)

##UNIVARITE

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/github/skew/"
}

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced

source(paste0(wd,"R/00_functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))


#names(makeSV(modT, stan_data_ped_noRep, ME=TRUE, delta=TRUE))

###-----------------------------###
###-----All models--------------###
###-----------------------------###

for (trait in c("tarsus_mm","headbill_mm","wing_mm","weight_g")){
# trait="tarsus_mm"
	if(trait=="weight_g"){
		stan_data_trait <- stan_data_ped_noRep
		stanModel_pedN <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_pedN_delta.stan"))
	}else{
		stan_data_trait <- c(list(y=THBW[,trait]),stan_data_ped)
		stanModel_pedN <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_ME_pedN_delta.stan"))
	}

	asreml_mod <- if(trait=="weight_g"){modM}else if(trait=="tarsus_mm"){modT}else if(trait=="headbill_mm"){modHB}else if(trait=="wing_mm"){modW}

	mod_stan_pedN <- stanModel_pedN$sample(
	    data = stan_data_trait,
	    chains = 4,
		iter_sampling = 6000, 
		iter_warmup = 4000,
		thin= 24,
		refresh = 1000,
		max_treedepth = 12,
		init=function() makeSV(asreml_mod, stan_data_trait, ME=trait!="weight_g", delta=TRUE )
	)

	model_z <- as.mcmc.list(lapply(mod_stan_pedN$output_files(), function(x){
		x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
		x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|gamma|sigma|omega|xi|delta|alpha|^nu",names(x1)))])
		return(x2)
	}))

	save(model_z, file= paste0(wd,"Data/Intermediate/stanMod_pedN_",if(reduced)"reduced_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

 	rm(stan_data_trait,stanModel_pedN,asreml_mod,mod_stan_pedN,model_z)

}






###-----------------------------###
###---Gaussian animal models----###
###-----------------------------###


for (trait in c("tarsus_mm","headbill_mm","wing_mm","weight_g")){
# trait="tarsus_mm"
	if(trait=="weight_g"){
		stan_data_trait <- stan_data_ped_reduced_noRep
		stanModel_pedN <- cmdstan_model(paste0(wd,"stan/normal_RE_day15_ped2.stan"))
	}else{
		stan_data_trait <- c(list(y=THBW[,trait]),stan_data_ped_reduced)
		stanModel_pedN <- cmdstan_model(paste0(wd,"stan/normal_RE_day15_ME_ped2.stan"))
	}

	# asreml_mod <- if(trait=="weight_g"){modM}else if(trait=="tarsus_mm"){modT}else if(trait=="headbill_mm"){modHB}else if(trait=="wing_mm"){modW}

	mod_stan_pedN <- stanModel_pedN$sample(
	    data = stan_data_trait,
	    chains = 4,
		iter_sampling = 6000, 
		iter_warmup = 4000,
		thin= 24,
		refresh = 1000,
		max_treedepth = 10
		# , init=function() makeSV(asreml_mod, stan_data_trait, ME=trait!="weight_g", delta=TRUE )
	)

	model_z <- as.mcmc.list(lapply(mod_stan_pedN$output_files(), function(x){
		x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
		x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|gamma|sigma|omega|xi|delta|alpha|^nu",names(x1)))])
		return(x2)
	}))

	save(model_z, file= paste0(wd,"Data/Intermediate/stanModNormal2_pedN_",if(reduced)"reduced_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

 	rm(stan_data_trait,stanModel_pedN,mod_stan_pedN,model_z)

}