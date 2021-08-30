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
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced

source(paste0(wd,"R/00_functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))


###-----------------------------###
###-----All models--------------###
###-----------------------------###

for (trait in c("tarsus_mm","headbill_mm","weight_g","wing_mm")){
	## t hb and m run for 4000 warmup and 6000 after
	## w run for 5000 warmup and 10000 after
# trait="tarsus_mm"
	if(trait=="weight_g"){
		stan_data_trait <- stan_data_DS_noRep
		stanModel_DS <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_DamSire_delta.stan"))
	}else{
		stan_data_trait <- c(list(y=THBW[,trait]),stan_data_DS)
		stanModel_DS <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_ME_DamSire_delta.stan"))
	}

	asreml_mod <- if(trait=="weight_g"){modM}else if(trait=="tarsus_mm"){modT}else if(trait=="headbill_mm"){modHB}else if(trait=="wing_mm"){modW}

	mod_stan_DS <- stanModel_DS$sample(
	    data = stan_data_trait,
	    chains = 4,
		iter_sampling = 10000, 
		iter_warmup = 5000,
		thin= 40,
		refresh = 1000,
		max_treedepth = 12,
		init=function() makeSV(asreml_mod, stan_data_trait, ME=trait!="weight_g", delta=TRUE ,animal=FALSE)
	)

	model_z <- as.mcmc.list(lapply(mod_stan_DS$output_files(), function(x){
		x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
		x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|gamma|sigma|omega|xi|delta|alpha|^nu",names(x1)))])
		return(x2)
	}))

	save(model_z, file= paste0(wd,"Data/Intermediate/stanMod_DS_",if(reduced)"reduced_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

 	rm(stan_data_trait,stanModel_DS,asreml_mod,mod_stan_DS,model_z)

}










###-----------------------------
###  alpha normal(0,10) prior
###-----------------------------

rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

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
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced

source(paste0(wd,"R/00_functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))


for (trait in c("wing_mm")){

	if(trait=="weight_g"){
		stan_data_trait <- stan_data_DS_noRep
		stanModel_DS <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_DamSire_alpha10.stan"))
	}else{
		stan_data_trait <- c(list(y=THBW[,trait]),stan_data_DS)
		stanModel_DS <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_ME_DamSire_alpha10.stan"))
	}

	asreml_mod <- if(trait=="weight_g"){modM}else if(trait=="tarsus_mm"){modT}else if(trait=="headbill_mm"){modHB}else if(trait=="wing_mm"){modW}

	mod_stan_DS <- stanModel_DS$sample(
	    data = stan_data_trait,
	    chains = 4,
		iter_sampling = 10000, 
		iter_warmup = 5000,
		thin= 40,
		refresh = 1000,
		max_treedepth = 12,
		init=function() makeSV(asreml_mod, stan_data_trait, ME=trait!="weight_g", delta=FALSE ,animal=FALSE)
	)

	model_z <- as.mcmc.list(lapply(mod_stan_DS$output_files(), function(x){
		x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
		x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|gamma|sigma|omega|xi|delta|alpha|^nu",names(x1)))])
		return(x2)
	}))

	save(model_z, file= paste0(wd,"Data/Intermediate/stanMod_DS_",if(reduced)"reduced_",trait,"_alpha10_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

 	rm(stan_data_trait,stanModel_DS,asreml_mod,mod_stan_DS,model_z)

}









###-----------------------------
###  alpha normal(0,1) prior
###-----------------------------

rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

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
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced

source(paste0(wd,"R/00_functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))


for (trait in c("wing_mm")){
	if(trait=="weight_g"){
		stan_data_trait <- stan_data_DS_noRep
		stanModel_DS <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_DamSire_alpha10.stan"))
	}else{
		stan_data_trait <- c(list(y=THBW[,trait]),stan_data_DS)
		stanModel_DS <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_ME_DamSire_alpha1.stan"))
	}

	asreml_mod <- if(trait=="weight_g"){modM}else if(trait=="tarsus_mm"){modT}else if(trait=="headbill_mm"){modHB}else if(trait=="wing_mm"){modW}

	mod_stan_DS <- stanModel_DS$sample(
	    data = stan_data_trait,
	    chains = 4,
		iter_sampling = 10000, 
		iter_warmup = 5000,
		thin= 40,
		refresh = 1000,
		max_treedepth = 12,
		init=function() makeSV(asreml_mod, stan_data_trait, ME=trait!="weight_g", delta=FALSE ,animal=FALSE)
	)

	model_z <- as.mcmc.list(lapply(mod_stan_DS$output_files(), function(x){
		x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
		x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|gamma|sigma|omega|xi|delta|alpha|^nu",names(x1)))])
		return(x2)
	}))

	save(model_z, file= paste0(wd,"Data/Intermediate/stanMod_DS_",if(reduced)"reduced_",trait,"_alpha1_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

 	rm(stan_data_trait,stanModel_DS,asreml_mod,mod_stan_DS,model_z)

}
