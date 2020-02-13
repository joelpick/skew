rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library("bayesplot")

dtPlot <- function(model, pars=character(), transform=list() ){
	available_mcmc(pattern = "_nuts_")
	color_scheme_set("darkgray")
	posterior <- as.array(model)
	np <- nuts_params(model)
	mcmc_parcoord(posterior ,transform = transform ,np = np, pars=pars)
# transform=function(x) {(x - mean(x)) / sd(x)}
}

##UNIVARITE

div_trans <- function(model){
	np <- nuts_params(model)
	div<-subset(np,Parameter=="divergent__")
	with(div, tapply(Value,Chain,sum))
}

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

load(file= paste0(wd,"Data/Intermediate/dam_sire_egg.Rdata"))


THBW_egg<- subset(THBW_egg,!is.na(sex) & !year%in%c(2019)& !is.na(tarsus_mm)& !is.na(headbill_mm)& !is.na(wing_mm)& !is.na(weight_g))
# assign individuals with hatch day 6, hatch day 3
THBW_egg[THBW_egg$hatch_day==6,"hatch_day"]<-3

THBW_egg$male_present <- as.factor(THBW_egg$male_present)
THBW_egg$year <- as.factor(THBW_egg$year)
THBW_egg$sex <- as.factor(THBW_egg$sex)
THBW_egg$hatch_day <- as.factor(THBW_egg$hatch_day)
THBW_egg$timeC <- scale(THBW_egg$time_hr, scale=FALSE)
THBW_egg$clutch_sizeC <- scale(THBW_egg$clutch_size, scale=FALSE)
THBW_egg$nest_hatch_dateC <- scale(THBW_egg$nest_hatch_date, scale=FALSE)
THBW_egg$egg_weightC <- scale(THBW_egg$egg_weight, scale=FALSE)

# THBW_egg$sireDummy <- THBW_egg$sire
# missing_sire <- any(is.na(THBW_egg$sireDummy))
# THBW_egg$sireDummy[missing_sire] <- paste0("dummy",1:sum(missing_sire))

dam_sire <- unique(c(THBW_egg$dam_P,THBW_egg$sire_P))
THBW_egg$dam <- as.numeric(factor(THBW_egg$dam_P, levels=dam_sire))
THBW_egg$sire <- as.numeric(factor(THBW_egg$sire_P, levels=dam_sire))

X <- model.matrix(~ male_present + clutch_sizeC + nest_hatch_dateC + hatch_day + year + timeC + sex + egg_weightC,THBW_egg)


stan_data_THBW <- list(
	N=nrow(THBW_egg),
	J=ncol(X), 
	N_nest=length(unique(THBW_egg$nest)),
	N_ind=length(unique(THBW_egg$bird_id)),
	# N_dam=length(unique(THBW_egg$dam_P)),
	N_dam_sire=length(dam_sire),
	X=X,
	nest_id=as.numeric(as.factor(THBW_egg$nest)), 
	ind_id=as.numeric(as.factor(THBW_egg$bird_id)), 
	#dam_id=as.numeric(as.factor(THBW_egg$dam_P))
	dam_id=THBW_egg$dam,
	sire_id=THBW_egg$sire
) 


stanModel_dam_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_DamSire2.stan"))
#stanModel_damN_ME <- stan_model(file = paste0(wd,"stan/normal_RE_day15_ME_DamSire2.stan"))


stan_data_weight <- list(
	N=nrow(THBW_egg_noRep),
	y=THBW_egg_noRep$weight_g,
	J=ncol(X_noRep ), 
	X=X_noRep ,
	N_nest=length(unique(THBW_egg_noRep$nest)),
	nest_id=as.numeric(as.factor(THBW_egg_noRep$nest)), 
	N_ind=length(unique(THBW_egg_noRep$bird_id)),
	ind_id=as.numeric(as.factor(THBW_egg_noRep$bird_id)), 
	N_ped = nrow(stan_ped2), 
	animal_id =animal_id_noRep+1, 
	dam = stan_ped2$dam,
	sire = stan_ped2$sire,
	MSV=c(1,MSV), 
	N_NoParents=length(NoParents), 
	N_ParentsOffspring=length(ParentsOffspring), 
	N_ParentsNoOffspring=length(ParentsNoOffspring), 
	NoParents=NoParents+1, 
	ParentsOffspring=ParentsOffspring+1,
	ParentsNoOffspring=ParentsNoOffspring+1
	)

#######-----------------------------------

THBW_egg_noRep<- subset(THBW_egg,!duplicated(bird_id))

X_noRep <- model.matrix(~ male_present + clutch_sizeC + nest_hatch_dateC  + hatch_day + year + timeC + sex + egg_weightC,THBW_egg_noRep)


stan_data_weight <- list(
	N=nrow(THBW_egg_noRep),
	y=THBW_egg_noRep$weight_g,
	J=ncol(X_noRep), 
	N_nest=length(unique(THBW_egg_noRep$nest)),
	N_ind=length(unique(THBW_egg_noRep$bird_id)),
	# N_dam=length(unique(THBW_egg_noRep$dam_P)),
	N_dam_sire=length(dam_sire),
	X=X_noRep,
	nest_id=as.numeric(as.factor(THBW_egg_noRep$nest)), 
	ind_id=as.numeric(as.factor(THBW_egg_noRep$bird_id)), 
	#dam_id=as.numeric(as.factor(THBW_egg_noRep$dam_P))
	dam_id=THBW_egg_noRep$dam,
	sire_id=THBW_egg_noRep$sire
) 

#######-----------------------------------


save(THBW_egg,THBW_egg_noRep,stan_data_THBW,stan_data_weight, file= paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))




###-----------------------------###
###-----Tarsus------------------###
###-----------------------------###



rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library("bayesplot")

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))

stanModel_dam_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_DamSire2.stan"))

stan_data_tarsus <- c(list(y=THBW_egg$tarsus_mm),stan_data_THBW)


load(paste0(wd,"Data/Intermediate/dam_sire_egg.Rdata"))


library(mvtnorm)
animal2DS <- function(x){
	out <- x[!names(x)%in%c("animal","R","nest_orig")]
	out["dam"] <- x["animal"]/4
	out["R"] <- x["R"] + out["dam"]*2
	return(out)
}

## function to make starting values from asreml models
makeSV <- function(asreml_data,stan_data) {
	sigma_vars <- colnames(asreml_data$mod_egg_DS_Ri)
	sigma_starting <- sqrt(rmvnorm(n = 1,
		mean = animal2DS(asreml_data$vars_mean["mod_egg_DS",])[sigma_vars],
		sigma = asreml_data$mod_egg_DS_Ri
		))
	beta_starting <- rmvnorm(n = 1,
		mean = asreml_data$fixed_mean["mod_egg_DS",],
		sigma = asreml_data$mod_egg_DS_Fi
		)
	## beta_tilde from QR decomposition
	beta_tilde_starting <- as.vector((qr.R(qr(stan_data$X))*1/sqrt(nrow(stan_data$X)-1)) %*% t(beta_starting))

	SV <- list(
		beta_tilde = beta_tilde_starting,

		sigma_nest = sigma_starting[,"nest"],
		alpha_nest = abs(rnorm(1,0,0.1))*-1,
		nu_nest = 30,
		nest_effects = rnorm(stan_data$N_nest,0,sigma_starting[,"nest"]),
		
		sigma_dam_sire=sigma_starting[,"dam"],
		alpha_dam_sire = abs(rnorm(1,0,0.1))*-1,
		nu_dam_sire = 30,
		dam_sire_effects = rnorm(stan_data$N_dam_sire,0,sigma_starting[,"dam"]),
		
		sigma_ind = sigma_starting[,"bird_id"],
		alpha_ind = abs(rnorm(1,0,0.1))*-1,
		nu_ind = 30,
		ind_effects = rnorm(stan_data$N_ind,0,sigma_starting[,"bird_id"]),

		sigma_E = sigma_starting[,"R"]
	)

	return(SV)
}


### scale everything so priors on fixed effects do better

n_chains=1
mod_stan_tarsus_dam <- sampling(
	stanModel_dam_ME, 
	data = stan_data_tarsus, 
	# pars =c("beta", 
	# 	"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
	# 	"xi_dam_sire", "omega_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_dam_sire", "gamma_dam_sire",
	# 	#"xi_dam", "omega_dam", "alpha_dam", "nu_dam", "sigma_dam", "gamma_dam",
	# 	"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
	# 	"sigma_E"
	# 	 ),
	chains=n_chains, 
	iter = 10000, 
	warmup = 4000,
	thin=3,
	, control=list(adapt_delta=0.97)#max_treedepth=12, 
	, init= lapply(1:n_chains,function(x) makeSV(modT, stan_data_tarsus))
	)
save(mod_stan_tarsus_dam, file= paste0(wd,"Data/Intermediate/stanModTarsus_DS",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
## 4 hours with default control settings
## nearly all transitions are divergent 




# mod_stan_tarsus_damN <- sampling(
# 	stanModel_damN_ME, 
# 	data = stan_data_tarsus, 
# 	pars =c("beta", "sigma_nest", "sigma_dam_sire", "sigma_ind","sigma_E"),
# 	chains=4, 
# 	iter = 5000, 
# 	warmup = 2000
# 	# , control=list(adapt_delta=0.95)#max_treedepth=12, 
# 	)
# save(mod_stan_tarsus_damN, file= paste0(wd,"Data/Intermediate/stanModTarsus_DS_N",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))


#load(paste0(wd,"Data/Intermediate/stanModTarsus_DS20191213_0005.Rdata"))
# summary(mod_stan_tarsus_dam)$summary[,c(4,1,8,9,10)]
# div_trans(mod_stan_tarsus_dam)
# pairs(mod_stan_tarsus_dam,pars =c("sigma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_ind", "alpha_ind", "nu_ind", "sigma_E","lp__"))
# rstan::traceplot(mod_stan_tarsus_dam, pars =c("sigma_nest", "gamma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "gamma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_ind", "gamma_ind", "alpha_ind", "nu_ind", "sigma_E","lp__"))
# plot(as.data.frame(mod_stan_tarsus_dam)[,"sigma_nest"],log(as.data.frame(mod_stan_tarsus_dam)[,"nu_nest"]))

# dtPlot(mod_stan_tarsus_dam,transform = function(x) {(x - mean(x)) / sd(x)})#,pars =c("sigma_nest", "alpha_nest", "nu_nest", "sigma_dam", "alpha_dam", "nu_dam", "sigma_ind", "alpha_ind", "nu_ind", "sigma_E"))



###-----------------------------###
###-----Head-Bill---------------###
###-----------------------------###


rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))

stanModel_dam_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_DamSire2.stan"))


stan_data_HB <- c(list(y=THBW_egg$headbill_mm),stan_data_THBW)

mod_stan_HB_dam <- sampling(
	stanModel_dam_ME, 
	data = stan_data_HB, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"xi_dam_sire", "omega_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_dam_sire", "gamma_dam_sire",
		#"xi_dam", "omega_dam", "alpha_dam", "nu_dam", "sigma_dam", "gamma_dam",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000
	#, control=list(adapt_delta=0.9)#max_treedepth=12, 
	)
save(mod_stan_HB_dam, file= paste0(wd,"Data/Intermediate/stanModHB_DS",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))


###-----------------------------###
###-----Wing--------------------###
###-----------------------------###
rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))

stanModel_dam_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_DamSire2.stan"))

stan_data_wing <- c(list(y=THBW_egg$wing_mm),stan_data_THBW)

mod_stan_wing_dam <- sampling(
	stanModel_dam_ME, 
	data = stan_data_wing, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"xi_dam_sire", "omega_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_dam_sire", "gamma_dam_sire",
		#"xi_dam", "omega_dam", "alpha_dam", "nu_dam", "sigma_dam", "gamma_dam",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000
	#, control=list(adapt_delta=0.9)#max_treedepth=12, 
	)
save(mod_stan_wing_dam, file= paste0(wd,"Data/Intermediate/stanModWing_DS",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))



###-----------------------------###
###-----Weight------------------###
###-----------------------------###


rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))


stanModel_dam <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_DamSire2.stan"))

mod_stan_weight_dam <- sampling(
	stanModel_dam, 
	data = stan_data_weight, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"xi_dam_sire", "omega_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_dam_sire", "gamma_dam_sire",
		"xi_E", "omega_E", "alpha_E", "nu_E", "sigma_E", "gamma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000)
save(X_noRep,mod_stan_weight_dam, file= paste0(wd,"Data/Intermediate/stanModWeight_DS",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))



#load(paste0(wd,"Data/Intermediate/stanModWeightDam20190816_2150.Rdata"))

# summary(mod_stan_weight_dam)$summary[,c(4,1,8,9,10)]
# div_trans(mod_stan_weight_dam)
# pairs(mod_stan_weight_dam,pars =c("sigma_nest", "gamma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "gamma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_E", "gamma_E", "alpha_E", "nu_E","lp__"))
# rstan::traceplot(mod_stan_weight_dam, pars =c("sigma_nest", "gamma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "gamma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_E", "gamma_E", "alpha_E", "nu_E","lp__"))
