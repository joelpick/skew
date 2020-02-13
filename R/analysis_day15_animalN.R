rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(MCMCglmm)
library(MasterBayes)


factorisePed <- function(pedigree, unknown=0){
    
    new_ped <- data.frame(
        1:nrow(pedigree), 
        ifelse(is.na(pedigree[,2]),unknown,match(pedigree[,2], pedigree[,1])), 
        ifelse(is.na(pedigree[,3]),unknown,match(pedigree[,3], pedigree[,1]))
        )
    colnames(new_ped) <- colnames(pedigree)[1:3]

    return(new_ped)
}

##UNIVARITE

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

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


X <- model.matrix(~ male_present + clutch_sizeC + nest_hatch_dateC  + hatch_day + year + timeC + sex + egg_weightC,THBW_egg)
# X <- model.matrix(~ 1,THBW_egg)
#any(!unique(THBW_tarsus$bird_id)%in%tPED$bird_id)

## get pruned pedigree, with dummy parents for birds with a single known parent


ped <- insertPed(prunePed(orderPed(tPED_use[1:3]), unique(THBW_egg$bird_id), make.base = TRUE))
names(ped) <- c("animal","dam", "sire")
MSV <- inverseA(ped)$dii

stan_ped <- factorisePed(ped[,1:3])

NoParents <- which(stan_ped[,2]==0 & stan_ped[,3]==0)
ParentsOffspring <- which(stan_ped[,1] %in% c(stan_ped[,2],stan_ped[,3]) & !stan_ped[,1] %in% NoParents)
ParentsNoOffspring <- which(!stan_ped[,1] %in% c(NoParents,ParentsOffspring))
length(NoParents); length(ParentsOffspring); length(ParentsNoOffspring)

## add first row that is base population
stan_ped2 <- rbind(c(0,-1,-1), stan_ped)+1

animal_id <- stan_ped[match(as.character(THBW_egg$bird_id), ped$animal),"animal"]

stan_data_THBW <- list(
	N=nrow(THBW_egg),
	J=ncol(X), 
	X=X,
	N_nest=length(unique(THBW_egg$nest)),
	nest_id=as.numeric(as.factor(THBW_egg$nest)), 
	N_ind=length(unique(THBW_egg$bird_id)),
	ind_id=as.numeric(as.factor(THBW_egg$bird_id)), 
	N_ped = nrow(stan_ped2), 
	animal_id =animal_id+1, 
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
# X_noRep <- model.matrix(~ 1,THBW_egg_noRep)

animal_id_noRep <- stan_ped[match(as.character(THBW_egg_noRep$bird_id), ped$animal),"animal"]


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


save(THBW_egg,stan_data_THBW,stan_data_weight, file= paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))

###-----------------------------###
###-----Tarsus------------------###
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

load(paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))

stan_data_tarsus <- c(list(y=THBW_egg$tarsus_mm),stan_data_THBW)

stanModel_pedN_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_pedN.stan"))

mod_stan_tarsus_pedN <- sampling(
	stanModel_pedN_ME, 
	data = stan_data_tarsus, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"sigma_A",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000
	#control=list(adapt_delta=0.9)#max_treedepth=12, 
	)
save(mod_stan_tarsus_pedN, file= paste0(wd,"Data/Intermediate/stanModTarsus_pedN",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))


###-----------------------------###
###-----headbill------------------###
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

load(paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))

stan_data_HB <- c(list(y=THBW_egg$headbill_mm),stan_data_THBW)

stanModel_pedN_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_pedN.stan"))

mod_stan_HB_pedN <- sampling(
	stanModel_pedN_ME, stan_data_HB, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"sigma_A",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000
	#control=list(adapt_delta=0.9)#max_treedepth=12, 
	)
save(mod_stan_HB_pedN, file= paste0(wd,"Data/Intermediate/stanModHB_pedN",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))


###-----------------------------###
###-----wing------------------###
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

load(paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))

stan_data_wing <- c(list(y=THBW_egg$wing_mm),stan_data_THBW)

stanModel_pedN_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_pedN.stan"))

mod_stan_wing_pedN <- sampling(
	stanModel_pedN_ME, stan_data_wing, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"sigma_A",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000
	#control=list(adapt_delta=0.9)#max_treedepth=12, 
	)
save(mod_stan_wing_pedN, file= paste0(wd,"Data/Intermediate/stanModWing_pedN",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))


###-----------------------------###
###-----Mass------------------###
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

load(paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))



stanModel_pedN <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_pedN.stan"))

mod_stan_weight_pedN <- sampling(
	stanModel_pedN, 
	data = stan_data_weight, 
	pars =c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"sigma_A",
		"xi_E", "omega_E", "alpha_E", "nu_E", "sigma_E", "gamma_E"
		 ),
	chains=4, 
	iter = 10000, 
	warmup = 4000)
save(mod_stan_weight_pedN, file= paste0(wd,"Data/Intermediate/stanModWeightPedN",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
#load(paste0(wd,"Data/Intermediate/stanModWeightPedN20190817_0218.Rdata"))
