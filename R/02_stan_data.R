rm(list=ls())

options(stringsAsFactors=FALSE)

#library(rstan)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

library(MasterBayes)
library(MCMCglmm)

##UNIVARITE

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/00_functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

## make stan data with both reduced and full fixed effect structures 
for(reduced in c(TRUE,FALSE)){
	
	#######
	## make design matrix
	#######

	if(reduced){
		X <- model.matrix(analysis_options$fixed_z_reduced,THBW)

		X_noRep <- model.matrix(analysis_options$fixed_z_reduced,THBW_noRep)

	}else{
		X <- model.matrix(analysis_options$fixed_z,THBW)

		X_noRep <- model.matrix(analysis_options$fixed_z,THBW_noRep)
	}

	#######
	## Data list for stan for Dam sire models 
	#######

	stan_data_DS <- list(
		N=nrow(THBW),
		J=ncol(X), 
		N_nest=length(unique(THBW$nest)),
		N_ind=length(unique(THBW$bird_id)),
		N_dam_sire=max(THBW$sire),
		X=X,
		nest_id=as.numeric(as.factor(THBW$nest)), 
		ind_id=as.numeric(as.factor(THBW$bird_id)), 
		dam_id=THBW$dam,
		sire_id=THBW$sire
	) 

	#######
	## Data list for stan for Dam sire models, with no repeat measures (i.e. for weight models)
	#######

	stan_data_DS_noRep <- list(
		N=nrow(THBW_noRep),
		y=THBW_noRep$weight_g,
		J=ncol(X_noRep), 
		N_nest=length(unique(THBW_noRep$nest)),
		N_ind=length(unique(THBW_noRep$bird_id)),
		N_dam_sire=max(THBW$sire),
		X=X_noRep,
		nest_id=as.numeric(as.factor(THBW_noRep$nest)), 
		ind_id=as.numeric(as.factor(THBW_noRep$bird_id)), 
		dam_id=THBW_noRep$dam,
		sire_id=THBW_noRep$sire
	) 



	#######
	## Prepare pedigree data
	#######

	## Mendelian sampling variance
	MSV <- inverseA(ped)$dii

	stan_ped <- factorisePed(ped[,1:3])

	NoParents <- which(stan_ped[,2]==0 & stan_ped[,3]==0)
	ParentsOffspring <- which(stan_ped[,1] %in% c(stan_ped[,2],stan_ped[,3]) & !stan_ped[,1] %in% NoParents)
	ParentsNoOffspring <- which(!stan_ped[,1] %in% c(NoParents,ParentsOffspring))
	length(NoParents); length(ParentsOffspring); length(ParentsNoOffspring)

	## add first row that is base population
	stan_ped2 <- rbind(c(0,-1,-1), stan_ped)+1


	animal_id <- stan_ped[match(as.character(THBW$bird_id), ped$animal),"animal"]

	animal_id_noRep <- stan_ped[match(as.character(THBW_noRep$bird_id), ped$animal),"animal"]


	#######
	## Data list for stan for animal models
	#######

	stan_data_ped <- list(
		N=nrow(THBW),
		J=ncol(X), 
		X=X,
		N_nest=length(unique(THBW$nest)),
		nest_id=as.numeric(as.factor(THBW$nest)), 
		N_ind=length(unique(THBW$bird_id)),
		ind_id=as.numeric(as.factor(THBW$bird_id)), 
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

	#######
	## Data list for stan for animal models, with no repeat measures (i.e. for weight models)
	#######

	stan_data_ped_noRep <- list(
		N=nrow(THBW_noRep),
		y=THBW_noRep$weight_g,
		J=ncol(X_noRep ), 
		X=X_noRep ,
		N_nest=length(unique(THBW_noRep$nest)),
		nest_id=as.numeric(as.factor(THBW_noRep$nest)), 
		N_ind=length(unique(THBW_noRep$bird_id)),
		ind_id=as.numeric(as.factor(THBW_noRep$bird_id)), 
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




	NoParentsOffspring <- which(stan_ped[,2]==0 & stan_ped[,3]==0 & stan_ped[,1] %in% c(stan_ped[,2],stan_ped[,3]))

	ParentsOffspring <- which( (stan_ped[,2]!=0 | stan_ped[,3]!=0) & stan_ped[,1] %in% c(stan_ped[,2],stan_ped[,3]))# & !stan_ped[,1] %in% NoParents
	#any(ParentsOffspring %in% NoParentsOffspring)

	NoParentsNoOffspring <- which(stan_ped[,2]==0 & stan_ped[,3]==0 & !stan_ped[,1] %in% c(stan_ped[,2],stan_ped[,3]))

	#any(NoParentsNoOffspring  %in% c(ParentsOffspring, NoParentsOffspring))
	ParentsNoOffspring <- which((stan_ped[,2]!=0 | stan_ped[,3]!=0) & !stan_ped[,1] %in% c(stan_ped[,2],stan_ped[,3]))


	stan_data_ped_reduced <- list( 
		N=nrow(THBW),
		J=ncol(X), 
		X=X,
		N_nest=length(unique(THBW$nest)),
		nest_id=as.numeric(as.factor(THBW$nest)), 
		N_ind=length(unique(THBW$bird_id)),
		ind_id=as.numeric(as.factor(THBW$bird_id)), 
		N_ped = nrow(stan_ped), 
		animal_id =animal_id, 
		dam = stan_ped$dam,
		sire = stan_ped$sire,
		MSV=MSV, 
		N_NoParentsOffspring=length(NoParentsOffspring), 
		N_ParentsOffspring=length(ParentsOffspring), 
		N_NoOffspring=length(c(NoParentsNoOffspring,ParentsNoOffspring)), 
		NoParentsOffspring=NoParentsOffspring, 
		ParentsOffspring=ParentsOffspring,
		NoOffspring=c(NoParentsNoOffspring,ParentsNoOffspring)
	)


	stan_data_ped_reduced_noRep <- list(
		N=nrow(THBW_noRep),
		y=THBW_noRep$weight_g,
		J=ncol(X_noRep ), 
		X=X_noRep ,
		N_nest=length(unique(THBW_noRep$nest)),
		nest_id=as.numeric(as.factor(THBW_noRep$nest)), 
		N_ind=length(unique(THBW_noRep$bird_id)),
		ind_id=as.numeric(as.factor(THBW_noRep$bird_id)), 
		N_ped = nrow(stan_ped), 
		animal_id =animal_id_noRep, 
		dam = stan_ped$dam,
		sire = stan_ped$sire,
		MSV=MSV, 
		N_NoParentsOffspring=length(NoParentsOffspring), 
		N_ParentsOffspring=length(ParentsOffspring), 
		N_NoOffspring=length(c(NoParentsNoOffspring,ParentsNoOffspring)), 
		NoParentsOffspring=NoParentsOffspring, 
		ParentsOffspring=ParentsOffspring,
		NoOffspring=c(NoParentsNoOffspring,ParentsNoOffspring)
	)


	save(stan_data_DS,stan_data_DS_noRep, stan_data_ped, stan_data_ped_noRep,stan_data_ped_reduced,stan_data_ped_reduced_noRep, file= paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))
}