rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MasterBayes)
library(asreml)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
source(paste0(wd,"R/00_functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))


## function to get variances out of asreml object
asreml_varcomp <- function(mod){
	x <- summary(mod)$varcomp[,c(1:3)]
	x <- as.matrix(x[rownames(x)!="units!R",])
	rownames(x) <- sub("!.+$", "", rownames(x))
	rownames(x) <- sub(".*\\((\\w*)\\,.+\\).*$", "\\1", rownames(x))
	return(x)
}
inv_hessian_varcomp <- function(mod){
	inv_hess <- mod$ai
	inv_hess <- as.matrix(inv_hess[rownames(inv_hess)!="units!R",colnames(inv_hess)!="units!R"])
	ih_names <- rownames(inv_hess)
	ih_names <- sub("!.+$", "", ih_names)
	ih_names <- sub(".*\\((\\w*)\\,.+\\).*$", "\\1", ih_names)
	colnames(inv_hess)<-rownames(inv_hess)<-ih_names
	return(inv_hess)
}

## function to get fixed effects out of asreml object
asreml_fixed <- function(formula,data,model) {
	X <- model.matrix(formula,data)
	beta <- summary(model,coef=TRUE)$coef.fixed
	rownames(beta) <- gsub("_","",rownames(beta))
	return(beta[colnames(X),])
}

inv_hessian_fixed <- function(formula,data,model) {
	X <- model.matrix(formula,data)
	inv_hess <- model$Cfixed
	colnames(inv_hess)<-rownames(inv_hess)<-gsub("_","",rownames(inv_hess))
	return(inv_hess[colnames(X),colnames(X)])

}

# trait="tarsus_mm";data=THBW;ainv=ainv;reduced=TRUE
asreml_mods <- function(trait, data, ainv, fixed_z){
	fixed <- as.formula(paste(c(trait,fixed_z) ,collapse=" "))
	
	# random_DS <- formula("~ nest + str(dam + and(sire), us(1):id(dam))")
	random_DS <- formula("~ nest + nest_orig + sire")
	random_A <- formula("~ nest + vm(animal,ainv)")
	if(trait!="weight_g") {
		random_DS <- update(random_DS,"~ . + bird_id")
		random_A <- update(random_A,"~ . + bird_id")
	}

	asreml.options(Cfixed=TRUE,keep.order=TRUE,
		trace=FALSE)

	mod_DS <- asreml(
		fixed=fixed,  
		random= random_DS,
		residual = ~idv(units),
		data=data)
	
	mod_DS_out <- list(
		V = asreml_varcomp(mod_DS)[,1],
		Vi = inv_hessian_varcomp(mod_DS),
		F = asreml_fixed(fixed,data,mod_DS)[,1],
		Fi = inv_hessian_fixed(fixed,data,mod_DS)
	)

	mod_A <- asreml(
		fixed=fixed,  
		random= random_A,
		residual = ~idv(units),
		data=data)

	mod_A_out <- list(
		V = asreml_varcomp(mod_A)[,1],
		Vi = inv_hessian_varcomp(mod_A),
		F = asreml_fixed(fixed,data,mod_A)[,1],
		Fi = inv_hessian_fixed(fixed,data,mod_A)
	)

	out <- list(mod_DS = mod_DS_out, mod_A = mod_A_out)
	return(out)
}



ainv<-ainverse(ped)


## make factors for asreml random effects
THBW$animal <- factor(THBW$bird_id, levels=ped[,1])
THBW$nest <- as.factor(THBW$nest)
THBW$nest_orig <- as.factor(THBW$nest_orig)
THBW$bird_id <- as.factor(THBW$bird_id)
THBW$dam<-as.factor(THBW$dam)
THBW$sire<-as.factor(THBW$sire)

THBW_noRep$animal <- factor(THBW_noRep$bird_id, levels=ped[,1])
THBW_noRep$nest <- as.factor(THBW_noRep$nest)
THBW_noRep$nest_orig <- as.factor(THBW_noRep$nest_orig)
THBW_noRep$bird_id <- as.factor(THBW_noRep$bird_id)
THBW_noRep$dam<-as.factor(THBW_noRep$dam)
THBW_noRep$sire<-as.factor(THBW_noRep$sire)


for(reduced in c(TRUE,FALSE)){

	fixed_z <- if(reduced){ analysis_options$fixed_z_reduced }else{ analysis_options$fixed_z }
	
	modT <- asreml_mods(trait="tarsus_mm", data=THBW, ainv=ainv, fixed_z=fixed_z)
	modHB <- asreml_mods(trait="headbill_mm", data=THBW, ainv=ainv, fixed_z=fixed_z)
	modW <- asreml_mods(trait="wing_mm", data=THBW, ainv=ainv, fixed_z=fixed_z)
	modM <- asreml_mods(trait="weight_g", data=THBW_noRep, ainv=ainv, fixed_z=fixed_z)

	save(modT,modHB,modW,modM, file= paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
}