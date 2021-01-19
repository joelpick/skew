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
source(paste0(wd,"R/functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced


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

ainv<-ainverse(ped)

THBW$animal <- factor(THBW$bird_id, levels=ped[,1])
THBW$nest <- as.factor(THBW$nest)
THBW$nest_orig <- as.factor(THBW$nest_orig)
THBW$bird_id <- as.factor(THBW$bird_id)
THBW$dam<-as.factor(THBW$dam)
THBW$sire<-as.factor(THBW$sire)
THBW$eggWeightC <- THBW$egg_weightC

THBW_noRep$animal <- factor(THBW_noRep$bird_id, levels=ped[,1])
THBW_noRep$nest <- as.factor(THBW_noRep$nest)
THBW_noRep$nest_orig <- as.factor(THBW_noRep$nest_orig)
THBW_noRep$bird_id <- as.factor(THBW_noRep$bird_id)
THBW_noRep$dam<-as.factor(THBW_noRep$dam)
THBW_noRep$sire<-as.factor(THBW_noRep$sire)
THBW_noRep$eggWeightC <- THBW_noRep$egg_weightC

THBW$malePresent <- THBW$male_present
THBW$hatchDay <- THBW$hatch_day
THBW$clutchSizeC <- THBW$clutch_sizeC
THBW$nestHatchDateC <- THBW$nest_hatch_dateC

THBW_noRep$malePresent <- THBW_noRep$male_present
THBW_noRep$hatchDay <- THBW_noRep$hatch_day
THBW_noRep$clutchSizeC <- THBW_noRep$clutch_sizeC
THBW_noRep$nestHatchDateC <- THBW_noRep$nest_hatch_dateC




# trait="tarsus_mm";data=THBW;ainv=ainv;reduced=TRUE
asreml_mods <- function(trait, data, ainv, reduced=TRUE){
	fixed <- if(reduced){ 
		formula(paste(trait," ~ year + timeC + sex + eggWeightC") )
	}else{
		formula(paste(trait," ~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex + eggWeightC") )
	}
	
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

modT <- asreml_mods(trait="tarsus_mm", data=THBW, ainv=ainv, reduced=reduced)
modHB <- asreml_mods(trait="headbill_mm", data=THBW, ainv=ainv, reduced=reduced)
modW <- asreml_mods(trait="wing_mm", data=THBW, ainv=ainv, reduced=reduced)
modM <- asreml_mods(trait="weight_g", data=THBW_noRep, ainv=ainv, reduced=reduced)

save(modT,modHB,modW,modM, file= paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
