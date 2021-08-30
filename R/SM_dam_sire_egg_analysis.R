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

## adapted pin function for asreml
pin<-function (object, transform) {
    pframe <- as.list(object$gammas)
    names(pframe) <- sub("!.*\\..*$", "", names(pframe))
	names(pframe) <- sub(".*\\((\\w*)\\).*$", "\\1", names(pframe))
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
        pframe)
    X <- as.vector(attr(tvalue, "gradient"))
    tname <- if (length(transform) == 3) 
        transform[[2]]
    else ""
    Vmat <- object$ai
    n <- length(pframe)
    i <- rep(1:n, 1:n)
    j <- sequence(1:n)
    k <- 1 + (i > j)
    se <- sqrt(sum(Vmat * X[i] * X[j] * k))
    data.frame(row.names = tname, Estimate = tvalue, SE = se)
}

## transforms dam/sire variances into Va and other variances
sire_to_animal_asreml <- function(mod){
	Va <- pin(mod, V_a~sire*4)
	Ve <- pin(mod, V_e~R-sire*2)
	Vo <- pin(mod, V_origin~dam-sire)
	Vr <- pin(mod, V_rear~nest)
	return( rbind(Va, Vr, Vo, Ve))
}
## transforms dam/sire multimembership variances into Va and other variances
## make this work for MM models
DSMM_to_animal_asreml <- function(mod, ME=TRUE, NO=TRUE){
	animal <- pin(mod, animal~dam*4)
	if(ME){
		bird_id <- pin(mod, bird_id~bird_id-dam*2)
		R <- pin(mod, R~R)
	}else{
		R <- pin(mod, R~R-dam*2)
	}
	if(NO) nest_orig <- pin(mod, nest_orig~nest_orig)
	nest <- pin(mod, nest~nest)
	out <- rbind(animal, nest)
	if(NO) out <- rbind(out,nest_orig)
	if(ME) out <- rbind(out,bird_id)
	out <- rbind(out,R)
	return(as.matrix(out))
}
## function to get variances out of asreml object
asreml_varcomp <- function(mod){
	x <- summary(mod)$varcomp[,c(2,3)]
	### add on standard errors
	rownames(x) <- sub("!.*\\..*$", "", rownames(x))
	rownames(x) <- sub(".*\\((\\w*)\\).*$", "\\1", rownames(x))
	as.matrix(x[rownames(x)!="R!variance",])
}
inv_hessian_varcomp <- function(mod){
	x <- summary(mod)$varcomp
	### add on standard errors
	rownames(x) <- sub("!.*\\..*$", "", rownames(x))
	rownames(x) <- sub(".*\\((\\w*)\\).*$", "\\1", rownames(x))
	exclude <- which(rownames(x)=="R!variance")
	inv_hess <- Tri2M(mod$ai,FALSE,TRUE)[-exclude,-exclude]
	colnames(inv_hess)<-rownames(inv_hess)<-rownames(x)[-exclude]
	return(inv_hess)
}

## function to get fixed effects out of asreml object
asreml_fixed <- function(formula,data,model) {
	X <- model.matrix(formula,data)
	beta <- summary(model,all=TRUE)$coef.fixed
	rownames(beta) <- gsub("_","",rownames(beta))
	return(beta[colnames(X),])
}

inv_hessian_fixed <- function(formula,data,model) {
	X <- model.matrix(formula,data)
	beta <- summary(model,all=TRUE)$coef.fixed
	rownames(beta) <- gsub("_","",rownames(beta))
	inv_hess <- Tri2M(model$Cfixed, lower.tri=FALSE, reverse=TRUE, diag=TRUE)
	colnames(inv_hess)<-rownames(inv_hess)<-rownames(beta)
	return(inv_hess[colnames(X),colnames(X)])

}

## fixed effect predictions for a given formula
predict_ASREML <- function(formula,data,model){
	X <- model.matrix(formula,data)
	beta <- asreml_fixed(formula,data,model)
	X %*% beta[,1]
}
rbind_match<-function(...){
  x <- list(...)
  names1 <- names(x[[1]])
  y <- lapply(x, function(y) y[names1] )
  do.call(rbind,y)
}

#trait="tarsus_mm";age=15; data=THBW_egg; ainv=ainv; measurement_error=TRUE
nest_orig_egg <- function(trait,age,data,ainv,measurement_error=TRUE){
	fixed <- formula(paste(trait," ~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex") )
	random_DS <- formula("~ nest_orig + nest + str(~dam+and(sire), ~us(1):id(dam))")
	random_A <- formula("~ nest_orig + nest + ped(animal)")
	if(trait!="weight_g" & measurement_error) {
		random_DS <- update(random_DS,"~ . + bird_id")
		random_A <- update(random_A,"~ . + bird_id")
	}

## models with nest_of_origin	
	## multi-membership model, where dam and sire variances are constrained to be equal
	mod_orig_DS <- asreml(
		fixed=fixed,  
		random= random_DS,
		rcov = ~idv(units) ,
		trace=FALSE,
		data=data,model.frame=TRUE)

	mod_orig_A <- asreml(
		fixed=fixed,  
		random= random_A,
		ginverse=list(animal=ainv),
		rcov = ~idv(units) ,
		trace=FALSE,
		data=data,keep.order=TRUE)

## models with egg size as a covariate
	fixed_egg <- update(fixed,"~ . +eggWeightC")
	mod_orig_egg_DS <- asreml(
		fixed=fixed_egg,  
		random= random_DS,
		rcov = ~idv(units) ,
		trace=FALSE,
		data=data)

	mod_orig_egg_A <- asreml(
		fixed=fixed_egg,  
		random= random_A,
		ginverse=list(animal=ainv),
		rcov = ~idv(units) ,
		trace=FALSE,
		data=data)

# models with egg size but no nest_of_origin
	random_DS_egg <- update(random_DS,"~ . -nest_orig")
	random_A_egg <- update(random_A,"~ . -nest_orig")
	mod_egg_DS <- asreml(
		fixed=fixed_egg,  
		random= random_DS_egg,
		rcov = ~idv(units) ,
		trace=FALSE,
		Cfixed=TRUE,
		data=data)

	mod_egg_A <- asreml(
		fixed=fixed_egg,  
		random= random_A_egg,
		ginverse=list(animal=ainv),
		rcov = ~idv(units) ,
		trace=FALSE,
		Cfixed=TRUE,
		data=data)

	out <- list()
	for(i in 1:2){
		out[[i]] <- round(rbind_match(
			mod_orig_DS=DSMM_to_animal_asreml(mod_orig_DS,ME=measurement_error)[,i],
			mod_orig_A=asreml_varcomp(mod_orig_A)[,i],
			mod_orig_egg_DS=DSMM_to_animal_asreml(mod_orig_egg_DS,ME=measurement_error)[,i],
			mod_orig_egg_A=asreml_varcomp(mod_orig_egg_A)[,i],
			mod_egg_DS=c(nest_orig=0,DSMM_to_animal_asreml(mod_egg_DS,NO=FALSE,ME=measurement_error)[,i]),
			mod_egg_A=c(nest_orig=0,asreml_varcomp(mod_egg_A)[,i])
			# c(nest_orig=0,asreml_varcomp(mod_egg_A))[c(1,3,2,4:(length(asreml_varcomp(mod_egg_A))+2))]
		),5)
	}

	for(i in 1:2){
	out[[i+2]] = rbind(
			mod_orig_DS=rbind(asreml_fixed(fixed,data,mod_orig_DS),eggWeightC=0)[,i],
			mod_orig_A=rbind(asreml_fixed(fixed,data,mod_orig_A),eggWeightC=0)[,i],
			mod_orig_egg_DS=asreml_fixed(fixed_egg,data,mod_orig_egg_DS)[,i],
			mod_orig_egg_A=asreml_fixed(fixed_egg,data,mod_orig_egg_A)[,i],
			mod_egg_DS=asreml_fixed(fixed_egg,data,mod_egg_DS)[,i],
			mod_egg_A=asreml_fixed(fixed_egg,data,mod_egg_A)[,i]
			)
	}
	names(out)<-c("vars_mean","vars_SE","fixed_mean","fixed_SE")

	fixed_time <- update(fixed,"~ . -timeC")
	fixed_egg_time <- update(fixed_egg,"~ . -timeC")

	out$fixed_var <- c(
		mod_orig_DS= var(predict_ASREML(fixed_time,data,mod_orig_DS)),
		mod_orig_A=var(predict_ASREML(fixed_time,data,mod_orig_A)),
		mod_orig_egg_DS=var(predict_ASREML(fixed_egg_time,data,mod_orig_egg_DS)),
		mod_orig_egg_A=var(predict_ASREML(fixed_egg_time,data,mod_orig_egg_A)),
		mod_egg_DS=var(predict_ASREML(fixed_egg_time,data,mod_egg_DS)),
		mod_egg_A=var(predict_ASREML(fixed_egg_time,data,mod_egg_A)))

		out$mod_egg_DS_Fi <- inv_hessian_fixed(fixed_egg,data,mod_egg_DS)
		out$mod_egg_A_Fi <- inv_hessian_fixed(fixed_egg,data,mod_egg_A)
		
		out$mod_egg_DS_Ri <- inv_hessian_varcomp(mod_egg_DS)
		out$mod_egg_A_Ri <- inv_hessian_varcomp(mod_egg_A)
		return(out)
}

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

THBW_egg<- subset(THBW_egg,!is.na(sex) & !year%in%c(2019)& !is.na(tarsus_mm)& !is.na(headbill_mm)& !is.na(wing_mm)& !is.na(weight_g)& !is.na(sire_P))

## pedigree
ped <- insertPed(prunePed(orderPed(tPED_use[1:3]), unique(THBW_egg$bird_id), make.base = TRUE))
names(ped) <- c("animal","dam", "sire")
ainv<-asreml.Ainverse(ped)$ginv

## random effects
THBW_egg$animal <- factor(THBW_egg$bird_id, levels=ped[,1])
THBW_egg$nest <- as.factor(THBW_egg$nest)
THBW_egg$nest_orig <- as.factor(THBW_egg$nest_orig)
THBW_egg$bird_id <- as.factor(THBW_egg$bird_id)

mm_levels <- unique(c(THBW_egg$sire_P,THBW_egg$dam_P))
THBW_egg$dam <- factor(THBW_egg$dam_P, levels=mm_levels)
THBW_egg$sire <- factor(THBW_egg$sire_P, levels=mm_levels)

## fixed effects
# assign individuals with hatch day 6, hatch day 3
THBW_egg[THBW_egg$hatch_day==6,"hatch_day"]<-3

THBW_egg$malePresent <- as.factor(THBW_egg$male_present)
THBW_egg$year <- as.factor(THBW_egg$year)
THBW_egg$sex <- as.factor(THBW_egg$sex)
THBW_egg$hatchDay <- as.factor(THBW_egg$hatch_day)
THBW_egg$timeC <- scale(THBW_egg$time_hr, scale=FALSE)
THBW_egg$clutchSizeC <- scale(THBW_egg$clutch_size, scale=FALSE)
THBW_egg$nestHatchDateC <- scale(THBW_egg$nest_hatch_date, scale=FALSE)
THBW_egg$eggWeightC <- scale(THBW_egg$egg_weight, scale=FALSE)


# models
modT <- nest_orig_egg(trait="tarsus_mm",age=15, data=THBW_egg, ainv=ainv_thbw)
modHB <- nest_orig_egg(trait="headbill_mm",age=15, data=THBW_egg, ainv=ainv_thbw)
modW <- nest_orig_egg(trait="wing_mm",age=15, data=THBW_egg, ainv=ainv_thbw)

THBW_egg_noRep<- subset(THBW_egg,!duplicated(bird_id))
modM <- nest_orig_egg(trait="weight_g",age=15, data=THBW_egg_noRep, ainv=ainv_thbw, measurement_error=FALSE)

save(modT,modHB,modW,modM,file=paste0(wd,"Data/Intermediate/dam_sire_egg.Rdata"))


### monogamy

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

## pedigree
THBW_egg<- subset(THBW_egg,!is.na(sex) & !year%in%c(2019)& !is.na(tarsus_mm)& !is.na(headbill_mm)& !is.na(wing_mm)& !is.na(weight_g)& !is.na(sire_P))

ped <- insertPed(prunePed(orderPed(tPED_use[c(1,4,5)]), unique(THBW_egg$bird_id), make.base = TRUE))
names(ped) <- c("animal","dam", "sire")
ainv<-asreml.Ainverse(ped)$ginv

## random effects
THBW_egg$animal <- factor(THBW_egg$bird_id, levels=ped[,1])
THBW_egg$nest <- as.factor(THBW_egg$nest)
THBW_egg$nest_orig <- as.factor(THBW_egg$nest_orig)
THBW_egg$bird_id <- as.factor(THBW_egg$bird_id)

mm_levels_M <- unique(c(THBW_egg$sire_M,THBW_egg$dam_M))
THBW_egg$dam <- factor(THBW_egg$dam_M, levels=mm_levels_M)
THBW_egg$sire <- factor(THBW_egg$sire_M, levels=mm_levels_M)

## fixed effects
# assign individuals with hatch day 6, hatch day 3
THBW_egg[THBW_egg$hatch_day==6,"hatch_day"]<-3

THBW_egg$malePresent <- as.factor(THBW_egg$male_present)
THBW_egg$year <- as.factor(THBW_egg$year)
THBW_egg$sex <- as.factor(THBW_egg$sex)
THBW_egg$hatchDay <- as.factor(THBW_egg$hatch_day)
THBW_egg$timeC <- scale(THBW_egg$time_hr, scale=FALSE)
THBW_egg$clutchSizeC <- scale(THBW_egg$clutch_size, scale=FALSE)
THBW_egg$nestHatchDateC <- scale(THBW_egg$nest_hatch_date, scale=FALSE)
THBW_egg$eggWeightC <- scale(THBW_egg$egg_weight, scale=FALSE)


modT_M <- nest_orig_egg(trait="tarsus_mm",age=15, data=THBW_egg, ainv=ainv_thbw_M)
modHB_M <- nest_orig_egg(trait="headbill_mm",age=15, data=THBW_egg, ainv=ainv_thbw_M)
modW_M <- nest_orig_egg(trait="wing_mm",age=15, data=THBW_egg, ainv=ainv_thbw_M)

THBW_egg_noRep<- subset(THBW_egg,!duplicated(bird_id))
modM_M <- nest_orig_egg(trait="weight_g",age=15, data=THBW_egg_noRep, ainv=ainv_thbw_M, measurement_error=FALSE)

save(modT_M,modHB_M,modW_M,modM_M,file=paste0(wd,"Data/Intermediate/dam_sire_egg_M.Rdata"))
