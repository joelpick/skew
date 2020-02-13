rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
library(MCMCglmm)
library(sn)
library(scales)
library(pbapply)
library(asreml)

wd <- "~/Dropbox/0_blue_tits"
data_wd <- paste0(wd,"/skew/Data/Intermediate/")
source(paste0(wd,"/functions.R"))

load(paste0(data_wd,"stanModTarsus_pedN20191217_1716.Rdata"))
load(paste0(data_wd,"stanModWing_DS20191210_1850.Rdata"))
load(paste0(data_wd,"stanModWing_pedN20191220_0127.Rdata"))

load(paste0(data_wd,"stan_pedN_data.Rdata"))
names(stan_data_THBW)

predict_ASREML <- function(formula,data,model){
	X <- model.matrix(formula,data)
	beta <- asreml_fixed(formula,data,model)
	X %*% beta[,1]
}
asreml_fixed <- function(formula,data,model) {
	X <- model.matrix(formula,data)
	beta <- summary(model,all=TRUE)$coef.fixed
	rownames(beta) <- gsub("_","",rownames(beta))
	return(beta[colnames(X),])
}
simSTRE <- function(levels,dp){
  dMat <- model.matrix(~ as.character(levels)-1)
  N <- length(unique(levels)) 
  pMat <- rst(N,dp=dp)
  return(dMat %*% pMat)
}
simRE <- function(levels,var){
  dMat <- model.matrix(~ as.character(levels)-1)
  N <- length(unique(levels)) 
  pMat <- rnorm(N,0,sqrt(var))
  return(dMat %*% pMat)
}
pars_ST <- function(model, variable){
	summary(model)$summary[paste(c("xi","omega", "alpha", "nu"),variable, sep="_"),c(1,4,8)]
}
pars <- function(model, par){
	out <- summary(model)$summary
	out[grep(par, rownames(out)),c(1,4,8)]
}
mean_CI <- function(x) c(mean(x),quantile(x, c(0.025,0.975)))



stan_mod=mod_stan_wing_pedN
stan_data=stan_data_THBW
ME=TRUE

# sim_preds_ped <- function(stan_mod,stan_data,ME=TRUE){
	
# ped <- cbind(animal=1:stan_data$N_ped,dam=stan_data$dam,sire=stan_data$sire)[-1,]-1
# ped[ped==0]<-NA
	
THBW_egg$malePresent <- as.factor(THBW_egg$male_present)
THBW_egg$hatchDay <- as.factor(THBW_egg$hatch_day)
THBW_egg$clutchSizeC <- scale(THBW_egg$clutch_size, scale=FALSE)
THBW_egg$nestHatchDateC <- scale(THBW_egg$nest_hatch_date, scale=FALSE)
THBW_egg$eggWeightC <- scale(THBW_egg$egg_weight, scale=FALSE)

THBW_egg$animal <- as.factor(stan_data$animal-1)
THBW_egg$nest <- as.factor(THBW_egg$nest)
THBW_egg$bird_id <- as.factor(THBW_egg$bird_id)
#yT_ped <- sim_preds_ped(mod_stan_tarsus_pedN,stan_data_THBW)
# ainv<-asreml.Ainverse(ped)$ginv


ped <- cbind(animal=1:stan_data$N_ped,dam=stan_data$dam,sire=stan_data$sire)[-1,]-1
ped[ped==0]<-NA
ainv<-asreml.Ainverse(ped)$ginv

fixed <- pars(stan_mod,"beta")
nest <- pars_ST(stan_mod, "nest")
animal <- pars(stan_mod, "sigma_A")
resid <- pars_ST(stan_mod, "ind")
me <- pars(stan_mod, "sigma_E")



fixed_form <- formula("~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex +eggWeightC")
X<- model.matrix(fixed_form,THBW_egg)
out<-pbreplicate(100,{
	y <- stan_data$X%*%fixed[,1] + 
	 simSTRE(levels=stan_data$nest,dp=nest[,1]) +
	 rbv(pedigree=ped,G=animal[1]^2)[stan_data$animal-1,]+
	 simSTRE(levels=stan_data$ind_id,dp=resid[,1]) +
	 rnorm(nrow(THBW_egg),0,me[1])

	# y <- X%*%fixed[,1] + 
	# 	 simSTRE(levels=THBW_egg$nest,dp=nest[,1]) +
	# 	 simSTRE(levels=THBW_egg$dam_P,dp=dam_sire[,1]) + 
	# 	 simSTRE(levels=THBW_egg$sire_P,dp=dam_sire[,1]) + 
	# 	 simSTRE(levels=THBW_egg$bird_id,dp=resid[,1]) +
	# 	 rnorm(nrow(THBW_egg),0,me[1])
	}
# stand_skew(y)
# stand_skew(THBW_egg$wing_mm)

THBW_egg$wing_sim <- y

mod <- asreml(
		fixed= update(fixed_form,"wing_sim ~."), 
		random= ~ nest + ped(animal) + bird_id,
		ginverse=list(animal=ainv),
		rcov = ~idv(units) ,
		trace=FALSE,
		data=THBW_egg)
return(c(
	var(predict_ASREML(fixed_form,THBW_egg,mod))
	,
	asreml_fixed(fixed_form,THBW_egg,mod)[2:3,1]
))
})

var(stan_data$X%*%fixed[,1])
fixed[2:3,1]

apply(out,1,mean_CI)
apply(out,1,sd)



mod_simple <- asreml(
		fixed= wing_mm ~ malePresent-1, 
		random= ~ nest + ped(animal) + bird_id,
		ginverse=list(animal=ainv),
		rcov = ~idv(units) ,
		trace=FALSE,
		data=THBW_egg)

asreml_fixed(formula("wing_mm ~ malePresent-1"),THBW_egg,mod_simple)
tapply(THBW_egg$wing_mm,THBW_egg$male_present,mean)
modW



mod_full <- asreml(
		fixed= wing_mm ~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex +eggWeightC, 
		random= ~ nest + ped(animal) + bird_id,
		ginverse=list(animal=ainv),
		rcov = ~idv(units) ,
		trace=FALSE,
		data=THBW_egg)
asreml_fixed(formula("wing_mm ~malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex +eggWeightC"),THBW_egg,mod_full)
beta_asreml <- asreml_fixed(formula("wing_mm ~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex +eggWeightC"),THBW_egg,mod_full)[,1]
var(X%*%beta_asreml)
summary(mod_full,all=TRUE)$coef.fixed

X%*%beta_asreml

plot(THBW_egg$wing_mm~I(X%*%fixed[,1]))
summary(lm(THBW_egg$wing_mm~I(X%*%fixed[,1])))
summary(lm(THBW_egg$wing_mm~I(X%*%beta_asreml)))
plot(I(X%*%fixed[,1])~I(X%*%beta_asreml))








## simulate skewed nest effects
## and a 


rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(lme4)
library(pbapply)

wd <- "~/Dropbox/0_blue_tits"
data_wd <- paste0(wd,"/skew/Data/Intermediate/")

source(paste0(wd,"/functions.R"))

load(paste0(data_wd,"chick_data.Rdata"))
load(paste0(data_wd,"stan_summary_data.Rdata"))

simSTRE <- function(levels,dp){
  dMat <- model.matrix(~ as.character(levels)-1)
  N <- length(unique(levels)) 
  pMat <- rst(N,dp=dp)
  return(dMat %*% pMat)
}

THBW_egg<- subset(THBW,!is.na(sex) & !year%in%c(2019)& !is.na(tarsus_mm)& !is.na(headbill_mm)& !is.na(wing_mm)& !is.na(weight_g))
THBW_egg[THBW_egg$hatch_day==6,"hatch_day"]<-3
THBW_egg_noRep<- subset(THBW_egg,!duplicated(bird_id))

THBW_egg_noRep$malePresent <- as.factor(THBW_egg_noRep$male_present)

X <- model.matrix(~1+malePresent,THBW_egg_noRep)

beta <- skewModT$fixed[1:3,1]

out<-pbreplicate(1000,{
	y <- X%*%beta +
	simSTRE(THBW_egg_noRep$nest, skew_parsT[grep("nest",rownames(skew_parsT)),1]) +
	rst(nrow(THBW_egg_noRep),dp=skew_parsT[grep("ind",rownames(skew_parsT)),1])

	return(summary(lmer(y~malePresent+(1|nest),THBW_egg_noRep))$coef[,1])
})

rowMeans(out)