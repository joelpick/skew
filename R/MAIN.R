rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)


if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}


analysis_options <- list()

analysis_options$reduced <- FALSE
# whether full set of covariates are included in trait models

analysis_options$fixed_z <- as.formula("~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex + eggWeightC")
#full set of fixed effects for trait models

analysis_options$fixed_z_reduced <- as.formula("~ year + timeC + sex + eggWeightC")
#reduced set of fixed effects for trait models

analysis_options$fixed_w <- "sex + malePresent + year+hatchDay + clutchSizeC + nestHatchDateC"

analysis_options$cond_term <- c("year", "sex") 
# terms to condition on in trait and survival models
## we assume that selection occurs within years ans within sexes

analysis_options$cont_term <- c("timeC")  
# terms to control for in trait analysis 

save(analysis_options, file= paste0(wd,"Data/Intermediate/analysis_options.Rdata"))



########################################
## extract day 15 body sizes data
########################################

source(paste0(wd,"R/01_data_prep.R"))
## Produces "Data/Intermediate/chick_data.Rdata" which contains data.frames:
## - 'TBHW' - all day 15 data
## - 'TBHW_noRep' - day 15 data without repeated measures
## - 'ped' - pedigree


########################################
## make data for trait models in stan
########################################

source(paste0(wd,"R/02_stan_data.R"))
## Produces "Data/Intermediate/stan_data_reduced.Rdata" and "Data/Intermediate/stan_data.Rdata" with reduced and filll fixed effect structures respectively. These contain the lists:
## - 'stan_data_DS' - data for stan skew t dam sire analysis with repeated measures
## - 'stan_data_DS_noRep' - data for stan skew t dam sire analysis without repeated measures
## - 'stan_data_ped' - data for stan skew t animal model analysis with repeated measures
## - 'stan_data_ped_noRep' - data for stan skew t animal model analysis without repeated measures


########################################
## generate starting values for trait models in stan using asreml
########################################

source(paste0(wd,"R/03_generate_starting_values.R"))
## Runs similarly structured asreml models assuming normality
## Produces "Data/Intermediate/starting_values.Rdata", which contains lists:
## - modT - tarsus asreml model results
## - modHB - head-bill asreml model results
## - modW - wing asreml model results
## - modM - weight asreml model results


########################################
## Run animal models
########################################

source(paste0(wd,"R/05_animalModel_analysis.R"))
## Runs animal models with normal breeding values and skew t between and within nest effects and normal residuals for models with repeated measures
## Produces


########################################
## Run survival models
########################################

source(paste0(wd,"R/06_survival_analysis.R"))
## Runs
## Produces


3_jarrod_script1.R
3_jarrod_script2.R

3_sim_ped.R



1_MA_data.R