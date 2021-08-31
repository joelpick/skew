rm(list=ls())

options(stringsAsFactors=FALSE)


if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/github/skew/"
}


analysis_options <- list()

analysis_options$reduced <- TRUE
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
