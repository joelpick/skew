rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)

# wd <- "~/Work"
if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

#load(paste0(wd,"Data/Intermediate/day15_survival_models.Rdata"))
load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))
THBW$year_sex <- paste(THBW$year,THBW$sex)

## THBW is the dataset used for the survival analysis
## mod_tarsus_fledge is model of survival from day 15 to fledging for tarsus
## mod_tarsus_recruit is model of survival from fledging to recruitment for tarsus

## function to apply a function to a subset in a dataframe and return a dataframe
groupFunc_dataFrame <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, c(y, make.row.names=FALSE))
}


## function to simulated normally distributed trait conditioned on other traits
rcmvnorm<-function(n,  
				   df, # dataframe of traits to be simulated from 
				   keep,  # variable to be retained
				   cond=NULL){
	mean <- colMeans(df)
	V <- cov(df)
	cV <- V[keep,keep]-V[keep, -keep]%*%solve(V[-keep, -keep])%*%V[-keep, keep]
 	cM <- sapply(1:nrow(df), function(x) mean[keep]+V[keep, -keep]%*%solve(V[-keep, -keep])%*%(cond[x,]-mean)[-keep])
 	return(rnorm(n=length(cM), mean=cM, sd=sqrt(as.numeric(cV))))
}

N <- 1000
n_it <- 1000

## assign trait

# trait = "tarsus_mmC"
# mod_biv<-mod_tarsus_bv

trait <- "weight_gC"
mod_biv<-mod_weight_bv


# take n_it iterations from posterior of survival models
sample_it <- sample(1:nrow(mod_biv$Sol),n_it, replace=FALSE)

betaS1_full <- mod_biv$Sol[sample_it,grep("at\\.level\\(age,\\ 1\\)", colnames(mod_biv$Sol))]
colnames(betaS1_full)<-gsub("at\\.level\\(age,\\ 1\\):", "", colnames(betaS1_full))
colnames(betaS1_full)[which(colnames(betaS1_full)=="age15")]<-"(Intercept)"

betaS2_full <- mod_biv$Sol[sample_it,grep("at\\.level\\(age,\\ 2\\)", colnames(mod_biv$Sol))]
colnames(betaS2_full)<-gsub(":at\\.level\\(age,\\ 2\\)", "", colnames(betaS2_full))
colnames(betaS2_full)[which(colnames(betaS2_full)=="age25")]<-"(Intercept)"

Vn<-mod_biv$VCV[sample_it,grep("nest", colnames(mod_biv$VCV))]

## empty vectors
S_skew <- vector(length=n_it)
S_normal <- vector(length=n_it)

##  partial derivative of fitness function

w_func<-expression(
  pnorm(i1 + s1 * x + q1*(x^2))*pnorm(i2 + s2 * x + q2*(x^2))
)

wD_func <- D(w_func, "x")

## fixed effects shared between morph and survival analyses
fixedEffects <- formula(~ male_present + hatch_day + clutch_sizeC + nest_hatch_dateC + sex + year)


for(j in 1:n_it){

	## for every year-sex combination, sample with replacement from fixed effects and trait, to keep dependencies between them
	DM_permute <- groupFunc_dataFrame(THBW, "year_sex", function(x){
		samp <- sample(x=1:nrow(x), size=N, replace=TRUE) 		# permuted rows
		trait_skew <- x[samp,trait] 							# permute trait
		DM_year_permute <- model.matrix(fixedEffects,x)[samp,]	# permute design matrix
		i1 <- DM_year_permute %*% betaS1_full[j,colnames(DM_year_permute)]	# intercept for fledging model
		i2 <- DM_year_permute %*% betaS2_full[j,colnames(DM_year_permute)]	# intercepts for recruitment model
		u<-mvrnorm(N, c(0,0), matrix(Vn[j,],2,2))
		i1<-i1+u[,1]
		i2<-i2+u[,2]
		trait_normal <- rcmvnorm(df= cbind(trait_skew,i1,i2), keep = 1, cond = cbind(NA,i1,i2))	# simulate 'normally distributed' trait values with same dependencies among covariates
		return(data.frame(trait_skew=trait_skew,trait_normal=trait_normal,i1=i1,i2=i2,year_sex=unique(x$year_sex)))
	})

	## data for skewed trait
	pars_skew <- data.frame(
		x = DM_permute[,"trait_skew"],
		i1 = DM_permute[,"i1"],
		s1 = betaS1_full[j,trait],
		q1 = betaS1_full[j,paste0(trait,2)],
		i2 = DM_permute[,"i2"],
		s2 = betaS2_full[j,trait],
		q2 = betaS2_full[j,paste0(trait,2)])
	
	## data for normal trait - just replace x
	pars_normal <- pars_skew 
	pars_normal$x <- DM_permute[,"trait_normal"]

	## partial derivatives over skewed trait distribution
	wD_skew <- eval(wD_func, pars_skew)

	## expected fitness over skewed trait distribution
	mW_skew <- eval(w_func, pars_skew)

	## selection gradient for skewed trait
	S_skew[j] <- mean(tapply(wD_skew,DM_permute$year_sex,mean) / tapply(mW_skew,DM_permute$year_sex,mean))

	## partial derivatives over normal trait distribution
	wD_normal <- eval(wD_func, pars_normal)
	
	## expected fitness over normal trait distribution
	mW_normal <- eval(w_func, pars_normal)
	
	## selection gradient for normal trait
	S_normal[j] <- mean(tapply(wD_normal,DM_permute$year_sex,mean) / tapply(mW_normal,DM_permute$year_sex,mean))

	## progress bar
	if((j/n_it) %in% ((1:10)/10)) cat(paste0((j/n_it) * 100,"% "))
}

save(S_normal,S_skew,file=paste0(wd,"Data/Intermediate/selection_gradient_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

par(mfrow=c(3,1))
hist(S_normal, xlim=range(c(S_normal,S_skew)), breaks=50)
hist(S_skew, xlim=range(c(S_normal,S_skew)), breaks=50)
hist((S_skew-S_normal), breaks=50)

sum(S_skew>S_normal)/n_it
sum(S_normal<0)/n_it
sum(S_skew<0)/n_it