rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

# library(rstan)
# rstan_options(auto_write = TRUE)
 options(mc.cores = parallel::detectCores())

library(cmdstanr)
library(coda)
library(mvtnorm)

##UNIVARITE

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced

source(paste0(wd,"R/functions.R"))
load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))
load(paste0(wd,"Data/Intermediate/stan_data",if(reduced)"_reduced",".Rdata"))




###-----------------------------###
###-----All models--------------###
###-----------------------------###

for (trait in c("tarsus_mm","wing_mm")){
# trait="tarsus_mm"
	if(trait=="weight_g"){
		stan_data_trait <- stan_data_ped_noRep
		stanModel_pedN <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_pedN2.stan"))
	}else{
		stan_data_trait <- c(list(y=THBW[,trait]),stan_data_ped)
		stanModel_pedN <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_ME_pedN2.stan"))
	}

	asreml_mod <- if(trait=="weight_g"){modM}else if(trait=="tarsus_mm"){modT}else if(trait=="headbill_mm"){modHB}else if(trait=="wing_mm"){modW}

	mod_stan_pedN <- stanModel_pedN$sample(
	    data = stan_data_trait,
	    chains = 4,
		iter_sampling = 6000, 
		iter_warmup = 4000,
		thin= 24,
		refresh = 1000,
		max_treedepth = 12,
		init=function() makeSV(asreml_mod, stan_data_trait, ME=if(trait=="weight_g"){FALSE}else{TRUE} )
	)

	model_z <- as.mcmc.list(lapply(mod_stan_pedN$output_files(), function(x){
		x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
		x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|gamma|sigma|omega|xi|alpha|^nu",names(x1)))])
		return(x2)
	}))

	save(model_z, file= paste0(wd,"Data/Intermediate/stanMod_pedN_",if(reduced)"reduced_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

 	rm(stan_data_trait,stanModel_pedN,asreml_mod,mod_stan_pedN,model_z)

}

load(paste0(wd,"Data/Intermediate/stanMod_pedN_tarsus_mm_20210117_0722.Rdata"))


SV <- makeSV(asreml_mod, stan_data_trait, ME=if(trait=="weight_g"){FALSE}else{TRUE})


sapply(SV,class)
alpha_nest <- rnorm(100,0,1)
sigma_nest <- abs(rcauchy(100,0,10))
nu_nest <- runif(100,4,40)
# alpha_nest<-SV$alpha_nest
# sigma_nest<-SV$sigma_nest
# nu_nest<-SV$nu_nest

b_nu_delta_nest = (alpha_nest / sqrt(1 + alpha_nest^2)) * sqrt(nu_nest/pi) * gamma((nu_nest-1)/2)/gamma(nu_nest/2)
omega_nest = sigma_nest/sqrt(nu_nest/(nu_nest-2) - (b_nu_delta_nest)^2)
xi_nest = 0 - omega_nest*b_nu_delta_nest



###-----------------------------###
###-----Tarsus------------------###
###-----------------------------###


stan_data_tarsus <- c(list(y=THBW$tarsus_mm),stan_data_ped)


# stanModel_pedN_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_pedN2.stan"))

# n_chains <- 1
# mod_stan_tarsus_pedN <- sampling(
# 	stanModel_pedN_ME, 
# 	data = stan_data_tarsus, 
# 	pars =c("beta", 
# 		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
# 		"sigma_A",
# 		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
# 		"sigma_E"
# 		 ),
# 	chains=n_chains, 
# 	iter = 15000, 
# 	warmup = 5000,
# 	thin=20
# 	#control=list(adapt_delta=0.9)#max_treedepth=12,
# 	, init= lapply(1:n_chains,function(x) makeSV(modT, stan_data_tarsus)) 
# 	)
# save(mod_stan_tarsus_pedN, file= paste0(wd,"Data/Intermediate/stanModTarsus_pedN",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
# rm(mod_stan_tarsus_pedN)



stanModel_pedN_ME <- cmdstan_model(paste0(wd,"stan/skew_t_RE_day15_ME_pedN2.stan"))
n_chains <- 4
mod_stan_tarsus_pedN <- stanModel_pedN_ME$sample(
  data = stan_data_tarsus,
  chains = n_chains,
	iter_sampling = 6000, 
	iter_warmup = 4000,
	thin= 24,
	init=function() makeSV(modT, stan_data_tarsus)

)

mod_chains <- as.mcmc.list(lapply(mod_stan_tarsus_pedN$output_files(), function(x){
	x1 <- read.table(x,sep=",", header=TRUE)
	x2 <- as.mcmc(x1[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(x1)))])
	return(x2)
}))

mod_summary <- as.data.frame(mod_stan_tarsus_pedN$summary(c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"sigma_A",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E")
		 ))
save(tarsus_mod_pedN, file= paste0(wd,"Data/Intermediate/stanModTarsus_pedN",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
save(tarsus_summary, file= paste0(wd,"Data/Intermediate/stanModTarsus_pedN_summary",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
tarsus_summary

# rm(stan_data_tarsus,mod_stan_tarsus_pedN,tarsus_mod_pedN,tarsus_summary)


init<- makeSV(modT, stan_data_tarsus)
names(init)
sapply(init[c("sigma_nest","alpha_nest","nu_nest"
	  ,"sigma_ind","alpha_ind","nu_ind"
	  ,"sigma_A"
	  ,"sigma_E"
	)],class)
alpha_nest <- rnorm(100,0,1)
sigma_nest <- abs(rcauchy(100,0,10))
nu_nest <- runif(100,4,40)
# alpha_nest<-init$alpha_nest
# sigma_nest<-init$sigma_nest
# nu_nest<-init$nu_nest

b_nu_delta_nest = (alpha_nest / sqrt(1 + alpha_nest^2)) * sqrt(nu_nest/pi) * gamma((nu_nest-1)/2)/gamma(nu_nest/2)
omega_nest = sigma_nest/sqrt(nu_nest/(nu_nest-2) - (b_nu_delta_nest)^2)
xi_nest = 0 - omega_nest*b_nu_delta_nest




# load("/Users/jpick/Dropbox/0_blue_tits/skew/Data/Intermediate/stanModTarsus_pedN20201212_0952.Rdata") 


# mod_stan_tarsus_pedN2 <- rstan::read_stan_csv(mod_stan_tarsus_pedN$output_files())
# save(rstan::read_stan_csv(mod_stan_tarsus_pedN$output_files()), file= paste0(wd,"Data/Intermediate/stanModTarsus_pedN20201212_0952_2",".Rdata"))

# mod_stan_tarsus_pedN$save_object(file = paste0(wd,"Data/Intermediate/stanModTarsus_pedN20201212_0952_2.RDS"))

# draws <- mod_stan_tarsus_pedN$draws(variables=c("beta", 
# 		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
# 		"sigma_A",
# 		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
# 		"sigma_E","lp__", "treedepth__", "divergent__"
# 		 ))
# mod_stan_tarsus_pedN$return_codes()

# as.data.frame(mod_stan_tarsus_pedN$summary(variables=c("beta", 
# 		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
# 		"sigma_A",
# 		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
# 		"sigma_E","lp__")))

# stanfit <- rstan::read_stan_csv(fit$output_files())

mod_stan_tarsus_pedN$cmdstan_diagnose()

as.data.frame(mod_stan_tarsus_pedN$summary(c("beta", 
		"xi_nest", "omega_nest", "alpha_nest", "nu_nest", "sigma_nest", "gamma_nest",
		"sigma_A",
		"xi_ind", "omega_ind", "alpha_ind", "nu_ind", "sigma_ind", "gamma_ind",
		"sigma_E")
		 ))

dd<-read.csv(mod_stan_tarsus_pedN$output_files()[1],skip=39)
names(dd)[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(dd)))]
read.csv(mod_stan_tarsus_pedN$output_files()[1],skip=39)
/var/folders/15/pdzngbj95lq4sdycdnbpp_k401rlyn/T/RtmpaTAU3d/skew_t_RE_day15_ME_pedN2-202012231149-1-7ebca0.csv

dd1 <- read.csv(mod_stan_tarsus_pedN$output_files()[1],nrows=38)
dd2 <-read.csv(mod_stan_tarsus_pedN$output_files()[1],skip=39)

dd1 <- read.csv("~/Desktop/skew_t_RE_day15_ME_pedN2-202012231149-1-7ebca0.csv",nrows=38, header=FALSE)

dd2 <- read.table(mod_stan_tarsus_pedN$output_files()[1],sep=",", header=TRUE)

dd3 <- as.mcmc(dd2[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(dd2)))])

mod_stan_tarsus_pedN_out <- lapply(mod_stan_tarsus_pedN$output_files(), function(x){
	x1 <- read.table(x,sep=",", header=TRUE)
	x2 <- as.mcmc(x1[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(x1)))])
	return(x2)
})

plot(dd3)

head(dd2[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(dd2)))])
nrow(dd2)

write.csv(dd1, file="~/Desktop/skew_t_RE_day15_ME_pedN2-202012231149-1.csv", row.names=FALSE)
write.table(dd2[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(dd2)))]
, file="~/Desktop/skew_t_RE_day15_ME_pedN2-202012231149-1.csv",append=TRUE, sep=",", row.names=FALSE)

head(dd2[c(1:7,grep("beta//.|gamma|sigma|omega|xi|alpha|^nu",names(dd2)))])


dd <- read.csv("~/Desktop/skew_t_RE_day15_ME_pedN2-202012231149-1-7ebca0.csv",skip=39,nrows=1, header=FALSE)
dd2 <- read.csv("~/Desktop/skew_t_RE_day15_ME_pedN2-202012231149-1-7ebca0.csv",skip=39)
head(dd2)

###-----------------------------###
###-----headbill------------------###
###-----------------------------###


# rm(list=ls())

# options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

# if(Sys.info()["user"]=="jhadfiel"){
# 	wd <- "..."
# }else{
# 	wd <- "~/Dropbox/0_blue_tits/skew/"
# }

# load(paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))

# stanModel_pedN_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_pedN.stan"))


stan_data_HB <- c(list(y=THBW$headbill_mm),stan_data_ped)

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
rm(mod_stan_HB_pedN)

###-----------------------------###
###-----wing------------------###
###-----------------------------###


# rm(list=ls())

# options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

# if(Sys.info()["user"]=="jhadfiel"){
# 	wd <- "..."
# }else{
# 	wd <- "~/Dropbox/0_blue_tits/skew/"
# }

# load(paste0(wd,"Data/Intermediate/stan_pedN_data.Rdata"))

# stanModel_pedN_ME <- stan_model(file = paste0(wd,"stan/skew_t_RE_day15_ME_pedN.stan"))

stan_data_wing <- c(list(y=THBW$wing_mm),stan_data_THBW)


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
rm(mod_stan_wing_pedN)

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
