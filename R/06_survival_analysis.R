rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(rstan)
library(parallel)
library(cmdstanr)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))
fixed_w <- analysis_options$fixed_w

model_ME <- TRUE ## whether to account for measurement error in body size traits

## make data in long format for survival analysis
fledge_dat<-THBW_noRep
colnames(fledge_dat)[which(colnames(fledge_dat)=="fledge")]<-"survival"
colnames(fledge_dat)[which(colnames(fledge_dat)=="recruit")]<-"age"
fledge_dat$age<-15

recruit_dat<-THBW_noRep[which(THBW_noRep$fledge),]
colnames(recruit_dat)[which(colnames(recruit_dat)=="recruit")]<-"survival"
colnames(recruit_dat)[which(colnames(recruit_dat)=="fledge")]<-"age"
recruit_dat$age<-25

surv_dat<-rbind(fledge_dat, recruit_dat)
surv_dat$age<-as.factor(surv_dat$age)



if(!model_ME){
	prior_bv <- list( R = list(V = 1, fix = 1), 
				   G = list(
				   	G1 = list(V = diag(2), nu = 2, alpha.mu=c(0,0), alpha.V=diag(2)*100))
				 )
	it_scale <- 100

	mclapply(c("tarsus_mm","headbill_mm","wing_mm","weight_g"), function(trait){

		traitC <- paste0(trait,"C")
		traitC2 <- paste0(trait,"C2")

		fixed <- as.formula(paste0("survival ~ -1+ at.level(age,1):(age + ", traitC," + ", traitC2, " + ",fixed_w,")+ at.level(age,2):(age + ", traitC," + ", traitC2, " + ",fixed_w,")"))


		model_w <- MCMCglmm(fixed=fixed, random=~us(age):nest, data=surv_dat, prior=prior_bv, family="threshold", nitt=13000*it_scale, thin=10*it_scale, burnin=3000*it_scale)

		save(model_w, file = paste0(wd,"Data/Intermediate/day15_survival_model_",trait,".Rdata"))

	}, mc.cores = 4)

}





if(model_ME){
	for (trait in c("tarsus_mm","headbill_mm","weight_g","wing_mm")){

	# trait="tarsus_mm"
		if(trait=="weight_g"){

			traitC <- paste0(trait,"C")
			traitC2 <- paste0(trait,"C2")

			fixed <- as.formula(paste0(" ~ -1+ at.level(age,1):(age + ",fixed_w,")+ at.level(age,2):(age + ",fixed_w,") +  age:", traitC," + age:", traitC2))

			X <- model.matrix(fixed,surv_dat)
			X <- X[,apply(X,2,function(x)sum(x==0))<nrow(X)]

			stan_dat <- list(
				N = nrow(surv_dat),
				J = ncol(X),
				K = 2,
				N_nest = n_unique(surv_dat$nest),
				Y = as.numeric(surv_dat$survival),
				X = X,
				nest_id = as.numeric(as.factor(surv_dat$nest)), 
				age = as.numeric(surv_dat$age)
				)

			stanModel_probit <- cmdstan_model(paste0(wd,"stan/survival_multi_probit.stan"))#stan_model(file = paste0(wd,"/stan/survival/survival_multi_probit.stan"))

			mod_stan_probit <- stanModel_probit$sample( 
				data = stan_dat, 
				chains=4, 
				iter_sampling = 2500, 
				iter_warmup = 2500,
				thin= 10,
				refresh = 1000
				,init = function() list(
					sigma_nest=array(abs(rnorm(2,0,0.1)),dim=2), 
			 	 	L_corr_nest=as.array(t(chol(matrix(c(1,0,0,1),2)))),
			 	 	beta_tilde=array(rnorm(stan_dat$J,0,0.1),dim=c(stan_dat$J)),
			 	 	nest_scaled=array(rnorm(stan_dat$N_nest*stan_dat$K,0,0.1),c(stan_dat$N_nest,stan_dat$K))
				 	)
			)
			model_w <- as.mcmc.list(lapply(mod_stan_probit$output_files(), function(x){
				x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
				x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|Sigma",names(x1)))])
				return(x2)
			}))

		}else{

			fixed_ME <- as.formula(paste0(" ~ -1+ at.level(age,1):(age + ",fixed_w,")+ at.level(age,2):(age + ",fixed_w,")"))
			X_ME <- model.matrix(fixed_ME,surv_dat)
			X_ME <- X_ME[,apply(X_ME,2,function(x)sum(x==0))<nrow(X_ME)]

			stan_dat <- list(
				N = nrow(surv_dat),
				J = ncol(X_ME),
				K = 2,
				N_nest = n_unique(surv_dat$nest),
				Y = as.numeric(surv_dat$survival),
				X = X_ME,
				nest_id = as.numeric(as.factor(surv_dat$nest)), 
				age = as.numeric(surv_dat$age),

				N_x1 = nrow(THBW),
				N_ind = n_unique(THBW$bird_id),
				x1 = THBW[,trait],
				
				x1_ind_id =  as.numeric(as.factor(THBW$bird_id)),
				ind_id =  as.numeric(as.factor(surv_dat$bird_id)),

				X1= model.matrix(~surv_dat$age-1)
				)


			stanModel_probit <- cmdstan_model(paste0(wd,"stan/survival_multi_probit_ME.stan"))

			mod_stan_probit <- stanModel_probit$sample( 
				data = stan_dat, 
				chains=4, 
				iter_sampling = 2500, 
				iter_warmup = 2500,
				thin= 10,
				refresh = 1000
				,init = function() list(
					sigma_nest=array(abs(rnorm(2,0,0.1)),dim=2), 
			 	 	L_corr_nest=as.array(t(chol(matrix(c(1,0,0,1),2)))),
			 	 	beta_tilde=array(rnorm(stan_dat$J,0,0.1),dim=c(stan_dat$J)),
			 	 	nest_scaled=array(rnorm(stan_dat$N_nest*stan_dat$K,0,0.1),c(stan_dat$N_nest,stan_dat$K)),
			 	 	beta0_x1 =rnorm(1,0,0.1),
					x1_hat_scaled = array(rnorm(stan_dat$N_ind,0,0.1),dim=stan_dat$N_ind),
					sigma_E_x1 = abs(rnorm(1,0,0.1)),
					sigma_ind_x1 = abs(rnorm(1,0,0.1)),
					beta_x1_hat = array(rnorm(2,0,0.1),dim=2),
					beta_x1_hat2 = array(rnorm(2,0,0.1),dim=2)
				 	)
			)



			model_w <- as.mcmc.list(lapply(mod_stan_probit$output_files(), function(x){
				x1 <- read.table(x,sep=",", header=TRUE) # using read table gets rid of commented lines
				x2 <- as.mcmc(x1[c(1:7,grep("beta\\.|beta_x1_hat|beta_x1_hat2|Sigma_nest|beta0_x1|sigma_ind_x1|sigma_E_x1|Sigma",names(x1)))])
				return(x2)
			}))
		}

		save(model_w, file= paste0(wd,"Data/Intermediate/day15_survival_ME_",trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))

	 	rm(stan_dat,stanModel_probit,mod_stan_probit,model_w)

	}
}




# library(MCMCvis)

# 	MCMCtrace(model_w,ind=TRUE,Rhat=TRUE,n.eff=TRUE,pdf=TRUE)







# trait="tarsus_mm"









# traceplot(mod_stan_probit_ME,pars=c("beta_x1_hat","beta_x1_hat2","Sigma_nest","beta0_x1","sigma_ind_x1","sigma_E_x1"))
# pairs(mod_stan_probit_ME,pars=c("beta_x1_hat","beta_x1_hat2","Sigma_nest","beta0_x1","sigma_ind_x1","sigma_E_x1"))


# summary(mod_stan_probit_ME)$summary
# 2421.23/60/60

# load(paste0(wd,"Data/Intermediate/day15_survival_model_",trait,".Rdata"))

# summary(model_w)


# par(mar=c(4,10,1,1))
# effectPlot(summary(mod_stan_probit_ME)$summary[c(1,31,33,2:16,32,34,17:30,35:38),c(1,4,8)])
# effectPlot(summary(mod_stan_probit)$summary[1:38,c(1,4,8)],add=TRUE,offset=-0.15,col="blue")
# effectPlot(rbind(summary(model_w)$sol[,1:3],summary(model_w)$Gc[,1:3]),add=TRUE,offset=-0.3,col="red")

# cbind(summary(mod_stan_probit_ME)$summary[c(31,33,32,34),c(1)],
# summary(mod_stan_probit)$summary[c(2:3,19:20),c(1)],
# summary(model_w)$sol[c(2:3,19:20),1]
# )

# summary(mod_stan_probit_ME)$summary[c(31,33,32,34),c(1)]/ summary(mod_stan_probit)$summary[c(2:3,19:20),c(1)]

# summary(mod_stan_probit_uni)$summary
# summary(mod_stan_probit)$summary

# #list(list(sigma_nest=array(c(0.5,0.5),dim=2), 
# # L_corr_nest=as.array(t(chol(matrix(c(1,0.5,0.5,1),2)))),
# # beta_tilde=array(rnorm(stan_dat$J),dim=c(stan_dat$J)),
# # nest_scaled=array(rnorm(stan_dat$N_nest*stan_dat$K),c(stan_dat$N_nest,stan_dat$K))) )


# X_15<- model.matrix( ~ 1 + weight_gC + weight_gC2 + sex + male_present + year + hatch_day + clutchSizeC + nestHatchDateC, data=fledge_dat)[as.numeric(as.factor(fledge_dat[,"bird_id"])),]

# stan_dat2 <- list(
# 	N = nrow(surv_dat),
# 	J = ncol(X_15),
# 	K = 2,
# 	N_nest = n_unique(surv_dat$nest),
# 	N_ind = n_unique(surv_dat$bird_id),
# 	Y = as.numeric(surv_dat$survival),
# 	X=X_15,
# 	nest_id = as.numeric(as.factor(surv_dat$nest)), 
# 	ind_id = as.numeric(as.factor(surv_dat$bird_id)), 
# 	age = as.numeric(surv_dat$age)
# 	)





# int_fun <- function(...) list(sigma_nest=sqrt(c(1.4,0.18)))
# stanModel_probit2 <- stan_model(file = paste0(wd,"/stan/survival/survival_multi_probit_test2.stan"))
# mod_stan_probit2 <- sampling(stanModel_probit, data = stan_dat, chains=1, iter = 5000, warmup = 2000,pars=c("beta","Sigma_nest")
# 	 , init=list(list(sigma_nest=array(c(0.5,0.5),dim=2), 
# 	 	 	L_corr_nest=as.array(t(chol(matrix(c(1,0.5,0.5,1),2)))),
# 	 	 	beta_tilde=array(rnorm(stan_dat$J*stan_dat$K),dim=c(stan_dat$J,stan_dat$K)),
# 	 	 	nest_scaled=array(rnorm(stan_dat$N_nest*stan_dat$K),c(stan_dat$N_nest,stan_dat$K))) )
# 	)



# summary(mod_stan_probit)$summary[,c(1,4,8,9)]
# hist(qnorm(x))

 
# me<-extract(mod_stan_probit,permute=FALSE)
# dim(me)
# me[,,"beta[1]"]

# X_15<- model.matrix( ~ 1 + weight_gC + weight_gC2 + sex + male_present + year + hatch_day + clutchSizeC + nestHatchDateC, data=fledge_dat)
# stan_dat_15 <- list(
# 	N = nrow(fledge_dat),
# 	J = ncol(X_15),
# 	N_nest = n_unique(fledge_dat$nest),
# 	Y = as.numeric(fledge_dat$survival),
# 	X=X_15,
# 	nest_id = as.numeric(as.factor(fledge_dat$nest))
# 	)
# stanModel_probit_U_RE <- stan_model(file = paste0(wd,"/stan/survival/survival_uni_probit.stan"))
# mod_15 <- sampling(stanModel_probit_U_RE, data = stan_dat_15, chains=1, iter = 2000, warmup = 1000,pars=c("beta", "sigma_nest"))
# summary(mod_15)$summary[,c(1,4,8,9)]

# prior <- list( R = list(V = 1, fix = 1),
# 			G = list (G1=list(V = 1, nu = 2, alpha.mu=0, alpha.V=100)))
# # it_scale=5
# mod_15_m<-MCMCglmm(survival~ 1 + weight_gC + weight_gC2 + sex + male_present + year + hatch_day + clutchSizeC + nestHatchDateC,random=~nest,family="threshold",prior=prior, data=fledge_dat,nitt=13000*it_scale, thin=10*it_scale, burnin=3000*it_scale)
# summary(mod_15_m)

# effectPlot(rbind(summary(mod_15_m)$sol[,1:3],sqrt(summary(mod_15_m)$Gc[,1:3])),add=TRUE,offset=0.1,col="red")


# effectPlot(summary(mod_15)$summary[1:17,c(1,4,8)])


# col="black"
# col[1:10]