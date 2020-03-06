rm(list=ls())

# source("~/Work/Skew/R/selection_gradients_skew_jarrod.R")

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)
library(rstan)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

trait<-"weight_g"
re_run<-FALSE
save<-FALSE
posterior_mean<-TRUE
save_plot<-FALSE
cond_term<-c("year", "sex") # terms to condition on 
cont_term<-c("timeC")   	# terms to control for
model_moments<-FALSE 		# when calculating selection gradients should model-based or sample moments be used

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))
# survival models (including 2010 when eggs weren't weight) and THBW is the data

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))
# stan_data_weight$X is the design matrix and THBW_egg_noRep the data
# excluding years in which eggs weren't measured and repeat measures

if(trait=="weight_g"){
 model_w<-mod_weight_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModWeightPedN20191220_0905.Rdata")))
}
if(trait=="tarsus_mm"){
 model_w<-mod_tarsus_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20191217_1716.Rdata"))) 
}
rm("mod_weight_bv", "mod_tarsus_bv")
# read in fitness models (model_w) and trait models (model_z): 
# Note the trait model excludes 2010 in which eggs were not weighed

pred_cond<-colnames(stan_data_weight$X)[unlist(sapply(cond_term, x=colnames(stan_data_weight$X), grep))] # terms to condition on 
pred_cont<-colnames(stan_data_weight$X)[unlist(sapply(cont_term, x=colnames(stan_data_weight$X), grep))] # terms to control for 
# gets the names of the effects associated with conditioning and controlling variables for the trait models 

pred_cond_pos<-match(pred_cond, colnames(stan_data_weight$X))
pred_cont_pos<-match(pred_cont, colnames(stan_data_weight$X))
# gets the position of those effects in the design matrix 

pred_pos<-match(setdiff(colnames(stan_data_weight$X), c(pred_cond, pred_cont_pos)), colnames(stan_data_weight$X))
# gets the position of the effects that are neither controlling or conditioning (marginal variables)

mu_pred_cond<-stan_data_weight$X[,pred_cond_pos,drop=FALSE]%*%pars(model_z, "beta")[pred_cond_pos,1]
mu_pred_cont<-stan_data_weight$X[,pred_cont_pos,drop=FALSE]%*%pars(model_z, "beta")[pred_cont_pos,1]
mu_pred<-stan_data_weight$X[,pred_pos]%*%pars(model_z, "beta")[pred_pos,1]
# gets the predictions for the trait models for conditioning, controlling and marginal variables. 

nestD <- pars_ST(model_z, "nest")
residD <- pars_ST(model_z, if(trait=="weight_g"){"E"}else{"ind"})
geneD <- pars(model_z,"_A")[1]
# gets the distributional parameters for the random effects 

fixedD<-doppelgangR::st.mle(y = mu_pred)$dp
names(fixedD)<-paste0(c("xi", "omega", "alpha", "nu"), "_fixed")
# approximates the distribution of predictions for marginal variables as skew-t

z_st<-list(
  n_st=nestD[,1],
  e_st=residD[,1],
  fixed_st=fixedD,
  g_st = c(0, geneD, 0, 1e+16)
)
# list of distribution parameters for the different components: xi, omega, alpha, nu

pmu<-post_mu(model_z, components=c("nest", if(trait=="weight_g"){"E"}else{"ind"}, "A"), X=stan_data_weight$X, pred_pos=pred_pos, standardise=TRUE)
# posterior distribution of moments

mup<-dp2cm(z_st, family="ST", standardise=TRUE)
# moments given posterior mean of parameters

smu<-scm(THBW_egg_noRep[[trait]]-mu_pred_cond-mu_pred_cont, standardise=TRUE)
# sample moments of deviation of trait from conditioning+controlling predictions

par(mfrow=c(2,2))
for(i in 1:4){
 hist(pmu[which(pmu[,i]<100),i], breaks=100, main=paste("Central Moment", i), xlab="")
 abline(v=smu[i], col="red")
 abline(v=mup[i], col="blue")
}
# posterior distribution of central moments with observed values in red

if(is.null(cond_term)){
	THBW$category <- as.factor(rep(1, nrow(THBW)))
	THBW_egg_noRep$category <- as.factor(rep(1, nrow(THBW_egg_noRep)))
}else{
	THBW$category <- as.factor(apply(THBW[,cond_term], 1, paste, collapse=""))
	THBW_egg_noRep$category <- as.factor(apply(THBW_egg_noRep[,cond_term], 1, paste, collapse=""))
}
n_comb<-nlevels(THBW$category)
# add a category column to the data.frames used in the fitness model (THBW) and trait model (THBW_egg_noRep)
# that designate data into the (n_comb) categories based on the conditioning variables

beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model_w$Sol))
# linear coefficients for trait values from the survival model

gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model_w$Sol))
# quadratic coefficients for trait values from the survival model

eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model_w$Sol))
eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))
# positions of non-trait terms in the survival model that pertain to Day15 to fledging. 

eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model_w$Sol))
eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))
# positions of non-trait terms in the survival model that pertain to fledging to recruitment. 

fixedEffects <- formula(~ sex+male_present + year+hatch_day +clutch_sizeC + nest_hatch_dateC)
# single trait version of the fixed effect formula used in the survival model

X_eta1<-model.matrix(fixedEffects, THBW)
X_eta2<-model.matrix(fixedEffects, THBW)
# design matrix for non-trait effects in survival model

if(any(!mapply(colnames(X_eta1), colnames(model_w$Sol)[eta1_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta1")}

if(any(!mapply(colnames(X_eta2), colnames(model_w$Sol)[eta2_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta2")}


z<-THBW[,paste0(trait, "C")]
# trait value in the survival model (not that it was globally mean centred)

if(posterior_mean){
	n_it<-1
}else{
   	n_it<-nrow(model_w$Sol)
}
# number of MCMC iterations to evaluate

if(re_run){

	beta1 <- matrix(NA, n_it, n_comb)  # linear term in linear best fit
	beta2 <- matrix(NA, n_it, n_comb)  # linear term in quadratic best fit
	beta3 <- matrix(NA, n_it, n_comb)  # average gradient
	h2a <- matrix(NA, n_it, n_comb)    # heritability after selection

	int_opt<- matrix(NA, n_it, n_comb)

    Wplot.points<-seq(min(z), max(z), length=100)
    # (mean-centred) trait values at which to evaluate the fitness function for visualisation purposes

    Wplot_tmp<-matrix(NA, nlevels(THBW$category), 100)   # fitness functions for each category at an iteration

    Wplot<-matrix(NA, nrow(model_w$Sol), 100)            # fitness function averaged over categories (i.e colMeans on Wplot_tmp) 

    if(posterior_mean){

        model_w$Sol[1,]<-colMeans(model_w$Sol)
        model_w$VCV[1,]<-colMeans(model_w$VCV)
        # replace first iteration with posterior means if not iterating over complete posterior
    }
       
	for(i in 1:n_it){


		eta1 <- X_eta1%*%model_w$Sol[i,eta1_pos]
		eta2 <- X_eta2%*%model_w$Sol[i,eta2_pos]
        # non-trait fixed effects predictions from survival model

	    beta<-model_w$Sol[i,beta_pos]
	    gamma<-model_w$Sol[i,gamma_pos]
        # trait coefficients from survival model

	    V_nest<-matrix(model_w$VCV[i,grep("nest", colnames(model_w$VCV))], 2,2)
	     # nest covariance matrix from survival model
	   

	    for(j in 1:n_comb){

	    	cat<-levels(THBW$category)[j]
	    	# name of the category defined by the conditioning variables 

	    	eta1_sub<-eta1[THBW$category==cat]
	    	eta2_sub<-eta2[THBW$category==cat]
	    	z_sub<-z[THBW$category==cat]
	    	# get (mean-centred) trait and linear predictors for the right category

		    mu_etaz<-c(mean(eta1_sub), mean(eta2_sub), mean(z_sub)) 
		    V_etaz<-cov(cbind(eta1_sub, eta2_sub, z_sub))
            # mean and covariance matrices of linear predictor and traits

		    W<-sapply(z_sub, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
		    # fitness evaluated for each (mean-centred) trait value

		    WD<-sapply(z_sub, wD_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
            # derivative of fitness with respect to  (mean-centred) trait value evaluated for each (mean-centred) trait value

	        WDmin<-wD_func(min(z), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
	        WDmax<-wD_func(max(z), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
            # slope of the fitness function at extreme (mean-centered) trait values

	        Wplot_tmp[j,]<-sapply(Wplot.points, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
 			# fitness evaluated for (mean-centred) trait values used for visualisation of the fitness function

			if(any(THBW_egg_noRep$category==cat) & model_moments){
			  mu<-dp2cm(z_st, family="ST")
			  # central moments of the deviation of the (non-mean-centred) trait values 
			  # from the conditional/controlling predictions as inferred from the trait model

			  mu[1]<-mu[1]+mean(mu_pred_cond[THBW_egg_noRep$category==cat])
			  # mean adjusted for conditional prediction for category j if that category was in the trait model 
			}else{
			  mu<-scm(z_sub)
			  # If category j was not in the trait model (i.e. 2010), or model_moments=FALSE use the observed moments
			}  

	        S<-mean(W*z_sub/mean(W))-mu[1]
	        # selection differential calculated conditional on observed (mean-centered) trait value

	        C<-mean(W*((z_sub-mu[1])^2)/mean(W))-mu[1]
            # quadratic selection differential calculated conditional on observed (mean-centered) trait value

		    beta1[i,j]<-S/mu[2]  
		    beta2[i,j]<-betaLA_2(mu, S=S, C=C, family="ST")
		    beta3[i,j]<-mean(WD)/mean(W)
			h2a[i,j]<-h2a(10000, g_st=z_st$g_st, e_st=z_st[-which(names(z_st)=="g_st")], adj_mean=mu[1]-attr(z, "scaled:center"), mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)

		    int_opt[i,j]<-(WDmax<0 & WDmin>0)  # is there an internal stationary point
		}

		Wplot[i,]<-colMeans(Wplot_tmp)  # average fitness function

		print(i)
	}

	THBW$traitCat <- bin(z,n=20)+attr(z, "scaled:center")

	par(mar=c(5,5,1,1),mfrow=c(3,1), bty="l")
	hist(z_sub+attr(z, "scaled:center"), col="grey",main="", xlab=paste(trait), breaks=30)

	binPlot(recruit~traitCat, subset(THBW, category==cat), xlab=paste(trait), ylab="Fitness", ylim=c(0,0.1), xlim=range(z)+attr(z, "scaled:center"),text.cex=0.7);
	abline(v=mean(z_sub)+attr(z, "scaled:center"), col="grey")
	lines(W[order(z_sub)]~I(z_sub[order(z_sub)]+attr(z, "scaled:center")), col="red", lwd=2)

	plot(WD[order(z_sub)]~I(z_sub[order(z_sub)]+attr(z, "scaled:center")), type="l", ylab="Fitness Derivative", xlab=paste(trait), lwd=2, col="red")
	abline(v=mean(z_sub)+attr(z, "scaled:center"), col="grey")

	if(save){
		save(beta1,beta2, beta3, int_opt, Wplot, Wplot.points, h2a, file=paste0(wd,"Data/Intermediate/selection_gradient_",if(is.null(cond_term)){""}else{"by_"}, trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
	}

}else{

	files<-list.files(paste0(wd,"Data/Intermediate/"))
	files<-files[grep(paste0("selection_gradient_",if(is.null(cond_term)){""}else{"by_"}, trait), files)]
	if(length(files)>1){warning("Multiple selection_gradient files, using the last one")}
	load(paste0(wd,"Data/Intermediate/", files[length(files)]))
}

if(save_plot){
  pdf(paste0(wd, "Tex/W_", trait, ".pdf"))
}



par(mar=c(5,5,1,1),mfrow=c(2,1))

hist(THBW_egg_noRep[[trait]], col="grey",main="", xlab=paste(trait), breaks=30)


THBW_egg_noRep$binned<-bin(THBW_egg_noRep[[trait]])
THBW_egg_noRep$w<-THBW$recruit[match(THBW_egg_noRep$bird_id, THBW$bird_id)]

post_it<-seq(1, nrow(model_w$Sol), 15)

plot(colMeans(Wplot)~I(Wplot.points+mean(THBW_egg_noRep[[trait]])), lwd=2, xlim=range(THBW_egg_noRep[[trait]]),  ylim=c(0,max(Wplot[post_it,])), xlab=paste(trait), ylab="Fitness", type="l")


for(i in post_it){

	lines(Wplot[i,]~I(Wplot.points+mean(THBW_egg_noRep[[trait]])), col="red", lwd=0.01)
}

lines(colMeans(Wplot)~I(Wplot.points+mean(THBW_egg_noRep[[trait]])), lwd=2)

binPoints(w~binned, THBW_egg_noRep)

if(save_plot){
dev.off()
}

beta1<-rowMeans(beta1)
beta2<-rowMeans(beta2)
beta3<-rowMeans(beta3)

par(mfrow=c(3,1))
hist(beta3, xlim=range(c(beta1,beta3)), breaks=50)
hist(beta1, xlim=range(c(beta1,beta3)), breaks=50)
hist(beta1-beta3, breaks=50)

par(mfrow=c(3,1))
hist(beta3, xlim=range(c(beta2,beta3)), breaks=50)
hist(beta2, xlim=range(c(beta2,beta3)), breaks=50)
hist(beta2-beta3, breaks=50)




