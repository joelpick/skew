rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
library(MCMCglmm)
library(sn)
library(mvtnorm)
library(cubature)
library(gam)
library(doppelgangR)

trait<-"weight_g"
simulate<-FALSE
nsim<-10000 # number of simulated parent-offspring pairs
pred_cond<-c("year2012", "year2013", "year2014", "year2015","year2016","year2017","year2018","sexM", "timeC") # terms to omit from the variance


if(simulate & nsim%%200!=0){
 stop("nsim must be divisible by 200")
}

if(Sys.info()["user"]=="jhadfiel"){
  wd <- "~/Work/Skew/"
}else{
  wd <- "~/Dropbox/0_blue_tits/skew/"
}

pars_ST <- function(model, variable){
  summary(model)$summary[paste(c("xi","omega", "alpha", "nu"),variable, sep="_"),c(1,4,8)]
}

pars <- function(model, par){
  out <- summary(model)$summary
  out[grep(par, rownames(out)),c(1,4,8)]
}

source(paste0(wd,"R/functions.R"))

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))

if(trait=="tarsus_mm"){
  model<-get(load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20191217_1716.Rdata")))
}
if(trait=="weight_g"){  
  model<-get(load(paste0(wd,"Data/Intermediate/stanModWeightPedN20191220_0905.Rdata")))
}

pred_cond_pos<-match(pred_cond, colnames(stan_data_weight$X))
pred_pos<-match(setdiff(colnames(stan_data_weight$X), pred_cond), colnames(stan_data_weight$X))

mu_pred_cond<-stan_data_weight$X[,pred_cond_pos]%*%pars(model, "beta")[pred_cond_pos,1]
mu_pred<-stan_data_weight$X[,pred_pos]%*%pars(model, "beta")[pred_pos,1]

nestD <- pars_ST(model, "nest")
residD <- pars_ST(model, if(trait=="weight_g"){"E"}else{"ind"})
geneD <- pars(model,"_A")[1]
fixedD<-st.mle(y = mu_pred)$dp
names(fixedD)<-paste0(c("xi", "omega", "alpha", "nu"), "_fixed")

mc<-max(hist(mu_pred, breaks=30)$counts)
dmu_pred<-dst(mu_pred, dp=fixedD)
lines(I(dmu_pred[order(mu_pred)]*mc/max(dmu_pred))~mu_pred[order(mu_pred)])
# check to make sure skew-t approximation for fixed effect predictors is good.

e_st<-list(
  n_st=nestD[,1],
  e_st=residD[,1],
  fixed_st=fixedD
)
# list of environmental distribution parameters: xi, omega, alpha, nu

g_st=c(0, geneD, 0, Inf)
# list of genetic distribution parameters: xi, omega, alpha, nu

Vp<-sum(unlist(lapply(e_st, function(x){dp2cp(x, family="ST")["s.d."]^2})))+geneD^2

h2<-(geneD^2)/Vp



conv<-function(par, z_p, g_st, e_st){

  ################################################################
  #   Function for obtaining the integrand in the convolution    #
  #             \int Pr(z_p-g_p)Pr(g_p) dg_p                     #
  #    the distribution of genetic and environmental effects     #
  #          are passed as parameters of the skew-t              #
  ################################################################

  if(length(par)!=length(e_st)){stop("par should be of length 2")}

  res<-dst(z_p-sum(par), dp=g_st)
  # Pr(z_p-e_p | \theta_g) = Pr(g_p | \theta_g)

  for(i in 1:length(e_st)){
     res<-res*dst(par[i], dp=e_st[[i]])
  }

  # Pr(g_p | \theta_g)*Pr(e_p| \theta_e) = Pr(g_p | \theta_g)*Pr(z_p-g_p| \theta_e)

  return(res)
}  


wconv<-function(par, z_p, g_st, e_st){

  #########################################################################
  #   Function for obtaining the integrand in the weighted convolution    #
  #                 \int g_p*Pr(z_p-g_p)Pr(g_p) dg_p                      #
  #      the distribution of genetic and environmental effects            #
  #          are passed as parameters of the skew-t                       #
  #########################################################################

  (z_p-sum(par))*conv(par, z_p, g_st, e_st)

}


POreg<-function(z_p, g_st, e_st, limit_prob=1e-4){

  ################################################################
  #           Function for calculating E[g_p|z_p]                #
  #  based on distribution of genetic and environmental effects  #
  #          (passed as parameters of the skew-t)                #
  ################################################################

  e_llimits<-unlist(lapply(e_st, function(x){qst(limit_prob, dp=x)}))
  e_ulimits<-unlist(lapply(e_st, function(x){qst(1-limit_prob, dp=x)}))
  # lower and upper limits when integrating over environmental effects

  Eg_p<-hcubature(wconv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
  Eg_p<-Eg_p/hcubature(conv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral

  return(Eg_p)
}


if(simulate){

  n<-10000

  g_p = rst(n=nsim, dp=g_st)
  e_p = lapply(e_st, function(x){rst(n=nsim, dp=x)})
  # simulate genetic and environmental effects for parents

  g_o = 0.5*g_p+0.5*rst(n=nsim, dp=g_st)+rnorm(nsim, 0, geneD/sqrt(2))
  e_o = lapply(e_st, function(x){rst(n=nsim, dp=x)})
  # simulate genetic and environmental effects for offspring
  # offspring genotype assumed equal to parental genotype: haploid model.

  z_p = g_p+rowSums(as.data.frame(e_p))
  z_o = g_o+rowSums(as.data.frame(e_o))
  z_n = rep(1, nsim)
  # obtain phenotypes by summing

}else{

  z_om<-tapply(THBW_egg_noRep[[trait]]-mu_pred_cond, THBW_egg_noRep$dam_P, mean, na.rm=TRUE)
  z_nm<-tapply(THBW_egg_noRep[[trait]], THBW_egg_noRep$dam_P, function(x){sum(!is.na(x))})
  z_m<-(THBW_egg_noRep[[trait]]-mu_pred_cond)[match(names(z_om), THBW_egg_noRep$bird_id)]

  z_of<-tapply(THBW_egg_noRep[[trait]]-mu_pred_cond, THBW_egg_noRep$sire_P, mean, na.rm=TRUE)
  z_nf<-tapply(THBW_egg_noRep[[trait]], THBW_egg_noRep$sire_P, function(x){sum(!is.na(x))})
  z_f<-(THBW_egg_noRep[[trait]]-mu_pred_cond)[match(names(z_of), THBW_egg_noRep$bird_id)]

  z_o<-c(z_om, z_of)
  z_p<-c(z_m, z_f)
  z_n<-c(z_nm, z_nf)

  complete<-which(!is.na(z_p) & !is.na(z_o))

  z_o<-z_o[complete]
  z_p<-z_p[complete]
  z_n<-z_n[complete]

}

npoints<-min(length(z_p),200)

pred_points<-sort(z_p)[seq(1, length(z_p), length(z_p)/npoints)]

Eg_p<-1:npoints

for(i in 1:npoints){
  Eg_p[i] <- POreg(z_p=pred_points[i], g_st=g_st, e_st=e_st)
} 
# E[g_p|z_p]

pdf(paste0(wd, "Tex/PO_", trait, ".pdf"))
par(mar=c(5,5,1,1),mfrow=c(2,1))
hist(THBW_egg_noRep[[trait]]-mu_pred_cond, col="grey",main="", xlab=paste(trait), breaks=30)

plot(I(0.5*Eg_p+mskt(dp=e_st$fixed_st))~pred_points, xlim=range(THBW_egg_noRep[[trait]]-mu_pred_cond), ylim=range(z_o), type="l", xlab="Parent", ylab="Offspring", lwd=2, col="red")

points(z_o~z_p, col=c("black", "grey")[(z_p==min(z_p))+1])

mt<-lm(z_o~offset((0.5*Eg_p+dp2cp(e_st$fixed_st, "ST")["mean"])[match(z_p, pred_points)]), weights=1/sqrt(z_n), subset=!is.na(z_p) & z_p>min(z_p))
m1<-lm(z_o~z_p, weights=1/sqrt(z_n), subset=!is.na(z_p) & z_p>min(z_p))           # & y_p>15 

z_pred<-predict(m1, newdata=list(z_p=pred_points), interval="confidence")

lines(z_pred[,1]~pred_points, lty=1)
lines(z_pred[,2]~pred_points, lty=2)
lines(z_pred[,3]~pred_points, lty=2)
dev.off()




