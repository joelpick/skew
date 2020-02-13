rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
library(MCMCglmm)
# library(MasterBayes)
library(sn)
library(mvtnorm)


if(Sys.info()["user"]=="jhadfiel"){
  wd <- "..."
}else{
  wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))
THBW_egg<- subset(THBW_egg,!is.na(sex) & !year%in%c(2018,2019)& !is.na(tarsus_mm)& !is.na(headbill_mm)& !is.na(wing_mm)& !is.na(weight_g))

load(paste0(wd,"Data/Intermediate/stanModWeightPedN20190813_0021.Rdata"))
load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20190812_2226.Rdata"))


pars_ST <- function(model, variable){
  summary(model)$summary[paste(c("xi","omega", "alpha", "nu"),variable, sep="_"),c(1,4,8)]
}

pars <- function(model, par){
  out <- summary(model)$summary
  out[grep(par, rownames(out)),c(1,4,8)]
}

nestT <- pars_ST(mod_stan_tarsus_pedN, "nest")
residT <- pars_ST(mod_stan_tarsus_pedN, "ind")
nestW <- pars_ST(mod_stan_weight_pedN, "nest")
residW <- pars_ST(mod_stan_weight_pedN, "E")

muT <- pars(mod_stan_tarsus_pedN,"beta")[1]
muW <- pars(mod_stan_weight_pedN,"beta")[1]
sigmaA_T <- pars(mod_stan_tarsus_pedN,"_A")[1]
sigmaA_W <- pars(mod_stan_weight_pedN,"_A")[1]

range(THBW_egg$weight_g)

dat_limT <-range(THBW_egg$tarsus_mm)- muT 
dat_limW <-range(THBW_egg$weight_g)- muW 

dat_lim <- dat_limT

e_st<-list(
  n_st=nestT[,1],
  e_st=residT[,1]
)
# list of environmental distribution parameters: xi, omega, alpha, nu

g_st=c(0, sigmaA_T,  0, Inf)
# list of genetic distribution parameters: xi, omega, alpha, nu
int = muT
n<-10000
# number of simulated parent-offspring pairs

g_p = rst(n=n, dp=g_st)
e_p = lapply(e_st, function(x){rst(n=n, dp=x)})
# simulate genetic and environmental effects for parents

g_o = g_p
e_o = lapply(e_st, function(x){rst(n=n, dp=x)})
# simulate genetic and environmental effects for offspring
# offspring genotype assumed equal to parental genotype: haploid model.

z_p = g_p+rowSums(as.data.frame(e_p))
z_o = g_o+rowSums(as.data.frame(e_o))
# obtain phenotypes by summing

conv<-function(par, z_p, g_st, e_st){
  if(length(par)!=length(e_st)){stop("par should be of length 2")}
  res<-dst(z_p-sum(par), dp=g_st)
  for(i in 1:length(e_st)){
     res<-res*dst(par[i], dp=e_st[[i]])
  }
  return(res)
}  
# Function for performing the convolution 

wconv<-function(par, z_p, g_st, e_st){

  (z_p-sum(par))*conv(par, z_p, g_st, e_st)
}
# Function for performing the convolution weighted by genotype effect

Ee <- sum(unlist(lapply(e_st, function(x){st.cumulants(dp=x, n=1)})))
# Expectation of the environmental component

e_llimits<-unlist(lapply(e_st, function(x){qst(1e-4, dp=x)}))
e_ulimits<-unlist(lapply(e_st, function(x){qst(1-1e-4, dp=x)}))
# integration limits

#x<-seq(quantile(z_p, 0.025), quantile(z_p, 1-0.025), length=100)
x<- seq(dat_lim[1],dat_lim[2],length=100)
# parental phenotypes (z_p) over which to calculate E[g_p|z_p]

Eg_p<-1:100
for(i in 1:100){
  Eg_p[i] <- hcubature(wconv, e_llimits, e_ulimits, z_p=x[i], g_st=g_st, e_st=e_st)$integral
  Eg_p[i] <- Eg_p[i]/hcubature(conv, e_llimits, e_ulimits, z_p=x[i], g_st=g_st, e_st=e_st)$integral
} 
# E[g_p|z_p]

cz_p<-cut(z_p, quantile(z_p, seq(0,1,length=20)))
# discretise z_p into 20 categories

mz_o<-tapply(z_o, cz_p, mean)
mz_p<-tapply(z_p, cz_p, mean)
# calculate expected parent and offspring phenotypes for each category.


# parent_T <- x + int
# offspring_T <- Eg_p+Ee + int
# parent_W <- x + int
# offspring_W <- Eg_p+Ee + int 
# save(parent_W,offspring_W,parent_T,offspring_T, file=paste0(wd,"/skew/Data/Intermediate/PO_regression.Rdata"))

par(mar=c(5,5,1,1),mfrow=c(2,1))
hist( THBW$tarsus_mm, col="grey",main="", xlab="Weight (g)")
#plot(I(mz_o + int)~I(mz_p + int), xlim=range(dat_limW  + int), ylim=range(dat_limW + int))
plot(I(Eg_p+Ee + int)~I(x + int), xlim=range(dat_lim  + int), ylim=range(dat_lim + int), type="l", xlab="Parent", ylab="Offspring", lwd=2)

