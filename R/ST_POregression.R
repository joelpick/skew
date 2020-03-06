rm(list=ls())

# source("~/Work/Skew/R/ST_POregression.R")

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
library(MCMCglmm)
library(sn)
library(mvtnorm)
library(cubature)

trait<-"weight_g"
re_run=FALSE
simulate<-FALSE
nsim<-10000 # number of simulated parent-offspring pairs
save<-TRUE
save_plot<-TRUE
pred_cond<-c("year2012", "year2013", "year2014", "year2015","year2016","year2017","year2018","sexM", "timeC") # terms to omit from the variance


if(simulate & nsim%%200!=0){
 stop("nsim must be divisible by 200")
}

if(Sys.info()["user"]=="jhadfiel"){
  wd <- "~/Work/Skew/"
}else{
  wd <- "~/Dropbox/0_blue_tits/skew/"
}

post_st<-function(model, it){

}


source(paste0(wd,"R/functions.R"))

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))

if(trait=="tarsus_mm"){
  model_z<-get(load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20191217_1716.Rdata")))
}
if(trait=="weight_g"){  
  model_z<-get(load(paste0(wd,"Data/Intermediate/stanModWeightPedN20191220_0905.Rdata")))
}

pred_cond_pos<-match(pred_cond, colnames(stan_data_weight$X))
pred_pos<-match(setdiff(colnames(stan_data_weight$X), pred_cond), colnames(stan_data_weight$X))

mu_pred_cond<-stan_data_weight$X[,pred_cond_pos]%*%pars(model_z, "beta")[pred_cond_pos,1]
mu_pred<-stan_data_weight$X[,pred_pos]%*%pars(model_z, "beta")[pred_pos,1]



nestD <- pars_ST(model_z, "nest")
residD <- pars_ST(model_z, if(trait=="weight_g"){"E"}else{"ind"})
geneD <- pars(model_z,"_A")[1]
fixedD<-doppelgangR::st.mle(y = mu_pred)$dp
names(fixedD)<-paste0(c("xi", "omega", "alpha", "nu"), "_fixed")


comp_skt(mu_pred, dp=fixedD, breaks=30)
# check to make sure skew-t approximation for fixed effect predictors is good.

e_st<-list(
  n_st=nestD[,1],
  e_st=residD[,1],
  fixed_st=fixedD
)
# list of distribution parameters: xi, omega, alpha, nu

g_st<-c(0, geneD, 0, 1e+16)



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

if(re_run){

  pred_points<-sort(z_p)[seq(1, length(z_p), length(z_p)/npoints)]

  Eg_p<-dEg_p<-dz_p<-1:npoints

  for(i in 1:npoints){
    dz_p[i] <- dz(z_p=pred_points[i], g_st=g_st, e_st=e_st)
    Eg_p[i] <- POreg(z_p=pred_points[i], g_st=g_st, e_st=e_st)
    dEg_p[i] <- dPOreg(z_p=pred_points[i], g_st=g_st, e_st=e_st, Eg_p=Eg_p[i])
    print(i)
  } 

  dz_p_norm<-dnorm(pred_points, sum(unlist(lapply(e_st, function(x){dp2cp(x, family="ST")[1]}))), g_st[2]^2+sum(unlist(lapply(e_st, function(x){dp2cp(x, family="ST")[2]^2}))))

  if(save){
  save(Eg_p,dEg_p, dz_p, dz_p_norm, pred_points, file=paste0(wd,"Data/Intermediate/po_regression_",if(is.null(pred_cond)){""}else{"by_"}, trait,"_",format(Sys.time(), "%Y%m%d_%H%M"),".Rdata"))
  }
}else{

  files<-list.files(paste0(wd,"Data/Intermediate/"))
  files<-files[grep(paste0("po_regression_",if(is.null(pred_cond)){""}else{"by_"}, trait), files)]
  if(length(files)>1){warning("Multiple po_regressions files, using the last one")}
  load(paste0(wd,"Data/Intermediate/", files[length(files)]))

}
# load(paste0(wd,"Data/Intermediate/po_regression_by_weight_g_20200228_1649.Rdata"))
# load(paste0(wd,"Data/Intermediate/po_regression_by_tarsus_mm_20200228_1708.Rdata"))

mu_eg<-0.5*dp2cp(g_st, family="ST")[1]+sum(unlist(lapply(e_st, function(x){dp2cp(x, family="ST")[1]})))

mt<-lm(z_o~offset(I(Eg_p)[match(z_p, pred_points)]), weights=z_n)

mu_eg<-coef(mt) # not sure why these differ by so much!

if(save_plot){
  pdf(paste0(wd, "Tex/PO_", trait, ".pdf"))
}

par(mar=c(5,5,1,1),mfrow=c(2,1))
hist(THBW_egg_noRep[[trait]]-mu_pred_cond, col="grey",main="", xlab=paste(trait), breaks=30)

plot(I(Eg_p+mu_eg)~pred_points, xlim=range(THBW_egg_noRep[[trait]]-mu_pred_cond), ylim=range(z_o), type="l", xlab="Parent", ylab="Offspring", lwd=2, col="red")

#for(i in 1:npoints){
#   arrows(pred_points[i], (Eg_p+mu_eg)[i], pred_points[i]+0.1, (Eg_p+mu_eg)[i]+0.1*dEg_p[i], length=0.1)
#}
# check to make sure derivative correct

points(z_o~z_p, cex=0.3*sqrt(z_n))


m1<-lm(z_o~z_p, weights=z_n)
mt<-lm(z_o~offset(I(Eg_p+mu_eg)[match(z_p, pred_points)])-1, weights=z_n)          

z_pred<-predict(m1, newdata=list(z_p=pred_points), interval="confidence")

lines(z_pred[,1]~pred_points, lty=1)
lines(z_pred[,2]~pred_points, lty=2)
lines(z_pred[,3]~pred_points, lty=2)

if(save_plot){
dev.off()
}



