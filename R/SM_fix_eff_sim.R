rm(list=ls())

options( stringsAsFactors=FALSE)

library(rstan)
library(sn)
library(scales)
library(pbapply)
library(lme4)

run<-FALSE

wd <- "~/Dropbox/0_blue_tits/skew"
data_wd <- paste0(wd,"/Data/Intermediate/")
source(paste0(wd,"/R/00_functions.R"))
load(paste0(data_wd,"chick_data.Rdata"))
load(paste0(data_wd,"stanModWing_DS20191210_1850.Rdata"))
load(paste0(data_wd,"stan_data.Rdata"))

simSTRE <- function(levels,dp){
  dMat <- model.matrix(~ as.character(levels)-1)
  N <- length(unique(levels)) 
  pMat <- rst(N,dp=dp)
  return(dMat %*% pMat)
}


fixed <- pars(mod_stan_wing_dam,"beta")
nest <- pars_ST(mod_stan_wing_dam, "nest")
dam_sire <- pars_ST(mod_stan_wing_dam, "dam_sire")
resid <- pars_ST(mod_stan_wing_dam, "ind")
me <- pars(mod_stan_wing_dam, "sigma_E")


if(run){
	out<-pbreplicate(1000,{
	
		y <- stan_data_DS$X%*%fixed[,1] + 
			 simSTRE(levels=THBW$nest,dp=nest[,1]) +
			 simSTRE(levels=THBW$dam_P,dp=dam_sire[,1]) + 
			 simSTRE(levels=THBW$sire_P,dp=dam_sire[,1]) + 
			 simSTRE(levels=THBW$bird_id,dp=resid[,1]) +
			 rnorm(nrow(THBW),0,me[1])
	
		THBW$wing_sim <- y

		res<-summary(lmer(wing_sim ~ malePresent + clutchSizeC + nestHatchDateC  + hatchDay + year + timeC + sex + eggWeightC + (1|nest) + (1|dam_P) + (1|sire_P) + (1|bird_id), data=THBW))$coef[,1]
		return(res)
	})


	save(out, file=paste0(data_wd,"lmer_sims.Rdata"))
}else{
	load(paste0(data_wd,"lmer_sims.Rdata"))	
}


mean_CI2 <- function(x){
	c(mean(x),mean(x)-se(x)*1.96,mean(x)+se(x)*1.96)
}

out[1,] <- out[1,]-35
fixed_plot <- fixed[,1]
fixed_plot[1] <- fixed_plot[1]-35
setEPS()
pdf(paste0(wd,"/R/plots/figure_SM_lmer_sim.pdf"), height=10, width=10)
{
par(mfrow=c(1,1),mar=c(4,12,1,1))
	effectPlot(t(apply(out,1,mean_CI2)), xlab="Estimate", cex.lab=1.75, arrow.length=0.2, xaxt="n")
arrows(fixed_plot,(17:1)-0.3,fixed_plot,(17:1)+0.3,code=0,col="red",lwd=3)
axis(1,seq(-5,8,1),c(seq(-5,5,1),40,41,42))
axis.break(1,5.5,style="slash")

}
dev.off()

