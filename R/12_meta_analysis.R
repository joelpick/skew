rm(list=ls())

options(stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/github/skew/"
}

source(paste0(wd,"R/00_functions.R"))

all<- read.csv(paste0(wd,"Data/Intermediate/MA_data.csv"))

library(scales)
library(MCMCglmm)

pMCMC <- function(x) 2*pmax(0.5/length(x), pmin(sum(x > 0)/length(x), 1 - sum(x > 0)/length(x)))



# Random effects meta-analysis in MCMCglmm.
# Population and species were modelled as random effects.
# Sex and age and their interaction were included as fixed effects.

## sampling error
# all$SE <-sqrt(6*all$n*(all$n-1)/((all$n-2)*(all$n+1)*(all$n+3)))
all$SE <- with(all,sqrt( ((n^2)/((n-2)*(n-1)))^3 *6*(n-2)/((n+1)*(n+3))))



prior <- list(R = list(V = diag(1), nu = 0.002),
               G = list(
               	 G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = 1000),
                 G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = 1000),
                 G3 = list(V = diag(1), fix = 1)))
a <- 5
MAmod<-MCMCglmm(skew~age-1,#all$
	random = ~ species + study_ref + idh(SE):units,
	data=all,
	prior = prior,
	nitt=13000*a,thin=10*a,burnin=3000*a)
summary(MAmod)
plot(MAmod)

save(MAmod, file= paste0(wd,"Data/Intermediate/meta_analysis.Rdata"), version=2)



setEPS()
pdf(paste0(wd,"Plots/figure2.pdf"), height=5, width=4)

{
	par(mar=c(3,4,1,2), cex.axis=1.1, cex.lab=1.25, mgp=c(2.5,1,0))

xA <- all[all$age=="A","skew"]
xJ <- all[all$age=="J","skew"]
# vioplot(xJ, xA, names=c("Juveniles","Adults"), col="grey", ylim=c(-7,1.7))
# mtext("Skew",side=2, line=2.5, cex=1.25)

boxplot(xJ, xA, names=c("Juveniles","Adults"), ylim=c(-5,1.7), ylab="Skew",xlab="", col="grey",pch=19)

abline(h=0)

text(c(2,1),c(1.5,1.5),paste0(table(all$age)," (",tapply(all$species,all$age,function(x)length(unique(x))),")"))
points(c(2,1),all$skew[all$study_ref==43], pch=19,col="red", cex=1.2)

}

dev.off()

