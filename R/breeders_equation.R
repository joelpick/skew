rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
# source(paste0(wd,"functions.R"))

data_wd <- paste0(wd,"Data/Intermediate/")
plot_wd <- "~/Dropbox/0_Presentations/images/"


load(paste0(data_wd,"chick_data.Rdata"))

load(paste0(data_wd,"dam_sire_egg.Rdata"))


tBIRDS <- read.csv(paste0(wd,"Data/Raw/tBIRDS.csv"))
head(tBIRDS)
adults <- unique(subset(read.csv(paste0(wd,"Data/Raw/tMORPH.csv")), agecode%in%c(5,6))$bird_id)

#remove repeated measures and those without a sex
THBW<- subset(THBW,!duplicated(bird_id) & !is.na(sex) & !is.na(tarsus_mm) & !is.na(headbill_mm) & !is.na(wing_mm) & !is.na(weight_g))

THBW$recruit <- THBW$bird_id %in% adults
THBW$fledge <- is.na(tBIRDS$max_death_date[match(THBW$bird_id , tBIRDS$bird_id)])

S_T <- mean(THBW$tarsus_mm[THBW$recruit]) - mean(THBW$tarsus_mm)

h2_T <- modT$vars_mean["mod_egg_A","animal"]/sum(c(modT$vars_mean["mod_egg_A",],modT$fixed_var["mod_egg_A"]))
#var(THBW$tarsus_mm)

S_W <- mean(THBW$weight_g[THBW$recruit]) - mean(THBW$weight_g)

h2_W <- modM$vars_mean["mod_egg_A","animal"]/sum(c(modM$vars_mean["mod_egg_A",],modM$fixed_var["mod_egg_A"]))


S*h2


skewModT

load(file=paste0(data_wd,"stan_summary_data.Rdata"))

load(file=paste0(data_wd,"selection_gradient_tarsus_mmC_20200110_1423.Rdata"))
S_normal_T <- S_normal
S_skew_T <- S_skew

load(file=paste0(data_wd,"selection_gradient_weight_gC_20200110_1433.Rdata"))
S_normal_W <- S_normal
S_skew_W <- S_skew

mean_CI <- function(x) c(mean(x),quantile(x, c(0.025,0.975)))

Tarsus <- rbind(Breeders = c(S_T*h2_T,NA,NA),rbind(Skew=mean_CI(S_skew_T),Normal=mean_CI(S_normal_T),Difference=mean_CI(S_skew_T-S_normal_T))* skewModPed_T$var["A",1])

Weight <- rbind(Breeders = c(S_W*h2_W,NA,NA),rbind(Skew=mean_CI(S_skew_W),Normal=mean_CI(S_normal_W),Difference=mean_CI(S_skew_W-S_normal_W)) * 0.2127457)

blankPlot <- function(ylim=c(-1,1),xlim=c(-1,1)) {
	op <- par(mar=c(0,0,0,0))
	plot(NULL,ylim=ylim,xlim=xlim, xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	par(op)
}
effectPlot <- function(x, y=NULL, xlab="", ylab= "", xlim=NULL, col="black", cex.axis=1,add=FALSE,...){ #want a matrix of mean, 0.025 and 0.975
	if(is.null(y)) y <- nrow(x):1
	if(is.null(xlim)){
		xlim<- if(min(x)>0 & max(x)>0) {c(0,max(x))
		}else if(min(x)<0&max(x)<0) {c(min(x),0)
		}else{	range(x)}}
	if(add){	
		points(x[,1], y, pch=19, col=col, cex=1.5)
	}else{
		plot(x[,1], y, xlab=xlab, ylab=ylab, xlim=xlim, pch=19, yaxt="n", ylim=c(0.5,nrow(x)+0.5),,col=col, cex.axis=cex.axis, cex=1.5,...)#c(min(y)-0.5,max(x)+0.5)
		axis(2,y,rownames(x),cex.axis=cex.axis*1.5, las=1)
		abline(v=0, col="grey")
	}
	arrows(x[,2], y, x[,3], y, code=3, angle=180, length=0.01,col=col, lwd=1.5)
}
cex.axis=1.5

#setEPS()
#pdf(paste0(plot_wd,"skew_day15_results6.pdf"), height=8, width=5)
{
layout(matrix(1:3,ncol=1,byrow=TRUE), height=c(2,rep(4,2)))
par(mar=c(0,0,0,0))
blankPlot(); text(0,0,"Selection\nResponse",cex=3)
par(mar=c(3,7,0,1),cex.axis=cex.axis, cex.lab=2)
effectPlot(Tarsus, xlim=c(-0.01,0.1))
effectPlot(Weight, xlim=c(-0.01,0.1), col=2)
}
#dev.off()
