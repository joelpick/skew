rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

wd <- "~/Dropbox/0_blue_tits/skew"

source(paste0(wd,"/R/functions.R"))
load(file= paste0(wd,"/Data/Intermediate/dam_sire_egg.Rdata"))
load(file= paste0(wd,"/Data/Intermediate/dam_sire_egg_M.Rdata"))


plot_results1 <- function(mod,ylim,main,rows, cols){
	position <- c(-0.2,0,0.2)
	col<-c(1,2,4)
	axis_names=c("Fixed","A","Rear", "Origin","Resid","ME")
	n_var <- length(cols)
	mod_P_means<-cbind(fixed=mod$fixed_var,mod$vars_mean)[rows,cols]
	mod_P_se<-mod$vars_SE[rows,cols[2:n_var]-1]
	
	plot(NA, ylim=ylim, xlim=c(0.5,n_var+0.5), xaxt="n", xlab="Variance component", ylab="Variance", pch=19, col=1,main=main)
	axis(1,c(1:(n_var)),axis_names[cols])
	abline(h=0)

	for(i in 1:length(rows)){
		points(mod_P_means[i,]~I(c(1:n_var)+position[i]), pch=19, col=col[i])
		arrows((2:n_var +position[i]),mod_P_means[i,2:n_var]+mod_P_se[i,]*1.96,(2:n_var +position[i]),mod_P_means[i,2:n_var]-mod_P_se[i,]*1.96, code=3,angle=90,length=0.05)
	}
}


setEPS()
pdf(paste0(wd,"/R/plots/figure_SM_eggSize.pdf"), height=10, width=10)
{
par(mfrow=c(2,2))
plot_results1(modT, ylim=c(0,0.2), main="Tarsus",rows=c(2,4,6), cols=1:6)
legend("topright",c("NO","NO-egg","egg"), pch=19, col=c(1,2,4))
plot_results1(modHB, ylim=c(0,0.2), main="Head-Bill",rows=c(2,4,6), cols=1:6)
plot_results1(modM, ylim=c(0,0.8), main="Mass",rows=c(2,4,6), cols=1:5)
plot_results1(modW, ylim=c(0,7.5), main="Wing",rows=c(2,4,6), cols=1:6)
}
dev.off()

setEPS()
pdf(paste0(wd,"/R/plots/figure_SM_DS_Animal.pdf"), height=10, width=10)
{
par(mfrow=c(2,2))
plot_results1(modT, ylim=c(0,0.2), main="Tarsus",rows=c(5,6), cols=c(1:3,5:6))
legend("topright",c("Dam-Sire","Animal"), pch=19, col=c(1,2))
plot_results1(modHB, ylim=c(0,0.2), main="Head-Bill",rows=c(5,6), cols=c(1:3,5:6))
plot_results1(modM, ylim=c(0,0.8), main="Mass",rows=c(5,6), cols=c(1:3,5))
plot_results1(modW, ylim=c(0,7.5), main="Wing",rows=c(5,6), cols=c(1:3,5:6))
}
dev.off()

setEPS()
pdf(paste0(wd,"/R/plots/figure_SM_pedigree.pdf"), height=10, width=10)
{
par(mfrow=c(2,2))
plot_results1(modT_all, ylim=c(0,0.2), main="Tarsus",rows=c(11,12), cols=c(1:3,5:6))
legend("topright",c("P","M"), pch=19, col=c(1,2))
plot_results1(modHB_all, ylim=c(0,0.2), main="Head-Bill",rows=c(11,12), cols=c(1:3,5:6))
plot_results1(modM_all, ylim=c(0,0.8), main="Mass",rows=c(11,12), cols=c(1:3,5))
plot_results1(modW_all, ylim=c(0,7.5), main="Wing",rows=c(11,12), cols=c(1:3,5:6))
}
dev.off()
