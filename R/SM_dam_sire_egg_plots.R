rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

#install.packages("rstan")
# library(MCMCglmm)
# library(MasterBayes)


wd <- "~/Dropbox/0_blue_tits/skew"

source(paste0(wd,"/R/functions.R"))
load(file= paste0(wd,"/Data/Intermediate/dam_sire_egg.Rdata"))
load(file= paste0(wd,"/Data/Intermediate/dam_sire_egg_M.Rdata"))

# plot(NA,ylim=c(-0.02,0.2), xlim=c(1,6), xaxt="n", xlab="Variance component", ylab="Variance")
# abline(h=0)
# axis(1,c(1:5)+0.25,colnames(modT$vars_mean))
# cols <- rep(c(1,2),3)
# pchs <- rep(15:17,each=2)

# for(i in 1:6){
# 	x <- (1:5) +(0.1*i-0.1)
# 	arrows(x,modT$vars_mean[i,]+modT$vars_SE[i,],x,modT$vars_mean[i,]-modT$vars_SE[i,], code=3,angle=90,length=0.05)
# 	points(modT$vars_mean[i,]~x,pch=pchs[i], col=cols[i])
# }
# n_mod=12
# 1+(1:n_mod)/10

# plot_results <- function(mod,ylim,n_var,axis_names=NULL,main){
# 	n_mod <- nrow(mod$vars_mean)
# 	cols <- rep(c(1,2),n_mod/2)
# 	pchs <- #rep(15:17,each=n_mod/3)
# 	 rep(c(15,22,16,21,17,24),each=n_mod/6)
# 	x1 <- 1+(1:n_mod)/20
# 	if(is.null(axis_names)) axis_names <- c("fixed",colnames(mod$vars_mean))
# 	plot(mod$fixed_var~x1,ylim=ylim, xlim=c(1,n_var+2), xaxt="n", xlab="Variance component", ylab="Variance",col=cols, pch=pchs,main=main)
# 	abline(h=0)
# 	axis(1,c(1:(n_var+1))+0.25,axis_names)

# 	for(i in 1:n_mod){
# 		x <- (1:n_var) +(0.075*i-0.075) +1
# 		# x <- x1 + i
# 		arrows(x,mod$vars_mean[i,]+mod$vars_SE[i,],x,mod$vars_mean[i,]-mod$vars_SE[i,], code=3,angle=90,length=0.05)
# 		points(mod$vars_mean[i,]~x,pch=pchs[i], col=cols[i],bg="white")
# 	}
# }

add_mod <- function(mod1,mod2){
	mod_all <- list()
	for(i in 1:6) mod_all$vars_mean <- rbind(mod_all$vars_mean , mod1$vars_mean[i,],mod2$vars_mean[i,])
	for(i in 1:6) mod_all$vars_SE <- rbind(mod_all$vars_SE , mod1$vars_SE[i,],mod2$vars_SE[i,])
	for(i in 1:6) mod_all$fixed_mean <- rbind(mod_all$fixed_mean , mod1$fixed_mean[i,],mod2$fixed_mean[i,])
	for(i in 1:6) mod_all$fixed_SE <- rbind(mod_all$fixed_SE , mod1$fixed_SE[i,],mod2$fixed_SE[i,])
	for(i in 1:6) mod_all$fixed_var <- c(mod_all$fixed_var , mod1$fixed_var[i],mod2$fixed_var[i])
	return(mod_all)
}
modT_all<-add_mod(modT,modT_M)
modHB_all<-add_mod(modHB,modHB_M)
modW_all<-add_mod(modW,modW_M)
modM_all<-add_mod(modM,modM_M)

# # par(mfrow=c(2,2))
# # plot_results(modT,c(-0.02,0.2),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"),"Tarsus")
# # legend("topright",c("DS - O","A - O","DS - O_E","A - O_E","DS - E","A - E"),col =rep(c(1,2),3), pch = rep(15:17,each=2))
# # plot_results(modHB,c(-0.02,0.2),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"), "Head-Bill")
# # plot_results(modW,c(-0.02,7),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"), "Wing")
# # plot_results(modM,c(-0.02,0.7),4,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid"), "Mass")
# setEPS()
# pdf(paste0(wd,"/skew/Plots/dam_sire_ped.pdf"), height=10, width=20)
# {
# par(mfrow=c(2,2))
# plot_results(modT_all,c(-0.02,0.2),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"),"Tarsus")
# legend("topright",c("pedP","pedM","","DS","A","","V_O","V_O + Egg","Egg"),col =c(1,2,0,1,1,0,1,1,1), pch = c(15,15,15,15,22,15,22,21,24),)
# plot_results(modHB_all,c(-0.02,0.2),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"),"Head-Bill")
# plot_results(modW_all,c(-0.02,7),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"),"Wing")
# plot_results(modM_all,c(-0.02,0.7),4,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid"),"Mass")
# }
# dev.off()


# plot_results(modT,c(-0.02,0.2),5,axis_names=c("Fixed","A","Nest_Rear", "Nest_Orig","Resid","ME"),"Tarsus")



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
plot_results1(modM, ylim=c(0,0.8), main="Weight",rows=c(2,4,6), cols=1:5)
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
plot_results1(modM, ylim=c(0,0.8), main="Weight",rows=c(5,6), cols=c(1:3,5))
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
plot_results1(modM_all, ylim=c(0,0.8), main="Weight",rows=c(11,12), cols=c(1:3,5))
plot_results1(modW_all, ylim=c(0,7.5), main="Wing",rows=c(11,12), cols=c(1:3,5:6))
}
dev.off()
