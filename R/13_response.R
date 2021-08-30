rm(list=ls())

options(stringsAsFactors=FALSE)

library(viridis)
library(scales)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")
source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"chick_data.Rdata"))

save_plot <- FALSE
reduced <- TRUE
traits <-  c("tarsus_mm","headbill_mm","weight_g","wing_mm")

response <- list()

for(trait in traits){
	load(paste0(data_wd,"selection_gradients_ME_",trait,".Rdata"))
	# load(paste0(data_wd,"h2_",trait,".Rdata"))
	model_files <- list.files(data_wd)[grep(paste0("stanModNormal2_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	V_A <- do.call(rbind,model_z)[,"sigma_A"]^2
	# pars(model_z,"sigma_A")
	response[[trait]] <- rowMeans(beta1) * V_A

}

save(response,file=paste0(data_wd,"response.Rdata"),version=2)



{if(save_plot){
	setEPS()
	pdf(paste0(wd,"R/plots/figure_SM_heywood.pdf"), , height=8.5, width=15)
}


layout(matrix(c(1:16),ncol=4),height=c(1,5,5,5))

# par(mfcol=c(3,4), mar=c(5,5,1,1), cex.lab=1.5)

traits_lab <- paste0(sub("_"," (",traits),")")
traits_lab <- gsub("weight","mass",traits_lab)
substr(traits_lab,1,1) <- LETTERS[match(substr(traits_lab,1,1),letters)]
cols <- inferno(5)

for(trait in traits){
	load(paste0(data_wd,"selection_gradients_ME_",trait,".Rdata"))
	load(paste0(data_wd,"nonLinearPO_",trait,".Rdata"))
	load(paste0(data_wd,"h2_",trait,".Rdata"))

	model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))

	mu_PO <- pars(model_z,"beta.1")[1,1]
	mu <- mean(THBW_noRep[,trait])


	POreg_out2 <- POreg_out[match(Wplot_points,POreg_out[,"z"]),]

	Wplot_pred <- colMeans(Wplot)

	beta1_pred <-  mean(THBW_noRep$recruit) + mean(THBW_noRep$recruit)*mean(beta1) * (Wplot_points-mu )

    trait_lab <- traits_lab[which(traits%in%trait)]
    trait_lab2 <- gsub(" .+","",trait_lab)
	par(mar=c(0,0,0,0), cex.lab=1.5)
	blankPlot(); text(0,0,trait_lab2,cex=2)

	par(mar=c(5,5,1,1), cex.lab=1.5)

	plot(Wplot_pred~Wplot_points, type="l", ylim=range(c(Wplot_pred,beta1_pred)), xlab=trait_lab, ylab="W", lwd=2)
	lines(beta1_pred~Wplot_points, col="red", lwd=2)

	W_resid <- Wplot_pred - beta1_pred

	h2 <- if(trait!="weight_g"){mean(h2bN_ME)}else{mean(h2bN)}
	nlpo_pred <- POreg_out2[,"Eg_p"] + mu
	density_pred <- POreg_out2[,"dz_p"]
	density_scaled <- density_pred/(max(density_pred))

	po_pred <- h2/2*(Wplot_points-mu) + mu

	plot(nlpo_pred~Wplot_points, type="l", ylim=range(c(nlpo_pred,po_pred)), xlab="Parent", ylab="Offspring", lwd=2)
	lines(po_pred~Wplot_points, col="red", lwd=2)
		
	PO_resid <- nlpo_pred-po_pred

	 plot(W_resid~PO_resid, type="l", col="grey", xlab="P-O residuals", ylab="W residuals", lwd=2)
	# points(W_resid~PO_resid, cex=density_scaled, col="red", pch=19)
	 points(W_resid~PO_resid, cex=1, col=alpha("red",density_scaled), pch=19)
	line_coef <- coef(lm(W_resid~PO_resid, weights=density_pred))
	# abline(line_coef, col="blue")
	# print(c(cov(W_resid,PO_resid),cov.wt(cbind(W_resid,PO_resid),wt=density_pred,cor=TRUE)$cor[1,2],line_coef[2]*var(PO_resid)))


	rm(Wplot_points)
}

if(save_plot) dev.off()

}
