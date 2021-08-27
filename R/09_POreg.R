rm(list=ls())

options(stringsAsFactors=FALSE)

# library(MCMCglmm)
# library(MASS)
# library(mvtnorm)
library(sn)
library(rstan)
library(cmdstanr)
library(coda)
library(cubature)
library(parallel)
library(viridis)


if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

se2ci <- function(x) c(x[1],x[1]+x[2]*qnorm(0.025),x[1]+x[2]*qnorm(0.975))

load(paste0(wd,"Data/Intermediate/analysis_options.Rdata"))

reduced <- analysis_options$reduced
cond_term<-analysis_options$cond_term  # terms to condition on 
cont_term<-analysis_options$cont_term  # terms to control for 

ncores <- 8
re_run_POreg<-FALSE
n_it <- 1000
nplot_points<- 100
save_plot<-FALSE

load(paste0(data_wd,"chick_data.Rdata"))
# THBW is the full data
# THBW_noRep excluding repeat measures

load(paste0(data_wd,"stan_data",if(reduced)"_reduced",".Rdata"))
# data used in stan models, stan_data_ped_noRep$X is the design matrix (without repeat measures)

load(paste0(wd,"Data/Intermediate/starting_values",if(reduced)"_reduced",".Rdata"))


h2_all <- list()
if(save_plot){
	setEPS()
	pdf(paste0(wd,"R/plots/figure5.pdf"), , height=6, width=12)
}
mat <- matrix(c(1:16),ncol=4)
mat2 <- cbind(c(0,0,17,0,0),rbind(mat[1:3,],c(0,18,18,0),mat[4,]))
layout(mat2,height=c(1,2,8,0.5,4),width=c(1,6,6,6,6))
traits <-  c("tarsus_mm","headbill_mm","weight_g","wing_mm")
traits_lab <- paste0(sub("_"," (",traits),")")
traits_lab <- sub("weight","mass",traits_lab)
substr(traits_lab,1,1) <- LETTERS[match(substr(traits_lab,1,1),letters)]
cex.heading = 1.75
cex.axis=1

#,"weight_g"
# if(re_run_POreg){

for (trait in traits){
	#trait="wing_mm"


	model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	## this gives the posterior distributions for stan model for the trait (model_z), as an mcmc.list object
	model_zpost<-do.call(rbind,model_z)[,-(1:7)]
	##combine into one matrix

	z<-THBW_noRep[[trait]]

	#design matrix for trait model
	X<-stan_data_ped_noRep$X

	model_zbeta<-pars(model_z, "beta")[,1]
	# posterior mean for location effects from trait model

	pred_pos <- which(!grepl((paste(c(cond_term,cont_term),collapse="|")),colnames(X)))
	#gets the position of the effects that are neither controlling or conditioning (marginal variables)
	
	predC_pos <- grep((paste(c(cond_term,cont_term),collapse="|")),colnames(X))
	#gets the position of the effects that are either controlling or conditioning 
	
	mu_pred <- X[,pred_pos,drop=FALSE]%*%model_zbeta[pred_pos]
	# gets the posterior mean predictions for the trait models for marginal variables. 
	
	mu_predC <- X[,predC_pos,drop=FALSE]%*%model_zbeta[predC_pos]
	# gets the posterior mean predictions for the trait models for marginal variables. 

	z_corrected <- z - mu_predC


	if(re_run_POreg){

		###################################
		##
		## estimate non-linear parent-offspring regression at the posterior mean of the trait model 
		##
		###################################
		
		Zplot_points <- sort(unique(c(seq(min(z),max(z),length.out=nplot_points), z_corrected[THBW_noRep$bird_id%in%c(THBW_noRep$dam_P, THBW_noRep$sire_P)])))
		## phenotypes used for estimating the parent offspring regression - equally spaced points across phenotypic distribution as well as parent phenotypes

		e_st<-list(n_st=pars_ST(model_z, "nest")[,1],
			       e_st=pars_ST(model_z, if(trait=="weight_g"){"E"}else{"ind"})[,1], 
			       fixed_st=st.mple(y = mu_pred)$dp)

		if(trait!="weight_g") e_st$me_st<-c(0, pars(model_z, "sigma_E")[,1], 0, 1e+16)

		# list of environmental distribution parameters: xi, omega, alpha, nu
		# fixed_st is a skew-t approximation for the marginal variable predictions 

		g_st<-c(0, pars(model_z, "sigma_A")[,1], 0, 1e+16)
		# list of genetic distribution parameters

		mu<-dp2cm(c(list(g_st), e_st), family="ST")

		limit_prob=1e-4
		e_llimits<-unlist(lapply(e_st, function(x){qst(limit_prob, dp=x)}))
		e_ulimits<-unlist(lapply(e_st, function(x){qst(1-limit_prob, dp=x)}))
  			# lower and upper limits when integrating over environmental effects



		cat("\n", trait," - Start: ", as.character(Sys.time()), " - Total points = ",length(Zplot_points),"\n")

# range(THBW[,"wing_mm"])
# mean(THBW[,"wing_mm"])
# range(z_corrected)

		# dz_p<-NULL    # phenotypic density function
		# Eg_p<-NULL    # parent-offspring regression 
		# dEg_p<-NULL   # derivative of parent-offspring regression 

		POreg_list <- mclapply(Zplot_points,function(j){
			# dz_p<- dz(j, g_st=g_st, e_st=e_st)
			# Eg_p <- POreg(j,  g_st=g_st, e_st=e_st)
			# dEg_p <- dPOreg(j, g_st=g_st, e_st=e_st, Eg_p=Eg_p)


			E_p<-hcubature(conv, e_llimits, e_ulimits, z_p=j, g_st=g_st, e_st=e_st)$integral
			Eg_p<-hcubature(wconv, e_llimits, e_ulimits, z_p=j, g_st=g_st, e_st=e_st)$integral
		    Eg_p2<-hcubature(wconv2, e_llimits, e_ulimits, z_p=j, g_st=g_st, e_st=e_st)$integral
		    Eg_p<-Eg_p/E_p
		    Eg_p2<-Eg_p2/E_p
		    dEg_p <- 0.5*(1+(Eg_p^2-Eg_p2)/(g_st[2]^2))
		    
			cat(which(Zplot_points==j)," ")
			return(c(z=j,dz_p=E_p,Eg_p=Eg_p*0.5,dEg_p=dEg_p))
		}, mc.cores = ncores)
		
		cat("\n")

		# phenotypic density function assuming normality
		dzn_p <- dnorm(Zplot_points, mu[1], sqrt(mu[2])) 

		POreg_out <- cbind(do.call(rbind,POreg_list),dzn_p=dzn_p)

		save(z_corrected,POreg_out, file=paste0(wd,"Data/Intermediate/nonLinearPO_", trait,".Rdata"))
	}else{
		load(paste0(data_wd,"nonLinearPO_",trait,".Rdata"))
	# if(trait!="wing_mm") load(paste0(data_wd,"nonLinearPO_",trait,".Rdata"))

	}
		# POreg_out <- cbind(POreg_out,dzn_p=dzn_p)
	trait_col <- inferno(5)[which(traits%in%trait)]
	trait_lab <- traits_lab[which(traits%in%trait)]


	load(paste0(data_wd,"h2_",trait,".Rdata"))

	z_om<-tapply(z_corrected, THBW_noRep$dam_P, mean, na.rm=TRUE)
	z_nm<-tapply(z, THBW_noRep$dam_P, function(x){sum(!is.na(x))})
	z_m<-(z_corrected)[match(names(z_om), THBW_noRep$bird_id)]

	z_of<-tapply(z_corrected, THBW_noRep$sire_P, mean, na.rm=TRUE)
	z_nf<-tapply(z, THBW_noRep$sire_P, function(x){sum(!is.na(x))})
	z_f<-(z_corrected)[match(names(z_of), THBW_noRep$bird_id)]

	z_o<-c(z_om, z_of)
	z_p<-c(z_m, z_f)
	z_n<-c(z_nm, z_nf)

	complete<-which(!is.na(z_p) & !is.na(z_o))

	z_o<-z_o[complete]
	z_p<-z_p[complete]
	z_n<-z_n[complete]


	## need to compare linear po, non-linear po and normal h2 with ME 
	m1<-lm(z_o~z_p, weights=z_n)
	
	m2<-lm(z_o~offset(I(POreg_out[,"Eg_p"])[match(z_p,POreg_out[,"z"])]), weights=z_n)
	m2b<-lm(z_o~offset(I(POreg_out[,"Eg_p"])[match(z_p,POreg_out[,"z"])])+z_p, weights=z_n)

	m3<-lm(z_o~offset(mean(h2bN_ME/2)*z_p), weights=z_n)
	m3b<-lm(z_o~offset(mean(h2bN_ME/2)*z_p)+z_p, weights=z_n)

	xlim <- range(z_corrected)
	# xlim <- range(POreg_out[,"z"])
	
	par(mar=c(0,0,0,0))
	blankPlot(); text(0,0,trait_lab,cex=cex.heading)
	
	par(mar=c(0,4,0,1))
	boxplot(list(z_p,z_corrected), horizontal=TRUE, ylim=xlim, yaxt="n", xaxt="n", frame=FALSE,col=alpha(c("blue","red"),0.5), border=c(1,1), names=if(trait=="tarsus_mm") c("survived","all") else c("","")) #, col=alpha(trait_col,c(0.3,0.7)

	par(mar=c(2,3,1,1))
	plot(z_o~z_p, cex=0.3*sqrt(z_n), pch=19, col=alpha(1,0.3), xlim=xlim, ylim=range(z_o),  xlab="", ylab= "")
	# if(trait!="wing_mm") 
		lines(I(POreg_out[,"Eg_p"]+ coef(m2))~POreg_out[,"z"], lty=1, col="red", lwd=2)
	z_pred1<-predict(m1, newdata=list(z_p=range(z_p)), interval="confidence")
	lines(z_pred1[,1]~range(z_p), lty=1, col="blue", lwd=2)
		
	# z_pred4<-coef(m4) + mean(h2b/2) * range(z_p)
	# lines(z_pred4~range(z_p), lty=1, col="orange", lwd=2)


	# h2_all <- rbind(
	# 	h2obs =se2ci(coef(summary(m1))[2,1:2]*2),
	# 	h2a=mean_CI(rowMeans(h2a)),
	# 	h2b=mean_CI(h2b))
	# 	#, diff=mean_CI(rowMeans(h2a)-h2b))
p_star <- function(x) ifelse(x>0.05,"NS", ifelse(x>0.005,"*","**"))

	h2s<- rbind(
		# h2obs =se2ci(coef(summary(m1))[2,1:2]*2),
		## need to correct this one for ME
		# h2bN=mean_CI(h2bN),
		h2a=mean_CI(rowMeans(h2a)),
		h2b=mean_CI(h2b))

	max_x <- min(max(h2s[,1])*2.2,1)
	x <- max_x * c(0.85)
	par(mar=c(3,4,2,1))
	# effectPlot(h2s, col= c("blue","red","blue","red"), names=c(expression(italic(h[PO]^2)),expression(italic(h[bN]^2)),expression(italic(h[a]^2)),expression(italic(h[b]^2))), xlim=c(0,max_x), cex.axis.scale=1)
	# abline(h=2.5, lty=2)
	
 	# arrows(x,3,x,4,length=0)
 	# text(x + max_x*0.075, c(3.5), p_star(anova(m1, m3)$`Pr(>F)`[2]),cex=1)

 	# arrows(x,1,x,2,length=0)
 	# text(x + max_x*0.075, c(1.5), p_star(pMCMC(rowMeans(h2a)-h2b)),cex=1)

	effectPlot(h2s, col= c("blue","red"), names=c(expression(italic(h[a]^2)),expression(italic(h[b]^2))), xlim=c(0,max_x), cex.axis.scale=1)
 	arrows(x,1,x,2,length=0)
 	text(x + max_x*0.075, c(1.5), p_star(pMCMC(rowMeans(h2a)-h2b)),cex=1)


	h2_all[[trait]] <- list(rbind(h2s,
		diff = mean_CI(rowMeans(h2a)-h2b),
		prop = mean_CI(rowMeans(h2a)/h2b)),
		prop = rowMeans(h2a)/h2b,
		po_vs_nlpo = anova(m2b, m2)$`Pr(>F)`[2], 
		po_vs_h2a=anova(m3b, m3)$`Pr(>F)`[2])

}
par(mar=c(0,0,0,0))
blankPlot()
text(0,0,"Offspring",srt=90,cex=1.5)
blankPlot()
text(0,0,"Parent",cex=1.5)
 
if(save_plot) dev.off()

save(h2_all, file=paste0(data_wd,"h2_results.Rdata"), version=2)


# load(paste0(wd,"Data/Intermediate/starting_values.Rdata"))

# var(X%*%model_zbeta)+sum(pars(model_z,"sigma")[,1]^2)
# var(z)
# sum(modT$mod_A$V) + var(X%*%modT$mod_A$F)

# rbind(c(fixed=var(X%*%modT$mod_A$F),modT$mod_A$V),c(var(X%*%model_zbeta),pars(model_z,"sigma")[,1]^2)

# 	)

		# m_summary<-list(m1.2=anova(m1, m2)$`Pr(>F)`[2], m1.3=anova(m1, m3)$`Pr(>F)`[2], m1.4=anova(m1, m4)$`Pr(>F)`[2], po_reg2=coef(summary(m1))[2,1:2]*2)

		# if(save){
		# 	save(m_summary, file=paste0(wd,"Data/Intermediate/", gsub("selection_gradient", "selection_gradient_msummary", files)), version=2)
		# 	save(m_summary, file=paste0(wd,"Data/Intermediate/", gsub(selection_gradient", "selection_gradient_msummary", files_date)), version=2)
	


# zpoints<-"cparents+0.1" 
# po-regression and fitness function to be evaluated at trait values: 
# even: spread evenly over the range    
# parents: equal to those that become parents
# parents+: equal to those that become parents + evenly spaced phenotypes up to the min/max
# cparents/cparents+: equal to those that become parents minus the posterior mean prediction from the cond/contr variables
# numerical value after even and parents+ determines the spacing

		# if(grepl("even", zpoints)){
		# 	Zplot_points<-seq(min(z), max(z), as.numeric(gsub("even", "", zpoints)))
		# }
		# if(grepl("parents", zpoints)){
		#     if(grepl("cparents", zpoints)){
		#     	z_corrected <- z - mu_predC	
		# 	}else{
		# 		z_p_corrected <- z
		# 	}  
		# 	Zplot_points <- sort(unique(z_p_corrected[THBW_noRep$bird_id%in%c(THBW_noRep$dam_P, THBW_noRep$sire_P)]))
		# }
		# if(grepl("parents\\+", zpoints)){	
		# 	Zplot_points<-sort(unique(c(Zplot_points, seq(min(z), max(z), as.numeric(gsub("c?parents\\+", "", zpoints))))))
		# }  
		# nplot_points<-length(Zplot_points)
