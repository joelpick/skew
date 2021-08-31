
rm(list=ls())

n_sims <- 2000
ncores <- 8
rerun <- FALSE
plot <- TRUE
save_plot <- TRUE
normal <- TRUE
xfoster <- TRUE
x_recip <- TRUE # reciprocal or round robin xfostering
trait<-"weight_g"

options(stringsAsFactors=FALSE)

library(parallel)
library(asreml) # asreml4
library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)
library(rstan)
library(cubature)


if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/github/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"analysis_options.Rdata"))

reduced <- analysis_options$reduced
cond_term<-analysis_options$cond_term  # terms to condition on 
cont_term<-analysis_options$cont_term  # terms to control for 
fixed_w <- analysis_options$fixed_w

############################
## load in data
############################

load(paste0(data_wd,"chick_data.Rdata"))

load(paste0(data_wd,"stan_data",if(reduced)"_reduced",".Rdata"))
# stan_data_weight$X is the design matrix and THBW_egg_noRep the data
# excluding years in which eggs weren't measured and repeat measures

model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
load(paste0(data_wd,model_files[length(model_files)]))
model_zpost<-do.call(rbind,model_z)[,-(1:7)]


model_files <- list.files(data_wd)[grep(paste0("day15_survival_ME_",trait),list.files(data_wd))]
load(paste0(data_wd,model_files[length(model_files)]))
model_wpost<-do.call(rbind,model_w)[,-(1:7)]

# rm("model_z", "model_w")



############################
## simulation parameters
############################

### from trait model

	X<-stan_data_ped_noRep$X

	model_zbeta<-pars(model_z, "beta")[,1]
	pred_pos <- which(!grepl((paste(c(cond_term,cont_term),collapse="|")),colnames(X)))
	mu_pred <- X[,pred_pos,drop=FALSE]%*%model_zbeta[pred_pos]
	
	fixed_st <- st.mple(y = mu_pred)$dp
	nest_st <- pars_ST(model_z, "nest")[,1]
	resid_st <- pars_ST(model_z, if(trait=="weight_g"){"E"}else{"ind"})[,1]
	g_st<-c(0, pars(model_z,"sigma_A")[1,1], 0, 1e+16)


nf<-1000
no<-10 ## needs to be even
ng<-3

ni<-nf*no	




if(rerun){


	model_w_out<- extract_w2(model_wpost, trait, fixedEffects=formula(paste("~", fixed_w)), data=THBW_noRep)
	attach(model_w_out)

	eta1 <- rowMeans(eta1_all)
	eta2 <- rowMeans(eta2_all)
	# non-trait fixed effects predictions from survival model

	beta<-colMeans(beta_all)
	gamma<-colMeans(gamma_all)
	# trait coefficients from survival model

	V_nest<-matrix(colMeans(V_nest_all), 2,2)

	mu_etaz<-c(mean(eta1), mean(eta2), mean(THBW_noRep[,paste0(trait, "C")])) 
	V_etaz<-cov(cbind(eta1, eta2, THBW_noRep[,paste0(trait, "C")]))
	# mean and covariance matrices of linear predictor and traits

	mu<-mean(THBW[,trait])


	selection<-TRUE
	for(normal in c(TRUE, FALSE)){
		for(xfoster in c(TRUE, FALSE))
		{
#selection=FALSE
# normal=FALSE
# xfoster=FALSE
	cat("\n","Starting: normal =",normal,"and xfoster =",xfoster," ... \n Time: ", as.character(Sys.time()),"\n")

	e_st <- if(normal){
		list(
			n_st= c(0, sqrt(dp2cm(nest_st,family="ST")[2]), 0, 1e+16),
			e_st= c(0, sqrt(dp2cm(resid_st,family="ST")[2]), 0, 1e+16), 
			fixed_st = c(dp2cm(fixed_st,family="ST")[1],sqrt(dp2cm(fixed_st,family="ST")[2]), 0, 1e+16)
			)
	}else{
		list(
			n_st=nest_st,
			e_st=resid_st, 
			fixed_st=fixed_st
		)
	}







	############################
	## simulation
	############################


	## nurse id - depends on cross fostering
	nurse_id <- 
		if(xfoster){ 
			if(x_recip){
				rep(as.vector(rbind(1:nf,(nf:1))),each=no/2) ## reciprocal xfoster
			}else{
				rep(1:nf, each=no)[c((no/2+1):(ni),1:(no/2))] ## round robin xfoster
			}
		}else{
			rep(1:nf, each=no)
		}


	sims <- mclapply(1:n_sims,function(j){
	#j=1
	# sims<-matrix(NA,5,n_sims)
	# for(j in 1:n_sims){

		ped<-matrix(NA, nf*no*ng, 6)
		## animal, dam, sire, phenotype(z), generation, nurse

		## generate breeding values, nest and residual effects
		g<-rz(ni, list(g_st))
		n<-rz(nf, list(e_st$n_st))[nurse_id]
		e<-rz(ni, e_st[-which(names(e_st)=="n_st")])

		ped[1:ni,1]<-1:ni
		ped[1:ni,4]<-g+n+e
		ped[1:ni,5]<-1
		ped[1:ni,6]<- as.numeric(paste0(1, nurse_id))

		for(i in 2:ng){
	# i=1
			W<-if(selection){
				sapply(g+n+e-mu, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
			}else{
				rep(1,ni)
			}

			# fitness evaluated for each (mean-centred) trait value

			## first half of previous gen pedigree are females/second half males
			## sample who survives - ni/2 of each sex, with their survival probability based in their size
			dam<-sample(1:(ni/2), nf, prob=W[1:(ni/2)], replace=FALSE)
			sire<-sample((ni/2)+1:(ni/2), nf, prob=W[(ni/2)+1:(ni/2)], replace=FALSE)

			## put surviving individual into pedigree as parents with no offspring each (random mating)
			ped[(i-1)*ni+1:ni,1]<-(i-1)*ni+1:ni
			ped[(i-1)*ni+1:ni,2]<-rep(dam+ni*(i-2), each=no)
			ped[(i-1)*ni+1:ni,3]<-rep(sire+ni*(i-2), each=no)

			## generate new breeding values, nest and residual effects
			g<-(g[rep(dam, each=no)]+g[rep(sire, each=no)])/2+rz(ni, list(g_st))/sqrt(2)
			n<-rz(nf, list(e_st$n_st))[nurse_id]
			e<-rz(ni, e_st[-which(names(e_st)=="n_st")])

			ped[(i-1)*ni+1:ni,4]<-g+n+e
			ped[(i-1)*ni+1:ni,5]<-i
			ped[(i-1)*ni+1:ni,6]<- as.numeric(paste0(i, nurse_id))
		}

		colnames(ped)<-c("id", "dam", "sire", "z", "generation","nurse")
		ped_df<-as.data.frame(ped)
		ped_df$id<-as.factor(ped_df$id)
		ped_df$dam<-as.factor(ped_df$dam)
		ped_df$sire<-as.factor(ped_df$sire)
		ped_df$z<-as.numeric(ped_df$z)
		ped_df$generation<-as.factor(ped_df$generation)
		ped_df$nurse<-as.factor(ped_df$nurse)

		#table(table(ped_df$nurse))
		#table(table(paste(ped_df$dam,ped_df$nurse)))

		## animal model
		# ped.ainv1 <- asreml.Ainverse(ped[,1:3])
		# assign("ped.ainv", ped.ainv1$ginv, envir = .GlobalEnv) 

		ped.ainv1<-ainverse(ped[,1:3])
		## make factors for asreml random effects
		ped_df$animal <- factor(ped_df$id, levels=ped[,1])

		# ped.ainv1 <- inverseA(ped[,1:3])$Ainv
		 assign("ped.ainv", ped.ainv1, envir = .GlobalEnv) 

		m1<-asreml(z~1, 
		random=~nurse+vm(animal,ped.ainv), 
		residual = ~idv(units),
		data=subset(ped_df, generation!=1), 
		trace=FALSE)

		m1_var <- summary(m1)$var[-4,1]

		## parent-offspring regression
		ped_df$pz <- (ped[ped[,"dam"],"z"] + ped[ped[,"sire"],"z"])/2
		po_df <- na.omit(aggregate(cbind(pz,z)~dam,ped_df,mean))
		
		out <- c(m1_var,h2_animal= m1_var[2]/sum(m1_var), h2_po =coef(lm(z~pz,po_df))[2])
		names(out)[1:3] <- c("Vnest","Va","Ve")
		cat(j," ")
		# sims[,j]<-out
		rm("ped.ainv")

		return(out)
	}
	, mc.cores = ncores)
	sims <- do.call(cbind,sims)

	save(sims, file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_", if(normal){"N"}else{"ST"}, if(xfoster){"_X"}else{""},if(x_recip){"R"}else{""},if(selection){""}else{"_noSel"},".Rdata"))

}
}
}
detach(model_w_out)

#####
if(plot){


obs_var <- c(
	dp2cm(pars_ST(model_z,"nest")[,1],family="ST")[2],
	pars(model_z,"sigma_A")[1,1]^2,
	var(mu_pred) + dp2cm(pars_ST(model_z,if(trait=="weight_g"){"E"}else{"ind"} )[,1],family="ST")[2]
	)

sim_h2 <- obs_var[2]/sum(obs_var)

mean_CI2 <- function(x) c(mean(x),mean(x)-se(x)*qnorm(0.975),mean(x)+se(x)*qnorm(0.975))
mean_se <- function(x) c(mean(x),se(x))

# pars(model_z,"sigma_nest")[1]^2
# dp2cm(pars_ST(model_z,"nest")[,1],family="ST")[2]

sim_data <- matrix(rep(c(obs_var,sim_h2,sim_h2),3),ncol=3)

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_STR.Rdata"))
simsST <- sims#t(apply(sims,1,mean_CI2))
sumST <- t(apply(simsST,1,mean_CI2))

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_NR.Rdata"))
#hist(sims[1,])
simsN <- sims#
sumN <- t(apply(simsN,1,mean_CI2))

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_ST_XR.Rdata"))
simsST_X <- sims#
sumST_X <- t(apply(simsST_X,1,mean_CI2))

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_N_XR.Rdata"))
simsN_X <- sims#
sumN_X <- t(apply(simsN_X,1,mean_CI2))




{
if(save_plot){
setEPS()
pdf(paste0(wd,"Plots/figure6.pdf"), height=5, width=6)
}

{par(mar=c(5,10,1,1))
effectPlot(sumN[5:4,], offset=0.3, xlim=c(0.1,0.3), names=c("P-O \nRegression","Animal Model"), xlab="Estimated Heritability", pch=24,bg=1, cex.axis=1, cex.lab=1.5, arrow.length=0.1)
effectPlot(sumN_X[5:4,], add=TRUE, offset=0.1, pch=25,bg=1, arrow.length=0.1)
effectPlot(sumST[5:4,], add=TRUE, offset=-0.1, pch=24,bg=0, arrow.length=0.1)
effectPlot(sumST_X[5:4,], add=TRUE, offset=-0.3, pch=25,bg=0, arrow.length=0.1)
abline(v=sim_data[4,1], col="grey")
legend("topright",c("Normal","Skewed","No x-foster","X-foster"), pch=c(21,21,24,25), pt.bg=c(1,0,"grey","grey"))
}



if(save_plot) dev.off()

}

}

