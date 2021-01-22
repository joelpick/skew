
# source("~/Work/Skew/R/selection_gradients_skew_jarrod.R")

## for running remotely
if(!exists("remote_run")){
	rm(list=ls())

	n_sims <- 2000
	ncores <- 8
	rerun <- FALSE
	plot <- TRUE
	normal <- TRUE
	xfoster <- TRUE
	x_recip <- TRUE # reciprocal or round robin xfostering
	trait<-"tarsus_mm"#"weight_g"

}

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(parallel)
library(asreml)
library(MCMCglmm)
library(MASS)
library(mvtnorm)
library(sn)
library(rstan)
library(cubature)


if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/00_functions.R"))

#posterior_mean<-TRUE
cond_term<-c("year", "sex") # terms to condition on 
cont_term<-c("timeC")   	# terms to control for

############################
## load in data
############################

load(paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))
# survival models (including 2010 when eggs weren't weight) and THBW is the data

load(paste0(wd,"Data/Intermediate/stan_dam_data.Rdata"))
# stan_data_weight$X is the design matrix and THBW_egg_noRep the data
# excluding years in which eggs weren't measured and repeat measures

load(paste0(wd,"Data/Intermediate/stan_summary_data.Rdata"))

if(trait=="weight_g"){
 model_w<-mod_weight_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModWeightPedN20191220_0905.Rdata")))
}
if(trait=="tarsus_mm"){
 model_w<-mod_tarsus_bv
 model_z<-get(load(paste0(wd,"Data/Intermediate/stanModTarsus_pedN20191217_1716.Rdata")))
}
rm("mod_weight_bv", "mod_tarsus_bv")



############################
## simulation parameters
############################

mu_pred<-stan_data_weight$X %*% pars(model_z,"beta")[,1]

if(rerun){

e_st <- if(normal){
	list(
		n_st= c(0, sqrt(dp2cm(pars_ST(model_z,"nest")[,1],family="ST")[2]), 0, 1e+16),
		e_st= c(0, sqrt(dp2cm(pars_ST(model_z,if(trait=="weight_g"){"E"}else{"ind"} )[,1],family="ST")[2]), 0, 1e+16), 
		fixed_st = c(mean(mu_pred), sd(mu_pred), 0, 1e+16)
		)
}else{
	list(
		n_st= pars_ST(model_z,"nest")[,1],
		e_st=pars_ST(model_z,if(trait=="weight_g"){"E"}else{"ind"} )[,1], 
		fixed_st=doppelgangR::st.mle(y = mu_pred)$dp
	)
}

g_st<-c(0, pars(model_z,"sigma_A")[1], 0, 1e+16)

beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model_w$Sol))
# linear coefficients for trait values from the survival model

gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model_w$Sol))
# quadratic coefficients for trait values from the survival model

eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model_w$Sol))
eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))
# positions of non-trait terms in the survival model that pertain to Day15 to fledging. 

eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model_w$Sol))
eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))
# positions of non-trait terms in the survival model that pertain to fledging to recruitment. 

fixedEffects <- formula(~ sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC)
# single trait version of the fixed effect formula used in the survival model

X_eta1<-model.matrix(fixedEffects, THBW)
X_eta2<-model.matrix(fixedEffects, THBW)
# design matrix for non-trait effects in survival model

eta1 <- X_eta1%*%colMeans(model_w$Sol[,eta1_pos])
eta2 <- X_eta2%*%colMeans(model_w$Sol[,eta2_pos])
# non-trait fixed effects predictions from survival model

beta<-colMeans(model_w$Sol[,beta_pos])
gamma<-colMeans(model_w$Sol[,gamma_pos])
# trait coefficients from survival model

V_nest<-matrix(colMeans(model_w$VCV[,grep("nest", colnames(model_w$VCV))]), 2,2)

mu_etaz<-c(mean(eta1), mean(eta2), mean(THBW[,paste0(trait, "C")])) 
V_etaz<-cov(cbind(eta1, eta2, THBW[,paste0(trait, "C")]))
# mean and covariance matrices of linear predictor and traits
mu<-attr(THBW[,paste0(trait, "C")], "scaled:center")



############################
## simulation
############################


nf<-1000
no<-10 ## needs to be even
ng<-3

ni<-nf*no	

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

		W<-sapply(g+n+e-mu, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
		# fitness evaluated for each (mean-centred) trait value

		## first half of previous gen pedigree are females/second half males
		## sample who survives - n1/2 of each sex, with their survival probability based in their
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

	ped.ainv1 <- inverseA(ped[,1:3])$Ainv
	assign("ped.ainv", ped.ainv1, envir = .GlobalEnv) 

	m1<-asreml(z~1, random=~nurse+giv(id), data=subset(ped_df, generation!=1), ginverse=list(id=sm2asreml(ped.ainv)),rcov = ~idv(units),trace=FALSE)
	summary(m1)

	## parent-offspring regression
	ped_df$pz <- (ped[ped[,"dam"],"z"] + ped[ped[,"sire"],"z"])/2
	po_df <- na.omit(aggregate(cbind(pz,z)~dam,ped_df,mean))
	
	out <- c(summary(m1)$var[c(1,2,4),2],h2_animal=pin(m1, h2_animal~id/(id+nurse+R))[1,1], h2_po =coef(lm(z~pz,po_df))[2])
	names(out)[1:3] <- c("Vnest","Va","Ve")
	cat(j," ")
	# sims[,j]<-out
	rm("ped.ainv")

	return(out)
}
, mc.cores = ncores)
sims <- do.call(cbind,sims)

save(sims, file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_", if(normal){"N"}else{"ST"}, if(xfoster){"_X"}else{""},if(x_recip){"R"}else{""},".Rdata"))

}


#####
if(plot){

obs_var <- c(
	dp2cm(pars_ST(model_z,"nest")[,1],family="ST")[2],
	pars(model_z,"sigma_A")[1]^2,
	var(mu_pred) + dp2cm(pars_ST(model_z,if(trait=="weight_g"){"E"}else{"ind"} )[,1],family="ST")[2]
	)

sim_h2 <- obs_var[2]/sum(obs_var)

mean_CI2 <- function(x) c(mean(x),mean(x)+se(x)*qnorm(0.975),mean(x)-se(x)*qnorm(0.975))

pars(model_z,"sigma_nest")[1]^2
dp2cm(pars_ST(model_z,"nest")[,1],family="ST")[2]

sim_data <- matrix(rep(c(obs_var,sim_h2,sim_h2),3),ncol=3)

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_ST.Rdata"))
simsST <- t(apply(sims,1,mean_CI2))

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_N.Rdata"))
#hist(sims[1,])
simsN <- t(apply(sims,1,mean_CI2))

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_ST_XR.Rdata"))
simsST_X <- t(apply(sims,1,mean_CI2))

load(file= paste0(wd,"Data/Intermediate/sim_ped_",trait, "_N_XR.Rdata"))
simsN_X <- t(apply(sims,1,mean_CI2))


rownames(simsN_X) <- rownames(simsST) <- rownames(simsST_X) <- rownames(simsN) <- rownames(sim_data) <- c("Vnest","Va","Ve","h2_animal","h2_po")


setEPS()
pdf(paste0(wd,"R/plots/figure_SM_pedSims_",trait, ".pdf"), , height=8, width=12)
{
par(mfrow=c(1,1),mar=c(4,8,1,1))

effectPlot(sim_data,xlim=c(0,0.625),offset=0.2,pch="|")
effectPlot(simsST,col=2,add=TRUE,offset=0.1,pch="|")
effectPlot(simsN,col=3,add=TRUE,offset=0,pch="|")
effectPlot(simsST_X,col=4,add=TRUE,offset=-0.1,pch="|")
effectPlot(simsN_X,col=5,add=TRUE,offset=-0.2,pch="|")

legend("bottomright",c("observed","simulated skew T","simulated normal","simulated skew T with xfoster","simulated normal with xfoster"), pch=19, col=1:5)
abline(h=2.5, lty=3)

}

#abline(v=c(0.132,0.16,0.182))
dev.off()




setEPS()
pdf(paste0(wd,"R/plots/figure_pedSims_",trait, ".pdf"), height=5, width=5)
{
par(mfrow=c(1,1),mar=c(3,5,1,1))

#points(sim_data[4:5,1]~c(0.8,1.8), pch=19, cex=2)
plot(simsST_X[5:4,1]~c(0.9,1.9), xlim=c(0.5,2.5), ylim=c(0.2,0.5), ylab="Estimated Heritability", pch=19, cex=2, col="red", xaxt="n",xlab="")
axis(1,c(1,2),c("P-O Regression","Animal Model"))
 points(simsST[5:4,1]~c(1.1,2.1), pch=19, cex=2, col="blue")
abline(h=sim_data[4,1], col="grey")
arrows(c(0.9,1.9),simsST_X[5:4,2],c(0.9,1.9),simsST_X[5:4,3], code=3, angle=90, length=0.15)
 arrows(c(1.1,2.1),simsST[5:4,2],c(1.1,2.1),simsST[5:4,3], code=3, angle=90, length=0.15)
legend("topleft",c("With x-foster","Without x-foster"), pch=19, col=c("red","blue"))}

dev.off()



}

simsN[5:4,1]/sim_data[4,1]
simsN_X[5:4,1]/sim_data[4,1]
simsST_X[5:4,1]/sim_data[4,1]
simsST[5:4,1]/sim_data[4,1]

y<-rbind(simsN[5:4,1],simsN_X[5:4,1],simsST[5:4,1],simsST_X[5:4,1])
yL<-rbind(simsN[5:4,2],simsN_X[5:4,2],simsST[5:4,2],simsST_X[5:4,2])
yU<-rbind(simsN[5:4,3],simsN_X[5:4,3],simsST[5:4,3],simsST_X[5:4,3])

x<-matrix(c(1,2),byrow=TRUE,ncol=2,nrow=4)+(matrix(1:4,ncol=2,nrow=4)-2.5)/10
{
plot(y~x, xlim=c(0.5,2.5), ylim=c(0.2,0.5), ylab="Estimated Heritability",  cex=2, pch=c(24,25),bg=c(1,1,0,0), xaxt="n",xlab="")
axis(1,c(1,2),c("P-O Regression","Animal Model"))
abline(h=sim_data[4,1], col="grey")
arrows(x,yU,x,yL, code=3, angle=90, length=0.15)
legend("topleft",c("Normal","Skewed","No x-foster","X-foster"), pch=c(21,21,24,25), pt.bg=c(1,0,"grey","grey"))
}
# setEPS()
# pdf("/Users/joelpick/Dropbox/0_Presentations/images/skew_h2_sim2.pdf", height=5, width=5)
# {
# par(mfrow=c(1,1),mar=c(3,5,1,1))

# #points(sim_data[4:5,1]~c(0.8,1.8), pch=19, cex=2)
# plot(simsST_X[5:4,1]~c(0.9,1.9), xlim=c(0.5,2.5), ylim=c(0.1,0.2), ylab="Estimated Heritability", pch=19, cex=2, col="red", xaxt="n",xlab="")
# axis(1,c(1,2),c("P-O Regression","Animal Model"))
# # points(simsST[5:4,1]~c(1.1,2.1), pch=19, cex=2, col="blue")
# abline(h=sim_data[4,1], col="grey")
# arrows(c(0.9,1.9),simsST_X[5:4,2],c(0.9,1.9),simsST_X[5:4,3], code=3, angle=90, length=0.15)
# # arrows(c(1.1,2.1),simsST[5:4,2],c(1.1,2.1),simsST[5:4,3], code=3, angle=90, length=0.15)
# legend("topleft",c("With x-foster","Without x-foster"), pch=19, col=c("red","blue"))
# }

# dev.off()

# simsST["Va",1] - sim_data["Va",1]
# simsST["Ve",1] - sim_data["Ve",1]
# simsST["Vnest",1] - sim_data["Vnest",1]
# (simsST["h2_animal",1] - sim_data["h2_animal",1])/sim_data["h2_animal",1]
# (simsST["h2_po",1] - sim_data["h2_po",1])/sim_data["h2_animal",1]

# (simsST["Va",1] - sim_data["Va",1])/sim_data["Va",1]


# rowMeans(sims)

# true_h2<-rep(obs_var["A",1]/sum(obs_var[,1]), n_sims)
# true_h2a<-rep(0.192137, n_sims)


# par(mfrow=c(2,2))
# hist(sims[2,],breaks=20, main="asreml Va");abline(v=obs_var["A",1], col="red")
# hist(sims[4,],breaks=20, main="asreml h2");abline(v=obs_var["A",1]/sum(obs_var[,1]), col="red");abline(v=0.192137, col="blue")
# hist(sims[5,],breaks=20, main="po-regression");abline(v=obs_var["A",1]/sum(obs_var[,1]), col="red");abline(v=0.192137, col="blue")

