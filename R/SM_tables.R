rm(list=ls())

options(stringsAsFactors=FALSE)

library(coda)
library(MCMCglmm)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
data_wd <- paste0(wd,"Data/Intermediate/")


load(paste0(data_wd,"table_data.Rdata"))

knitr::kable(model_w_table$tarsus_mm)

knitr::kable(model_z_table_A$tarsus_mm)
model_z_table_DS,
model_z_table_N,

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"analysis_options.Rdata"))

reduced <- analysis_options$reduced
cond_term <- analysis_options$cond_term  # terms to condition on 
cont_term <- analysis_options$cont_term  # terms to control for 
fixed_w <- analysis_options$fixed_w

load(paste0(data_wd,"chick_data.Rdata"))
load(paste0(data_wd,"stan_data",if(reduced)"_reduced",".Rdata"))

table_pci <- function(x, digits=3, pMCMC=TRUE){
	x<-as.mcmc(x)
	point_mean<-mean(x)
    point_mode<-MCMCglmm::posterior.mode(x)
    interval<-round(as.numeric(coda::HPDinterval(x)), digits)



c(paste0(formatC(point_mean,digits=digits, format="f"), " (", formatC(point_mode,digits=digits, format="f"),") [", formatC(interval[1],digits=digits,format="f"), ", ", formatC(interval[2], digits=digits,format="f"), "]"),formatC(pMCMC(x),digits=digits, format="f"))

	# c(PCI(x),formatC(pMCMC(x),digits=digits, format="f"))
}


model_w_table <- list()
model_z_table_A <- list()
model_z_table_DS <- list()
model_z_table_N <- list()

traits <-  c("tarsus_mm","headbill_mm","wing_mm","weight_g")
traits_lab <- sub("_.+","",traits)
traits_lab <- sub("weight","mass",traits_lab)
substr(traits_lab,1,1) <- LETTERS[match(substr(traits_lab,1,1),letters)]

X_eta<-model.matrix(formula(paste("~", fixed_w)), THBW_noRep)

for (trait in traits){
# trait="tarsus_mm"

	model_files <- list.files(data_wd)[grep(paste0("stanMod_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zA <- model_z
	model_zpostA<-do.call(rbind,model_zA)[,-(1:7)]

	model_files <- list.files(data_wd)[grep(paste0("stanMod_DS_",if(reduced)"reduced_",trait,"_\\d"),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zDS <- model_z
	model_zpostDS<-do.call(rbind,model_zDS)[,-(1:7)]

	model_files <- list.files(data_wd)[grep(paste0("stanModNormal2_pedN_",if(reduced)"reduced_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_zN <- model_z
	model_zpostN<-do.call(rbind,model_zN)[,-(1:7)]


	if(trait!="weight_g") {

	colnames(model_zpostA) <- gsub("_E","_ME",colnames(model_zpostA))
	colnames(model_zpostA) <- gsub("_ind","_E",colnames(model_zpostA))
	
	colnames(model_zpostDS) <- gsub("_E","_ME",colnames(model_zpostDS))
	colnames(model_zpostDS) <- gsub("_ind","_E",colnames(model_zpostDS))

	colnames(model_zpostN) <- gsub("_E","_ME",colnames(model_zpostN))
	colnames(model_zpostN) <- gsub("_ind","_E",colnames(model_zpostN))
	}



	model_files <- list.files(data_wd)[grep(paste0("day15_survival_ME_",trait),list.files(data_wd))]
	load(paste0(data_wd,model_files[length(model_files)]))
	model_wpost<-do.call(rbind,model_w)[,-(1:7)]



  ####	 survival models

	if(trait!="weight_g"){
		beta_pos<-grep("beta_x1_hat\\.", colnames(model_wpost))
		# linear coefficients for trait values from the survival model

		gamma_pos<-grep("beta_x1_hat2\\.", colnames(model_wpost))
		# quadratic coefficients for trait values from the survival model
	}else{
		beta_pos<-grep("beta\\.", colnames(model_wpost))[(2*ncol(X_eta))+c(1,2)]
		# linear coefficients for trait values from the survival model

		gamma_pos<-grep("beta\\.", colnames(model_wpost))[(2*ncol(X_eta))+c(3,4)]
		# quadratic coefficients for trait values from the survival model
	}

	eta2_pos <- grep("beta\\.", colnames(model_wpost))[(1:ncol(X_eta))+ncol(X_eta)]

	eta1_pos <- grep("beta\\.", colnames(model_wpost))[1:ncol(X_eta)]




	eta1_sum<-t(apply(model_wpost[,eta1_pos], 2,table_pci))
	eta2_sum<-t(apply(model_wpost[,eta2_pos], 2, table_pci))

	beta_sum<-t(apply(model_wpost[,beta_pos], 2, table_pci))
	gamma_sum<-t(apply(model_wpost[,gamma_pos], 2, table_pci))

	fledge_sum <- rbind(eta1_sum,beta_sum[1,], gamma_sum[1,])
	recruit_sum <- rbind(eta2_sum,beta_sum[2,], gamma_sum[2,])

	# rownames(fledge_sum) <- rownames(recruit_sum) <- c(colnames(X_eta),trait, paste0(trait,2))

	V_nest_sum <- t(apply(model_wpost[,grep("Sigma_nest", colnames(model_wpost))], 2, table_pci, pMCMC=FALSE))[c(1,2,4)]

	# rownames(V_nest_sum) <- c("Nest Fledge","Covariance","Nest Recruit")


	row_names<-c("Intercept","Sex (M)","Male Presence (1)","Male Presence (2)","Year (2012)","Year (2013)","Year (2014)","Year (2015)","Year (2016)","Year (2017)","Year (2018)","Hatch Day (1)", "Hatch Day (3)","Clutch Size","Nest Hatch Date",NA,NA, "Nest Variance")
	trait_lab <- traits_lab[which(traits%in%trait)]
	row_names[16:17] <- c(trait_lab,paste0(trait_lab,"$^2$"))


	table_out<-rbind(cbind(fledge_sum,"",recruit_sum),c(V_nest_sum[1],"",V_nest_sum[2],V_nest_sum[3],""))
	rownames(table_out)<-row_names
	colnames(table_out) <- c( "Fledging & & & Recruitment \\\\ & Posterior Mean (Mode) [CI]","pMCMC","Covariance", " Posterior Mean (Mode) [CI]","pMCMC")

	model_w_table[[trait]] <- table_out


	##### TRAIT MODELS
	X_z <-stan_data_ped_noRep$X

	# skew_var<-c("xi","omega","alpha","delta","nu")

	skew_var <- c("sigma","delta","nu")

	A_var <- c("sigma_A",
				paste0(skew_var,"_nest"),
				paste0(skew_var,"_E"),
				if(trait!="weight_g") "sigma_ME")
	DS_var <- c(paste0(skew_var,"_dam_sire"),
				paste0(skew_var,"_nest"),
				paste0(skew_var,"_E"),if(trait!="weight_g") "sigma_ME")
	N_var <- paste0("sigma_",c("A","nest","E",if(trait!="weight_g") "ME"))
	
	beta_names <- c( "Intercept", "Year (2012)","Year (2013)","Year (2014)","Year (2015)","Year (2016)","Year (2017)","Year (2018)","Time of Day","Sex (M)", "Egg Mass")

	col_names <- c("Posterior Mean (mode) [CI]","pMCMC")


	A_var2 <- gsub("_","_{",paste0("$\\",A_var,"}$"))
	A_var2[1] <- paste("\\hline ",A_var2[1])

	DS_var2 <- gsub("dam_sire","DS",DS_var)
	DS_var2 <- gsub("_E","_R",DS_var2)
	DS_var2 <- gsub("_","_{",paste0("$\\",DS_var2,"}$"))
	DS_var2[1] <- paste("\\hline ",DS_var2[1])

	N_var2 <- gsub("_","_{",paste0("$\\",N_var,"}$"))
	N_var2[1] <- paste("\\hline ",N_var2[1])


	table_A <- t(apply(model_zpostA[,c(paste0("beta.",1:ncol(X_z)),A_var)], 2,table_pci))
	rownames(table_A) <-  c(beta_names,A_var2)	
	table_A[,2]<-ifelse(grepl("sigma|nu",rownames(table_A)),"-",table_A[,2])

	table_DS <- t(apply(model_zpostDS[,c(paste0("beta.",1:ncol(X_z)),DS_var)], 2,table_pci))
	rownames(table_DS) <-  c(beta_names,DS_var2)	
	table_DS[,2]<-ifelse(grepl("sigma|nu",rownames(table_DS)),"-",table_DS[,2])

	table_N <- t(apply(model_zpostN[,c(paste0("beta.",1:ncol(X_z)),N_var)], 2,table_pci))
	rownames(table_N) <-  c(beta_names,N_var2)	
	table_N[,2]<-ifelse(grepl("sigma|nu",rownames(table_N)),"-",table_N[,2])

	colnames(table_A) <- colnames(table_DS) <- colnames(table_N) <- col_names

	model_z_table_A[[trait]] <- table_A
	model_z_table_DS[[trait]] <- table_DS
	model_z_table_N[[trait]] <- table_N

}

save(
model_w_table,
model_z_table_A,
model_z_table_DS,
model_z_table_N, file=paste0(data_wd,"table_data.Rdata"),version=2)


library(knitr)
knitr::kable(model_w_table)

# rownames(model_w_table$tarsus_mm[[1]])
knitr::kable(model_z_table_DS[[1]])




