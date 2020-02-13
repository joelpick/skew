rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
# source(paste0(wd,"R/functions.R"))

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

tBIRDS <- read.csv(paste0(wd,"Data/Raw/tBIRDS.csv"))
head(tBIRDS)
adults <- unique(subset(read.csv(paste0(wd,"Data/Raw/tMORPH.csv")), agecode%in%c(5,6))$bird_id)

#remove repeated measures and those without a sex
THBW<- subset(THBW,!duplicated(bird_id) & !is.na(sex) & !is.na(tarsus_mm) & !is.na(headbill_mm) & !is.na(wing_mm) & !is.na(weight_g))

THBW[THBW$hatch_day==6,"hatch_day"]<-3

#sex, day of hatching within the nest, year, clutch size, male presence, nest hatch date

THBW$male_present <- as.factor(THBW$male_present)
THBW$year <- as.factor(THBW$year)
THBW$sex <- as.factor(THBW$sex)
THBW$hatch_day <- as.factor(THBW$hatch_day)
THBW$clutch_sizeC <- scale(THBW$clutch_size, scale=FALSE)
THBW$nest_hatch_dateC <- scale(THBW$nest_hatch_date, scale=FALSE)

THBW$tarsus_mmC <- scale(THBW$tarsus_mm, scale=FALSE)
THBW$tarsus_mmC2 <- THBW$tarsus_mmC^2

THBW$headbill_mmC <- scale(THBW$headbill_mm, scale=FALSE)
THBW$headbill_mmC2 <- THBW$headbill_mmC^2

THBW$wing_mmC <- scale(THBW$wing_mm, scale=FALSE)
THBW$wing_mmC2 <- THBW$wing_mmC^2

THBW$weight_gC <- scale(THBW$weight_g, scale=FALSE)
THBW$weight_gC2 <- THBW$weight_gC^2

THBW$recruit <- THBW$bird_id %in% adults
THBW$fledge <- is.na(tBIRDS$max_death_date[match(THBW$bird_id , tBIRDS$bird_id)])

THBW2<-THBW
colnames(THBW2)[which(colnames(THBW2)=="fledge")]<-"survival"
colnames(THBW2)[which(colnames(THBW2)=="recruit")]<-"age"
THBW2$age<-15

THBW3<-THBW[which(THBW$fledge),]
colnames(THBW3)[which(colnames(THBW3)=="recruit")]<-"survival"
colnames(THBW3)[which(colnames(THBW3)=="fledge")]<-"age"
THBW3$age<-25

THBW2<-rbind(THBW2, THBW3)
THBW2$age<-as.factor(THBW2$age)


prior_bv <- list( R = list(V = 1, fix = 1), 
			   G = list(G1 = list(V = diag(2), nu = 2, alpha.mu=c(0,0), alpha.V=diag(2)*100), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002))
			 )
it_scale <- 50


mod_tarsus_bv <- MCMCglmm(survival ~ -1+at.level(age,1):(age+tarsus_mmC + tarsus_mmC2 + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC)+at.level(age,2):(age+tarsus_mmC + tarsus_mmC2 + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC), random=~us(age):nest, data=THBW2, prior=prior_bv, family="threshold", nitt=13000*it_scale, thin=5*it_scale, burnin=3000*it_scale)


mod_weight_bv <- MCMCglmm(survival ~ -1+at.level(age,1):(age+weight_gC + weight_gC2 + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC)+at.level(age,2):(age+weight_gC + weight_gC2 + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC), random=~us(age):nest, data=THBW2, prior=prior_bv, family="threshold", nitt=13000*it_scale, thin=5*it_scale, burnin=3000*it_scale)

save(THBW, mod_tarsus_bv, mod_weight_bv , file = paste0(wd,"Data/Intermediate/day15_survival_models_bv.Rdata"))




### splines

mod_tarsus_bv_splines <- MCMCglmm(survival ~ -1+at.level(age,1):(age+tarsus_mmC  + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC)+at.level(age,2):(age+tarsus_mmC + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC), random=~us(age):nest + idv(at.level(age,1):spl(tarsus_mmC,k=30))+idv(at.level(age,2):spl(tarsus_mmC,k=30)), data=THBW2,pr=TRUE, prior=prior_bv, family="threshold", nitt=13000*it_scale, thin=5*it_scale, burnin=3000*it_scale)

save(THBW2, mod_tarsus_bv_splines, file = paste0(wd,"Data/Intermediate/day15_survival_models_bv_spline.Rdata"))


mod_weight_bv_splines <- MCMCglmm(survival ~ -1+at.level(age,1):(age+weight_g  + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC)+at.level(age,2):(age+weight_g + sex + male_present + year + hatch_day + clutch_sizeC + nest_hatch_dateC), random=~us(age):nest + idv(at.level(age,1):spl(weight_g,k=30))+idv(at.level(age,2):spl(weight_g,k=30)), data=THBW2,pr=TRUE, prior=prior_bv, family="threshold", nitt=13000*it_scale, thin=5*it_scale, burnin=3000*it_scale)

save(THBW2, mod_weight_bv_splines, file = paste0(wd,"Data/Intermediate/day15_survival_models_bv_splineW.Rdata"))



##############################

