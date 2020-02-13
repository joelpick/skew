rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
# source(paste0(wd,"functions.R"))

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

load(paste0(wd,"Data/Intermediate/dam_sire_egg.Rdata"))


tBIRDS <- read.csv(paste0(wd,"Data/Raw/tBIRDS.csv"))
head(tBIRDS)
adults <- unique(subset(read.csv(paste0(wd,"Data/Raw/tMORPH.csv")), agecode%in%c(5,6))$bird_id)

#remove repeated measures and those without a sex
THBW<- subset(THBW,!duplicated(bird_id) & !is.na(sex) & !is.na(tarsus_mm) & !is.na(headbill_mm) & !is.na(wing_mm) & !is.na(weight_g))

THBW$recruit <- THBW$bird_id %in% adults
THBW$fledge <- is.na(tBIRDS$max_death_date[match(THBW$bird_id , tBIRDS$bird_id)])

S <- mean(THBW$tarsus_mm[THBW$recruit]) - mean(THBW$tarsus_mm)

h2 <- modT$vars_mean["mod_egg_A","animal"]/sum(c(modT$vars_mean["mod_egg_A",],modT$fixed_var["mod_egg_A"]))

S*h2