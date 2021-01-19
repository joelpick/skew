######################
#
#---- This script makes the chick body size dataset use in all subsequent analyses
#
######################

rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

library(MCMCglmm)
library(MasterBayes)

data_wd <- paste0(wd,"Data/Raw/")

source(paste0(wd,"R/functions.R"))

tBIRDS <- read.csv(paste0(data_wd,"tBIRDS.csv"))
tPED <- read.csv(paste0(data_wd,"tPED.csv"))
tMORPH <- read.csv(paste0(data_wd,"tMORPH.csv"))
tEGGS <- read.csv(paste0(data_wd,"tEGGS.csv"))
tNEST_SURVEY <- read.csv(paste0(data_wd,"tNEST_SURVEY.csv"))

### --------------------------
## Make pedigree
### --------------------------

## there are two pedigrees. They differ in how they they assign sires in colony, one assumes the females only mate with one unknown male ("M"), and one assumes that she can mate multiply ("P")

tPED_M <- subset(tPED, substr(nest_orig,1,2)=="09" | ped_type=="DSCM")[,c("bird_id","dam","sire","nest_orig")]
tPED_P <- subset(tPED, substr(nest_orig,1,2)=="09" | ped_type=="DSCP")[,c("bird_id","dam","sire")]
names(tPED_M)<-c("bird_id","dam_M","sire_M","nest_orig")
names(tPED_P)<-c("bird_id","dam_P","sire_P")
tPED_use<-merge(tPED_P,tPED_M)
n_unique(tPED_use$bird_id)==nrow(tPED_use)


### --------------------------
## Data exclusion 
### --------------------------


## extract year
tMORPH$year <- substr(tMORPH$date,1,4)
tail(tMORPH)

## exclude 2009,2010, 2019 and 2020, measures without a nest (winter catches), those not processed by Jarrod, and adult measurements
tMORPH1 <- subset(tMORPH, 
	year %in% 2011:2018
	& !is.na(nest) 
	& processed_by%in%c("JDH",NA)
	& is.na(agecode)
	)
nrow(tMORPH1)

## exclude weird nests (most found when chicks had already hatched) and those with no egg measurements
exclude <- nest_to_exclude(tMORPH1,tNEST_SURVEY)
meanEGGS <- aggregate(egg_weight~nest_orig, tEGGS, function(x) mean(x, na.rm=TRUE))
no_egg <- unique(tMORPH1$nest[!tMORPH1$nest%in%meanEGGS$nest_orig])

# exclude[!exclude %in% no_egg]
# subset(THBW,nest_orig%in%no_egg[!no_egg %in% exclude])

tMORPH2 <- subset(tMORPH1, !nest %in% unique(c(exclude,no_egg)))


### --------------------------
## add other variables in to data
### --------------------------

## convert time into hours
tMORPH2$time_hr <- as.numeric(strptime(
	ifelse(nchar(tMORPH2$time)==4, paste0(0,tMORPH2$time),
	ifelse(nchar(tMORPH2$time)==5,tMORPH2$time,
	ifelse(nchar(tMORPH2$time)>5,substr(tMORPH2$time,1,5),NA
	))),"%H:%M")- strptime(Sys.Date(), "%F"))

# add male presence, clutch size, origin/parentage and egg weight
mp <- male_presence(tNEST_SURVEY,tMORPH)
cs <- clutch_size(tEGGS)
tMORPH2$clutch_size <- cs[match(tMORPH2$nest,cs$nest),"clutch_size"]
tMORPH2$male_present <- mp[match(tMORPH2$nest,mp$nest),"male_present"]

tPED_use$egg_weight <- meanEGGS[match(tPED_use$nest_orig,meanEGGS$nest_orig),"egg_weight"] 

tMORPH3 <- merge(tMORPH2,tPED_use, all.x=TRUE)
nrow(tMORPH2)
nrow(tMORPH3)

head(tPED_use)

## add chick age and nest hatch date
tMORPH4 <- groupFunc_dataFrame(tMORPH3,"nest", function(x){
	x$nest_age <- days_from_minimum(x$date)
	x$nest_hatch_date <- april_hatch_date(x$date)
	return(x)
	})

## add chick hatch day
tMORPH5 <- groupFunc_dataFrame(tMORPH4,"bird_id", function(x){ 
	x$hatch_day <- as.character(min(x$nest_age))
	return(x)
	})

## dataframe for day 15 measurements (Tarsus, Head-Bill and Wing)
THBW <- subset(tMORPH5, nest_age==15)[,!colnames(tMORPH5) %in% c("agecode","ugcs","uas","morph_sex","processed_by","ringing_id","colour_id","time","date")]


### use information of chicks caught as breeding adults to update sexing information 
### - sexing adults in hand is more reliable than molecular sexing
adults <- subset(tMORPH, agecode%in%c(5,6) & !is.na(nest) )
## this selects only adults that were caught when breeding

## check any adults recorded with multiple sexes
# any(sapply(split(adults, adults$bird_id), function(x) length(unique(x$morph_sex))>1 ))

## unique adult sexes from morphological data
adult_sex <- t(sapply(split(adults, adults$bird_id), function(x) c(x$bird_id[1],unique(x$morph_sex), min(x$year)) ))

## if sexed morphologically as adult, take that over molecular sexing of chicks
THBW_adult_sex <- adult_sex[,2][match(THBW$bird_id,adult_sex[,1])]
THBW_chick_sex <- tBIRDS$sex[match(THBW$bird_id,tBIRDS$bird_id)]

THBW$sex <- ifelse(is.na(THBW_adult_sex),THBW_chick_sex,THBW_adult_sex)


## get recruitment data 
adult_ids <- unique(adults$bird_id)

# recruit are those chicks we have caught as adult *breeders*
THBW$recruit <- THBW$bird_id %in% adult_ids
THBW$fledge <- is.na(tBIRDS$max_death_date[match(THBW$bird_id , tBIRDS$bird_id)])


## work out how many individuals dont have sex or morphology measurement
THBW2 <- subset(THBW, !is.na(tarsus_mm)& !is.na(headbill_mm)& !is.na(wing_mm)& !is.na(weight_g))
n_unique(THBW$bird_id)-n_unique(THBW2$bird_id)
n_unique(subset(THBW2,is.na(sex))$bird_id) 

## exclude those with no sex or a missing morphological measurement
THBW<- subset(THBW,!is.na(sex) & !is.na(tarsus_mm) & !is.na(headbill_mm) & !is.na(wing_mm) & !is.na(weight_g))

# number of birds
n_unique(THBW$bird_id)

# number of repeated measurements
nrow(THBW)-n_unique(THBW$bird_id)

#number of nests
n_unique(THBW$nest)




# assign individuals with hatch day 6, hatch day 3
THBW[THBW$hatch_day==6,"hatch_day"]<-3

# make factors into factors
THBW$male_present <- as.factor(THBW$male_present)
THBW$year <- as.factor(THBW$year)
THBW$sex <- as.factor(THBW$sex)
THBW$hatch_day <- as.factor(THBW$hatch_day)


## continuous variables centered
THBW$timeC <- scale(THBW$time_hr, scale=FALSE)
THBW$clutch_sizeC <- scale(THBW$clutch_size, scale=FALSE)
THBW$nest_hatch_dateC <- scale(THBW$nest_hatch_date, scale=FALSE)
THBW$egg_weightC <- scale(THBW$egg_weight, scale=FALSE)

THBW$tarsus_mmC <- scale(THBW$tarsus_mm, scale=FALSE)
THBW$tarsus_mmC2 <- THBW$tarsus_mmC^2

THBW$headbill_mmC <- scale(THBW$headbill_mm, scale=FALSE)
THBW$headbill_mmC2 <- THBW$headbill_mmC^2

THBW$wing_mmC <- scale(THBW$wing_mm, scale=FALSE)
THBW$wing_mmC2 <- THBW$wing_mmC^2

THBW$weight_gC <- scale(THBW$weight_g, scale=FALSE)
THBW$weight_gC2 <- THBW$weight_gC^2

dam_sire <- unique(c(THBW$dam_P,THBW$sire_P))
THBW$dam <- as.numeric(factor(THBW$dam_P, levels=dam_sire))
THBW$sire <- as.numeric(factor(THBW$sire_P, levels=dam_sire))


## dataframe with no repeated measures
THBW_noRep<- subset(THBW,!duplicated(bird_id))


## prune pedigree
ped <- insertPed(prunePed(orderPed(tPED_use[1:3]), unique(THBW$bird_id), make.base = TRUE))
names(ped) <- c("animal","dam", "sire")

#save all the data 
save(THBW, THBW_noRep, ped, file= paste0(wd,"Data/Intermediate/chick_data.Rdata"))




