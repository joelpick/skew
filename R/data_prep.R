rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

data_wd <- paste0(wd,"Data/Raw/")

source(paste0(wd,"R/functions.R"))

#download_bluetit("tNEST_SURVEY",data_wd)

tBIRDS <- read.csv(paste0(data_wd,"tBIRDS.csv"))
tPED <- read.csv(paste0(data_wd,"tPED.csv"))
tMORPH <- read.csv(paste0(data_wd,"tMORPH.csv"))
tEGGS <- read.csv(paste0(data_wd,"tEGGS.csv"))
tNEST_SURVEY <- read.csv(paste0(data_wd,"tNEST_SURVEY.csv"))

tPED_M <- subset(tPED, substr(nest_orig,1,2)=="09" | ped_type=="DSCM")[,c("bird_id","dam","sire","nest_orig")]
tPED_P <- subset(tPED, substr(nest_orig,1,2)=="09" | ped_type=="DSCP")[,c("bird_id","dam","sire")]
names(tPED_M)<-c("bird_id","dam_M","sire_M","nest_orig")
names(tPED_P)<-c("bird_id","dam_P","sire_P")
tPED_use<-merge(tPED_P,tPED_M)
n_unique(tPED_use$bird_id)==nrow(tPED_use)

# identical(tPED_use$sire_M,tPED_use$sire_P)
# any(tPED_use$sire_M!=tPED_use$sire_P)

tMORPH$year <- substr(tMORPH$date,1,4)
tail(tMORPH)
## convert time into hours
tMORPH$time_hr <- as.numeric(strptime(
	ifelse(nchar(tMORPH$time)==4, paste0(0,tMORPH$time),
	ifelse(nchar(tMORPH$time)==5,tMORPH$time,
	ifelse(nchar(tMORPH$time)>5,substr(tMORPH$time,1,5),NA
	))),"%H:%M")- strptime(Sys.Date(), "%F"))

## exclude 2009, measures without a nest, and those not processed by jarrod
tMORPH1 <- subset(tMORPH, !year%in%c("2009") & !is.na(nest) & processed_by%in%c("JDH",NA))

## exclude weird nests and adult measurements
exclude <- nest_to_exclude(tMORPH1,tNEST_SURVEY)
tMORPH2 <- subset(tMORPH1, !nest %in% exclude & is.na(agecode))

#male presence, clutch size, origin/parentage and sex
mp <- male_presence(tNEST_SURVEY,tMORPH)
cs <- clutch_size(tEGGS)

tMORPH3 <- merge(merge(merge(merge(tMORPH2,mp),cs),tBIRDS[,c("bird_id","sex")]),tPED_use, all.x=TRUE)

## add and chick age nest hatch date
tMORPH4 <- groupFunc_dataFrame(tMORPH3,"nest", function(x){
	x$nest_age <- days_from_minimum(x$date)
	x$nest_hatch_date <- april_hatch_date(x$date)
	return(x)
	})

tMORPH4 <- tMORPH4[,!colnames(tMORPH4) %in% c("agecode","ugcs","uas","morph_sex","processed_by","ringing_id","colour_id","time","date")]

head(tMORPH4)

# when catching has been done
# tMORPH4b <- subset(groupFunc_dataFrame(tMORPH1,"nest", function(x){
# 	x$nest_age <- days_from_minimum(x$date)
# 	x$nest_hatch_date <- april_hatch_date(x$date)
# 	return(x)
# 	}), !is.na(agecode) & !nest %in% exclude)
# table(tMORPH4b$year,tMORPH4b$nest_age)
# par(mfrow=c(3,3))
# for(i in unique(tMORPH4b$year)) {hist(subset(tMORPH4b,year==i)$nest_age, xlim=range(tMORPH4b$nest_age), breaks=seq(-0.5,20.5,1), col="grey", main=i, xlab="nest age", ylim=c(0,max(table(tMORPH4b$year,tMORPH4b$nest_age)
# ))); abline(v=12, col=2)}

## add hatch day and remove repeated day 15 measurements
tMORPH5 <- groupFunc_dataFrame(tMORPH4,"bird_id", function(x){ 
	if(sum(x$nest_age==15)==2) x <- x[-which(x$nest_age==15)[2],]
	x$hatch_day <- as.character(min(x$nest_age))
	return(x)
	})[,!colnames(tMORPH4) %in% c("tarsus_mm","headbill_mm","wing_mm")]


THBW <- groupFunc_dataFrame(tMORPH4,"bird_id", function(x){ 
	x$hatch_day <- as.character(min(x$nest_age))
	x <- subset(x, nest_age==15)
	return(x)
	})

nrow(tMORPH)
nrow(THBW)
#THBW <- subset(tMORPH4, !is.na(tarsus_mm|headbill_mm|wing_mm))
head(THBW) 
n_unique(THBW$bird_id);nrow(THBW)

#nrow(tMORPH4)==nrow(tMORPH5)
#n_unique(tMORPH4$bird_id)==n_unique(tMORPH5$bird_id)
#head(tMORPH5)




## add in egg data

meanEGGS <- aggregate(egg_weight~nest_orig, tEGGS, function(x) mean(x, na.rm=TRUE))

# hist(meanEGGS$egg_weight, breaks=50)
# hist(aggregate(egg_weight~nest_orig,tEGGS,mean)$egg_weight)
# tail(meanEGGS)
# nrow(meanEGGS)
 
tMORPH6 <- subset(merge(tMORPH5, meanEGGS), egg_weight!="NaN")

THBW_egg <- subset(merge(THBW, meanEGGS), egg_weight!="NaN")

tPED_use<-tPED_use[,1:5]
save(tMORPH5, tMORPH6, THBW, THBW_egg, tPED_use, file= paste0(wd,"Data/Intermediate/chick_data.Rdata"))


### SEX
# sum(is.na(tMORPH6$sex))
# nrow(tMORPH6)
# tMORPH6_noSex <- subset(tMORPH6, is.na(tMORPH6$sex) &year!=2018)
# table(tMORPH6_noSex$year)
# length(unique(tMORPH6_noSex$bird_id))
# length(unique(tMORPH6_noSex$nest))
# nrow(tMORPH6_noSex)

###




table(tMORPH5$year)
table(tMORPH6$year)

tm <- merge(tMORPH5, meanEGGS)
nrow(tMORPH5)
nrow(tm)
nrow(tMORPH6)
n_unique(tMORPH5$bird_id)-n_unique(tMORPH6$bird_id)
subset(tMORPH5,year=="2018")$bird_id %in% tBIRDS$bird_id

subset(meanEGGS,grepl("^18", nest_orig))$bird_id %in% tBIRDS$bird_id





nrow(tMORPH6)
head(tMORPH5)
tail(tMORPH5)


hist(subset(tMORPH4,nest_age==9)$weight_g)


resid6 <- resid(lm(weight_g~as.factor(male_present)+hatch_day+nest_hatch_date+as.factor(year),tMORPH5,subset = nest_age==6))
stand_skew(resid6)
stand_skew(subset(tMORPH5,nest_age==6)$weight_g)
# par(mfrow=c(2,1))
# hist(resid6)
# hist(subset(tMORPH4,nest_age==6)$weight_g)


resid9 <- resid(lm(weight_g~as.factor(male_present)+hatch_day+nest_hatch_date+as.factor(year),tMORPH5,subset = nest_age==9))
# par(mfrow=c(2,1))
# hist(subset(tMORPH4,nest_age==9)$weight_g)
# hist(resid9)
stand_skew(resid9)
stand_skew(subset(tMORPH4,nest_age==9)$weight_g)



resid12 <- resid(lm(weight_g~as.factor(male_present)+hatch_day+nest_hatch_date+as.factor(year),tMORPH5,subset = nest_age==12))
stand_skew(resid12)
stand_skew(subset(tMORPH5,nest_age==12)$weight_g)

resid15 <- resid(lm(weight_g~as.factor(male_present)+hatch_day+nest_hatch_date+as.factor(year),tMORPH5,subset = nest_age==15))
#hist(resid15)
stand_skew(resid15)
stand_skew(subset(tMORPH4,nest_age==15)$weight_g)



par(mfrow=c(2,1))
hist(subset(tMORPH4,nest_age==12)$weight_g)
hist(resid12)

subset(tMORPH4,bird_id%in%subset(tMORPH4,nest_age==12 & weight_g<4)$bird_id)




