rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

library(MCMCglmm)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

load(paste0(wd,"Data/Intermediate/chick_data.Rdata"))

tMORPH <- read.csv(paste0(wd,"Data/Raw/tMORPH.csv"))
head(tMORPH)
tNESTBOX <- read.csv(paste0(wd,"Data/Raw/tNESTBOX.csv"))
head(tNESTBOX)
sexes <- read.csv(paste0(wd,"Data/Raw/tBIRDS.csv"))[,c("bird_id","sex")]
names(sexes)[2]<-"genetic_sex"
### are immigrants different from

tMORPH2 <- subset(tMORPH, !is.na(nest))
tMORPH2$year <- substr(tMORPH2$nest,1,2)
tMORPH2$nestbox <- substring(tMORPH2$nest,4)
tMORPH2$nestbox <- sub("b","",tMORPH2$nestbox)
tMORPH2 <- tMORPH2[order(tMORPH2$date),]

dd <- groupFunc_dataFrame(tMORPH2,"bird_id",function(x){
	ch <- subset(x, is.na(agecode))
	ad <- subset(x, !is.na(agecode))
	
	cols<-c("bird_id","year","nestbox","nest","agecode","tarsus_mm","morph_sex")
	out <- if(nrow(ad)>0) {
		if(nrow(ch)>0) {
			cbind(rbind(ch[nrow(ch),cols],ad[1,cols]),type = "recruit")
		}else{
			c(ad[1,cols],type = "immigrant")
		}
	}else{NULL}
	return(out)
})
head(dd,20)

adults <- subset(dd,!is.na(agecode) & !year%in%c("09") )
table(adults$year,adults$type)

boxplot(tarsus_mm~type+morph_sex,adults, ylab="Adult tarsus length (mm)")
library(lme4)
summary(lmer(tarsus_mm~type+morph_sex+(1|year),adults))

table(adults$agecode,adults$type)
table(adults$type,adults$morph_sex)
table(adults$type,adults%year,adults$morph_sex)


recruits <- subset(dd,type=="recruit")
recruits2 <- merge(recruits,tNESTBOX)
recruits2 <- recruits2[order(recruits2$year),]


# recruits <- subset(recruits, !bird_id%in%recruits[recruits$tarsus_mm<15,"bird_id"])
# recruits$age <- ifelse(is.na(recruits$agecode),"J",recruits$agecode)
# boxplot(tarsus_mm~age,recruits)
# t.test(subset(recruits,is.na(agecode))$tarsus_mm,subset(recruits,!is.na(agecode))$tarsus_mm,paired=TRUE)
# mean(subset(recruits,is.na(agecode))$tarsus_mm-subset(recruits,!is.na(agecode))$tarsus_mm,na.rm=TRUE)

recruits3 <- groupFunc_dataFrame(recruits2,"bird_id",function(x){
	diffY <-abs(x$lat[1]-x$lat[2])
	diffX <-abs(x$long[1]-x$long[2])
	x$distance <- sqrt(diffY^2 + diffX^2)
	cbind(x[2,c("bird_id","morph_sex","distance")],x[1,c("nest","tarsus_mm","year")])
})
head(recruits3)
table(table(recruits3$nest))
plot(table(table(recruits3$nest)))


## nest of rearing and ## year

plot(distance~tarsus_mm,recruits3, pch=19, col=c(1,2)[as.factor(recruits3$morph_sex)])
summary(lm(distance~tarsus_mm+morph_sex,recruits3))

subset(tMORPH,bird_id==recruits3[recruits3$tarsus_mm<15,"bird_id"])

par(mfrow=c(2,1))
hist(recruits3[recruits3$morph_sex=="M","distance"], breaks=20, xlim=c(0,0.06))
hist(recruits3[recruits3$morph_sex=="F","distance"], breaks=80, xlim=c(0,0.06))

site_moves <- subset(tMORPH, bird_id %in% recruits3[recruits3$distance>0.03,"bird_id"])
site_moves[order(site_moves$bird_id),]





dd2 <- groupFunc_dataFrame(tMORPH2,"bird_id",function(x){
	ch <- subset(x, is.na(agecode))
	ad <- subset(x, !is.na(agecode))
	
	cols<-c("bird_id","year","nestbox","nest","agecode","tarsus_mm","morph_sex")
	out <- if(nrow(ad)>0) {
		if(nrow(ch)>0) {
			cbind(ad[1,cols],type = "recruit")
		}else{
			cbind(ad[1,cols],type = "immigrant")
		}
	}else{
		cbind(ch[nrow(ch),cols],type = "juvenile")
	}
	return(out)
})
head(dd2)
dd3 <- merge(subset(dd2, !is.na(tarsus_mm)),sexes)
head(dd3)
#dd3$morph_sex <- ifelse(is.na(dd3$morph_sex),"U",dd3$morph_sex)
dd3$sex <- ifelse(is.na(dd3$morph_sex),dd3$genetic_sex,dd3$morph_sex)
boxplot(tarsus_mm~type+sex,dd3, ylab="Adult tarsus length (mm)")
aggregate(tarsus_mm~type+sex,dd3, mean)

dd3$morph_sex <- ifelse(is.na(dd3$morph_sex),"U",dd3$morph_sex)



dd2 <- groupFunc_dataFrame(tMORPH2,"bird_id",function(x){
	ch <- subset(x, is.na(agecode))
	ad <- subset(x, !is.na(agecode))
	
	cols<-c("bird_id","year","nestbox","nest","agecode","tarsus_mm","morph_sex")
	out <- if(nrow(ad)>0) {
		if(nrow(ch)>0) {
			cbind(ad[1,cols],type = "recruit")
		}else{
			cbind(ad[1,cols],type = "immigrant")
		}
	}else{
		cbind(ch[nrow(ch),cols],type = "juvenile")
	}
	return(out)
})