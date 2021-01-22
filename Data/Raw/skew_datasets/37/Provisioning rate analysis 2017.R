#Is feeding rate affected by daily windspeed maximum temperature or whether it had rained in the previous 3 days?
#(correcting for the number of nestlings)
library(tidyverse)

feeding <- read.csv("Provisioning rate data 2017.csv")
feeding_female <- feeding %>% filter(!is.na(ffv))
feeding_male <- feeding %>% filter(!is.na(mfv))


###females

#mean windspeed
mod <- lm(ffv~meanwindspeed*NestlingsAlive, data=feeding_female)
plot(mod)
hist(resid(mod))
shapiro.test(resid(mod))
plot(resid(mod)~feeding_female$meanwindspeed)
plot(resid(mod)~feeding_female$NestlingsAlive)

anova(mod)
#Nestlings alone is the better model

#Max temperature
mod <- lm(ffv~MaxTemp*NestlingsAlive, data=feeding_female)
plot(mod)
hist(resid(mod))
shapiro.test(resid(mod))
plot(resid(mod)~feeding_female$MaxTemp)
plot(resid(mod)~feeding_female$NestlingsAlive)

anova(mod) 
#Nestlings alone is the better model

#Rain in the previous 3 days
mod <- lm(ffv ~NestlingsAlive*TotalRain_3day2, data=feeding_female, na.action="na.fail")
plot(mod)
hist(resid(mod))
shapiro.test(resid(mod))
plot(resid(mod)~feeding_female$NestlingsAlive)
plot(resid(mod)~as.factor(feeding_female$TotalRain_3day2))
anova(mod)

mam_rain3day <- lm(ffv ~NestlingsAlive+TotalRain_3day2, data=feeding_female)
anova(mam_rain3day)
summary(mam_rain3day)


#males
#Mean windspeed
mod <- lm(mfv~meanwindspeed*NestlingsAlive, data=feeding_male, na.action="na.fail")
plot(mod)
hist(resid(mod))
plot(resid(mod)~feeding_male$meanwindspeed)
summary(mod)
anova(mod, test="F") 
#Nope not at all. 

#Max Temp
mod <- lm(mfv~MaxTemp*NestlingsAlive,  data=feeding_male, na.action="na.fail")
plot(mod)
hist(resid(mod))
plot(resid(mod)~feeding_male$MaxTemp)
anova(mod, test="F")



#Do females get to slack a bit if the males work harder? 
feeding2 <- feeding %>% filter(!is.na(mfv) & !is.na(ffv))
mod <- lm(mfv ~ffv, data=feeding2)
plot(mod)
hist(resid(mod))
plot(resid(mod)~feeding2$mfv)
anova(mod)
summary(mod)
#Nope females who have harder working males are also working harder!

