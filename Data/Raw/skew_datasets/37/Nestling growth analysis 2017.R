#Nestling growth based on local weather conditions analysis

library(tidyverse)
library(lme4) #for doing the mixed effect models
library(nlme)
library(lmerTest)

# ###############
# #Load in the data
dat <- read.csv("Nestling growth data 2017.csv")


#Calculate out the residual body mass (first must decide on best model to show mass growth across age)

M.lm <- gls(Mass ~ poly(Age, 3), data=dat)
M.gls1 <- gls(Mass ~ poly(Age, 3), weights = varFixed(~Age), data=dat)
M.gls2 <- gls(Mass ~ poly(Age, 3), weights = varIdent(~1|as.factor(Age)), data=dat)

AICc(M.lm, M.gls1, M.gls2)
anova(M.lm, M.gls1, M.gls2)


summary(M.gls1)
plot(M.lm)
plot(M.gls1)
plot(M.gls2)

#Should we use the polynomial or a linear? 
M.gls1_3 <- gls(Mass ~ poly(Age, 3), weights = varFixed(~Age), data=dat)
M.gls1_2 <- gls(Mass ~ poly(Age, 2), weights = varFixed(~Age), data=dat)
M.gls1_1 <- gls(Mass ~ poly(Age, 1), weights = varFixed(~Age), data=dat)
plot(M.gls1_3)
plot(M.gls1_2)
plot(M.gls1_1)

fit1 <- nls(Mass ~ SSlogis(Age, Asym, xmid, scal), data = dat)
fit2 <- nls(Mass ~ SSlogis(Age, Asym, xmid, scal), weights = Age, data = dat)
plot(fit1)
plot(fit2)

AICc(M.gls1_1, M.gls1_2, M.gls1_3, M.gls1_4, fit1, fit2)
#OK We will use the 3rd order with varFixed(Age)


mamMass <- gls(Mass ~ poly(Age, 3), weights = varFixed(~Age), data=dat)
anova(mamMass)
dat$ResidMass <- resid(mamMass)

### Scale residual body mass so each nesting is mean centered. 
dat <- dat %>% group_by(NestlingID) %>% mutate(ResidMass_scaled=ResidMass-mean(ResidMass))


########################################################################
###How do local weather conditions influence nestling body mass?

##Max temperature the prior 3 days
mod1 <- lmer(ResidMass_scaled ~ MaxTemp3day*ThermoReg + (1|NestID/NestlingID), data=dat, REML=FALSE)
plot(mod1) #this is OK
hist(resid(mod1)) 
shapiro.test(resid(mod1))
plot(resid(mod1)~dat$MaxTemp3day)
plot(resid(mod1)~dat$ThermoReg)
plot(resid(mod1)~dat$NestID)

summary(mod1)
#Don't need nestling ID as a random effect

mod2 <- lmer(ResidMass_scaled ~ MaxTemp3day*ThermoReg + (1|NestID), data=dat, REML=FALSE)
summary(mod2)
#Do not need nest ID as a random effect


mod3 <- lm(ResidMass_scaled ~ MaxTemp3day*ThermoReg, data=dat, na.action = "na.fail")

anova(mod3)
dredge(mod3)

mam_maxtemp <- lm(ResidMass_scaled ~ MaxTemp3day*ThermoReg, data=dat)
summary(mam_maxtemp)
anova(mam_maxtemp)


##########################################################
#############Does whether is rains or not predict residual mass?
mod1 <- lmer(ResidMass_scaled ~ TotalRainFall3day*ThermoReg + (1|NestID/NestlingID), data=dat, REML=FALSE)
plot(mod1) #this is OK
hist(resid(mod1)) 
shapiro.test(resid(mod1)) 
plot(resid(mod1)~dat$TotalRainFall3day)
plot(resid(mod1)~dat$ThermoReg)
plot(resid(mod1)~dat$NestID)

summary(mod1)
#Can drop the random nestling ID effect

mod2 <- lmer(ResidMass_scaled ~ TotalRainFall3day*ThermoReg + (1|NestID), data=dat, REML=FALSE)
summary(mod2)
#Should not keep the nestID effect

mod3 <- lm(ResidMass_scaled ~ TotalRainFall3day*ThermoReg, data=dat)


dredge(mod3)
anova(mod3)
#Best model is the full model. However, with delta=1.97 is the much simpler just total precipitation

mam_rain <- lm(ResidMass_scaled ~ TotalRainFall3day*ThermoReg, data=dat)
summary(mam_rain)
anova(mam_rain)

##########################################################
#############Does meanwindspeed predict residual mass?
mod1 <- lmer(ResidMass_scaled ~ MeanWindspeed3day*ThermoReg + (1|NestID/NestlingID), data=dat, REML=FALSE)
plot(mod1) 
hist(resid(mod1)) . 
shapiro.test(resid(mod1)) 
plot(resid(mod1)~dat3$MeanWindspeed3day)
plot(resid(mod1)~dat3$ThermoReg)
plot(resid(mod1)~dat3$NestID)

summary(mod1)
#Good to drop the nestling ID random effect

mod2 <- lmer(ResidMass_scaled ~ MeanWindspeed3day*ThermoReg + (1|NestID), data=dat, REML=FALSE)
summary(mod2)
#don't need to keep the NestID random effect

mod3 <- lm(ResidMass_scaled ~ MeanWindspeed3day*ThermoReg, data=dat)


dredge(mod3)
anova(mod3)

mam_windspeed <- lm(ResidMass_scaled ~ MeanWindspeed3day*ThermoReg , data=dat )
summary(mam_windspeed)


##Rank all 3 of the models based on AICc and R^2


A <- AICc(mam_maxtemp, mam_rain, mam_windspeed)
A$delta <- A$AICc-min(A$AICc)
A
#Based on AICc, rainfall is EASILLY the best predictor. 
anova(mam_maxtemp)
summary(mam_rain)
summary(mam_maxtemp)