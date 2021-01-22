#Is adult body size changing through time? 

#Packages Needed
library(lme4)
library(tidyverse)
library(lmerTest)

adult <- read.csv("Long term adult body size data.csv")


#Has female body mass changed throught time? 
adult_Fm <- adult %>% filter(sex=="F" & !is.na(mass))
adult_Mm <- adult %>% filter(sex=="M" & diff>17 & !is.na(mass))

#Testing first through third order polynomials
mod <- lmer(mass ~ poly(diff, 1)*year2 + (1|band), data=adult_Fm, na.action="na.fail")
summary(mod)
mod_2 <- lmer(mass ~ poly(diff, 2)*year2 + (1|band), data=adult_Fm, na.action="na.fail")
summary(mod_2)
mod_3 <- lmer(mass ~ poly(diff, 3)*year2 + (1|band), data=adult_Fm, na.action="na.fail")
AICc(mod, mod_2, mod_3)

#Definitely need to retain the random effect of bird. 
#Chose to use the 2rd order polynomial. 

mod_F <- lmer(mass ~ poly(diff, 2)*year2 + (1|band), data=adult_Fm, na.action="na.fail")
plot(mod_F)
hist(resid(mod_F))
plot(resid(mod_F)~adult_F$diff)
plot(resid(mod_F)~adult_F$year)
plot(resid(mod_F)~adult_F$Period)

anova(mod_F, test="F")

mam_F <- lmer(mass ~ poly(diff, 2)*year2 + (1|band), data=adult_Fm, na.action="na.fail")
summary(mam_F)
anova(mam_F, test="F")


##################
#Has male body mass changed throught time? 
#Testing first through third order polynomials
mod <- lmer(mass ~ poly(diff, 1)*year2 + (1|band), data=adult_Mm, na.action="na.fail")
summary(mod)
mod_2 <- lmer(mass ~ poly(diff, 2)*year2 + (1|band), data=adult_Mm, na.action="na.fail")
summary(mod_2)
mod_3 <- lmer(mass ~ poly(diff, 3)*year2 + (1|band), data=adult_Mm, na.action="na.fail")
AICc(mod, mod_2, mod_3)
#SHould retain the random effect and go with a linear
#relationship.


mod_M <- lmer(mass ~ diff*year2+ (1|band), data=adult_Mm, na.action="na.fail")
plot(mod_M)
hist(resid(mod_M))
shapiro.test(resid(mod_M))
plot(resid(mod_M)~adult_Mm$diff)
plot(resid(mod_M)~adult_Mm$year)

anova(mod_M)

mam_M <- lmer(mass ~ diff*year2 + (1|band), data=adult3_Mm, na.action="na.fail")
summary(mam_M)
anova(mam_M)




#Has wing chord also decreased for either sex? 

adult_Fw <- adult %>% filter(sex=="F" & !is.na(wingChord)& wingChord>80 & wingChord<135 & diff>-5 & diff<30)
adult_Mw<- adult %>% filter(sex=="M" &  !is.na(wingChord)& wingChord>80 & wingChord<135 & diff>-5 & diff<30)

#Females
mod_fw <- lmer(wingChord ~ year + (1|band), data=adult_Fw, na.action="na.fail")
plot(mod_fw)
shapiro.test(resid(mod_fw))
hist(resid(mod_fw))
plot(resid(mod_fw)~adult_Fw$year)
summary(mod_fw) #Need to retain band ID-- you are measured multiple times and you'll be about the same
anova(mod_fw)
#global is best model

#Males
mod_mw <- lmer(wingChord ~ year + (1|band), data=adult_Mw, na.action="na.fail")
plot(mod_mw)
shapiro.test(resid(mod_mw))
hist(resid(mod_mw))
plot(resid(mod_mw)~adult_Mw$year)
summary(mod_mw) #Need to retain band ID-- you are measured multiple times and you'll be about the same
anova(mod_mw)
#global is best model


