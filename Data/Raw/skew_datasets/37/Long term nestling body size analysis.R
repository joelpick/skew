#Are nestlings fledging at a lower weight than they did in the past? 


library(tidyverse)
library(lme4)
library(lmerTest)


dat <- read.csv("Long term nestling body size data.csv")

#Has nestling body mass growth (10-16 day old nestlings) changed over time?
dat2 <- dat %>% filter(mass <30 ) %>% group_by(nestlingID) %>% slice(which.min(abs(age)-12))

mod <- lmer(mass~age*year2 + (1|nestID), data=dat2, na.action="na.fail")

plot(mod)
hist(resid(mod))
plot(resid(mod)~dat2$age)
plot(resid(mod)~dat2$year2)
summary(mod)

anova(mod, test="F")

mam <- lmer(mass~age*year2 + (1|nestID), data=dat2, na.action="na.fail")
summary(mam)





#Has nestling wing chord growth (10-16 day old nestlings) changed over time?
dat3 <- dat %>% filter( ninprim<80 ) %>% group_by(nestlingID) %>% slice(which.min(abs(age)-12))

wmod <- lmer(ninprim~age*year2 + (1|nestID), data=dat3, na.action="na.fail")
plot(wmod)
hist(resid(wmod))
plot(resid(wmod)~dat3$year)
plot(resid(wmod)~dat3$age)
#looks good. Although there is a fairly large chunk of time where we are missing measurements
summary(wmod)
dredge(wmod)
anova(wmod)

wmam <- lmer(ninprim~age*year2 + (1|nestID), data=dat3, na.action="na.fail")
summary(wmam)
anova(wmam)


