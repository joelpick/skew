
load("/Users/jpick/Desktop/stanModTarsus_DS20200206_1944.Rdata")
summary(mod_stan_tarsus_dam)$summary[c("sigma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_ind", "alpha_ind", "nu_ind", "sigma_E","lp__"),c(4,1,8,9,10)]
div_trans(mod_stan_tarsus_dam)
pairs(mod_stan_tarsus_dam,pars =c("beta"))

pairs(mod_stan_tarsus_dam,pars =c("sigma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_ind", "alpha_ind", "nu_ind", "sigma_E","lp__"))
rstan::traceplot(mod_stan_tarsus_dam, pars =c("sigma_nest", "gamma_nest", "alpha_nest", "nu_nest", "sigma_dam_sire", "gamma_dam_sire", "alpha_dam_sire", "nu_dam_sire", "sigma_ind", "gamma_ind", "alpha_ind", "nu_ind", "sigma_E","lp__"))

mod<-summary(mod_stan_tarsus_dam)$summary
mod_names<-rownames(mod)
nest_effects <- mod[grep("nest_effects",mod_names),1]
dam_sire_effects <- mod[grep("dam_sire_effects",mod_names),1]
ind_effects <- mod[grep("ind_effects",mod_names),1]
beta <- mod[grep("beta\\[",mod_names),1]

preds_stan <- stan_data_tarsus$X %*% beta + 
	nest_effects[stan_data_tarsus$nest_id] +
	dam_sire_effects[stan_data_tarsus$dam_id] +
	dam_sire_effects[stan_data_tarsus$sire_id] +
	ind_effects[stan_data_tarsus$ind_id]

resid_stan <- THBW_egg$tarsus_mm - preds_stan
hist(preds,breaks=200)
c(nest=mean(nest_effects),dam_sire=mean(dam_sire_effects),resid=mean(ind_effects))
sd(nest_effects)
THBW_egg$tarsus_mm

THBW_egg$nest <- as.factor(THBW_egg$nest)
THBW_egg$bird_id <- as.factor(THBW_egg$bird_id)
mm_levels <- unique(c(THBW_egg$sire_P,THBW_egg$dam_P))
THBW_egg$dam <- factor(THBW_egg$dam_P, levels=mm_levels)
THBW_egg$sire <- factor(THBW_egg$sire_P, levels=mm_levels)

mod_egg_DS <- asreml(
	fixed= tarsus_mm ~ male_present + clutch_sizeC + nest_hatch_dateC  + hatch_day + year + timeC + sex + egg_weightC,  
	random= ~ nest + str(~dam+and(sire), ~us(1):id(dam)) +bird_id,
	rcov = ~idv(units) ,
	trace=FALSE,
	Cfixed=TRUE,
	data=THBW_egg)


resid_asreml <- residuals(mod_egg_DS)
preds_asreml <- THBW_egg$tarsus_mm - resid_asreml
plot(resid_asreml~resid_stan)
plot(preds_asreml~preds_stan)

nest_effects_asreml <- coef(mod_egg_DS, list=TRUE)$nest[,1]
ind_effects_asreml <- coef(mod_egg_DS, list=TRUE)$bird_id[,1]
plot(nest_effects_asreml,nest_effects)
plot(ind_effects_asreml,ind_effects)


tapply(ind_effects_asreml[stan_data_tarsus$ind_id], THBW_egg$hatch_day, mean)
	tapply(ind_effects[stan_data_tarsus$ind_id], THBW_egg$hatch_day, mean)

## predictions
rbind(stan=tapply(preds_stan, THBW_egg$male_present, mean),
asreml=tapply(preds_asreml, THBW_egg$male_present, mean),
observed=tapply(THBW_egg$tarsus_mm, THBW_egg$male_present, mean))

## nest effects
rbind(asreml=tapply(nest_effects_asreml[stan_data_tarsus$nest_id], THBW_egg$male_present, mean),
stan=tapply(nest_effects[stan_data_tarsus$nest_id], THBW_egg$male_present, mean))

## residual effects
rbind(asreml=tapply(ind_effects_asreml[stan_data_tarsus$ind_id], THBW_egg$male_present, mean),
stan=tapply(ind_effects[stan_data_tarsus$ind_id], THBW_egg$male_present, mean))

rbind(asreml=tapply(dam_sire_effects_asreml[stan_data_tarsus$dam_id], THBW_egg$male_present, mean),
stan=tapply(dam_sire_effects[stan_data_tarsus$dam_id], THBW_egg$male_present, mean))

nest_meansT <- aggregate(tarsus_mm~male_present+nest,THBW_egg,mean)
round(rbind(
asreml=with(nest_meansT,tapply( nest_effects_asreml,male_present,mean)),
stan=with(nest_meansT,tapply( nest_effects,male_present,mean))),3)

round(rbind(observed=with(nest_meansT,tapply( tarsus_mm,male_present,stand_skew)),
asreml=with(nest_meansT,tapply( nest_effects_asreml,male_present,stand_skew)),
stan=with(nest_meansT,tapply( nest_effects,male_present,stand_skew))),3)

round(rbind(observed=with(nest_meansT,tapply( tarsus_mm,male_present,var)),
asreml=with(nest_meansT,tapply( nest_effects_asreml,male_present,var)),
stan=with(nest_meansT,tapply( nest_effects,male_present,var))),3)

with(THBW_egg, tapply(tarsus_mm,hatch_day,stand_skew))

nest_meansW <- aggregate(wing_mm~nest+male_present,THBW_egg,mean)
with(nest_meansW,tapply( wing_mm,male_present,stand_skew))
with(nest_meansW,tapply( wing_mm,male_present,var))

nest_meansM <- aggregate(weight_g~nest+male_present,THBW_egg,mean)
with(nest_meansM,tapply( weight_g,male_present,stand_skew))
with(nest_meansM,tapply( weight_g,male_present,var))


plot(nest_effects_asreml~nest_meansT$tarsus_mm)
plot(nest_effects~nest_meansT$tarsus_mm)


par(mfrow=c(3,1))
for(i in unique(nest_meansT$male_present)) hist(nest_meansT[nest_meansT$male_present==i,"tarsus_mm"], breaks=75, xlim=c(14,18.5), ylim=c(0,40))
for(i in unique(nest_meansW$male_present)) hist(nest_meansW[nest_meansW$male_present==i,"wing_mm"], breaks=50, xlim=c(20,50))
for(i in unique(nest_meansM$male_present)) hist(nest_meansM[nest_meansM$male_present==i,"weight_g"], breaks=50, xlim=c(5,13))


table(nest_means$male_present)
table(THBW_egg$hatch_day)

# tapply(resid_stan, THBW_egg$male_present, mean)
# tapply(resid_asreml, THBW_egg$male_present, mean)




summary(mod_egg_DS)