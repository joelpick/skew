
rm(list=ls())

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

source(paste0(wd,"R/00_functions.R"))

dryad_wd <- paste0(wd,"Data/Raw/skew_datasets/")


sum_func <- function(study_ref, data, metadata=NULL, sex=FALSE){
	formula <- if(sex){ tarsus~sex+age }else{ tarsus~age}
	index <- if(sex){ 3 }else{ 2 }
	out <- data.frame(study_ref=study_ref,
 		age = aggregate(formula,data,mean)[,index-1],
 		sex = if(sex){aggregate(formula,data,mean)[,1]}else{"U"},
		mean= aggregate(formula,data,mean, na.rm=TRUE)[,index],
		variance= aggregate(formula,data,var, na.rm=TRUE)[,index],
		skew= aggregate(formula,data,stand_skew,sample=TRUE)[,index],
		n=aggregate(formula,data,length)[,index])
	if(!is.null(metadata)){out <- cbind(out,species=metadata[metadata$study_ref==study_ref,"species"])}
	if(!sex) colnames(out)[2] <- "age"
	return(out)
}

refs <- read.csv(paste0(dryad_wd,"refs.csv"))



###------------------------------------------------------
# 1. Krause, E. Tobias; Krüger, Oliver; Schielzeth, Holger (2018), Data from: Long-term effects of early nutrition and environmental matching on developmental and personality traits in zebra finches, Dryad, Dataset, https://doi.org/10.5061/dryad.6j700
# https://doi.org/10.1016/j.anbehav.2017.04.003
###------------------------------------------------------

dd1 <- subset(read.table(paste0(dryad_wd,"01/Krause\ Krüger\ Schielzeth_Morph.txt"),header=TRUE),!is.na(Sex))
dd1$sex <- dd1$Sex

dd1a <- na.omit(rbind_notAnnoying(cbind(dd1[,c("sex","Tarsus_d17")],age="J"),cbind(dd1[,c("sex","Tarsus_d400")],age="A")))
dd1a$tarsus <- dd1a$Tarsus_d17

s01 <- sum_func("1",dd1a, metadata=refs, sex=FALSE)


###------------------------------------------------------
# 2. Pap, Peter Laszlo et al. (2019), Data from: Selection on multiple sexual signals in two Central- and Eastern-European populations of the barn swallow, v2, Dryad, Dataset, https://doi.org/10.5061/dryad.64p7k2f
# https://doi.org/10.1002/ece3.5629
###------------------------------------------------------
dd2 <- read.csv(paste0(dryad_wd,"02/BS_sexual_selection_CZ_RO_FULL_data.csv"))
dd2$age <- "A"
dd2$sex <- ifelse(dd2$sex==1,"F","M")
dd2$tarsus <- dd2$ltarsus

dd2a <- subset(dd2, country=="RO")
dd2b <- subset(dd2, country=="CZ")

s02 <- rbind(sum_func("2",dd2a, metadata=refs, sex=FALSE) ,sum_func("2",dd2b, metadata=refs, sex=FALSE) )


###------------------------------------------------------
# 3. Jacob, Staffan et al. (2015), Data from: Microbiome affects egg carotenoid investment, nestling development and adult oxidative costs of reproduction in Great tits, Dryad, Dataset, https://doi.org/10.5061/dryad.9n741
# https://doi.org/10.1111/1365-2435.12404
###------------------------------------------------------

dd3a <- read.table(paste0(dryad_wd,"03/Data_oxi.txt"),header=TRUE)
dd3a$age <-"A"
dd3a$sex <- dd3a$Sex

dd3b <- read.table(paste0(dryad_wd,"03/Data_nestlings.txt"),header=TRUE)
dd3b$age <-"J"
dd3b$sex <-"U"


dd3 <- na.omit(rbind_notAnnoying(dd3a[,c("sex","age","Tarsus_length")],dd3b[,c("sex","age","Tarsus_length_j14")]) )
dd3$tarsus <- dd3$Tarsus_length

s03 <- sum_func("3",dd3, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 4. Forstmeier, Wolfgang et al. (2017), Data from: Testing the phenotype-linked fertility hypothesis in the presence and absence of inbreeding, Dryad, Dataset, https://doi.org/10.5061/dryad.515c6
# https://doi.org/10.1111/jeb.13062
###------------------------------------------------------

# better sample in 13b
# dd4 <- read.csv(paste0(dryad_wd,"04/Dryad_data_Forstmeieretal2017_SpermPhenoCorrel.csv"))
# dd4$age <- "A"
# dd4$sex <- "M"
# dd4$tarsus <- dd4$Tarsus

# s04 <- sum_func("04",dd4, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 5. Cornell, Allison; Gibson, Kate F.; Williams, Tony D. (2017), Data from: Physiological maturity at a critical life-history transition and flight ability at fledging, Dryad, Dataset, https://doi.org/10.5061/dryad.c2n66
# https://doi.org/10.1111/1365-2435.12777
###------------------------------------------------------

dd5 <- read.csv(paste0(dryad_wd,"05/chickdatadryad.csv"))
dd5$age <- "J"

s05 <- sum_func("5",dd5, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 6. Vermeulen, Anke; Müller, Wendt; Eens, Marcel (2016), Data from: Vitally important – does early innate immunity predict recruitment and adult innate immunity?, Dryad, Dataset, https://doi.org/10.5061/dryad.p0s3g
# https://doi.org/10.1002/ece3.1939
###------------------------------------------------------

dd6 <- read.csv(paste0(dryad_wd,"06/Vermeulen_et_al_DryadData.csv"))
dd6$age <- "J"
dd6$sex <- "U"
dd6$tarsus <- dd6$Tarsus
s06 <- sum_func("6",dd6, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 7. Song, Zitan et al. (2018), Data from: Silver spoon effects of hatching order in an asynchronous hatching bird, Dryad, Dataset, https://doi.org/10.5061/dryad.184c1dj
# https://doi.org/10.1093/beheco/ary191
###------------------------------------------------------

dd7 <- read.csv(paste0(dryad_wd,"07/extracted.csv"))
dd7$age <- "J"
dd7$sex <- ifelse(dd7$sex==0,"F","M")

s07 <- sum_func("7",dd7, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 8. Simpson, Richard K.; McGraw, Kevin J. (2017), Data from: Multiple signaling in a variable environment: expression of song and color traits as a function of ambient sound and light, Dryad, Dataset, https://doi.org/10.5061/dryad.1j81k
# https://doi.org/10.1111/btp.12528
###------------------------------------------------------

dd8 <- read.csv(paste0(dryad_wd,"08/RTAT_data.csv"))
dd8$age <- "A"
dd8$sex <- "M"
dd8$tarsus <- dd8$Tarsus
s08 <- sum_func("8",dd8, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 9. Torres, Roxana; Chin, Eunice; Rampton, Rowan; Williams, Tony D. (2019), Data from: Are there synergistic or antagonistic effects of multiple maternally-derived egg components (antibodies and testosterone) on offspring phenotype?, Dryad, Dataset, https://doi.org/10.5061/dryad.j348s75
# https://doi.org/10.1242/jeb.196956
###------------------------------------------------------

dd9a <- subset(read.table(paste0(dryad_wd,"09/Chick_growth.prn"),header=TRUE),age==15)
dd9a$age <- "J"
dd9a$tarsus <- as.numeric(dd9a$tarsus)

dd9b <- read.table(paste0(dryad_wd,"09/Mom_repro.prn"),header=TRUE)
dd9b$age <- "A"
dd9b$sex <- "F"
dd9b$tarsus <- as.numeric(dd9b$ftars)

dd9<-rbind(dd9a[,c("sex","age","tarsus")],dd9b[,c("sex","age","tarsus")])

s09 <- sum_func("9",dd9, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 10. Camacho, Carlos; Canal, David; Potti, Jaime (2014), Data from: Nonrandom dispersal drives phenotypic divergence within a bird population, Dryad, Dataset, https://doi.org/10.5061/dryad.h22n9
# https://doi.org/10.1002/ece3.563
###------------------------------------------------------

dd10 <- read.csv(paste0(dryad_wd,"10/INDIV_DATA.csv"))
dd10$age <- "A"
dd10$sex <- dd10$SEX
dd10$tarsus <- dd10$T.LENGTH

s10 <- sum_func("10",dd10, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 11. Poissant, Jocelyn et al. (2016), Data from: Multivariate selection and intersexual genetic constraints in a wild bird population, Dryad, Dataset, https://doi.org/10.5061/dryad.qt745
# https://doi.org/10.1111/jeb.12925
###------------------------------------------------------

dd11 <- read.csv(paste0(dryad_wd,"11/selection_analysis_data.csv"))

dd11$age <- "A"
dd11$tarsus <- dd11$tars

s11 <- sum_func("11",dd11, metadata=refs, sex=FALSE) 





###------------------------------------------------------
# 12. Dubuc-Messier, Gabrielle et al. (2018), Data from: Gene flow does not prevent personality and morphological differentiation between two blue tit populations, Dryad, Dataset, https://doi.org/10.5061/dryad.31tc3s8
# https://doi.org/10.1111/jeb.13291
###------------------------------------------------------

dd12 <- read.csv(paste0(dryad_wd,"12/Tarsus_Dryad.csv"))
dd12$age <- "A"
dd12$sex <- ifelse(dd12$Sex=="female","F","M")
dd12$tarsus <- dd12$Tarsus

dd12a <- subset(dd12, Population=="Evergreen")
dd12b <- subset(dd12, Population=="Deciduous")

s12 <- rbind(sum_func("12",dd12a, metadata=refs, sex=FALSE) ,sum_func("12",dd12b, metadata=refs, sex=FALSE) )



###------------------------------------------------------
# 13. Husby, Arild et al. (2012), Data from: Sex chromosome linked genetic variance and the evolution of sexual dimorphism of quantitative traits, Dryad, Dataset, https://doi.org/10.5061/dryad.451n7
# https://doi.org/10.1111/j.1558-5646.2012.01806.x
###------------------------------------------------------

dd13a <- read.csv(paste0(dryad_wd,"13/data_file_dryad_collared_flycatcher.csv"), na.strings="*")
dd13a$age <- "A"
dd13a$tarsus <- dd13a$tars

dd13aa<-aggregate(tarsus~age+sex+ring,dd13a,mean)

s13a <- sum_func("13a",dd13aa, metadata=refs, sex=FALSE) 



dd13b <- read.csv(paste0(dryad_wd,"13/ZebraFinch_TarsusLength.csv"))
dd13b$sex <- ifelse(dd13b$Sex==0,"F","M")
dd13b$age <- "A"
dd13b$tarsus <- dd13b$Tarsus

s13b <- sum_func("13b",dd13b, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# Nishida, Yuusuke; Takagi, Masaoki (2018), Data from: Song performance is a condition-dependent dynamic trait honestly indicating the quality of paternal care in the Bull-headed Shrike, Dryad, Dataset, https://doi.org/10.5061/dryad.c84f7c4
#  https://doi.org/10.1111/jav.01794

###------------------------------------------------------

dd14 <- read.csv(paste0(dryad_wd,"14/extracted.csv"))

s14 <- sum_func("14",dd14, metadata=refs, sex=FALSE) 




###------------------------------------------------------
# 15. Moiron, Maria et al. (2018), Data from: Functional relations between body mass and risk-taking behavior in wild great tits, Dryad, Dataset, https://doi.org/10.5061/dryad.14cn58v
# https://doi.org/10.1093/beheco/ary199
###------------------------------------------------------

dd15 <- read.table(paste0(dryad_wd,"15/Moiron\ et\ al_Dataset\ Integration\ MS.txt"),header=TRUE)

dd15$age <- "A"
dd15$sex <- "M"
dd15$tarsus <- dd15$Tarsus

s15 <- sum_func("15",aggregate(tarsus~age+sex+IndividualIdentity,dd15,mean), metadata=refs, sex=FALSE) 




###------------------------------------------------------
# 16. Krist, Miloš et al. (2010), Data from: Egg size and offspring performance in the collared flycatcher (Ficedula albicollis): a within-clutch approach, Dryad, Dataset, https://doi.org/10.5061/dryad.1758
# https://doi.org/10.1007/s00442-004-1568-5

###------------------------------------------------------

dd16 <- read.csv(paste0(dryad_wd,"16/dataset.csv"))

dd16$tarsus <- dd16$tarsus.length
dd16$sex <- ifelse(dd16$sex=="female","F","M")
dd16$age <- "J"

s16 <- sum_func("16",dd16, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 17. Krist, Miloš; Janča, Martin; Edme, Anaïs; Dzuro, Rudolf (2016), Data from: Are prenatal maternal resources more important in competitive than in benign postnatal environments?, Dryad, Dataset, https://doi.org/10.5061/dryad.823f0
# https://doi.org/10.1642/AUK-14-236.1
###------------------------------------------------------

dd17 <- read.table(paste0(dryad_wd,"17/dataset.csv"),sep=";",header=TRUE)
dd17$tarsus <- dd17$Fledging_tarsus
dd17$sex <- "U"
dd17$age <- "J"

s17 <- sum_func("17",dd17, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 18. Kvalnes, Thomas et al. (2017), Data from: Reversal of response to artificial selection on body size in a wild passerine, Dryad, Dataset, https://doi.org/10.5061/dryad.v50r8
# https://doi.org/10.1111/evo.13277
###------------------------------------------------------

dd18 <- read.table(paste0(dryad_wd,"18/sparrow_experiment_2002-2012_data_2017-05-12.txt"),sep=";",header=TRUE)

dd18$tarsus <- dd18$tars
dd18$sex <- ifelse(dd18$sex=="f","F","M")
dd18$age <- "A"

s18 <- sum_func("18",aggregate(tarsus~age+sex+ringnr,dd18,mean), metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 19. Nord, Andreas; Nilsson, Jan-Åke (2011), Data from: Incubation temperature affects growth and energy metabolism in blue tit nestlings, Dryad, Dataset, https://doi.org/10.5061/dryad.jb314
# https://doi.org/10.1086/662172
###------------------------------------------------------

dd19 <- read.table(paste0(dryad_wd,"19/nestling\ biometric\ data.txt"),header=TRUE)

dd19$tarsus <- dd19$tars14
dd19$sex <- "U"
dd19$age <- "J"

s19 <- sum_func("19", dd19, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 20. Bebbington, Kat et al. (2016), Data from: Consequences of sibling rivalry vary across life in a passerine bird, Dryad, Dataset, https://doi.org/10.5061/dryad.12np0
# https://doi.org/10.1093/beheco/arw167
###------------------------------------------------------

dd20 <- read.csv(paste0(dryad_wd,"20/chicks.csv"))

dd20$tarsus <- dd20$TarsusChick
dd20$sex <- ifelse(dd20$Sex=="Female","F","M")
dd20$age <- "J"

s20 <- sum_func("20", dd20, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 21. Arct, Aneta et al. (2019), Data from: Parental genetic similarity and offspring performance in blue tits in relation to brood size manipulation, v2, Dryad, Dataset, https://doi.org/10.5061/dryad.v6r0758
# https://doi.org/10.1002/ece3.5367
###------------------------------------------------------

dd21 <- read.csv(paste0(dryad_wd,"21/gs_arct_ecolevol2019.csv"))

dd21$tarsus <- dd21$RTARS14
dd21$sex <- dd21$SEX
dd21$age <- "J"

s21 <- sum_func("21", dd21, metadata=refs, sex=FALSE) 





###------------------------------------------------------
# 22. Kvalnes, Thomas et al. (2018), Data from: Offspring fitness and the optimal propagule size in a fluctuating environment, Dryad, Dataset, https://doi.org/10.5061/dryad.m74c7m9
# https://doi.org/10.1111/jav.01786
###------------------------------------------------------

dd22 <- read.table(paste0(dryad_wd,"22/house_sparrow_ind_data.txt"),sep=";",header=TRUE)

dd22$tarsus <- dd22$tars11
dd22$sex <- "U"
dd22$age <- "J"

s22 <- sum_func("22", dd22, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 23. DuRant, Sarah E. et al. (2014), Data from: Are thyroid hormones mediators of incubation temperature-induced phenotypes in birds?, Dryad, Dataset, https://doi.org/10.5061/dryad.hb87k
# https://doi.org/10.1098/rsbl.2013.0950
###------------------------------------------------------

#dd23 <- read.table(paste0(dryad_wd,"23/DuRant\ et\ al.\ Dataset\ 2.txt"),header=TRUE)

### data a real mess



###------------------------------------------------------
# 24. Grunst, Melissa L. et al. (2019), Data from: Artificial light at night does not affect telomere shortening in a developing free-living songbird: a field experiment, Dryad, Dataset, https://doi.org/10.5061/dryad.8216g63
# https://doi.org/10.1016/j.scitotenv.2018.12.469
###------------------------------------------------------

dd24 <- subset(read.csv(paste0(dryad_wd,"24/ALAN\ and\ telomere\ shortening_DRYAD_2118.csv")), Sample.pt ==15 & !is.na(Mean.Cp.T) & Sex!="")

dd24$tarsus <- dd24$Tarsus
dd24$sex <- dd24$Sex
dd24$age <- "J"

s24 <- sum_func("24", dd24, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 25. Silva, Catarina N. S. et al. (2017), Data from: Insights into the genetic architecture of morphological traits in two passerine bird species, Dryad, Dataset, https://doi.org/10.5061/dryad.786m4
# https://doi.org/10.1038/hdy.2017.29
###------------------------------------------------------

# better sample in 18
# dd25a <- read.table(paste0(dryad_wd,"25/pheno_rep_Sparrows.txt"),header=TRUE)

# dd25a$sex <- ifelse(dd25a$sex==0,"M","F")
# dd25a$age <- "A"

# s25a <- sum_func("25a",aggregate(tarsus~sex+age+id,dd25a,mean), metadata=refs, sex=FALSE) 


# better sample in 40
# dd25b <- subset(read.table(paste0(dryad_wd,"25/meanphenos_flycatchers.txt"),header=TRUE), tarsus>16)
# ##excluded one individual with very extreme and unlikely phenotype

# dd25b$sex <- ifelse(dd25b$sex_cat=="m","M","F")
# dd25b$age <- "A"

# s25b <- sum_func("25b", dd25b, metadata=refs, sex=FALSE) 





###------------------------------------------------------
# 26. Perrier, Charles et al. (2018), Data from: Heritability estimates from genome wide relatedness matrices in wild populations: application to a passerine, using a small sample size, Dryad, Dataset, https://doi.org/10.5061/dryad.k6r1mk8
# https://doi.org/10.1111/1755-0998.12886
###------------------------------------------------------

dd26 <- read.table(paste0(dryad_wd,"26/traits.txt"),header=TRUE)

dd26$tarsus <- as.numeric(sub(",",".",dd26$tarsed))
dd26$sex <- "U"
dd26$age <- "A"

s26 <- sum_func("26", dd26, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 27. Björklund, Mats; Gustafsson, Lars (2017), Data from: Subtle but ubiquitous selection on body size in a natural population of collared flycatchers over 33 years, Dryad, Dataset, https://doi.org/10.5061/dryad.24tm0
#  https://doi.org/10.1111/jeb.13117
###------------------------------------------------------

# repeated from 13a
# dd27f <- read.csv(paste0(dryad_wd,"27/FemaleSurvival_Dryad_Tot.csv"))[,1:4]
# dd27f$sex <- "F"

# dd27m <- read.csv(paste0(dryad_wd,"27/MaleSurv_DryadTot.csv"))[,1:4]
# dd27m$sex <- "M"

# #table(rbind(dd27m,dd27f)$Age,rbind(dd27m,dd27f)$sex)
# # both sexes have most observations at age 1, so will just include those to avoid repeated measurements

# dd27 <- subset(rbind(dd27m,dd27f), Age==1)
# dd27$age <- "A"
# dd27$tarsus <- dd27$Tarsus
# s27 <- sum_func("27", dd27, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 28. Becker, Philipp J. J. et al. (2015), Data from: Mother-offspring and nest mate resemblance but no heritability in early-life telomere length in white-throated dippers, Dryad, Dataset, https://doi.org/10.5061/dryad.b2v37
#  https://doi.
###------------------------------------------------------

dd28 <- read.table(paste0(dryad_wd,"28/BeckerReichert-et-al_data_telomeres.txt"), skip=21, header=TRUE)

dd28$age <- "J"
dd28$tarsus <- dd28$TARSUS
dd28$sex <- dd28$SEX

s28 <- sum_func("28", dd28, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 29. Weber, Beth M. et al. (2019), Data from: Pre- and post-natal effects of experimentally manipulated maternal corticosterone on growth, stress reactivity, and survival of nestling house wrens, Dryad, Dataset, https://doi.org/10.5061/dryad.16049f4
# https://doi.org/10.1111/1365-2435.13126
###------------------------------------------------------

dd29a <- read.csv(paste0(dryad_wd,"29/Dryad\ data\ Weber\ et\ al1.csv"), na.strings=".")
dd29b <- read.csv(paste0(dryad_wd,"29/Dryad\ data\ Weber\ et\ al2.csv"), na.strings=".")

dd29 <- rbind(
	na.omit(data.frame(age="J",sex="U",tarsus=dd29a$TARSUS11)),
	na.omit(data.frame(age="J",sex="U",tarsus=dd29b$BD11Tarsus)),
	aggregate(tarsus~id+age+sex,data.frame(id=c(dd29a$BIOL_MOTHER, dd29a$FBAND),tarsus=c(dd29a$BIOL_MOTHER_TARSUS, dd29a$FEM_TARSUS),age="A",sex="F"),mean)[,2:4]
	)

s29 <- sum_func("29", dd29, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 30. Poorboy, Dylan et al. (2018), Data from: Experimental cross-fostering of eggs reveals effects of territory quality on reproductive allocation, Dryad, Dataset, https://doi.org/10.5061/dryad.h8v8157
# https://doi.org/10.1093/beheco/ary098
###------------------------------------------------------

dd30m <- read.csv(paste0(dryad_wd,"30/MALES.csv"))
dd30m$sex <- "M"

dd30f <- read.csv(paste0(dryad_wd,"30/FEMALES.csv"))
dd30f$sex <- "F"

dd30a <- aggregate(TARSUS~BAND.+sex,rbind(dd30m,dd30f),mean)
dd30a$age <- "A"
dd30a$tarsus <- dd30a$TARSUS

dd30c1 <- read.csv(paste0(dryad_wd,"30/CHICKS_PRE.csv"))
dd30c2 <- read.csv(paste0(dryad_wd,"30/CHICKS_POST.csv"))


dd30 <- rbind(dd30a[,c(2,4,5)],data.frame(sex="U",age="J",tarsus=c(dd30c1$Tarsus,dd30c2$Tarsus))
	)

s30 <- sum_func("30", dd30, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 31. Sakaluk, Scott K. et al. (2014), Data from: Genetic and environmental variation in condition, cutaneous immunity, and haematocrit in house wrens, Dryad, Dataset, https://doi.org/10.5061/dryad.jk2m0
# https://doi.org/10.1186/s12862-014-0242-8
###------------------------------------------------------

dd31a <- subset(read.csv(paste0(dryad_wd,"31/Sakaluk\ house\ wren\ \ pedigree\ and\ trait\ values.csv"), na.strings="."),TARSUS>14)

##excluded one individual with very extreme and unlikely phenotype
# subset(dd31a,TARSUS<14)
# plot(dd31a$MASS,dd31a$TARSUS)

dd31 <- na.omit(data.frame(
	tarsus = dd31a$TARSUS, 
	sex = ifelse(dd31a$SEX=="MALE","M","F"),
		#ifelse(dd31a$SEX=="FEMALE","F",NA)), 
	age = "J"))

s31 <- sum_func("31", dd31, metadata=refs, sex=FALSE) 




### check ifelse statements for sex to make sure there aren't NAs that are assigned to a sex



###------------------------------------------------------
# 32. Berzins, Lisha L.; Gilchrist, H. Grant; Burness, Gary (2015), Data from: No assortative mating based on size in black guillemots breeding in the Canadian Arctic, Dryad, Dataset, https://doi.org/10.5061/dryad.1bm5t
# https://doi.org/10.1675/063.032.0313
###------------------------------------------------------

dd32 <- read.csv(paste0(dryad_wd,"32/Black\ guillemot\ discriminant\ function\ analysis.csv"))

head(dd32)

dd32$sex <- dd32$SEX
dd32$tarsus <- dd32$TARSUS
dd32$age <- "A"

s32 <- sum_func("32", dd32, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 33. Bebbington, Kat et al. (2016), Data from: Telomere length reveals cumulative individual and transgenerational inbreeding effects in a passerine bird, v2, Dryad, Dataset, https://doi.org/10.5061/dryad.52fp4
# https://doi.org/10.1111/mec.13670
###------------------------------------------------------

dd33 <- read.csv(paste0(dryad_wd,"/33/Adults\ Hz.csv"))

dd33$age <- "A"
dd33$sex <- ifelse(dd33$Sex=="Females","F","M")
dd33$tarsus <- dd33$Tarsus

s33 <- sum_func("33", aggregate(tarsus~BirdID+age+sex,dd33,mean), metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 34. Caizergues, Aude E.; Gregoire, Arnaud; Charmantier, Anne (2018), Data from: Urban versus forest ecotypes are not explained by divergent reproductive selection, Dryad, Dataset, https://doi.org/10.5061/dryad.tv45802
# https://doi.org/10.1098/rspb.2018.0261
###------------------------------------------------------


dd34 <- read.csv(paste0(dryad_wd,"34/morpho_data_open.csv"))

head(dd34)

dd34$age <- "A"
dd34$sex <- ifelse(dd34$sex==2,"F","M")
dd34$tarsus <- dd34$TD

dd34a <- subset(dd34, habitat=="urb")
dd34b <- subset(dd34, habitat=="rur")

s34 <- rbind(sum_func("34",dd34a, metadata=refs, sex=FALSE) ,sum_func("34",dd34b, metadata=refs, sex=FALSE) )



###------------------------------------------------------
# 35. Podofillini, Stefano et al. (2020), Data from: Benefits of extra food to reproduction depend on maternal condition, Dryad, Dataset, https://doi.org/10.5061/dryad.5db0168
# https://doi.org/10.1111/oik.06067
###------------------------------------------------------

dd35 <- read.csv(paste0(dryad_wd,"35/nestling_data.csv"))

dd35$age <- "J"
dd35$sex <- ifelse(dd35$SE==0,"F","M")
dd35$tarsus <- dd35$TA

s35 <- sum_func("35",dd35, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 36. Rioux Paquette, Sébastien et al. (2014), Data from: Severe recent decrease of adult body mass in a declining insectivorous bird population, Dryad, Dataset, https://doi.org/10.5061/dryad.67t23
# https://doi.org/10.1098/rspb.2014.0649
###------------------------------------------------------

dd36f <- read.csv(paste0(dryad_wd,"36/femaledata.csv"),skip=14)
dd36f$sex<-"F"
dd36m <- read.csv(paste0(dryad_wd,"36/maledata.csv"),skip=13)
dd36m$sex<-"M"

dd36 <- rbind(aggregate(tarsus~femaleID+sex,dd36f,mean)[,2:3],aggregate(tarsus~maleID+sex,dd36m,mean)[,2:3])
dd36$age <- "A"

s36 <- sum_func("36",dd36, metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 37. Cox, Amelia R. et al. (2019), Data from: Rainy springs linked to poor nestling growth in a declining avian aerial insectivore (Tachycineta bicolor), Dryad, Dataset, https://doi.org/10.5061/dryad.7m41jd8
# https://doi.org/10.1098/rspb.2019.0018
###------------------------------------------------------

dd37a <- subset(read.csv(paste0(dryad_wd,"37/Long\ term\ adult\ body\ size\ data.csv")),!is.na(tarsus))
dd37a$age <- "A"
## clearly outlier tarsi in one year - decimal in incorrect place
dd37a[dd37a$tarsus<5,"tarsus"] <- dd37a[dd37a$tarsus<5,"tarsus"]*10

## something weird with this data!!
# hist(subset(dd37a,sex=="M")$tarsus,breaks=100)
# hist(subset(dd37a,sex=="F")$tarsus,breaks=100)


## very strange things going on in long term chick data, so just use the detailed data from 2017
dd37j <- subset(read.csv(paste0(dryad_wd,"37/Nestling\ growth\ data\ 2017.csv")), Age==12)

dd37j$age <- "J"
dd37j$sex <- "U"
dd37j$tarsus <- dd37j$Tarsus

dd37 <- rbind(aggregate(tarsus~band+sex+age,dd37a,mean)[,2:4],dd37j[,c("sex","age","tarsus")])

s37 <- sum_func("37",dd37, metadata=refs, sex=FALSE) 




###------------------------------------------------------
# 38. DeSimone, Joely G.; Clotfelter, Ethan D.; Black, Elizabeth C.; Knutie, Sarah A. (2017), Data from: Avoidance, tolerance, and resistance to ectoparasites in nestling and adult tree swallows, Dryad, Dataset, https://doi.org/10.5061/dryad.9bb60
# https://doi.org/10.1111/jav.01641
###------------------------------------------------------
dd38a <- read.csv(paste0(dryad_wd,"38/DeSimone\ et\ al\ dataset\ 2014-2016\ grouped\ by\ nest.csv"))

dd38j <- read.csv(paste0(dryad_wd,"38/DeSimone\ et\ al\ dataset\ 2014-2016\ grouped\ by\ nestling.csv"))

dd38 <- na.omit(rbind(
	data.frame(age="A",sex="F",tarsus=rowMeans(dd38a[,c("Mom_firstcap_tarsus","Mom_secondcap_tarsus")],na.rm=TRUE)),
	data.frame(age="J",sex="U",tarsus=dd38j$d13_tarsus)
	))

s38 <- sum_func("38",dd38, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 39. Class, Barbara; Brommer, Jon (2020), Data from: Can dominance genetic variance be ignored in evolutionary quantitative genetic analyses of wild populations?, v3, Dryad, Dataset, https://doi.org/10.5061/dryad.zpc866t6d
#https://doi.org/10.5061/dryad.zpc866t6d
###------------------------------------------------------

dd39 <- read.csv(paste0(dryad_wd,"39/Young_pheno.csv"))

head(dd39)
dd39$tarsus <- dd39$Tarsus
dd39$sex <- ifelse(dd39$Sex=="m","M",ifelse(dd39$Sex=="f","F","U"))
dd39$age <- "J"
s39 <- sum_func("39",dd39, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 40. Eryns data
###------------------------------------------------------

dd40 <- read.csv(paste0(dryad_wd,"40/eryn_skew.csv"))

s40a <- cbind(study_ref="40a",dd40[grep("pied",dd40$X),2:7], species="Ficedula hypoleuca")
s40b <- cbind(study_ref="40b",dd40[grep("collared",dd40$X),2:7],species="Ficedula albicollis")



###------------------------------------------------------
# 41. Sparrows - Ihle et al. 2019
###------------------------------------------------------

dd41 <- read.csv(paste0(dryad_wd,"41/R_MY_TABLE_perChick.csv"))

dd41$tarsus <- dd41$AvgOfTarsus
dd41$age <- "J"
dd41$sex <- "U"

s41 <- sum_func("41",dd41, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 42. Dunnocks - from Shinichi
###------------------------------------------------------

dd42j <- subset(read.csv(paste0(dryad_wd,"42/chicks.csv")),!is.na(tarsus)&Age=="8")

# day 8 gives most tarsus measures
# table(dd42j$Age) 

dd42j$sex <- ifelse(dd42j$Sex=="Male","M",ifelse(dd42j$Sex=="Female","F","U"))
dd42j$age <- "J"

dd42a <- subset(read.csv(paste0(dryad_wd,"42/adults.csv")), grepl("\\+",Age) & tarsus>0)
dd42a$sex <- ifelse(dd42a$Sex=="Male","M",ifelse(dd42a$Sex=="Female","F","U"))
dd42a$age <- "A"

dd42 <- rbind(dd42j[,c("sex","age","tarsus")],aggregate(tarsus~bird_id+sex+age,dd42a,mean)[2:4])

s42 <- sum_func("42",dd42, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 43. Our data
###------------------------------------------------------

dd43 <- subset(read.csv(paste0(wd,"Data/Raw/tMORPH.csv")), substr(date,1,4) %in% 2010:2018)
dd43s <- read.csv(paste0(wd,"Data/Raw/tBIRDS.csv"))[,c("bird_id","sex")]
dd43$tarsus <- dd43$tarsus_mm
dd43$sex <- dd43$morph_sex

dd43j <- merge(aggregate(tarsus~bird_id,subset(dd43, is.na(agecode) & !is.na(tarsus)),mean),dd43s)
dd43j$age <- "J"

dd43a <- aggregate(tarsus~bird_id+sex,subset(dd43, !is.na(agecode) & !is.na(tarsus)),mean)
dd43a$age <- "A"


s43 <- sum_func("43",rbind(dd43a[,2:4],dd43j[,2:4]), metadata=refs, sex=FALSE) 



###------------------------------------------------------
# 44. 
###------------------------------------------------------

dd44 <- read.csv(paste0(dryad_wd,"44/Moller.csv"), na.string=".")
dd44$tarsus <- dd44$Tarsus..mm.x.100/100
dd44$sex <- ifelse(dd44$Sex=="Male","M",ifelse(dd44$Sex=="Female","F","U"))
dd44$age <- "A" ##ifelse(dd44$Age..years.>1,"A","J")

s44 <- sum_func("44",dd44, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 45. 
###------------------------------------------------------

dd45 <- read.csv(paste0(dryad_wd,"45/Santiago.csv"))
str(dd45)
dd45$tarsus <- dd45$Tarsus
dd45$sex <- ifelse(dd45$Sex=="male","M",ifelse(dd45$Sex=="female","F",NA))
dd45$age <- ifelse(dd45$Age=="adult","A",ifelse(dd45$Age=="juvenile","J",NA)) 

s45 <- sum_func("45",dd45, metadata=refs, sex=FALSE) 


###------------------------------------------------------
# 46. 
###------------------------------------------------------

dd46 <- read.csv(paste0(dryad_wd,"46/jovani.csv"), na.string="")

dd46$sex <- ifelse(is.na(dd46$SEX),"U",dd46$SEX)
dd46$age <- ifelse(dd46$AGE%in%c("3","3J"),"J",ifelse(dd46$AGE==2,NA,"A"))
dd46$tarsus <- dd46$TARS
dd46$SPECIES <- ifelse(dd46$SPECIES=="Parus caeruleus","Cyanistes caeruleus",dd46$SPECIES)
Cyanistes
## clearly outlier tarsi in one species
dd46 <- subset(dd46, !(SPECIES=="Erithacus rubecula" & tarsus<15))

str(dd46)

dd46_split <- split(dd46,dd46$SPECIES)
s46 <- lapply_dataFrame(1:length(dd46_split), function(x){
	cbind(sum_func(paste0("46",letters[x]),dd46_split[[x]]),species=dd46_split[[x]]$SPECIES[1])
})
s46 <- subset(s46, n>2)


###------------------------------------------------------
# Others. 
###------------------------------------------------------




###------------------------------------------------------
# Summary
###------------------------------------------------------

all <- do.call(rbind,lapply(ls()[grep("s\\d",ls())],get))

write.csv(all, file=paste0(wd,"Data/Intermediate/MA_data.csv"), row.names=FALSE)


