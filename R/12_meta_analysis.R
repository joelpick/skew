rm(list=ls())

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}
dryad_wd <- paste0(wd,"Data/Raw/skew_datasets/")

all<- read.csv(paste0(wd,"Data/Intermediate/MA_data.csv"))


length(unique(all$species))

table(all$sex,all$age)
boxplot(skew~age,all)
boxplot(skew~age+sex,all, col=c(2,5)); abline(h=0)
#points(c(1,3,2,4),s43$skew,pch=19, cex=2)
with(all,tapply(skew,list(age),mean))
with(all,tapply(skew,list(age,sex),mean))

library(scales)
all$group <- as.numeric(as.factor(paste(all$sex,all$age)))
plot(skew~group,all,cex=0)
abline(h=0)
symbols(all$group,all$skew, circle= sqrt( (all$n)/pi )/200, inches=FALSE, bg=alpha(c(2,1),0.5)[as.factor(all$age)], add=TRUE)
symbols(c(1,3,2,4),s43$skew,, circle= sqrt( (s43$n)/pi )/200, inches=FALSE,bg=3, add=TRUE)


# all$SE <-sqrt(6*all$n*(all$n-1)/((all$n-2)*(all$n+1)*(all$n+3)))
all$SE <- with(all,sqrt( ((n^2)/((n-2)*(n-1)))^3 *6*(n-2)/((n+1)*(n+3))))

plot(skew~sqrt(n),all, col=as.factor(all$age), pch=19)
plot(skew~n,all, col=as.factor(all$age), pch=19, log="x")
with(all,arrows(n,skew+SE,n,skew-SE,length=0.01, angle=90, code=3))
abline(h=0, col="grey")

# Random effects meta-analysis in MCMCglmm.
# Population and species were modelled as random effects.
# Sex and age and their interaction were included as fixed effects.

## sampling error

library(MCMCglmm)

prior <- list(R = list(V = diag(1), nu = 0.002),
               G = list(
               	 G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = 1000),
                 G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = 1000),
                 G3 = list(V = diag(1), fix = 1)))
a <- 5
MAmod<-MCMCglmm(skew~age,#all$
	random = ~ species + study_ref + idh(SE):units,
	data=all,
	prior = prior,
	nitt=13000*a,thin=10*a,burnin=3000*a)
summary(MAmod)
plot(MAmod)

#aod::wald.test(cov(MAmod$Sol[,2:3]), colMeans(MAmod$Sol[,2:3]), Terms=1:2)$result$chi2["P"]
## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2017q3/025933.html

resid <- all$skew - predict(MAmod) 
plot(sqrt(all$n)~resid, pch=19, col=as.factor(all$age))
abline(v=0, col="grey")

n<-1:10000
SE<-sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3)))
lines(sqrt(n)~I(-1*SE*1.96), col="blue")
lines(sqrt(n)~I(SE*1.96), col="blue")

save(MAmod, file= paste0(wd,"Data/Intermediate/meta_analysis.Rdata"))

setEPS()
pdf(paste0(wd,"R/plots/figure_meta_analysis.pdf"), height=6, width=6)

par(mfrow=c(1,1),mar=c(5,5,0,0),cex.lab=1.5)
boxplot(skew~age,all, ylim=c(-5,1.7), ylab="Skew", xlab="Age",names=c("Adult","Juvenile"));abline(h=0, col="darkgrey")
text(c(1,2),c(1.5,1.5),paste0(table(all$age)," (",tapply(all$species,all$age,function(x)length(unique(x))),")"))
points(c(1,2),all$skew[all$study_ref==43], pch=19,col="blue", cex=1.2)


dev.off()




rbind(summary(MAmod)$Gcov,summary(MAmod)$Rcov)
summary(MAmod)$sol

head(MAmod$Sol)

