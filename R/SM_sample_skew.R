rm(list=ls())

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

options(width=Sys.getenv("COLUMNS"), stringsAsFactors=FALSE)

source(paste0(wd,"R/functions.R"))

library(sn)


n<-11

(n*(n-1))/((n-2)^2)
(n^2-n)/(n^2-4*n+4)
(n^2-n)/((n-3)*(n-1))
1/(1-4/n)

1+4/n

skew_summary <- function(y){ c(mean=mean(y),
			sd=sqrt(sum((y-mean(y))^2)/length(y)), 
			sdC=sd(y),
			ssk=stand_skew(y),
			sskC1=stand_skew(y,sample=TRUE),
			sskC2=stand_skew(y,sample=TRUE,correction=2),
			sk=skew(y),
			skC=skew(y,sample=TRUE)#,jarrod_skew(y)
			)}


simSkew <- -0.99
n <- c(5,10,20,50,100,200,500,1000)
sims<-lapply(n, function(x){
	samp_V<- replicate(1000, {
		y <- rskn(x, mean=0, sd=1, skew=simSkew)
		skew_summary(y)
		})
	})



# ((n^2)/((n-2)*(n-1)))^1.5
# #((n)/(n-1))^1.5
# (n^2.5)/(((n-2)^2)*((n-1)^0.5))
# (sqrt(n*(n-1))/(n-2) )* ((n^2)/((n-2)*(n-1)))
# (sqrt(n*(n-1))/(n-2) )+ ((n^2)/((n-2)*(n-1)))
#(sqrt(n*(n-1))/(n-2) )

sim_means <- sapply(sims,apply,1,mean)
sim_sd <- sapply(sims,apply,1,sd)


sim_means

sim_means["skC",]/sim_means["sdC",]^3

  sapply(sims,function(x) cov(1/(x["sdC",]^3),x["skC",]))

  sapply(sims_exp,function(x) cor(1/(x[2,]^1.5),x[6,]))
lines(sim_exp_means[3,] - sapply(sims_exp,function(x) cov(1/(x[2,]^3),x[6,]))~n, col="purple")

{
	plot(sim_means["ssk",]~n, type="l", ylim=c(-1,0),log="x",ylab="simulation mean", xlab="log N")
	lines(sim_means["sskC1",]~n, col="red")
	lines(sim_means["sskC2",]~n, col="orange")
	lines(sim_means["sk",]~n, col="blue")
	lines(sim_means["skC",]~n, col="purple")
	legend("topright", c("stand skew","corrected stand skew 1","corrected stand skew 2","skew","corrected skew"),lwd=2, col=c(1,2,"orange",4,"purple"))

	abline(h=simSkew, col="grey")
}

 
 sapply(sims,function(x) cor(1/(x[2,]^1.5),x[7,]))
lines(sim_means[3,] - sapply(sims,function(x) cov(1/x[2,]^1.5,x[6,]))~n, col="purple")


sims_exp<-lapply(n, function(x){
	samp_V<- replicate(5000, {
		y <- rexp(x)
		skew_summary(y)
		})
	})

sim_exp_means <- sapply(sims_exp,apply,1,mean)
sim_exp_sd <- sapply(sims_exp,apply,1,sd)

{
	plot(sim_exp_means["ssk",]~n, type="l", ylim=c(-1,3),log="x",ylab="simulation mean", xlab="log N")
	lines(sim_exp_means["sskC1",]~n, col="red")
	lines(sim_exp_means["sskC2",]~n, col="orange")
	lines(sim_exp_means["sk",]~n, col="blue")
	lines(sim_exp_means["skC",]~n, col="purple")
	legend("topright", c("stand skew","corrected stand skew 1","corrected stand skew 2","skew","corrected skew"),lwd=2, col=c(1,2,"orange",4,"purple"))

	abline(h=2, col="grey")
}

sim_exp_means

sim_exp_means["skC",]/sim_exp_means["sdC",]^3

  sapply(sims_exp,function(x) cov(1/(x[2,]^3),x[6,]))
  sapply(sims_exp,function(x) var(1/(x[2,]^3)))
  sapply(sims_exp,function(x) var(x[6,]))

  sapply(sims_exp,function(x) cor(1/(x[2,]^1.5),x[6,]))
lines(sim_exp_means[3,] - sapply(sims_exp,function(x) cov(1/(x[2,]^3),x[6,]))~n, col="purple")

sims_exp

mean(1/(sims_exp[[1]][2,]^3)) * mean(sims_exp[[1]][6,]) + sapply(sims_exp,function(x) cov(1/(x[2,]^3),x[6,]))[1]

mean(sims_exp[[1]][6,]/(sims_exp[[1]][2,]^3) )

SE<-sqrt( ((n^5)/(((n-2)^4)*(n-1))) *6*(n-2)/((n+1)*(n+3)))
#SE<-sqrt( ((n^2)/((n-2)*(n-1)))^3 *6*(n-2)/((n+1)*(n+3)))
#((n^2)/((n-2)*(n-1)))^2 * 6*(n*(n-1)/((n-2)*(n+1)*(n+3)))

{
plot(sim_means[5,]~n, col="red", type="l", ylim=c(-2,1),log="x",ylab="simulation mean", xlab="log N")
lines(sim_means[5,]+sim_sd[5,]~n, col="red")
lines(sim_means[5,]-sim_sd[5,]~n, col="red")
lines(simSkew-SE~n, col="grey")
lines(simSkew+SE~n, col="grey")
	abline(h=simSkew, col="grey")
}

lines(sim_means[6,]+sim_sd[6,]~n, col="purple")
lines(sim_means[6,]-sim_sd[6,]~n, col="purple")

# sd skew covariance
# sapply(sims,function(x) cov(1/x[2,],x[6,]))
 sapply(sims,function(x) cov(1/x[2,]^1.5,x[6,]))
 sapply(sims,function(x) cor(1/(x[2,]^1.5),x[6,]))
  sapply(sims_exp,function(x) cov(1/(x[2,]^1.5),x[6,]))
  sapply(sims_exp,function(x) cor(1/(x[2,]^1.5),x[6,]))
# sim_means[3,] - sapply(sims,function(x) cov(1/x[2,],x[6,]))
# sim_means[4,] - sapply(sims,function(x) cov(1/x[2,],x[6,]))


sim_CI <- sapply(sims,apply,1,mean_CI,simplify = "array")

sim_means <- sapply(sims,apply,1,mean)
sim_sd <- sapply(sims,apply,1,sd)

par(mfrow=c(1,1))
plot(sim_CI[1,3,]~n, type="l", ylim=c(-2,1),log="x")
lines(sim_CI[2,3,]~n)
lines(sim_CI[3,3,]~n)

lines(sim_CI[1,4,]~n, col="red")
lines(sim_CI[2,4,]~n, col="red")
lines(sim_CI[3,4,]~n, col="red")

plot(sim_means[3,]~n, type="l", ylim=c(-2,1),log="x")
lines(sim_means[3,]+sim_sd[3,]~n)

plot(sim_means[3,]~n, type="l", ylim=c(-2,1),log="x")
lines(sim_means[3,]+sim_sd[3,]~n)
lines(sim_means[3,]-sim_sd[3,]~n)

lines(sim_means[4,]~n, col="red")
lines(sim_means[4,]+sim_sd[4,]~n, col="red")
lines(sim_means[4,]-sim_sd[4,]~n, col="red")


plot(sim_sd[4,]~SE)
abline(0,1)
points(sim_sd[3,]~SE)

sim_sd

SE<-sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3)))
lines(0.9-SE~n, col="blue")
lines(0.9+SE~n, col="blue")

# lines(sim_means[1,5,]~n, col="blue")
# lines(sim_means[2,5,]~n, col="blue")
# lines(sim_means[3,5,]~n, col="blue")

abline(h=-0.9, col="grey")