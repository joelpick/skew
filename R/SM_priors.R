

rm(list=ls())

options( stringsAsFactors=FALSE)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "~/Work/Skew/"
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/00_functions.R"))

get_gamma<-function(alpha=NULL, delta=NULL, nu){
	if(is.null(delta)) delta = alpha / sqrt(1 + alpha^2); 
	b_nu = sqrt(nu/pi) * gamma((nu-1)/2)/gamma(nu/2);
	sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
	skew = (b_nu*delta)/sigma_z^3 * ( (nu*(3-delta^2))/(nu-3) - (3*nu)/(nu-2)  + 2*(b_nu*delta)^2 );
	return(skew);
}


density_plot <- function(x,col="grey", xlim=range(x)){
	d <- density(x, adjust=2.5 )
	plot(d$x,d$y, type="l", bty="n", xlim=xlim, xlab="", ylab="")#, xaxt="n", yaxt="n"
	polygon(d$x,d$y, col=col, border=NA)
	abline(h=0, col="grey")
}


plot_priors <- function(nu,alpha=NULL,delta=NULL){
	hist(nu,breaks=1000, ylim=c(0,6000), xlab="", main="")
	mtext(expression(nu),1, line=3, cex=1.5)
	if(!is.null(alpha)){
		delta = alpha / sqrt(1 + alpha^2)
	}else{
		alpha = delta/(sqrt(1-delta^2))
	}
	hist(alpha,breaks=1000, ylim=c(0,20000), xlab="", main="")
	mtext(expression(alpha),1, line=3, cex=1.5)
	hist(delta,breaks=1000, ylim=c(0,20000), xlab="", main="")
	mtext(expression(delta),1, line=3, cex=1.5)
	hist(get_gamma(alpha=alpha,delta=delta,nu=nu),breaks=1000, xlab="", main="", xlim=c(-4,4),ylim=c(0,40000))
	mtext(expression(gamma),1, line=3, cex=1.5)
}


n <- 1000000

setEPS()
pdf(paste0(wd,"R/plots/figure_SM_skew_priors.pdf"), , height=7, width=13)
{
set.seed(255)
par(mfrow=c(3,4), mar=c(5,5,1,1))
plot_priors(nu = runif(n,4,40), alpha = rnorm(n,0,10))
plot_priors(nu = runif(n,4,40),alpha = rnorm(n,0,1))
plot_priors(nu = runif(n,4,40),delta = runif(n,-1,1))

}
dev.off()


hist(rbeta(10000,0.5,0.5)*2-1, breaks=100)

par(mfrow=c(4,4), mar=c(5,5,1,1))
plot_priors(nu = runif(n,4,40), alpha = rnorm(n,0,10))
plot_priors(nu = runif(n,4,40),alpha = rnorm(n,0,1))
plot_priors(nu = runif(n,4,40),delta = runif(n,-1,1))
plot_priors(nu = runif(n,4,40),delta = rbeta(n,0.6,0.6)*2-1)


par(mfrow=c(5,4))
plot_priors(nu = runif(n,4,40), alpha = rnorm(n,0,10))
plot_priors(nu = runif(n,4,40),alpha = rnorm(n,0,1))
plot_priors(nu = exp(runif(n,log(4),log(40))),alpha = rnorm(n,0,1))
plot_priors(nu = runif(n,4,40),delta = runif(n,-1,1))
plot_priors(nu = exp(runif(n,log(4),log(40))),delta = runif(n,-1,1))

nu = runif(n,4,40)
delta = runif(n,-1,1)
alpha = delta/(sqrt(1-delta^2))
	hist(alpha,breaks=10000, main="alpha",xlim=c(-10,10))


plot_priors(nu = runif(n,4,40), alpha = rnorm(n,0,40))


hist(nu,breaks=1000, ylim=c(0,6000))
hist(alpha,breaks=1000)
hist(get_gamma(alpha,nu),breaks=1000)

nu <- exp(runif(n,log(4),log(40)));
delta <- runif(n,-1,1);
hist(nu,breaks=1000, ylim=c(0,6000))
hist(delta,breaks=1000)
hist(get_gammaD(delta,nu),breaks=1000)



# nu <- abs(rnorm(n,0,10))*-1+40;
# nu[nu<4]<-4
hist(get_gamma(alpha,nu),breaks=300)
density_plot(get_gamma(alpha,nu),col="blue")



nu <- runif(n,4,40);
#delta <- runif(n,-1,1);
 delta <- rbeta(n,1,1)*2-1;
#hist(get_gammaD(delta,nu),breaks=300, ylim=c(0,80000))
plot(density(get_gammaD(delta,nu)))

nu <- exp(runif(n,log(4),log(40)));
delta <- runif(n,-1,1);
#hist(get_gammaD(delta,nu),breaks=300, ylim=c(0,80000))
lines(density(get_gammaD(delta,nu)),col="red")
density_plot(get_gamma(delta,nu),col="blue")
lines(density(get_gamma(alpha,nu)),col="red")

hist(

density_plot(get_gammaD(delta,nu))
hist(delta,breaks=300)
hist(nu,breaks=300)

abline(v=4,col="red")
abline(h=0,col="grey")


https://discourse.mc-stan.org/t/ideas-for-modeling-systematically-skewed-outliers/636/6