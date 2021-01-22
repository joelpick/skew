

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
	hist(nu,breaks=1000, ylim=c(0,6000), main="nu")
	if(!is.null(alpha)){
		delta = alpha / sqrt(1 + alpha^2)
	}else{
		alpha = delta/(sqrt(1-delta^2))
	}
	hist(alpha,breaks=1000, ylim=c(0,6000), main="alpha")
	hist(delta,breaks=1000, ylim=c(0,6000), main="delta")
	hist(get_gamma(alpha=alpha,delta=delta,nu=nu),breaks=1000, main="gamma", xlim=c(-4,4),ylim=c(0,40000))
}


n <- 1000000

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