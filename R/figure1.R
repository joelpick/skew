rm(list=ls())


library(scales)
library(sn)

if(Sys.info()["user"]=="jhadfiel"){
	wd <- "..."
}else{
	wd <- "~/Dropbox/0_blue_tits/skew/"
}

source(paste0(wd,"R/functions.R"))

blankPlot <- function(ylim=c(-1,1),xlim=c(-1,1)) {
	op <- par(mar=c(0,0,0,0))
	plot(NULL,ylim=ylim,xlim=xlim, xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	par(op)
}

nlh2_plot <- function(G_skew,E_skew,h2,colG=NULL,colE=NULL){
	if(is.null(colG)) colG<-"grey"
	if(is.null(colE)) colE<-"grey"
	
	set.seed(24)
	G <- rskn(10000, mean=0,sd=sqrt(h2),skew=G_skew)
	E <- rskn(10000, mean=0,sd=sqrt(1-h2),skew=E_skew)

	P <- G + E
	G_d <- density(G, adjust=2.5 )
	E_d <- density(E, adjust=2.5 )
	P_d <- density(P, adjust=2.5 )
	P_skew <- round(add_skew(means=c(0,0),sds=sqrt(c(h2,1-h2)),skews=c(G_skew,E_skew))[3],2)
	ss <- smooth.spline(P,G)

	x<-seq(-4,4,length.out=1000)
	lims <- c(-4,4)
	ylims1 <- c(0,0.65)
	ylims2 <- c(0,2250)
	txt_pos <- 0.1
	txt_cex <- 2.5

	plot(G_d$x,G_d$y, type="l", bty="n", xlim=lims, ylim=ylims1, xaxt="n", yaxt="n", xlab="", ylab="", cex.lab=txt_cex);polygon(G_d$x,G_d$y, col=colG, border=NA); text(0,txt_pos,"G", cex=txt_cex); mtext(G_skew, side=1, line=0.8, at=0, cex=txt_cex/1.5);abline(h=0, col="grey")
	plot(E_d$x,E_d$y, type="l", bty="n", xlim=lims, ylim=ylims1, xaxt="n", yaxt="n", xlab="", ylab="+", cex.lab=txt_cex); polygon(E_d$x,E_d$y, col=colE, border=NA); text(0,txt_pos,"E", cex=txt_cex); mtext(E_skew, side=1, line=0.8, at=0, cex=txt_cex/1.5);abline(h=0, col="grey")
	plot(P_d$x,P_d$y, type="l", bty="n", xlim=lims, ylim=ylims1, xaxt="n", yaxt="n", xlab="", ylab="II", cex.lab=txt_cex); polygon(P_d$x,P_d$y, col="grey", border=NA); text(0,txt_pos,"P", cex=txt_cex); mtext(P_skew, side=1, line=0.8, at=0, cex=txt_cex/1.5);abline(h=0, col="grey")
	plot(ss, bty="L", xlim=lims,ylim=lims, type="l", lwd=3, xaxt="n", yaxt="n", xlab="Parent", ylab="Offspring", cex.lab=1.5)
}


setEPS()
pdf(paste0(wd,"R/plots/figure1.pdf"), height=8, width=12)
{

nPlot <- 5
# layout(matrix(c(1:(5*(nPlot+1))),ncol=5,byrow=TRUE), height=c(1,rep(3,nPlot+1)),width=c(1,rep(3,4)))

# par(mar=c(1,0,0,0))
# blankPlot()
# blankPlot(); text(0,0,"Genetic",cex=2.5)
# blankPlot(); text(c(-1,0),c(0,0),c("+","Environ."),cex=2.5)
# blankPlot(); text(c(-1,0),c(0,0),c("=","Phenotype"),cex=2.5)
# blankPlot()#; text(0,0,"h2 (0.5)",cex=3)

layout(matrix(c(1:(5*(nPlot))),ncol=5,byrow=TRUE), width=c(1,rep(3,4)))
par(mar=c(2,2,1,0),mgp=c(0.5,1,0)) 
blankPlot(); text(0.5,0,c("A)"),cex=3);	mtext("Skew =", side=1, at=0.5, line=0.7, cex=1.5)
nlh2_plot(G_skew=0.0,E_skew=0.0,h2=0.5, colG=alpha("darkgreen",0.5), colE=alpha("yellow",0.5))
blankPlot(); text(0.5,0,c("B)"),cex=3);	mtext("Skew =", side=1, at=0.5, line=0.7, cex=1.5)
nlh2_plot(G_skew=-0.99,E_skew=0.0,h2=0.5, colG=alpha("darkgreen",0.5), colE=alpha("yellow",0.5))
blankPlot(); text(0.5,0,c("C)"),cex=3);	mtext("Skew =", side=1, at=0.5, line=0.7, cex=1.5)
nlh2_plot(G_skew=0.0,E_skew=-0.99,h2=0.5, colG=alpha("darkgreen",0.5), colE=alpha("yellow",0.5))
blankPlot(); text(0.5,0,c("D)"),cex=3);	mtext("Skew =", side=1, at=0.5, line=0.7, cex=1.5)
nlh2_plot(G_skew=-0.99,E_skew=-0.99,h2=0.5, colG=alpha("darkgreen",0.5), colE=alpha("yellow",0.5))
blankPlot(); text(0.5,0,c("E)"),cex=3);	mtext("Skew =", side=1, at=0.5, line=0.7, cex=1.5)
nlh2_plot(G_skew=-0.99,E_skew=0.99,h2=0.5, colG=alpha("darkgreen",0.5), colE=alpha("yellow",0.5))

}
dev.off()
