
########################################
########################################
####    General data wrangling functions
########################################
########################################

#################################################
#  
#################################################
rbind_notAnnoying <- function(..., names=NULL){
  x <- list(...)
  y <- lapply(x, function(y){
    names(y) <- if(is.null(names)) names(x[[1]]) else names
    return(y)  
  } )
  do.call(rbind,y)
}

#################################################
#  
#################################################
lapply_dataFrame <- function(X,FUN) {
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, c(y, make.row.names=FALSE))
}

#################################################
#  
#################################################
groupFunc_dataFrame <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, c(y, make.row.names=FALSE))
}

#################################################
#  
#################################################
groupFunc_vector <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, y)
}

#################################################
#  
#################################################
groupFunc_vectorC <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(c, y)
}


#################################################
#  Number of unique levels in vector 			  
#################################################
n_unique <- function(x) length(unique(x))

#################################################
#  Makes vector of Mean and 95% CI from vector  
#################################################
mean_CI <- function(x) c(mean(x),quantile(x, c(0.025,0.975)))

#################################################
#  SE of vector 								  
#################################################
se <- function(x) sd(x)/sqrt(length(x))


#################################################
#  divide vector into bins								  
#################################################
bin <- function(x, n=10){
	breaks<-seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE)*1.0001,length.out=n+1)
	xCat <- min(breaks) + (abs(breaks[2]-breaks[1]))*sapply(x,function(y) sum(y>=breaks)-0.5)
	xCat
}



########################################
# calculate days from minimum date in vector
########################################
days_from_minimum <- function(dates,format="%F") as.numeric(strptime(dates,format=format) - min(strptime(dates,format=format)))/c(60*60*24)

########################################
# calculate days after 1st april
########################################
april_day <- function(date, format="%F") (as.numeric(strptime(date,format=format)) - as.numeric(strptime(paste0(format(strptime(date,format=format),"%Y"),"-04-01"),format="%F")))/c(60*60*24) + 1

########################################
# turns april day back into date
########################################
un_april_day <- function(april_date, year) as.character(strptime(paste0(year,"-04-01"),format="%F") + (april_date-1)*(60*60*24))







########################################
########################################
####    Blue tit database functions
########################################
########################################



########################################
# works out which nests to exclude from blue tit database
########################################
nest_to_exclude <- function(tMORPH, tNEST_SURVEY){
## nests where eggs not found before chicks
	chicks_no_eggs <- sapply(split(tNEST_SURVEY, tNEST_SURVEY$nest), function(x) "C" %in% x$state & !"E" %in% x$state)
## exclude nest where ages are not in c(0,1,3,6,9,12,15)
	chicks <- subset(tMORPH, is.na(agecode))
	ages <- c(0,1,3,6,9,12,15)
	wrong_ages <- sapply(split(chicks, chicks$nest), function(x){
		x$nest_age <- days_from_minimum(x$date)
		nest_ages <- unique(x$nest_age)
		incorrect_ages <- any(!nest_ages %in% ages)
		incorrect_seq <- if(!incorrect_ages){
			!length(nest_ages)==which(max(nest_ages)==ages)
		}else{TRUE}
		return(incorrect_seq)
	})
## excluded
	unique(c(names(chicks_no_eggs)[chicks_no_eggs], names(wrong_ages)[wrong_ages]))	
}

########################################
# works out hatch date of nest from blue tit database
########################################
april_hatch_date <- function(dates, format="%F") april_day( min(strptime(dates,format=format)))

########################################
# works out  male presence codes from blue tit database
########################################
male_presence <- function(tNEST_SURVEY,tMORPH){
	## assign whether a male was seen or not at each nest, from nest surveys, during chick phase
	males <- groupFunc_dataFrame(tNEST_SURVEY,"nest", function(x){
			x$nBirds <- 
			ifelse(is.na(x$observation) & is.na(x$birds), 0,
			ifelse(!is.na(x$observation) & is.na(x$birds), 1,
			ifelse(!is.na(x$observation) & nchar(x$birds)>=2, 2,
			nchar(x$birds)/2)))
			x$male_seen <- grepl("M",x$birds) | x$nBirds>1
			chicks <- subset(x, state=="C")
			male_seen <-if(nrow(chicks)>=1){ data.frame(nest=unique(chicks$nest), male_seen=sum(chicks$male_seen)>0)}else{NULL}
			return(male_seen)
		})

	## male presence
	tMORPH_M <- subset(tMORPH, agecode %in% c(5,6) & !grepl("_",bird_id) & morph_sex=="M")

	## assign whether a male was caught or not at each nest
	males$male_caught <- males$nest %in% unique(tMORPH_M$nest)

	males$male_present <- 
		ifelse(males$male_caught,2,
		ifelse(males$male_seen & !males$male_caught, 1, 0
	))

	return(males[,c("nest","male_present")])
}

########################################
# works out clutch size for each nest from blue tit database
########################################
clutch_size <- function(tEGGS){
	# tEGGS$broken<-tEGGS$treatment=="B"
	clutch_size <- aggregate(egg_id ~ nest_rear,tEGGS,length)
	names(clutch_size)<-c("nest","clutch_size")
	return(clutch_size)
}






########################################
########################################
####    Skew functions
########################################
########################################


#################################################
# calculate skew as third central moment
#################################################
skew <- function(x, na.rm=TRUE, sample=FALSE){
	if(na.rm) x <- na.omit(x) 
	n <- length(x)
	mu <- mean(x)
	sk <- mean((x-mu)^3)
	if(sample) sk <- sk * (n^2)/((n-2)*(n-1))
	return(sk)
	}

#################################################
# calculate skew as third standardised moment
#################################################
stand_skew <- function(x, na.rm=TRUE, sample=FALSE, correction=1){ 
	if(na.rm) x <- na.omit(x) 
	n <- length(x)
	mu <- mean(x)
    sigma <- sqrt(sum((x-mu)^2)/n)
    g1 <- mean(((x-mu)/sigma)^3)
	if(sample) {
		if(correction==1) g1 <- g1 * (sqrt(n*(n-1)))/(n-2) 
		if(correction==2) g1 <- g1 * ((sqrt(n*(n-1))/(n-2) )* ((n^2)/((n-2)*(n-1))))
	}
	#(n^2.5)/(((n-2)^2)*((n-1)^0.5))
	#((n^2)/((n-2)*(n-1)))^1.5
	#(sqrt(n*(n-1)))/(n-2) 
		#(n^2*g1)/((n-1)*(n-2)) }
	return(g1)
	}


#################################################
# mean, sd and skew tha result from adding variables 
# of different distributions
#################################################
add_skew <- function(means, sds, skews){
	new_mean <- sum(means)
	new_sd <- sqrt(sum(sds^2))
	new_skew <- sum(skews*sds^3)/new_sd^3
	return(c(new_mean,new_sd,new_skew))
}



skew_param <- function(xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){

	b <- sqrt(2/pi)
 	b2 <- (b^2 -1)
 	a2 <- (2*b^3-b)^(2/3)

	if((is.null(xi)&is.null(mean))) stop("Location parameter (xi or mean) needed")
	if((is.null(omega)&is.null(sd))) stop("Scale parameter (omega or sd) needed")
	if(	(is.null(alpha)&is.null(delta)&is.null(skew)))	stop("Skew parameter (alpha, delta or skew) needed")

	if((!is.null(xi)&!is.null(mean))) warning("2 location parameters given, using xi")
	if((!is.null(omega)&!is.null(sd))) warning("2 scale parameters given, using omega")
	if((!is.null(alpha)&!is.null(delta)&!is.null(skew))) warning("Multiple skew parameters given; alpha is used over delta, and delta over skew")

	## work out alpha, if not given
	## from delta
	if(is.null(alpha)&!is.null(delta)){
		if(delta>1 | delta< -1) stop("Delta must be between -1 and 1")
		alpha <- delta / sqrt(1 - delta^2)
	}

	## from skew
	if(is.null(alpha)&is.null(delta)){
		if(skew>0.999 | skew< -0.999) stop("Skew must be between -0.999 and 0.999")
		alpha <- sign(skew)*abs(skew)^(1/3)/sqrt(a2 + b2*abs(skew)^(2/3))
	}

	## calculate delta
	delta <- alpha / sqrt(1 + alpha^2)

	## work out omega if not given
	if(is.null(omega)) omega <- sd / sqrt(1-(b*delta)^2)
	
	## work out xi if not given
	if(is.null(xi)) xi <- mean - omega * delta * b	

	# b <- sqrt(2/pi)
	# b2 <- (b^2 -1)
	# a2 <- (2*b^3-b)^(2/3)
	# alpha <- sign(skew)*abs(skew)^(1/3)/sqrt(a2 + b2*abs(skew)^(2/3))
	# delta <- alpha / sqrt(1 + alpha^2)
	# omega <- sd / sqrt(1-(b*delta)^2)
	# xi <- mean - omega * delta * b
  return(list(xi=xi,omega=omega,alpha=alpha))
}







#################################################
#  Simulate vector of random effects            
#################################################
simRE <- function(levels,var){

  dMat <- model.matrix(~ as.character(levels)-1)
  N <- length(unique(levels)) 
  pMat <- rnorm(N,0,sqrt(var))
  return(dMat %*% pMat)
}


#################################################
#  Simulate skew normal             
#################################################
rskn <- function(n, xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){
	do.call(rsn,c(list(n=n), skew_param(xi=xi, omega=omega, alpha=alpha, delta=delta, mean=mean, sd=sd, skew=skew)))
}

dskn <- function(x, xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){
	do.call(dsn,c(list(x=x), skew_param(xi=xi, omega=omega, alpha=alpha, delta=delta, mean=mean, sd=sd, skew=skew)))
}

#################################################
#  Simulate vector of random effects from skew normal distribution            
#################################################
simRE_skew <- function(levels, xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){
	dMat <- model.matrix(~ as.character(levels)-1)
	pMat <- matrix(rskn(n=length(unique(levels)), xi=xi, omega=omega, alpha=alpha, delta=delta, mean=mean, sd=sd, skew=skew), ncol=1) 
	return(dMat %*% pMat)
	}


mskn<-function(xi, omega, alpha){
delta<-alpha/sqrt(1+alpha^2)
xi+omega*delta*sqrt(2/pi)
}

vskn<-function(omega, alpha){
delta<-alpha/sqrt(1+alpha^2)
(omega^2)*(1-(2*delta^2)/pi)
}

sskn<-function(omega, alpha=NULL, delta=NULL, standardised=TRUE){
  if(is.null(delta)) delta<-alpha/sqrt(1+alpha^2)
  sk<-((4-pi)/2)*((delta*sqrt(2/pi))^3)
  if(!standardised){
    sk<-sk*omega^3
  }else{
    sk<-sk/(1-2*(delta^2)/pi)^(3/2)
  }  
  return(sk)  
}

mskt<-function(xi=0, omega=1, alpha=0, nu=Inf, dp=NULL){
	# from Azzalini & Capitanio 203
	if(!is.null(dp)){
		if(length(dp)!=4){stop("dp should be of length 4")}
		xi<-dp[1]
		omega<-dp[2]
		alpha<-dp[3]
		nu<-dp[4]
	}
    delta<-alpha/sqrt(1+alpha^2)
	mu<-delta*((nu/pi)^(0.5))*gamma(0.5*(nu-1))/gamma(0.5*nu)
	xi+omega*mu
}



########################################
########################################
####    Stan functions
########################################
########################################


###############################################
#    extract parameters from rstan model       
###############################################
pars <- function(model, par){
	if(class(model)=="stanfit"){
		x <- summary(model)$summary
		x <- x[grep(par, rownames(x)),c(1,4,8), drop=FALSE]
	}else if(class(model)=="mcmc.list"){
		x <- as.data.frame(posterior::summarize_draws(posterior::as_draws_array(model),"mean",~quantile(.x, probs = c(0.025, 0.975))))[-(1:7),]
		rownames(x) <- x$variable
		x <- x[grep(par, rownames(x)),colnames(x)!="variable"]
  	}else{	
  		stop("model is not class stanfit or mcmc.list")
	}
	return(x)
}


###############################################
#  extract skew T parameters from rstan model  
###############################################
pars_ST <- function(model, variable){
	if(class(model)=="stanfit"){
		x <- summary(model)$summary[paste(c("xi","omega", "alpha", "nu"),variable, sep="_"),c(1,4,8)]
	}else if(class(model)=="mcmc.list"){
		x <- as.data.frame(posterior::summarize_draws(posterior::as_draws_array(model_z),"mean",~quantile(.x, probs = c(0.025, 0.975))))#[-(1:7),]
		rownames(x) <- x$variable
		x <- x[paste(c("xi","omega", "alpha", "nu"),variable, sep="_"),colnames(x)!="variable"]
	}else{
		stop("model is not class stanfit or mcmc.list")
	}
	return(x)
}


###############################################
#  factorise pedigree, IDs refer to row numbers
###############################################
factorisePed <- function(pedigree, unknown=0){
    new_ped <- data.frame(
        1:nrow(pedigree), 
        ifelse(is.na(pedigree[,2]),unknown,match(pedigree[,2], pedigree[,1])), 
        ifelse(is.na(pedigree[,3]),unknown,match(pedigree[,3], pedigree[,1]))
        )
    colnames(new_ped) <- colnames(pedigree)[1:3]

    return(new_ped)
}

###############################################
#  make starting values based on equivalent models in asreml
###############################################
makeSV <- function(asreml_data,stan_data, animal=TRUE, ME=TRUE, delta=FALSE) {

	mod <- if(animal){
		asreml_data$mod_A
	}else{
		asreml_data$mod_DS
	}

	sigma_vars <- colnames(mod$Vi)
	sigma_starting <- sqrt(rmvnorm(n = 1,
		mean = mod$V,
		sigma = as.matrix(mod$Vi)
		))
	beta_starting <- rmvnorm(n = 1,
		mean = mod$F,
		sigma = as.matrix(mod$Fi)
		)
	## beta_tilde from QR decomposition
	beta_tilde_starting <- as.vector((qr.R(qr(stan_data$X))*1/sqrt(nrow(stan_data$X)-1)) %*% t(beta_starting))

	SV <- list(
		beta_tilde = beta_tilde_starting,

		sigma_nest = sigma_starting[,"nest"],
		nu_nest = rnorm(1,30,0.1),
		#log_nu_nest = log(30),
		nest_effects = rnorm(stan_data$N_nest,0,sigma_starting[,"nest"])
	)	

	if(delta){
		SV$delta_nest = rnorm(1,0,0.1)
	}else{
		SV$alpha_nest = rnorm(1,0,0.1)
	}

	if(animal){
		SV$sigma_A <- sigma_starting[,"animal"]
		SV$A_NoParents <- rnorm(stan_data$N_NoParents,0,1)
		SV$A_ParentsOffspring <- rnorm(stan_data$N_ParentsOffspring,0,1)
		SV$A_ParentsNoOffspring <- rnorm(stan_data$N_ParentsNoOffspring,0,1)
	}else{
		SV$sigma_dam_sire <- sigma_starting[,"dam"]
		#SV$log_nu_dam_sire <- log(30)
		SV$nu_dam_sire <- rnorm(1,30,0.1)
		SV$dam_sire_effects <- rnorm(stan_data$N_dam_sire,0,sigma_starting[,"dam"])
		if(delta){
			SV$delta_dam_sire <- rnorm(1,0,0.1)
		}else{
			SV$alpha_dam_sire <- rnorm(1,0,0.1)
		}
	}
	
	if(ME){
		SV$sigma_ind <- sigma_starting[,"bird_id"]
		#SV$log_nu_ind <- log(30)
		SV$nu_ind <- rnorm(1,30,0.1)
		SV$ind_effects <- rnorm(stan_data$N_ind,0,sigma_starting[,"bird_id"])
		
		if(delta){
			SV$delta_ind <- rnorm(1,0,0.1)
		}else{
			SV$alpha_ind <- rnorm(1,0,0.1)
		}

		SV$sigma_E <- sigma_starting[,"units"]


	}else{
		SV$sigma_E <- sigma_starting[,"units"]
		#SV$log_nu_E <- log(30)
		SV$nu_E <- rnorm(1,30,0.1)

		if(delta){
			SV$delta_E <- rnorm(1,0,0.1)
		}else{
			SV$alpha_E <- rnorm(1,0,0.1)
		}
	}

	return(SV)
}




make_stan_dat <- function(fixed, random, pedigree=NULL, data){
	F_vars <- all.vars(fixed)
	R_vars <- all.vars(random, unique = FALSE)

	if(!is.null(pedigree) ){
		R_terms <- attr(x = terms(random), which = "term.labels")
		mm_terms <- grepl("mm\\(",R_terms)
		mm <- any(grepl("mm\\(",R_terms))
		animal_terms <- NULL
		for(i in 1:length(R_terms)){
			 y<- if(mm_terms[i]) {rep(R_terms[i],2)} else {R_terms[i]}
			 animal_terms <- c(animal_terms,grepl("animal\\(",y))
			}
		animal_term <- R_vars[animal_terms]
	}
	
	stan_dat <- list(N= nrow(data), y=as.numeric(data[,F_vars[1]]))
	
	for(i in which(!animal_terms)){
		stan_dat[[paste0("N_",R_vars[i])]] <- length(unique(data[,R_vars[i]]))
		stan_dat[[paste0(R_vars[i],"_id")]] <- as.numeric(as.factor(as.character(data[,R_vars[i]])))
	}

	if(!is.null(pedigree)){
		stan_ped <- factorisePed(pedigree)
		stan_dat[["N_ped"]] <- nrow(stan_ped)
		stan_dat[[paste0(animal_term,"_Ped_id")]] <- stan_ped[match(data[,animal_term], pedigree[,1]),1]
		stan_dat[["dam_Ped_id"]] <- stan_ped[,2]
		stan_dat[["sire_Ped_id"]] <- stan_ped[,3]

		invA <- inverseA(pedigree)
		stan_dat[["MSV"]] <- invA$dii
	}
	return(stan_dat)
}




########################################
########################################
####    ASReml functions
########################################
########################################


#################################################
## function to get variances out of asreml object
#################################################
var_out <- function(mod){
	x <- mod$gammas
	names(x) <- sub("!.*\\..*$", "", names(x))
	names(x) <- sub(".*\\((\\w*)\\).*$", "\\1", names(x))
	x[names(x)!="R!variance"]
}


pin<-function (object, transform) {
    pframe <- as.list(object$gammas)
    names(pframe) <- sub("!.*\\..*$", "", names(pframe))
	names(pframe) <- sub(".*\\((\\w*)\\).*$", "\\1", names(pframe))
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
        pframe)
    X <- as.vector(attr(tvalue, "gradient"))
    tname <- if (length(transform) == 3) 
        transform[[2]]
    else ""
    Vmat <- object$ai
    n <- length(pframe)
    i <- rep(1:n, 1:n)
    j <- sequence(1:n)
    k <- 1 + (i > j)
    se <- sqrt(sum(Vmat * X[i] * X[j] * k))
    data.frame(row.names = tname, Estimate = tvalue, SE = se)
}





########################################
########################################
####    Plotting functions
########################################
########################################

#################################################
#   Makes a forest plot.                        #
#   Takes matrix of mean, 0.025 and 0.975 	  #
#################################################
effectPlot <- function(x, y=NULL,add=FALSE,offset=NULL, fixed=NULL, xlab="", ylab= "", xlim=NULL,yaxis=2, col="black", cex.axis=1,pch=19,names=NULL,...){ 

	if(is.null(y)) y <- nrow(x):1
	y_offset <- if(is.null(offset)){ y }else{ y + offset }
	if(is.null(xlim)){
		xlim<- if(min(x)>0 & max(x)>0) {c(0,max(x))
		}else if(min(x)<0&max(x)<0) {c(min(x),0)
		}else{	range(x)}}
	if(is.null(names)) names <- rownames(x)
	if(!add){	
		plot(x[,1], y, xlab=xlab, ylab=ylab, xlim=xlim, col=0, yaxt="n", ylim=c(0.5,nrow(x)+0.5), cex.axis=cex.axis,...)#c(min(y)-0.5,max(x)+0.5)
		axis(yaxis,y,names,cex.axis=cex.axis*1.5, las=1)
		abline(v=0, col="grey")
	}
	points(x[,1], y_offset, pch=pch, col=col, cex=1.5)
	arrows(x[,2], y_offset, x[,3], y_offset, code=3, angle=180, length=0.01,col=col, lwd=1.5)
	if(!is.null(fixed)){
		points(x[fixed,1], y_offset[fixed], col="red", cex=1.75)

	}
}

#################################################
## function to make a blank plot
#################################################
blankPlot <- function(ylim=c(-1,1),xlim=c(-1,1)) {
	op <- par(mar=c(0,0,0,0))
	plot(NULL,ylim=ylim,xlim=xlim, xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	par(op)
}

#################################################
## function to make a blank plot
#################################################
binPlot <- function(formula,data=NULL,text=TRUE,text.cex=1,...){
	bindedMeans<- aggregate(formula,data,mean)
	bindedCounts<- aggregate(formula,data,length)
	plot(formula,bindedMeans, pch=19,...)
	if(text) text(formula,bindedMeans, bindedCounts[,2], pos=3,cex=text.cex)
}

#################################################
## function to make a blank plot
#################################################
binPoints <- function(formula,data,text.cex=1,...){
	bindedMeans<- aggregate(formula,data,mean)
	bindedCounts<- aggregate(formula,data,length)
	points(formula,bindedMeans, pch=19,...)
	text(formula,bindedMeans, bindedCounts[,2], pos=3,cex=text.cex)
}



#################################################
## function to extract posterior distributions relevant parameters from survival model for calculation of selection gradients
## the fixedEffects argument single trait version of the fixed effect formula used in the survival model (excluding the trait)
##Creates a list with:
#### beta_all and gamma_all - linear and quadratic trait coefficients from survival model, respectively
#### eta2_all and eta2_all - non-trait fixed effects predictions from survival model
#### V_nest_all - nest covariance matrix from survival model
#################################################

extract_w <- function(model_w, trait, fixedEffects, data){

	beta_pos<-grep(paste0(trait, "C:|", trait, "C$"), colnames(model_w$Sol))
	# linear coefficients for trait values from the survival model

	gamma_pos<-grep(paste0(trait, "C2:|", trait, "C2$"), colnames(model_w$Sol))
	# quadratic coefficients for trait values from the survival model

	eta1_pos<-grep("at\\.level\\(age,\\ 1\\)", colnames(model_w$Sol))
	eta1_pos<-setdiff(eta1_pos, c(beta_pos, gamma_pos))
	# positions of non-trait terms in the survival model that pertain to Day15 to fledging. 

	eta2_pos<-grep("at\\.level\\(age,\\ 2\\)", colnames(model_w$Sol))
	eta2_pos<-setdiff(eta2_pos, c(beta_pos, gamma_pos))
	# positions of non-trait terms in the survival model that pertain to fledging to recruitment. 

	X_eta<-model.matrix(fixedEffects, THBW_noRep)
	# design matrix for non-trait effects in survival model


	if(any(!mapply(colnames(X_eta), colnames(model_w$Sol)[eta1_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta1")}

	if(any(!mapply(colnames(X_eta), colnames(model_w$Sol)[eta2_pos], FUN=function(x,y){grepl(x,y)})[-1])){stop("design matrices and model effects don't match for eta2")}

	eta1_all<-apply(model_w$Sol[,eta1_pos], 1, function(x) X_eta%*%x)
	eta2_all<-apply(model_w$Sol[,eta2_pos], 1, function(x) X_eta%*%x)
    # non-trait fixed effects predictions from survival model

	beta_all <- model_w$Sol[,beta_pos]
	gamma_all <- model_w$Sol[,gamma_pos]

    V_nest_all <- model_w$VCV[,grep("nest", colnames(model_w$VCV))]

	list(beta_all=beta_all, gamma_all=gamma_all, eta1_all=eta1_all, eta2_all=eta2_all, V_nest_all=V_nest_all)
}








###############################################
#    Sample mean and 2nd-4th central moments  #
###############################################
scm<-function(x, standardise=FALSE){

  mu<-1:4
  mu[1]<-mean(x)
  for(i in 2:4){
    mu[i]<-mean((x-mu[1])^i)
  }
  if(standardise){
  	for(i in 3:4){
       mu[i]<-mu[i]/sqrt(mu[2])^i
    }   
  }	
  return(mu)
}  


################################################################
# Computes mean and 2nd-4th central moments
# if dp is a list it is the moments of the sum
################################################################
dp2cm<-function(dp, family, object = NULL, cp.type = "proper", upto = NULL, standardise=FALSE){

  if(!is.list(dp)){
  	dp<-list(dp)
  }

  stdm<-lapply(dp, dp2cp, family=family, object = object, cp.type = cp.type, upto = upto)

  cm<-lapply(stdm, function(x){c(x[1], x[2]^2, x[3]*x[2]^3, (x[4]+3)*x[2]^4)})

  mu = 1:4

  mu[1]<-sum(unlist(lapply(cm, function(x){x[1]})))
  mu[2]<-sum(unlist(lapply(cm, function(x){x[2]})))
  mu[3]<-sum(unlist(lapply(cm, function(x){x[3]})))
  # mean and 2nd + 3rd central moments are additive

  # mu[4] = k[4] + 3*mu[2]^2 where k[4] is the fourth cumulamnt of the sum. 
  # The 4th cumulant of component i is k_i[4] = cm_i[4]-3*cm_i[2]^2
  # Cumulants are additive so k[4] = sum(cm[4]-3*cm[2]^2)

  mu[4] = sum(unlist(lapply(cm, function(x){x[4]-3*x[2]^2})))+3*mu[2]^2

  if(standardise){
  	for(i in 3:4){
       mu[i]<-mu[i]/sqrt(mu[2])^i
    }   
  }	

  return(mu)
}


###############################################################################
# Computes the linear regression coefficient from a quadratic fit
# of fitness on the trait, where mu are the mean and central moments (2:4),
# S is the linear and C the quadratic selection differential 			
###############################################################################
betaLA_2<-function(mu, S, C, family, object = NULL, cp.type = "proper", upto = NULL){

	return(((mu[4]-mu[2]^2)*S-mu[3]*C)/(mu[2]*(mu[4]-mu[2]^2)-mu[3]^2))

}


################################################################
# Function for plotting a histogram of x and overlaying 
# a skew-t density function with parameters dp 
################################################################
comp_skt<-function(x, dp, breaks="Sturges", ...){

  plot.points<-hist(x, breaks=breaks)

  plot.points$breaks[1]<--Inf
  plot.points$breaks[length(plot.points$breaks)]<-Inf

  pz<-diff(pst(plot.points$breaks, dp=dp))

  lines(pz*sum(!is.na(x))~plot.points$mids, ...)

} 


################################################################
# Function for obtaining the integrand in the convolution
# \int Pr(z_p-e_p | \theta_g)*Pr(e_p| \theta_e) de_p 
# the distribution of genetic and environmental effects
# are passed as parameters of the skew-t
################################################################
conv<-function(par, z_p, g_st, e_st){

  if(length(par)!=length(e_st)){stop("par should be of length 2")}

  res<-dst(z_p-sum(par), dp=g_st)
  # Pr(z_p-e_p | \theta_g)

  for(i in 1:length(e_st)){
     res<-res*dst(par[i], dp=e_st[[i]])
  }

  # Pr(z_p-e_p | \theta_g)*Pr(e_p| \theta_e)

  return(res)
}  


wconv<-function(par, z_p, g_st, e_st){

  #########################################################################
  #   Function for obtaining the integrand in the weighted convolution    #
  #     \int (z_p-e_p)* Pr(z_p-e_p | \theta_g)*Pr(e_p| \theta_e) de_p     #
  #      the distribution of genetic and environmental effects            #
  #          are passed as parameters of the skew-t                       #
  #########################################################################

  (z_p-sum(par))*conv(par, z_p, g_st, e_st)

}

wconv2<-function(par, z_p, g_st, e_st, limit_prob=1e-4){

  #########################################################################
  #   Function for obtaining the integrand in the weighted convolution    #
  #    \int (z_p-e_p)^2* Pr(z_p-e_p | \theta_g)*Pr(e_p| \theta_e) de_p    #
  #      the distribution of genetic and environmental effects            #
  #          are passed as parameters of the skew-t                       #
  #########################################################################

  ((z_p-sum(par))^2)*conv(par, z_p, g_st, e_st)

}


dz<-function(z_p, g_st, e_st, limit_prob=1e-4){

  ################################################################
  #               Function for calculating Pr(z)                 #
  #  based on distribution of genetic and environmental effects  #
  #          (passed as parameters of the skew-t)                #
  ################################################################

  e_llimits<-unlist(lapply(e_st, function(x){qst(limit_prob, dp=x)}))
  e_ulimits<-unlist(lapply(e_st, function(x){qst(1-limit_prob, dp=x)}))
  # lower and upper limits when integrating over environmental effects

  dz<-hcubature(conv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral

  return(dz)
}

rz<-function(n, z_st){

  ################################################################
  #               Function for simulating z                      #
  #  based on distribution of genetic and environmental effects  #
  #          (passed as parameters of the skew-t)                #
  ################################################################

  u<-data.frame(c(lapply(z_st, function(x){rst(n, dp=x)})))

  return(rowSums(u))
}


################################################################
#  Function for calculating E[g|z] based on distribution of 
#  genetic and environmental effects (passed as parameters of the skew-t)
################################################################

POreg<-function(z_p, g_st, e_st, limit_prob=1e-4){

  e_llimits<-unlist(lapply(e_st, function(x){qst(limit_prob, dp=x)}))
  e_ulimits<-unlist(lapply(e_st, function(x){qst(1-limit_prob, dp=x)}))
  # lower and upper limits when integrating over environmental effects

  Eg_p<-hcubature(wconv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
  Eg_p<-0.5*Eg_p/hcubature(conv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral

  return(Eg_p)
}



################################################################
#  Function for calculating \partial E[g|z] \partial z     
#  based on distribution of genetic and environmental effects  
#  (passed as parameters of the skew-t)                
################################################################

dPOreg<-function(z_p=NULL, g_st=NULL, e_st=NULL, Eg_p=NULL, limit_prob=1e-4){

  e_llimits<-unlist(lapply(e_st, function(x){qst(limit_prob, dp=x)}))
  e_ulimits<-unlist(lapply(e_st, function(x){qst(1-limit_prob, dp=x)}))
  # lower and upper limits when integrating over environmental effects

  if(is.null(Eg_p)){
     Eg_p<-hcubature(wconv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
     Eg_p2<-hcubature(wconv2, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
     E_p<-hcubature(conv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
     Eg_p<-Eg_p/E_p
     Eg_p2<-Eg_p2/E_p
  }else{
    Eg_p<-2*Eg_p
    Eg_p2<-hcubature(wconv2, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
    Eg_p2<-Eg_p2/hcubature(conv, e_llimits, e_ulimits, z_p=z_p, g_st=g_st, e_st=e_st)$integral
  }
  dEg_p <- 0.5*(1+(Eg_p^2-Eg_p2)/(g_st[2]^2))

  return(dEg_p)
}


#######################################################################################
#  Function for calculating the conditional means and (co)variances
#  for the set of variates in keep_var conditional on the set of variates in cond_var 
#  cond are the values of the cond_var variates to be conditioned on
#  mean and sigma are the unconditional means and (co)variances
#######################################################################################

cmvnorm<-function(mean=NULL, sigma=NULL, cond=NULL, keep_var, cond_var=NULL){
		
	cV <- sigma[keep_var,keep_var]-sigma[keep_var, cond_var]%*%solve(sigma[cond_var, cond_var])%*%sigma[cond_var, keep_var]

	if(length(dim(cond))==2){
		cM <- sapply(1:nrow(cond), function(x) mean[keep_var]+sigma[keep_var, cond_var]%*%solve(sigma[cond_var, cond_var])%*%(cond[x,]-mean)[cond_var])
	}else{
		cM <- mean[keep_var]+sigma[keep_var, cond_var]%*%solve(sigma[cond_var, cond_var])%*%(cond-mean)[cond_var]
	}

	return(list(cM=cM, cV=cV))
}


############################################################################
# Fitness as a function of z for a probit two-event history analysis where 
## beta and gamma are vectors for the linear and quadratic coefficient for the two events
## mu_etaz and V_etaz are the means and covariances of the linear predictors (excluding z) for the two events followed by z
## V_nest is the between nest covariance matrix and the error structure is assumed to be diag(2)
############################################################################

w_func<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	c_s <- cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)
	g_s<-c_s$cM+beta*z+gamma*z^2
    V_s<-c_s$cV+V_nest+diag(2)

	return(pmvnorm(lower=c(0,0), mean=c(g_s), sigma=V_s))
}


#################################################################################################
# Derivative of fitness with respect to z for a probit two-event history analysis 
#  where beta and gamma are vectors for the linear and quardatic coefficient for the two events 
# mu_etaz and V_etaz are the means and covariances of the linear predictors (excluding z) for the two events followed by z
# V_nest is the between nest covrainace matrix and the error structure is assumed to be diag(2) 
#################################################################################################
wD_func<-function(z, mu_etaz, V_etaz, beta, gamma,  V_nest){

	c_s <- cmvnorm(mean=mu_etaz, sigma=V_etaz, cond=c(NA,NA, z), keep_var=1:2, cond_var=3)
	g_s<-c_s$cM+beta*z+gamma*z^2
    V_s<-c_s$cV+V_nest+diag(2)

    c_sc <- cmvnorm(mean=c(0,0), sigma=V_s, cond=c(g_s), keep_var=1, cond_var=2)
    g_sc<-g_s[1]-c_sc$cM

	V_sc<-c_sc$cV

    ch1<-dnorm(g_sc[1], 0, sqrt(V_sc))*pnorm(g_s[2], 0, sqrt(V_s[2,2]))*(V_etaz[1,3]/V_etaz[3,3]+beta[1]+2*gamma[1]*z-(V_s[1,2]/V_s[2,2])*(V_etaz[2,3]/V_etaz[3,3]+beta[2]+2*gamma[2]*z))
    ch2<-pnorm(g_sc[1], 0, sqrt(V_sc))*dnorm(g_s[2], 0, sqrt(V_s[2,2]))*(V_etaz[2,3]/V_etaz[3,3]+beta[2]+2*gamma[2]*z)

    return(ch1+ch2)
}


post_mu<-function(model, components=NULL, X=NULL, pred_pos=NULL, standardise=FALSE){

  ################################################################################################
  # Posterior distribution of mean and 2nd-4th central moments from a stan fit
  # for the sum of effects listed in components	and the linear predictor defined by X[,pred_pos] 
  # component terms are assumed to be either Gaussain or Skew-t
  ################################################################################################

    post<-as.data.frame(extract(model))

    post_mu<-matrix(NA, nrow(post), 4)
    z_st<-pos_comp<-rep(list(c(0,0,0, 1e+16)), length(components))

    ST<-paste0("xi", "_", components)%in%names(post)

    for(j in 1:length(components)){
      if(ST[j]){
        pos_comp[[j]]<-match(paste0(c("xi", "omega", "alpha", "nu"), "_", components[j]), names(post))
      }else{
        pos_comp[[j]]<-match(paste0("sigma_", components[j]), names(post))
      }  
    }

    for(i in 1:nrow(post)){
      for(j in 1:length(components)){
        if(ST[j]){
          z_st[[j]]<-unlist(post[i,pos_comp[[j]]])
        }else{
          z_st[[j]][2]<-post[i,pos_comp[[j]]]
        }  
      }
      post_mu[i,]<-dp2cm(z_st, family="ST")
    }

    if(!is.null(X)){
	    Xpred<-X[,pred_pos]

	    pred_pos<-grep("beta\\.", names(post))[pred_pos]
	 
	    post_mu<-post_mu+t(apply(post[,pred_pos], 1, function(x){scm(Xpred%*%x)}))
    }
    
    if(standardise){
	  	for(i in 3:4){
	       post_mu[,i]<-post_mu[,i]/sqrt(post_mu[,2])^i
	    }   
  }	


    return(post_mu)
}




###############################################################
#  h2 - Function for obtaining (2X) the single-parent-offspring regression after selection using simulation where n is the number of simulated values                                
#  *From the Trait Model*                                            
#  - g_st and e_st are the distributional parameters for the trait                           
#  - if adj_mean is not NULL the mean defined by g_st and e_st is replaced with the value specified  	 
#  *From the Survival Model*                                          
#  - mu_etaz and V_etaz are the mean and covariance matrix of trait values and linear predictors         
#	 - beta and gamma are the linear/quadratic regression terms for the (mean-centred) trait            
#	 - V_nest is the between nest covariance matrix                                
###############################################################
h2<-function(n, g_st=NULL, e_st=NULL, adj_mean=NULL, mu_etaz=NULL, V_etaz=NULL, beta=NULL, gamma=NULL,  V_nest=NULL, after=TRUE){

  mu<-dp2cm(c(list(g_st), e_st), family="ST")

  if(after){
	  g<-rst(n, dp=g_st)
	  zp<-g+rz(n, e_st)
	  zo<-0.5*g+rz(n, e_st)+rst(n, dp=g_st)*sqrt(0.75)

	  if(!is.null(adj_mean)){
	    adj_mean<-adj_mean-mu[1]
	    mu[1]<-mu[1]+adj_mean
	  }else{
	  	adj_mean<-0
	  }	
	 
	  zp<-zp+adj_mean
	  zo<-zo+adj_mean

	  wz<-sapply(zp, w_func, mu_etaz=mu_etaz, V_etaz=V_etaz, beta=beta, gamma=gamma,  V_nest=V_nest)
	  wz<-wz/mean(wz)

      h2<-2*(mean(wz*zp*zo)-mean(wz*zp)*mean(wz*zo))/mean((wz*zp^2)-mean(wz*zp)^2)

   }else{
   	  h2<-dp2cm(g_st, family="ST")[2]/mu[2]
   }	  

  return(h2)

}

