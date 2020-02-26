rbind_notAnnoying <- function(..., names=NULL){
  x <- list(...)
  y <- lapply(x, function(y){
    names(y) <- if(is.null(names)) names(x[[1]]) else names
    return(y)  
  } )
  do.call(rbind,y)
}

download_bluetit <- function(tables, dir){
	links <- list(
		tMORPH = "3gcvx",
		tNEST_SURVEY = "7dm4v",
		tEGGS = "bvfkg",
		tBIRDS = "958x7",
		tEXTRACTIONS = "68fz9",
		tGENOTYPES = "bmke2",
		tIB = "s75pj",
		tIBSTAT = "4dcgb"
		
	)
	links_to_save <- links[tables]
	for (i in 1:length(links_to_save)){
		Link <- paste0("https://osf.io/",links_to_save[[i]],"/download")
		FullPath <- paste0(dir,"/",names(links_to_save)[i],".csv")
		httr::GET(Link,httr::write_disk(FullPath, overwrite = TRUE))
	}
}


days_from_minimum <- function(dates,format="%F") as.numeric(strptime(dates,format=format) - min(strptime(dates,format=format)))/c(60*60*24)


april_day <- function(date, format="%F") (as.numeric(strptime(date,format=format)) - as.numeric(strptime(paste0(format(strptime(date,format=format),"%Y"),"-04-01"),format="%F")))/c(60*60*24) + 1

un_april_day <- function(april_date, year) as.character(strptime(paste0(year,"-04-01"),format="%F") + (april_date-1)*(60*60*24))



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

april_hatch_date <- function(dates, format="%F") april_day( min(strptime(dates,format=format)))


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


clutch_size <- function(tEGGS){
	# tEGGS$broken<-tEGGS$treatment=="B"
	clutch_size <- aggregate(egg_id ~ nest_rear,tEGGS,length)
	names(clutch_size)<-c("nest","clutch_size")
	return(clutch_size)
}


lapply_dataFrame <- function(X,FUN) {
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, c(y, make.row.names=FALSE))
}

groupFunc_dataFrame <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, c(y, make.row.names=FALSE))
}
groupFunc_vector <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(rbind, y)
}
groupFunc_vectorC <- function(dat, var, FUN) {
	X <- split(dat, dat[,var])
	y <- lapply(X=X,FUN=FUN)
	if(any(!sapply(y,is.null))) do.call(c, y)
}

n_unique <- function(x) length(unique(x))

mean_CI <- function(x) c(mean(x),quantile(x, c(0.025,0.975)))

## function to get variances out of asreml object
var_out <- function(mod){
	x <- mod$gammas
	names(x) <- sub("!.*\\..*$", "", names(x))
	names(x) <- sub(".*\\((\\w*)\\).*$", "\\1", names(x))
	x[names(x)!="R!variance"]
}


simRE <- function(levels, var){
	dMat <- model.matrix(~ as.character(levels)-1)
	pMat <- matrix(rnorm(length(unique(levels)), 0, sqrt(var)), ncol=1) 
	return(dMat %*% pMat)
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


rskn <- function(n, xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){
	do.call(rsn,c(list(n=n), skew_param(xi=xi, omega=omega, alpha=alpha, delta=delta, mean=mean, sd=sd, skew=skew)))
}

dskn <- function(x, xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){
	do.call(dsn,c(list(x=x), skew_param(xi=xi, omega=omega, alpha=alpha, delta=delta, mean=mean, sd=sd, skew=skew)))
}

simRE_skew <- function(levels, xi=NULL, omega=NULL, alpha=NULL, delta=NULL, mean=NULL, sd=NULL, skew=NULL){
	dMat <- model.matrix(~ as.character(levels)-1)
	pMat <- matrix(rskn(n=length(unique(levels)), xi=xi, omega=omega, alpha=alpha, delta=delta, mean=mean, sd=sd, skew=skew), ncol=1) 
	return(dMat %*% pMat)
	}


# calculate skew as third central moment
skew <- function(x, na.rm=TRUE){
	if(na.rm) x <- na.omit(x) 
	mu <- mean(x)
	mean((x-mu)^3)
	}

# calculate skew as third standardised moment
stand_skew <- function(x, na.rm=TRUE){ 
	if(na.rm) x <- na.omit(x) 
	mu <- mean(x)
	sigma <- sd(x)
	mean(((x-mu)/sigma)^3)
	}


# mean, sd and skew tha result from adding variables of different distributions
add_skew <- function(means, sds, skews){
	new_mean <- sum(means)
	new_sd <- sqrt(sum(sds^2))
	new_skew <- sum(skews*sds^3)/new_sd^3
	return(c(new_mean,new_sd,new_skew))
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

factorisePed <- function(pedigree, unknown=0){
    new_ped <- data.frame(
        1:nrow(pedigree), 
        ifelse(is.na(pedigree[,2]),unknown,match(pedigree[,2], pedigree[,1])), 
        ifelse(is.na(pedigree[,3]),unknown,match(pedigree[,3], pedigree[,1]))
        )
    colnames(new_ped) <- colnames(pedigree)[1:3]

    return(new_ped)
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



mean_CI <- function(x) c(mean(x),quantile(x, c(0.025,0.975)))

