# Extra functions for data visualization and nonparametric analysis, mostly based on `rogme` functions
# Zack Williams
# Started 06/07/2018
# Last edit 04/18/2021

# Package loading using pacman
if(!require(pacman)){
  install.packages("pacman")
  require(pacman)
}
if(!require(rogme)){
  if(!require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("GRousselet/rogme")
}
if(!require(orddom)){
  devtools::install_github("cran/orddom")
}
# Required packages: 
pacman::p_load(rogme,orddom,tidyverse,easystats,psych,optimParallel,orddom,parallel,norm,mipfp,
               lawstat,cocor,pROC,Partiallyoverlapping,patchwork,qgraph,eva,ggokabeito,ggthemes)

# Utility functions: rep.row, rep.col
rep.row<-function(x,n){ 
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# %!in% - opposite of the %in% operator (credit to Brenton Wiernik); equivalent to Hmisc::%nin%
`%!in%` <- Negate(`%in%`)

# Skew-Normal Statistics
skewnorm <- function(x,na.rm=TRUE,type=3){
  mx <- mean(x)
  sdx <- sd(x)
  skx <- psych::skew(x,na.rm=na.rm,type=type)
  delx <- sqrt(pi/2 * skx^(2/3)/(skx^(2/3) + ((4-pi)/2)^(2/3)))
  
  scalex <- sdx/sqrt(1-2/pi * delx^2)
  locx <- mx - sqrt(2/pi)*delx*scalex
  return(c("Mean"=mx,"SD"=sdx,"Skew"=skx,"Location"=locx,"Scale"=scalex))
}

#  Compute the Harrell-Davis estimate of the qth quantile (q can now be a vector)
#  The vector x contains the data,
#  and the desired quantile is q
#  The default value for q is .5.
hd <- function(x,qs=.5,na.rm=TRUE){

  elimna<-function(m){
    #
    # remove any rows of data having missing values
    #
    DONE=FALSE
    if(is.list(m) && is.matrix(m)){
      z=pool.a.list(m)
      m=matrix(z,ncol=ncol(m))
      DONE=TRUE
    }
    if(!DONE){
      if(is.list(m) && is.matrix(m[[1]])){
        for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
        e=m
        DONE=TRUE
      }}
    if(!DONE){
      if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
        for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
        e=m
        DONE=TRUE
      }}
    if(!DONE){
      #if(!is.list(m)){
      #if(is.null(dim(m)))
      m<-as.matrix(m)
      ikeep<-c(1:nrow(m))
      for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
      e<-m[ikeep[ikeep>=1],]
      #}
    }
    e
  }
  if(na.rm){x=elimna(x)}
  
  n<-length(x)
  hd <- NULL
  for(q in qs){
    m1<-(n+1)*q
    m2<-(n+1)*(1-q)
    vec<-seq(along=x)
    w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
    y<-sort(x)
    hd<-c(hd,sum(w*y))
  }
  hd
}

p.transform <- function(p,PProb.H0=NULL,digits=3){
  # Assuming 1:1 prior odds
  BFB <- 1/(-exp(1)*p*log(p))
  pH1 <- BFB/(1+BFB)
  # Incorporating prior odds
  if(!is.null(PProb.H0)){
    PProb.H1 <- 1 - PProb.H0
    prior.odds <- PProb.H1/PProb.H0
    BFB.p <- 1/(-exp(1)*p*log(p))*prior.odds
    pH1.p <- BFB.p/(1+BFB.p)
  } else{
    BFB.p <- pH1.p <- NULL
    prior.odds <- 0.5
  }
  
  # S-value (# coin flips all coming up heads = p value)
  S=-log2(p)
  
  # Intrinsic Credibility P-value (Held, 2019: https://doi.org/10.1098/rsos.181534)
  p.IC <- 2*(1 - pnorm(qnorm(1 - p/2)/sqrt(2)))
  p.rep <- 1 - p.IC/2
  
  return(round(c("BFB.equal"=BFB,"P.H1.equal"=pH1,
                 "BFB.prior"=BFB.p,"P.H1.equal"=pH1.p,
                 "Prior(H1/H0)"=prior.odds,
                 "S.value"=S,"p.IC"=p.IC,"p.rep"=p.rep),digits=digits))
}

# Mean and SD
MSD <- function(x,digits=2){
  return(round(c("M"=mean(x,na.rm=T),"SD"=sd(x,na.rm = T)),digits))
}

# empty.rows
# Goes through the rows of a df and determines which are empty
empty.rows <- function(df){
  return(apply(df,1,function(X){all(is.na(X))}))
}


# Spearman Correlation test
Spearman.test <- function(x,y,digits=3,alpha=0.05,H0=c(-0.1,0.1),s2="cc",pandeR=F){
  require(Hmisc)
  n <- sum(complete.cases(cbind(x,y)))
  rho <- cor(x,y,use="pairwise",method="spearman")
  switch(s2,
         f = s2rho <- 1.06/(n-3), # Fieller, Hartley, & Pearson, 1957; Biometrika
         bw = s2rho <- (1 + rho^2/2)/(n-3), # Bonett & Wright, 2000; Psychometrika
         cc = s2rho <- 1/(n-2) + abs(atanh(rho))/(6*n + 4*sqrt(n)) # Caruso & Cliff, 1997; EPM
         )
  
  delta <- qnorm(1-alpha)*sqrt(s2rho)
  rho.low <- tanh(atanh(rho) - delta)
  rho.high <- tanh(atanh(rho) + delta)
  
  # Z and p levels
  prho <- c(2*pnorm(-abs(atanh(rho)/sqrt(s2rho))), # p.null
            pnorm((atanh(H0[1]) - atanh(rho))/sqrt(s2rho)), # p.lower bound
            pnorm((atanh(rho) - atanh(H0[2]))/sqrt(s2rho))) # p.higher bound
  
  p.equiv <- max(prho[2:3])
  
  result <- signif(c(rho,rho.low,rho.high,prho[1],p.equiv),digits)
  names(result) <- c("rho.s",paste0("rho.ci",100*(1-2*alpha),".low"),paste0("rho.ci",100*(1-2*alpha),".high"),"p.val","p.equiv")
  
  if(pandeR==T) result <- pander(result)
  
  return(result)
}

#### Signal Detection Metrics: classic and non-parametric
### Zack Williams, 10/28/19
# based on Stanislaw, H., & Todorov, N. (1999). 
# Calculation of signal detection theory measures. 
# Behavior research methods, instruments, & computers, 31(1), 137-149.
#### Includes A and b from Zhang & Mueller (2005)
## Take a hit rate and false alarm rate and calculate these indices
## http://doi.org/10.1007/s11336-003-1119-8
## Code adapted from Shane Mueller (https://sites.google.com/a/mtu.edu/whynotaprime/)
## Note that b = 1 is unbiased. log(b) is a symmetric bias measure (< 0 = more liberal, > 0 = more conservative)
# When calculating d', implements the log-linear correction (all Hits/Htrials/FA/FAtrials corrected)
SDTmeasures <-function(Hits,FAs,Htrials,FAtrials,h=NULL,f=NULL,digits=3,loglinear.correct=F){
  if(all(is.null(c(h,f)))){
    h <- Hits/Htrials
    f <- FAs/FAtrials
  }
  if(!is.null(Htrials) & !is.null(FAtrials)){
    Pc <- (Hits + FAtrials - FAs)/(Htrials + FAtrials)
  } else{
    Pc <- NULL
  }
  # Classic measures: calculate d' (sensitivity), beta (bias), and c (criterion)
  # For c, negative values = more liberal, positive = more conservative
  if(h==1 | loglinear.correct){ # Correction in d' for no misses 
    if(f==0 | loglinear.correct){ # Correction for no FAs
      dprime <- qnorm((Hits+0.5)/(Htrials+1))-qnorm((FAs+0.5)/(FAtrials+1))
      beta <- exp((qnorm((FAs+0.5)/(FAtrials+1))^2 - qnorm((Hits+0.5)/(Htrials+1))^2)/2)
      c <- -(qnorm((Hits+0.5)/(Htrials+1))+qnorm((FAs+0.5)/(FAtrials+1)))/2
    } else{ 
      dprime <- qnorm((Hits+0.5)/(Htrials+1))-qnorm(FAs)
      beta <- exp((qnorm(f)^2 - qnorm((Hits+0.5)/(Htrials+1))^2)/2)
      c <- -(qnorm((Hits+0.5)/(Htrials+1))+qnorm(f))/2
    }
  } else{ 
    if(f==0){ # correction for no FAs
      dprime <- qnorm(h)-qnorm((FAs+0.5)/(FAtrials+1))
      beta <- exp((qnorm((FAs+0.5)/(FAtrials+1))^2 - qnorm(h)^2)/2)
      c <- -(qnorm(h)+qnorm((FAs+0.5)/(FAtrials+1)))/2
    } else{
      dprime <- qnorm(h)-qnorm(f)
      beta <- exp((qnorm(f)^2-qnorm(h)^2)/2)
      c <- -(qnorm(h)+qnorm(f))/2
    }
  }
  # Calculation of Az, which is a transformation of d'
  Az <- pnorm(dprime/(sqrt(2)))
  # Calculation of c' (bias as a proportion of sensitivity [c/d'])
  cprime <- c/dprime
  # Calculate A' and B''
  Aprime <- 0.5 + sign(h - f)*((h - f)^2 + abs(h - f))/(4*max(h,f)-4*h*f)
  Bprime <- sign(h-f)*(h*(1-h) - f*(1-f))/(h*(1-h) + f*(1-f))
  # A calculation (Zhang & Mueller, 2005)
  if(f<=.5 & h>=.5){
    a <- .75 + (h-f)/4 - f*(1-h)
    b <-(5-4*h)/(1+4*f)
  } else if(f<=h & h<=.5){
    a <- .75 + (h-f)/4 - f/(4*h)
    b <-(h^2+h)/(h^2+f)
  } else {
    a <- .75 + (h-f)/4 - (1-h)/(4 * (1-f))
    b <- ((1-f)^2 + (1-h))/((1-f)^2 + (1-f))
  }
  return(round(c("Pc"=Pc,"d'"=dprime,"Az"=Az,"Beta"=beta,"log(Beta)"=log(beta),"c"=c,"c'"=cprime,
                 "A'"=Aprime,"B''"=Bprime,"exp(B'')"=exp(Bprime),
                 "A"=a,"b"=b,"log(b)"=log(b)),digits))
}

# calcSMD - calculates Cohen's d and confidence interval
# Can take means and variances, as well as raw data
# Zack Williams, 11/27/2019
calcSMD <- function(M1=0,M2=0,V1=1,V2=1,N1=NULL,N2=NULL,x=NULL,y=NULL,ci.width=0.95,digits=3){
  # if x is given
  if(!is.null(x)){
    M1 <- mean(x,na.rm=T)
    V1 <- var(x,na.rm=T)
    N1 <- sum(!is.na(x) & !is.nan(x),na.rm=T)
  }
  # if y is given
  if(!is.null(y)){
    M2 <- mean(y,na.rm=T)
    V2 <- var(y,na.rm=T)
    N2 <- sum(!is.na(y) & !is.nan(y),na.rm=T)
  }
  
  alpha <- 1 - ci.width
  pooledSD <- sqrt((V1*(N1-1) + V2*(N2-1))/(N1 + N2 - 2))
  d <- (M2 - M1)/pooledSD
  var.d <- ((N1 + N2)/(N1*N2) + d^2/(2*(N1 + N2 - 2))) * (N1 + N2)/(N1 + N2 - 2)
  se.d <- sqrt(var.d)
  ci.lo <- qnorm(alpha/2,mean=d,sd=se.d)
  ci.hi <- qnorm(1-alpha/2,mean=d,sd=se.d)
  return(round(c("d"=d,"d.ci.lo"=ci.lo,"d.ci.hi"=ci.hi,"CI"=ci.width,"n1"=N1,"n2"=N2),digits))
}

# var.test.partover
# Tests the equality of two groups that have at least some dependence 
# Can work with paired samples with any number of independent samples
# See Derrick et al. (2018): Tests for Equality of Variances between 
# Two Samples which Contain Both Paired Observations and Independent Observations
# http://www.jaqm.ro/issues/volume-13,issue-2/3_BEANDEPA.PHP
# x1 and x2 are two samples of values, including both paired and unpaired observations
var.test.partover <- function(x1 = NULL, x2 = NULL,var.equal = FALSE, mu = 0, 
                             alternative = "two.sided"){
  X1 <- abs(x1 - median(x1))
  X2 <- abs(x2 - median(x2))
  return(Partiallyoverlapping::Partover.test(x1=X1,x2=X2,var.equal=var.equal,
          mu=mu,alternative=alternative,stacked = TRUE))
}

# Second-generation P value calculation
# Zack Williams
# Updated 3/31/2019
# See also `sgpv` package (https://github.com/weltybiostat/sgpv)
# Blume, J. D., Greevy, R. A., Welty, V. F., Smith, J. R., & Dupont, W. D. (2019).
# An Introduction to Second-Generation p-Values. Am Stat, 73(sup1), 157-167.
# https://doi.org/10.1080/00031305.2018.1537893
# Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). 
# Second-generation p-values: Improved rigor, reproducibility, & transparency in statistical analyses. 
# PLoS ONE 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299
pval.2g <- function(CI, H0_l, H0_u,digits=3,delta.gap=F){
  if(length(CI) != 2) stop("ERROR: 'CI' must contain two values")
  if(CI[1] > CI[2]) stop("ERROR: 'CI' must include the lower bound first")
  if(H0_l >= H0_u) stop("ERROR: 'H0_l' must be less than 'H0_u'")
  CI_l <- CI[1]
  CI_u <- CI[2]
  I <- diff(CI)
  H0 <- diff(c(H0_l, H0_u))
  # Calculate overlap of two CIs
  intersect <- max(diff(c(min(c(max(c(H0_l,CI_l)),H0_u)),min(c(H0_u,CI_u)))),0)
  if(I==0){ # CI bounds are the same
    if(CI_l > H0_l & CI_l < H0_u){
      p.2g <- 1 # Return 1 if point CI within H0 bounds
    } else{
      p.2g <- 0 # Return 0 if point CI out of H0 bounds
    }
  } else{
    p.2g <- signif(intersect/I*max(I/(2*H0),1),digits)
  }
  # Added in delta.gap, a proposed effect size metric for p.2g
  # Delta gap is the number of H0 interval lengths between the upper part of H0 and the lower part of CI
  if(delta.gap==T){
    delta <- H0/2
    if(H0_u < CI_l){ # CI is above interval null
      d.gap <- (CI_l - H0_u)/delta
    } else if(CI_u < H0_l){ # CI is below interval null
      d.gap <- (CI_u - H0_l)/delta
    } else{
      d.gap <- 0
    }
    return(c("pval.2g"=p.2g,"delta.gap"=d.gap))
  } else{
    return(p.2g)
  }
}

# Outlier detection using Tukey's G and H distribution, from below paper:
# Xu, Y., Iglewicz, B., & Chervoneva, I. (2014). 
# Robust estimation of the parameters of g-and-h distributions, with applications
# to outlier detection. Computational statistics & data analysis, 75, 66-80.
gh_outliers_robust <- function(x,m=10,alpha=0.05,skewed=NULL){
  x <- na.omit(x)
  # Random starting values
  theta0 <- c(mean(x),rnorm(1),rnorm(1),0.5)
  # 10 quantile method. i = 1:10, m = 10; quantiles are (i-1/3)/(m+1/3)
  i <- 1:m
  q <- (i-1/3)/(m+1/3)
  
  # Function to calculate the QLS quantile discrepancy function for the g-h dist
  qls_gh <- function(pars,rawdat,q){
    emp_q <- hd(rawdat,q)
    est_q <- qgh(q,pars[1],exp(pars[2]),pars[3],pars[4])
    return(sum((emp_q-est_q)^2))
  }
  
  biweight_qs <- function(pars,rawdat,q,c=5){
    emp_q <- hd(rawdat,q)
    est_q <- qgh(q,pars[1],exp(pars[2]),pars[3],pars[4])
    outliers <- abs(emp_q - est_q) > c
    wts <- 1-((emp_q - est_q)/c)^2
    wts[outliers] <- 0
    return(list(q=wts*emp_q,i=which(!outliers),w=wts))
  }
  
  # Adaptive calculation of the c parameter for biweight
  get_c <- function(pars,rawdat,q){
    # C is best found in the range (a,b/2), with a and b defined as below
    # a is the min value of c where >50% data retained
    # b is the min value of c where 100% data retained
    get_ab <- function(pars,rawdat,q){
      emp_q <- hd(rawdat,q)
      est_q <- qgh(q,pars[1],exp(pars[2]),pars[3],pars[4])
      diffs <- abs(emp_q - est_q)
      return(c(median(diffs),max(diffs)))
    }
    ab <- get_ab(pars,rawdat,q)
    v <- 10^(log(ab[2]/2,10)-1) # step size
    if(ab[2]/2 < ab[1]) return(ab[1])
    c <- c(seq(ab[2]/2,ab[1],-v),ab[1])
    # Move through the weights looking for stopping criteria
    for(ci in c){
      bwi <- biweight_qs(pars,rawdat,q,c=ci)
      wts <- bwi$w[bwi$i]
      if(all(wts > 0.8)){
        return(ci)
      }
      n <- length(wts)
      s1 <- max(c(which(wts[seq(1,n/2,1)]<0.7),0))
      s2 <- min(c(which(wts[seq(1,n/2,1)]>0.8),n/2))
      s3 <- max(c(which(wts[seq(ceiling(n/2),n,1)]>0.8),0))
      s4 <- min(c(which(wts[seq(ceiling(n/2),n,1)]<0.7),n))
      if(s1 < s2 & which.min(wts[seq(1,n/2,1)]) == 1 & 
         s3 < s4 & which.min(wts[seq(1,n/2,1)]) == length(wts[seq(1,n/2,1)])){
        return(ci)
      }
    }
    return(ci)
  }
  
  # Function to calculate the QLS quantile discrepancy function for the g-h dist
  rqls_gh <- function(pars,rawdat,q,c=NULL){
    if(is.null(c)){
      c <- get_c(pars,rawdat,q)
    }
    bw <- biweight_qs(thetaQLS,rawdat,q,c=c)
    emp_q <- bw$q[bw$i]
    est_q <- qgh(q[bw$i],pars[1],exp(pars[2]),pars[3],pars[4])
    return(sum((emp_q-est_q)^2))
  }

  # Ordinary QLS estimation first
  thetaQLS <- suppressWarnings(optimParallel(theta0,qls_gh,rawdat=x,q=q,
              parallel = list(cluster=makeCluster(detectCores()-1,type="FORK")))$par)
  # Robust QLS estimation using ordinary QLS parameter estimates
  # Repeat until convergence
  # try(stopCluster(cl))
  cl <- makeCluster(detectCores()-1,type="FORK"); setDefaultCluster(cl = cl)
  theta_converge <- function(theta,tol=0.01,max.iter=50){
    iter=1
    theta_new <- suppressWarnings(optimParallel(theta,qls_gh,rawdat=x,q=q,
                 parallel = list(cluster=makeCluster(detectCores()-1,type="FORK")))$par)
    while(mean(abs(theta_new - theta)) > tol & iter < max.iter){
      theta <- theta_new
      theta_new <- suppressWarnings(optimParallel(theta,qls_gh,rawdat=x,q=q,
                   parallel = list(cluster=makeCluster(detectCores()-1,type="FORK")))$par)
      iter <- iter + 1
    }
    return(theta_new)
  }
  thetaQLSr <- theta_converge(thetaQLS) # Run iterative algorithm
  
  names(thetaQLSr) <- c("A","log(B)","g","h")
  
  gh_q <- function(q){
    pgh(q, thetaQLSr[1], exp(thetaQLSr[2]), thetaQLSr[3], thetaQLSr[4])
  }
  # Estimate K for boxplot outliers
  if(is.null(skewed)){
    sk <- psych::skew(x)
    if(sk >= 1) skewed = "Right"
    else if(sk <= -1) skewed = "Left"
    else skewed = "Symmetric"
  }
  if(skewed %in% c("symmetric","none","Symmetric","None")){
    k <- diff(gh_q(c(0.75,(1-alpha/2)^(1/length(x)))))/diff(gh_q(c(0.25,0.75)))
    outliers <- x > hd(x,0.75) + k*diff(hd(x,c(0.25,0.75))) | x < hd(x,0.25) - k*diff(hd(x,c(0.25,0.75)))
  }else if(skewed %in% c("right","Right")){
    k <- diff(gh_q(c(0.75,(1-alpha)^(1/length(x)))))/diff(gh_q(c(0.5,0.75)))
    outliers <- x > hd(x,0.75) + k*diff(hd(x,c(0.5,0.75)))
  }else if(skewed %in% c("left","Left")){
    k <- diff(gh_q(c(0.25,1-(1-alpha)^(1/length(x)))))/diff(gh_q(c(0.5,0.25)))
    outliers <- x < hd(x,0.25) - k*diff(hd(x,c(0.25,0.5)))
  }
  
  # Return MOM, GH-identified outliers, their indices, and their values
  MOM <- mean(x[!outliers])
  result <- list(nOutliers=sum(outliers),Indices=which(outliers),Vals=x[outliers],
                 MOM=MOM,Outlier.mean=mean(x[outliers]),"GH.pars"=thetaQLSr,k=unname(k),
                 skewed=skewed)
  return(result)
}

# Faster gh outliers function without all the robust stuff
gh_outliers <- function(x,m=10,alpha=0.05,skewed=NULL){
  x <- na.omit(x)
  # Random starting values
  theta0 <- c(mean(x),rnorm(1),rnorm(1),0.25)
  # 10 quantile method. i = 1:10, m = 10; quantiles are (i-1/3)/(m+1/3)
  i <- 1:m
  q <- (i-1/3)/(m+1/3)
  
  # Function to calculate the QLS quantile discrepancy function for the g-h dist
  qls_gh <- function(pars,rawdat,q){
    emp_q <- hd(rawdat,q)
    est_q <- qgh(q,pars[1],exp(pars[2]),pars[3],pars[4])
    return(sum((emp_q-est_q)^2))
  }
  
  # Parallel processing for faster optimization
  # Ordinary QLS estimation first
  thetaQLSr <- optimParallel(theta0,qls_gh,rawdat=x,q=q,
               parallel = list(cluster=makeCluster(detectCores()-1,type="FORK")))$par
  
  names(thetaQLSr) <- c("A","log(B)","g","h")
  
  gh_q <- function(q){
    pgh(q, thetaQLSr[1], exp(thetaQLSr[2]), thetaQLSr[3], thetaQLSr[4])
  }
  # Estimate K for boxplot outliers
  if(is.null(skewed)){
    sk <- psych::skew(x)
    if(sk >= 1) skewed = "Right"
    else if(sk <= -1) skewed = "Left"
    else skewed = "Symmetric"
  }
  if(skewed %in% c("symmetric","none","Symmetric","None")){
    k <- diff(gh_q(c(0.75,(1-alpha/2)^(1/length(x)))))/diff(gh_q(c(0.25,0.75)))
    outliers <- x > hd(x,0.75) + k*diff(hd(x,c(0.25,0.75))) | x < hd(x,0.25) - k*diff(hd(x,c(0.25,0.75)))
  }else if(skewed %in% c("right","Right")){
    k <- diff(gh_q(c(0.75,(1-alpha)^(1/length(x)))))/diff(gh_q(c(0.5,0.75)))
    outliers <- x > hd(x,0.75) + k*diff(hd(x,c(0.5,0.75)))
  }else if(skewed %in% c("left","Left")){
    k <- diff(gh_q(c(0.25,1-(1-alpha)^(1/length(x)))))/diff(gh_q(c(0.5,0.25)))
    outliers <- x < hd(x,0.25) - k*diff(hd(x,c(0.25,0.5)))
  }
  
  # Return MOM, GH-identified outliers, their indices, and their values
  MOM <- mean(x[!outliers])
  result <- list(nOutliers=sum(outliers),Indices=which(outliers),Vals=x[outliers],
                 MOM=MOM,Outlier.mean=mean(x[outliers]),"GH.pars"=thetaQLSr,k=unname(k),
                 skewed=skewed)
  return(result)
}
# Outliers using the boxplot rule (1.5 * IQR away from Q1/Q3)
boxplot_outliers <- function(x,na.rm=TRUE,coef=1.5){
  if(na.rm)x<-x[!is.na(x)] #Remove missing values
  q13 <- hd(x,c(0.25,0.75)) # Quartiles using HD values
  iqr <- diff(q13)
  is_outlier <- x < q13[1]-coef*iqr | x > q13[2]+coef*iqr # Q1 - 1.5IQR or Q3 + 1.5*IQR
  MOM <- mean(x[!is_outlier])
  result <- list(nOutliers=sum(is_outlier),Indices=which(is_outlier),Vals=x[is_outlier],MOM=MOM,Outlier.mean=mean(x[is_outlier]))
  return(result)
}

# Coefficient of Quartile Variation (Q3 - Q1)/(Q3 + Q1)
cqv <- function(x){
  q13 <- hd(x,c(0.25,0.75))
  return(diff(q13)/sum(q13))
}
# Modified One-step Estimator (MOM) from Rand Wilcox's WRS2 Package
# Also outlier detection based on this
# Zack Williams
# 06/07/2018
mad_outliers <- function(x,bend=2.24,na.rm=TRUE){
  if(na.rm)x<-x[!is.na(x)] #Remove missing values
  is_outlier <- abs((x-hd(x)))/(1.4826*hd(abs(x-hd(x)))) > bend
  MOM <- mean(x[!is_outlier])
  result <- list(nOutliers=sum(is_outlier),Indices=which(is_outlier),Vals=x[is_outlier],MOM=MOM,Outlier.mean=mean(x[is_outlier]))
  return(result)
}

mom<-function(x,bend=2.24,na.rm=TRUE){
  #
  #  Compute MOM-estimator of location.
  #  The default bending constant is 2.24
  #
  if(na.rm)x<-x[!is.na(x)] #Remove missing values
  flag1<-(x>median(x)+bend*mad(x))
  flag2<-(x<median(x)-bend*mad(x))
  flag<-rep(T,length(x))
  flag[flag1]<-F
  flag[flag2]<-F
  mom<-mean(x[flag])
  mom
}

# Summary function: Mean, SD, Median (HD), IQR (HD), MAD, Min, Max
summ <- function(x,name=NULL,hdQ=F,NAs=F,dig=2){
  if(hdQ==TRUE){
    result <- c("M" = mean(x,na.rm=T),
                "SD" = sd(x,na.rm=T),
                "Mdn" = hd(x),
                "Q1" = hd(x,0.25),
                "Q3" = hd(x,0.75),
                "IQR" = diff(hd(x,c(0.25,0.75))),
                "MAD" = 1.4826*hd(abs(x-hd(x))),
                "Min" = min(x,na.rm=T),
                "Max" = max(x,na.rm=T))
  } else{
    result <- c("M" = mean(x,na.rm=T),
                "SD" = sd(x,na.rm=T),
                "Mdn" = median(x,na.rm = T),
                "Q1" = unname(quantile(x,0.25,na.rm = T)),
                "Q3" = unname(quantile(x,0.75,na.rm = T)),
                "IQR" = IQR(x,na.rm = T),
                "MAD" = 1.4826*hd(abs(x-hd(x))),
                "Min" = min(x,na.rm=T),
                "Max" = max(x,na.rm=T))
  }
  if(NAs){
    result <- c(result,"NAs"=sum(is.na(x)))
  }
  
  if(!is.null(name)){
    names(result) <- paste(names(result),name,sep=".")
  }
  
  return(round(result,digits=dig))
}

# Median Absolute Deviation using the Harrell-Davis Quantile Estimator
madhd <- function(x,constant=1.4826,na.rm=T){
  if(na.rm)x<-x[!is.na(x)] 
  return(constant*hd(abs(x-hd(x))))
}

# Kernel Density Plots 
# Zack Williams
# 06/07/2018
# Using code from `rogme` package that I debugged
# 1 group
kdeplot.1g <- function(data=df,fill.colour="grey30",fill.alpha=.3){
  data <- mkt1(data)
  cdat <- plyr::ddply(data, "gr", plyr::summarise, deciles=hd(obs,c(0.1,0.4,0.6,0.9)))
  hd05 <- plyr::ddply(data, "gr", plyr::summarise, hd=hd(obs,0.5))
  cc <- "grey80" # colour to plot deciles
  p <- ggplot(data, aes(x=obs)) +
    geom_density(alpha=fill.alpha,fill=fill.colour,colour="black") +
    geom_vline(xintercept=hd05$hd, colour="black", linetype="solid",
               size=2, alpha=0.5) + # thicker median
    geom_vline(data=cdat, aes(xintercept=deciles),
               linetype="solid", size=1, alpha=.5, colour="black") +
    geom_rug() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold")) +
    ylab("Density")
  p
}

# Kernel Density Plots 
# Zack Williams
# 06/07/2018
# Using code from `rogme` package that I debugged
# 2 groups
kdeplot.2g <- function(x,y=NULL,data=NULL,labs=NULL,dvlab="obs",denslab="Density"){
  
  if(is.null(y) & class(x)=="formula"){
    # Function to make tibble from formula
    makeTibble <- function(formula,data){
      if(missing(formula)
         || (length(formula) != 3L)
         || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
      m <- match.call(expand.dots = FALSE)
      if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
      ## need stats:: for non-standard evaluation
      m[[1L]] <- quote(stats::model.frame)
      mf <- eval(m, parent.frame())
      names(mf) <- NULL
      response <- attr(attr(mf, "terms"), "response")
      g <- factor(mf[[-response]])
      if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
      DATA <- setNames(split(mf[[response]], g), c("x", "y"))
      y <- do.call("mkt2", DATA)
      levels(y$gr) <- levels(g)
      y
    }
    tib <- makeTibble(x,data=data)
  } else{
    tib <- mkt2(c(na.omit(x)),c(na.omit(y)))
  }
  
  if(!is.null(labs)){levels(tib$gr) <- labs} # Replace group labels
  
  cdat <- plyr::ddply(tib, "gr", plyr::summarise, qs = hd(obs,c(0.1,0.25,0.75,0.9)))
  hd05 <- plyr::ddply(tib, "gr", plyr::summarise, hd = hd(obs,0.5))
  #cc <- "grey80" # colour to plot deciles
  p <- ggplot(tib, aes(x=obs, fill=gr)) + geom_density(alpha=.3) +
    facet_grid(gr ~ .) +
    geom_vline(data=hd05, aes(xintercept=hd,  colour=gr),
               linetype="solid", size=2, alpha=.5) + # thicker median
    geom_vline(data=cdat, aes(xintercept=qs,  colour=gr),
               linetype="solid", size=1, alpha=.5) +
    geom_rug() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold"),
          strip.text.y = element_text(size = 20, colour = "white"),
          strip.background = element_rect(colour="darkgrey", fill="darkgrey")) +
    ylab(eval(denslab)) + xlab(eval(dvlab))
  p
}

# 2-group scatterplot with link fn (with modifications to make it like a boxplot)
quantile_scatterplot <- function(x,y=NULL,data=NULL,labs=NULL,q = c(.1,.25,.5,.75,.9),shift=TRUE,
                        nboot=1000,flip=FALSE,dvlab="Scores (a.u.)",outliers="Boxplot",stats=T,
                        returnPlot=T){
  if(is.null(y) & class(x)=="formula"){
    # Function to make tibble from formula
    makeTibble <- function(formula,data){
      if(missing(formula)
         || (length(formula) != 3L)
         || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
      m <- match.call(expand.dots = FALSE)
      if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
      ## need stats:: for non-standard evaluation
      m[[1L]] <- quote(stats::model.frame)
      mf <- eval(m, parent.frame())
      names(mf) <- NULL
      response <- attr(attr(mf, "terms"), "response")
      g <- factor(mf[[-response]])
      if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
      DATA <- setNames(split(mf[[response]], g), c("x", "y"))
      y <- do.call("mkt2", DATA)
      levels(y$gr) <- levels(g)
      return(y)
    }
    tib <- makeTibble(x,data=data)
  } else{
    tib <- mkt2(c(na.omit(x)),c(na.omit(y))) 
  }
  if(!is.null(labs)){levels(tib$gr) <- labs} # Replace group labels
  
  # Outlier detection for point re-coloring
  grps <- levels(tib$gr)
  g1 <- unlist(tib[which(tib$gr==grps[1]),2])
  g2 <- unlist(tib[which(tib$gr==grps[2]),2])
  ng1 <- length(g1)
  ng2 <- length(g2)
  if(outliers %in% c("MAD","mad")){
    ol <- c(mad_outliers(g1)$Indices,ng1+mad_outliers(g2)$Indices)
  } else if(outliers %in% c("boxplot","Boxplot")){
    ol <- c(boxplot_outliers(g1)$Indices,ng1+boxplot_outliers(g2)$Indices)
  } else if(outliers %in% c("gh","GH")){
    ol <- c(gh_outliers(g1)$Indices,ng1+gh_outliers(g2)$Indices)
  } else if(outliers %in% c("ghr","GHR","robust","Robust")){
    ol <- c(gh_outliers_robust(g1)$Indices,ng1+gh_outliers_robust(g2)$Indices)
  } else{
    ol <- NULL
  }
  pointcols <- c(rep("#66A61E",ng1),rep("#E7298A",ng2))
  pointcols[ol] <- "white" # Color outlier points with white
  
  # Shift function
  sf <- shifthd_pbci(data = tib, formula = obs ~ gr, nboot = nboot, q = q)
  # Scatter plot by group
  p <- plot_scat2(tib,
                  xlabel = "",
                  ylabel = dvlab,
                  alpha = .7,
                  shape = 21,
                  colour = "grey10",
                  fill = pointcols)    
  if(shift==TRUE){
    p <- plot_hd_links(p, sf[[1]],
                       q_size = 0.75,
                       md_size = 1.5,
                       # add_rect = TRUE,
                       # rect_alpha = 0.1,
                       rect_col = "grey50",
                       add_lab = TRUE,
                       text_size = 5) # superimposed deciles + rectangle 
    
    p <- p + annotate("label", # Label for median difference
                      x = 1.5,
                      y = min(sf[[1]][3,2],sf[[1]][3,3]) + abs(sf[[1]][3,2] - sf[[1]][3,3]) / 2,
                      label = round(sf[[1]][3,4],2),
                      fill = c("darkviolet","darkorange2")[sign(sf[[1]][3,4]>0)+1],
                      colour = "white",
                      fontface = "bold",
                      alpha = 1)
  }
  
  p <- p + annotate("rect", xmin = 0.75, xmax = 1.25, ymin = sf[[1]][2,2], 
               ymax = sf[[1]][4,2], alpha = 0.1,color="grey21",size=0.75) + 
           annotate("rect", xmin = 1.75, xmax = 2.25, ymin = sf[[1]][2,3], 
               ymax = sf[[1]][4,3], alpha = 0.1,color="grey21",size=0.75) +
           # Vertical boxplot lines
           annotate("segment",x = 1, xend = 1, y = sf[[1]][1,2], yend = sf[[1]][2,2],
                    color="grey21",size=0.75)  +
           annotate("segment",x = 1, xend = 1, y = sf[[1]][4,2], yend = sf[[1]][5,2],
                    color="grey21",size=0.75)  +
           annotate("segment",x = 2, xend = 2, y = sf[[1]][1,3], yend = sf[[1]][2,3],
                    color="grey21",size=0.75)  +
           annotate("segment",x = 2, xend = 2, y = sf[[1]][4,3], yend = sf[[1]][5,3],
                    color="grey21",size=0.75)
  
  if(flip==TRUE){
    p <- p + coord_flip() # flip axes  
  }
  if(returnPlot==T){
    return(p)
  }else{
    print(p)
  }
  if(stats==T & shift==T){return(round(sf[[1]],4))}
}

# Paired Differences plot w/ lines
pairdiff_plot <- function(df,groupvar,outcome,subjectvar,
                          xlabs=c("Condition 1","Condition 2"),ylab="Scores (a.u.)"){
  ggplot(df, aes(x=df[,groupvar], y=df[,outcome], group=df[,subjectvar])) +
    geom_line(aes(colour=df[,subjectvar]),size=1, alpha=.5,
              position=position_dodge(width = 0.1)) +
    geom_point(aes(colour=df[,subjectvar]),position=position_dodge(width = 0.1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16,face="bold"),
          legend.position="none") +
    labs(title="Paired observations") +
    ylab(label=ylab) +
    scale_x_discrete(labels = c("condition1" = xlabs[1],"condition2" = xlabs[2])) 
  #scale_y_continuous(limits=c(0, 1.1*max(df[,outcome])))#,breaks=seq(0,30,5))
}

# Conversion between Cohen's U1 (overlap) and cohen's d - used as delta/CohD conversion in `orddom`
u12cohd <-function(d){ #delta (Overlap) to Cohen's d
  # sign(d)*qnorm(-1/(abs(d)-2))*2
  return (sign(d)*qnorm(-1/(abs(d)-2))*2)
}
cohd2delta <- function(d){ #U1 in Table 2.2.1 in Cohen's Statistical Power Analysis, p. 22
  return(sign(d)*((2*pnorm(abs(d)/2))-1)/pnorm(abs(d)/2))
}

# Conversion between cliff delta and cohen's d based on AUC
delta2cohd <-function(d){ #delta (2*AUC-1) to Cohen's d
  # note that this is different from orddom package, which finds u12cohd
  return (qnorm((d + 1)/2)*sqrt(2))
}
cohd2delta <- function(d){ #Cohen's d to Cliff's delta (based on AUC)
  return(2*pnorm(d/sqrt(2))-1) # (2*AUC-1)
}
#Based on AUC to coh D conversion
cohd2vda <- cohd2auc <- function(d){ 
  return(pnorm(d/sqrt(2)))
}
vda2cohd <- auc2cohd <- function(A){ 
  return(qnorm(A)*sqrt(2))
}
# Convert Cliff's delta to AUC/Vargha-Delaney A
delta2vda <- delta2auc <- function(d){return((d+1)/2)}
# Convert AUC/Vargha-Delaney A to Cliff's delta
vda2delta <- auc2delta <- function(A){return(2*A-1)}

# Cliff's delta effect size (with pretty printing)
cliff.d <- function(X, Y, alpha=0.05,paired=F,SGPV=F,H0=c(cohd2delta(-0.2),cohd2delta(0.2)),
                    coh.d=F,VDA = F,...){
  # Remove NAs
  X <- c(na.omit(X))
  Y <- c(na.omit(Y))
  
  if(paired){
    if(length(X) != length(Y)){stop("ERROR: Paired Samples do not have equal length. Try again.")}
    ns <- list(n1 = length(X),n2 = length(Y))
    cliff <- as.list(as.numeric(orddom::orddom(x=X,y=Y,
                                alpha=alpha,paired=TRUE)[c(4:6,11:14,19,29),1],...))
    
    names(cliff) <- c("N(Y>X)","N(Y=X)","N(Y<X)","delta","ci.size","delta.ci.low","delta.ci.high",
                      "p.val","df")
    
    cliff <- c(ns,cliff)
    cliff$"cohen.d" <- delta2cohd(cliff$delta)
    cliff$"d.ci.low" <- delta2cohd(cliff$delta.ci.low)
    cliff$"d.ci.high" <- delta2cohd(cliff$delta.ci.high)
    cliff$n.comp <- cliff$n1
    cliff$"P(Y>X)" <- cliff$"N(Y>X)"/cliff$n.comp
    cliff$"P(Y<X)" <- cliff$"N(Y<X)"/cliff$n.comp
    cliff$"P(Y=X)" <- cliff$"N(Y=X)"/cliff$n.comp
    # Vargha and Delaney's A statistic (equal to a transformation of delta)
    if(VDA==T){
      cliff$VD.A <- (cliff$delta+1)/2
      cliff$A.ci.low <- (cliff$delta.ci.low+1)/2
      cliff$A.ci.high <- (cliff$delta.ci.high+1)/2 
    }
    cliff$NNT <- abs(1/cliff$delta)
  } else{ # if independent samples:
    cliff <- as.list(as.numeric(orddom::orddom(x=X,y=Y,
                                alpha=alpha,paired=FALSE)[c(4:8,13:16,22,29),1],...))
    names(cliff) <- c("n1","n2","N(Y>X)","N(Y=X)","N(Y<X)","delta","ci.size","delta.ci.low","delta.ci.high",
                      "p.val","df")
    if(coh.d==T){
      cliff$"cohen.d" <- delta2cohd(cliff$delta)
      cliff$"d.ci.low" <- delta2cohd(cliff$delta.ci.low)
      cliff$"d.ci.high" <- delta2cohd(cliff$delta.ci.high)      
    }
    cliff$n.comp <- cliff$n1 * cliff$n2
    cliff$"P(Y>X)" <- cliff$"N(Y>X)"/cliff$n.comp
    cliff$"P(Y<X)" <- cliff$"N(Y<X)"/cliff$n.comp
    cliff$"P(Y=X)" <- cliff$"N(Y=X)"/cliff$n.comp
    # Vargha and Delaney's A statistic (equal to a transformation of delta)
    if(VDA==T){
      cliff$VD.A <- (cliff$delta+1)/2
      cliff$A.ci.low <- (cliff$delta.ci.low+1)/2
      cliff$A.ci.high <- (cliff$delta.ci.high+1)/2 
    }
    cliff$NNT <- abs(1/cliff$delta)
  }
  cliff$paired <- paired
  cliff$Call <- match.call()
  class(cliff) <- "cliff.delta"
  cliff$summ.X <- summ(X,"X")
  cliff$summ.Y <- summ(Y,"Y")
  
  # 2nd generation p-value
  if(SGPV==T){
    cliff$p.val.2g <- pval.2g(c(cliff$delta.ci.low,cliff$delta.ci.high),H0[1],H0[2])
    cliff$H0.low <- signif(H0[1],3)
    cliff$H0.high <- signif(H0[2],3)    
  }
  
  return(cliff)
}

# Pretty printing for cliff's delta
print.cliff.delta <- function(x,digits=3,summary=T){
  cat("Call: ")
  print(x$Call)
  cat("\n")
  if(summary == T){
    cat("Nx =", x$n1,"| Ny =",x$n2,"| Total comparisons:",x$n.comp,"\n")
    cat(paste0(" P(Y>X): ", round(x$"P(Y>X)",3)," (N = ",x$"N(Y>X)",")\n"))
    cat(paste0(" P(Y<X): ", round(x$"P(Y<X)",3)," (N = ",x$"N(Y<X)",")\n"))
    cat(paste0(" P(Y=X): ", format(x$"P(Y=X)",digits=3)," (N = ",x$"N(Y=X)",")\n"))
    cat("\n")
    print(x$summ.X)
    print(x$summ.Y)
    cat("\n") 
  }
  if(x$paired==TRUE){cat("Paired ")}
  cat(paste0("Cliff's Delta (Y>X) = ",round(x$delta,digits)," (",x$ci.size,"% CI:",
      paste(format(round(x$delta.ci.low,digits=digits),nsmall=1),
            format(round(x$delta.ci.high,digits=digits),nsmall=1),sep="—"),")\n"))
  cat(paste0(" df = ",x$df,", Pr(>|d|) = ",signif(x$p.val,digits=3)))
  if(x$p.val<0.05){cat("*")}
  if(x$p.val<0.01){cat("*")}
  if(x$p.val<0.005){cat("*")}
  cat("\n")
  print(p.transform(x$p.val,digits = digits)[-3])
  if(!is.null(x$p.val.2g)){
    cat("\n")
    cat(paste0(" 2nd-generation p-value = ",x$p.val.2g," with interval H0 [",x$H0.low,", ",x$H0.high,"]"))
    if(x$p.val.2g==0){cat(" (reject H0)")}
    else if(x$p.val.2g == 0){cat(" (reject |d| > H0)")}
    else {cat(" (inconclusive)")}
  }
  cat("\n\nQualitative Effect Size: ") # From romano et al. (2006)
  # J. Romano, J. D. Kromrey, J. Coraggio, J. Skowronek, 
  # Appropriate statistics for ordinal level data: Should we really be using 
  # t-test and cohen's d for evaluating group differences on the NSSE and other surveys?, 
  # in: Annual meeting of the Florida Association of Institutional Research, 2006.
  if(abs(x$delta) < cohd2delta(0.2)){
    cat("negligible (|d|<0.112)")
  } else if(abs(x$delta) < cohd2delta(0.5)){
    cat("small (|d|<0.276)")
  } else if(abs(x$delta) < cohd2delta(0.8)){
    cat("medium (|d|<0.428)")
  } else{
    cat("large (|d|≥0.429)")
  }
  if(!is.null(x$cohen.d)){
    cat("\nApproximate Cohen's d:",round(x$cohen.d,digits))
    cat(paste0(" (",x$ci.size,"% CI: ",
               paste(format(round(x$d.ci.low,digits=digits),nsmall=1),
                     format(round(x$d.ci.high,digits=digits),nsmall=1),sep="—"),")\n")) 
  }
  if(!is.null(x$p.val.high)){ # equivalence test results
    cat(paste0("\n\nEquivalence tests (TOST) with equivalence bounds [",x$H0.low,", ",x$H0.high,"]"))
    cat(paste0("\n Pr(equivalence) = ",signif(max(x$p.val.low,x$p.val.high),digits=3)))
    if(x$p.val.low < 0.05 & x$p.val.high < 0.05){
      cat("\n Equivalence accepted at p < 0.05 level")
    } else{
      cat("\n Equivalence not accepted (p > 0.05)")
    }
  }
  if(!is.null(x$VD.A)){ # equivalence test results
    cat("\nProbability of stochastic superiority of Y over X: (0<A<1, 0.5 = no dominance)")
    cat("\n Vargha & Delaney's A =",round(x$VD.A,digits))
    cat(paste0(" (",x$ci.size,"% CI: ",
               paste(format(round(x$A.ci.low,digits=digits),nsmall=1),
                     format(round(x$A.ci.high,digits=digits),nsmall=1),sep="—"),")\n"))    
  }
}

# Support functions for Vargha-Delaney A calculations
# Dominance Matrix (from Orddom package)
dm <- function (x,y) { #produces dominance matrix
  dx <- matrix(nrow=length(x), ncol=length(y))
  for (i in 1:length(x)){
    for (j in 1:length(y)){
      dx[i,j]<-sign(y[j]-x[i])
    }
  } 
  rownames(dx) <- x
  colnames(dx) <- y
  return(dx)
}
# Calculation of Vargha & Delaney A
vdA <- function(x,y,ROC=T){
  require(pROC)
  x <- c(na.omit(x))
  y <- c(na.omit(y))
  
  # if not using pROC calculations, uses orddom calculations (required for analytical CI)
  if(!ROC){
    dom <- dm(x,y)
    nxny <- length(dom)
    # Calculate actual statistics A and d
    Ac <- sum(dom>0)/nxny + 0.5*sum(dom==0)/nxny
    d <- 2*Ac - 1
    
    # Cliff's delta variance calculation
    # Based on code from 'orddom' function
    n_x <- length(x)
    n_y <- length(y)
    di_<- matrix (nrow =n_x)
    dj_<- matrix (nrow =n_y)
    dw <- mean(dom)
    for (i in 1:n_x) {
      di_[i] <- mean(dom[i,]) 
    }
    for (j in 1:n_y) {
      dj_[j]<-mean(dom[,j]) 
    }
    # Variance of Cliff's delta
    delta_var <- ((((n_y*n_y)*sum((di_-dw)^2))+(((n_x*n_x)*sum((dj_-dw)^2)))-(sum((dom-dw)^2)))/(n_x*n_y*(n_x-1)*(n_y-1)))
    delta_sd <- sqrt(delta_var)
    a_sd <- delta_sd/2
    return(list("A"=Ac,"Asd"=a_sd,"D"=d,"Dsd"=delta_sd))
  } else{
    Ac <- suppressMessages(roc(cases=x,controls=y)$auc)
    d <- 2*Ac - 1
    return(list("A"=Ac,"D"=d))
  }
}

# Multi-group Vargha-delaney A
# Function by Zack Williams, 2/16,2019
# Calculate the multiple group extensions of Vargha and Delaney's A 
# Ruscio, J., & Gera, B. L. (2013). Generalizations and extensions of the probability 
# of superiority effect size estimator. Multivariate behavioral research, 48(2), 208-219.
# http://doi.org/10.1080/00273171.2012.738184
vdA.multi <- function(x,method="AAD",grpnames=NULL){
  if(!is.list(x)){
    stop("ERROR: First argument must be a list containing 3 or more sets of values.")
  } else if(length(x) < 3){
    stop("ERROR: First argument must be a list containing 3 or more sets of values.")
  }
  
  x <- lapply(x,function(vec){c(na.omit(vec))}) # Remove NA values from list entries
  ngroups <- length(x)
  # Name entries of list (Default is Grp1, Grp2, Grp3, etc.)
  if(is.null(names(x))){
    if(is.null(grpnames) | length(grpnames) != ngroups){
      names(x) <- paste0("Grp",1:ngroups)
    } else{
      names(x) <- grpnames
    }
  }
  nx <- unlist(lapply(x,length)) # Number of non-NA entries in each group
  
  if(method %in% c("AAD","apd","groupwise","Groupwise","1")){
    # A for average absolute group deviation from all other groups
    method <- "AAD"
    # A for each group
    As <- sapply(1:ngroups,function(i){
      suppressMessages(abs(roc(cases=x[[i]],controls=unlist(x[-i]))$auc - 0.5) + 0.5)
    })
    names(As) <- names(x)
    A <- mean(As)
  } else if(method %in% c("AAPD","aapd","pairwise","Pairwise","2")){
    # A for average pairwise group deviation 
    method <- "AAPD"
    # A for each pair
    As <- unlist(sapply(1:(ngroups-1),function(i){
      sapply((i+1):ngroups,function(j){
        Aij <- suppressMessages(abs(roc(cases=x[[i]],controls=x[[j]])$auc - 0.5) + 0.5)
        names(Aij) <- paste0(names(x)[i],"-",names(x)[j])
        return(Aij)
      })
    }))
    A <- mean(As)
  } else if(method %in% c("ORD","ord","Ordinal","ordinal","3")){
    # A for ordinal trend (each group larger than last)
    method <- "ORD"
    As <- sapply(1:(ngroups-1),function(i){
      Aij <- suppressMessages(cases=roc(x[[i]],controls=x[[i+1]])$auc)
      names(Aij) <- paste0(names(x)[i],"-",names(x)[i+1])
      return(Aij)
    })
    A <- mean(As)
  }else{
    stop("ERROR: Method argument not recognized. Valid options include 'AAD' (group deviation),
         'AAPD' (pairwise deviation) or 'ORD' (ordinal comparisons).")
  }
  
  obj <- list("A"=A,"A.group"=As,"N.group"=nx,"Method"=method)
  return(obj)
}

# Vargha-delaney A BCA
# Function by Zack Williams, 12/1/2018
### Calculate the SE and construct a CI for Vargha & Delaney's A using bootstrap methods
### Based on Ruscio & Mullen's (2012) Bootstrap.SE.CI.A function (see orddom package Manual)
vdA.boot <- function(x,y,B=1999,alpha=.05,A0=0.5,H0=c(cohd2vda(-0.2),cohd2vda(0.2)),p.2g=TRUE,boot.p=F,
                     seed=1){  # initialize variables
  x <- c(na.omit(x))
  y <- c(na.omit(y))
  
  set.seed(seed)
  nx <- length(x)
  ny <- length(y)
  if(boot.p==T){
    Vals <- vdA(x,y,ROC = F)
    Dsd.obs <- Vals$Dsd
  } else{
    Vals <- vdA(x,y)
  }
  A.obs <- c(Vals$A)
  D.obs <- Vals$D
  CI.Lower <- CI.Upper <- A.obs
  
  # perform bootstrap to generate B values of A
  BS <- mcmapply(function(i){
    xs <- sample(x,replace=TRUE)
    ys <- sample(y,replace=TRUE)
    if(boot.p==T){
      bootvals <- vdA(xs,ys,ROC = F)
      # Calculate ratio of (boot.delta - delta)/sd(boot.delta)
      bootvals$t <- (bootvals$D - D.obs)/bootvals$Dsd
    } else{
      bootvals <- vdA(xs,ys)
    }
    return(unlist(bootvals))
  },1:B)
  
  BS.Values <- sort(BS[1,])
  # if all bootstrap samples yield same value for A, use it for both ends of CI
  if (min(BS.Values) == max(BS.Values)){
    CI.Lower <- CI.Upper <- BS.Values[1]
    if(A0 != A.obs){
      pval = 0
    } else{
      pval = 1
    }
  } else {
    # bootstrap percentile bias corrected and accelerated
    # cf. Efron & Tibshirani (1993, Ch. 14) 
    # calculate bias-correction and acceleration parameters (z0 and a)
    z0 <- qnorm(((sum(BS.Values < A.obs) + (sum(BS.Values == A.obs)+1)/2)/(B+1)))
    # zo calculation: cf. Efron & Tibshirani (1993, p. 186, formula 14.14) 
    # However, if ties are present, original formula overestimates bias, thus tie term added
    # See https://www.medcalc.org/manual/note-bcabootstrap.php for more details
    jk <- jk_t <- rep(0, (nx + ny))
    for (i in 1:nx){
      tempjk <- vdA(x[-i], y)
      jk[i] <- tempjk$A #jackknife
      if(boot.p){
        jk_t[i] <- (tempjk$D - D.obs)/tempjk$Dsd #jackknife
      }
      for (i in 1:ny){
        tempjk <- vdA(x, y[-i])
        jk[nx+i] <- tempjk$A
        if(boot.p){
          jk_t[nx+i] <- (tempjk$D - D.obs)/tempjk$Dsd
        }
      }
    }
    Diff <- mean(jk)-jk
    a <- sum(Diff^3)/(6*(sum(Diff^2))^1.5) # cf. Efron & Tibshirani (1993, p.186/15.15 and p. 328/22.29) 
    # adjust location of endpoints, cf. Efron & Tibshirani (1993, p.185/14.10) 
    alpha1 <- pnorm(z0+((z0+qnorm(alpha/2))/(1-(a*(z0+qnorm(alpha/2))))))
    alpha2 <- pnorm(z0+((z0-qnorm(alpha/2))/(1-(a*(z0-qnorm(alpha/2))))))
    # if either endpoint undefined, replace it with value for percentile CI
    if (is.na(alpha1)) {alpha1 <- (alpha/2)}
    if (is.na(alpha2)) {alpha2 <- (1-(alpha/2))}
    if (round(alpha1*B)<1) {
      CI.Lower <- BS.Values[1]
    } else {
      CI.Lower <- BS.Values[round(alpha1 * B)]
    }
    CI.Upper <- BS.Values[round(alpha2 * B)]	
  }
  
  if(boot.p==T){
    # P value calculation for H0: A = 0.5 (added in by Zack Williams)
    # Based on guidelines from Hall, P., & Wilson, S. R. (1991). 
    # Two guidelines for bootstrap hypothesis testing. Biometrics, 757-762.
    # Basically, in order to ensure the bootstrapped distribution approximates the null,
    # we're bootstrapping the distribution of t-statistics based on the following statistic:
    #### t = (Dboot - D.obs)/Dsdboot
    # Compare to t statistic from (D.obs - D0)/Dsd.obs to get p-value for data given null
    D0 <- 2*A0 - 1
    t <- (D.obs - D0)/Dsd.obs
    BS.ts <- sort(BS[4,])
    # Calculate acceleration parameter (a; z0 is the same) for bootstrap t dist
    Diff_t <- mean(jk_t)-jk_t
    a_t <- sum(Diff_t^3)/(6*(sum(Diff_t^2))^1.5) # cf. Efron & Tibshirani (1993, p.186/15.15 and p. 328/22.29) 
    # adjust location of endpoints, cf. Efron & Tibshirani (1993, p.185/14.10) 
    quant <- (sum((BS.ts < t))+(sum((BS.ts == t))+1)/2)/(B+1) # Quantile of bootstrap dist
    b0 <- qnorm(quant) # z-score
    c0 <- ((2+a_t*b0-a_t*z0)*z0-b0)/(1+a_t*b0-a_t*z0) # transform z-score using BCa parameters
    p0 <- pnorm(c0) # P-value for a left-tailed test
    pc <- 2*min(p0,1-p0) # 2-tailed p value
    p <- signif(pc,digits=3)    
  } else{
    p <- NULL
  }

  # return A, SE of A, lower limit of CI, upper limit of CI
  obj <- list(A=c(A.obs),D.obs=D.obs,A.CI.Lower=CI.Lower,A.CI.Upper=CI.Upper,
              D.CI.Lower=2*CI.Lower-1,D.CI.Upper=2*CI.Upper-1,A0=A0,alpha=alpha,
              p.val=p,n1=nx,n2=ny,nboot=B,call=match.call())
  
  # 2nd generation p-value
  if(p.2g==T){
    obj$p.val.2g <- pval.2g(c(CI.Lower,CI.Upper),H0[1],H0[2])
    obj$H0.low <- signif(H0[1],3)
    obj$H0.high <- signif(H0[2],3)
  }
  
  obj$n.comp <- obj$n1 * obj$n2
  class(obj) <- "A.boot"
  return(obj)
}

# Print Function for bootstrap Vargha-Delaney A
print.A.boot <- function(x,digits=3){
  cat("Call: ")
  print(x$call)
  cat("\n")
  if(!is.null(x$method)){
    cat(paste0("Multi-group Vargha-Delaney A (",x$method," method; k = ",length(x$N.group),")\n"))
    cat(paste0(" A = ",round(x$A,digits)," (",100*(1-x$alpha),"% CI:",
               paste(format(round(x$A.CI.Lower,digits=digits),nsmall=1),
                     format(round(x$A.CI.Upper,digits=digits),nsmall=1),sep="—"),")\n"))
  } else{
    cat(paste0("Vargha-Delaney A (Y>X) = ",round(x$A,digits)," (",100*(1-x$alpha),"% CI:",
               paste(format(round(x$A.CI.Lower,digits=digits),nsmall=1),
                     format(round(x$A.CI.Upper,digits=digits),nsmall=1),sep="—"),")\n"))
  }
  cat(" BCa confidence based on",x$nboot,"bootstrapped samples\n")
  if(!is.null(x$p.val)){
    cat(paste0(" Pr(>|A-",x$A0,"|) = ",signif(x$p.val,digits=3)))
    cat("\n")    
  }
  if(is.null(x$method)){
    cat(paste0("Approximate Cohen's d = ",round(auc2cohd(x$A),digits)," (",100*(1-x$alpha),"% CI:",
               paste(format(round(auc2cohd(x$A.CI.Lower),digits=digits),nsmall=1),
                     format(round(auc2cohd(x$A.CI.Upper),digits=digits),nsmall=1),sep="—"),")\n"))
  }
  if(!is.null(x$p.val.2g)){
    cat(paste0(" 2nd-generation p-value = ",x$p.val.2g," with interval H0 [",x$H0.low,", ",x$H0.high,"]"))
    if(x$p.val.2g==0){cat(" (reject interval null)")}
    else if(x$p.val.2g == 0){cat(" (accept equivalence to interval null)")}
    else {cat(" (inconclusive)")}    
  }
}

# Multi-group Vargha-delaney A with BCa bootstrapped CI
# Function by Zack Williams, 2/16,2019
# Calculate the multiple group extensions of Vargha and Delaney's A 
# Ruscio, J., & Gera, B. L. (2013). Generalizations and extensions of the probability 
# of superiority effect size estimator. Multivariate behavioral research, 48(2), 208-219.
# http://doi.org/10.1080/00273171.2012.738184
vdA.multi.boot <- function(x,method="AAD",grpnames=NULL,B=1999,alpha=.05,seed=1){
  set.seed(seed)
  Vals <- vdA.multi(x,method=method,grpnames=grpnames)
  nx <- Vals$N.group
  A.obs <- Vals$A
  A.group <- Vals$A.group
  CI.Lower <- CI.Upper <- A.obs
  
  ngroups <- length(x)
  # Name entries of list (Default is Grp1, Grp2, Grp3, etc.)
  if(is.null(names(x))){
    if(is.null(grpnames) | length(grpnames) != ngroups){
      names(x) <- paste0("Grp",1:ngroups)
    } else{
      names(x) <- grpnames
    }
  }
  
  # perform bootstrap to generate B values of A
  BS <- mcmapply(function(i){
    xb <- lapply(x,sample,replace=TRUE)
    bootvals <- vdA.multi(xb,method=method,grpnames=grpnames)$A
    return(unlist(bootvals))
  },1:B)
  
  BS.Values <- sort(BS)
  # if all bootstrap samples yield same value for A, use it for both ends of CI
  if (min(BS.Values) == max(BS.Values)){
    CI.Lower <- CI.Upper <- BS.Values[1]
  } else {
    # bootstrap percentile bias corrected and accelerated
    # cf. Efron & Tibshirani (1993, Ch. 14) 
    # calculate bias-correction and acceleration parameters (z0 and a)
    z0 <- qnorm(((sum(BS.Values < A.obs) + (sum(BS.Values == A.obs)+1)/2)/(B+1)))
    # zo calculation: cf. Efron & Tibshirani (1993, p. 186, formula 14.14) 
    # However, if ties are present, original formula overestimates bias, thus tie term added
    # See https://www.medcalc.org/manual/note-bcabootstrap.php for more details
    jk <- unlist(mapply(function(grp){
      mcmapply(function(i){
        if(grp != ngroups){
          jklist <- c(x[0:(grp-1)],list(x[[grp]][-i]),x[(grp+1):(ngroups)])
        }else{
          jklist <- c(x[0:(grp-1)],list(x[[grp]][-i]))
        }
        return(vdA.multi(jklist,method=method)$A)
      },1:nx[grp])
    },1:ngroups))
    
    Diff <- mean(jk)-jk
    a <- sum(Diff^3)/(6*(sum(Diff^2))^1.5) # cf. Efron & Tibshirani (1993, p.186/15.15 and p. 328/22.29) 
    # adjust location of endpoints, cf. Efron & Tibshirani (1993, p.185/14.10) 
    alpha1 <- pnorm(z0+((z0+qnorm(alpha/2))/(1-(a*(z0+qnorm(alpha/2))))))
    alpha2 <- pnorm(z0+((z0-qnorm(alpha/2))/(1-(a*(z0-qnorm(alpha/2))))))
    # if either endpoint undefined, replace it with value for percentile CI
    if (is.na(alpha1)) {alpha1 <- (alpha/2)}
    if (is.na(alpha2)) {alpha2 <- (1-(alpha/2))}
    if (round(alpha1*B)<1) {
      CI.Lower <- BS.Values[1]
    } else {
      CI.Lower <- BS.Values[round(alpha1 * B)]
    }
    CI.Upper <- BS.Values[round(alpha2 * B)]	
  }
  
  # return A, SE of A, lower limit of CI, upper limit of CI
  obj <- list(A=A.obs,A.CI.Lower=CI.Lower,A.CI.Upper=CI.Upper,alpha=alpha,
              A.group,N.group=nx,method=method,nboot=B,call=match.call())
  class(obj) <- "A.boot"
  return(obj)
}


# Equivalence testing using Cliff's Delta 
# Two One-sided Tests (TOST) method
# Lakens, D. (2017). Equivalence tests: a practical primer for t tests, correlations, 
# and meta-analyses. Social Psychological and Personality Science, 8(4), 355-362.
# H0 consists of equivalence bounds (in cliff delta units) representing smallest ES of interest
cliff.d.TOST <- function(x,y,H0=c(cohd2delta(-0.2),cohd2delta(0.2)),alpha=0.05){
  x <- c(na.omit(x))
  y <- c(na.omit(y))
  signlev <- alpha
  
  # Dominance Matrix (from Orddom package)
  dm <- function (x,y) { #produces dominance matrix
    dx <- matrix(nrow=length(x), ncol=length(y))
    for (i in 1:length(x)){
      for (j in 1:length(y)){
        dx[i,j] <- sign(y[j]-x[i])
      }
    } 
    rownames(dx) <- x
    colnames(dx) <- y
    return(dx)
  }
  dom <- dm(x,y)
  nxny <- length(dom)
  
  # Calculate cliff's delta (dw in this code)
  dw <- sum(dom>0)/nxny - sum(dom<0)/nxny
  
  # Cliff's delta variance calculation
  # Based on code from 'orddom' function
  n_x <- length(x)
  n_y <- length(y)
  di_<- matrix (nrow =n_x)
  dj_<- matrix (nrow =n_y)
  for (i in 1:n_x) {
    di_[i] <- mean(dom[i,]) 
  }
  for (j in 1:n_y) {
    dj_[j]<-mean(dom[,j]) 
  }
  # Variance of Cliff's delta
  delta_var <- ((((n_y*n_y)*sum((di_-dw)^2))+(((n_x*n_x)*sum((dj_-dw)^2)))-(sum((dom-dw)^2)))/(n_x*n_y*(n_x-1)*(n_y-1)))

  #Confidence Intervals
  #cf. Feng & Cliff (2004), Journal of Modern Applied Statistical Methods, 3(2), 322-332 and
  #cf. Feng (2007), in:  Shlomo S. Sawilowsky (Ed.), Real Data Analysis (pp. 163-183).
  #here derived fm t-test as qt((1-signlev),degrees of freedom).
  t_level <- qt(1-(signlev),n_x+n_y-2) #Student t-Distribution
  
  #asymmetric CI
  ci_dw_lo <- (dw-(dw^3)-(t_level*sqrt(delta_var)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*delta_var)))) / (1-(dw*dw)+(t_level*t_level*delta_var))
  ci_dw_hi <- (dw-(dw^3)+(t_level*sqrt(delta_var)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*delta_var)))) / (1-(dw*dw)+(t_level*t_level*delta_var))
  #special case Feng & Cliff (2004), eq. 8, p.324: d=+-1
  if ((dw*dw)==1) {n_b<-n_x
  if (n_y<n_x){n_b<-n_y}}
  if (dw==1)  {ci_dw_lo <-((n_x-(t_level*t_level))/(n_x+(t_level*t_level)))
  if (delta_var==0){ci_dw_hi<-Inf}}
  if (dw==-1) {ci_dw_hi <--((n_x-(t_level*t_level))/(n_x+(t_level*t_level)))
  if (delta_var==0){ci_dw_lo<--Inf}}
  if (ci_dw_lo < -1) {ci_dw_lo<--1}
  if (ci_dw_hi > 1) {ci_dw_hi<-1}
  #t and p levels
  # Zero point (two-tailed)
  tdw <- dw/sqrt(delta_var)
  ptdw <- 2*pt(-abs(tdw),n_x+n_y-2)
  
  # Functions to optimize to get one-tailed p-vals with asymmetric CI
  # Lower bound p-val (one tailed)
  findp_lo <- function(p){
    t_level <- qt(1-(p),n_x+n_y-2) #Student t-Distribution
    #asymmetric CI
    ci_dw_lo <- (dw-(dw^3)-(t_level*sqrt(delta_var)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*delta_var)))) / (1-(dw*dw)+(t_level*t_level*delta_var))
    return(abs(ci_dw_lo-H0[1]))  
  }
  # Upper bound p-val (one tailed)
  findp_hi <- function(p){
    t_level <- qt(1-(p),n_x+n_y-2) #Student t-Distribution
    #asymmetric CI
    ci_dw_hi <- (dw-(dw^3)+(t_level*sqrt(delta_var)*sqrt(1-(2*dw*dw)+(dw^4)+(t_level*t_level*delta_var)))) / (1-(dw*dw)+(t_level*t_level*delta_var))
    return(abs(ci_dw_hi-H0[2]))  
  }
  
  pdw <- c(ptdw,optimize(findp_lo,c(0,1))$minimum,optimize(findp_hi,c(0,1))$minimum)
  
  cliff <- as.list(c(n_x,n_y,sum(dom==1),sum(dom==-1),sum(dom==0),dw,100*(1-2*alpha),ci_dw_lo,ci_dw_hi,
                   n_x+n_y-2,H0[1],H0[2],pdw[1],pdw[2],pdw[3]))
  names(cliff) <- c("n1","n2","N(Y>X)","N(Y=X)","N(Y<X)","delta","ci.size","delta.ci.low","delta.ci.high",
                    "df","ESbound.low","ESbound.high","p.val","p.val.low","p.val.high")
  cliff$n.comp <- cliff$n1 * cliff$n2
  cliff$"P(Y>X)" <- cliff$"N(Y>X)"/cliff$n.comp
  cliff$"P(Y<X)" <- cliff$"N(Y<X)"/cliff$n.comp
  cliff$"P(Y=X)" <- cliff$"N(Y=X)"/cliff$n.comp
  cliff$paired <- FALSE
  cliff$Call <- match.call()
  class(cliff) <- "cliff.delta"
  cliff$summ.X <- summ(x,"X")
  cliff$summ.Y <- summ(y,"Y")
  
  # 2nd generation p-value
  cliff$H0.low <- signif(H0[1],3)
  cliff$H0.high <- signif(H0[2],3)
  
  return(cliff)
}

cliff.d.power <- function(coh.d=NULL,n1,n2,TOST = F,n.iter=10000,alpha=0.05,
                          H0=c(-0.1475576,0.1475576),delta=0,H0.coh.d=NULL,seed=12345){
  set.seed(seed)
  # Note that H0 (equivalence bounds) are in Cliff delta units
  # Bounds in Cohen's d units can be set with H0.coh.d parameter
  if(!is.null(H0.coh.d)){ # Group difference (in cohen's d by default can be set in delta units with delta)
    H0 <- cohd2delta(H0.coh.d)
  }
  if(is.null(coh.d)) coh.d <- delta2cohd(delta)
  if(TOST == T){
    ps <- mcmapply(function(x){
      g1 <- rnorm(n1,0,1)
      g2 <- rnorm(n2,coh.d,1)
      tost <- cliff.d.TOST(g1,g2,H0=H0,alpha=alpha)
      return(max(tost$p.val.low,tost$p.val.high))
    },1:n.iter)
  } else{ # Power analysis for one-sided hypothesis test
    ps <- mcmapply(function(x){
      g1 <- rnorm(n1,0,1)
      g2 <- rnorm(n2,coh.d,1)
      return(cliff.d.TOST(g1,g2,alpha=alpha)$p.val)
    },1:n.iter)
  }
  return(signif(sum(ps < alpha)/n.iter,3))
}

spearman.power <- function(rs,n,TOST = F,n.iter=10000,alpha=0.05,s2="cc",H0=c(-0.1,0.1)){
  # Note that H0 (equivalence bounds) are in Spearman correlation units
  require(MASS)
  if(is.null(rs)) rs <- rp <- 0 
  else{
    rp <- 2*sin(rs*pi/6) # get pearson from spearman corr
  }
  S <- matrix(c(1,rp,rp,1),ncol=2)
  if(TOST == T){
    ps <- mcmapply(function(x){
      s <- mvrnorm(n, c(0,0), S)
      rs_s <- Spearman.test(s[,1],s[,2],alpha=alpha,s2=s2,H0=H0)
      return(rs_s[5])
    },1:n.iter)
  } else{ # Power analysis for one-sided hypothesis test
    ps <- mcmapply(function(x){
      s <- mvrnorm(n, c(0,0), S)
      rs_s <- Spearman.test(s[,1],s[,2],alpha=alpha,s2=s2,H0=H0)
      return(rs_s[4])
    },1:n.iter)
  }
  return(signif(sum(ps < alpha)/n.iter,3))
}

### Permutation rank test based on maximum of Brunner-Munzel/Rank Welch test statistics
## Zack Williams
## Last updated 02.24.2019
# Based on code from https://datadryad.org/bitstream/handle/10255/dryad.169771/MaximumTestExact.r?sequence=1
# Test described in Welz, Ruxton & Neuhäuser (2018) 
# A non-parametric maximum test for the Behrens–Fisher problem, 
# Journal of Statistical Computation and Simulation, 88:7, 1336-1347
# https://doi.org/10.1080/00949655.2018.1431236
# If # permuted values greater than test >= 10, can use typical permutation p-value method
# Else approximate p-value using generalized Pareto distribution
# See Knijnenburg et al., (2009) for more details
# http://doi.org/10.1093/bioinformatics/btp211

## Function to calculate the maximum standardized test statistic
maxtest <- function(group1,group2){
  n <- length(group1)
  m <- length(group2)
  N <- n + m
  
  allvals <- c(group1,group2)
  ranks <- rank(allvals,ties.method="average")
  rank1 <- ranks[1:n]
  rank2 <- ranks[(n+1):N]
  
  ### Calculation of test statistics on observed data
  #---------Welch-Test on ranks----------------
  Trank <- t.test(rank2,rank1,var.equal=F)    #Welch-Test on ranks
  df_w <- Trank$parameter
  p_Trankexact <- Trank$p.value
  Trank_results <- c("Rank.t"=unname(Trank$statistic),df_w,"p.val"=p_Trankexact)
  # When df < 2, implement corrections to not break variance calc
  if(Trank$parameter<2){df_w=2.01}
  # SD of test statistic
  varTrank <- df_w/(df_w-2)
  sdTrank <- sqrt(varTrank)
  
  # Standardized statistic
  Trankstd <- Trank$statistic/sdTrank
  
  #----------Brunner-Munzel-Test-------
  BM <- lawstat::brunner.munzel.test(group1, group2) #from Package lawstat
  BMstat <- BM$statistic
  df <- BM$parameter
  p_BMexact <- BM$p.value
  BM_results <- c("BM.Stat"=unname(BMstat),df,"p.val"=p_BMexact,BM$estimate,
                  "ci.lo"=unname(BM$conf.int[1]),"ci.hi"=unname(BM$conf.int[2]))
  # When BM is infinite or df < 2, implement corrections to not break variance calc
  if(BM$statistic!=Inf & BM$statistic!= -Inf){if(df<2){df=2.01}}
  if(BM$statistic==Inf | BM$statistic== -Inf){
    df=2*(n-1)*(m-1)/(N-2)
    if(df<2){df=2.01}
    Rmean1=mean(rank1)
    Rmean2=mean(rank2)
    sigma2=N/(2*n*m)
    BMstat= sqrt((n*m)/N)*((Rmean1-Rmean2)/sqrt(sigma2))
    p_BMexact=0
  }
  varW=df/(df-2)
  sdW=sqrt(varW)
  
  # Standardized statistic
  BMstd <- BMstat/sdW
  
  # Maximum test statistic
  MAXstat <- max(abs(Trankstd),abs(BMstd))
  max_test <- unname(ifelse(abs(Trankstd)>abs(BMstd),"Trank","BM"))
  
  return(list(MAXstat=MAXstat,MAX.test=max_test,"BM.test"=BM_results,"Rank.Welch.test"=Trank_results))
}

# Calculate delta (for the effect size A described in below paper)
# A = 1 - delta/(mean of permuted deltas)
# See Berry, K. J., & Mielke Jr, P. W. (1992). 
# A family of multivariate measures of association for nominal independent variables. 
# Educational and Psychological Measurement, 52(1), 41-55.
# https://doi.org/10.1177/001316449205200104
calc_delta <- function(x,y){
  nx <- length(x)
  ny <- length(y)
  Nxy <- nx + ny
  delta <- weighted.mean(c(GiniMd(x),GiniMd(y)),c(nx,ny))
}

# Do the permutation test, make a nice object, etc.
MAX.perm.test=function(x,y,Nperm=2000,plot=T,fill="gray80",
                       linecol="blue",grpnames=c("Group 1","Group 2"),cores=detectCores()-1){
  require(lawstat)
  require(eva)
  X <- c(na.omit(x))
  Y <- c(na.omit(y))
  
  nx <- length(X)
  ny <- length(Y)
  Nxy <- nx + ny
  
  # Maximum statistic for observed data
  MAX <- maxtest(X,Y)
  MAXstat <- MAX$MAXstat
  # Delta from observed data
  delta <- calc_delta(X,Y)
  
  # Permutation test
  perm_vals <- pbsapply(1:Nperm,function(x){
    
    perm <- base::sample(c(X,Y))
    perm_1 <- perm[1:nx]
    perm_2 <- perm[(nx+1):Nxy]
    
    c(maxtest(perm_1,perm_2)$MAXstat,calc_delta(perm_1,perm_2))
  },cl=cores)
  
  perms <- sort(perm_vals[1,]) # Test statistics
  deltas <- perm_vals[2,]

  # M: Number of permutations greater than MAXstat
  M <- sum(perms>MAXstat)
  # If M >= 10, can use typical permutation p-value method
  # Else approximate p-value using generalized Pareto distribution
  # See Knijnenburg et al., (2009) for more details
  # http://doi.org/10.1093/bioinformatics/btp211
  if(M>=10){
    pval <- M/length(perms)
    sd_pval <- sqrt(pval * (1-pval)/length(perms))
    pval.ci.95 <- pval + qnorm(c(0.025,0.975)) * sd_pval
  } else{ # Find p-value with tail approximation from generalized Pareto distribution
    permEVs <- perms[(Nperm-249):Nperm]
    AD <- gpdAd(permEVs) # Anderson-darling GoF test for pareto
    # If poor fit to distribution, cut off 10 values and try again
    while(length(permEVs) > 20 & AD$p.value < 0.05){
      permEVs <- permEVs[-1:-10]
      AD <- try(gpdAd(permEVs),silent = T)
      # If need to bootstrap, do so (500 reps)
      if(class(AD)=="try-error"){
        message("Bootstrapping AD test. May take a bit longer... ")
        if(grepl("bootstrap",attr(AD,"condition")$message)){
          AD <- gpdAd(permEVs,bootstrap = T,allowParallel = T,bootnum=500,numCores=cores)
        }
      }
    }
    # If still a problem after all that, use different p-value method
    if(class(AD)=="try-error" | AD$p.value < 0.05){
      warning("Generalized Pareto distribution did not fit well. Estimating p-value using (M+1)/(NPerm+1).")
      pval <- (M+1)/(Nperm+1)
    } else{
      # Number of Exceedences
      nExc <- length(permEVs)
      # Exceedences Threshold
      tExc <- mean(perms[(Nperm-nExc):(Nperm-nExc+1)])
      # Calculate P-value from pareto distribution
      pval <- unname(nExc/Nperm * (1 - eva::pgpd(MAXstat,
                                                 loc=tExc,scale=AD$theta[1],shape=AD$theta[2])))
    }
  }

  # Effect size: Berry & Mielke's A (https://doi.org/10.1177/001316449205200104)
  ES <- 1 - delta/mean(deltas)
  
  if(plot==T){
    p <- ggplot(data=as.data.frame(perms),aes(x=perms)) +
      geom_density(colour="black",fill=fill) +
      geom_vline(xintercept=MAXstat,linetype="dashed",color=linecol,size=1) +
      geom_rug() +
      theme(legend.position="none",
            plot.title = element_text(size=14, face="bold", hjust=0.5),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            axis.title.x = element_text(size=16,face="bold"),
            axis.title.y = element_text(size=16,face="bold"),
            strip.text.y = element_text(size = 20, colour = "white"),
            strip.background = element_rect(colour="darkgrey", fill="darkgrey")) +
      ylab("Density") + xlab("Permuted Test Statistics") + 
      ggtitle(paste0("Permuted Rank Test: ",grpnames[1]," vs. ",grpnames[2]))
    # Visualization of permutation test
    p <- p + annotate("label", size=5, x = MAXstat, y = quantile(density(perms)$y,0.925), 
                      label = paste0("T = ",round(MAXstat,3),"\nP = ",round(pval,3)))
    print(p)
  }
  
  result <- list("Nx"=nx,"Ny"=ny,"MAXstat"=MAXstat,"p.val"=pval,"ES"=ES,"n.perm"=Nperm,
                 "BM.test"=MAX$BM.test,"Rank.Welch.test"=MAX$Rank.Welch.test)
  result$Call <- match.call()
  result$summ.X <- summ(X,"X")
  result$summ.Y <- summ(Y,"Y")
  
  class(result) <- "BM.perm"
  return(result)
}

# Pretty printing for permutation test
print.BM.perm <- function(x,digits=3,summary=T){
  cat("Call: ")
  print(x$Call)
  cat("\n")
  if(summary == T){
    cat("Nx =", x$Nx,"| Ny =",x$Ny,"| Total N =",x$Nx + x$Ny,"\n")
    print(x$summ.X)
    print(x$summ.Y)
    cat("\n") 
  }
  cat(paste0("Maximum permutation test (B = ",x$n.perm,") of Brunner Munzel/Rank Welch test statistics:\n"))
  cat(paste0(" T = ",round(x$MAXstat,digits)," | Permutation Effect Size (A) = ",round(x$ES,digits),"\n"))
  cat(paste0(" Pr(>|T|) = ",signif(x$p.val,digits=3)))
  if(x$p.val<0.05){cat("*")}
  if(x$p.val<0.01){cat("*")}
  if(x$p.val<0.005){cat("*")}
  
  cat("\n\nProbability of stochastic superiority of Y over X: (0<A<1, 0.5 = no dominance)")
  cat("\n Vargha & Delaney's A =",round(x$BM.test[4],digits))
  cat(paste0(" (95% CI: ",
             paste(format(round(x$BM.test[5],digits=digits),nsmall=1),
                   format(round(x$BM.test[6],digits=digits),nsmall=1),sep="—"),")\n")) 
  cat(" Approximate Cohen's d =",round(auc2cohd(x$BM.test[4]),digits))
  cat(paste0(" (95% CI: ",
             paste(format(round(auc2cohd(x$BM.test[5]),digits=digits),nsmall=1),
                   format(round(auc2cohd(x$BM.test[6]),digits=digits),nsmall=1),sep="—"),")\n")) 
  
  cat("\nAdditional Significance Tests (non-permuted):")
  cat("\n Brunner Munzel Test:")
  cat(paste0("  t = ",round(x$BM.test[1] ,digits),","))
  cat(paste0(" df = ",round(x$BM.test[2],2),", Pr(>|t|) = ",signif(x$BM.test[3],digits=3)))
  if(x$p.val<0.05){cat("*")}
  if(x$p.val<0.01){cat("*")}
  if(x$p.val<0.005){cat("*")}
  
  cat("\n Rank Welch Test:")
  cat(paste0("  t = ",round(x$Rank.Welch.test[1] ,digits),","))
  cat(paste0(" df = ",round(x$Rank.Welch.test[2],2),", Pr(>|t|) = ",signif(x$Rank.Welch.test[3],digits=3)))
  if(x$p.val<0.05){cat("*")}
  if(x$p.val<0.01){cat("*")}
  if(x$p.val<0.005){cat("*")}
}

# Brunner Munzel Equivalence Testing
# Zack Williams
# March 7th, 2019
# Brunnzer Munzel Test code from `lawstat` package
BM.TOST <-function (x, y, alpha = 0.05, H0 = c(cohd2vda(-0.2),cohd2vda(0.2))){
  x <- na.omit(x)
  y <- na.omit(y)
  n1 = length(x)
  n2 = length(y)
  N <- n1 + n2
  r1 = rank(x)
  r2 = rank(y)
  r = rank(c(x, y))
  # Mean Ranks
  m1 = mean(r[1:n1])
  m2 = mean(r[n1 + 1:n2])
  pst = (m2 - (n2 + 1)/2)/n1 # Estimate of P hat (Vargha/Delaney A)
  # Variances
  v1 = sum((r[1:n1] - r1 - m1 + (n1 + 1)/2)^2)/(n1 - 1) 
  v2 = sum((r[n1 + 1:n2] - r2 - m2 + (n2 + 1)/2)^2)/(n2 - 1)
  vn <- N * (v1/n2 + v2/n1) # Pooled variance
  sn <- sqrt(vn) # Pooled SD
  # Variance in the metric of P hat
  vp <- vn/(N * n1 * n2)
  sp <- sqrt(vp)
  #  Brunner-Munzel Test Statistic
  statistic <- (m2 - m1)/sn * sqrt(n1*n2/N) # Equal to (P hat - 0.5)/sp
  # Estimated Degrees of Freedom
  dfbm <- ((n1 * v1 + n2 * v2)^2)/(((n1 * v1)^2)/(n1 - 1) + 
                                     ((n2 * v2)^2)/(n2 - 1))
  # 2-tailed P-value
  p.val = 2 * min(pt(abs(statistic), dfbm), (1 - pt(abs(statistic),dfbm)))
  # 90% Confidence interval: P hat ± tcrit * Vn/(N * n1 * n2)
  conf.int = c(pst - qt(1 - alpha, dfbm) * sqrt(v1/(n1 * n2^2) + v2/(n2 * n1^2)),
               pst + qt(1 - alpha, dfbm) * sqrt(v1/(n1 * n2^2) + v2/(n2 * n1^2)))
  # Equivalence Testing (TOST)
  p.val.lower <- 1 - pt((pst - H0[1])/sp,dfbm)
  p.val.upper <- pt((pst - H0[2])/sp,dfbm)
  p.equiv <- max(p.val.lower,p.val.upper)
  
  
  BM <- structure(as.list(c("Nx"=n1,"Ny"=n2,"BM.t"=statistic,"df"=dfbm,"ci.size"=100*(1-2*alpha),
                            "H0.lo"=H0[1],"H0.hi"=H0[2],"p.val"=p.val,"p.equiv"=p.equiv,
                            "p.equiv.lo"=p.val.lower,"p.equiv.hi"=p.val.upper,
                            "VD.A"=pst,"A.ci.lo"=conf.int[1],"A.ci.hi"=conf.int[2],"Call"=match.call())),
                  class = "BM.test")
  BM$summ.X <- summ(x)
  BM$summ.Y <- summ(y)
  
  return(BM)
}

# Pretty printing for Brunner-Munzel Test
print.BM.test <- function(x,digits=3,summary=T){
  cat("Call: ")
  print(x$Call)
  cat("\n")
  if(summary == T){
    cat("Nx =", x$Nx,"| Ny =",x$Ny,"| Total N =",x$Nx + x$Ny,"\n")
    print(x$summ.X)
    print(x$summ.Y)
    cat("\n") 
  }
  
  cat(paste0("Brunner-Munzel Test (Y>X):\n"))
  cat(paste0(" T = ",round(x$BM.t,digits)," df = ",round(x$df,2),"\n"))
  cat(paste0(" Pr(>|T|) = ",signif(x$p.val,digits=3)))
  if(x$p.val<0.05){cat("*")}
  if(x$p.val<0.01){cat("*")}
  if(x$p.val<0.005){cat("*")}
  cat("\n")
  print(p.transform(x$p.val,digits = digits)[-3])
  
  if(!is.null(x$VD.A)){ # equivalence test results
    cat("\nProbability of stochastic superiority of Y over X: (0<A<1, 0.5 = no dominance)")
    cat("\n Vargha & Delaney's A =",round(x$VD.A,digits))
    cat(paste0(" (",x$ci.size,"% CI: ",
               paste(format(round(x$A.ci.lo,digits=digits),nsmall=1),
                     format(round(x$A.ci.hi,digits=digits),nsmall=1),sep="—"),")\n"))    
  }
  
  cat(" Qualitative Effect Size: ") # Based on Cohen's guidelines and conversion
  if(abs(x$VD.A - 0.5) < cohd2vda(0.2)-0.5){
    cat("negligible (|A-0.5|<0.056)")
  } else if(abs(x$VD.A - 0.5) < cohd2vda(0.5)-0.5){
    cat("small (|A-0.5|<0.138)")
  } else if(abs(x$VD.A - 0.5) < cohd2vda(0.8)-0.5){
    cat("medium (|A-0.5|<0.214)")
  } else{
    cat("large (|A-0.5|≥0.214)")
  }
  if(!is.null(x$p.equiv)){ # equivalence test results
    cat(paste0("\n\nEquivalence tests (TOST) with equivalence bounds [",round(x$H0.lo,3),", ",round(x$H0.hi,3),"]"))
    cat(paste0("\n Pr(equivalence) = ",signif(x$p.equiv,digits=3)))
    if(x$p.equiv.lo < 0.05 & x$p.equiv.hi < 0.05){
      cat("\n Equivalence accepted at p < 0.05 level")
    } else{
      cat("\n Equivalence not accepted (p > 0.05)")
    }
  }
}

###############################################################################
# Analysis of Credibility (AnCred) Function
# Zack Williams
# 3/25/2019
# See https://doi.org/10.1098/rsos.171047; https://doi.org/10.1098/rsos.181534; https://doi.org/10.1080/00031305.2018.1543136
AnCred <- function(CI,M=NULL,type="SL",OR=F){
  L <- CI[1]
  U <- CI[2]
  
  # Intrinsic credibility - if M provided
  if(!is.null(M)){
    if(OR == T){
      skep <- exp(log(U/L)^2/(4*sqrt(log(U)*log(L))))
      SL <- c("SL.low"=1/skep,"SL.high"=skep)
    } else{
      skep <- (U-L)^2/(4*sqrt(U*L))
      SL <- c("SL.low"=-skep,"SL.high"=skep)
    }
    # Assess Intrinsic Credibility
    int.cred <- ifelse(unname(M < SL[1] | M > SL[2]),"YES","NO")
    return(c(round(c(SL,"M"=M),3),"Int.Cred"=int.cred))
  }
  
  # For statistically significant findings - Skepticism Level (SL)
  if(type == "SL"){
    if(OR == T){
      skep <- exp(log(U/L)^2/(4*sqrt(log(U)*log(L))))
      SL <- c("SL.low"=1/skep,"SL.high"=skep)
    } else{
      skep <- (U-L)^2/(4*sqrt(U*L))
      SL <- c("SL.low"=-skep,"SL.high"=skep)
    }
    return(SL)
  }
  # For nonsignificant findings - Advocacy Limit (AL)
  if(type == "AL"){
    if(OR == T){
      adv <- exp(log(U*L)*log(U/L)^2/(2*log(U)*log(L)))
      AL <- c("AL.low"=1,"AL.high"=adv)
    } else{
      adv <- -(U+L)/(2*U*L)*(U-L)^2
      AL <- c("AL.low"=0,"AL.high"=adv)
    }
    return(AL)
  }
}

### Bootstrap Correlation Difference Test
# Zack Williams
# 9/21/2019
# Little function to calculate the difference between two independent correlation coefficients and bootstrap it
# Uses BCa bootstrap intervals
# Can handle multiple different input types
#   - x and y are two different columns containing variables to be correlated
#   - x is matrix/df containing two variables to be correlated
# Can accept point null hypothesis (default 0) or interval null
cor.diff.boot <- function(x,y=NULL,groupvar,method="pearson",B=10000,CI.width=0.95,H0=0,
                          nCores=detectCores()-1,seed=256,varnames=NULL){
  set.seed(seed)
  call <- match.call()
  options("mc.cores"=nCores)
  # Handle different data types
  if(!is.null(y)){
    data <- cbind(x,y)
  } else{
    if(ncol(x)!=2){
      warning("'x' has more than two columns. Assuming first two columns contain variables of interest.\n")
      data <- x[,1:2]
    } else{
      data <- x
    }
  }
  if(!is.null(varnames)){colnames(data) <- varnames}
  # Make sure groupvar is the right length and has only two levels
  if(length(groupvar)!=nrow(data)){stop("ERROR: 'groupvar' is not the same length as supplied variables.")}
  if(length(unique(as.character(groupvar)))!=2){
    stop("ERROR: 'groupvar' should only have two levels.")
  } else{
    grp1 <- unique(as.character(groupvar))[1]
    grp2 <- unique(as.character(groupvar))[2]
  }
  # Delete cases with NAs
  which_cases <- complete.cases(data)
  cor_data <- data[which_cases,]
  cor_grps <- groupvar[which_cases]
  # Declare some useful variables
  alpha <- 1- CI.width
  N <- nrow(cor_data)
  grp1_data <- cor_data[which(cor_grps==grp1),]
  grp2_data <- cor_data[which(cor_grps==grp2),]
  n1 <- nrow(grp1_data)
  n2 <- nrow(grp2_data)
  
  # Point estimates of correlations and difference
  grp1_cor <- cor(grp1_data[,1],grp1_data[,2],method=method)
  grp2_cor <- cor(grp2_data[,1],grp2_data[,2],method=method)
  r_diff <- grp2_cor - grp1_cor
  
  # perform bootstrap to generate B values of r_diff
  BS <- mcmapply(function(i){
    b1 <- grp1_data[sample(1:n1,replace=TRUE),]
    b2 <- grp2_data[sample(1:n2,replace=TRUE),]
    b_diff <- cor(b2[,1],b2[,2],method=method) - cor(b1[,1],b1[,2],method=method)
    return(b_diff)
  },1:B)
  
  BS.Values <- sort(BS)
  # if all bootstrap samples yield same value for A, use it for both ends of CI
  if (min(BS.Values) == max(BS.Values)){
    CI.Lower <- CI.Upper <- BS.Values[1]
  } else {
    # bootstrap percentile bias corrected and accelerated
    # cf. Efron & Tibshirani (1993, Ch. 14) 
    # calculate bias-correction and acceleration parameters (z0 and a)
    z0 <- qnorm(((sum(BS.Values < r_diff) + (sum(BS.Values == r_diff)+1)/2)/(B+1)))
    # zo calculation: cf. Efron & Tibshirani (1993, p. 186, formula 14.14) 
    # However, if ties are present, original formula overestimates bias, thus tie term added
    # See https://www.medcalc.org/manual/note-bcabootstrap.php for more details
    jk <- rep(0,N)
    jk <- mcmapply(function(i){
      data_jk <- cor_data[-i,]
      grps_jk <- cor_grps[-i]
      return(cor(data_jk[grps_jk==grp2,],method=method)[1,2] - cor(data_jk[grps_jk==grp1,],method=method)[1,2])
    },1:N)
    
    Diff <- mean(jk)-jk
    a <- sum(Diff^3)/(6*(sum(Diff^2))^1.5) # cf. Efron & Tibshirani (1993, p.186/15.15 and p. 328/22.29) 
    # adjust location of endpoints, cf. Efron & Tibshirani (1993, p.185/14.10) 
    alpha1 <- pnorm(z0+((z0+qnorm(alpha/2))/(1-(a*(z0+qnorm(alpha/2))))))
    alpha2 <- pnorm(z0+((z0-qnorm(alpha/2))/(1-(a*(z0-qnorm(alpha/2))))))
    # if either endpoint undefined, replace it with value for percentile CI
    if (is.na(alpha1)) {alpha1 <- (alpha/2)}
    if (is.na(alpha2)) {alpha2 <- (1-(alpha/2))}
    if (round(alpha1*B)<1) {
      CI.Lower <- BS.Values[1]
    } else {
      CI.Lower <- BS.Values[round(alpha1 * B)]
    }
    CI.Upper <- BS.Values[round(alpha2 * B)]	
  }
  
  # Include cocor for analytical results and Fisher Z test
  cc <- cocor.indep.groups(grp2_cor,grp1_cor,n2,n1,alpha = alpha,conf.level = CI.width,null.value = 0,
                           return.htest = T)
  
  result <- list("Call"=call,"Group.1"=grp1,"Group.2"=grp2,"Var1"=colnames(data)[1],"Var2"=colnames(data)[2],
                 "n1"=n1,"n2"=n2,"r1"=grp1_cor,"r2"=grp2_cor,"r.diff"=r_diff,"CI.boot"=c(CI.Lower,CI.Upper),
                 "Fisher.Z"=cc$fisher1925,"Zou2007"=cc$zou2007,method=method,"CI.level"=CI.width,
                 "n.boot"=B,"H0.val"=H0,"boot.vals"=BS.Values)
  class(result) <- "boot.cor"
  # 2nd generation p-value - if interval null hypothesis provided
  if(length(H0)>1){
    result$p.val.2g <- pval.2g(c(CI.Lower,CI.Upper),H0[1],H0[2])
    result$H0.low <- signif(H0[1],3)
    result$H0.high <- signif(H0[2],3)
    result$H0 <- NULL
  }
  
  return(result)
}

# print.boot.cor function
# Pretty-printing the above function
print.boot.cor <- function(x,digits=3){
  cat("Call: ")
  print(x$Call)
  cat("\n")
  cat(paste0("Bootstrapped correlation difference for ",x$Var1," and ",x$Var2,":"))
  cat(paste0("\n Group 1: ",x$Group.1,", r1 = ",round(x$r1,digits)))
  cat(paste0("\n Group 2: ",x$Group.2,", r2 = ",round(x$r2,digits)))
  cat("\n")
  cat("\n")

  cat(paste0("Correlation Difference (r2 - r1) = ",round(x$r.diff,digits)))
  cat(paste0("\n BCa bootstrapped (B = ",x$n.boot,") ",100*(x$CI.level),"% CI: [",
             paste(format(round(x$CI.boot[1],digits=digits),nsmall=1),
                   format(round(x$CI.boot[2],digits=digits),nsmall=1),sep=","),"]"))
  cat(paste0("\n Zou (2007) analytic ",100*(x$CI.level),"% CI: [",
             paste(format(round(x$Zou2007$conf.int[1],digits=digits),nsmall=1),
                   format(round(x$Zou2007$conf.int[2],digits=digits),nsmall=1),sep=","),"]"))
  cat(paste0("\n Fisher (1925) Z-test: Z = ",round(x$Fisher.Z$statistic,2),", p = ",
             signif(x$Fisher.Z$p.value,digits=3),"\n\n"))
  
  if(is.null(x$p.val.2g)){
    cat("Null Hypothesis: r2 - r1 =",x$H0.val,"\n")
    cat(" Bootstrap: ")
    if(x$CI.boot[1] < x$H0.val & x$H0.val < x$CI.boot[2]){
      cat(paste0("Null hypothesis retained (interval includes ",x$H0.val,")\n"))
    } else{
      cat(paste0("Null hypothesis rejected (interval does not include ",x$H0.val,")\n"))
    }
    cat(" Zou (2007): ")
    if(x$Zou2007$conf.int[1] < x$H0.val & x$H0.val < x$Zou2007$conf.int[2]){
      cat(paste0("Null hypothesis retained (interval includes ",x$H0.val,")\n"))
    } else{
      cat(paste0("Null hypothesis rejected (interval does not include ",x$H0.val,")\n"))
    }
    cat(" Fisher (1925): ")
    if(x$Fisher.Z$statistic > 1-x$CI.level){
      cat(paste0("Null hypothesis retained (p > ",1-x$CI.level,")\n"))
    } else{
      cat(paste0("Null hypothesis rejected (p < ",1-x$CI.level,")\n"))
    }
  } else {
    cat(paste0(" 2nd-generation p-value = ",x$p.val.2g," with interval H0 [",x$H0.low,", ",x$H0.high,"]"))
    if(x$p.val.2g==0){cat(" (reject interval null)")}
    else if(x$p.val.2g == 0){cat(" (accept equivalence to interval null)")}
    else {cat(" (inconclusive)")}    
  }
}

# Quick nice-looking Histogram with vector as only needed argument
quick.hist <- function(X,title="Histogram",x.lab="Value",y.lab="Density",density.fill="#FF6666",density.alpha=0.2,bins=30,
                       adjust=1){
  ggplot(data.frame(X),aes(x=X)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=bins) +
    geom_density(alpha=density.alpha,adjust=adjust, fill=density.fill) + xlab(x.lab) + ylab(y.lab) +
    ggtitle(title) + 
    theme(legend.position="right",
          plot.title = element_text(size=14, face="bold", hjust=0.5),
          legend.title = element_text(size=10, face="bold", hjust=0.5),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold"),
          strip.text.y = element_text(size = 20, colour = "white"),
          strip.background = element_rect(colour="darkgrey", fill="darkgrey"))
}

# Histogram for multiple groups
# Modified from https://stackoverflow.com/questions/6957549/overlaying-histograms-with-ggplot2-in-r
plot_multi_histogram <- function(df, feature, label_column,x.lab=NULL,y.lab="Density",
                                 bins=30,adjust=1,hist.alpha=0.7,density.alpha=0.35) {
  if(is.null(x.lab)) x.lab <- feature
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=hist.alpha, position="identity", aes(y = ..density..), color="black",bins=bins) +
    geom_density(alpha=density.alpha,adjust=adjust) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=x.lab, y = y.lab)
  plt + guides(fill=guide_legend(title=label_column))
}

# Get Sensitivity, Specificity, PPV, NPV from table
# Can now input prevalence and get adjusted PPV/NPV
dxVals <- function(tab,prev_adj=NULL){
  sens <- tab[1]/sum(tab[,1])
  spec <- tab[2,2]/sum(tab[,2])
  ppv <- tab[1]/sum(tab[1,])
  npv <- tab[2,2]/sum(tab[2,])
  posLR <- sens/(1-spec)
  negLR <- (1-sens)/spec
  if(!is.null(prev_adj)){
    ppv_adj <- (sens*prev_adj)/(sens*prev_adj + (1-spec)*(1-prev_adj))
    npv_adj <- (spec*(1-prev_adj))/((1-sens)*prev_adj + spec*(1-prev_adj))
    names(ppv_adj) <- paste0("PPV[",prev_adj,"]")
    names(npv_adj) <- paste0("NPV[",prev_adj,"]")
  } else{
    ppv_adj <- npv_adj <- NULL
  }
  
  return(c("Sensitivity"=sens,"Specificity"=spec,
           "PPV"=ppv,"NPV"=npv,"LR+"=posLR,"LR-"=negLR,ppv_adj,npv_adj))
}

# (Stratified) Percentile Bootstrap of sensitivity/specificity/PPV/NPV/LR
# Zack Williams
# 05/05/2020
# X = predictor (dichotomous)
# Y = outcome
sens.spec <- function(X,Y=NULL,B=10000,ci.width=0.95,rev=F,revpars=c(1,2),
                           stratified = T,seed=12345,prev_adj=NULL,digits=3){
  set.seed(seed)
  alpha <- 1-ci.width
  # If table, convert to matrix
  if(class(X)=="table" & length(dim(X))== 2 & all(dim(X)==2)){
    data <- rbind(cbind(rep(c("F"),X[1,1]),rep(c("F"),X[1,1])),
                  cbind(rep(c("F"),X[1,2]),rep(c("T"),X[1,2])),
                  cbind(rep(c("T"),X[2,1]),rep(c("F"),X[2,1])),
                  cbind(rep(c("T"),X[2,2]),rep(c("T"),X[2,2])))
    X <- data[,1]
    Y <- data[,2]
  }
  
  XY <- cbind(X,Y)
  X <- X[complete.cases(XY)]
  Y <- Y[complete.cases(XY)]
  XY <- XY[complete.cases(XY),]
  if(rev){
    tab.obs <- DescTools::Rev(table(X,Y),revpars)
  } else{
    tab.obs <- table(X,Y)
  }
  if(all(c(dim(tab.obs),length(dim(tab.obs)))!=2)){
    stop("ERROR: Both X and Y must have two levels.")
  }
  
  vals.obs <- dxVals(tab.obs,prev_adj=prev_adj)
  
  # If B==0 (no bootstrap), return sens/spec alone
  if(B==0){
    return(round(vals.obs,digits))
  }
  # Simple percentile bootstrap
  BS <- t(mcmapply(function(i){
    if(stratified){ # Stratified bootstrap (same prevalence of outcome [Y] each iteration)
      lvls <- unique(Y)
      y1_rows <- sample(which(Y==lvls[1]),replace=T)
      y2_rows <- sample(which(Y==lvls[2]),replace=T)
      b <- XY[c(y1_rows,y2_rows),] # Column marginals are the same as original table
    } else{ # Simple percentile bootstrap
      rows <- sample(nrow(XY),replace = T)
      b <- XY[rows,]
    }
    # Reverse table if rev=T (for putting TRUE/TRUE before FALSE/FALSE in logicals)
    if(rev){
      tab <- DescTools::Rev(table(b[,1],b[,2]),revpars)
    } else{
      tab <- table(b[,1],b[,2])
    }
    dxVals(tab,prev_adj=prev_adj)
  },1:B))
  
  result <- t(rbind("Est"=vals.obs,apply(BS,2,quantile,c(alpha/2,1-alpha/2))))

  return(round(result,digits))
}

# Method to return the optimal cut point of a ROC object (from pROC) based on the
# concordance probability method (maximize Se * Sp)
# Liu, X. (2012). Classification accuracy and cut point selection. Statistics in medicine, 31(23), 2676-2686.
coords.CP <- function(ROC){
  c <- coords(ROC,transpose=F)
  CPs <- apply(c,1,function(X){X[2]*X[3]})
  res <- c[which.max(CPs),]
  return(rbind("best.CP"=res))
}

# Method to return the optimal cut point of a ROC object (from pROC) based on the
# Index of Union (IU) method (minimize abs(Se - AUC) + abs(Sp - AUC))
# Unal, I. (2017). Defining an optimal cut-point value in ROC analysis: an alternative approach. 
# Computational and mathematical methods in medicine, 2017.
# https://doi.org/10.1155/2017/3762651
coords.IU <- function(ROC){
  a <- c(ROC$auc)
  c <- coords(ROC,transpose=F)
  IUs <- apply(c,1,function(X){abs(X[2] - a) + abs(X[3] - a)})
  res <- c[which.min(IUs),]
  return(rbind("best.IU"=res))
}

# Calculate I^2 from a list of study weights and tau^2
calcI2 <- function(tau2,wi,digits=2){
  k <- length(wi)
  H2 <- ((sum(wi) - sum(wi^2)/sum(wi)) * tau2 / (k-1) + 1)
  I2 <- (H2 - 1)/H2 * 100
  H <- sqrt(H2)
  return(round(c("I^2"=I2,"H"=H),digits))
}

### Inverse Normal Rank Transformation (INT)
# Popular in Behavioral Genetics
# Can run a normal LMEM on these data and it's basically a Wilcoxon
# cf. https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html
INT <- function(u,k=3/8){
  Ru <- rank(u,na.last = "keep")
  n <- sum(!is.na(Ru))
  return(qnorm((Ru-k)/(n + 1 - 2*k)))
}

cor_OR <- function(x,use="pairwise"){
  if(any(x > 1)){
    nb_rows <- apply(x,2,function(X){length(levels(factor(X)))>2})
    message(paste0("Removing the following row(s) with non-binary items: ",paste(which(nb_rows),collapse=", ")))
    x <- x[!nb_rows,!nb_rows]
  }
  or_mat <- mipfp::Corr2Odds(cor(x,use=use),marg.probs = apply(x,2,mean,na.rm=T))$odds
  diag(or_mat) <- NA
  return(or_mat)
}

# Little's MCAR Test (mcar)
# Zack Williams, 02/09/21
# Taken from: https://stats-bayes.com/post/2020/08/14/r-function-for-little-s-test-for-data-missing-completely-at-random/
mcar <- function(x){ 
  if(!require(norm)) {
    stop("You must have norm installed to use LittleMCAR") 
  } 
  
  # if(!require(data.table)) {
  #   stop("Please install the R-package data.table to use mcar")
  # }
  
  if(!(is.matrix(x) | is.data.frame(x))) {
    stop("Data should be a matrix or dataframe")
  }
  
  if (is.data.frame(x)){
    x <- data.matrix(x)
  }
  
  # delete rows of complete missingness
  foo <- function(x) return(any(!is.na(x)))
  dd <- apply(X = x, MARGIN = 1L, FUN = foo)
  dd <- which(!dd, arr.ind = TRUE)
  if(length(dd) > 0) 
    x <- x[-dd,]
  
  # define variables        
  n.var <- ncol(x) # number of variables
  n <- nrow(x)  #number of respondents
  var.names <- colnames(x)
  r <- 1 * is.na(x)
  
  nmis <- as.integer(apply(r, 2, sum))  #number of missing data for each variable REWRITE
  mdp <- (r %*% (2^((1:n.var - 1)))) + 1  #missing data patterns
  x.mp <- data.frame(cbind(x,mdp)) # add column indicating pattern
  colnames(x.mp) <- c(var.names,"MisPat") # set name of new column to MisPat
  n.mis.pat <- length(unique(x.mp$MisPat)) # number of missing data patterns
  p <- n.mis.pat-1 # number of Missing Data patterns minus 1 (complete data row)
  
  
  s <- prelim.norm(x)
  ll <- em.norm(s)
  fit <- getparam.norm(s = s, theta = ll)
  
  # gmean<-mlest(x)$muhat #ML estimate of grand mean (assumes Normal dist)
  gmean <- fit$mu
  # gcov<-mlest(x)$sigmahat #ML estimate of grand covariance (assumes Normal dist)
  gcov <- fit$sigma
  colnames(gcov) <- rownames(gcov) <- colnames(x)
  
  #recode MisPat variable to go from 1 through n.mis.pat
  x.mp$MisPat2 <- rep(NA,n)
  for (i in 1:n.mis.pat){ 
    x.mp$MisPat2[x.mp$MisPat == sort(unique(x.mp$MisPat), partial=(i))[i]]<- i 
  }
  
  x.mp$MisPat<-x.mp$MisPat2
  x.mp<-x.mp[ , -which(names(x.mp) %in% "MisPat2")]
  
  #make list of datasets for each pattern of missing data
  datasets <- list() 
  for (i in 1:n.mis.pat){
    datasets[[paste("DataSet",i,sep="")]]<-x.mp[which(x.mp$MisPat==i),1:n.var]
  }
  
  #degrees of freedom
  kj<-0
  for (i in 1:n.mis.pat){ 
    no.na<-as.matrix(1* !is.na(colSums(datasets[[i]]))) 
    kj<-kj+colSums(no.na) 
  }
  
  df<-kj -n.var
  
  #Little's chi-square
  d2<-0
  cat("this could take a while")
  
  # this crashes at the missingness pattern where every column is missing
  # this for-loop can be handled faster with plyr-function
  for (i in 1:n.mis.pat){ 
    mean <- (colMeans(datasets[[i]])-gmean) 
    mean <- mean[!is.na(mean)] 
    keep <- 1* !is.na(colSums(datasets[[i]])) 
    keep <- keep[which(keep[1:n.var]!=0)] 
    cov <- gcov 
    cov <- cov[which(rownames(cov) %in% names(keep)) , which(colnames(cov) %in% names(keep))] 
    d2 <- as.numeric(d2+(sum(x.mp$MisPat==i)*(t(mean)%*%solve(cov)%*%mean)))
  }
  
  #p-value for chi-square
  p.value<-1-pchisq(d2,df)
  
  #descriptives of missing data
  amount.missing <- matrix(nmis, 1, length(nmis))
  percent.missing <- amount.missing/n
  amount.missing <- rbind(amount.missing,percent.missing)
  colnames(amount.missing) <- var.names
  rownames(amount.missing) <- c("Number Missing", "Percent Missing")
  
  list(chi.square = d2, 
       df = df, 
       p.value = p.value, 
       missing.patterns = n.mis.pat, 
       amount.missing = amount.missing, 
       data = datasets)
}

### leftjoin_closest_val
## Zack Williams, 04/18/2021
# Useful when merging two dataframes where the second df may have more than one value for each case
# Allows you to pick a case to retain for each based on the row of origin in df1 and the difference between selected variables 
# max.difference - only add in a row if the difference in variables is < this value (e.g., don't merge two measures administered >2 years apart)
leftjoin_closest_val <- function(df1,df2,comp1,comp2,merge.by=NULL,max.difference=NULL){
  if(!is.data.frame(df1) | !is.data.frame(df2)){stop("Both 'df1' and 'df2' arguments should be data frames.")}
  if(!comp1 %in% names(df1)){stop("Unidentified column selected for 'comp1'—check the variable names in df1.")}
  if(!comp2 %in% names(df2)){stop("Unidentified column selected for 'comp2'—check the variable names in df2.")}
  df1$UNIQUE_VALUE <- paste0("UNIQ",1:nrow(df1)) # Assumes you want every row in df to be preserved
  
  # Identify columns that are repeated
  same_cols <- names(df1)[names(df1) %in% names(df2)]
  # If one of difference cols is repeated, cut the one in the other df
  if(comp1 %in% same_cols){
    message(paste0("Column '",comp1,"' is present in both dataframes. Variable excluded from df2 to preserve uniqueness."))
    df2[,comp1] <- NULL
  }
  if(comp2 %in% same_cols){
    message(paste0("Column '",comp2,"' is present in both dataframes. Variable excluded from df1 to preserve uniqueness."))
    df1[,comp2] <- NULL
  }
  same_cols <- same_cols[!same_cols %in% c(merge.by,comp1,comp2)]
  
  # Put df and df2 together (creating many duplicate rows)
  df_merged <- left_join(df1,df2,by=merge.by,suffix=c(".FIRSTCOL",".SECONDCOL"))
  # If comp1 or comp2 in both dataframes, change values of comp1
  df_merged$DIFF_SCORE_VAR <- abs(df_merged[,comp1] - df_merged[,comp2])
  df_merged <- df_merged[order(df_merged$DIFF_SCORE_VAR),]
  # Now cut duplicated rows (leaving only the one with the smallest value of DIFF_SCORE_VAR)
  df_merged <- df_merged[!duplicated(df_merged$UNIQUE_VALUE),]
  # If max.difference provided, remove values with DIFF_SCORE_VAR > max.difference
  if(is.numeric(max.difference)){
    df_merged[which(df_merged$DIFF_SCORE_VAR > max.difference),(ncol(df1)+1):ncol(df_merged)] <- NA
  }
  # Remove columns used to merge from same_cols (don't need to be transferred)
  same_cols <- same_cols[!(same_cols %in% names(df_merged))]
  # Now transfer over information from new columns
  for(c in same_cols){
    first.col <- paste0(c,".FIRSTCOL")
    second.col <- paste0(c,".SECONDCOL")
    df_merged[is.na(df_merged[,first.col]),first.col] <- df_merged[is.na(df_merged[,first.col]),second.col]
    df_merged[,second.col] <- NULL
  }
  names(df_merged) <- gsub("\\.FIRSTCOL$","",names(df_merged))
  df_merged$UNIQUE_VALUE <- df_merged$DIFF_SCORE_VAR <- NULL
  return(df_merged)
}