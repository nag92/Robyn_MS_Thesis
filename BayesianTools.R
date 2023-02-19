################################################################################
###### BayesianTools.R - a library of functions to do Bayesian estimation ######
################################################################################
### Zack Williams, created 02/23/2020, last updated 11/06/2022
# Most recent update: Added BBcor wrapper for Bayesian bootstrap correlations

# Package loading using pacman
if(!require(pacman)){
  install.packages("pacman")
  require(pacman)
}
if(!require(BBcor)){
  if(!require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("donaldRwilliams/BBcor")
}
if(!require(cmdstanr)){
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}
# if(!require(tidyverse)){
#   install.packages("tidyverse")
#   require(tidyverse)
# }

# Required packages: 
pacman::p_load(LaplacesDemon,brms,tidyverse,easystats,BayesFactor,parallel,pbapply,
               logspline,emmeans,formula.tools,mvnfast,coda,wCorr,
               tidybayes,DescTools,missForest,ggridges,mice,posterior,cmdstanr)

options("brms.backend" = "cmdstanr")

# Utility functions: rep.row, rep.col
rep.row<-function(x,n){ 
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# Summarize MCMC output (using posterior package) - no longer uses coda package
#### No longer 
MCMC_summary <- function(MCMClist,ci=0.95,ci.type="eti"){
  ci.lo <- (1-ci)/2
  ci.hi <- 1-ci.lo
  summ <- summary(MCMClist)
  rawlist <- suppressWarnings(MCMClist[,1:(ncol(MCMClist)-3)])
  CI <- t(apply(rawlist,2,function(temp){
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(temp,c(ci.lo,ci.hi))
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3])
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    return(CI)
  }))
  out <- cbind(summ[,c(2,4,3)],apply(rawlist,2,LaplacesDemon::Mode),
               CI,summ[,8:10])
  colnames(out)[3:4] <- c("Mdn","MAP")
  return(out)
}

### post_summary: Summarize a posterior distribution (Mean, SD, Mdn, MAP, 95% CrI, Pd, Bayes Factor vs H0, Bayes Factor vs. ROPE)
# Can handle posteriors that are transformed or not (H0/ROPE should also be specified on scale of the posterior)
# For instance, this would be specifying H0 = 0, ROPE = c(log(1/1.1),log(1.1)) for a log-OR with transform = exp
post_summ <- function(post,ci=0.95,ci.type="hdi",transform=NULL,prior=NULL,H0=0,ROPE=c(-0.1,0.1),logBF=TRUE,seed=12345,na.rm=TRUE){
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  set.seed(seed)
  if(is.null(transform)){ # No function specified to transform the posterior distribution
    M <- mean(post,na.rm=na.rm)
    SD <- sd(post,na.rm=na.rm)
    MAP <- LaplacesDemon::Mode(post)
    Mdn <- median(post,na.rm=na.rm)
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(post,c(ci.lo,ci.hi),na.rm=na.rm)
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- unname(unlist(bayestestR::hdi(post,ci=ci))[2:3])
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    pd <- mean(sign(post-H0)==sign(Mdn-H0),na.rm=na.rm) # Probability of direction
    p.exceeds.rope <- max(mean(post < ROPE[1],na.rm=na.rm),mean(post > ROPE[2],na.rm=na.rm))
    p.rope <- 1 - mean(post < ROPE[1],na.rm=na.rm) - mean(post > ROPE[2],na.rm=na.rm) # Probability posterior ES within ROPE
    # Calculate Bayes factors if prior is supplied
    if(!is.null(prior)){
      bf.rope <- exp(bf_rope(post,prior=prior,null = ROPE)$log_BF)
      bf.sd <- exp(bf_parameters(post,prior=prior,null=H0)$log_BF)
      # Log-transform  BF vals if desired
      if(logBF){
        bf.rope <- log(bf.rope)
        bf.sd <- log(bf.sd)
      }
    } else{
      bf.rope <- bf.sd <- NA
    }
  } else{ # transform is specified (should be a function)
    if(!any(class(transform) %in% c("standardGeneric","function"))){
      stop("To transform posterior, 'transform' argument must a function name (e.g., transform = exp for log Odds Ratio input).")
    }
    # Now calculate the same indices, but on untransformed scale, then transform
    M <- transform(mean(post,na.rm=na.rm)) # This is the geometric mean if transform = exp
    SD <- sd(transform(post,na.rm=na.rm))
    MAP <- transform(LaplacesDemon::Mode(post))
    Mdn <- transform(median(post,na.rm=na.rm))
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(post,c(ci.lo,ci.hi),na.rm=na.rm)
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- unname(unlist(bayestestR::hdi(post,ci=ci))[2:3])
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    CI[1] <- transform(CI[1])
    CI[2] <- transform(CI[2])
    # Calculate Bayes factors if prior is supplied
    if(!is.null(prior)){
      bf.rope <- exp(bf_rope(post,prior=prior,null = ROPE)$log_BF)
      bf.sd <- exp(bf_parameters(post,prior=prior,null=H0)$log_BF)
      # Log-transform  BF vals if desired
      if(logBF){
        bf.rope <- log(bf.rope)
        bf.sd <- log(bf.sd)
      }
    } else{
      bf.rope <- bf.sd <- NA
    }
    pd <- mean(sign(post-H0)==sign(median(post)-H0),na.rm=na.rm) # Probability of direction
    
    p.exceeds.rope <- max(mean(post < ROPE[1],na.rm=na.rm),mean(post > ROPE[2],na.rm=na.rm))
    p.rope <- 1 - mean(post < ROPE[1],na.rm=na.rm) - mean(post > ROPE[2],na.rm=na.rm) # Probability posterior ES within ROPE
    
    # Transforms the null values to the scale of the transformed parameter (tests done on untransformed values)
    H0 <- transform(H0)
    ROPE[1] <- transform(ROPE[1])
    ROPE[2] <- transform(ROPE[2])
  }
  
  # Summary of values for the standardized posterior distribution
  if(logBF){
    return(c("M"=M,"SD"=SD,"MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"logBF.10"=bf.sd,"H0"=H0,
             "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
             "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
             "logBF.ROPE"=bf.rope))
  } else{
    return(c("M"=M,"SD"=SD,"MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"BF.10"=bf.sd,"H0"=H0,
             "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
             "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
             "BF.ROPE"=bf.rope))
  }
}


# Jeffreys's (1961) BF Guidelines 
# > 100	Extreme evidence for H1
# 30 – 100	Very strong evidence for H1
# 10 – 30	Strong evidence for H1
# 3 – 10	Moderate evidence for H1
# 1 – 3	Anecdotal evidence for H1
# 1	No evidence
# 1/3 – 1	Anecdotal evidence for H1
# 1/3 – 1/10	Moderate evidence for H1
# 1/10 – 1/30	Strong evidence for H1
# 1/30 – 1/100	Very strong evidence for H1
# < 1/100	Extreme evidence for H1
bf_evidence <- function(BF,model=T){
  BF <- as.numeric(BF)
  if(model){
    as.character(cut(BF,c(-1,1/100,1/30,1/10,1/3,0.999,1.001,3,10,30,100,Inf),
                     c("Extreme Evidence Against Model","Very Strong Evidence Against Model",
                       "Strong Evidence Against Model",
                       "Moderate Evidence Against Model","Anecdotal Evidence Against Model","Reference Model",
                       "Anecdotal Evidence For Model","Moderate Evidence For","Strong Evidence For Model",
                       "Very Strong Evidence For Model","Extreme Evidence For Model")))
  } else{ # regarding null and alternative hypothesis
    as.character(cut(BF,c(-1,1/100,1/30,1/10,1/3,0.999,1.001,3,10,30,100,Inf),
                     c("Extreme Evidence For H0","Very Strong Evidence For H0",
                       "Strong Evidence For H0",
                       "Moderate Evidence For H0","Anecdotal Evidence For H0","Equal Evidence For H0 and H1",
                       "Anecdotal Evidence For H1","Moderate Evidence For H1","Strong Evidence For H1",
                       "Very Strong Evidence For H1","Extreme Evidence For H1")))
  }
  
}

### Summary for logistic regression brmsfit models
# Zack Williams, 08/15/2020
# Tests Odds Ratios against H0 of OR = 1 and ROPE of OR = c(0.833,1.2)
# Prior is a normal(0,1) distribution on logit scale
LR.summary <- function(BRM,ci=0.95,ci.type="hdi",ROPE=c(1/1.2,1.2),
                       priorfun=rnorm,seed=12345,digits=3,...){
  
  set.seed(seed)
  if(any(class(BRM)=="list")){
    if(any(class(BRM[[1]])=="brmsfit")){
      BRM <- pool.brm(BRM)
    } else{
      stop("ERROR: model must be of class 'brmsfit'")
    }
  }
  post <- as_draws_df(BRM,grep("^b_|^bs_|^bsp_",names(BRM$fit),value = T)) # Just get fixed effects
  n <- nrow(post)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  
  outlist <- mclapply(1:(ncol(post)-3),function(i){
    temp <- unlist(post[,i])
    MAP <- exp(LaplacesDemon::Mode(temp))
    Mdn <- exp(median(temp))
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(exp(temp),c(ci.lo,ci.hi))
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- exp(unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3]))
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    pd <- c(p_direction(temp)) # Probability of direction
    p.exceeds.rope <- max(mean(exp(temp) < ROPE[1]),mean(exp(temp) > ROPE[2]))
    p.rope <- 1 - mean(exp(temp) < ROPE[1]) - mean(exp(temp) > ROPE[2])
    bf.rope <- exp(bf_rope(exp(temp),prior=exp(priorfun(length(temp),...)),null = ROPE)$log_BF)
    bf.sd <- exp(bf_parameters(exp(temp),prior=exp(priorfun(length(temp),...)),null=1)$log_BF)
    return(c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"BF.0"=bf.sd,
             "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
             "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
             "BF.ROPE"=bf.rope))
  })
  names(outlist) <- gsub("b_","OR_",colnames(post)[1:(ncol(post)-3)])
  return(round(do.call(rbind,outlist),digits))
}

### Summary for negative binomial regression brmsfit models
# Zack Williams, 03/23/2021
# Tests Incidence Rate Ratios against H0 of IRR = 1 and ROPE of IRR = c(0.833,1.2)
# Prior is a normal(0,1) distribution on logit scale
NB.summary <- function(BRM,ci=0.95,ci.type="hdi",ROPE=c(1/1.2,1.2),
                       priorfun=rnorm,seed=12345,digits=3,...){
  
  set.seed(seed)
  if(any(class(BRM)=="list")){
    if(any(class(BRM[[1]])=="brmsfit")){
      BRM <- pool.brm(BRM)
    } else{
      stop("ERROR: model must be of class 'brmsfit'")
    }
  }
  post <- as_draws_df(BRM,grep("^b_|^bs_|^bsp_",names(BRM$fit),value = T)) # Just get fixed effects
  n <- nrow(post)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  
  outlist <- mclapply(1:(ncol(post)-3),function(i){
    temp <- unlist(post[,i])
    MAP <- exp(LaplacesDemon::Mode(temp))
    Mdn <- exp(median(temp))
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(exp(temp),c(ci.lo,ci.hi))
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- exp(unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3]))
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    pd <- c(p_direction(temp)) # Probability of direction
    p.exceeds.rope <- max(mean(exp(temp) < ROPE[1]),mean(exp(temp) > ROPE[2]))
    p.rope <- 1 - mean(exp(temp) < ROPE[1]) - mean(exp(temp) > ROPE[2])
    bf.rope <- exp(bf_rope(exp(temp),prior=exp(priorfun(length(temp),...)),null = ROPE)$log_BF)
    bf.sd <- exp(bf_parameters(exp(temp),prior=exp(priorfun(length(temp),...)),null=1)$log_BF)
    return(c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"BF.0"=bf.sd,
             "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
             "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
             "BF.ROPE"=bf.rope))
  })
  names(outlist) <- gsub("b_","IRR_",colnames(post)[1:(ncol(post)-3)])
  return(round(do.call(rbind,outlist),digits))
}

### Summary for linear regression brmsfit models
# Zack Williams, 08/15/2020
# UPDATED 03/23/2021 to now recognize LR/NBR and switch to those functions
# Tests betas against H0 of b = 0 and ROPE of b = c(-0.1,0.1)
# Prior is a normal(0,1) distribution on logit scale
# ... = arguments to prior function
regression.summary <- function(BRM,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                               priorfun=rnorm,seed=12345,digits=3,...){
  # Wrap LR.summary and NB.summary
  if(BRM$family$family %in% c("bernoulli","cumulative") & BRM$family$link == "logit"){
    if(all(ROPE %in% c(-0.1,0.1))){ROPE <- c(1/1.2,1.2)}
    return(LR.summary(BRM,ci=ci,ci.type=ci.type,ROPE=ROPE,priorfun=priorfun,seed=seed,digits=digits,...))
  } else if(BRM$family$family == "negbinomial"){
    if(all(ROPE %in% c(-0.1,0.1))){ROPE <- c(1/1.2,1.2)}
    return(NB.summary(BRM,ci=ci,ci.type=ci.type,ROPE=ROPE,priorfun=priorfun,seed=seed,digits=digits,...))
  }
  set.seed(seed)
  if(any(class(BRM)=="list")){
    if(any(class(BRM[[1]])=="brmsfit")){
      BRM <- pool.brm(BRM)
    } else{
      stop("ERROR: model must be of class 'brmsfit'")
    }
  }
  post <- as_draws_df(BRM,grep("^b_|^bs_|^bsp_",names(BRM$fit),value = T)) # Just get fixed effects
  n <- nrow(post)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  
  outlist <- mclapply(1:(ncol(post)-3),function(i){
    temp <- unlist(post[,i])
    MAP <- LaplacesDemon::Mode(temp)
    Mdn <- median(temp)
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(temp,c(ci.lo,ci.hi))
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3])
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    pd <- c(p_direction(temp)) # Probability of direction
    p.exceeds.rope <- max(mean(temp < ROPE[1]),mean(temp > ROPE[2]))
    p.rope <- 1 - mean(temp < ROPE[1]) - mean(temp > ROPE[2])
    bf.rope <- exp(bf_rope(temp,prior=priorfun(length(temp),...),null = ROPE)$log_BF)
    bf.sd <- exp(bf_parameters(temp,prior=priorfun(length(temp),...),null=0)$log_BF)
    return(c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"BF.0"=bf.sd,
             "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
             "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
             "BF.ROPE"=bf.rope))
  })
  names(outlist) <- gsub("b_","",colnames(post)[1:(ncol(post)-3)])
  return(round(do.call(rbind,outlist),digits))
}

### Bayesian Power Analysis for a t-test
## Zack Williams, 02.23.2020
# Goal: find power to, given a certain effect and sample size, find a "significant" Bayes Factor >3 or <1/3
# Uses the Bayes factor vs. ROPE calculated by BayesTestR (among others that can be found in the BEST package)
# Based on blog posts at: https://solomonkurz.netlify.com/post/bayesian-power-analysis-part-i/ 
# and https://vuorre.netlify.com/post/2017/01/02/how-to-compare-two-groups-with-robust-bayesian-estimation-using-r-stan-and-brms/#robust-bayesian-estimation
## Note that default priors are not quite what Kruschke does (which is sort of an empirical bayes type thing)
## Default priors: N(0,1) on intercept/effect size
# Now allows for unequal sample sizes (n1 and n2); only n1 specified = equal sample sizes
pwr.BEST.ttest <- function(n1,d=0.5,s1=1,n.sim=200,n2=n1,s2=s1,ROPE=c(-0.1,0.1),ci=0.95,BF.sig=3,
                           Prior=NULL,H0.test=F,dig=3,seed=1234){
  if(H0.test){d <- 0}
  ## Function to simulate (normal) data with chosen Effect Size
  sim_d <- function(seed,n1,n2,s1,s2,es=d) {
    mu_t <- es * sqrt(mean(c(s1^2,s2^2)))
    mu_c <- 0
    set.seed(seed)
    tibble(group = c(rep("control",n1),rep("treatment",n2))) %>% 
      mutate(treatment = ifelse(group == "control", 0, 1),
             y         = ifelse(group == "control", 
                                rnorm(n1, mean = mu_c, sd = s1),
                                rnorm(n2, mean = mu_t, sd = s2)) %>% 
               scale() %>% c() # now scales outcome vars to accommodate different SDs
             )
  }
  
  dat <- sim_d(seed,n1,n2,s1,s2,es=d)
  
  if(is.null(Prior)){
    pr <- c(prior(normal(0, 1), class = Intercept),prior(normal(0,1), class = b))
  } else{
    pr <- Prior
  }
  message("Compiling initial Stan model...")
  mod_robust <- suppressMessages(brm(
    bf(y ~ group, sigma ~ group),
    prior=pr,
    family=student,
    data = dat, 
    cores=detectCores()-1,
    seed=seed,
    silent=T,refresh=0,iter=3000,warmup = 1000
  ))
  message(" Done!\n")
  
  sims <- pbsapply(1:n.sim,function(it){
    dat_s <- sim_d(seed + it,n1,n2,s1,s2,es=d)
    mod_s <- suppressMessages(update(mod_robust,newdata=dat_s,seed=seed+it,silent=T,refresh=0))
    mc <- do.call(rbind,as.mcmc(mod_s))
    sd1 <- exp(mc[,2])
    sd2 <- exp(mc[,2] + mc[,4])
    sdpooled <- sqrt(rowMeans(cbind(sd1,sd2)^2))
    effSize <- mc[,3]/sdpooled
    return(unlist(c(exp(bf_rope(effSize,prior = rnorm(length(effSize),0,1),null=ROPE)$log_BF),
             "SD"=exp(bayesfactor_parameters(effSize,prior = rnorm(length(effSize),0,1),null=0)$log_BF),
             "Pd"=p_direction(effSize),quantile(effSize,c((1-ci)/2,1-(1-ci)/2)))))
  },cl=detectCores()-1)
  
  BF.ROPE <- sims[1,]
  BF.SD <- sims[2,]
  Pd <- sims[3,]
  ci.lo <- sims[4,]
  ci.hi <- sims[5,]
  if(H0.test){
    return(rbind("Value"=round(c("Pwr.BF.ROPE.H0"=mean(BF.ROPE < 1/BF.sig),"Pwr.BF.SD.H0"=mean(BF.SD < 1/BF.sig),
             "Pwr.CrI.ROPE.H0"=mean(ci.lo < ROPE[1] & ci.hi < ROPE[2]),
           "ROPE.l"=ROPE[1],"ROPE.u"=ROPE[2],"ES"=d,"N1"=n1,"N2"=n2,"SD1"=s1,"SD2"=s2,
           "BF.ROPE.Mdn"=median(BF.ROPE),"BF.SD.Mdn"=median(BF.SD)),dig)))
  } else{
    return(rbind("Value"=round(c("Pwr.BF.ROPE"=mean(BF.ROPE > BF.sig),"Pwr.BF.0"=mean(BF.SD > BF.sig),
           "Pwr.CrI.ROPE"=mean(ci.lo > ROPE[2]),"Pwr.CrI.0"=mean(Pd > 1-((1-ci)/2)),
           "ROPE.l"=ROPE[1],"ROPE.u"=ROPE[2],"ES"=d,"N1"=n1,"N2"=n2,"SD1"=s1,"SD2"=s2,
           "BF.ROPE.Mdn"=median(BF.ROPE),"BF.SD.Mdn"=median(BF.SD),"Pd.Mdn"=median(Pd)),dig)))
  }
}

### BEST.ttest.brm - use a similar procedure to BEST to conduct a robust Bayesian t-test
### Zack Williams, 08/08/2020
### Scales the DV to M = 0, SD = 1 and puts a prior on the SMD between groups
## Default priors include:
# N(0,1) for b (mean difference), intercept, and log(sigma [estimated separately for each grp])
# gamma (2,0.1) for nu
# Prior on ES (cohen's d) similar to a t(2,0,0.5) distribution, 95% CrI within ±2.6
# Tests ES (cohen's d) against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
BEST.ttest.brm <- function(form,data=NULL,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                           chains=5,iter=10000,warmup = 2000,seed=12345,cores=detectCores()-1,
                           file=NULL,overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  if(length(form) != 3 || length(all.vars(form)) != 2 || op(form) != "~" || is.null(lhs(form)) || is.null(rhs(form))){
    stop("Formula is incorrectly specified! Please enter formula of the type 'y ~ x'.")
  }
  if(!is.data.frame(data)){stop("Variable 'data' is empty and/or not a data frame.")}
  y <- lhs(form)
  x <- rhs(form)
  if(!deparse(x) %in% names(data)){stop(paste0("Variable '",deparse(x),"' not found! Undefined columns selected."))}
  if(!deparse(y) %in% names(data)){stop(paste0("Variable '",deparse(y),"' not found! Undefined columns selected."))}
  # Get M and SD (for rescaling results)
  My <- mean(data[,deparse(y)],na.rm=T)
  SDy <- sd(data[,deparse(y)],na.rm=T)
  # make sure x is a factor
  data[,deparse(x)] <- as.factor(data[,deparse(x)])
  lvls_x <- levels(data[,deparse(x)])
  if(length(lvls_x) != 2){stop("Grouping variable does not have only two levels! Check your data/formula.")}
  ns <- rowSums(table(data[,deparse(x)],data[,deparse(y)]))
  
  # Construct formula
  bform <- bf(paste0("scale(",y,") ~ ",x),paste0("sigma ~ 0 +",x))
  
  # Use default priors unless otherwise specified (N(0,1) for ES, intercept, slope between sigmas, and log(sigma))
  if(is.null(Prior)){
    pr <- c(prior(normal(0,1), class = b),
            prior(normal(0,1), class = Intercept),
            prior(normal(0,1), class = b, dpar="sigma"))
  } else{
    pr <- Prior
  }
  # Actually fit brm model
  mod_robust <- suppressWarnings(brm(
    bform,
    prior=pr,
    family=student,
    data = data, 
    cores=cores,
    chains=chains,
    seed=seed,sample_prior = T,
    iter=iter,warmup = warmup,
    ... # additional arguments can go here
  ))
  
  post <- as_draws_df(mod_robust)
  # Now modify the posterior to include the group parameters values/CrIs:
  sigma1 <- suppressWarnings(unlist(exp(post[,grep("^b_sigma_",colnames(post))[1]])))
  sigma2 <- suppressWarnings(unlist(exp(post[,grep("^b_sigma_",colnames(post))[2]])))
  
  post$Mu.GROUP1 <- suppressWarnings(unlist(post[,1])) * SDy + My
  post$Mu.GROUP2 <- suppressWarnings(rowSums(post[,c(1,2)])) * SDy + My
  
  post$Sigma.GROUP1 <- sigma1 * SDy # On scale of original data
  post$Sigma.GROUP2 <- sigma2 * SDy # On scale of original data
  post$Nu <- suppressWarnings(unlist(post[,5]))
  post$Mu.Diff <- suppressWarnings(unlist(post[,2])) * SDy
  post$Cohen.D <- temp <- suppressWarnings(unlist(post[,2])/sqrt((sigma1^2 + sigma2^2)/2))
  
  # Look at effect size samples, derive useful Bayesian indices of effect significance/existence
  prior_cols <- suppressWarnings(post[,grep("^prior_",colnames(post))])
  prior_sig1 <- exp(prior_cols[,grep("b_sigma$",colnames(prior_cols))])
  prior_sig2 <- sample(prior_sig1)
  es_prior <- unlist(prior_cols[,grep("prior_b$",colnames(prior_cols))]/sqrt(rowMeans(cbind(prior_sig1^2,prior_sig2^2))))
  
  post <- post[,c(grep("^Mu.GROUP1$",names(post)):grep("^Cohen.D$",names(post)),grep("^\\.chain|^\\.iteration|^\\.draw",names(post)))]
  
  colnames(post) <- gsub("GROUP1",lvls_x[1],colnames(post))
  colnames(post) <- gsub("GROUP2",lvls_x[2],colnames(post))
  
  post_summ <- MCMC_summary(post,ci=ci,ci.type = ci.type)
  
  set.seed(seed)
  n <- length(temp)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  MAP <- LaplacesDemon::Mode(temp)
  Mdn <- median(temp)
  if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
    CI <- quantile(temp,c(ci.lo,ci.hi))
    names(CI) <- paste0("Q.",names(CI))
  } else{ # HDI
    CI <- unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3])
    names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
  }
  pd <- c(p_direction(temp)) # Probability of direction
  p.exceeds.rope <- max(mean(temp < ROPE[1]),mean(temp > ROPE[2]))
  pmap <- c(p_map(temp))
  p.rope <- 1 - mean(temp < ROPE[1]) - mean(temp > ROPE[2]) # Probability posterior ES within ROPE
  bf.rope <- exp(bf_rope(temp,prior=es_prior,null = ROPE)$log_BF)
  bf.sd <- exp(bf_parameters(temp,prior=es_prior)$log_BF)
  
  # Summary of values for the standardized effect size parameter
  ES.summary <- c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"P.MAP"=pmap,"BF.0"=bf.sd,
           "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
           "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
           "BF.ROPE"=bf.rope)
  
  # Now make the function output
  result <- list("Summary"=post_summ,"ES"=ES.summary,"BRM"=mod_robust,"posterior.samples"=post,
                 "ROPE"=c("ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2]),
                 "BF.10"=bf.sd,"BF.01"=1/bf.sd,"BF.ROPE"=bf.rope,"BF.inROPE"=1/bf.rope,
                 "ci"=ci,"ci.type"=ifelse(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile"),"eti","hdi"),
                 "samples"=n,"N"=ns,"Formula"=form,"Call"=match.call())
  class(result) <- "BEST.ttest"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Pretty print function for Bayesian t-test
print.BEST.ttest <- function(x,verbose=F,dig=3){
  cat("Call:")
  print(x$Call)
  cat(paste0("\nBayesian two-sample t-test (BEST; Kruschke, 2013) with Region of Practical Equivalence (ROPE) d = [",paste(x$ROPE,collapse=", "),"]:"))
  cat(paste0("\n data: group ",names(x$N)[1]," (n = ",x$N[1],") and group ",names(x$N)[2]," (n = ",x$N[2],")"))
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  if(verbose){
    cat(paste0("\n\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):\n\n"))
    print(round(as.matrix(x$Summary[1:6,]),dig))
  } else{
    cat("\n")
  }
  cat(paste0("\n Effect Size (Standardized Mean Difference; ",names(x$N)[2]," — ",names(x$N)[1],"):"))
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$ES[2],dig),
             " [",round(x$ES[3],dig),", ",round(x$ES[4],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$ES[5]>0.999,"> 0.999",paste0("= ",signif(x$ES[5],dig))),
             " (",ifelse(x$ES[5]>0.999,">99.9",signif(x$ES[5]*100,dig)),
             "% chance that d ",ifelse(sign(x$ES[2])==1,">","<")," 0)"))
  cat(paste0("\n  ",ifelse(x$ES[10]<0.001,"<0.1",signif(x$ES[10]*100,dig)),"% of posterior density inside [",
             paste(x$ROPE,collapse=", "),"] (",ifelse(x$ES[11]>0.999,">99.9",signif(x$ES[11]*100,dig)),"% exceeds ROPE)"))
  cat(paste0("\n\n Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$BF.10>999 | x$BF.10 < 0.001,
                    format(x$BF.10,dig=3,scientific=T),
                    signif(x$BF.10,dig=3)),
             " (BF01 = ",ifelse(x$BF.01>999 | x$BF.01 < 0.001,
                                format(x$BF.01,dig=3,scientific=T),
                                signif(x$BF.01,dig=3)),") — ",
             bf_evidence(x$BF.10,model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$BF.ROPE>999 | x$BF.ROPE < 0.001,
                                                           format(x$BF.ROPE,dig=3,scientific=T),
                                                           signif(x$BF.ROPE,dig=3)),
             " (BF.inROPE = ",ifelse(x$BF.inROPE>999 | x$BF.inROPE < 0.001,
                                     format(x$BF.inROPE,dig=3,scientific=T),
                                     signif(x$BF.inROPE,dig=3)),") — ",
             bf_evidence(x$BF.ROPE,model = F))) # Bayes Factor Qualitative Interpretation
}


### brm_mi - a "fixed" version of brm_multiple that works with cmdstanr (default backend) and also calculates Rhat individually on every model to examine convergence
# Zack Williams, 11/06/2022 - updated 02/04/2023
# data should be a list of completed dataframes (output of multiple imputation)
# Also includes "better" calculation of MCMC diagnostics for multiply-imputed data using posterior package Rhat and bulk/tail ESS metrics than the brm_multiple defaults
brm_mi <- function(formula,data=NULL,family=gaussian(),prior=NULL,...,seed=1234,file=NULL,backend="cmdstanr",save_all_pars=TRUE,overwrite.file=FALSE,cores=detectCores()-1){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  
  # Error Catching for data input
  if(is.null(data) | !is.list(data) | !is.data.frame(data[[1]]) | !all(dim(data[[1]] == dim(data[[2]])))){
    stop("ERROR: 'data' must be a list of dataframes with identical columns (i.e., the output of multiple imputation).")
  }
  
  ### Meat of the function: Iteratively fit the function using the non-multiple brm function
  MOD_LIST <- lapply(1:length(data),function(imp){
    message(paste0("***Fitting model ",imp," of ",length(data),":***")) # Nice little message to say where you are in the fitting (I've always wanted this in brm_multiple)
    MOD <- brm(formula=formula,data=data[[imp]],family=family,prior=prior,seed=seed+imp,backend=backend,cores=cores,save_pars=brms::save_pars(all = save_all_pars),...,)
  })
  
  result <- combine_models(mlist=MOD_LIST,check_data = FALSE) # result - brmsfit object created from combining MOD_LIST together
  
  # Calculate Rhats and ESS from each imputed dataset using the posterior package (includes only fixed effects, random effects, distributional parameters)
  draws_list <- lapply(MOD_LIST,as_draws_rvars,variable = c("^b_","^bs_","^bsp_","^sd_","^simo_","^sigma","^alpha","^cor_","^disc_","^shape","^nu","^phi",
                                                            "^kappa","^beta","^zi","^hu","^zoi","^coi","^ndt","^bias","^xi","^quantile","^Intercept"),regex=T)
  message("Calculating marginal likelihoods...")
  set.seed(seed) # Set seed for reproducibility
  log_marginal_lik_list <- suppressMessages(pblapply(MOD_LIST,function(X){brms::bridge_sampler(X,cores=cores,silent=TRUE)$logml}))
  if(length(log_marginal_lik_list)>1){names(log_marginal_lik_list) <- 1:length(log_marginal_lik_list)}
  message("Done!")
  # Calculate bulk/tail ESS using posterior package
  ess_list <- lapply(draws_list,function(X){
    lapply(X,function(draws){
      c("ESS_Bulk"=mean(posterior::ess_bulk(draws)),"ESS_Tail"=mean(posterior::ess_tail(draws)))
    })
  })
  # Calculate rhat using posterior package
  rhat_list <- sapply(1:length(draws_list),function(i){
    X <- draws_list[[i]]
    rhats <- lapply(mclapply(X,posterior::rhat,mc.cores=detectCores()-1),mean)
    if(max(unlist(rhats)) > 1.05){
      warning("Warning: Imputation ",i," did not converge (max rhat = ",round(max(unlist(rhats)),2),"). More MCMC iterations may be needed.\n")
    }
    return(unlist(rhats))
  })
  # Combine ESS and Rhats together into one dataframe
  total_ess <- t(rowSums(array(unlist(ess_list),dim=c(length(ess_list[[1]][[1]]),length(ess_list[[1]]),length(ess_list))),dims = 2))
  rownames(total_ess) <- names(ess_list[[1]])
  colnames(total_ess) <- c("ESS_Bulk","ESS_Tail")
  
  # Throw warnings if bulk/tail ESS is low (<1000 in either case)
  if(min(total_ess[,1]) < 1000){warning("Warning: Smallest Bulk ESS is",round(min(total_ess[,1])),"). More MCMC iterations may be needed.\n")}
  if(min(total_ess[,2]) < 1000){warning("Warning: Smallest Tail ESS is",round(min(total_ess[,2])),"). More MCMC iterations may be needed.\n")}
  
  result$logml <- unlist(log_marginal_lik_list)
  result$MCMC_diagnostics <- cbind(total_ess,"Rhat_max"=apply(rhat_list,1,max)) # Add this to the model output
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

### BEST.ordtest.brm - use a similar procedure to Liddell & Kruschke, 2018
### Zack Williams, 08/28/2020
# Based on tutorial paper: Bürkner, P. C., & Vuorre, M. (2019). Ordinal regression models in psychology: A tutorial. 
# Advances in Methods and Practices in Psychological Science, 2(1), 77-101. https://doi.org/10.1177/2515245918823199
### Ordinal DV is turned into a latent normal var (sd = 1 for group 1, sd variable for grp 2)
### Prior on the MD between groups (in normal distribution units, i.e., d)
## Default priors include:
# N(0,1) for b (mean diff) and log(1/sd) [group 2 only], student_t(3, 0, 2.5) for intercepts
# sd [group 1] fixed at 1 to set scale of latent variable (as per Bürkner et al. paper)
# ES prior is similar to t(5.4,0,0.8) - ~95% of density between -2 and 2
# Tests ES (cohen's d) against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
BEST.ordtest.brm <- function(form,data=NULL,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                           chains=5,iter=10000,warmup = 2000,seed=12345,inits="0",family=cumulative("probit"),
                           cores=detectCores()-1,file=NULL,overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  if(length(form) != 3 || length(all.vars(form)) != 2 || op(form) != "~" || is.null(lhs(form)) || is.null(rhs(form))){
    stop("Formula is incorrectly specified! Please enter formula of the type 'y ~ x'.")
  }
  if(!is.data.frame(data)){stop("Variable 'data' is empty and/or not a data frame.")}
  y <- lhs(form)
  x <- rhs(form)
  if(!deparse(x) %in% names(data)){stop(paste0("Variable '",deparse(x),"' not found! Undefined columns selected."))}
  if(!deparse(y) %in% names(data)){stop(paste0("Variable '",deparse(y),"' not found! Undefined columns selected."))}
  # Convert y variable to numeric to make sure levels are correct
  data[,deparse(y)] <- ordered(data[,deparse(y)])
  
  lvls_x <- levels(as.factor(data[,deparse(x)]))
  ns <- rowSums(table(data[,deparse(x)],data[,deparse(y)]))
  
  if(length(levels(as.factor(data[,deparse(x)]))) != 2){stop("Grouping variable does not have only two levels! Check your data/formula.")}
  # Construct formula
  bform <- bf(paste0(y," ~ ",x),paste0("disc ~ 0 +",x))
  # Make it so that cmc=F (lf() method shown in https://journals.sagepub.com/doi/full/10.1177/2515245918823199 is broken)
  attr(bform$pforms$disc,"cmc") <- FALSE
  
  # Use default priors unless otherwise specified (N(0,1) for ES, b_disc, and log(disc intercept), t(3,0,2.5) for intercepts)
  if(is.null(Prior)){
    pr <- c(prior(normal(0,1), class = b),
            prior(normal(0,1), class = b,dpar="disc"),
            prior(student_t(3,0,2.5), class = Intercept)) # default prior for intercepts, but spelled out here anyways
  } else{
    pr <- Prior
  }
  # Actually fit brm model
  mod_robust <- suppressWarnings(brm(
    bform,
    prior=pr,
    family=family,
    data = data, 
    cores=cores,
    chains=chains,
    inits = inits,
    seed=seed,sample_prior = T,
    iter=iter,warmup = warmup,
    ... # additional arguments can go here
  ))
  

  post <- as_draws_df(mod_robust)
  # Now modify the posterior to include the group parameters values/CrIs:
  sigma1 <- 1
  sigma2 <- suppressWarnings(unlist(1/exp(post[,grep("^b_disc_",colnames(post))])))
  
  post$Mu.Diff <- suppressWarnings(unlist(post[,grep(paste0("^b_",x),colnames(post))]))
  post$Sigma.GROUP1 <- sigma1 # Fixed to 1
  post$Sigma.GROUP2 <- sigma2 # Allowed to vary
  temp <- post$Cohen.D <- suppressWarnings(unlist(post$Mu.Diff)/sqrt((sigma1^2 + sigma2^2)/2))
  
  # Look at effect size samples, derive useful Bayesian indices of effect significance/existence
  prior_cols <- suppressWarnings(post[,grep("^prior_",colnames(post))])
  prior_sig1 <- 1
  prior_sig2 <- 1/exp(prior_cols[,grep("b_disc$",colnames(prior_cols))])
  es_prior <- unlist(prior_cols[,grep("prior_b$",colnames(prior_cols))]/sqrt(rowMeans(cbind(prior_sig1^2,prior_sig2^2))))
  
  post <- post[,c(grep("^Mu.Diff$",names(post)):grep("^Cohen.D$",names(post)),grep("^\\.chain|^\\.iteration|^\\.draw",names(post)))]
  
  colnames(post) <- gsub("GROUP1",lvls_x[1],colnames(post))
  colnames(post) <- gsub("GROUP2",lvls_x[2],colnames(post))
  
  post_summ <- MCMC_summary(post,ci=ci,ci.type = ci.type)
  
  set.seed(seed)
  n <- length(temp)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  MAP <- LaplacesDemon::Mode(temp)
  Mdn <- median(temp)
  if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
    CI <- quantile(temp,c(ci.lo,ci.hi))
    names(CI) <- paste0("Q.",names(CI))
  } else{ # HDI
    CI <- unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3])
    names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
  }
  pd <- c(p_direction(temp)) # Probability of direction
  p.exceeds.rope <- max(mean(temp < ROPE[1]),mean(temp > ROPE[2]))
  pmap <- c(p_map(temp))
  p.rope <- 1 - mean(temp < ROPE[1]) - mean(temp > ROPE[2]) # Probability posterior ES within ROPE
  bf.rope <- exp(bf_rope(temp,prior=es_prior,null = ROPE)$log_BF)
  bf.sd <- exp(bf_parameters(temp,prior=es_prior)$log_BF)
  
  # Summary of values for the standardized effect size parameter
  ES.summary <- c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"P.MAP"=pmap,"BF.0"=bf.sd,
                  "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
                  "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
                  "BF.ROPE"=bf.rope)
  
  # Now make the function output
  result <- list("Summary"=post_summ,"ES"=ES.summary,"BRM"=mod_robust,"posterior.samples"=post,
                 "prior.samples"=cbind(prior_cols,"prior_Cohen.D"=es_prior),
                 "ROPE"=c("ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2]),
                 "BF.10"=bf.sd,"BF.01"=1/bf.sd,"BF.ROPE"=bf.rope,"BF.inROPE"=1/bf.rope,
                 "ci"=ci,"ci.type"=ifelse(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile"),"eti","hdi"),
                 "samples"=n,"N"=ns,"Formula"=form,"Call"=match.call())
  class(result) <- "BEST.ordtest"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Pretty print function for Bayesian ordinal test
print.BEST.ordtest <- function(x,verbose=T,dig=3){
  cat("Call:")
  print(x$Call)
  cat(paste0("\nBayesian two-sample mean difference test for ordinal variables with Region of Practical Equivalence (ROPE) d = [",paste(x$ROPE,collapse=", "),"]:"))
  cat(paste0("\n data: group ",names(x$N)[1]," (n = ",x$N[1],") and group ",names(x$N)[2]," (n = ",x$N[2],")"))
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  if(verbose){
    cat(paste0("\n\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):\n\n"))
    print(round(x$Summary[1:4,],dig))
  } else{
    cat("\n")
  }
  cat(paste0("\n Effect Size (Standardized Mean Difference; ",names(x$N)[2]," — ",names(x$N)[1],"):"))
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$ES[2],dig),
             " [",round(x$ES[3],dig),", ",round(x$ES[4],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$ES[5]>0.999,"> 0.999",paste0("= ",signif(x$ES[5],dig))),
             " (",ifelse(x$ES[5]>0.999,">99.9",signif(x$ES[5]*100,dig)),
             "% chance that d ",ifelse(sign(x$ES[2])==1,">","<")," 0)"))
  cat(paste0("\n  ",ifelse(x$ES[10]<0.001,"<0.1",signif(x$ES[10]*100,dig)),"% of posterior density inside [",
             paste(x$ROPE,collapse=", "),"] (",ifelse(x$ES[11]>0.999,">99.9",signif(x$ES[11]*100,dig)),"% exceeds ROPE)"))
  cat(paste0("\n\n Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$BF.10>999 | x$BF.10 < 0.001,
                    format(x$BF.10,dig=3,scientific=T),
                    signif(x$BF.10,dig=3)),
             " (BF01 = ",ifelse(x$BF.01>999 | x$BF.01 < 0.001,
                                format(x$BF.01,dig=3,scientific=T),
                                signif(x$BF.01,dig=3)),") — ",
             bf_evidence(x$BF.10,model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$BF.ROPE>999 | x$BF.ROPE < 0.001,
                                                           format(x$BF.ROPE,dig=3,scientific=T),
                                                           signif(x$BF.ROPE,dig=3)),
             " (BF.inROPE = ",ifelse(x$BF.inROPE>999 | x$BF.inROPE < 0.001,
                                     format(x$BF.inROPE,dig=3,scientific=T),
                                     signif(x$BF.inROPE,dig=3)),") — ",
             bf_evidence(x$BF.ROPE,model = F))) # Bayes Factor Qualitative Interpretation
}

### bayes.cor.brm - use a similar procedure to the baysianFirstAid bayes.cor.test, but in BRMS
### Based on https://solomonkurz.netlify.app/post/bayesian-correlations-let-s-talk-options/
### Zack Williams, 08/08/2020
### Does not scale variables, estimates both intercepts and sigma parameters 
# (i.e., assumes population mean/SD are not the same as sample mean/SD for tested variables)
## Default priors include:
# lkj(2) for correlation (more weight around zero with relatively heavy tails)
# Student's t(3,0,MAD) for log(sigma) [BRMS default]
# Student's t(3,Mdn,MAD) for intercepts [BRMS default]
# gamma (2,0.1) for nu (same for both X and Y) [BRMS default]
# Tests r against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
bayes.cor.brm <- function(x,y,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                           chains=5,iter=5000,warmup = 1000,seed=12345,cores=detectCores()-1,file=NULL,
                          overwrite.file=FALSE,...){
  xy <- tibble("x"=x,"y"=y)
  xy <- xy[complete.cases(xy),]
  N <- nrow(xy)
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  if(is.null(Prior)){
    pr <- prior(lkj(2), class = rescor)
  } else{
    pr <- Prior
  }
  
  cor_brm <- brm(bf(mvbind(x, y) ~ 1) + set_rescor(TRUE),
                 family=student,
                 data=xy,
                 prior = pr,
                 sample_prior = T,
                 iter = iter, warmup = warmup, 
                 chains = chains, cores = cores, 
                 seed = seed,
                 ... # additional arguments
                 )
  # Get posterior samples of model parameters
  post <- as_draws_df(cor_brm)
  model_pars <- MCMC_summary(post,ci=ci,ci.type = ci.type)
  model_pars <- model_pars[c(1:5,8),]
  row.names(model_pars) <- c("Mu.X","Mu.Y","Sigma.X","Sigma.Y","Nu","Rho.XY")
  # Get prior/posterior correlation distributions
  post_cor <- suppressWarnings(unlist(post[,8]))
  prior_cor <- suppressWarnings(unlist(post[,13]))
  nsamp <- length(post_cor)
  set.seed(seed)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  MAP <- LaplacesDemon::Mode(post_cor)
  Mdn <- median(post_cor)
  if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
    CI <- quantile(post_cor,c(ci.lo,ci.hi))
    names(CI) <- paste0("Q.",names(CI))
  } else{ # HDI
    CI <- unname(unlist(bayestestR::hdi(post_cor,ci=ci))[2:3])
    names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
  }
  pd <- c(p_direction(post_cor)) # Probability of direction
  p.exceeds.rope <- max(mean(post_cor < ROPE[1]),mean(post_cor > ROPE[2]))
  pmap <- c(p_map(post_cor))
  p.rope <- 1 - mean(post_cor < ROPE[1]) - mean(post_cor > ROPE[2]) # Probability posterior ES within ROPE
  bf.rope <- exp(bf_rope(post_cor,prior=prior_cor,null = ROPE)$log_BF)
  bf.sd <- exp(bf_parameters(post_cor,prior=prior_cor)$log_BF)
  
  # Summary of values for the correlation parameter
  ES.summary <- c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"P.MAP"=pmap,"BF.0"=bf.sd,
                  "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
                  "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
                  "BF.ROPE"=bf.rope)
  
  # Now make the function output
  result <- list("Par.Summary"=model_pars,"ES.summary"=ES.summary,"BRM"=cor_brm,"posterior.samples"=post,
                 "ROPE"=c("ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2]),
                 "BF.10"=bf.sd,"BF.01"=1/bf.sd,"BF.ROPE"=bf.rope,"BF.inROPE"=1/bf.rope,ordinal=FALSE,
                 "ci"=ci,"ci.type"=ifelse(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile"),"eti","hdi"),
                 "samples"=nsamp,"N"=N,name.x=deparse(x),name.y=deparse(y),"Call"=match.call())
  class(result) <- "bayes.cortest"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Pretty print function for Bayesian correlations
print.bayes.cortest <- function(x,verbose=F,dig=3){
  cat("Call:")
  print(x$Call)
  if(x$ordinal==T){
    cat(paste0("\nBayesian polyserial correlation with Region of Practical Equivalence (ROPE) rho = [",paste(x$ROPE,collapse=", "),"]:"))
  } else{
    cat(paste0("\nRobust Bayesian correlation with Region of Practical Equivalence (ROPE) rho = [",paste(x$ROPE,collapse=", "),"]:"))
  }
  cat(paste0("\n n = ",x$N," complete cases"))
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  if(verbose){
    cat(paste0("\n\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):\n\n"))
    print(round(x$Par.Summary,dig))
  } else{
    cat("\n")
  }
  if(x$ordinal==T){
    cat(paste0("\n Polyserial Correlation Coefficient (derived from ordered-probit model; Breen et al., 2014):"))
  } else{
    cat(paste0("\n Robust Correlation Coefficient (based on multivariate Student-t distribution):")) 
  }
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$ES[2],dig),
             " [",round(x$ES[3],dig),", ",round(x$ES[4],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$ES[5]>0.999,"> 0.999",paste0("= ",signif(x$ES[5],dig))),
             " (",ifelse(x$ES[5]>0.999,">99.9",signif(x$ES[5]*100,dig)),
             "% chance that r ",ifelse(sign(x$ES[2])==1,">","<")," 0)"))
  cat(paste0("\n  ",ifelse(x$ES[10]<0.001,"<0.1",signif(x$ES[10]*100,dig)),"% of posterior density inside [",
             paste(x$ROPE,collapse=", "),"] (",ifelse(x$ES[11]>0.999,">99.9",signif(x$ES[11]*100,dig)),"% exceeds ROPE)"))
  cat(paste0("\n\n Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$BF.10>999 | x$BF.10 < 0.001,
                    format(x$BF.10,dig=3,scientific=T),
                    signif(x$BF.10,dig=3)),
             " (BF01 = ",ifelse(x$BF.01>999 | x$BF.01 < 0.001,
                                format(x$BF.01,dig=3,scientific=T),
                                signif(x$BF.01,dig=3)),") — ",
             bf_evidence(x$BF.10,model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$BF.ROPE>999 | x$BF.ROPE < 0.001,
                                                           format(x$BF.ROPE,dig=3,scientific=T),
                                                           signif(x$BF.ROPE,dig=3)),
             " (BF.inROPE = ",ifelse(x$BF.inROPE>999 | x$BF.inROPE < 0.001,
                                     format(x$BF.inROPE,dig=3,scientific=T),
                                     signif(x$BF.inROPE,dig=3)),") — ",
             bf_evidence(x$BF.ROPE,model = F))) # Bayes Factor Qualitative Interpretation
}

## bayes.pcor: Partial correlation estimated in a Bayesian manner based on robust Student-t correlation
## Only works in trivariate case at the moment (unlike polyserial pcor, which handles as many control vars as you want)
### Zack Williams, 09/25/2020
### Does not scale variables, estimates both intercepts and sigma parameters
### Partial correlations estimated in the typical frequentist manner for each posterior correlation matrix separately
## Default priors include:
# lkj(2) for correlation (same prior on partial cor)
# Student's t(3,0,MAD) for log(sigma) [BRMS default]
# Student's t(3,Mdn,MAD) for intercepts [BRMS default]
# gamma (2,0.1) for nu (same for both X and Y) [BRMS default]
# Tests r against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
bayes.pcor.brm <- function(x,y,z,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                          chains=5,iter=5000,warmup = 1000,seed=12345,cores=detectCores()-1,file=NULL,
                          overwrite.file=FALSE,...){
  ### Partial correlation function based on 3 correlations (can take either 3 correlations separately or input)
  pcor <- function(rs){
    rxy <- rs[1]
    rxz <- rs[2]
    ryz <- rs[3]
    rxyz <- (rxy - (rxz * ryz))/(sqrt(1 - rxz^2)*sqrt(1 - ryz^2))
    rxyz[which(abs(rxyz)>1)] <- 0 # Turn improper values to 0
    return(rxyz)
  }
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  xyz <- tibble("x"=x,"y"=y,"z"=z)
  xyz <- xyz[complete.cases(xyz),]
  N <- nrow(xyz)
  
  if(is.null(Prior)){
    pr <- prior(lkj(2), class = rescor)
  } else{
    pr <- Prior
  }
  
  cor_brm <- brm(bf(mvbind(x, y, z) ~ 1) + set_rescor(TRUE),
                 family=student,
                 data=xyz,
                 prior = pr,
                 iter = iter, warmup = warmup, 
                 chains = chains, cores = cores, 
                 seed = seed,
                 sample_prior = T,
                 ... # additional arguments
  )
  # Get posterior samples of model parameters
  post <- as_draws_df(cor_brm)
  post$"Rho.xy,z" = suppressWarnings(apply(post[,c(11,12,13)],1,pcor))
  post$"Rho.xz,y" = suppressWarnings(apply(post[,c(11,13,12)],1,pcor))
  post$"Rho.yz,x" = suppressWarnings(apply(post[,c(12,13,11)],1,pcor))
  post$"Rho.xy - Rho.xz"= suppressWarnings(unlist(post[,11] - post[,12])) # On scale of original data
  post$"Rho.xy - Rho.yz"= suppressWarnings(unlist(post[,11] - post[,13]))
  post$"Rho.xz - Rho.yz"= suppressWarnings(unlist(post[,12] - post[,13]))
  
  # Get prior/posterior correlation distributions
  post_pcor <- unname(unlist(post$`Rho.xy,z`))
  nsamp <- length(post_pcor)
  prior_pcor <- post$prior_rescor
  
  post <- post[,c(1:7,11:13,23:28,20:22)]
  
  model_pars <- MCMC_summary(post,ci=ci,ci.type = ci.type)
  row.names(model_pars)[1:10] <- c("Mu.X","Mu.Y","Mu.Z","Sigma.X","Sigma.Y","Sigma.Z","Nu","Rho.xy","Rho.xz","Rho.yz")
  set.seed(seed)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  MAP <- LaplacesDemon::Mode(post_pcor)
  Mdn <- median(post_pcor)
  if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
    CI <- quantile(post_pcor,c(ci.lo,ci.hi))
    names(CI) <- paste0("Q.",names(CI))
  } else{ # HDI
    CI <- unname(unlist(bayestestR::hdi(post_pcor,ci=ci))[2:3])
    names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
  }
  pd <- c(p_direction(post_pcor)) # Probability of direction
  p.exceeds.rope <- max(mean(post_pcor < ROPE[1]),mean(post_pcor > ROPE[2]))
  pmap <- c(p_map(post_pcor))
  p.rope <- 1 - mean(post_pcor < ROPE[1]) - mean(post_pcor > ROPE[2]) # Probability posterior ES within ROPE
  bf.rope <- exp(bf_rope(post_pcor,prior=prior_pcor,null = ROPE)$log_BF)
  bf.sd <- exp(bf_parameters(post_pcor,prior=prior_pcor)$log_BF)
  
  # Summary of values for the correlation parameter
  ES.summary <- c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"P.MAP"=pmap,"BF.0"=bf.sd,
                  "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
                  "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
                  "BF.ROPE"=bf.rope)
  
  Call <- match.call()
  # Now make the function output
  result <- list("Par.Summary"=model_pars,"ES.summary"=ES.summary,"BRM"=cor_brm,"posterior.samples"=post,
                 "ROPE"=c("ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2]),
                 "BF.10"=bf.sd,"BF.01"=1/bf.sd,"BF.ROPE"=bf.rope,"BF.inROPE"=1/bf.rope,ordinal=F,
                 "ci"=ci,"ci.type"=ifelse(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile"),"eti","hdi"),
                 "samples"=nsamp,"N"=N,name.x=as.character(Call)[2],name.y=as.character(Call)[3],
                 name.z=as.character(Call)[4],"Call"=Call)
  class(result) <- "bayes.pcor"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Pretty print function for Bayesian partial correlations
print.bayes.pcor <- function(x,verbose=F,dig=3){
  cat("Call:")
  print(x$Call)
  if(x$ordinal==T){
    cat(paste0("\nBayesian partial polyserial correlation with Region of Practical Equivalence (ROPE) rho = [",paste(x$ROPE,collapse=", "),"]:"))
    cat(paste0("\n Correlation between ",x$name.x," and ",x$name.y," after controlling for ",paste(x$name.z,collapse=" + ")))
  } else{
    cat(paste0("\nRobust Bayesian partial correlation with Region of Practical Equivalence (ROPE) rho = [",paste(x$ROPE,collapse=", "),"]:"))
    cat(paste0("\n Correlation between ",x$name.x," and ",x$name.y," after controlling for ",x$name.z))
  }
  cat(paste0("\n n = ",x$N," complete cases"))
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  if(verbose){
    cat(paste0("\n\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):\n\n"))
    print(round(x$Par.Summary,dig))
  } else{
    cat("\n")
  }
  
  if(x$ordinal==T){
    cat(paste0("\n Partial Polyserial Correlation Coefficient (derived from ordered-probit model; Breen et al., 2014):"))
  } else{
    cat(paste0("\n Robust Partial Correlation Coefficient (based on multivariate Student-t distribution):"))
  }
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$ES[2],dig),
             " [",round(x$ES[3],dig),", ",round(x$ES[4],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$ES[5]>0.999,"> 0.999",paste0("= ",signif(x$ES[5],dig))),
             " (",ifelse(x$ES[5]>0.999,">99.9",signif(x$ES[5]*100,dig)),
             "% chance that r ",ifelse(sign(x$ES[2])==1,">","<")," 0)"))
  cat(paste0("\n  ",ifelse(x$ES[10]<0.001,"<0.1",signif(x$ES[10]*100,dig)),"% of posterior density inside [",
             paste(x$ROPE,collapse=", "),"] (",ifelse(x$ES[11]>0.999,">99.9",signif(x$ES[11]*100,dig)),"% exceeds ROPE)"))
  cat(paste0("\n\n Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$BF.10>999 | x$BF.10 < 0.001,
                    format(x$BF.10,dig=3,scientific=T),
                    signif(x$BF.10,dig=3)),
             " (BF01 = ",ifelse(x$BF.01>999 | x$BF.01 < 0.001,
                                format(x$BF.01,dig=3,scientific=T),
                                signif(x$BF.01,dig=3)),") — ",
             bf_evidence(x$BF.10,model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$BF.ROPE>999 | x$BF.ROPE < 0.001,
                                                           format(x$BF.ROPE,dig=3,scientific=T),
                                                           signif(x$BF.ROPE,dig=3)),
             " (BF.inROPE = ",ifelse(x$BF.inROPE>999 | x$BF.inROPE < 0.001,
                                     format(x$BF.inROPE,dig=3,scientific=T),
                                     signif(x$BF.inROPE,dig=3)),") — ",
             bf_evidence(x$BF.ROPE,model = F))) # Bayes Factor Qualitative Interpretation
}

### bayes.IDAcor.brm - uses a multilevel model to calculate summary correlation from IDA
### Zack Williams, 12/28/21 
## Default priors include:
# Student's t(3,0,2.5) for intercepts [BRMS default]
# N(0,0.5) for b_yx (slope parameter) [r_yx = b_yx/sqrt(b_yx^2 + 1) if x is standardized]
# This puts ~50% of prior weight between ±0.33, 95% between ±0.7
# Tests r against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
bayes.IDAcor.brm <- function(form,data=NULL,family=brms::student,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                             chains=10,iter=2500,warmup = 1000,seed=12345,cores=detectCores()-1,file=NULL,
                             overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  ## Accommodate imputed data in the form of DF lists
  if(!is.data.frame(data) & !is.data.frame(data[[1]])){stop("Variable 'data' is empty and/or invalid. Please input a data frame or list of data frame objects.")}
  
  mod_form <- as.formula(form)
  vars <- all.vars(mod_form)
  n.covar <- length(vars) - 2
  
  for(v in vars){
    if(is.data.frame(data)){
      if(!v %in% names(data)){stop(paste0("Variable '",v,"' not found! Undefined columns selected."))}
    } else{ # List of imputed data frames
      if(!v %in% names(data[[1]])){stop(paste0("Variable '",v,"' not found! Undefined columns selected."))}
    }
  }
  # Outcome data
  DV <- vars[1]
  IV <- vars[2]
  covars <- vars[-c(1:2)]
  
  # Standardize data (done on every multiple imputation if list)
  # If single data frame, listwise deletion (complete case analysis) is used
  if(is.data.frame(data)){
    incomplete <- which(!complete.cases(data[,vars]))
    if(length(incomplete) > 0){
      message(paste0("Removed ",length(incomplete)," cases with one or more NA values."))
      data <- data[-incomplete,vars]
    } else{
      data <- data[,vars]
    }
    N <- nrow(data)
    data <- standardize(data)
  } else{ # List of data frames - assumed to be complete (imputed) data
    N <- nrow(data[[1]])
    data <- lappply(data,standardize)
  }
  
  if(is.null(Prior)){
    pr <- c(prior(normal(0,0.5),class = b),prior(student_t(3, 0, 2.5),class = Intercept),prior(lkj(2), class = cor))
  } else{
    pr <- Prior
  }
  
  # Actually fit brm model
  if(is.data.frame(data)){
    mod_cor <- suppressWarnings(brm(
      form,
      prior=pr,
      family=family,
      data = data, 
      cores=cores,
      chains=chains,
      seed=seed,sample_prior = T,
      iter=iter,warmup = warmup,
      ... # additional arguments can go here
    ))
  } else{ # List of data frames
    mod_cor <- suppressWarnings(brm_multiple(
      form,
      prior=pr,
      family=family,
      data = data, 
      cores=cores,
      chains=chains,
      seed=seed,sample_prior = T,
      iter=iter,warmup = warmup,
      ... # additional arguments can go here
    ))
  }
  
  post <- as_draws_df(mod_cor)
  # Now extract correlation
  post$COR <- suppressWarnings(unname(unlist(beta_to_pcor(post[,grep(paste0("^b_",IV,"$"),names(post))]))))
  post$prior_COR <- unlist(beta_to_pcor(post$prior_b))
  set.seed(seed)
  post$pred_cor <- beta_to_pcor(sapply(1:nrow(post),function(i){
    set.seed(seed + i)
    return(suppressWarnings(rnorm(1,unlist(post[i,grep(paste0("^b_",IV,"$"),names(post))]),unlist(post[i,grep(paste0("^sd_.*__",IV,"$"),names(post))]))))
  }))
  
  summ_cor <- post_summ(post$COR,ci = ci,ci.type=ci.type,ROPE = ROPE,prior=post$prior_COR)
  summ_pred <- post_summ(post$pred_cor,ci = ci,ci.type=ci.type,ROPE = ROPE)[c(1:7,10:13)]
  summ_pred <- c(summ_pred,"Pdiff"=1-unname(summ_pred)[7])
  # Heterogeneity parameters (tau^2 on outcome scale, I^2 is standardized)
  post$tau2 <- unname(unlist(suppressWarnings(post[,grep(paste0("^sd_.*_",IV),names(post))])))^2
  var_b <- suppressWarnings(var(unlist(post[,2])))
  post$I2 <- 100*(post$tau2/(var_b + post$tau2))
  post$ICC <- 1 - post$sigma^2/rowSums(suppressWarnings(post[,grep("^sd_|^sigma",names(post))])^2)
  
  
  # Make result object
  result <- list(ES=summ_cor,ES_Pred=summ_pred,tau2=c("tau2"=post_summ(post$tau2,ci = ci,ci.type=ci.type)[4:6]),I2=c("I2"=post_summ(post$I2,ci = ci,ci.type=ci.type)[4:6]),
                 ICC=c("ICC"=post_summ(post$ICC,ci = ci,ci.type=ci.type)[4:6]),ci=ci,ci.type=ci.type,N=N,ROPE=ROPE,
                 posterior=post,form=form,vars=vars,BRM=mod_cor,"mi"=ifelse(is.data.frame(data),0,length(data)),samples=nrow(post),Call=match.call())
  class(result) <- "IDAcor"
  result$BF.10 <- unname(exp(result$ES[8]))
  result$BF.01 <- unname(exp(-result$ES[8]))
  result$BF.ROPE <- unname(exp(result$ES[14]))
  result$BF.inROPE <- unname(exp(result$ES[14]))
  # Probabilities of the summary effect being (a) very small (<0.1), (b) small (0.1-0.3), medium (0.3-0.5), and large (>0.5)
  result$P.vsmall <- mean(abs(post$COR) < 0.1)
  result$P.small <- mean(abs(post$COR) >= 0.1 & abs(post$COR) < 0.3)
  result$P.med <- mean(abs(post$COR) >= 0.3 & abs(post$COR) < 0.5)
  result$P.large <- mean(abs(post$COR) >= 0.5)
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Pretty-print function
print.IDAcor <- function(x,dig=3){
  cat("Call:")
  print(x$Call)
  cat(paste0("\nRobust Bayesian partial correlation with Region of Practical Equivalence (ROPE) rho = [",paste(x$ROPE,collapse=", "),"]:"))
  cat(paste0("\n Correlation between ",x$vars[1]," and ",x$vars[2]," after controlling for ",x$vars[3]))
  if(x$mi > 0){
    cat(paste0("\n n = ",x$N," cases [multiply imputed ",x$mi,"-fold]"))
  } else{
    cat(paste0("\n n = ",x$N," complete cases"))
  }
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  
  cat(paste0("\n\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):"))
  
  cat(paste0("\n  Mdn r [",x$ci*100,"% CrI] = ",format(x$ES[4],digits=1,nsmall=dig),
             " [",format(x$ES[5],digits=1,nsmall=dig),", ",format(x$ES[6],digits=1,nsmall=dig),"]"))
  cat(paste0("\n   Probability of Direction (Pd) ",ifelse(x$ES[7]>0.999,"> 0.999",paste0("= ",signif(x$ES[7],dig))),
             " (",ifelse(x$ES[7]>0.999,">99.9",signif(x$ES[7]*100,dig)),
             "% chance that r ",ifelse(sign(x$ES[4])==1,">","<")," 0)"))
  cat(paste0("\n   ",ifelse(x$ES[12]<0.001,"<0.1",signif(x$ES[12]*100,dig)),"% of posterior density inside [",
             paste(x$ROPE,collapse=", "),"] (",ifelse(x$ES[11]>0.999,">99.9",signif(x$ES[13]*100,dig)),"% exceeds ROPE)"))
  cat(paste0("\n P.vsmall (<0.1) = ",round(x$P.vsmall,3)," | P.small (0.1-0.3) = ",round(x$P.small,3)," | P.med (0.3-0.5) = ",
             round(x$P.med,3)," | P.large (>0.5) = ",round(x$P.large,3)))
  
  cat(paste0("\n\n Heterogeneity Metrics:"))
  cat(paste0("\n  ",x$ci*100,"% Prediction Interval = ",
             "[",format(x$ES_Pred[5],digits=1,nsmall=dig),", ",format(x$ES_Pred[6],digits=1,nsmall=dig),"], P.diff ",(ifelse(x$ES_Pred[12] < 0.001,"< 0.001",paste0("= ",round(x$ES_Pred[12],dig)))),
             " (",(ifelse(x$ES_Pred[12] < 0.001,"<0.1",100*round(x$ES_Pred[12],dig))),"% chance of r.pred having opposite sign)"
             ))
  tau2 <- format(x$tau2,nsmall=1,digits=2)
  I2 <- format(x$I2,nsmall=1,digits=1)
  ICC <- format(round(x$ICC,digits=3),digits=3)
  cat(paste0("\n  tau^2 (raw heterogeneity) = ",tau2[1],", ",x$ci*100,"% CrI [",tau2[2],", ",tau2[3],"]"))
  cat(paste0("\n  I^2 (standardized heterogeneity) = ",I2[1],"%, ",x$ci*100,"% CrI [",I2[2],", ",I2[3],"]"))
  cat(paste0("\n  ICC (random effect variance/total variance) = ",ICC[1],", ",x$ci*100,"% CrI [",ICC[2],", ",ICC[3],"]"))
  
  cat(paste0("\n\n Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$BF.10>999 | x$BF.10 < 0.001,
                    format(x$BF.10,dig=3,scientific=T),
                    signif(x$BF.10,dig=3)),
             " (BF01 = ",ifelse(x$BF.01>999 | x$BF.01 < 0.001,
                                format(x$BF.01,dig=3,scientific=T),
                                signif(x$BF.01,dig=3)),") — ",
             bf_evidence(x$BF.10,model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$BF.ROPE>999 | x$BF.ROPE < 0.001,
                                                           format(x$BF.ROPE,dig=3,scientific=T),
                                                           signif(x$BF.ROPE,dig=3)),
             " (BF.inROPE = ",ifelse(x$BF.inROPE>999 | x$BF.inROPE < 0.001,
                                     format(x$BF.inROPE,dig=3,scientific=T),
                                     signif(x$BF.inROPE,dig=3)),") — ",
             bf_evidence(x$BF.ROPE,model = F))) # Bayes Factor Qualitative Interpretation
}

### bayes.cor.polyserial.brm - uses a probit model to calculate polyserial (also biserial) correlation
### Based on equations from Breen et al., 2014 (https://doi.org/10.1177/0049124114544224)
### Zack Williams, 09/26/2020 - updated 06/21/21 to standardize x by default
## Default priors include:
# Student's t(3,0,2.5) for intercepts [BRMS default]
# N(0,1) for b_yx (slope parameter) [r_yx = b_yx/sqrt(b_yx^2 + 1) if x is standardized]
# This puts ~50% of prior weight between ±0.55, 95% between ±0.9
# Tests r against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
bayes.cor.polyserial.brm <- function(form,data=NULL,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                             chains=5,iter=10000,warmup = 2000,seed=12345,std=TRUE,
                             cores=detectCores()-1,file=NULL,overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  if(length(form) != 3 || length(all.vars(form)) > 2 || op(form) != "~" || is.null(lhs(form)) || is.null(rhs(form))){
    stop("Formula is incorrectly specified! Please enter formula of the type 'y ~ x'.")
  }
  if(!is.data.frame(data)){stop("Variable 'data' is empty and/or not a data frame.")}
  y <- lhs(form)
  x <- rhs(form)
  if(!deparse(x) %in% names(data)){stop(paste0("Variable '",deparse(x),"' not found! Undefined columns selected."))}
  if(!deparse(y) %in% names(data)){stop(paste0("Variable '",deparse(y),"' not found! Undefined columns selected."))}
  if(!is.numeric(data[,deparse(x)])){stop("Predictor variable does not appear to be continuous! Check your data/formula.")}
  # Data
  ncat_y <- length(levels(as.factor(data[,deparse(y)])))
  if(ncat_y == 2){
    fam <- bernoulli("probit") # if only 2 categories do non-ordinal regression
  } else{
    fam <- cumulative("probit")
  }
  
  if(ncat_y > 10){message("Outcome has more than 10 categories. Polyserial correlation may be unnecessary. ")}
  data <- data[complete.cases(data[,all.vars(form)]),all.vars(form)] # Remove NAs
  # Convert y variable to ordered if factor
  if(!is.ordered(data[,deparse(y)])){data[,deparse(y)] <- ordered(data[,deparse(y)])}
  N <- nrow(data)
  
  if(std){
    data <- effectsize::standardize(data)
  }
  sdx <- sd(data[,deparse(x)])
  
  # Construct formula
  bform <- bf(form)
  
  # Use default priors unless otherwise specified (N(0,1) for b, t(3,0,2.5) for intercepts)
  if(is.null(Prior)){
    pr <- c(prior(normal(0,1), class = b),
            prior(student_t(3,0,2.5), class = Intercept)) # default prior for intercepts, but spelled out here anyways
  } else{
    pr <- Prior
  }
  
  mod_poly <- brm(bform,
                  family=fam,
                  data=data,
                  prior = pr,
                  sample_prior = T,
                  iter = iter, warmup = warmup, 
                  chains = chains, cores = cores, 
                  seed = seed,
                  ... # additional arguments
  )

  # Get posterior samples of model parameters
  post <- as_draws_df(mod_poly)
  model_pars <- MCMC_summary(post,ci=ci,ci.type = ci.type)
  model_pars <- model_pars[c(1:ncat_y),]
  # Get prior/posterior correlation distributions
  post_byx <- suppressWarnings(unlist(post[,ncat_y]))
  # Posterior correlation (based on formula from Breen et al., 2014)
  post$r_yx <- post_cor <- post_byx * sdx/sqrt(post_byx^2 * sdx^2 + 1)
  # Prior correlation
  post$prior_r <- prior_cor <- post$prior_b * sdx/sqrt(post$prior_b^2 * sdx^2 + 1)
  nsamp <- length(post_cor)
  set.seed(seed)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  MAP <- LaplacesDemon::Mode(post_cor)
  Mdn <- median(post_cor)
  if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
    CI <- quantile(post_cor,c(ci.lo,ci.hi))
    names(CI) <- paste0("Q.",names(CI))
  } else{ # HDI
    CI <- unname(unlist(bayestestR::hdi(post_cor,ci=ci))[2:3])
    names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
  }
  pd <- c(p_direction(post_cor)) # Probability of direction
  p.exceeds.rope <- max(mean(post_cor < ROPE[1]),mean(post_cor > ROPE[2]))
  pmap <- c(p_map(post_cor))
  p.rope <- 1 - mean(post_cor < ROPE[1]) - mean(post_cor > ROPE[2]) # Probability posterior ES within ROPE
  bf.rope <- exp(bf_rope(post_cor,prior=prior_cor,null = ROPE)$log_BF)
  bf.sd <- exp(bf_parameters(post_cor,prior=prior_cor)$log_BF)
  
  # Summary of values for the correlation parameter
  ES.summary <- c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"P.MAP"=pmap,"BF.0"=bf.sd,
                  "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
                  "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
                  "BF.ROPE"=bf.rope)
  
  # Now make the function output
  result <- list("Par.Summary"=model_pars,"ES.summary"=ES.summary,"BRM"=mod_poly,"posterior.samples"=post,
                 "ROPE"=c("ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2]),
                 "BF.10"=bf.sd,"BF.01"=1/bf.sd,"BF.ROPE"=bf.rope,"BF.inROPE"=1/bf.rope,ordinal=TRUE,
                 "ci"=ci,"ci.type"=ifelse(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile"),"eti","hdi"),
                 "samples"=nsamp,"N"=N,name.x=deparse(x),name.y=deparse(y),"Call"=match.call())
  class(result) <- "bayes.cortest"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

### bayes.pcor.polyserial.brm - Bayesian partial polyserial correlations using a probit model to calculate polyserial (also biserial) correlations
### Based on equations from Breen et al., 2014 (https://doi.org/10.1177/0049124114544224)
### Zack Williams, 09/26/2020 - updated 06/21/21 to standardize x by default
## Default priors include:
# Student's t(3,0,2.5) for intercepts [BRMS default]
# N(0,1) for b_ (slope parameter), same for all predictors
# This puts ~50% of prior weight between ±0.55, 95% between ±0.9
# Tests r against 0 (Savage-Dickey ratio) and a ROPE of -0.1,0.1
bayes.pcor.polyserial.brm <- function(form,data=NULL,Prior=NULL,ci=0.95,ci.type="hdi",ROPE=c(-0.1,0.1),
                                     chains=5,iter=10000,warmup = 2000,seed=12345,
                                     cores=detectCores()-1,file=NULL,overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  if(length(form) != 3 || length(all.vars(form)) < 3 || op(form) != "~" || is.null(lhs(form)) || is.null(rhs(form))){
    stop("Formula is incorrectly specified! Please enter formula of the type 'y ~ x + z1 + ... + zn'.")
  }
  if(!is.data.frame(data)){stop("Variable 'data' is empty and/or not a data frame.")}
  y <- lhs(form)
  preds <- rhs.vars(form) # All predictors (first one is "x" others are "z")
  x <- as.name(preds[1])
  z <- preds[-1]
  # Check for variables in data
  if(!deparse(x) %in% names(data)){stop(paste0("Variable '",deparse(x),"' not found! Undefined columns selected."))}
  if(!deparse(y) %in% names(data)){stop(paste0("Variable '",deparse(y),"' not found! Undefined columns selected."))}
  for(v in z){ # Go through control vars and make sure they're in the data
    if(!v %in% names(data)){stop(paste0("Variable '",v,"' not found! Undefined columns selected."))}
  }
  if(!is.numeric(data[,deparse(x)])){stop(paste0("Predictor variable '",deparse(x),"' does not appear to be continuous! Check your data/formula."))}
  # Data
  ncat_y <- length(levels(as.factor(data[,deparse(y)])))
  if(ncat_y == 2){
    fam <- bernoulli("probit") # if only 2 categories do non-ordinal regression
  } else{
    fam <- cumulative("probit")
  }
  
  if(ncat_y > 10){message("Outcome has more than 10 categories. Polyserial correlation may be unnecessary. ")}
  data <- data[complete.cases(data[,all.vars(form)]),all.vars(form)] # Remove NAs
  # Convert y variable to ordered if factor
  if(!is.ordered(data[,deparse(y)])){data[,deparse(y)] <- ordered(data[,deparse(y)])}
  N <- nrow(data)
  
  if(std){
    data <- effectsize::standardize(data)
  }
  sdx <- sd(data[,deparse(x)])
  
  
  N <- nrow(data)
  # Get sd(x|z)
  sdx.z <- sd(residuals(lm(paste0(x," ~ ",paste0(z,collapse=" + ")),data=data)))
  # Convert y variable to numeric if factor
  if(is.factor(data[,deparse(y)])){
    data[,deparse(y)] <- as.numeric(data[,deparse(y)])
  }
  # Get rid of 0 values if any
  if(any(data[,deparse(y)]==0)){
    data[,deparse(y)] <- as.numeric(as.factor(data[,deparse(y)]))
  }
  
  # Construct formula
  bform <- bf(form)
  
  # Use default priors unless otherwise specified (N(0,1) for b, t(3,0,2.5) for intercepts)
  if(is.null(Prior)){
    pr <- c(prior(normal(0,1), class = b),
            prior(student_t(3,0,2.5), class = Intercept)) # default prior for intercepts, but spelled out here anyways
  } else{
    pr <- Prior
  }
  
  mod_poly <- brm(bform,
                  family=fam,
                  data=data,
                  prior = pr,
                  sample_prior = T,
                  iter = iter, warmup = warmup, 
                  chains = chains, cores = cores, 
                  seed = seed,
                  ... # additional arguments
  )
  
  # Get posterior samples of model parameters
  post <- as_draws_df(mod_poly)
  model_pars <- MCMC_summary(post,ci=ci,ci.type = ci.type)
  model_pars <- model_pars[1:(ncat_y+length(z)),]
  # Get prior/posterior correlation distributions
  post_byx.z <- suppressWarnings(unlist(post[,ncat_y]))
  # Posterior partial correlation (based on formula from Breen et al., 2014)
  post$r_yx.z <- post_pcor <- post_byx.z * sdx.z/sqrt(post_byx.z^2 * sdx.z^2 + 1)
  # Prior partialcorrelation
  post$prior_r <- prior_pcor <- post$prior_b * sdx.z/sqrt(post$prior_b^2 * sdx.z^2 + 1)
  nsamp <- length(post_pcor)
  set.seed(seed)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  MAP <- LaplacesDemon::Mode(post_pcor)
  Mdn <- median(post_pcor)
  if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
    CI <- quantile(post_pcor,c(ci.lo,ci.hi))
    names(CI) <- paste0("Q.",names(CI))
  } else{ # HDI
    CI <- unname(unlist(bayestestR::hdi(post_pcor,ci=ci))[2:3])
    names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
  }
  pd <- c(p_direction(post_pcor)) # Probability of direction
  p.exceeds.rope <- max(mean(post_pcor < ROPE[1]),mean(post_pcor > ROPE[2]))
  pmap <- c(p_map(post_pcor))
  p.rope <- 1 - mean(post_pcor < ROPE[1]) - mean(post_pcor > ROPE[2]) # Probability posterior ES within ROPE
  bf.rope <- exp(bf_rope(post_pcor,prior=prior_pcor,null = ROPE)$log_BF)
  bf.sd <- exp(bf_parameters(post_pcor,prior=prior_pcor)$log_BF)
  
  # Summary of values for the correlation parameter
  ES.summary <- c("MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=pd,"P.MAP"=pmap,"BF.0"=bf.sd,
                  "ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2],
                  "P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,
                  "BF.ROPE"=bf.rope)
  
  Call <- match.call()
  # Now make the function output
  result <- list("Par.Summary"=model_pars,"ES.summary"=ES.summary,"BRM"=mod_poly,"posterior.samples"=post,
                 "ROPE"=c("ROPE.lo"=ROPE[1],"ROPE.hi"=ROPE[2]),
                 "BF.10"=bf.sd,"BF.01"=1/bf.sd,"BF.ROPE"=bf.rope,"BF.inROPE"=1/bf.rope,ordinal=T,
                 "ci"=ci,"ci.type"=ifelse(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile"),"eti","hdi"),
                 "samples"=nsamp,"N"=N,name.x=deparse(x),name.y=deparse(y),
                 name.z=z,"Call"=Call)
  class(result) <- "bayes.pcor"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Convert standardized beta into partial correlation coefficient
beta_to_pcor <- function(b,sdx = 1,sdy = 1){
  return(b * (sdx/sdy)/sqrt(b^2 * sdx^2 + 1))
}


## bayes.proptest: Test proportions using a Bayes factor (from Bayesfactor package) and independent multinomial sampling
## Zack Williams, 08/29/20, updated 01/10/2021
# Also calculate out the proportion of interest for each group and its 95% CrI (HDI by default)
# Lastly, if 2 groups, calculate the (a) difference in proportions for 1st column, (b) odds ratio for first column
# ROPE updated to be OR = [0.833,1.2], (approximately equivalent to d = [-0.1, 0.1])
bayes.proptest <- function(X, Y=NULL,priorConcentration=1,which.prop=1,H0.prop=NULL,ROPE.OR=c(1/1.2,1.2),
                           iter=15000,ci.type="hdi",ci=0.95,rev=F,revpars=c(1,2),seed=12345,digits=3,
                           file=NULL,overwrite.file=FALSE){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  set.seed(seed)
  if(any(class(X)=="table")){
    tbl <- X
  } else{
    tbl <- table(X,Y)
  }
  dimtbl <- dim(tbl)
  if(rev){
    tbl <- Rev(tbl,revpars)
  }
  # Get CI quantiles
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  
  tbl_rownames <- rownames(tbl)
  tbl_colnames <- colnames(tbl)
  ngrps <- dimtbl[1]
  # Sample prior distribution (dirichlet (a,a,...))
  priors <- lapply(1:ngrps,function(i){
    brms::rdirichlet(iter,alpha = rep(priorConcentration,ncol(tbl)))
  })
  set.seed(seed)
  # Calculate the Bayes factor analytically for the point null
  BF <- contingencyTableBF(tbl,sampleType = "indepMulti",fixedMargin = "rows",priorConcentration=priorConcentration)
  chains <- posterior(BF, iterations = iter)
  post <- as.matrix(chains)[,grep("^pi",colnames(chains))]
  
  
  props <- sapply(1:ngrps,function(i){
    which.cols <- grep(paste0("^pi\\[",i),colnames(post),value = T) # Finds one proportion and calculates that
    post[,which.cols[which.prop]]/post[,grep("\\*",which.cols,value=T)] # returns prob of interest divided by row prob
  })
  prior_props <- sapply(1:ngrps,function(i){priors[[i]][,which.prop]})
  colnames(props) <- colnames(prior_props) <- sapply(1:ngrps,function(i){paste0("p(",tbl_colnames[which.prop],"|",tbl_rownames[i],")")})
  
  # Null point for Pd: NA for proportions (unless specified H0s), 1 for OR, and 0 for difference
  if(length(H0.prop)==length(ngrps)){ # Full vector of proportion H0s specified
    H0 <- c(H0.prop,0) 
  } else if(!is.null(H0.prop)){
    # if single value specified, repeat that value
    if(length(H0.prop)==1){
      H0 <- c(rep(H0.prop,ngrps),0)
    } else{
      message(paste0("H0.prop not the correct length. Ignoring arguments.\n"))
      H0 <- c(rep(NA,ngrps),0)
    }
  } else {
    H0 <- c(rep(NA,ngrps),0)
  }
  
  # If 2 groups being compared only, calculate proportion difference and OR
  if(ncol(props)==2){ 
    p1 <- props[,1] # proportion for group 1
    p2 <- props[,2] # proportion for group 2
    propDiff <- p1 - p2
    logOR <- log((p1/(1 - p1))/(p2/(1 - p2)))
    props <- cbind(props,propDiff,logOR)
    # Now do the same for the prior
    pp1 <- prior_props[,1] # proportion for group 1
    pp2 <- prior_props[,2] # proportion for group 2
    prior_propDiff <- pp1 - pp2
    logpOR <- log((pp1/(1 - pp1))/(pp2/(1 - pp2)))
    prior_props <- cbind(prior_props,prior_propDiff,logpOR)
    # Rename difference column
    colnames(props)[3] <- colnames(prior_props)[3] <- paste(tbl_rownames,collapse=" - ") # Rename to show difference direction
  }
  
  post_summary <- t(sapply(which(colnames(props)!="logOR"),function(i){ # for everything except OR
    temp <- props[,i]
    pr <- prior_props[,i]
    M <- mean(temp)
    SD <- sd(temp)
    MAP <- LaplacesDemon::Mode(temp)
    Mdn <- median(temp)
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- quantile(temp,c(ci.lo,ci.hi))
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- unname(unlist(bayestestR::hdi(temp,ci=ci))[2:3])
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    if(!is.na(H0[i])){
      Pd <- mean(sign(temp-H0[i])==sign(Mdn-H0[i]))
      BF.10 <- exp(bayesfactor_parameters(temp,prior=pr,null=H0[i])$log_BF)
      BF.01 <- 1/BF.10
    } else{
      Pd <- BF.10 <- BF.01 <- NA
    }
    # if Odds ratio, treat it differently 
    
    return(c("M"=M,"SD"=SD,"MAP"=MAP,"Mdn"=Mdn,CI,"H0"=H0[i],"Pd"=Pd,"BF.10"=BF.10,"BF.01"=BF.01))
  }))
  rownames(post_summary) <- colnames(props)[which(colnames(props)!="logOR")]
  
  post <- cbind(post,props)
  
  result <- list("Table"=tbl,"BF.10"=exp(unlist(BF@bayesFactor$bf)),"BF.01"=1/exp(unlist(BF@bayesFactor$bf)),
                 "Summary"=post_summary,"samples"=iter,"a"=priorConcentration,N=sum(tbl),
                 "ROPE"=c("ROPE.lo"=ROPE.OR[1],"ROPE.hi"=ROPE.OR[2]),
                 "Prior"=priors,"Posterior"=post,"ci.type"=ci.type,"ci"=ci,Call=match.call())
  # If odds ratio is calculated, add it to the 
  if(any(colnames(props)=="logOR")){
    which.OR <- which(colnames(props)=="logOR")
    # Now calculate ROPE Bayes Factor
    prior_or <- prior_props[,which.OR]
    post_or <- props[,which.OR]
    M <- exp(mean(post_or))
    SD <- sd(exp(post_or))
    MAP <- exp(LaplacesDemon::Mode(post_or))
    Mdn <- exp(median(post_or))
    Pd <- c(p_direction(post_or))
    if(ci.type %in% c("ETI","eti","QI","qi","quantile","percentile")){ # ETI
      CI <- exp(quantile(post_or,c(ci.lo,ci.hi)))
      names(CI) <- paste0("Q.",names(CI))
    } else{ # HDI
      CI <- exp(unname(unlist(bayestestR::hdi(post_or,ci=ci))[2:3]))
      names(CI) <- c(paste0("HDI.",ci.lo*100,"%"),paste0("HDI.",ci.hi*100,"%"))
    }
    p.exceeds.rope <- max(mean(post_or < log(ROPE.OR[1])),mean(post_or > log(ROPE.OR[2])))
    p.rope <- 1 - mean(post_or < log(ROPE.OR[1])) - mean(post_or > log(ROPE.OR[2])) # Probability posterior ES within ROPE
    result$BF.ROPE <- bf.rope <- exp(bf_rope(post_or,prior=prior_or,null = log(ROPE.OR))$log_BF)
    result$BF.inROPE <- 1/bf.rope
    result$OR <- c("M"=M,"SD"=SD,"MAP"=MAP,"Mdn"=Mdn,CI,"Pd"=Pd,"P.ROPE"=p.rope,"P>ROPE"=p.exceeds.rope,"BF.ROPE"=bf.rope)
  } else{
    result$OR <- NULL
  }
  class(result) <- "bayes.proptest"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  return(result)
}

# Pretty printing for bayes.proptest objects
print.bayes.proptest <- function(x,dig=3,summary=FALSE){
  cat("Call:")
  print(x$Call)
  cat(paste0("\nBayesian contingency table test (Independent Multinomial [fixed row] method, a = ",x$a,"):"))
  cat(paste0("\n n = ",x$N," complete cases: "))
  cat(paste(rowSums(x$Table),rownames(x$Table),collapse = ", "))
  cat("\n")
  cat("\nObserved percentages in each group:\n")
  print(100*round(prop.table(x$Table,1),dig))
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  
  if(summary){
    cat(paste0("\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):\n\n"))
    print(round(x$Summary,dig))
  }
  
  if(!is.null(x$OR)){
    cat(paste0("\n Odds Ratio (",paste(rownames(x$Table),collapse="/"),"):"))
    cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$OR[4],dig),
               " [",round(x$OR[5],dig),", ",round(x$OR[6],dig),"]"))
    cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$OR[7]>0.999,"> 0.999",paste0("= ",signif(x$OR[7],dig))),
               " (",ifelse(x$OR[7]>0.999,">99.9",signif(x$OR[7]*100,dig)),
               "% chance that OR ",ifelse(x$OR[4]>=1,">","<")," 1)"))
    cat(paste0("\n  ",ifelse(x$OR[8]<0.001,"<0.1",signif(x$OR[8]*100,dig)),"% of posterior density inside [",
               paste(round(x$ROPE,3),collapse=", "),"] (",ifelse(x$OR[9]>0.999,">99.9",signif(x$OR[9]*100,dig)),"% exceeds ROPE)"))
  }
  cat(paste0("\n\n Bayes Factor vs. Independence (BF10): ",
             ifelse(x$BF.10>999 | x$BF.10 < 0.001,
                    format(x$BF.10,dig=3,scientific=T),
                    signif(x$BF.10,dig=3)),
             " (BF01 = ",ifelse(x$BF.01>999 | x$BF.01 < 0.001,
                                format(x$BF.01,dig=3,scientific=T),
                                signif(x$BF.01,dig=3)),") — ",
             bf_evidence(x$BF.10,model = F))) # Bayes Factor Qualitative Interpretation
  if(!is.null(x$OR)){
    cat(paste0("\n Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$BF.ROPE>999 | x$BF.ROPE < 0.001,
                                                             format(x$BF.ROPE,dig=3,scientific=T),
                                                             signif(x$BF.ROPE,dig=3)),
               " (BF.inROPE = ",ifelse(x$BF.inROPE>999 | x$BF.inROPE < 0.001,
                                       format(x$BF.inROPE,dig=3,scientific=T),
                                       signif(x$BF.inROPE,dig=3)),") — ",
               bf_evidence(x$BF.ROPE,model = F))) # Bayes Factor Qualitative Interpretation
  }
}


### mi.missForest: Multiple Imputation using MissForest
# Takes a dataframe, uses missForest to impute values, returns list of dataframes
# default 10 multiple imputations
# Automatically takes characters and turns them to factors. 
# ordered = vector of column names/numbers that should be turned into ordered categories 
# nominal = vector of column names/numbers that should be turned into unordered categories
# ... - input to actual missForest function
mi.missForest <- function(DF,m=10,ordered=NULL,nominal=NULL,cores=detectCores()-1,seed=12345,file=NULL,overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  # If given "ordered" as column names, convert to column numbers
  if(class(ordered)=="character"){
    ordered <- sapply(ordered,function(X){grep(paste0("^",X,"$"),names(DF))})
  }
  # If given "nominal" as column names, convert to column numbers
  if(class(nominal)=="character"){
    nominal <- sapply(nominal,function(X){grep(paste0("^",X,"$"),names(DF))})
  }
  # Search for character variables, turn to factors
  for(i in 1:ncol(DF)){
    if(i %in% ordered){
      # Since missForest can't handle categorical predictors with >53 categories, skip those columns
      DF[,i] <- factor(DF[,i],ordered=T)
    } else if(i %in% nominal | class(DF[,i])[1]=="character" & length(unique(DF[,i])) <= 53 | class(DF[,i])[1]=="logical"){
      DF[,i] <- factor(DF[,i],ordered=F)
    }
    if(length(levels(DF[,i])) > 53){DF[,i] <- as.factor(as.character(DF[,i]))}
  }
  # exclude columns that can't go into missForest
  exclude <- unlist(sapply(1:ncol(DF),function(X){
    if(is.numeric(DF[,X])){
      return(NULL)
    } else if(length(unique(DF[,X])) > 53){
      return (X)
    } else {
      return(NULL)
    }
  }))
  # Save these columns for incorporation back into full DF
  if(!is.null(exclude)){
    DF_excluded <- DF[,exclude,drop=FALSE]
    DF <- DF[,-exclude]
  } else{
    DF_excluded <- NULL
  }
  
  # Now run a loop with missForests
  result <- pblapply(1:m,function(it){
    set.seed(seed+it)
    imp <- missForest(DF,...)
    return(bind_cols(DF_excluded,imp$ximp))
  },cl=cores)
  
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file = file)
  }
  if(m==1){result <- result[[1]]} # Return DF rather than list
  return(result)
}

### mi.mice: Multiple Imputation using mice
# Takes a dataframe, uses mice to impute values, returns list of dataframes
# default 5 multiple imputations (same as 'mice' function)
# Automatically takes characters and turns them to factors. 
# ordered = vector of column names/numbers that should be turned into ordered categories 
# nominal = vector of column names/numbers that should be turned into unordered categories
# ... - input to actual mice function
mi.mice <- function(DF,m=5,ordered=NULL,nominal=NULL,seed=12345,nnet.MaxNWts=2000,return_mids=TRUE,
                    file=NULL,overwrite.file=FALSE,...){
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  # If given "ordered" as column names, convert to column numbers
  if(class(ordered)=="character"){
    ordered <- sapply(ordered,function(X){grep(paste0("^",X,"$"),names(DF))})
  }
  # If given "nominal" as column names, convert to column numbers
  if(class(nominal)=="character"){
    nominal <- sapply(nominal,function(X){grep(paste0("^",X,"$"),names(DF))})
  }
  # Search for character variables, turn to factors
  for(i in 1:ncol(DF)){
    if(i %in% ordered){
      DF[,i] <- factor(DF[,i],ordered=T)
    } else if(i %in% nominal | class(DF[,i])[1]=="character"){
      DF[,i] <- factor(DF[,i],ordered=F)
    }
  }
  
  # Now actually impute datasets
  result <- mice(DF,m=m,...,cluster.seed=seed,nnet.MaxNWts=nnet.MaxNWts)
  
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file = file)
  }
  if(!return_mids){
    result <- lapply(1:m,function(i){
      return(mice::complete(result,i))
    })
  }
  
  return(result)
}

### bf_brm_mi
# Calculate bayes factor from multiply imputed models using the average of BFs calculated from each imputed model
bf_brm_mi <- bayesfactor_brm_mi <- bayesfactor_brm_multiple <- bf_brm_multiple <- bayes_factor_brm_mi <- 
  function(modlist,blmod,cores=detectCores()-1,seed=12345){
    set.seed(seed)
    # If there's already logml calculated
    if("brmsfit" %in% class(modlist) & !is.null(modlist$logml)){
      logml_mods <- modlist$logml
      if("brmsfit" %in% class(blmod)){
        if(is.null(blmod$logml)){
          logml_bl <- bridge_sampler(blmod,cores=cores,silent=T)$logml
        } else{
          logml_bl <- blmod$logml
          if(length(logml_bl) != length(logml_mods)){
            stop("Models contain different numbers of imputations—check your inputs!")
          }
        }
      }
      logBFs <- logml_mods - logml_bl
    } else if(any(class(blmod)=="list")){ # Multiply-imputed baseline model can also be used (has to be BRM list) 
      if(length(blmod)!=length(modlist)){stop("Length of baseline and augmented model lists are not the same.")}
      # Calculate Bayes factors for each imputed dataset compared to the baseline model calculated with that same imputed dataset
      logBFs <- pbsapply(1:length(modlist),function(it){
        mod_it <- modlist[[it]]
        bl_it <- blmod[[it]]
        bayes_factor(mod_it,bl_it,silent=T,cores=cores,log=T)$bf
      })
    } else{
      logBFs <- pbsapply(1:length(modlist),function(it){
        mod_it <- modlist[[it]]
        bayes_factor(mod_it,blmod,silent=T,cores=cores,log=T)$bf
      })
    }
    
    result <- list("BF"=exp(mean(logBFs)),"BF.raw"=exp(logBFs),"logBF.raw"=logBFs,"k"=length(logBFs),"m1.name"=deparse(substitute(modlist)),"m0.name"=deparse(substitute(blmod)))
    class(result) <- "BF.mi"
    return(result)
  }
# Print results from previous function
print.BF.mi <- function(X,round=3){
  cat(paste0("Multiply imputed Bayes factor (k = ",X$k,") for ",X$m1.name,
             " over ",X$m0.name,": ",round(X$BF,digits = round),"\n\nImputed Bayes Factor Summary:\n"))
  bf.summ <- as.numeric(format(summ(X$BF.raw,dig = round+1,hdQ = T)[c(8,4,3,1,5,9)],nsmall=1,digits=round))
  names(bf.summ) <- paste0("BF.",c("Min","Q1","Mdn","Mean","Q3","Max"))
  if(any(bf.summ > 1e5)){bf.summ <- scales::scientific(bf.summ)}
  print(bf.summ)
}

### Pool BRM models (wrapper for combine_models)
pool.brm <- function(modlist){
  return(combine_models(mlist=modlist,check_data=FALSE))
}

### bayes_R2_diff: Compare Bayesian R2 of two different models
bayes_R2_diff <- function(fullmod,blmod,ci=0.95,round=3){
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  r2_full <- unlist(bayes_R2(fullmod,summary=FALSE))
  r2_bl <- unlist(bayes_R2(blmod,summary=FALSE))
  r2_diff <- r2_full - r2_bl
  result <- rbind(c(mean(r2_full),LaplacesDemon::Mode(r2_full),quantile(r2_full,c(0.5,ci.lo,ci.hi))),
                  c(mean(r2_bl),LaplacesDemon::Mode(r2_bl),quantile(r2_bl,c(0.5,ci.lo,ci.hi))),
                  c(mean(r2_diff),LaplacesDemon::Mode(r2_diff),quantile(r2_diff,c(0.5,ci.lo,ci.hi))))
  rownames(result) <- c("Full","BL","Diff")
  colnames(result)[1:3] <- c("Mean","MAP","Mdn")
  return(round(result,round))
}

## mvMode - uses mean-shift algorithm implemented in 'mvnfast::ms' to calculate the multivariate mode of a dataset
## default settings include initializing values at univariate means, using 0.1*covariance matrix as bandwidth
mvMode <- function(X,startfun=mean,H=0.1*cov(X),cores=detectCores()-1,...){
  pacman::p_load(mvnfast) # requires mvnfast package
  cnames <- colnames(X)
  # Meat of the function - actually use the mean-shift function
  result <- mvnfast::ms(X,init=apply(X,2,startfun),H=H,ncores=cores,...)$final
  names(result) <- cnames # name variables
  
  return(result)
}

### fitmodels_brm Function to fit BRMS models with all selected predictors
## Zack Williams, 08/09/20
## Input: a baseline fitted model with terms to be included in all models (can be just intercept)
## Also a vector of the exact names of candidate predictors as characters (e.g., c("foo","bar"))
## If you want interactions, they should be specifically specified (e.g., preds = c("foo","bar","foo:bar"))
# Options for the fit function can be programmed into the initial fit and will be applied to all models
fitmodels_brm <- function(bl_mod,preds,newprior=NULL){
  
  # Make list with all combinations of all numbers of predictors
  pred_combos <- unlist(unlist(lapply(1:length(preds),
                                      function(i){
                                        apply(combn(preds,i),2,function(x) list(x))
                                      }),recursive = F),recursive = F)
  
  # New addition for including interaction terms only if both predictors are in the model
  blpreds <- rownames(fixef(brm_CEtrials_bl))
  
  
  # Go through pred_combos, find models with interaction terms, flag if both predictors not included
  which_mods <- sapply(pred_combos,
                       function(varnames){
                         interact_terms <- unique(unlist(strsplit(grep(":",varnames,value=T),":")))
                         mod_terms <- c(blpreds,varnames) # all terms that are in model
                         if(is.null(interact_terms)){
                           return(TRUE)
                         } else{
                           return(all(interact_terms %in% mod_terms))
                         }
                       })
  # Include only models that have all terms within an included interaction
  pred_combos <- pred_combos[which_mods]
  
  cat(paste0("\nFitting ",length(pred_combos)," models. May take a while..."))
  
  modlist <- pblapply(pred_combos,function(preds_i){
    f <- formula(paste(deparse(formula(bl_mod)),
                       paste(preds_i,collapse="+"),sep = "+"))
    if(!is.null(newprior)){
      suppressWarnings(update(bl_mod,f,prior=newprior))
    } else{
      suppressWarnings(update(bl_mod,f))
    }
    
  })
  cat(" Done!\n")
  names(modlist) <- sapply(pred_combos,function(x) paste0("BL+",paste(x,collapse = "+")))
  modlist$BL <- bl_mod
  
  return(modlist)
}

### bayes_factor_BMA: uses Bayesian Model Averaging to derive posterior probabilities for each model
## Zack Williams, 08/13/2020, updated on 02/13/2023 to account for multiple imputations (mean of posterior model probabilities across imputations)
## Calculates inclusion Bayes factors for every variable, as well as Bayes factors for each model
bayes_factor_BMA <- bf_BMA <- function(modlist,mod_names=NULL,prior_prob=NULL,match.models=T,
                                       method="normal",which.BL=NULL,logml=NULL,
                                       cores=detectCores()-1,seed=123){
  set.seed(seed)
  # Check to make sure data is OK
  if(!any(class(modlist)=="list") | !all(sapply(modlist,function(mod){any(class(mod)=="brmsfit")}))){
    stop("Invalid input! 'modlist' must be a list of brmsfit objects.")
  }
  
  # If no prior prob set then equiprobable models
  if(is.null(prior_prob)){prior_prob <- rep(1/length(modlist), length(modlist))}
  # Throw some errors if input is wrong
  if(length(modlist) != length(prior_prob)){
    stop("Number of objects needs to match number of elements in prior_prob.", call. = FALSE)
  }
  if(!isTRUE(all.equal(sum(prior_prob), 1))){stop("Prior model probabilities do not sum to one.", call. = FALSE)}
  
  # Extract fixed effects from models
  modFEs <- mclapply(modlist,function(X){
    FE.names <- rownames(fixef(X))
    # If interaction orders are switched, turn them around to be consistent between models
    unname(sapply(FE.names,function(eff){
      paste(sort(unlist(str_split(eff,":"))),collapse=":")
    }))
  })
  allFEs <- unique(unlist(modFEs))
  which.FEs <- lapply(modFEs,function(X){which(allFEs %in% X)})
  which.FEs.interact <- grep(":",allFEs)
  
  # Same for random effects
  modREs <- mclapply(modlist,function(X){
    if(nrow(X$ranef)==0){return(NULL)}
    # If interaction orders are switched, turn them around to be consistent between models
    coefs <- unname(sapply(X$ranef$coef,function(eff){
      paste(sort(unlist(str_split(eff,":"))),collapse=":")
    }))
    paste0(coefs,"|",X$ranef$group)
  })
  allREs <- unique(unlist(modREs))
  if(!is.null(allREs)){
    which.REs <- lapply(modREs,function(X){which(allREs %in% X)})
    which.REs.interact <- grep(":",allREs) # Determine which are interactions
  } else{
    which.REs <- NULL
  }
  
  #### Uses bridgesampling package to calculate marginal likelihoods (unless logml provided)
  if(is.null(logml)){
    message("Calculating Model Likelihoods. May take a while...")
    mod_logmls <- pblapply(modlist,function(X){
      result <- X$logml # if logml already calculated, use that instead of recalculating
      if(is.null(result)){
        result <- bridge_sampler(X,cores=cores,method=method,silent=T)$logml # otherwise use bridgesampling to calculate
      }
      return(result)
    })
    if(any(is.na(mod_logmls))) {
      post_prob <- rep(NA_real_, length(modlist))
      warning("NAs in marginal likelihood values. No posterior probabilities calculated.", call. = FALSE)
    } 
    logml <- bind_cols(mod_logmls)
  } 
  
  # Code lifted from bridgesampling::.post_prob_calc
  e <- Brobdingnag::as.brob(exp(1))
  
  if(nrow(logml)==1){
    post_probs <- as.numeric(e^logml*prior_prob / sum(e^logml*prior_prob))
  } else{
    post_probs <- colMeans(t(apply(logml,1,function(logml_i){
      post_prob <- as.numeric(e^logml_i*prior_prob / sum(e^logml_i*prior_prob))
    })))
  }
  
  
  if(!isTRUE(all.equal(sum(post_probs), 1))){
    warning("Posterior model probabilities do not sum to one.", call. = FALSE)
  }
  # Now that we have the posterior model probabilities, we can calculate model odds and Bayes factors
  post_odds <- post_probs/(1-post_probs)
  prior_odds <- prior_prob/(1-prior_prob)
  bf_modelodds <- post_odds/prior_odds
  
  # Fixed Effect Bayes Factor calculation
  out_fes <- t(sapply(allFEs,function(eff){
    i <- which(allFEs==eff)
    mods_i <- sapply(which.FEs,function(FE){any(FE==i)})
    prior_prob_i <- sum(mods_i * prior_prob)
    
    if(prior_prob_i==1){ # If term in all models, then no Bayes Factor
      incbf_i <- NA
      post_prob_i <- 1
    } else{
      post_prob_i <- sum(mods_i * post_probs)
      # Check if interaction; if so, only look at probability of effect in models with component main effects
      if(match.models & eff %in% allFEs[which.FEs.interact]){
        component_effs <- which(allFEs %in% unlist(strsplit(eff,":")))
        eligible_mods <- which(sapply(which.FEs,function(mod){all(component_effs %in% mod)}))
        prior_prob_i <- prior_prob_i/sum(prior_prob[eligible_mods]) # divide prior prob by proportion of models interaction could be in
        post_prob_i <- post_prob_i/sum(post_probs[eligible_mods]) # recalibrate post prob in same way
      }
      incbf_i <- (post_prob_i/(1-post_prob_i))/(prior_prob_i/(1-prior_prob_i)) # posterior odds/prior odds (bayes factor)
    }
    return(c("Prior.Prob"=prior_prob_i,"Post.Prob"=post_prob_i,"BF.inclusion"=incbf_i,"BF.exclusion"=1/incbf_i))
  }))
  
  # Random effect Bayes Factor calculation
  out_res <- t(sapply(allREs,function(eff){
    i <- which(allREs==eff)
    mods_i <- sapply(which.REs,function(RE){any(RE==i)})
    prior_prob_i <- sum(mods_i * prior_prob)
    
    if(prior_prob_i==1){ # If term in all models, then no Bayes Factor
      incbf_i <- NA
      post_prob_i <- 1
    } else{
      post_prob_i <- sum(mods_i * post_probs)
      # Check if interaction; if so, only look at probability of effect in models with component random effects
      if(match.models & eff %in% allREs[which.REs.interact]){
        component_effs <- which(allREs %in% unlist(strsplit(gsub(" | .*$","",eff),":")))
        eligible_mods <- which(sapply(which.REs,function(mod){all(component_effs %in% mod)}))
        prior_prob_i <- prior_prob_i/sum(prior_prob[eligible_mods]) # divide prior prob by proportion of models interaction could be in
        post_prob_i <- post_prob_i/sum(post_probs[eligible_mods]) # recalibrate post prob in same way
      }
      incbf_i <- (post_prob_i/(1-post_prob_i))/(prior_prob_i/(1-prior_prob_i)) # posterior odds/prior odds (bayes factor)
    }
    return(c("Prior.Prob"=prior_prob_i,"Post.Prob"=post_prob_i,"BF.inclusion"=incbf_i,"BF.exclusion"=1/incbf_i))
  }))
  # Create Model names if there are none
  if(is.null(mod_names)){
    baseline_FEs <- which(out_fes[,1]==1)
    if(length(out_fes)>0){
      baseline_FEs <- which(out_fes[,1]==1)
    } else{
      baseline_FEs <- NULL
    }
    if(length(out_res)>0){
      baseline_REs <- which(out_res[,1]==1)
    } else{
      baseline_REs <- NULL
    }
    mod_names <- sapply(1:length(modlist),function(i){
      FE_i <- which.FEs[[i]]
      RE_i <- which.REs[[i]]
      if(all(FE_i %in% baseline_FEs) & all(RE_i %in% baseline_REs)){
        name_i <- "BL"
      } else if(all(1:length(allFEs) %in% FE_i) & all(1:length(allREs) %in% RE_i)){
        "Full Model"
      } else{
        FE_noBL <- FE_i[-which(FE_i %in% baseline_FEs)]
        RE_noBL <- RE_i[-which(RE_i %in% baseline_REs)]
        name_i <- gsub(" \\+ $","",paste("BL",paste(allFEs[FE_noBL],collapse = " + "),
                                         paste(allREs[RE_noBL],collapse = " + "),sep = " + ")) # BL + FEs + REs
      }
    })
  }
  # BL model used as reference for Bayes Factors
  if(isTRUE(which.BL %in% 1:length(modlist))){ # which.BL declared in arguments
    logML_BL <- logml[which.BL]
    if(nrow(logml)>1){
      logML_BL <- rep.col(unlist(logml[,which.BL]),ncol(logml))
      BF.vsBL <- exp(colMeans(logml - logML_BL))
    } else{
      logML_BL <- logml[which.BL]
      BF.vsBL <- as.numeric(e^logml / e^logML_BL)
    }
  } else if("BL" %in% mod_names){ # Else use model named BL
    which.BL <- grep("^BL$",mod_names)[1]
    if(nrow(logml)>1){
      logML_BL <- rep.col(unlist(logml[,which.BL]),ncol(logml))
      BF.vsBL <- exp(colMeans(logml - logML_BL))
    } else{
      logML_BL <- logml[which.BL]
      BF.vsBL <- as.numeric(e^logml / e^logML_BL)
    }
  } else { # Else just use the first model
    message("No baseline model detected. Arbitrarily selecting first model as baseline.\n")
    which.BL <- 1
    if(nrow(logml)>1){
      logML_BL <- rep.col(unlist(logml[,which.BL]),ncol(logml))
      BF.vsBL <- exp(colMeans(logml - logML_BL))
    } else{
      logML_BL <- logml[which.BL]
      BF.vsBL <- as.numeric(e^logml / e^logML_BL)
    }
  }
  
  # Now create output
  model_out <- cbind("Prior.Prob"=prior_prob,"Post.Prob"=post_probs,"BF.model"=bf_modelodds,"1/BF"=1/bf_modelodds,
                       "BF.BL"=BF.vsBL) # Original output
  

  rownames(model_out) <- mod_names
  # Find first and second best models
  mod_order <- order(model_out[,2],decreasing = T)
  best.mod <- mod_names[mod_order[1]]
  second.best.mod <- mod_names[mod_order[2]]
  if(nrow(logml)>1){
    BF.1vs2 <- unname(exp(colMeans(logml[,mod_order[1]] - logml[,mod_order[2]])))
  } else{
    BF.1vs2 <- as.numeric(e^logml[mod_order[1]] / e^logml[mod_order[2]])
  }
  
  # Abridge the names of model_out if too long
  if(max(nchar(mod_names)) > 60){
    short_names <- sapply(mod_names,function(NAME){
      effs <- unlist(strsplit(NAME," \\+ "))
      NEWNAME <- paste(sapply(effs,function(EFF){
        if(str_detect(EFF,"\\|")){ # Random effect term
          RE_split <- unlist(str_split(EFF," \\| "))
          # Determine whether slope is an interaction, trim accordingly
          if(str_detect(RE_split[1],":")){
            RE1_split <- unlist(str_split(RE_split[1],":"))
            RE1 <- paste(strtrim(RE1_split,3),collapse=":")
          } else{
            RE1 <- strtrim(RE_split[1],3)
          }
          # Trim random effect name
          RE2 <- strtrim(RE_split[2],3)
          return(paste(RE1,RE2,sep="|"))
        } else if(str_detect(EFF,":")){ # Interaction term but not RE
          EFF_split <- unlist(str_split(EFF,":"))
          return(paste(strtrim(EFF_split,3),collapse=":"))
        } else{ # Not an interaction term (simply trim)
          return(strtrim(EFF,3))
        }
      }),collapse="/")
      return(NEWNAME)
    })
    short_names <- unname(sapply(short_names,function(X){
      if(nchar(X)>60){X <- str_trunc(X,60)}
      return(X)
    }))
    rownames(model_out) <- short_names
  }
  
  result <- list("Models"=model_out,
                 "Models.Ordered"=model_out[mod_order,],"logml"=logml,
                 "model.weights"=model_out[,2],"method"=method,
                 "Call"=match.call(),"match.models"=match.models,"Model.Names"=mod_names,
                 "Baseline.FEs"=paste(modFEs[[which.BL]],collapse=" + "),
                 "Baseline.REs"=paste(modREs[[which.BL]],collapse=" + "),
                 "Model.FEs"=modFEs,"Model.REs"=modREs,
                 "BestMod"=best.mod,"SecondBestMod"=second.best.mod,
                 "BF.BestvsBL"=BF.vsBL[mod_order[[1]]],"BF.1vs2"=BF.1vs2)
  if(length(out_fes)>0){ result$Fixed.Effects <- out_fes }
  if(length(out_res)>0){ result$Random.Effects <- out_res }
  names(result$logml) <- mod_names
  class(result) <- "BayesFactor.BMA"
  return(result)
}

print.BayesFactor.BMA <- function(x,dig=3,suppressBLVars=T){
  cat("Call:")
  print(x$Call)
  cat(paste0("\nBayes Factors calculated using bridge sampling (",x$method," method):"))
  cat(paste0("\n\n Baseline Fixed Effects: ",x$Baseline.FEs))
  if(!is.null(x$Random.Effects)){cat(paste0("\n Baseline Random Effects: ",x$Baseline.REs))}
  cat("\n\n Model Probabilities and Bayes Factors") 
  if(nrow(x$Models)>10){cat(" (top 10 most likely)")}
  cat(":\n")
  print(signif(x$Models.Ordered[1:(min(10,nrow(x$Models))),]),dig)
  cat(paste0("\nBest Model: ",x$BestMod," (PP = ",round(x$Models.Ordered[1,2],dig),
             "; Bayes factor vs. BL = ",signif(x$BF.BestvsBL,dig),")"))
  cat(paste0("\n Bayes factor vs. second best model (",x$SecondBestMod,") = ",signif(x$BF.1vs2,dig)))
  cat("\n")
  if(!is.null(x$Fixed.Effects)){
    cat("\n Fixed Effects")
    if(suppressBLVars & any(x$Fixed.Effects[,1]==1)){
      cat(paste0(" (omitted ",sum(x$Fixed.Effects[,1]==1)," variables present in all models)"))
      cat(":\n")
      print(signif(x$Fixed.Effects[which(x$Fixed.Effects[,1]!=1),],dig))
    } else{
      print(signif(x$Fixed.Effects,dig)) # no suppression of baseline variables (BF is NA)
    }
  }
  
  if(!is.null(x$Random.Effects)){
    cat("\n Random Effects")
    if(suppressBLVars & any(x$Random.Effects[,1]==1)){
      cat(paste0(" (omitted ",sum(x$Random.Effects[,1]==1)," variables present in all models)"))
      cat(":\n")
      print(signif(x$Random.Effects[which(x$Random.Effects[,1]!=1),,drop=FALSE],dig))
    } else{
      print(signif(x$Random.Effects,dig)) # no suppression of baseline variables (BF is NA)
    }
  }
}


### bayes.mediate - Bayesian Mediation analysis based on bayestestR::mediation()
### Zack Williams - 06/14/2021
### Largely based on cannibalized code from https://rdrr.io/cran/bayestestR/src/R/mediation.R
#################################################################################
## Can take input of either a multivariate brmsfit model (two regressions) or a formula to fit said model
## Additionally includes mediation effect sizes (proportion mediated, upsilon)
# TODO: Recode to accommodate multiple mediation
# TODO: Add Moderated Mediation
bayes.mediate <- function(form, treatment, mediator, moderator = NULL, response = NULL,data=NULL,Prior=NULL,ROPE=c(-0.1,0.1),
                          chains=5,iter=5000,warmup = 1000,seed=12345,cores=detectCores()-1,imputeMissing = TRUE,std=TRUE,
                          file=NULL,overwrite.file=FALSE,ci = 0.95, ci.type = "hdi", pattern = "b_%s_%s",...){
  
  # If file exists, just load that instead
  if(!is.null(file) & !overwrite.file){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    if(file.exists(file)){
      return(readRDS(file=file))
    }
  }
  
  # Check to see what kind of input is provided (formula vs. brmsfit vs. other)
  if(any(class(form) == "mvbrmsformula")){ # Formula provided - fit model internally
    if(!is.data.frame(data)){stop("Variable 'data' is empty and/or not a data frame.")} # throw error if no data provided
    # Get names of all variables (IVs and DVs)
    reg_vars <- unique(unlist(lapply(form$forms,function(f){
      f <- as.formula(f)
      vars <- all.vars(f)
    })))
    # Check to make sure variables are in data:
    if(!all(reg_vars %in% names(data))){
      stop(paste0("Undefined columns selected: ",paste(reg_vars[!reg_vars %in% names(data)],collapse = ", ")))
    }
    
    data <- data[,reg_vars] # Use only defined variables for imputation
    temp_data <- data.frame(lapply(data,function(v){
      if(is.numeric(v) & length(unique(v))==2){
        v <- as.factor(v)
      } else if(is.numeric(v) & length(unique(v))<=11){
        v <- ordered(v)
      }
      return(v)
    }))
    # use missForest to impute any missing values (can be turned off using "imputeMissing = FALSE")
    if(imputeMissing & any(is.na(data))){
      message("Missing values detected. Imputing using 'missForest' function (uses only variables in model).")
      # Make sure categorical variables are defined as such (≤11 vals = ordinal)
      set.seed(seed)
      mf <- missForest(temp_data)
      impute.OOBerror <- mf$OOBerror
      data <- mf$ximp
    } else{
      data <- temp_data
      impute.OOBerror <- NULL
    }
    # Now standardize all (numeric) predictors
    if(std){
      data <- effectsize::standardize(data)
    }
    
    # Use default priors unless otherwise specified (also N(0,1) for beta)
    if(is.null(Prior)){
      if(!std) message("Variables are not standardized but default priors are selected. Interpret results with caution!")
      pr <- c(brms::prior(normal(0,1), class = b))
    } else{
      pr <- Prior
    }
    # Actually fit brm model
    model <- brm(
      form + set_rescor(FALSE),
      prior=pr,
      data = data, 
      cores=cores,
      chains=chains,
      seed=seed,sample_prior = T,
      iter=iter,warmup = warmup,
      ... # additional arguments can go here
    )
  } else if(any(class(form)=="brmsfit")){ # Model provided
    # Use existing model (ignores any arguments for fitting model and imputation)
    impute.OOBerror <- NULL
    model <- form
    form <- model$formula
  } else { # Unrecognized argument
    stop("Function must be passed either a 'mvbrmsformula' object or fitted 'brmsfit' object as the first argument.")
  }
  
  ##################### Begin the code I cannibalized from bayestestR::mediation()
  .fix_factor_name <- function(model, variable) { # Helper function to fix brms-changed names
    # check for categorical. if user has not specified a treatment variable
    # and this variable is categorical, the posterior samples contain the
    # samples from each category of the treatment variable - so we need to
    # fix the variable name
    
    mf <- insight::get_data(model)
    if (variable %in% colnames(mf)) {
      check_fac <- mf[[variable]]
      if (is.factor(check_fac)) {
        variable <- sprintf("%s%s", variable, levels(check_fac)[nlevels(check_fac)])
      } else if (is.logical(check_fac)) {
        variable <- sprintf("%sTRUE", variable)
      }
    }
    variable
  }
  
  
  # only one HDI interval
  if (length(ci) > 1) ci <- ci[1]
  
  # model responses
  if(is.null(response)){
    response <- insight::find_response(model, combine = TRUE)
  }
  fix_mediator <- FALSE
  
  # find mediator, if not specified
  if(missing(mediator)){
    predictors <- insight::find_predictors(model, flatten = TRUE)
    mediator <- predictors[predictors %in% response]
    fix_mediator <- TRUE
  }
  
  # find treatment, if not specified
  if(missing(treatment)){
    predictors <- lapply(
      insight::find_predictors(model),
      function(.f) .f$conditional
    )
    
    treatment <- predictors[[1]][predictors[[1]] %in% predictors[[2]]][1]
    treatment <- .fix_factor_name(model, treatment)
  }
  
  mediator.model <- which(response == mediator)
  treatment.model <- which(response != mediator)
  
  if(fix_mediator) mediator <- .fix_factor_name(model, mediator)
  
  if (inherits(model, "brmsfit")) {
    response_name <- names(response)
  } else {
    response_name <- unname(response)
  }
  response_name_fixed <- grep(paste0("^",paste(unlist(strsplit(response[treatment.model],NULL)),
                                               collapse=".*"),"$"),names(model$data),value=T)
  
  # brms removes underscores from variable names when naming estimates
  # so we need to fix variable names here
  response <- names(response)
  
  # All posteriors:
  post_coef <- posterior_samples(model)
  
  # Direct effect: coef(treatment) from model_y_treatment
  coef_treatment <- sprintf(pattern, response[treatment.model], treatment)
  effect_direct <- post_coef[[coef_treatment]]
  
  # Mediator effect: coef(mediator) from model_y_treatment
  coef_mediator <- sprintf(pattern, response[treatment.model], mediator)
  effect_mediator <- post_coef[[coef_mediator]]
  
  # Indirect effect: coef(treament) from model_m_mediator * coef(mediator) from model_y_treatment
  coef_indirect <- sprintf(pattern, response[mediator.model], treatment)
  tmp.indirect <- post_coef[c(coef_indirect, coef_mediator)]
  effect_indirect <- tmp.indirect[[coef_indirect]] * tmp.indirect[[coef_mediator]]
  
  # Total effect
  effect_total <- effect_indirect + effect_direct
  # proportion mediated: indirect effect / total effect
  proportion_mediated <- effect_indirect/effect_total
  # Upsilon (from Lachowicz et al. 2018; http://dx.doi.org/10.1037/met0000165) - squared standardized effect
  # beta(mediator~x)^2 * beta(dv~mediator|x)^2
  upsilon <- effect_indirect^2
  
  # Add all of these into the posterior dataframe
  post <- data.frame(post_coef[,grep("^b_|^bsp",names(post_coef))],effect_direct,effect_indirect,effect_total,proportion_mediated,upsilon)
  n <- nrow(post)
  ci.lo <- (1-ci)/2
  ci.hi <- 1-(1-ci)/2
  
  # Get prior
  prior_slopes <- post_coef[,grep("^prior_b_",names(post_coef))]
  prior_slopes$prior_indirect <- prior_slopes[,1] * prior_slopes[,2]
  prior_slopes$prior_total <- prior_slopes$prior_indirect + prior_slopes[,2]
  prior_slopes$prior_proportion <- prior_slopes$prior_indirect/prior_slopes$prior_total
  # Remove the extreme values from the prior proportion distribution to stop density from not converging
  large_props <- which(abs(prior_slopes$prior_proportion) > 3)
  prior_slopes$prior_proportion[large_props] <- sample(prior_slopes$prior_proportion[-large_props],length(large_props),replace = T)
  
  outlist <- mclapply(1:ncol(post),function(i){
    name <- names(post)[i]
    if(grepl("_Intercept",name)){ # Intercepts
      post_summ(post[,i],ci=0.95,ci.type=ci.type,ROPE=c(NA,NA))
    } else if(any(c(grepl("^b_",name),grepl("^effect_direct",name)))){ # Betas (and direct effect)
      post_summ(post[,i],ci=0.95,ci.type=ci.type,ROPE=ROPE,prior = prior_slopes[,1])
    } else if(grepl("^effect_indirect",name)){ # Indirect Effect
      post_summ(post[,i],ci=0.95,ci.type=ci.type,ROPE=ROPE,prior = prior_slopes$prior_indirect)
    } else if(grepl("^effect_total",name)){ # Indirect Effect
      post_summ(post[,i],ci=0.95,ci.type=ci.type,ROPE=ROPE,prior = prior_slopes$prior_total)
    } else if(grepl("^proportion_mediated",name)){
      post_summ(post[,i],ci=0.95,ci.type=ci.type,ROPE=c(-Inf,0.5), prior = prior_slopes$prior_proportion)
    } else{ # Anything else
      post_summ(post[,i],ci=0.95,ci.type=ci.type,ROPE=c(NA,NA))
    }
  })
  names(outlist) <- gsub("b_","",colnames(post))
  out_table <- do.call(rbind,outlist)
  
  result <- list(model=model,formula=form,coefs=out_table,posterior=post,ci=ci,ci.type=ci.type,ROPE=ROPE,
                 samples=nrow(post),name.x=treatment,name.m=mediator,name.y=response_name_fixed,
                 imputed=imputeMissing,call=match.call())
  if(!is.null(impute.OOBerror)){
    result$impute.OOBerror <- impute.OOBerror
  }
  class(result) <- "bayes.mediation"
  
  # If desired, save result as a file
  if(!is.null(file)){
    if(!grepl("\\.[Rr][Dd][Ss]$",file)){ # whoops forgot extension!
      file <- paste0(file,".RDS")
    }
    saveRDS(result,file=file) # Saves an RDS file with chosen name
  }
  
  return(result)
}

# Pretty print function for bayes.mediate class
print.bayes.mediation <- function(x,dig=3){
  cat("Call:")
  print(x$call)
  cat(paste0("\nBayesian mediation analysis with Region of Practical Equivalence (ROPE) beta = [",paste(x$ROPE,collapse=", "),"]:"))
  cat(paste0("\n n = ",nrow(x$model$data)," complete cases"))
  cat(paste0("\n Formula:"))
  cat(paste0("\n  ",lapply(x$formula$forms,function(X){as.character(as.formula(X))})))
  
  ci.type.label <- ifelse(x$ci.type=="eti","Equal-tailed","Highest Density")
  cat(paste0("\n\n Estimates [",x$ci*100,"% ",ci.type.label," CrI] (based on ",x$samples," posterior samples):\n\n"))
  # Now re-format the coefficient matrix
  coef_mat <- data.frame(round(x$coefs[-grep("_Intercept",rownames(x$coefs)),c(4:14)],digits=dig))
  coef_mat[,5] <- as.numeric(scales::scientific(coef_mat[,5],digits = dig))
  coef_mat[,11] <- as.numeric(scales::scientific(coef_mat[,11],digits = dig))
  colnames(coef_mat)[10] <- "P>ROPE"
  print(as.matrix(coef_mat))
  
  cat("\n")
  direct_row <- grep("^effect_direct",rownames(x$coefs))
  indirect_row <- grep("^effect_indirect",rownames(x$coefs))
  total_row <- grep("^effect_total",rownames(x$coefs))
  prop_row <- grep("^proportion_mediated",rownames(x$coefs))
  ups_row <- grep("^upsilon",rownames(x$coefs))
  
  cat(paste0(" Direct effect of ",x$name.x," on ",x$name.y,":"))
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$coefs[direct_row,4],dig)," [",round(x$coefs[direct_row,5],dig),", ",
             round(x$coefs[direct_row,6],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$coefs[direct_row,7]>0.999,"> 0.999",
                                                         paste0("= ",signif(x$coefs[direct_row,7],dig))),
             " (",ifelse(x$coefs[direct_row,7]>0.999,">99.9",signif(x$coefs[direct_row,7]*100,dig)),
             "% chance that Direct Effect (c') ",ifelse(sign(x$coefs[direct_row,4])==1,">","<")," 0)"))
  cat(paste0("\n  Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$coefs[direct_row,8]>999 | x$coefs[direct_row,8] < 0.001,
                    format(x$coefs[direct_row,8],dig=3,scientific=T),
                    signif(x$coefs[direct_row,8],dig=3)),
             " (BF01 = ",ifelse(1/x$coefs[direct_row,8]>999 | 1/x$coefs[direct_row,8] < 0.001,
                                format(1/x$coefs[direct_row,8],dig=3,scientific=T),
                                signif(1/x$coefs[direct_row,8],dig=3)),") — ",
             bf_evidence(x$coefs[direct_row,8],model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n  Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$coefs[direct_row,14]>999 | x$coefs[direct_row,14] < 0.001,
                                                            format(x$coefs[direct_row,14],dig=3,scientific=T),
                                                            signif(x$coefs[direct_row,14],dig=3)),
             " (BF.inROPE = ",ifelse(1/x$coefs[direct_row,14]>999 | 1/x$coefs[direct_row,14] < 0.001,
                                     format(1/x$coefs[direct_row,14],dig=3,scientific=T),
                                     signif(1/x$coefs[direct_row,14],dig=3)),") — ",
             bf_evidence(x$coefs[direct_row,14],model = F))) # Bayes Factor Qualitative Interpretation
  
  cat(paste0("\n\n Indirect effect of ",x$name.x," on ",x$name.y," (via ",x$name.m,"):"))
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$coefs[indirect_row,4],dig)," [",round(x$coefs[indirect_row,5],dig),", ",
             round(x$coefs[indirect_row,6],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$coefs[indirect_row,7]>0.999,"> 0.999",
                                                         paste0("= ",signif(x$coefs[indirect_row,7],dig))),
             " (",ifelse(x$coefs[indirect_row,7]>0.999,">99.9",signif(x$coefs[indirect_row,7]*100,dig)),
             "% chance that Indirect Effect (a*b) ",ifelse(sign(x$coefs[indirect_row,4])==1,">","<")," 0)"))
  cat(paste0("\n  Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$coefs[indirect_row,8]>999 | x$coefs[indirect_row,8] < 0.001,
                    format(x$coefs[indirect_row,8],dig=3,scientific=T),
                    signif(x$coefs[indirect_row,8],dig=3)),
             " (BF01 = ",ifelse(1/x$coefs[indirect_row,8]>999 | 1/x$coefs[indirect_row,8] < 0.001,
                                format(1/x$coefs[indirect_row,8],dig=3,scientific=T),
                                signif(1/x$coefs[indirect_row,8],dig=3)),") — ",
             bf_evidence(x$coefs[indirect_row,8],model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n  Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$coefs[indirect_row,14]>999 | x$coefs[indirect_row,14] < 0.001,
                                                            format(x$coefs[indirect_row,14],dig=3,scientific=T),
                                                            signif(x$coefs[indirect_row,14],dig=3)),
             " (BF.inROPE = ",ifelse(1/x$coefs[indirect_row,14]>999 | 1/x$coefs[indirect_row,14] < 0.001,
                                     format(1/x$coefs[indirect_row,14],dig=3,scientific=T),
                                     signif(1/x$coefs[indirect_row,14],dig=3)),") — ",
             bf_evidence(x$coefs[indirect_row,14],model = F))) # Bayes Factor Qualitative Interpretation
  
  
  cat(paste0("\n\n Total effect of ",x$name.x," on ",x$name.y," (both direct and via ",x$name.m,"):"))
  cat(paste0("\n  Mdn [",x$ci*100,"% CrI] = ",round(x$coefs[total_row,4],dig)," [",round(x$coefs[total_row,5],dig),", ",
             round(x$coefs[total_row,6],dig),"]"))
  cat(paste0("\n  Probability of Direction (Pd) ",ifelse(x$coefs[total_row,7]>0.999,"> 0.999",
                                                         paste0("= ",signif(x$coefs[total_row,7],dig))),
             " (",ifelse(x$coefs[total_row,7]>0.999,">99.9",signif(x$coefs[total_row,7]*100,dig)),
             "% chance that Total Effect (c') ",ifelse(sign(x$coefs[total_row,4])==1,">","<")," 0)"))
  cat(paste0("\n  Bayes Factor vs. Point Null (BF10): ",
             ifelse(x$coefs[total_row,8]>999 | x$coefs[total_row,8] < 0.001,
                    format(x$coefs[total_row,8],dig=3,scientific=T),
                    signif(x$coefs[total_row,8],dig=3)),
             " (BF01 = ",ifelse(1/x$coefs[total_row,8]>999 | 1/x$coefs[total_row,8] < 0.001,
                                format(1/x$coefs[total_row,8],dig=3,scientific=T),
                                signif(1/x$coefs[total_row,8],dig=3)),") — ",
             bf_evidence(x$coefs[total_row,8],model = F))) # Bayes Factor Qualitative Interpretation
  cat(paste0("\n  Bayes Factor vs. ROPE (BF.ROPE): ",ifelse(x$coefs[total_row,14]>999 | x$coefs[total_row,14] < 0.001,
                                                            format(x$coefs[total_row,14],dig=3,scientific=T),
                                                            signif(x$coefs[total_row,14],dig=3)),
             " (BF.inROPE = ",ifelse(1/x$coefs[total_row,14]>999 | 1/x$coefs[total_row,14] < 0.001,
                                     format(1/x$coefs[total_row,14],dig=3,scientific=T),
                                     signif(1/x$coefs[total_row,14],dig=3)),") — ",
             bf_evidence(x$coefs[total_row,14],model = F))) # Bayes Factor Qualitative Interpretation
  
  cat(paste0("\n\n Percentage Mediated: ",100*round(x$coefs[prop_row,4],dig),"% [",100*round(x$coefs[prop_row,5],dig),"%, ",
             100*round(x$coefs[prop_row,6],dig),"%]"))
  cat(paste0("  |  Upsilon (R^2 analog): ",round(x$coefs[ups_row,4],dig)," [",round(x$coefs[ups_row,5],dig),", ",
             round(x$coefs[ups_row,6],dig),"]"))
}

# cor.bayesboot
# Zack Williams 04/18/2022, updated 1/16/23
# wrapper function for BBcor::bbcor to easily calculate correlations between 2 vectors
# Now includes polyserial/mixed correlations with the inclusion of polycor_auto
# Also tests correlations against a ROPE (default ±0.2) and presents median + 95% HDI
cor.bb <- cor.BB <- cor.bayesboot <- function(x,y=NULL,type="auto",B=5000,ROPE=c(-0.2,0.2),HDI=TRUE,conf=0.95,seed=12345,cores=detectCores()-1,ordinalLevelMax=11,...){
  # Take vector input and turn into data frame by combining x/y
  if(is.vector(x)){
    # Throw errors if data input incorrect
    if(is.null(y)){stop("Vector input detected for 'x' but 'y' is null! Provide a vector of equal length to 'x' for the 'y' argument.")}
    if(!is.vector(y) | length(y) != length(x)){stop("Invalid input for 'y' - must be a vector of equal length to 'x'.")}
    # Bind x and y together
    cordat <- cbind("VAR1"=x,"VAR2"=y)
  } else { # x is the wide-form data matrix
    if(!(is.matrix(x) | is.data.frame(x)) | ncol(x) < 2){stop("Invalid format for 'x'! Please provide a *wide-form* dataframe with each variable in a different column.")}
    cordat <- x
  }
  
  ## Cannibalizing BBcor::bbcor function to allow for ... input (e.g., additional input to psych::polychoric)
  
  ## Utility Functions
  bb_weights <- function(n){
    wts <- stats::rgamma(n, 1, 1)
    wts <- wts / sum(wts)
    return(wts)
  }
  normalize <- function(data, n) {
    x_ranked <- rank(data, ties.method = "random")
    x <- x_ranked / (n + 1)
    normalized_data <- stats::qnorm(x)
    return(normalized_data)
  }
  
  # polycor.auto - similar to psych mixed cor and lavaan cor auto; allows for sampling weights (for use in Bayesian bootstrap)
  # Written by Zack Williams 01/16/2023
  # Defaults to polychoric correlations for 10 or fewer categories
  # Uses wCorr for weighted polyserials, pearsons, polychorics
  polycor_auto <- function(data,ordinalLevelMax=11,weight=NULL){
    data <- apply(data,2,as.numeric) # convert to matrix and turn ordered/logical/factor vars to numerics
    data <- data[,which(!apply(data,2,function(x){all(is.na(x))}))] # Drop columns that couldn't be converted to numeric (e.g., characters)
    vnames <- colnames(data)
    nLevels <- apply(data,2,function(X){length(unique(na.omit(X)))})
    cont <- which(nLevels > ordinalLevelMax)
    ord <- which(nLevels >= 2 & nLevels <= ordinalLevelMax)
    ncase <- nrow(data)
    nvar <- ncol(data)
    
    # Calculate correlation matrix
    cormat <- sapply(1:nvar,function(i){
      sapply(1:nvar,function(j){
        temp <- data[,c(i,j)] 
        pairwise.complete <- complete.cases(temp)
        temp <- temp[pairwise.complete,]
        wts <- ncase*weight[pairwise.complete]
        if(j > i){ # Don't calculate half of the values to save time
          return(0)
        } else if(i==j){ # Same variable
          return(0.5) # Diagonal correlations will be doubled when matrix is folded
        } else if(i %in% ord & j %in% ord){ # Both ordinal
          return(wCorr::weightedCorr(temp[,1],temp[,2],method="Polychoric",weights=wts))
        } else if((i %in% ord & j %in% cont) | (i %in% cont & j %in% ord)){ # One is continuous
          # Figure out which one is continuous and put that as "y"
          if(i %in% cont){return(wCorr::weightedCorr(temp[,1],temp[,2],method="Polyserial",weights=wts))} # i continuous
          if(j %in% cont){return(wCorr::weightedCorr(temp[,2],temp[,1],method="Polyserial",weights=wts))} # j continuous
        } else{ # Both continuous
          return(wCorr::weightedCorr(temp[,1],temp[,2],method="Pearson",weights=wts))
        }
      })
    })
    # Fold over correlation matrix (double)
    cormat <- cormat + t(cormat)
    # Smooth correlation matrix if necessary
    if(smooth & !LaplacesDemon::is.positive.definite(cormat)){
      cormat <- as.matrix(Matrix::nearPD(cormat,corr = T)$mat)
    }
    
    rownames(cormat) <- colnames(cormat) <- vnames
    
    return(cormat)
  }
  
  bbcor <- function(x, 
                    method = "auto", 
                    iter = 1000, 
                    ordinalLevelMax=11,
                    cores = 2){
    
    # data matrix
    # (for return)
    Y <- x
    # na.omit
    x <- stats::na.omit(x)
    # variables
    p <- ncol(x)
    # observations
    n <- nrow(x)
    # parallel computing
    # cl <- parallel::makeCluster(cores)
    # redundancy for efficiency
    if ( method == "pearson" ) {
      # scale
      x <- scale(x)
      
      # draw from posterior
      samps  <- pbapply::pbreplicate(n = iter,
                                     stats::cov.wt(x,
                                                   wt = bb_weights(n),
                                                   cor = TRUE)$cor,
                                     cl = cores)
      
    } else if ( method == "spearman" ) {
      
      # ranks (makes sampling faster)
      x <- sapply(1:p, function(i) rank(x[, i]))
      
      # draw from posterior
      samps  <- pbapply::pbreplicate(n = iter,
                                     stats::cov.wt(x,
                                                   wt = bb_weights(n),
                                                   cor = TRUE)$cor,
                                     cl = cores)
      
    } else if ( method %in% c("polychoric","polyserial","poly","auto","mixed")) {
      # draw from posterior
      samps  <- pbapply::pbreplicate(n = iter,
                                     polycor_auto(x, weight = bb_weights(n),ordinalLevelMax=ordinalLevelMax),
                                     cl = cores)
    } else if ( method == "gaussian_rank" ) {
      
      # normalized ranks
      x <- sapply(1:p, function(i) normalize(x[, i], n))
      
      # draw from posterior
      samps  <- pbapply::pbreplicate(n = iter,
                                     stats::cov.wt(x,
                                                   wt = bb_weights(n),
                                                   cor = TRUE)$cor,
                                     cl = cores)
    } else {
      
      # draw from posterior  
      samps <- pbapply::pbreplicate(n = iter,
                                    wdm::wdm(x, method = method, 
                                             weights = bb_weights(n)), 
                                    cl = cores)
    }
    
    # Calculated quantities from correlation posteriors (note they're caught with capture.output so the progress bars don't show up)
    trash <- utils::capture.output(cor_mean <- pbapply::pbapply(samps, MARGIN = 1:2,cl = cores, FUN = mean)) # compute mean
    trash <- utils::capture.output(cor_median <- pbapply::pbapply(samps, MARGIN = 1:2,cl = cores, FUN = median)) # compute median
    
    if(HDI){
      trash <- utils::capture.output(cor_CrI_low <- pbapply::pbapply(samps, MARGIN = 1:2,cl = cores, FUN = function(X){ggdist::hdi(X,.width = conf)[1]})) # HDI.low
      trash <- utils::capture.output(cor_CrI_high <- pbapply::pbapply(samps, MARGIN = 1:2,cl = cores, FUN = function(X){ggdist::hdi(X,.width = conf)[2]})) # HDI.low
    } else{
      trash <- utils::capture.output(cor_CrI_low <- pbapply::pbapply(samps, MARGIN = 1:2,cl = cores, FUN = function(X){quantile(X,(1-conf)/2)})) # ETI.low
      trash <- utils::capture.output(cor_CrI_high <- pbapply::pbapply(samps, MARGIN = 1:2,cl = cores, FUN = function(X){quantile(X,1-(1-conf)/2)})) # ETI.low
    }
    
    # Smooth the summary correlation matrices to allow them to all be positive definite
    if(!is.positive.definite(cor_mean)){cor_mean <- as.matrix(Matrix::nearPD(cor_mean,corr = T)$mat)}
    if(!is.positive.definite(cor_median)){cor_median <- as.matrix(Matrix::nearPD(cor_median,corr = T)$mat)}
    if(!is.positive.definite(cor_CrI_low)){cor_CrI_low <- as.matrix(Matrix::nearPD(cor_CrI_low,corr = T)$mat)}
    if(!is.positive.definite(cor_CrI_high)){cor_CrI_high <- as.matrix(Matrix::nearPD(cor_CrI_high,corr = T)$mat)}
    
    # returned object
    returned_object <- list(cor_mean = cor_mean, 
                            cor_median = cor_median, 
                            cor_CrI = list("CrI_low"=cor_CrI_low,"CrI_high"=cor_CrI_high),
                            samps = samps, 
                            method = method, 
                            iter = iter, 
                            Y = Y)
    
    class(returned_object) <- c("bbcor", 
                                "default")
    return(returned_object)
  }
  
  # Back to function at hand
  set.seed(seed) # Set seed for reproducibility
  cor_out <- bbcor(cordat,method=type,iter=B,cores=cores,ordinalLevelMax=ordinalLevelMax,...) # bbcor function actually calculates the posterior samples
  
  post <- BBcor::posterior_samples(cor_out)
  cor_summary <- apply(post,2,post_summ,ci=conf,ci.type = ifelse(HDI,"hdi","eti"),ROPE = ROPE)[c(4:7,10:13),]
  
  result <- list("Summary"=cor_summary,"cor.type"=type,"n.subjects"=nrow(na.omit(cordat)),"n.var"=ncol(cordat),"BBcor"=cor_out,"Mean.Cor"=cor_out$cor_mean,"Median.Cor"=cor_out$cor_median,
                 "Cormat.CrI"=cor_out$cor_CrI,"posterior"=post,"ci.size"=conf,"ci.type"=ifelse(HDI,"hdi","eti"),"ROPE"=ROPE,"n.boot"=B,"data"=cordat,"Call"=match.call())
  class(result) <- "cor.bayesboot"
  return(result)
}

print.cor.bayesboot <- function(x,digits=3){
  cat("Call: ")
  print(x$Call)
  cat("\n")
  cat("Number of subjects =", x$n.subjects,"| Correlation type =",x$cor.type,"|",x$n.boot,"bootstrapped posterior samples\n\n")
  cat("Posterior Mean Correlation:\n")
  print(round(x$Mean.Cor,digits))
  cat("\n")
  cat(paste0(100*x$ci.size,"% ",ifelse(x$ci.type=="hdi","Highest-density","Equal-tailed")," Credible Interval (CrI):\n Lower Bound:\n"))
  print(round(x$Cormat.CrI[[1]],digits))
  cat("\n Upper Bound:\n")
  print(round(x$Cormat.CrI[[2]],digits))
  cat(paste0("\nSummary Indices (","based on ROPE of [",paste(x$ROPE,collapse=", "),"])"))
  cat("\n")
  print(round(x$Summary[-c(1:3,5,6),],digits))
}


### postprob.benefit.harm - Function to calculate the probability that a given heterogeneous effect will have a practically significant effect for a randomly-selected individual in the sample
### Zack Williams, 11/23/2022
# MOD: a brms model with a fixed effect and random slope ('varname')
# benefit.positive: logical that determines whether positive ES should be interpreted as a benefit or not
# ROPE region (default is [-0.2, 0.2])
# MCID = minimal clinically important difference (default = ±0.5 depending on benefit.positive)
# cores = number of cores for parallel computing
postprob.benefit.harm <- function(MOD,varname,benefit.positive=TRUE,ROPE=c(-0.2,0.2),MCID=ifelse(benefit.positive,0.5,-0.5),cores=1){
  
  # Correct mismatch between benefit.positive and MCID (sign of MCID must be in direction of benefit)
  if(MCID < 0 & benefit.positive){
    benefit.positive <- FALSE
  } else if(MCID > 0 & !benefit.positive){
    benefit.positive <- TRUE
  }
  
  # use posterior package to get the posterior distribution of MOD
  post <- as_draws_df(MOD)
  # Throw errors for incompatible models
  ncases <- length(grep(paste0("r_.*\\[.*,",varname),names(post)))
  if(ncases==0){stop(paste0("Model does not have a random slope defined for variable '",varname,"'. Check the spelling of varname and model specification."))}
  if(length(grep(paste0("^b_",varname),names(post)))==0){stop(paste0("Model does not have a fixed effect defined for variable '",varname,"'. Check the spelling of varname and model specification."))}
  
  # make "indivs," a dataframe of the overall effect plus the disturbances from each individual's random slope term (at each posterior draw)
  indivs <- suppressWarnings(post[,grep(paste0("r_.*\\[.*,",varname),names(post))]) + suppressWarnings(rep.col(post[,grep(paste0("^b_",varname),names(post))[1]],ncases))
  
  suppressWarnings(post[,grep(paste0("^b_",varname),names(post))[1]])
  
  # Posterior summary of fixed effect
  FE_summary <- suppressWarnings(post_summ(unlist(post[,grep(paste0("^b_",varname),names(post))[1]]),ROPE=ROPE)[c(4:7,12:13)])
  # Add P(exceeds MCID) to FE_summary
  p.exceeds.mcid <- suppressWarnings(ifelse(benefit.positive,
                                            mean(unlist(post[,grep(paste0("^b_",varname),names(post))[1]]) > MCID),
                                            mean(unlist(post[,grep(paste0("^b_",varname),names(post))[1]]) < MCID)))
  FE_summary <- c(FE_summary,"P>MCID"=p.exceeds.mcid)
  
  # Variance component of random effect
  RE_summary <- post_summ(apply(indivs,1,sd))[4:6]
  
  # Define the cutoffs from the ROPE as the sign that matches the MCID (default is positive)
  ben_cut <- ROPE[which(sign(ROPE)==sign(MCID))]
  harm_cut <- ROPE[which(sign(ROPE)!=sign(MCID))]
  
  # Meat of function: Calculate individual-level probabilities of benefit/harm/null-effect/above MCID for each case
  if(benefit.positive){ # If positive benefit
    p_benefit <- pbapply(indivs,1,function(X){mean(X > ben_cut)},cl=cores)
    p_harm <- pbapply(indivs,1,function(X){mean(X < harm_cut)},cl=cores)
    p_null <- pbapply(indivs,1,function(X){mean(X > harm_cut & X < ben_cut)},cl=cores)
    p_mcid <- pbapply(indivs,1,function(X){mean(X > MCID)},cl=cores)
  } else{ # Benefit is negative
    p_benefit <- pbapply(indivs,1,function(X){mean(X < ben_cut)},cl=cores)
    p_harm <- pbapply(indivs,1,function(X){mean(X > harm_cut)},cl=cores)
    p_null <- pbapply(indivs,1,function(X){mean(X < harm_cut & X > ben_cut)},cl=cores)
    p_mcid <- pbapply(indivs,1,function(X){mean(X < MCID)},cl=cores)
  }
  
  # Utility function: summarize posterior distribution of probability (Mean, SD, quantiles)
  SUMMARY <- function(X){c("M"=mean(X),"SD"=sd(X),quantile(X,c(0.025,0.25,0.5,0.75,0.95,0.975)),"Min"=min(X),"Max"=max(X))}
  
  result <- list("P.benefit"=SUMMARY(p_benefit),"P.harm"=SUMMARY(p_harm),"P.null"=SUMMARY(p_null),"P.MCID"=SUMMARY(p_mcid),
                 "Summary.FixEff"=FE_summary,"Summary.RanEff"=RE_summary,"PARAM"=varname,"OUTCOME"=lhs(as.formula(MOD$formula)),"formula"=as.formula(MOD$formula),
                 "N.obs"=nrow(MOD$data),"N.cases"=ncases,"ROPE"=ROPE,"MCID"=MCID,"posterior.effects"=indivs,"posterior"=post,"Call"=match.call())
  
  class(result) <- "bayes.BHpredict"
  return(result)
}

### Pretty-print function for output of postprob.benefit.harm
print.bayes.BHpredict <- function(x,dig=3){
  cat("Call:")
  print(x$Call)
  cat(paste0("\nModel formula: ",x$formula))
  cat(paste0("\n\nBayesian analysis of posterior probabilities of benefit (P.ben), harm (P.harm), trivially small effect (P.null), and a clinically significant effect (P.mcid)"))
  cat(paste0("\n n = ",x$N.obs," observations from ",x$N.cases," cases/groups"))
  cat(paste0("\n Based on ROPE bounds of [",paste(x$ROPE,collapse=", "),"] and Minimal Clinically Important Difference (MCID) of ",x$MCID,":"))
  
  cat(paste0("\n\n  Fixed effect of ",x$PARAM," on ",x$OUTCOME,":\n"))
  print(round(x$Summary.FixEff,dig))
  cat(paste0("\n  Random effect (SD) of ",x$PARAM," across cases/groups:\n"))
  print(round(x$Summary.RanEff,dig))
  
  cat(paste0("\nPosterior probabilities (M±SE) and 95% HDIs (based on ",nrow(x$posterior)," posterior samples):"))
  
  if(x$MCID < 0){ # Negative benefit
    cat(paste0("\n Probability of Benefit [P(D < ",x$ROPE[1],"|data)]: ", format(x$P.benefit[1],nsmall=dig,digits=dig)," ± ",format(x$P.benefit[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.benefit[3],nsmall=dig,digits=dig),", ",format(x$P.benefit[8],nsmall=dig,digits=dig),"]"))
    cat(paste0("\n Probability of Harm [P(D > ",x$ROPE[2],"|data)]: ", format(x$P.harm[1],nsmall=dig,digits=dig)," ± ",format(x$P.harm[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.harm[3],nsmall=dig,digits=dig),", ",format(x$P.harm[8],nsmall=dig,digits=dig),"]"))
    cat(paste0("\n Probability of Null Effect [P(D in [",x$ROPE[1],",",x$ROPE[2],"]","|data)]: ", format(x$P.null[1],nsmall=dig,digits=dig)," ± ",format(x$P.null[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.null[3],nsmall=dig,digits=dig),", ",format(x$P.null[8],nsmall=dig,digits=dig),"]"))
    cat(paste0("\n Probability of Clinically Significant Effect [P(D < ",x$MCID,"|data)]: ", format(x$P.MCID[1],nsmall=1,digits=dig)," ± ",format(x$P.MCID[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.MCID[3],nsmall=dig,digits=dig),", ",format(x$P.MCID[8],nsmall=dig,digits=dig),"]"))
  } else{ # Positive benefit (MCID > 0)
    cat(paste0("\n Probability of Benefit [P(D > ",x$ROPE[2],"|data)]: ", format(x$P.benefit[1],nsmall=dig,digits=dig)," ± ",format(x$P.benefit[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.benefit[3],nsmall=dig,digits=dig),", ",format(x$P.benefit[8],nsmall=dig,digits=dig),"]"))
    cat(paste0("\n Probability of Harm [P(D < ",x$ROPE[1],"|data)]: ", format(x$P.harm[1],nsmall=dig,digits=dig)," ± ",format(x$P.harm[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.harm[3],nsmall=dig,digits=dig),", ",format(x$P.harm[8],nsmall=dig,digits=dig),"]"))
    cat(paste0("\n Probability of Null Effect [P(D in [",x$ROPE[1],",",x$ROPE[2],"]","|data)]: ", format(x$P.null[1],nsmall=dig,digits=dig)," ± ",format(x$P.null[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.null[3],nsmall=dig,digits=dig),", ",format(x$P.null[8],nsmall=dig,digits=dig),"]"))
    cat(paste0("\n Probability of Clinically Significant Effect [P(D > ",x$MCID,"|data)]: ", format(x$P.MCID[1],nsmall=1,digits=dig)," ± ",format(x$P.MCID[2],nsmall=dig,digits=dig),", 95% HDI [", format(x$P.MCID[3],nsmall=dig,digits=dig),", ",format(x$P.MCID[8],nsmall=dig,digits=dig),"]"))
  }
}

# variance decomposition (from performance package) using a model fit to multiply-imputed data
# Zack Williams - 11/24/22
bayes.ICC.mi <- bayes.ICC <- function(MOD,data=NULL,re_formula=NULL,...,thin=10,digits=3,cores=detectCores()-1){
  # Determines if data is a multiply imputed dataframe and processes appropriiately
  if(family(MOD)$link == "identity"){ # If identity link, use the actual predictions
    if(is.data.frame(data) | is.null(data)){
      post_FE <- pbapply(posterior_predict(MOD,re_formula=NA,draw_ids=seq(1,ndraws(MOD),thin),newdata = data,cores=cores),1,var)
      post_RE <- pbapply(posterior_predict(MOD,re_formula=re_formula,draw_ids=seq(1,ndraws(MOD),thin),newdata = data,cores=cores),1,var)
    } else{ # List of dataframes (from multiple imputatation)
      post_FE <- unlist(pblapply(data,function(X){apply(posterior_predict(MOD,re_formula=NA,draw_ids=seq(1,ndraws(MOD),thin),newdata = X,cores=cores),1,var)}))
      post_RE <- unlist(pblapply(data,function(X){apply(posterior_predict(MOD,re_formula=re_formula,draw_ids=seq(1,ndraws(MOD),thin),newdata = X,cores=cores),1,var)}))
    }
  } else{ # Otherwise use the linear predictors
    if(is.data.frame(data) | is.null(data)){
      post_FE <- pbapply(posterior_linpred(MOD,re_formula=NA,draw_ids=seq(1,ndraws(MOD),thin),newdata = data,cores=cores),1,var)
      post_RE <- pbapply(posterior_linpred(MOD,re_formula=re_formula,draw_ids=seq(1,ndraws(MOD),thin),newdata = data,cores=cores),1,var)
    } else{ # List of dataframes (from multiple imputatation)
      post_FE <- unlist(pblapply(data,function(X){apply(posterior_linpred(MOD,re_formula=NA,draw_ids=seq(1,ndraws(MOD),thin),newdata = X,cores=cores),1,var)}))
      post_RE <- unlist(pblapply(data,function(X){apply(posterior_linpred(MOD,re_formula=re_formula,draw_ids=seq(1,ndraws(MOD),thin),newdata = X,cores=cores),1,var)}))
    }
  }
  
  # Variance Decomposition (Similar to ICC): (Total (RE + FE) variance - FE only variance) / Total variance
  vardecomp <- (post_RE - post_FE)/post_RE
  
  # Summarize the posterior of the ICC
  return(round(c("ICC"=post_summ(vardecomp,...)[4:7]),digits))
}

