##### Extra IRT functions for the mirt package
##### by Zack Williams, Vanderbilt University
##### Updated 05/12/2021
##### No guarantee that these functions will work in any way. Use at your own risk.
##### Feedback, comments, etc.? Contact Zack at zachary.j.williams@vanderbilt.edu
# Package loading using pacman
if(!require(pacman)){
  install.packages("pacman")
  require(pacman)
}
# Required packages: mirt, gtools, mvtnorm, parallel, tidyverse, bayestestR, glue, missForest, cubature, mvnfast
pacman::p_load(mirt,gtools,mvtnorm,parallel,bayestestR,glue,missForest,cubature,mvnfast)
require(tidyverse,quietly = T)

# Utility functions: rep.row, rep.col
rep.row<-function(x,n){ # used by DIF.wABC
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
# Quiet utility function - quiet return of object
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# loads() or loadings() function
# Utility function to get standardized loadings from lavaan or other object
# Now able to handle "lslx" objects (regularized SEM from lslx package)
loadings <- loads <- function(obj,bifactor=NULL,add.h2=T,cut=NULL,selector="bic",...){
  if("lavaan" %in% class(obj)){
    loadingmat <- lavInspect(obj,"std",...)$lambda
    Phi <- lavInspect(obj,"cor.lv",...)
    h2 <- diag(loadingmat %*% Phi %*% t(loadingmat))
    # figure out whether bifactor model is appropriate
    if(is.null(bifactor)){ 
      bifactor <- (all(abs(loadingmat[,1])>0.1) & ncol(loadingmat) > 1 & all(Phi == diag(nrow(Phi))))
    }
  } else if(any(attr(class(obj),"package")=="mirt")){ # mirt object
    loadingmat <- extract.mirt(obj,"F")
    h2 <- extract.mirt(obj,"h2")
    Phi <- diag(ncol(loadingmat))
    # Switch signs of loadings in matrix
    signed <- sign(colSums(loadingmat))
    signed[signed==0] <- 1
    if(ncol(loadingmat) > 1){
      cnames <- colnames(loadingmat)
      loadingmat <- loadingmat %*% diag(signed)  # flips factors to be in positive direction but loses the colnames
      colnames(loadingmat) <- cnames # put the names back
    } else{ # for 1-factor case
      loadingmat <- loadingmat * signed
    }
    if(is.null(bifactor)){ 
      bifactor <- (all(abs(loadingmat[,1])>0.1) & ncol(loadingmat) > 1 & all(Phi == diag(nrow(Phi))))
    }
  } else if("lslx" %in% class(obj)){
    loadingmat <- obj$extract_coefficient_matrix(selector=selector,block = "y<-f")$g
    Phi <- obj$extract_coefficient_matrix(selector=selector,block = "f<->f")$g
    h2 <- diag(loadingmat %*% Phi %*% t(loadingmat))
  } else{
    loadingmat <- stats::loadings(obj,...)
    Phi <- obj$Phi
    if(!is.null(Phi)){
      h2 <- diag(loadingmat %*% Phi %*% t(loadingmat))
    } else{
      h2 <- diag(loadingmat %*% t(loadingmat))
    }
  }
  # Add communalities
  if(add.h2){
    # figure out whether bifactor model is appropriate
    if(is.null(bifactor)){
      bifactor <- (all(abs(loadingmat[,1])>0.1)  & ncol(loadingmat) > 1 & all(Phi == diag(nrow(Phi))))
    }
    if(bifactor){IECV <- apply(loadingmat,1,function(x){x[1]^2/sum(x^2)})} else{IECV <- NULL}
    if(!is.null(cut)){loadingmat[which(abs(loadingmat) < cut)] <- NA}
    result <- cbind(loadingmat,"h2"=h2,"I-ECV"=IECV)
  } else{
    if(!is.null(cut)){loadingmat[which(abs(loadingmat) < cut)] <- NA}
    result <- loadingmat
  }
  
  class(result) <- c("lavaan.matrix","matrix")
  return(result)
}
# Wrapper function for factor correlation matrix extraction
Phi <- function(obj,plot=F,selector="bic",...){
  if("lavaan" %in% class(obj)){
    phi <- lavInspect(obj,"cor.lv")
  } else if(any(attr(class(obj),"package")=="mirt")){ # mirt model
    if("MultipleGroupClass" %in% class(obj)){ # Use phi from first group
      message("MultipleGroupClass detected. Using only factor correlations from first group.")
      phi <- quiet(summary(obj)[[1]]$fcor)
    } else{
      phi <- quiet(summary(obj)$fcor)
    }
    # if phi is a symmetric matrix, return it. Otherwise, use lower tri to make it symmetric
    if(isSymmetric(phi)){
      class(phi) <- c("lavaan.matrix.symmetric","matrix")
    } else{
      phi_lower <- phi[lower.tri(phi)] # take values below diagonal
      phi <- t(phi)
      phi[lower.tri(phi)] <- phi_lower
      class(phi) <- c("lavaan.matrix.symmetric","matrix")
    }
  } else if("lslx" %in% class(obj)){
    phi <- obj$extract_coefficient_matrix(selector=selector,block = "f<->f")$g
    class(phi) <- c("lavaan.matrix.symmetric","matrix")
  } else{
    if(is.matrix(obj) & is.null(obj$Phi)){return(NULL)}
    phi <- obj$Phi # psych::fa object or GPArotation object
    if(is.null(phi)){ # If no phi found, assume orthogonal
      ls <- loads(obj,add.h2 = F)
      phi <- diag(ncol(ls))
      rownames(phi) <- colnames(phi) <- colnames(ls)
    }
    class(phi) <- c("lavaan.matrix.symmetric","matrix")
  }
  # Plot correlations
  if(plot==T){
    psych::corPlot(phi,xlas=2,...)
  }
  return(phi)
}

# Harrell-Davis Quantile
hd <- function(x,qs=.5,na.rm=TRUE){
  
  elimna<-function(m){
    # remove any rows of data having missing values
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

# Wrapper function for factor correlation matrix extraction
Phi <- function(obj,plot=F,...){
  if("lavaan" %in% class(obj)){
    phi <- lavInspect(obj,"cor.lv")
  } else if(any(attr(class(obj),"package")=="mirt")){ # mirt model
    phi <- quiet(summary(obj)$fcor)
    # if phi is a symmetric matrix, return it. Otherwise, use lower tri to make it symmetric
    if(isSymmetric(phi)){
      class(phi) <- c("lavaan.matrix.symmetric","matrix")
    } else{
      phi_lower <- phi[lower.tri(phi)] # take values below diagonal
      phi <- t(phi)
      phi[lower.tri(phi)] <- phi_lower
      class(phi) <- c("lavaan.matrix.symmetric","matrix")
    }
  } else{
    phi <- obj$Phi # psych::fa object
    class(phi) <- c("lavaan.matrix.symmetric","matrix")
  }
  # Plot correlations
  if(plot==T){
    psych::corPlot(phi,xlas=2,...)
  }
  return(phi)
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

# Get Q3 residuals and examine relative to cutoff of 0.2 (wrapper for mirt::residuals)
# Zack Williams 04/08/2021
# Added option to return Q3* (Q3 values minus mean of Q3) using normalize=T
LD.Q3 <- function(object,Theta=fscores(object,QMC=extract.mirt(object,"nfact")>3),Q3.cut=0.2,round=3,pos.only=F,normalize=F){
  res <- quiet(residuals(object,type="Q3",Theta=Theta,QMC=extract.mirt(object,"nfact")>3))
  if(normalize){res <- res - mean(res[lower.tri(res)])}
  large_res <- corCheck(res,Q3.cut)
  
  inames <- extract.mirt(object,"itemnames")
  name_mat <- sapply(1:length(inames),function(i){
    sapply(1:length(inames),function(j){
      paste0(inames[i],"|",inames[j])
    })
  })
  Q3s <- res[lower.tri(res)]
  names(Q3s) <- name_mat[lower.tri(res)]
  if(pos.only){
    large_Q3s <- Q3s[which(Q3s>Q3.cut)]
  } else{
    large_Q3s <- Q3s[which(abs(Q3s)>Q3.cut)]
  }
  
  diag(res) <- NA
  Q3.Range <- range(Q3s)
  names(Q3.Range) <- c(paste0("Min:",names(Q3s)[which.min(Q3s)]),
                       paste0("Max:",names(Q3s)[which.max(Q3s)]))
  if(normalize){
    return(list("Q3*"=round(res,round),"Large.Q3*s"=round(large_Q3s,round),"Cutoff"=round(Q3.cut,round),"Q3*.Range"=round(Q3.Range,round)))
  } else{
    return(list("Q3"=round(res,round),"Large.Q3s"=round(large_Q3s,round),"Cutoff"=round(Q3.cut,round),"Q3.Range"=round(Q3.Range,round)))
  }
}

### LD.Q3.boot
### Q3-based LD detection using parametric bootstrap of maximum Q3 value proposed by Christensen et al. (2017)
### Zack Williams, 04/24/2021
## Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values for Yen’s Q3: 
## Identification of local dependence in the Rasch model using residual correlations. 
## Applied psychological measurement, 41(3), 178-194. https://doi.org/10.1177/0146621616677520
### obj = a mirt object
### q = quantile of the empirical Q3 distribution to use as critical value (default is 0.95, only pos LD flagged)
### n.boot = number of simulated datasets to generate (default is 1000)
### Theta = factor scores used to calculate expected item probabilities (default EAP scores from object with QMC used for nfactors > 3)
### cores = number of cores used for multi-core processing (default = 1 fewer than maximum number)
### seed = random seed for reproducibility (default = 12345)
### Nnormalize = should the Q3 values be tested using absolute deviations from mean Q3? (default = false)
LD.Q3.boot <- LD.Q3.pb <- function(obj,q=0.95,n.boot=1000,Theta=fscores(obj,QMC=extract.mirt(obj,"nfact")>3),normalize=F,
                                   round=3,cores=1,seed=12345){
  if(attr(class(obj),"package")!="mirt"){stop("ERROR: Input must be a mirt model!")}
  # Breaks if mirtCluster is defined: remove MirtCluster
  if(cores > 1){
    suppressMessages(mirtCluster(remove=T))
  }
  message("Simulating datasets...")
  sim <- pblapply(1:n.boot,function(i){
    set.seed(seed+i)
    simdata(model=obj,Theta = Theta)
  },cl=cores)
  ### Run bootstrap
  message("Bootstrapping Q3 distribution...")
  q3_boot <- pbsapply(sim,function(s){
    # Calculate item residuals for each item
    res <- sapply(1:extract.mirt(obj,"nitems"),function(i){
      ei <- extract.item(obj, item=i)
      EI <- expected.item(ei, Theta=Theta)
      res_i <- s[ ,i] - EI
    })
    # Q3 mat is matrix of correlations between residuals
    Q3s <- cor(res,use="pairwise",method="pearson")[lower.tri(cor(res))]
    if(normalize){Q3s <- abs(Q3s - mean(Q3s))}
    return(max(Q3s)) # return maximum value of Q3 mat for each iteration
  },cl=cores)
  result <- LD.Q3(obj,Theta = Theta,Q3.cut = quantile(q3_boot,q),pos.only = !normalize,round = round,normalize = normalize)
  result$Q3.bootdist <- round(quantile(q3_boot,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)),round)
  return(result)
}

# Standardized LD X2 (in the method of IRTPRO)
# Values > ±10 indicate LD, values between 5 and 10 are borderline
# Standardized X2 = (X2 - df)/sqrt(2*df) - see IRTPRO User Guide for more info
# Updated to now include Q3 and Q3* values for significant X2 values (sort of like an effect size)
LD.standard <- function(object,flag=10,round=2,QMC=extract.mirt(object,"nfact")>3,pos.only=F,add.Q3=T,...){
  inames <- extract.mirt(object,"itemnames")
  res <- quiet(residuals(object, type="LD",df.p=T,QMC=QMC,...))
  if(extract.mirt(object,"ngroups")==1){
    if(add.Q3){
      Q3_mat <- quiet(residuals(object, type="Q3",Theta=fscores(object,QMC=QMC),...))
    }
    x2s <- abs(res$LD[lower.tri(res$LD)])
    dfs <- res$df.p[lower.tri(res$df.p)]
    signs <- sign(t(res$LD)[lower.tri(res$LD)])
    std_x2s <- round(signs*(x2s-dfs)/sqrt(2*dfs),round)
    
    result <- round(res$LD,round)
    class(result) <- "matrix"
    result[lower.tri(result)] <- std_x2s
  } else{ # Multiple Groups
    if(add.Q3){
      Q3_mat <- quiet(residuals(object, type="Q3",...))
    }
    result <- lapply(res,function(resi){
      x2s <- abs(resi$LD[lower.tri(resi$LD)])
      dfs <- resi$df.p[lower.tri(resi$df.p)]
      signs <- sign(t(resi$LD)[lower.tri(resi$LD)])
      std_x2s <- round(signs*(x2s-dfs)/sqrt(2*dfs),round)
      
      result <- round(resi$LD,round)
      class(result) <- "matrix"
      result[lower.tri(result)] <- std_x2s
      return(result)
    })
  }
  
  # Matrix of names (for vector of item pairs)
  name_mat <- sapply(1:length(inames),function(i){
    sapply(1:length(inames),function(j){
      paste0(inames[i],"|",inames[j])
    })
  })
  
  
  if(extract.mirt(object,"ngroups")==1){
    X2s <- result[lower.tri(result)]
    if(add.Q3){
      Q3s <- round(Q3_mat[lower.tri(result)],round)
      Q3stars <- round(Q3_mat[lower.tri(result)] - mean(Q3_mat[lower.tri(result)]),round)
      names(Q3s) <- names(Q3stars) <- name_mat[lower.tri(result)]
    }
    CramerV <- t(result)[lower.tri(result)]
    resids <- residCheck(object,cut = 0)$Large.Resids
    names(X2s) <- names(CramerV) <- name_mat[lower.tri(result)]
    if(pos.only){ # only flag positive LD
      large_X2 <- which(X2s > flag) # Which X2s > 10
    } else{
      large_X2 <- which(abs(X2s) > flag) # Which X2s > ±10 (default)
    }
    
    if(length(large_X2)>0){
      LDvals <- cbind("Std.X2"=X2s[large_X2],"Cramer.V"=CramerV[large_X2],"Resid"=resids[large_X2])
      if(add.Q3){
        LDvals <- cbind(LDvals,"Q3*"=Q3stars[large_X2],"Q3"=Q3s[large_X2])
      }
    } else{
      LDvals <- integer(0)
    }
    
    X2.Range <- range(X2s)
    names(X2.Range) <- c(paste0("Min:",names(X2s)[which.min(X2s)]),
                         paste0("Max:",names(X2s)[which.max(X2s)]))
    
  } else{ # Multiple Groups
    LDvals <- lapply(1:length(result),function(i){
      resulti <- result[[i]]
      X2s <- resulti[lower.tri(resulti)]
      CramerV <- t(resulti)[lower.tri(resulti)]
      if(add.Q3){
        Q3i <- round(Q3_mat[[i]],round)
        Q3s <- Q3i[lower.tri(resulti)]
        Q3stars <- round(Q3i[lower.tri(resulti)] - mean(Q3i[lower.tri(resulti)]),round)
        names(Q3s) <- names(Q3stars) <- name_mat[lower.tri(resulti)]
      }
      names(X2s) <- names(CramerV) <- name_mat[lower.tri(resulti)]
      if(pos.only){ # only flag positive LD
        large_X2 <- which(X2s > flag) # Which X2s > 10
      } else{
        large_X2 <- which(abs(X2s) > flag) # Which X2s > ±10 (default)
      }
      
      if(length(large_X2)>0){
        LDvals <- cbind("Std.X2"=X2s[large_X2],"Cramer.V"=CramerV[large_X2])
        if(add.Q3){
          LDvals <- cbind(LDvals,"Q3*"=Q3stars[large_X2],"Q3"=Q3s[large_X2])
        }
      } else{
        LDvals <- integer(0)
      }
      return(LDvals)
    })
    
    X2.Range <- lapply(result,function(resulti){
      X2s <- resulti[lower.tri(resulti)]
      names(X2s) <- name_mat[lower.tri(resulti)]
      X2.Range <- range(X2s)
      names(X2.Range) <- c(paste0("Min:",names(X2s)[which.min(X2s)]),
                           paste0("Max:",names(X2s)[which.max(X2s)]))
      return(X2.Range)
    })
  }
  return(list("LD"=result,"LDvals"=LDvals,"X2.Range"=X2.Range))
}

# C2 - assess the C2 fit index of Cai & Monroe. Wrapper for mirt::M2 with different options.
# na.rm is always on, QMC is chosen by number of factors, and C2 is default type
C2 <- function(object,type="C2",QMC=extract.mirt(object,"nfact")>3,round=3,...){
  if(!type %in% c("C2","M2","M2*")){stop("'type' not recognized. Valid values are C2, M2, or M2*")}
  if(type=="C2"){suppressMessages(mirtCluster(detectCores()-1))} # Run C2 in parallel
  
  out <- try(mirt::M2(object,na.rm=T,QMC=QMC,type=type,residmat=F,...),silent = T)
  if(any(class(out)=="try-error")){stop(out)} # Catch errors including too few DF for M2*
  if(!is.na(round)){
    out <- round(out,round) # can return unrounded by round = NA
  }
  
  colnames(out)[1] <- type
  colnames(out)[-c(1:3)] <- paste0(colnames(out)[-c(1:3)],".",type)
  colnames(out) <- gsub("SRMSR","SRMR",colnames(out))
  return(unlist(out))
}

# itemfit function that overrides the mirt itemfit, adding adjusted p-vals and having na.rm by default
itemfit <- function(object,p.adj="BH",cutoff=0.05,QMC=extract.mirt(object,"nfact")>3,na.rm=T,...){
  ifit <- mirt::itemfit(object,QMC=QMC,na.rm=na.rm,...)
  # If there's a p-value, adjust it.
  if(length(grep("^p\\.",names(ifit)))>0){
    ifit$p.adj <- p.adjust(ifit[,grep("^p\\.",names(ifit))[1]],method = p.adj)
    ifit$sig <- ""
    ifit$sig[which(ifit[,grep("^p\\.",names(ifit))[1]] < cutoff)] <- "."
    ifit$sig[which(ifit$p.adj < cutoff)] <- "*"
  }
  return(ifit)
}

#### mirt.2step / multipleGroup.2step
### Zack Williams, 05/07/2021
### Estimates mirt model in two steps
## First step uses MHRM algorithm to get close to MLE, then second step uses params from that as input to ML fitting (default QMCEM)
mirt.2step <- function(...,method="QMCEM",SE.type="Oakes",SE=T,TOL=NULL,TOL.mhrm=0.001){
  mod_mhrm <- mirt(...,method="MHRM",SE=F,TOL=TOL.mhrm)
  mhrm_vals <- mod2values(mod_mhrm)
  mod_ml <- mirt(...,method = method,pars=mhrm_vals,SE=T,SE.type = SE.type,TOL=TOL)
  return(mod_ml)
}
# Same thing for multipleGroup
multipleGroup.2step <- function(...,method="QMCEM",SE.type="Oakes",SE=T,TOL=NULL,TOL.mhrm=0.001){
  mod_mhrm <- multipleGroup(...,method="MHRM",SE=F,TOL=TOL.mhrm)
  mhrm_vals <- mod2values(mod_mhrm)
  mod_ml <- multipleGroup(...,method = method,pars=mhrm_vals,SE=T,SE.type = SE.type,TOL=TOL)
  return(mod_ml)
}

#### Effect sizes for differential item functioning
#### Zack Williams
#### Updated 12/18/2021
# Based on the taxonomy proposed by Meade, 2010 (http://doi.org/10.1037/a0018966)
# Also note the erratum to the above paper which presents a different standard deviation term for ESSD (http://doi.org/10.1037/a0020897)
# Works for any number of groups, although the effect sizes are based on each 
# focal group compared to a reference group (default reference being group 1)
# Note that this function is similar to mirt's empirical_ES function (but works for multidimensional models)
# Now allows for individual-level ESs as of December 2021
DIF.ES <- function(difMod,ref.grp=1,which.items=1:extract.mirt(difMod,"nitems"),digits=3){
  # Item names
  inames <- extract.mirt(difMod,"itemnames")
  which.item.names <- inames[which.items]
  # Define the groups and which group is reference group (default, first one by factor level)
  grps <- extract.mirt(difMod,"groupNames")
  ngrps <- length(grps)
  grp_assign <- extract.mirt(difMod,"group")
  if(class(ref.grp)=="character"){
    ref.grp <- which(grps==ref.grp) # Change to numeric
    if(length(ref.grp)==0){stop("ERROR: Invalid reference group value")}
  }
  # Focal groups are all groups other than reference
  focal_grps <- which(grps!=grps[ref.grp])
  # Get reference group model
  mod_Ref <- suppressMessages(extract.group(difMod,ref.grp))
  # Theta scores from DIF model
  theta_mat <- fscores(difMod,QMC=extract.mirt(difMod,"nfact")>3)
  # Probability traces for all groups
  prob.traces <- lapply(seq_len(length(grps)),function(grp){
    mod <- suppressMessages(extract.group(difMod,grp))
    return(probtrace(mod,theta_mat))
  })
  # What values are the items taking in each column? (To calculate estimated score)
  item_scores <- as.numeric(gsub(".*P\\.(\\d)$","\\1",colnames(prob.traces[[1]])))
  # Multiply the probability * the item score
  Trace.Vals <- mclapply(prob.traces,function(tr){
    sapply(1:ncol(tr),function(i){
      tr[,i] <- tr[,i]*item_scores[i]
    })
  })
  # Expected scores based on probability traces for each group
  Escores <- mclapply(which.item.names,function(i){
    which.cols <- grep(paste0("^",i,"\\.P"),colnames(prob.traces[[1]]))
    # Bind together the item scores across groups 
    grp.iscores <- sapply(Trace.Vals,function(val){
      rowSums(val[,which.cols])
    })
    colnames(grp.iscores) <- grps
    return(grp.iscores)
  })
  names(Escores) <- which.item.names
  
  ES_list <- lapply(focal_grps,function(focal.grp){
    focal_Ss <- which(grp_assign==grps[focal.grp])
    Nf <- length(focal_Ss) # Number of cases in focal group
    # Expected test scores for focal group based on ref/foc parameters
    ETscores <- cbind("ET.ref"=rowSums(Trace.Vals[[ref.grp]]),
                      "ET.foc"=rowSums(Trace.Vals[[focal.grp]]))[focal_Ss,]
    
    # Item-level indices
    # SIDS: signed difference in expected item score for the sample
    # UIDS: unsigned (absolute) difference in expected item score for the sample
    # D-Max: maximum difference in expected item score for any one person in sample
    # ESSD: Standardized mean item score difference (Cohen's d) between groups
    item.indices <- t(sapply(Escores,function(ES){
      ES2 <- ES[focal_Ss,focal.grp]
      ES1 <- ES[focal_Ss,ref.grp]
      SDipool <- sqrt(((Nf-1)*var(ES2) + (Nf-1)*var(ES1))/(2*Nf-2)) # Item SD pooled - from Meade erratum (doi: 10.1037/a0020897)
      ESSD <- (mean(ES2) - mean(ES1))/SDipool
      ES_diffs <- ES2 - ES1
      ABS_diffs <- abs(ES_diffs)
      round(c("SIDS"=mean(ES_diffs),"UIDS"=mean(ABS_diffs),
              "Dmax"=ES_diffs[which.max(ABS_diffs)],"ESSD"=ESSD),digits)
    }))
    
    ### Individual-level item indices
    # SIDi: signed difference in expected item score for the individual
    # UIDi: unsigned (absolute) difference in expected item score for the individual
    # ESSDi: Standardized item score difference for an individual in SD (Cohen's d) units
    item.indices.individual <- lapply(Escores,function(ES){
      ES2 <- ES[focal_Ss,focal.grp]
      ES1 <- ES[focal_Ss,ref.grp]
      SDipool <- sqrt(((Nf-1)*var(ES2) + (Nf-1)*var(ES1))/(2*Nf-2)) # Item SD pooled - from Meade erratum (doi: 10.1037/a0020897)
      ESSDi <- (ES2 - ES1)/SDipool
      ES_diffs <- ES2 - ES1
      ABS_diffs <- abs(ES_diffs)
      round(cbind("SIDi"=ES_diffs,"UIDi"=ABS_diffs,"ESSDi"=ESSDi),digits)
    })
    
    # Test-level indices
    ## UETSDS: unsigned expected test score difference in the sample
    UETSDS <- mean(abs(apply(ETscores,1,diff)))
    ## Test D-Max: maximum difference in expected test score for any one person in sample
    TestDMax <- apply(ETscores,1,diff)[which.max(abs(apply(ETscores,1,diff)))]
    ##  STDS: signed expected test score difference in the sample (sum of SIDS)
    STDS <- sum(item.indices[,1])
    ## UTDS: unsigned expected test score difference in the sample (sum of UIDS)
    UTDS <- sum(item.indices[,2])
    ## ETSSD: expected test score standardized difference (Cohen's d metric)
    SDtpool <- sqrt(((Nf-1)*var(ETscores[,2]) + (Nf-1)*var(ETscores[,1]))/(2*Nf-2)) # Pooled test score SD 
    ETSSD <- unname(diff(colMeans(ETscores))/SDtpool)
    
    test.indices <- round(c("STDS"=STDS,"UTDS"=UTDS,"UETSDS"=UETSDS,"TestDmax"=TestDMax,"ETSSD"=ETSSD),digits)
    
    # Test-level indices for individuals
    ## UETSD: unsigned expected test score difference in an individual
    UETSDi <- abs(apply(ETscores,1,diff))
    ##  STDi: signed expected test score difference in the individual (sum of SIDS)
    STDi <- Reduce("+",item.indices.individual)[,1]
    ## UTDi: unsigned expected test score difference in the sample (sum of UIDS)
    UTDi <- Reduce("+",item.indices.individual)[,2]
    ## ETSSD: expected test score standardized difference (Cohen's d metric)
    ETSSDi <- unname(UETSDi/SDtpool)
    
    test.indices.individual <- round(cbind("STDi"=STDi,"UTDi"=UTDi,"UETSDi"=UETSDi,"ETSSDi"=ETSSDi),digits)
    
    return(list("items"=item.indices,"test"=test.indices,"item.individual"=item.indices.individual,"test.individual"=test.indices.individual))
  })
  names(ES_list) <- paste0(grps[ref.grp],":",grps[-ref.grp])
  return(ES_list)
}

### ggtraceplot - Traceline plot for items in 1-factor model
## Zack Williams
## 09/06/2020
ggtraceplot <- ggitemplot <- function(mod,which.items=1:extract.mirt(mod,"nitems"),theta_lim=c(-6,6),itemnames=NULL){
  nitems <- extract.mirt(mod,"nitems")
  if(extract.mirt(mod,"nfact") > 1){stop("Error: Function only works on unidimensional models. Use 'itemplot.bf' instead.")}
  which.item.names <- extract.mirt(mod,"itemnames")[which.items]
  # All of the points
  eval_pts <- seq(theta_lim[1],theta_lim[2],length.out = 1000)
  prob <- probtrace(mod,eval_pts)
  
  trace_list <- lapply(which.item.names,function(i){
    which.cols <- grep(paste0("^",i,"\\.P"),colnames(prob))
    as.data.frame(prob[,which.cols])
  })
  
  if(is.null(itemnames)){
    names(trace_list) <- which.item.names # use default item names
  } else{
    if(length(which.items)!=length(itemnames)){
      message("\nItemnames parameter is not the same length as number of items selected. Using default names.")
      names(trace_list) <- which.item.names # use default item names
    } else{
      names(trace_list) <- itemnames
    }
  }
  
  plot_list <- mclapply(1:length(trace_list),function(i){
    itrace <- trace_list[[i]]
    iname <- names(trace_list)[i]
    respOptions <- ncol(itrace)
    colnames(itrace) <- 1:respOptions
    itrace$Theta <- eval_pts
    tracedf <- tidyr::pivot_longer(itrace,1:respOptions,names_to = "Response",values_to = "Probability")
    ggplot(tracedf, aes(Theta,Probability)) + geom_line(aes(colour = Response)) + ggtitle(paste0("Trace Lines for ",iname)) +
      theme(legend.position="right",
            plot.title = element_text(size=14, face="bold", hjust=0.5),
            legend.title = element_text(size=10, face="bold", hjust=0.5),
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11),
            axis.title.x = element_text(size=12,face="bold"),
            axis.title.y = element_text(size=12,face="bold"),
            strip.text.y = element_text(size = 20, colour = "white"),
            strip.background = element_rect(colour="darkgrey", fill="darkgrey"))
  })
  return(plot_list)
}

#### itemplot.bf Marginal Traceline plot for items in bifactor model
### Zack Williams
## 05/08/2021
itemplot.bf.marginal <- function(mod,which.items=1:extract.mirt(mod,"nitems"),
                                 theta_lim=c(-6,6),quadpts=1000,cores=detectCores()-1){
  nitems <- extract.mirt(mod,"nitems")
  nfact <- extract.mirt(mod,"nfact")
  if(nfact < 2){stop("Error: Function only works on bifactor models. Use 'ggitemplot' instead.")}
  which.item.names <- extract.mirt(mod,"itemnames")[which.items] # Item names for plots
  item.factors <- t(apply(loads(mod,add.h2 = F),1,function(X){X != 0})) # Which dimensions to integrate over
  # Cut items that don't load on selected factor
  which.items <- which.items[(item.factors[,1])[which.items]]
  mod_coefs <- coef(mod,simplify=T) # Coefficients of model (needed for integral)
  
  theta_vals <- seq(theta_lim[1],theta_lim[2],length.out = quadpts)
  
  ### Function to integrate over nuisance dimensions
  fn <- function(theta_mat,theta_val,which.resp,item=item,mod=mod,which.factors=which.factors){
    th <- matrix(0,nfact,ncol(theta_mat))
    th[1,] <- theta_val
    th[which.factors,] <- theta_mat
    
    prob <- probtrace(item, t(th))
    sapply(1:nrow(prob),function(r){ # MVN density times item endorsement probability
      return(prob[r,] * mvnfast::dmvn(th[-1,r],mu=mod_coefs$mean[-1],sigma = mod_coefs$cov[-1,-1]))
    })[which.resp,,drop=FALSE]
  }
  if(cores > 1){
    suppressMessages(mirtCluster(remove=T)) # mirtClusters may screw this up, so remove any existing ones
  }
  # Now create list of probability traces for each item
  trace_list <- lapply(which.items,function(i){
    message(paste0("Calculating marginal trace lines for item #",i,"..."))
    item <- extract.item(mod,i) # extract specific item
    which.factors <- which(item.factors[i,-1]) + 1 # All specific factors with nonzero loadings
    n.resp <- ncol(probtrace(item,t(rep(0,nfact)))) # How many response categories?
    
    if(length(which.factors)==0){ # Only loads on general factor
      item_trace <- probtrace(item,cbind(theta_vals,rep.col(rep(0,quadpts),nfact-1)))
    } else{
      ### Generate marginal tracelines (the part that takes some time)
      item_trace <- t(pbsapply(theta_vals,function(theta){ # Parallelize over general factor theta values
        sapply(1:n.resp,function(resp){
          hcubature(fn, lower = matrix(rep(theta_lim[1],length(which.factors))), upper = matrix(rep(theta_lim[2],length(which.factors))),
                    mod=mod, absError = 1e-4, vectorInterface = T, theta_val = theta, item = item, which.factors = which.factors,
                    which.resp = resp)$integral
        })
      },cl=cores))
    }
    if(length(which.factors)>0){ # No need to smooth if no integration
      ### Item trace line smoothing
      item_trace <- apply(item_trace,2,function(X){
        X <- smooth.spline(x=theta_vals,y=X,spar=0.75)$y
        return(X)
      })
    }
    return(item_trace)
  })
  
  plot_list <- mclapply(1:length(trace_list),function(i){
    itrace <- data.frame(trace_list[[i]])
    iname <- which.item.names[i]
    respOptions <- ncol(itrace)
    # Rename columns of trace to be 1:k or 0:k-1
    if(any(mod@Data$data[,i]==0)){ # Item starts with 0 or 1?
      colnames(itrace) <- 0:(respOptions-1)
    } else{
      colnames(itrace) <- 1:respOptions
    }
    itrace$Theta <- theta_vals
    # Turn to 3-column df (Theta, Response, Probability)
    tracedf <- tidyr::pivot_longer(itrace,1:respOptions,names_to = "Response",values_to = "Probability")
    ggplot(tracedf, aes(Theta,Probability)) + geom_line(aes(colour = Response)) + ggtitle(paste0("Marginal Trace Lines for ",iname)) +
      ylab("Marginal Response Probability") + xlab(expression("Latent Trait Score ("*theta[G]*")")) +
      theme(legend.position="right",
            plot.title = element_text(size=14, face="bold", hjust=0.5),
            legend.title = element_text(size=10, face="bold", hjust=0.5),
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11),
            axis.title.x = element_text(size=12,face="bold"),
            axis.title.y = element_text(size=12,face="bold"),
            strip.text.y = element_text(size = 20, colour = "white"),
            strip.background = element_rect(colour="darkgrey", fill="darkgrey"))
  })
  return(plot_list)
}

# Conditional trace lines for bifactor models (theta exactly at zero)
# Zack Williams
# 11/29/19
itemplot.bf.conditional <- function(mod,which.items=1:extract.mirt(mod,"nitems"),theta_lim=c(-6,6),
                                    theta_step=0.1,cond.theta=rep(0,extract.mirt(mod, "nfact")-1),zero.start=TRUE){
  nitems <- extract.mirt(mod,"nitems")
  nfact <- extract.mirt(mod,"nfact")
  if(nfact < 2){stop("Error: Function only works on bifactor models. Use 'itemplot' instead.")}
  which.item.names <- extract.mirt(mod,"itemnames")[which.items]
  quadpts <- length(seq(theta_lim[1],theta_lim[2],theta_step))
  # All of the points to be evaluated
  theta_mat <- cbind(seq(theta_lim[1],theta_lim[2],length.out = quadpts),
                     rep.row(cond.theta,quadpts))
  
  prob <- probtrace(mod,theta_mat)

  # Probability traces for selected items at all theta points on grid
  trace_list <- lapply(which.item.names,function(i){
    which.cols <- grep(paste0("^",i,"\\.P"),colnames(prob))
    as.data.frame(prob[,which.cols])
  })
  names(trace_list) <- which.item.names
  
  plot_list <- mclapply(trace_list,function(itrace){
    iname <- gsub("(^.*)\\.P\\.\\d+$","\\1",colnames(itrace)[1])
    respOptions <- ncol(itrace)
    if(zero.start){
      colnames(itrace) <- 0:(respOptions-1)
    } else{
      colnames(itrace) <- 1:respOptions
    }
    itrace$Theta <- theta_mat[,1]
    tracedf <- tidyr::pivot_longer(itrace,1:respOptions,names_to = "Response",values_to = "Probability")
    ggplot(tracedf, aes(Theta,Probability)) + geom_line(aes(colour = Response)) + ggtitle(paste0("Conditional Trace Lines for ",iname)) +
      ylab("Conditional Response Probability") + xlab(expression("Latent Trait Score ("*theta[G]*")")) +
      theme(legend.position="right",
            plot.title = element_text(size=14, face="bold", hjust=0.5),
            legend.title = element_text(size=10, face="bold", hjust=0.5),
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11),
            axis.title.x = element_text(size=12,face="bold"),
            axis.title.y = element_text(size=12,face="bold"),
            strip.text.y = element_text(size = 20, colour = "white"),
            strip.background = element_rect(colour="darkgrey", fill="darkgrey"))
  })
  return(plot_list)
}

#### wABC - Weighted area between curves (or volumes between surfaces in case of multidimensional model)
#### Zack Williams
#### 9/22/2019
# This is an effect size for DIF. For unidimensional models, values < 0.06*number of response options
# are suggestive of negligible DIF.
# Edelen et al. (2015) http://doi.org/10.1007/s11136-013-0540-4
# Note that the multidimensional, marginal, and conditional versions of this index were made up by me. No references published. 
# I would caution against the use of non-unidimensional wABC metrics for more than descriptive reasons.
DIF.WABC <- dif.WABC <- dif.wabc <- dif.wABC <- DIF.wABC <- 
  function(dif_model,theta_lim=c(-4,4),which.items=1:extract.mirt(dif_model, "nitems"),
                 quadpts=NULL,theta_step=NULL,conditional=F,marginal=F,
                 cond.theta=rep(0,extract.mirt(dif_model, "nfact")-1),digits=3){
  options("mc.cores"=detectCores()-1) # get that multicore on
  # Group models
  if(extract.mirt(dif_model,"ngroups")!=2){stop("ERROR: wABC only defined for 2-group models")}
  g1 <- extract.group(dif_model,1)
  g2 <- extract.group(dif_model,2)
  nitems <- extract.mirt(dif_model,"nitems")
  nfact <- extract.mirt(dif_model,"nfact")
  which.item.names <- extract.mirt(dif_model,"itemnames")[which.items]
  if((conditional | marginal) & nfact==1){
    conditional <- marginal <- F # Can't calculate conditional wABC for one factor model
  }
  if(is.null(quadpts) & is.null(theta_step)){ # Neither quadpts nor theta_step defined
    if(conditional){
      theta_step <- 0.1
      quadpts <- length(seq(theta_lim[1],theta_lim[2],theta_step))
    } else{
      quadpts <- switch(as.character(nfact), '1'=81, '2'=31, '3'=15, '4'=9, '5'=7, 3) # Default: same number of quadpts as model
      theta_step <- diff(seq(theta_lim[1],theta_lim[2],length.out = quadpts))[1]
    }
  } else if(!is.null(quadpts)){ # quadpts only defined
    theta_step <- diff(seq(theta_lim[1],theta_lim[2],length.out = quadpts))[1]
  } else if(!is.null(theta_step)){ # theta_step defined
    quadpts <- length(seq(theta_lim[1],theta_lim[2],theta_step))
  } else{ # neither defined
    quadpts <- switch(as.character(nfact), '1'=81, '2'=31, '3'=15, '4'=9, '5'=7, 3)
    theta_step <- diff(seq(theta_lim[1],theta_lim[2],length.out = quadpts))[1]
  }
  
  if(marginal){
    # All of the points
    eval_pts <- gtools::permutations(quadpts,nfact-1,seq(theta_lim[1],theta_lim[2],length.out = quadpts))
    
    # All of the conditional point values
    theta_mat <- unlist(apply(eval_pts,1,function(pt){
      list(cbind(seq(theta_lim[1],theta_lim[2],length.out = quadpts),
                 matrix(pt,nrow = quadpts,ncol = nfact-1)))
    }),recursive=F)
  } else if(!conditional){
    theta_mat <- gtools::permutations(quadpts,nfact,seq(theta_lim[1],theta_lim[2],length.out = quadpts))
  } else{
    theta_mat <- cbind(seq(theta_lim[1],theta_lim[2],length.out = quadpts),
                       rep.row(cond.theta,quadpts))
  }
  
  # If marginal, do all this many times
  if("list" %in% class(theta_mat)){
    prob.1 <- mclapply(theta_mat,function(X) probtrace(g1,X))
    prob.2 <- mclapply(theta_mat,function(X) probtrace(g2,X))
    item_scores <- as.numeric(gsub(".*P\\.(\\d)$","\\1",colnames(prob.1[[1]])))
    for(i in 1:length(prob.1)){
      for(j in 1:ncol(prob.1[[1]])){
        prob.1[[i]][,j] <- prob.1[[i]][,j]*item_scores[j]
        prob.2[[i]][,j] <- prob.2[[i]][,j]*item_scores[j]
      }
    }
    # Estimated scores for each group on selected items at all theta points on grid
    Escores <- lapply(which.item.names,function(i){
      which.cols <- grep(paste0("^",i,"\\.P"),colnames(prob.1[[1]]))
      mclapply(1:length(prob.1),function(j){
        cbind("E1"=rowSums(prob.1[[j]][,which.cols]),"E2"=rowSums(prob.2[[j]][,which.cols]))
      })
    })
    names(Escores) <- which.item.names
  } else{ # Not marginal
    prob.1 <- probtrace(g1,theta_mat)
    prob.2 <- probtrace(g2,theta_mat)
    item_scores <- as.numeric(gsub(".*P\\.(\\d)$","\\1",colnames(prob.1)))
    for(i in 1:ncol(prob.1)){
      prob.1[,i] <- prob.1[,i]*item_scores[i]
      prob.2[,i] <- prob.2[,i]*item_scores[i]
    }
    # Estimated scores for each group on selected items at all theta points on grid
    Escores <- mclapply(which.item.names,function(i){
      which.cols <- grep(paste0("^",i,"\\.P"),colnames(prob.1))
      cbind("E1"=rowSums(prob.1[,which.cols]),"E2"=rowSums(prob.2[,which.cols]))
    })
    names(Escores) <- which.item.names
  }
  
  # Generate values for group 1 weights
  mu_g1 <- coef(g1)$GroupPars[1:nfact]
  sigma_g1 <- diag(nfact)
  sigma_g1[lower.tri(sigma_g1,diag=T)] <- coef(g1)$GroupPars[-(1:nfact)]
  sigma_g1 <- t(sigma_g1)
  sigma_g1[lower.tri(sigma_g1,diag=T)] <- coef(g1)$GroupPars[-(1:nfact)]
  # Generate values for group 2 weights
  mu_g2 <- coef(g2)$GroupPars[1:nfact]
  sigma_g2 <- diag(nfact)
  sigma_g2[lower.tri(sigma_g2,diag=T)] <- coef(g2)$GroupPars[-(1:nfact)]
  sigma_g2 <- t(sigma_g2)
  sigma_g2[lower.tri(sigma_g2,diag=T)] <- coef(g2)$GroupPars[-(1:nfact)]
  
  # Weights derived from multivariate normal density
  if("list" %in% class(theta_mat)){ # marginal
    # Group 1 weights
    g1_weights <- mclapply(theta_mat,function(MAT){
      apply(MAT,1,function(x){
        dmvnorm(x,mu_g1,sigma_g1)
      })
    })
    # Group 2 weights
    g2_weights <- mclapply(theta_mat,function(MAT){
      apply(MAT,1,function(x){
        dmvnorm(x,mu_g2,sigma_g2)
      })
    })
  } else{ # not marginal
    # Group 1 weights
    g1_weights <- apply(theta_mat,1,function(x){
      dmvnorm(x,mu_g1,sigma_g1)
    })
    # Group 2 weights
    g2_weights <- apply(theta_mat,1,function(x){
      dmvnorm(x,mu_g2,sigma_g2)
    })
  }
  
  # wABC is a weighted mean of the areas between the curves using weights from both groups.
  if(marginal){
    # Marginal wABC: weighted mean of all conditionals across all lists
    wABC_list <- mclapply(Escores,function(item){
      wABCs_Weights <- sapply(seq_len(length(item)),function(qpt){
        wABC_Mean <- wABC_1 <- wABC_2 <- rep(NA,quadpts)
        # Get (multivariate) normal density weights for each group
        eval_point <- theta_mat[[qpt]][1,-1] # group factor values for this conditional wABC
        if(nfact==2){ # 2-factor case
          wt_1 <- dnorm(eval_point,mu_g1[-1],sigma_g1[-1,-1])
          wt_2 <- dnorm(eval_point,mu_g2[-1],sigma_g2[-1,-1])
        } else{ # more than 2 factors - mvn density
          wt_1 <- dmvnorm(eval_point,mu_g1[-1],sigma_g1[-1,-1])
          wt_2 <- dmvnorm(eval_point,mu_g2[-1],sigma_g2[-1,-1])
        }
        # These weights are the same if group factor means/covs are the same
        # But the weighted weight handles cases where they aren't (do these matter?)
        weighted_weight <- weighted.mean(c(wt_1,wt_2),c(g1@Data$N,g2@Data$N))
        
        for(i in 1:quadpts){
          wABC_1[i] <- sum(abs(apply(item[[qpt]],1,diff))*g1_weights[[i]]*theta_step)
          wABC_2[i] <- sum(abs(apply(item[[qpt]],1,diff))*g2_weights[[i]]*theta_step)
          wABC_Mean[i] <- weighted.mean(c(wABC_1[i],wABC_2[i]),c(g1@Data$N,g2@Data$N))
        }
        return(c(mean(wABC_Mean),weighted_weight)) # Equivalent to a single conditional wABC * the weight (p(GrpFactors))
      })
      # Return weighted mean of wABCs at all quadrature points (weighted by mean of group factor densities)
      return(weighted.mean(wABCs_Weights[1,],wABCs_Weights[2,]))
    })
    return(cbind("wABC.margin"=signif(unlist(lapply(wABC_list,mean)),digits)))
  } else if(!conditional){
    # wABC_1 = sum(abs(E1-E2)*g1_weights*theta_step^nfact)
    # wABC_2 = sum(abs(E1-E2)*g2_weights*theta_step^nfact)
    # wABC = weighted mean of wABC_1 and wABC_2
    wABC <- unlist(lapply(Escores,function(X){
      wABC_1 <- sum(abs(apply(X,1,diff))*g1_weights*theta_step^nfact)
      wABC_2 <- sum(abs(apply(X,1,diff))*g2_weights*theta_step^nfact)
      return(weighted.mean(c(wABC_1,wABC_2),c(g1@Data$N,g2@Data$N)))
    }))
    return(cbind("wABC"=signif(wABC,digits)))
  } else{ # If conditional, just looking at area, not volume
    # Same as wABC, but without theta_step being raised to a power
    wABC <- unlist(lapply(Escores,function(X){
      wABC_1 <- sum(abs(apply(X,1,diff))*g1_weights*theta_step)
      wABC_2 <- sum(abs(apply(X,1,diff))*g2_weights*theta_step)
      return(weighted.mean(c(wABC_1,wABC_2),c(g1@Data$N,g2@Data$N)))
    }))
    return(cbind("wABC.cond"=signif(wABC,digits)))
  }
}

# Wald Test
# Cannibalized from mirt 'DIF' code
WALD.DIF <- function(MGmodel, which.par=NULL, items2test = 1:extract.mirt(MGmodel, 'nitems'),
                 p.adj = 'fdr'){
  if(!extract.mirt(MGmodel,"secondordertest")){
    stop("Second-order test: indicates information matrix is too inaccurate. DIF tests cannot be computed.")
  }
  coefs <- coef(MGmodel,simplify=T)[[1]]$items
  coefs[apply(coefs,2,function(X){X==0 | duplicated(X)})] <- NA # if constrained to 0 or set equal to another loading, ignore par
  all.pars <- colnames(coefs)
  if(is.null(which.par)){
    which.par <- all.pars
  } else if(any(which.par %in% c("slopes","slope","Slopes","Slope"))){
    which.par <- all.pars[grep("^a",all.pars)]
  } else if(any(which.par %in% c("intercepts","intercept","Intercepts","Intercept"))){
    which.par <- all.pars[grep("^[cd]",all.pars)]
  }
  # Extract item names, values
  itemnames <- colnames(MGmodel@Data$data)
  values <- mod2values(MGmodel)
  parnum <- list()
  
  # Code cannibalized from DIF loop_test sub-function, put into a sapply loop
  Waldmat <- t(sapply(items2test,function(item){
    # Check to see which parameters item has
    which.par.i <- all.pars[!is.na(as.data.frame(coefs)[item,])]
    which.par.i <- which.par.i[which(which.par.i %in% which.par)]
          
    if(!length(which.par.i)){
      stop('Item ', item, ' does not contain any of the parameters defined in which.par.
                 Consider removing it from the item2test input or adding relevant parameters
                 to which.par', call.=FALSE)
    }
    # Which parameter numbers
    for(i in seq_len(length(which.par.i))){
      parnum[[i]] <- values$parnum[values$name == which.par.i[i] &
                                       values$item == itemnames[item]]
    }
    wv <- mirt::wald(MGmodel)
    infoname <- names(wv)
    # Collect only values of infoname that contain which.par.i
    # infoname <- infoname[sort(unlist(sapply(which.par.i,function(par){grep(paste0("^",par,"\\.\\d+$"),infoname)})))]
    L <- matrix(0, length(parnum), length(infoname))
    for(i in seq_len(length(which.par.i))){
      L[i, paste0(which.par.i[i], '.', parnum[[i]][1L]) == infoname] <- 1
      L[i, paste0(which.par.i[i], '.', parnum[[i]][2L]) == infoname] <- -1
    }
    res <- unlist(suppressWarnings(mirt::wald(MGmodel, L)))
    return(res)
  }))
  rownames(Waldmat) <- itemnames[items2test]
  apply(Waldmat,2,round,3)
  
  Waldmat <- as.data.frame(round(cbind(Waldmat,"p.adj"=p.adjust(Waldmat[,3],p.adj)),3))
  return(Waldmat)
}


### DIF.wald.AATA function
### Zack Williams
### Updated 05/07/21 to allow for two-step implementation
# Based on methods of Wang & Woods (2017)
# http://doi.org/10.1177/0146621616668014
# Note: Cited paper uses SEM standard errors, which I find do not converge much of the time
# The default is thus to use SE.type="Oakes" for more stable SE estimates. See ?mirt::mirt for more info
###########################################################
# Performs the following steps:
# Uses Wald Anchor-All-Test-All (AATA) algorithm in first step to find anchor items
# In second step, uses anchor items (determined using largest slope parameter items that didn't show DIF) to
# calculate a Wald test for DIF using all of the remaining items.
#   - Up to 25% items selected as anchors
#   - Default result is to FDR adjust the p-values
DIF.wald.AATA <- function(data,model,model2 = paste0("G = 1-", ncol(data)),group,test_pars=NULL,
                          bifactor=F,include.ES=T,method=NULL,p.adj="BH",seed=12345,technical=NULL,TOL=NULL,method.2step="QMCEM",
                          NCYCLES=5000,MHRM_SE_draws=5000,SEtol = 1e-4,p.cutoff=0.05,verbose=T,fix_MV=NULL,
                          SE.type="Oakes",ref.group=NULL,cores=detectCores()-1,...){
  suppressMessages(mirtCluster(spec=cores)) 
  # If any NA values in group vector, remove those rows from data
  if(any(is.na(group))){
    message(paste0("Removing ",sum(is.na(group))," rows due to missing values of grouping variable.\n"))
    data <- data[which(!is.na(group)),]
    group <- group[which(!is.na(group))]
  }
  # If any rows in data are all NA, remove those too
  empty_rows <- apply(data,1,function(X){all(is.na(X))})
  if(any(empty_rows)){
    message(paste0("Removing ",sum(empty_rows)," rows due to completely missing item data.\n"))
    data <- data[-which(empty_rows),]
    group <- group[-which(empty_rows)]
  }
  # Fix reference group to be first level
  group <- as.factor(group)
  grp_levels <- levels(group)
  ngrps <- length(grp_levels)
  if(is.null(ref.group)){
    ref.group <- 1
  } else if(class(ref.group)=="character"){
    ref.group <- which(levels(group)==ref.group) # Change to numeric
    if(length(ref.group)==0){stop("ERROR: Invalid reference group value")}
  }
  grp_levels <- levels(group) <- c(levels(group)[ref.group],levels(group)[-ref.group])
  
  # Fit a dummy model to extract information
  if(bifactor==T){
    mod_dummy <- suppressMessages(bfactor(data=data,model=model,model2=model2,group=group,SE=F,TOL=NA,
                                          verbose=F,technical=c(technical,list(NCYCLES=1))))
  } else{
    mod_dummy <- suppressMessages(multipleGroup(data=data,model=model,group=group,SE=F,TOL=NA,
                                                verbose=F,technical=c(technical,list(NCYCLES=1))))
    # Also carries over additional model information if input (e.g., constraints)
    if(length(model)==1 & (is.character(model)) | "mirt.model" %in% class(model)){
      model2 <- model
    } else{
      model2 <- NULL
    }
  }
  # Extract info from dummy model
  loads <- extract.mirt(mod_dummy,"F")
  nfac <- ncol(loads)
  if(!is.null(fix_MV)){
    which_factors <- (1:nfac)[-fix_MV]
  } else{
    which_factors <- 1:nfac
  }
  nitems <- nrow(loads)
  fac_names <- colnames(loads)
  item_names <- rownames(loads)
  fac_items <- apply(loads,2,function(X){which(abs(X)>0)})
  # Turn into list for 1-factor case or multi-factor case where same # of items per factor
  # Thanks to Thomas McCauley for bringing the latter problem to my attention - fixed as of 04/10/20
  if(!"list" %in% class(fac_items)){
    if(nrow(fac_items)==ncol(fac_items)){ # same number of factors per item
      fac_items <- lapply(fac_names,function(X){
        which(loads[,X] > 0)
      })
    } else{ # Unidimensional case
      fac_items <- list(fac_items[,1])
    }
  }
  
  # Create mirt model with all parameters fixed between groups, except group mean/var free to vary
  # Unfortunately, free_means and free_vars only do this for factor 1, so custom model needs to be used (updated 03/06/21)
  ### Update 05/11/21 - Pasting together mirt.model input and text does not preserve constraints. Now converts mirt.model back to string
  if(class(model2)[1]=="mirt.model"){
    base_model <- model_freeMV <- paste(apply(model2$x,1,paste,collapse=" = "),collapse ="\n ")
  } else{
    base_model <- model_freeMV <- paste(c(sapply(1:nfac,function(i){
      paste0(fac_names[i],"=",paste(fac_items[[i]],collapse=","))
    }),model2),collapse="\n")
  }
  
  for(i in 2:ngrps){
    free_pars <- paste(sapply(which_factors,function(f){
      paste0("(GROUP, MEAN_",f,"), (GROUP, COV_",f,f,")")
    }),collapse=", ")
    model_iter <- model_freeMV <- paste0(model_freeMV,"\n",
                                         "FREE [",grp_levels[i],"] = ",free_pars,"\n")
    # If fix_MV is true, then M/V are still freed for WALD-2 analysis; they're just fixed in future analyses at the WALD-2 values.
    if(!is.null(fix_MV)){
      fixed_pars <- paste(sapply(fix_MV,function(f){
        paste0("(GROUP, MEAN_",f,"), (GROUP, COV_",f,f,")")
      }),collapse=", ")
      model_iter <- model_freeMV # model_iter is used for the iterations—keeps fixed_pars fixed 
      model_freeMV <- paste0(model_freeMV,"\n",
                             "FREE [",grp_levels[i],"] = ",fixed_pars,"\n")
    }
  }
  
  # Start by fitting model with all parameters fixed between groups, with group mean/var free to vary
  if(bifactor==T){
    if(!is.null(method)){
      if(verbose==T){message("Only EM method can be used for bifactor models. Ignoring 'method' argument.")}
    }
    mod_fixed <- suppressMessages(bfactor(data=data,model=model,model2=model_freeMV,group=group,
                                          invariance = c("slopes","intercepts"),SE=F,SE.type=SE.type,TOL=TOL,
                                          verbose=verbose,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
  } else{
    # Default method: EM (3 or fewer factors) or 2step (>3 factors)
    if(is.null(method)){method <- ifelse(nfac > 3,"2step","EM")}
    # if stochastic method, switch to MHRM standard error
    if(method %in% c("SEM","MHRM","MCEM") & SE.type=="Oakes"){
      SE.type <- "FMHRM"
    } else if(method %in% c("2step","twostep","2STEP","TWOSTEP")){
      # If 2step is method, use MHRM
      mod_fixed <- suppressMessages(multipleGroup(data=data,model=model_freeMV,group=group,method="MHRM",
                                                  invariance = c("slopes","intercepts"),SE=F,TOL=NULL,
                                                  verbose=verbose,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
    } else{
      mod_fixed <- suppressMessages(multipleGroup(data=data,model=model_freeMV,group=group,method=method,
                                                  invariance = c("slopes","intercepts"),SE=F,TOL=TOL,
                                                  verbose=verbose,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
    }
}
  
  # Extract relevant information from model
  fixed_coefs <- coef(mod_fixed,simplify=T)
  mod_means <- lapply(fixed_coefs,function(GRP){GRP$means})
  mod_cov <- lapply(fixed_coefs,function(GRP){GRP$cov})
  
  # Now make model where M/V are fixed at earlier values
  # Start with base model + fixed means/variances for identification
  if(!is.null(fix_MV)){ # If there are some fixed model pars, fix both them and the free pars
    model_MV <- paste0(base_model,"\nFIXED = ",free_pars,", ",fixed_pars) # if fixed pars exist, have them fixed here as well
  } else {
    model_MV <- paste0(base_model,"\nFIXED = ",free_pars) 
  }
  
  
  for(i in 2:ngrps){
    for(f in 1:nfac){
      model_MV <- paste0(model_MV,"\n",
                         "START [",grp_levels[i],"] = (GROUP, COV_",f,f,",",as.matrix(unlist(mod_cov[[i]]))[f,f],")\n",
                         "START [",grp_levels[i],"] = (GROUP, MEAN_",f,",",(mod_means[[i]])[f],")\n")
      # Also provide these start values for model_iter
      model_iter <- paste0(model_iter,"\n",
                           "START [",grp_levels[i],"] = (GROUP, COV_",f,f,",",as.matrix(unlist(mod_cov[[i]]))[f,f],")\n",
                           "START [",grp_levels[i],"] = (GROUP, MEAN_",f,",",(mod_means[[i]])[f],")\n")
    }
  }

  # Fix mean and variance of focal group, then re-fit model allowing all parameters to vary (Wald-2)
  if(bifactor==T){
    mod_fixedMV <- suppressMessages(bfactor(data=data,model=model,model2=model_iter,group=group,verbose=verbose,SE.type=SE.type,
                                            SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
    
  } else{
    if(method %in% c("2step","twostep","2STEP","TWOSTEP")){
      # If 2step is method, use multipleGroup.2step function
      mod_fixedMV <- suppressMessages(multipleGroup.2step(data=data,model=model_MV,group=group,SE.type=SE.type,verbose=verbose,method=method.2step,TOL=TOL,
                                                    SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,
                                                                                    SEtol = SEtol)),...))
    } else{
      mod_fixedMV <- suppressMessages(multipleGroup(data=data,model=model_MV,group=group,SE.type=SE.type,verbose=verbose,method = method,TOL=TOL,
                                                    SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,MHRM_SE_draws=MHRM_SE_draws,
                                                                                    SEtol = SEtol)),...)) 
    }
  }
  
  fixedMV_pars <- coef(mod_fixedMV,simplify=T)[[2]]$items

  # Calculate Wald-2 DIF and store items that exhibit no DIF as potential anchors
  dif_wald_2 <- suppressMessages(WALD.DIF(mod_fixedMV,test_pars,p.adj = p.adj))
  if(p.adj != "none"){
    non_dif_items <- which(dif_wald_2[,4] > p.cutoff) # 
  } else{
    non_dif_items <- which(dif_wald_2[,3] > p.cutoff) # non-adjusted p.value
  }
  # Add effect sizes to Wald-2 output (if there are 2 groups; otherwise, find in $DIF.EffectSize call)
  if(include.ES==T){
    ES.wald2 <- DIF.ES(mod_fixedMV)
    if(ngrps==2){
      if(nfac==1){ # 1 factor: use wABC, SIDS, UIDS, ESSD
        wABC <- DIF.wABC(mod_fixedMV)
        itemES <- ES.wald2[[1]]$items
        dif_wald_2 <- cbind(dif_wald_2,wABC,itemES[,c(1,2,4)])
      } else{ # more than 1 factor
        itemES <- ES.wald2[[1]]$items
        dif_wald_2 <- cbind(dif_wald_2,itemES[,c(1,2,4)])
      }
    }
  } else{
    ES.wald2 <- NULL
  }
  # If all items show DIF, stop here.
  if(length(non_dif_items)==0){ # all items show DIF on Wald-2
    message("No anchor items! Look over data carefully.\n")
    # If 2 groups, add conditional wABC to DIF output
    result <- list("Iterations"=0,"DIF.model"=mod_fixedMV,"coef.DIF.model"=coef(mod_fixedMV,simplify=T),
                   "Results"=dif_wald_2,"DIF.items"=item_names,"DIF.EffectSize"= ES.wald2,
                   "p.adjust.method"=p.adj,"p.cutoff"=p.cutoff,bifactor=bifactor,"which.pars"=test_pars,
                   "Call"=match.call(),"Type"="wald.2step","DIF.model.wald2"=mod_fixedMV,"ES.anchors"=FALSE)
    class(result) <- "DIF.wald"
    
    if(verbose==T){ # if verbose, print final model (in this case, the Wald-2 model).
      cat("\n")
      print(result$DIF.model)
      cat("\n")
    }
    return(result)
  }
  # Assuming there are eligible anchor items, select items as anchors:
  # No more than 25% of items should be selected as anchors per Wang & Woods paper
  n_anchors <- ifelse(length(non_dif_items)>floor(nitems*0.25),
                      max(floor(nitems*0.25),1), # minimum of 1 anchor selected
                      length(non_dif_items)) # use all non-DIF items if less than 25% of items

  unadjusted_nonDIF <- which(dif_wald_2[,3] > p.cutoff)
  # Now pick the n_anchors with the highest discrimination parameter a1, provided they didn't show DIF
  which_anchors <- (non_dif_items[order(fixedMV_pars[non_dif_items,1],decreasing = T)])[1:n_anchors]
  
  # Now refit model with selected anchor items and test DIF
  if(bifactor==T){
    mod_anchor <- suppressMessages(bfactor(data=data,model=model,model2=model_iter,group=group,verbose=verbose,SE.type=SE.type,
                                           invariance = c(item_names[which_anchors]),TOL=TOL,
                                           SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
  } else{
    if(method %in% c("2step","twostep","2STEP","TWOSTEP")){
      # If 2step is method, use multipleGroup.2step function
      mod_anchor <- suppressMessages(multipleGroup.2step(data=data,model=model_iter,group=group,verbose=verbose,SE.type=SE.type,
                                                   invariance = c(item_names[which_anchors]),TOL=TOL,method=method.2step,
                                                   SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,MHRM_SE_draws=MHRM_SE_draws,
                                                                                   SEtol = SEtol)),...))
    } else{
      mod_anchor <- suppressMessages(multipleGroup(data=data,model=model_iter,group=group,verbose=verbose,SE.type=SE.type,method = method,
                                                   invariance = c(item_names[which_anchors]),TOL=TOL,
                                                   SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,MHRM_SE_draws=MHRM_SE_draws,
                                                                                   SEtol = SEtol)),...)) 
    }
  }
  
  dif <- WALD.DIF(mod_anchor,test_pars,items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
  
  # Post-hoc tests for just slopes and just intercepts
  if(is.null(test_pars)){
    dif_slope <- WALD.DIF(mod_anchor,which.par = "slopes",
                          items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
    dif_intercept <- WALD.DIF(mod_anchor,which.par = "intercepts",
                              items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
  } else{
    dif_slope <- WALD.DIF(mod_anchor,which.par = test_pars[grep("^a",test_pars)],
                          items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
    dif_intercept <- WALD.DIF(mod_anchor,which.par = test_pars[grep("^[dc]",test_pars)],
                              items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
  }
  # Find which items were tested
  tested_items <- item_names[-which_anchors]
  # Evaluate which tested items showed DIF or no DIF
  if(p.adj != "none"){
    non_dif_items <- which(item_names %in% tested_items[which(dif[,4] > p.cutoff)])
    dif_items <- item_names[which(item_names %in% tested_items[which(dif[,4] < p.cutoff)])] 
  } else{
    non_dif_items <- which(item_names %in% tested_items[which(dif[,3] > p.cutoff)])
    dif_items <- item_names[which(item_names %in% tested_items[which(dif[,3] < p.cutoff)])]
  }
  if(include.ES){
    # Add Effect sizes
    ES <- DIF.ES(mod_anchor,which.items = which(item_names %in% tested_items))
    if(ngrps==2){
      if(nfac==1){ # 1 factor: use wABC, SIDS, UIDS, ESSD
        wABC <- DIF.wABC(mod_anchor,which.items = which(item_names %in% tested_items))
        itemES <- ES[[1]]$items
        dif <- cbind(dif,wABC,itemES[,c(1,2,4)])
      } else{ # more than 1 factor
        itemES <- ES[[1]]$items
        dif <- cbind(dif,itemES[,c(1,2,4)])
      }
    }
  } else{
    ES <- NULL
  }
  result <- list("DIF.model"=mod_anchor,"coef.DIF.model"=coef(mod_anchor,simplify=T),
                 "Results"=dif,"DIF.items"=dif_items,"DIF.EffectSize"= ES,"p.adjust.method"=p.adj,
                 "p.cutoff"=p.cutoff,bifactor=bifactor,"Call"=match.call(),"Type"="wald.2step",
                 "which.pars"=test_pars,"DIF.model.wald2"=mod_fixedMV,
                 "Posthoc.Slope"=dif_slope,"Posthoc.Intercept"=dif_intercept,"ES.anchors"=FALSE)
  class(result) <- "DIF.wald"
  if(verbose==T){ # if verbose, print final model after fitting all models.
    cat("\n")
    print(result$DIF.model)
    cat("\n")
  }
  return(result)
}

### DIF.wald.iter function
### Zack Williams
### Updated 05/07/21 to allow for 2step integration (MHRM them QMC)
# Based on methods of Cao, Tay, & Liu (2017)
# https://doi.org/10.1177/0013164416637104
# Note: Cited paper uses SEM standard errors, which I find do not converge much of the time
# The default is thus to use SE.type="Oakes" for more stable SE estimates. See ?mirt::mirt for more info
###########################################################
# Performs the following steps:
# Uses Wald Anchor-All-Test-All (AATA) algorithm in first step to test all items for DIF
# In second step, uses all non-DIF (FDR p > 0.05 by default) items as anchors
#   - If no DIF items in step 1, returns that data frame, stating there is no DIF
# Iteratively test DIF using anchor items. 
#   - Any items not showing DIF used as anchor items in the next iteration
#   - Test until all non-anchor items show DIF or all items are selected as anchors (i.e., no DIF)
#   - Default result is to FDR adjust the p-values
DIF.wald.iter <- function(data,model,model2 = paste0("G = 1-", ncol(data)),group,
                          test_pars=NULL,bifactor=F,ref.group=NULL,include.ES=T,method=NULL,p.adj="BH",TOL=NULL,
                          seed=12345,NCYCLES=5000,MHRM_SE_draws=5000,SEtol = 1e-4,p.cutoff=0.05,verbose=T,fix_MV=NULL,
                          ES_anchors=FALSE,ES.cutoff=0.2,method.2step="QMCEM",technical=NULL,SE.type="Oakes",cores=detectCores()-1,...){
  suppressMessages(mirtCluster(spec=cores)) # Define mirtCluster (does nothing if already defined)
  # If any NA values in group vector, remove those rows from data
  if(any(is.na(group))){
    message(paste0("Removing ",sum(is.na(group))," rows due to missing values of grouping variable.\n"))
    data <- data[which(!is.na(group)),]
    group <- group[which(!is.na(group))]
  }
  # If any rows in data are all NA, remove those too
  empty_rows <- apply(data,1,function(X){all(is.na(X))})
  if(any(empty_rows)){
    message(paste0("Removing ",sum(empty_rows)," rows due to completely missing item data.\n"))
    data <- data[-which(empty_rows),]
    group <- group[-which(empty_rows)]
  }
  # Fix reference group to be first level
  group <- as.factor(group)
  grp_levels <- levels(group)
  ngrps <- length(grp_levels)
  if(is.null(ref.group)){
    ref.group <- 1
  } else if(class(ref.group)=="character"){
    ref.group <- which(levels(group)==ref.group) # Change to numeric
    if(length(ref.group)==0){stop("ERROR: Invalid reference group value")}
  }
  grp_levels <- levels(group) <- c(levels(group)[ref.group],levels(group)[-ref.group])
  
  # Fit a dummy model to extract information
  if(bifactor==T){
    mod_dummy <- suppressMessages(bfactor(data=data,model=model,model2=model2,group=group,SE=F,TOL=NA,
                                          verbose=F,technical=c(technical,list(NCYCLES=1))))
  } else{
    mod_dummy <- suppressMessages(multipleGroup(data=data,model=model,group=group,SE=F,TOL=NA,
                                                verbose=F,technical=c(technical,list(NCYCLES=1))))
    # Also carries over additional model information if input (e.g., constraints)
    if(length(model)==1 & (is.character(model)) | "mirt.model" %in% class(model)){
      model2 <- model
    } else{
      model2 <- NULL
    }
  }
  # Extract info from dummy model
  loads <- extract.mirt(mod_dummy,"F")
  nfac <- ncol(loads)
  if(!is.null(fix_MV)){
    which_factors <- (1:nfac)[-fix_MV]
  } else{
    which_factors <- 1:nfac
  }
  nitems <- nrow(loads)
  fac_names <- colnames(loads)
  item_names <- rownames(loads)
  fac_items <- apply(loads,2,function(X){which(abs(X)>0)})
  # Turn into list for 1-factor case or multi-factor case where same # of items per factor
  # Thanks to Thomas McCauley for bringing the latter problem to my attention - fixed as of 04/10/20
  if(!"list" %in% class(fac_items)){
    if(nrow(fac_items)==ncol(fac_items)){ # same number of factors per item
      fac_items <- lapply(fac_names,function(X){
        which(loads[,X] > 0)
      })
    } else{ # Unidimensional case
      fac_items <- list(fac_items[,1])
    }
  }
  
  # Create mirt model with all parameters fixed between groups, except group mean/var free to vary
  # Unfortunately, free_means and free_vars only do this for factor 1, so custom model needs to be used (updated 03/06/21)
  ### Update 05/11/21 - Pasting together mirt.model input and text does not preserve constraints. Now converts mirt.model back to string
  if(class(model2)[1]=="mirt.model"){
    base_model <- model_freeMV <- paste(apply(model2$x,1,paste,collapse=" = "),collapse ="\n ")
  } else{
    base_model <- model_freeMV <- paste(c(sapply(1:nfac,function(i){
      paste0(fac_names[i],"=",paste(fac_items[[i]],collapse=","))
    }),model2),collapse="\n")
  }

  for(i in 2:ngrps){
    free_pars <- paste(sapply(which_factors,function(f){
      paste0("(GROUP, MEAN_",f,"), (GROUP, COV_",f,f,")")
    }),collapse=", ")
    model_iter <- model_freeMV <- paste0(model_freeMV,"\n",
                           "FREE [",grp_levels[i],"] = ",free_pars,"\n")
    # If fix_MV is true, then M/V are still freed for WALD-2 analysis; they're just fixed in future analyses at the WALD-2 values.
    if(!is.null(fix_MV)){
      fixed_pars <- paste(sapply(fix_MV,function(f){
        paste0("(GROUP, MEAN_",f,"), (GROUP, COV_",f,f,")")
      }),collapse=", ")
      model_iter <- model_freeMV # model_iter is used for the iterations—keeps fixed_pars fixed 
      model_freeMV <- paste0(model_freeMV,"\n",
                             "FREE [",grp_levels[i],"] = ",fixed_pars,"\n")
    }
  }
  
  # Start by fitting model with all parameters fixed between groups, with group mean/var free to vary
  if(bifactor==T){
    if(!is.null(method)){
      if(verbose==T){message("Only EM method can be used for bifactor models. Ignoring 'method' argument.")}
    }
    mod_fixed <- suppressMessages(bfactor(data=data,model=model,model2=model_freeMV,group=group,
                                          invariance = c("slopes","intercepts"),SE=F,SE.type=SE.type,TOL=TOL,
                                          verbose=verbose,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
  } else{
    # Default method: EM (3 or fewer factors) or 2step (>3 factors)
    if(is.null(method)){method <- ifelse(nfac > 3,"2step","EM")}
    # if stochastic method, switch to MHRM standard error
    if(method %in% c("SEM","MHRM","MCEM") & SE.type=="Oakes"){
      SE.type <- "FMHRM"
    } else if(method %in% c("2step","twostep","2STEP","TWOSTEP")){
      # If 2step is method, use MHRM method
      mod_fixed <- suppressMessages(multipleGroup(data=data,model=model_freeMV,group=group,method="MHRM",
                                                  invariance = c("slopes","intercepts"),SE=F,TOL=NULL,
                                                  verbose=verbose,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
    } else{
      mod_fixed <- suppressMessages(multipleGroup(data=data,model=model_freeMV,group=group,method=method,
                                                  invariance = c("slopes","intercepts"),SE=F,TOL=TOL,
                                                  verbose=verbose,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
    }
  }
  
  # Extract relevant information from model
  fixed_coefs <- coef(mod_fixed,simplify=T)
  mod_means <- lapply(fixed_coefs,function(GRP){GRP$means})
  mod_cov <- lapply(fixed_coefs,function(GRP){GRP$cov})
  
  # Now make model where M/V are fixed at earlier values
  # Start with base model + fixed means/variances for identification
  if(!is.null(fix_MV)){ # If there are some fixed model pars, fix both them and the free pars
    model_MV <- paste0(base_model,"\nFIXED = ",free_pars,", ",fixed_pars) # if fixed pars exist, have them fixed here as well
  } else {
    model_MV <- paste0(base_model,"\nFIXED = ",free_pars) 
  }
  
  
  for(i in 2:ngrps){
    for(f in 1:nfac){
      model_MV <- paste0(model_MV,"\n",
                         "START [",grp_levels[i],"] = (GROUP, COV_",f,f,",",as.matrix(unlist(mod_cov[[i]]))[f,f],")\n",
                         "START [",grp_levels[i],"] = (GROUP, MEAN_",f,",",(mod_means[[i]])[f],")\n")
      # Also provide these start values for model_iter
      model_iter <- paste0(model_iter,"\n",
                         "START [",grp_levels[i],"] = (GROUP, COV_",f,f,",",as.matrix(unlist(mod_cov[[i]]))[f,f],")\n",
                         "START [",grp_levels[i],"] = (GROUP, MEAN_",f,",",(mod_means[[i]])[f],")\n")
    }
  }
  
  # Fix mean and variance of focal group, then re-fit model allowing all parameters to vary (Wald-2)
  if(bifactor==T){
    mod_fixedMV <- suppressMessages(bfactor(data=data,model=model,model2=model_MV,group=group,verbose=verbose,SE.type=SE.type,TOL=TOL,
                                            SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
    
  } else{
    if(method %in% c("2step","twostep","2STEP","TWOSTEP")){
      # If 2step is method, use multipleGroup.2step function
      mod_fixedMV <- suppressMessages(multipleGroup.2step(data=data,model=model_MV,group=group,SE.type=SE.type,verbose=verbose,method=method.2step,TOL=TOL,
                                                    SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,
                                                                                    SEtol = SEtol)),...))
    } else{
      mod_fixedMV <- suppressMessages(multipleGroup(data=data,model=model_MV,group=group,SE.type=SE.type,verbose=verbose,method = method,TOL=TOL,
                                                    SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,MHRM_SE_draws=MHRM_SE_draws,
                                                                                    SEtol = SEtol)),...))  
    }  
  }
  
  fixedMV_pars <- coef(mod_fixedMV,simplify=T)[[2]]$items
  
  # Calculate Wald-2 DIF and store items that exhibit no DIF as potential anchors
  dif_wald_2 <- suppressMessages(WALD.DIF(mod_fixedMV,test_pars,p.adj = p.adj))
  # Add effect sizes to Wald-2 output (if there are 2 groups; otherwise, find in $DIF.EffectSize call)
  if(include.ES==T){
    ES.wald2 <- DIF.ES(mod_fixedMV)
    if(ngrps==2){
      if(nfac==1){ # 1 factor: use wABC, SIDS, UIDS, ESSD
        wABC <- DIF.wABC(mod_fixedMV)
        itemES <- ES.wald2[[1]]$items
        dif_wald_2 <- cbind(dif_wald_2,wABC,itemES[,c(1,2,4)])
      } else{ # more than 1 factor
        itemES <- ES.wald2[[1]]$items
        dif_wald_2 <- cbind(dif_wald_2,itemES[,c(1,2,4)])
      }
    }
  } else{
    ES.wald2 <- NULL
  }
  if(p.adj != "none"){
    non_dif_items <- which(dif_wald_2[,4] > p.cutoff) #
    unadjusted_nonDIF <- which(dif_wald_2[,3] > p.cutoff)
  } else{
    non_dif_items <- which(dif_wald_2[,3] > p.cutoff) # non-adjusted p.value
  }
  # If ES_anchors are turned on, add those to "non-dif" items
  if(ES_anchors & include.ES & ngrps==2){
    non_dif_items <- sort(unique(c(non_dif_items,which(abs(ES.wald2[[1]]$items[,4]) < ES.cutoff))))
  }
  # If all or no items show DIF, stop here.
  if(length(non_dif_items)==0){ # all items show DIF on Wald-2
    message("No anchor items! Look over data carefully.\n")
    result <- list("Iterations"=0,"DIF.model"=mod_fixedMV,"coef.DIF.model"=coef(mod_fixedMV,simplify=T),
                   "Results"=dif_wald_2,"DIF.items"=item_names,"DIF.EffectSize"=ES.wald2,
                   "p.adjust.method"=p.adj,"p.cutoff"=p.cutoff,"ES.anchors"=ES_anchors,"ES.cutoff"=ES.cutoff,"bifactor"=bifactor,"which.pars"=test_pars,
                   "Call"=match.call(),"Type"="wald.iter","DIF.model.wald2"=mod_fixedMV)
    class(result) <- "DIF.wald"

    if(verbose==T){ # if verbose, print final model (in this case, the Wald-2 model).
      cat("\n")
      print(result$DIF.model)
      cat("\n")
    }
    return(result)
  } else if(length(non_dif_items)==nitems){ # No DIF on Wald-2. Stop here
      result <- list("Iterations"=0,"DIF.model"=mod_fixedMV,"coef.DIF.model"=coef(mod_fixedMV,simplify=T),
                     "Results"=dif_wald_2,"DIF.items"=integer(0),"DIF.EffectSize"=ES.wald2,
                     "p.adjust.method"=p.adj,"p.cutoff"=p.cutoff,"ES.anchors"=ES_anchors,"ES.cutoff"=ES.cutoff,"bifactor"=bifactor,"which.pars"=test_pars,
                     "Call"=match.call(),"Type"="wald.iter","DIF.model.wald2"=mod_fixedMV)
      class(result) <- "DIF.wald"

      if(verbose==T){ # if verbose, print final model (in this case, the Wald-2 model).
        cat("\n")
        print(result$DIF.model)
        cat("\n")
      }
      return(result)
    } else { # Some items do show DIF on this. Test those using all other items as anchors.
    # Now pick the items not showing DIF
    unadjusted_nonDIF <- which(dif_wald_2[,3] > p.cutoff)
    # If ES_anchors included, add to list of anchors
    if(ES_anchors & include.ES & ngrps==2){
      unadjusted_nonDIF <- sort(unique(c(unadjusted_nonDIF,which(abs(ES.wald2[[1]]$items[,4]) < ES.cutoff))))
    }
    which_anchors <- unadjusted_nonDIF
    if(verbose){ # Print results of wald-2 before moving on
      cat(paste0("\nDIF detected in ",nitems - length(non_dif_items)," items on Wald-2 testing:\n"))
      print(as.matrix(dif_wald_2))
      cat("\n")
    }
  }
  
  # Set up variables for iterative procedure
  n.iter <- 0
  iter_done <- FALSE
  # Iterate until all items are either deemed to be anchors or show DIF
  while(iter_done==FALSE){
    # Now refit model with selected anchor items and test DIF again
    
    if(bifactor==T){
      mod_anchor <- suppressMessages(bfactor(data=data,model=model,model2=model_iter,group=group,verbose=verbose,SE.type=SE.type,
                                             invariance = c(item_names[which_anchors]),TOL=TOL,
                                             SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES)),...))
      
    } else{
      if(method %in% c("2step","twostep","2STEP","TWOSTEP")){
        # If 2step is method, use method.2step (default "QMCEM") starting with previous model's pars
        mod_anchor <- suppressMessages(multipleGroup.2step(data=data,model=model_iter,group=group,verbose=verbose,TOL=TOL,
                                                     invariance = c(item_names[which_anchors]),method=method.2step,
                                                     SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,
                                                                                     SEtol = SEtol)),...))
      } else{
        mod_anchor <- suppressMessages(multipleGroup(data=data,model=model_iter,group=group,verbose=verbose,method=method,TOL=TOL,
                                                     invariance = c(item_names[which_anchors]),
                                                     SE=T,technical=c(technical,list(set.seed=seed,NCYCLES=NCYCLES,MHRM_SE_draws=MHRM_SE_draws,
                                                                                     SEtol = SEtol)),...))
      }  
    }
    
    dif <- WALD.DIF(mod_anchor,test_pars,items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
    
    # Find which items were tested
    tested_items <- item_names[-which_anchors]
    # Add Effect sizes
    if(include.ES){
      ES <- DIF.ES(mod_anchor,which.items = (1:nitems)[-which_anchors])
      if(ngrps==2){
        if(nfac==1){ # 1 factor: use wABC, SIDS, UIDS, ESSD
          wABC <- DIF.wABC(mod_anchor,which.items = which(item_names %in% tested_items))
          itemES <- ES[[1]]$items
          dif <- cbind(dif,wABC,itemES[,c(1,2,4),drop=FALSE])
        } else{ # more than 1 factor
          itemES <- ES[[1]]$items
          dif <- cbind(dif,itemES[,c(1,2,4),drop=FALSE])
        }
      }
    } else{
      ES <- NULL
    }
    
    # Evaluate which tested items showed DIF or no DIF (uncorrected)
    if(p.adj != "none"){
      if(ES_anchors & include.ES & ngrps==2){
        non_dif_items <- which(item_names %in% tested_items[which(dif[,3] > p.cutoff | abs(ES[[1]]$items[,4]) < ES.cutoff)])
      } else{
        non_dif_items <- which(item_names %in% tested_items[which(dif[,3] > p.cutoff)])
      }
      dif_items <- item_names[which(item_names %in% tested_items[which(dif[,4] < p.cutoff)])]
    } else{
      if(ES_anchors & include.ES & ngrps==2){
        non_dif_items <- which(item_names %in% tested_items[which(dif[,3] > p.cutoff | abs(ES[[1]]$items[,3]) < ES.cutoff)])
      } else{
        non_dif_items <- which(item_names %in% tested_items[which(dif[,3] > p.cutoff)])
      }
      dif_items <- item_names[which(item_names %in% tested_items[which(dif[,3] < p.cutoff)])]
    }
    # Now add the non-DIF items to the list of anchors
    which_anchors <- c(which_anchors,non_dif_items)
    # Increase list of iterations
    n.iter <- n.iter + 1
    # If all items tested show DIF (without p-value correction), end loop
    if(ES_anchors & include.ES & ngrps==2){
      iter_done <- all(dif[,3] <= p.cutoff) & all(abs(ES[[1]]$items[,4]) >= ES.cutoff)
    } else{
      iter_done <- all(dif[,3] <= p.cutoff)
    }
    
    # Also end loop if no items show DIF 
    if(all(item_names[non_dif_items]==tested_items)){
      iter_done <- TRUE
    }
  }
  # Post-hoc tests for just slopes and just intercepts
  if(is.null(test_pars)){
    dif_slope <- WALD.DIF(mod_anchor,which.par = "slopes",
                          items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
    dif_intercept <- WALD.DIF(mod_anchor,which.par = "intercepts",
                              items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
  } else{
    dif_slope <- WALD.DIF(mod_anchor,which.par = test_pars[grep("^a",test_pars)],
                          items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
    dif_intercept <- WALD.DIF(mod_anchor,which.par = test_pars[grep("^[dc]",test_pars)],
                              items2test = (1:nitems)[-which_anchors],p.adj = p.adj)
  }
  result <- list("Iterations"=n.iter,"DIF.model"=mod_anchor,"coef.DIF.model"=coef(mod_anchor,simplify=T),
                 "Results"=dif,"DIF.items"=dif_items,"DIF.EffectSize"=ES,"p.adjust.method"=p.adj,"p.cutoff"=p.cutoff,
                 "ES.anchors"=ES_anchors,"ES.cutoff"=ES.cutoff,"bifactor"=bifactor,"Call"=match.call(),"Type"="wald.iter","which.pars"=test_pars,
                 "DIF.model.wald2"=mod_fixedMV,"Posthoc.Slope"=dif_slope,"Posthoc.Intercept"=dif_intercept)
  class(result) <- "DIF.wald"
  if(verbose==T){ # if verbose, print final model after fitting all models.
    cat("\n")
    print(result$DIF.model)
    cat("\n")
  }
  return(result)
}

print.DIF.wald <- function(x){
  cat("Call: ")
  print(x$Call)
  cat("\n")
  cat("Differential Item Functioning (DIF) using ")
  if(x$Type=="wald.iter"){cat("iterative Wald test procedure (Cao et al., 2017):\n")}
  if(x$Type=="wald.2step"){cat("two-step Wald test procedure (Wang & Woods, 2017):\n")}
  if(!is.null(x$Iterations)){
    cat(paste0(" Solution converged after ",x$Iterations," iterations"))
    if(x$Iterations==0){cat(paste0(" ('Wald-2' results presented)"))}
    cat("\n\n")
  }
  print(x$Results) # meat of the results: DIF table
  cat("\n")
  # List of items showing DIF (with adjustment if adjustment used)
  if(x$p.adjust.method!="none"){
    if(length(x$DIF.items)==0){
      cat(paste0("----------No significant DIF detected (ps > ",x$p.cutoff))
      if(x$p.adjust.method %in% c("BH","fdr")){cat(" [FDR corrected]")}
      if(x$p.adjust.method == c("bonferroni")){cat(" [Bonferroni corrected]")}
      if(x$p.adjust.method == c("holm")){cat(" [Bonferroni-Holm corrected]")}
      if(x$p.adjust.method == c("hommel")){cat(" [Hommel corrected]")}
      if(x$p.adjust.method == c("hochberg")){cat(" [Hochberg corrected]")}
      if(x$p.adjust.method == "BY"){cat(" [Benjamini-Yekutieli corrected]")}
      if(!is.null(x$ES.anchors) && x$ES.anchors){cat(paste0(" OR ESSDs < ",x$ES.cutoff))}
      cat(")----------")
    } else{
      cat(paste0("--------**Significant DIF detected in ",length(x$DIF.items)))
      if(length(x$DIF.items)==1){cat(" item")} else{cat(" items")}
      cat(paste0(" (p < ",x$p.cutoff))
      if(x$p.adjust.method %in% c("BH","fdr")){cat(" [FDR corrected]")}
      if(x$p.adjust.method == c("bonferroni")){cat(" [Bonferroni corrected]")}
      if(x$p.adjust.method == c("holm")){cat(" [Bonferroni-Holm corrected]")}
      if(x$p.adjust.method == c("hommel")){cat(" [Hommel corrected]")}
      if(x$p.adjust.method == c("hochberg")){cat(" [Hochberg corrected]")}
      if(x$p.adjust.method == "BY"){cat(" [Benjamini-Yekutieli corrected]")}
      if(!is.null(x$ES.anchors) && x$ES.anchors){cat(paste0(" AND all ESSDs > ",x$ES.cutoff))}
      cat(")**--------")
    }
  } else{ # No adjustment
    if(length(x$DIF.items)==0){
      if(!is.null(x$ES.anchors) & x$ES.anchors){
        cat(paste0("----------No significant DIF detected (ps > ",x$p.cutoff," [uncorrected] or ESSDs < ",x$ES.cutoff,")----------"))
      } else{
        cat(paste0("----------No significant DIF detected (ps > ",x$p.cutoff," [uncorrected])----------"))
      }
    } else{
      cat(paste0("--------**Significant DIF detected in ",length(x$DIF.items)))
      if(length(x$DIF.items)==1){cat(" item")} else{cat(" items")}
      if(!is.null(x$ES.anchors) & x$ES.anchors){
        cat(paste0(" (p < ",x$p.cutoff)," [uncorrected] AND all ESSDs > ",x$ES.cutoff,")**--------")
      } else{
        cat(paste0(" (p < ",x$p.cutoff)," [uncorrected])**--------")
      }
    }
  }
  # Print test-level effect sizes
  if(!is.null(x$DIF.EffectSize)){
    cat("\n\nTest-level DIF Effect Sizes:\n")
    n.contrasts <- length(x$DIF.EffectSize)
    print(t(sapply(x$DIF.EffectSize,function(ES){
      ES$test
    })))
  }
}

# wald.posthoc
# Does post-hoc single-df wald tests for each parameter after omnibus Wald test is significant
# Uses same p-value correction as original model (default is BH)
# Updated 03/01/21 to now allow for testing of specific factor slopes
wald.posthoc <- function(mod,test_pars=NULL,sig.only=F){
  if(!("DIF.wald" %in% class(mod))){stop("ERROR: model must be a 'DIF.wald' object!")}
  difmod <- mod$DIF.model
  inames <- extract.mirt(difmod,"itemnames")
  grpnames <- names(mod$coef.DIF.model)
  # Get IRT parameters for each group (and eliminate columns not estimated)
  coefs_1 <- rbind(coef(difmod,simplify=T)[[1]]$items[which(inames %in% mod$DIF.items),])
  coefs_2 <- rbind(coef(difmod,simplify=T)[[2]]$items[which(inames %in% mod$DIF.items),])
  coefs_1[which(coefs_1==0)] <- coefs_2[which(coefs_2==0)] <- NA 
  coefs_1 <- coefs_1[,!apply(coefs_1,2,function(X){all(is.na(X))}),drop=FALSE]
  coefs_2 <- coefs_2[,!apply(coefs_2,2,function(X){all(is.na(X))}),drop=FALSE]
  all_pars <- colnames(coefs_1)
  # If test_pars not provided, test what model tested
  if(is.null(test_pars)){
    if(!is.null(mod$which.pars)){
      test_pars <- mod$which.pars
    } else{
      test_pars <- all_pars
    }
  }
  cnames <- c(paste0("par.",grpnames),"par.diff","p",paste0("p.",mod$p.adjust.method),"higher.grp")
  
  dif_mat <- bind_rows(lapply(test_pars,function(p){
    if(!p %in% all_pars){return(NULL)}
    which.items <- which(!is.na(as.data.frame(coefs_1)[,p]))
    ps <- WALD.DIF(difmod,which.par=p,
                   items2test = which(inames %in% mod$DIF.items)[which.items],
                   p.adj = mod$p.adjust.method)[,3:4]
    pars_1 <- coefs_1[which.items,which(all_pars==p)]
    pars_2 <- coefs_2[which.items,which(all_pars==p)]
    
    out <- data.frame(pars_1,pars_2,pars_2 - pars_1,ps,ifelse(pars_2 - pars_1 > 0,grpnames[2],grpnames[1])) %>% 
      mutate(across(where(is.numeric),round,3))
    rownames(out) <- paste0(p,".",mod$DIF.items[which.items])
    colnames(out) <- cnames
    out$sig <- ifelse(out[,5] < 0.05, "*",ifelse(out[,4] < 0.05, ".",""))
    out
  }))
  if(sig.only){
    dif_mat <- dif_mat[which(dif_mat$sig=="*"),]
  }
  return(dif_mat)
}

# BCa Bootstrapped Confidence Intervals for the 'empirical_rxx' function
irt.rxx.boot <- function(fs,B=9999,CI.width=0.95,digits=3,seed=12345){
  rxx.obs <- empirical_rxx(fs)
  n <- nrow(fs)
  nf <- length(rxx.obs)
  alpha <- 1-CI.width
  CI.Lower <- CI.Upper <- z0 <- rep(NA,nf)
  flag <- rep(FALSE,nf) # flag factor with TRUE if CI happens to be a point estimate
  set.seed(seed)
  # perform bootstrap to generate B reliability values for each factor
  BS <- mcmapply(function(i){
    fs.b <- fs[sample(n,replace=TRUE),]
    return(empirical_rxx(fs.b))
  },1:B)
  
  # If not a matrix, make it a matrix
  if(!"matrix" %in% class(BS)){
    BS <- t(matrix(BS))
  }
  
  BS.Values <- apply(BS,1,sort)
  
  # if all bootstrap samples yield same value for rxx, use it for both ends of CI
  for(i in 1:nf){
    # bootstrap percentile bias corrected and accelerated
    # cf. Efron & Tibshirani (1993, Ch. 14) 
    # calculate bias-correction and acceleration parameters (z0 and a)
    z0[i] <- qnorm(((sum(BS.Values[,i] < rxx.obs[i]) + (sum(BS.Values[,i] == rxx.obs[i])+1)/2)/(B+1)))
    # zo calculation: cf. Efron & Tibshirani (1993, p. 186, formula 14.14) 
    # However, if ties are present, original formula overestimates bias, thus tie term added
    # See https://www.medcalc.org/manual/note-bcabootstrap.php for more details
  }
  # Jackknife (only one time for all factors)
  jk <- mcmapply(function(i){empirical_rxx(fs[-i,])},1:n)

  # If not a matrix, make it a matrix
  if(!"matrix" %in% class(jk)){
    jk <- t(matrix(jk))
  }
  
  Diff <- rowMeans(jk)-jk
  a <- rowSums(Diff^3)/(6*(sum(Diff^2))^1.5) # cf. Efron & Tibshirani (1993, p.186/15.15 and p. 328/22.29) 
  # adjust location of endpoints, cf. Efron & Tibshirani (1993, p.185/14.10) 
  alpha1 <- pnorm(z0+((z0+qnorm(alpha/2))/(1-(a*(z0+qnorm(alpha/2))))))
  alpha2 <- pnorm(z0+((z0-qnorm(alpha/2))/(1-(a*(z0-qnorm(alpha/2))))))
  # if either endpoint undefined, replace it with value for percentile CI
  alpha1[is.na(alpha1)] <- alpha/2
  alpha2[is.na(alpha2)] <- 1-(alpha/2)
  
  for(i in 1:nf){
    if(round(alpha1[i]*B)<1){
      CI.Lower[i] <- BS.Values[1,i]
    } else{
      CI.Lower[i] <- BS.Values[round(alpha1[i] * B),i]
    }
    CI.Upper[i] <- BS.Values[round(alpha2[i] * B),i]
  }
  
  result <- round(rbind(rxx.obs,CI.Lower,CI.Upper),digits)
  rownames(result) <- c("rxx.obs",paste0("CI.",alpha/2*100,"%"),paste0("CI.",(1-alpha/2)*100,"%"))
  return(result)
}
  
# Jackknife Slope Index
# Based on the work of Edwards, Houts, & Cai (2018), 
# http://dx.doi.org/10.1037/met0000121
# Zack Williams, 10/15/19
# Right now only works for unidimensional models - note doesn't match mirt implementation (I'd trust that one more)
LD.JSI <- function(obj,digits=3,bend=2.24,...){
  inames <- extract.mirt(obj, "itemnames") # item names
  # Jackknife your estimates of the model minus each item
  modk_list <- pblapply(1:extract.mirt(obj, "nitems"),function(k){
    # bifactor model - use bfactor instead
    if(obj@Options$dentype=="bfactor"){
      spec_k <- attr(extract.mirt(obj, "model"),"specific")[-k]
      mod_k <- bfactor(data=extract.mirt(obj, "data")[,-k],model=spec_k,
                       SE = T,itemtype = extract.mirt(obj, "itemtype")[-k],verbose=F,...)
    } else{
      mod_k <- mirt(data=extract.mirt(obj, "data")[,-k],model = extract.mirt(obj, "model"),
                    method = extract.mirt(obj, "method"),SE = T,
                    itemtype = extract.mirt(obj, "itemtype")[-k],verbose=F,...)
    }
    return(mod_k)
  })
  
  # Return matrix of calculated JSI values
  JSI_mat <- sapply(1:extract.mirt(obj, "nitems"),function(j){
    it_j <- extract.item(obj, j)
    aj <- it_j@par[1]
    j_name <- inames[j]
    
    # Jackknife by comparing slope with resampled slope
    JSI <- sapply(1:extract.mirt(obj, "nitems"),function(k){
      if(j==k){return(NA)}
      mod_k <- modk_list[[k]]
      it_jk <- extract.item(mod_k, which(extract.mirt(mod_k, "itemnames")==j_name))
      ajk <- it_jk@par[1]
      se_ajk <- it_jk@SEpar[1]
      return(unname((aj-ajk)/se_ajk))
    })
    
    return(JSI)
  })
  # Name matrix for better printing
  rownames(JSI_mat) <- colnames(JSI_mat) <- inames
  JSI_folded <- JSI_mat
  # Matrix of names (for vector of item pairs)
  name_mat <- sapply(1:length(inames),function(i){
    sapply(1:length(inames),function(j){
      paste0(inames[i],"|",inames[j])
    })
  })
  
  JSI_vals <- c(JSI_mat)
  # "Fold" matrix over the midline and sum JSI(j,k) and JSI(k,j)
  JSI_fold <- JSI_mat[lower.tri(JSI_mat)] + t(JSI_mat)[lower.tri(JSI_mat)]
  JSI_folded[lower.tri(JSI_mat)] <- JSI_mat[lower.tri(JSI_mat)] 
  JSI_folded <- t(JSI_folded)
  JSI_folded[lower.tri(JSI_mat)] <- JSI_mat[lower.tri(JSI_mat)] 
  names(JSI_fold) <- name_mat[lower.tri(JSI_mat)]
  # Calculate outlier JSI values using Mdn-MAD rule
  JSI_outlier <- mad_outliers(JSI_fold,bend=bend)$Vals
  
  return(list("JSI"=round(JSI_folded,digits),"Outliers"=JSI_outlier))
}

### residCheck- now for mirt and Lavaan objects (and EFA from psych)
# Takes object, prints out residuals greater than value (default 0.1, though 0.2 is another decent cutoff)
# Updated 04/08/21 to now print out average RMS residual for each item and check for outliers within that set
# Updated 08/30/21 to now average residuals from multiple groups
residCheck <- function(obj,cut=0.1,dig=3,n.large=T,pos.only=F,item.resids=F){
  if(is.null(attr(class(obj),"package"))){
    resmat <- residuals(obj)
  } else if(attr(class(obj),"package")=="lavaan"){
    resmat <- residuals(obj)$cov
  } else if(attr(class(obj),"package")=="mirt"){
    QMC <- extract.mirt(obj,"nfact")>3
    resmat <- M2(obj,na.rm = T,residmat = T,QMC=QMC)
    if(length(obj@Data$groupNames) > 1){
      group_props <- prop.table(table(obj@Data$group))
      resmat <- Reduce("+",lapply(1:length(resmat),function(i){
        return(resmat[[i]] * group_props[i])
      }))
    }
    resmat[upper.tri(resmat)] <- t(resmat)[upper.tri(resmat)]
  } else{
    resmat <- residuals(obj)
  }
  
  inames <- colnames(resmat)
  
  # Matrix of names (for vector of item pairs)
  name_mat <- sapply(1:length(inames),function(i){
    sapply(1:length(inames),function(j){
      paste0(inames[i],"|",inames[j])
    })
  })
  
  res <- resmat[lower.tri(resmat)]
  names(res) <- name_mat[lower.tri(resmat)]
  RMSR <- sqrt(mean(res^2,na.rm=T))
  
  if(pos.only){ # positive residuals only
    largeResids <- res[which(res > cut)]
  } else{ # positive and negative residuals
    largeResids <- res[which(abs(res) > cut)]
  }
  
  if(n.large){
    ## Look at how many high/low outliers from each item (only items with 2+)
    n_high_res <- sort(sapply(inames,function(i){
      length(grep(paste0("\\|",i,"$"),names(largeResids))) + length(grep(paste0("^",i,"\\|"),names(largeResids)))
    }),decreasing=T)
    n_high_res <- n_high_res[n_high_res > 1]
    if(length(n_high_res)==0){n_high_res <- NULL}
  } else{
    n_high_res <- NULL
  }
  
  itemResids <- apply(resmat,1,function(X){sqrt(mean(X^2,na.rm=T))})

  result <- list("Residuals"=round(resmat,dig),"Large.Resids"=round(largeResids,dig),"RMSR"=round(RMSR,dig),"N.Large.Resids"=n_high_res)
  if(item.resids){ # Include item-level RMSR in output
    result$Item.RMSR <- round(itemResids,dig)
    item_outliers <- which((itemResids-hd(itemResids))/(1.4826*hd(abs(itemResids-hd(itemResids)))) > 2.24)
    if(length(item_outliers)>0){
      result$RMSR.Outliers <- round(itemResids,dig)[item_outliers]
    }
    result$Cutoff <- cut
  } else{
    item_outliers <- which((itemResids-hd(itemResids))/(1.4826*hd(abs(itemResids-hd(itemResids)))) > 2.24)
    if(length(item_outliers)>0){
      result$RMSR.Outliers <- round(itemResids,dig)[item_outliers]
    }
    result$Cutoff <- cut
  }
  return(result)
}

### resid_cor - get correlation residuals from mirt object
resid_cor <- function(obj,cut=0.1,dig=3,n.large=T,pos.only=F,item.resids=F){
  if(attr(class(obj),"package")!="mirt"){
    stop("Object must be a 'mirt' model!")
  }
  
  if(obj@Data$ngroups > 1){
    resmat <- lapply(1:obj@Data$ngroups,function(i){
      grp_i <- extract.group(obj,i)
      cormat <- cor_auto(grp_i@Data$data,ordinalLevelMax = max(obj@Data$K)+2,verbose = F)
      cormat[cormat==0] <- NA
      rescors <- cormat - grp_i@Fit$F %*% Phi(grp_i) %*% t(grp_i@Fit$F)
      diag(rescors) <- NA
      return(rescors)
    })
    
    group_props <- prop.table(table(obj@Data$group))
    resmat <- Reduce("+",lapply(1:length(resmat),function(i){
      return(resmat[[i]] * group_props[i])
    }))
  } else { # Just one group
    cormat <- cor_auto(obj@Data$data,ordinalLevelMax = max(obj@Data$K)+2,verbose = F)
    cormat[cormat==0] <- NA
    resmat <-  cormat - obj@Fit$F %*% Phi(obj) %*% t(obj@Fit$F)
    diag(resmat) <- NA
  }
  
  inames <- colnames(resmat)
  
  # Matrix of names (for vector of item pairs)
  name_mat <- sapply(1:length(inames),function(i){
    sapply(1:length(inames),function(j){
      paste0(inames[i],"|",inames[j])
    })
  })
  
  res <- resmat[lower.tri(resmat)]
  names(res) <- name_mat[lower.tri(resmat)]
  CRMR <- sqrt(mean(res^2,na.rm=T))
  
  if(pos.only){ # positive residuals only
    largeResids <- res[which(res > cut)]
  } else{ # positive and negative residuals
    largeResids <- res[which(abs(res) > cut)]
  }
  
  if(n.large){
    ## Look at how many high/low outliers from each item (only items with 2+)
    n_high_res <- sort(sapply(inames,function(i){
      length(grep(paste0("\\|",i,"$"),names(largeResids))) + length(grep(paste0("^",i,"\\|"),names(largeResids)))
    }),decreasing=T)
    n_high_res <- n_high_res[n_high_res > 1]
    if(length(n_high_res)==0){n_high_res <- NULL}
  } else{
    n_high_res <- NULL
  }
  
  itemResids <- apply(resmat,1,function(X){sqrt(mean(X^2,na.rm=T))})
  
  result <- list("Residuals"=round(resmat,dig),"Large.Resids"=round(largeResids,dig),"CRMR"=round(CRMR,dig),"N.Large.Resids"=n_high_res)
  if(item.resids){ # Include item-level CRMR in output
    result$Item.CRMR <- round(itemResids,dig)
    item_outliers <- which((itemResids-hd(itemResids))/(1.4826*hd(abs(itemResids-hd(itemResids)))) > 2.24)
    if(length(item_outliers)>0){
      result$CRMR.Outliers <- round(itemResids,dig)[item_outliers]
    }
  }
  result$Cutoff <- cut
  return(result)
}

### Infoplot - plot test information for a test (as mirt plot for "info" not working at the moment)
# Zack Williams updated 05/08/21 - now allows for fully marginal reliability (doesn't marginalize information in plot though)
# Note: only tests one factor at a time for multidimensional models (default: factor 1)
# Now also includes reliability calculations (use "reliability=T" or 'reliabilityplot()' wrapper) 
# Can change latent density using density argument and ... for inputs to density
infoplot <- function(mod,Theta.lim=c(-6,6),which.items=1:extract.mirt(mod,"nitems"),conditional=TRUE,
                     linecolor="darkblue",linewidth=1,areafill="dodgerblue3",areaalpha=0.4,
                     which.factor=1,testname=deparse(substitute(mod)),
                     factorname=extract.mirt(mod,"factorNames")[which.factor],density=dnorm,
                     reliability=FALSE,rel.cutoff=0.7,rel.cutoff.color="darkred",rel.cutoff.display=TRUE,...){
  nfact <- extract.mirt(mod,"nfact")
  # Creates thetas in one dimension with all other theta vals set to 0
  thetas <- do.call(cbind,c(rep(list(0),which.factor-1),
                            list(seq(Theta.lim[1],Theta.lim[2],length.out = 1000)),
                            rep(list(0),extract.mirt(mod,"nfact")-which.factor)))
  if(nfact==1){
    degrees <- NULL
    plot_title <- bquote(bold("Test Reliability for"~.(testname)))
  } else{
    degrees <- rep(90,nfact)
    degrees[which.factor] <- 0
    plot_title <- bquote(bold("Test Reliability for"~.(testname)["["*.(factorname)*"]"]))
  }
  info <- testinfo(mod,Theta=thetas,degrees=degrees,which.items = which.items)
  
  if(reliability==T){
    rel <- info/(info + 1)
    maxrel <- which.max(rel) # index where maximum reliability occurs
    if(rel[maxrel] <= rel.cutoff){ # no reliability values are high enough to be over cutoff
      rel_roots <- NULL
    } else{ # Reliability cutoff intersects curve
      rel_roots <- c(thetas[which.min(abs(rel[1:maxrel] - rel.cutoff)),which.factor],
                     thetas[maxrel + which.min(abs(rel[(maxrel+1):length(rel)] - rel.cutoff)),which.factor])
    }
    
    if(conditional==T | nfact==1){ # Use 1-dimensional integral (default; much faster)
      # Information evaluation function
      fn <- function(theta, mod, den, ...){
        th <- do.call(cbind,c(rep(list(0),which.factor-1),
                              list(theta),
                              rep(list(0),extract.mirt(mod,"nfact")-which.factor)))
        th[,which.factor] <- theta
        TI <- testinfo(mod, th,degrees=degrees, which.items=which.items)
        TI / (TI + 1) * den(theta, ...)
      }
      marg_rxx <- round(integrate(fn, lower = -Inf, upper=Inf, mod=mod, den=density)$value,3) # Marginal Reliability 
    } else { # Fully marginal reliability (integrates over all dimensions)
      # Note that this ignores the density function, as only multivariate normal is supported
      fn <- function(theta_mat,mod){
        mod_coefs <- coef(mod,simplify=T)
        TI <- testinfo(mod, t(theta_mat),degrees=degrees)
        matrix(TI / (TI + 1) * pbapply(theta_mat,2,function(th){mvnfast::dmvn(th,mu=mod_coefs$mean,sigma = mod_coefs$cov)}))
      }
      message("Performing multidimensional integration. May take a while...")
      marg_rxx <- round(hcubature(fn, lower = matrix(rep(Theta.lim[1],nfact)), upper = matrix(rep(Theta.lim[2],nfact)),
                                  mod=mod, absError = 1e-3, vectorInterface = T)$integral,3) # Marginal Reliability 
    }
    
    marg_label <- bquote("Marginal"~rho[xx]~"="~.(marg_rxx))

    p <- ggplot(data=data.frame("Theta"=thetas[,which.factor],"Reliability"=rel),aes(x=Theta,y=Reliability)) + 
      geom_line(color=linecolor,size=linewidth) + geom_area(alpha=areaalpha,fill=areafill) +
      ggtitle(plot_title) +
      theme_bw() + xlab(expression(bold("Latent Trait Score"~(theta)))) + ylab(expression(bold(r[xx](theta)))) + 
      annotate("label",x = Theta.lim[1],y=1,label=deparse(marg_label),size=3,parse=TRUE,hjust="left") + 
      theme(legend.position="right",
            plot.title = element_text(size=14, face="bold", hjust=0.5),
            legend.title = element_text(size=10, face="bold", hjust=0.5),
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11),
            axis.title.x = element_text(size=12,face="bold"),
            axis.title.y = element_text(size=12,face="bold"),
            strip.text.y = element_text(size = 20, colour = "white"),
            strip.background = element_rect(colour="darkgrey", fill="darkgrey")) + 
      coord_cartesian(xlim=Theta.lim,ylim=c(0,1))
    
    if(rel.cutoff.display==T){
      p <- p + geom_hline(data=NULL,yintercept=rel.cutoff,color=rel.cutoff.color,linetype="dashed",alpha=0.7) + 
        annotate("label",x = Theta.lim[1],y=rel.cutoff,label=deparse(bquote(r[xx]~"="~.(rel.cutoff))),parse=TRUE,
                 size=3,color=rel.cutoff.color,hjust=0.3)
      if(length(rel_roots)==2){
        p <- p + annotate("label",x = rel_roots[1],y=rel.cutoff,label=paste0(round(rel_roots[1],2)),size=3) + 
          annotate("label",x = rel_roots[2],y=rel.cutoff,label=paste0(round(rel_roots[2],2)),size=3)
      }
    }
    
  } else{
    p <- ggplot(data=data.frame("Theta"=thetas[,which.factor],"Info"=info),aes(x=Theta,y=Info)) + 
      geom_line(color=linecolor,size=linewidth) + geom_area(alpha=areaalpha,fill=areafill) +
      ggtitle(paste0("Test Information for ",testname)) +
      theme_bw() + xlab(expression(bold("Latent Trait Score ("*theta*")"))) + ylab("Test Information") + 
      theme(legend.position="right",
            plot.title = element_text(size=14, face="bold", hjust=0.5),
            legend.title = element_text(size=10, face="bold", hjust=0.5),
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11),
            axis.title.x = element_text(size=12,face="bold"),
            axis.title.y = element_text(size=12,face="bold"),
            strip.text.y = element_text(size = 20, colour = "white"),
            strip.background = element_rect(colour="darkgrey", fill="darkgrey"))
  }
  return(p)
}
### reliabilityplot - alias for infoplot (see above function for more details)
reliabilityplot <- function(mod,testname=deparse(substitute(mod)),...){
  infoplot(mod = mod,reliability=TRUE,testname=testname,...)
}
### Iteminfoplot
## Zack Williams, 06/26/20
## Like infoplot, except superimposes all items on top of each other
## Has the option to superimpose itemplot
## If testinfo desired but you don't want the shading, use test.areafill=NA
## By default looks at all items and uses loadings along first factor. 
## Items with 0 loadings need to be excluded manually (though the information should be a pretty flat line at 0)
iteminfoplot <- function(mod,Theta.lim=c(-6,6),which.items=1:extract.mirt(mod,"nitems"),
                         itemnames=extract.mirt(mod,"itemnames")[which.items],which.factor=1,
                         linewidth=1,linealpha=0.8,title.label="Item Information Plot",
                         test.info=F,test.linecolor="darkblue",test.linewidth=1.1,
                         which.items.testinfo=which.items,
                         test.areafill="dodgerblue3",test.areaalpha=0.2){
  nfact <- extract.mirt(mod,"nfact")
  
  # Creates thetas in one dimension with all other theta vals set to 0
  thetas <- do.call(cbind,c(rep(list(0),which.factor-1),
                            list(seq(Theta.lim[1],Theta.lim[2],length.out = 1000)),
                            rep(list(0),extract.mirt(mod,"nfact")-which.factor)))
  if(nfact==1){
    degrees <- NULL
  } else{
    degrees <- rep(90,nfact)
    degrees[which.factor] <- 0
  }
  
  infomat <- cbind(sapply(which.items,function(i){
    testinfo(mod,Theta=thetas,degrees=degrees,which.items = i)
  }))
  colnames(infomat) <- itemnames
  imat_wide <- data.frame(cbind("Theta"=thetas[,which.factor],infomat))
  # Make long-form usinv tidyr::pivot_longer
  imat_long <- pivot_longer(imat_wide,cols=2:ncol(imat_wide),names_to="Item",values_to = "Info")
  
  p <- ggplot(data=imat_long,aes(x=Theta,y=Info,color=Item,fill=Item)) + 
    ggtitle(title.label) +
    theme_bw() + xlab(expression(bold("Latent Trait Score"~(theta)))) + ylab("Information") + 
    theme(legend.position="right",
          plot.title = element_text(size=14, face="bold", hjust=0.5),
          legend.title = element_text(size=10, face="bold", hjust=0.5),
          axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=11),
          axis.title.x = element_text(size=12,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          strip.text.y = element_text(size = 20, colour = "white"),
          strip.background = element_rect(colour="darkgrey", fill="darkgrey"))
  
  # Add test information curve if desired
  if(test.info){
    test_info <- testinfo(mod,Theta=thetas,degrees=degrees,which.items = which.items.testinfo)
    p <- p + geom_line(data=data.frame("Theta"=thetas[,which.factor],"Info"=test_info),
                       aes(x=Theta,y=Info),inherit.aes = F,color=test.linecolor,
                       size=test.linewidth) + 
      geom_area(data=data.frame("Theta"=thetas[,which.factor],"Info"=test_info),
                aes(x=Theta,y=Info),inherit.aes = F,alpha=test.areaalpha,fill=test.areafill) +
      geom_line(size=linewidth,alpha=linealpha)
  }
  
  ## Add the lines for each item on top of the shading
  p <- p + geom_line(size=linewidth,alpha=linealpha)
  
  # Now return plot
  return(p)
}


### fscores.pv - computes summary statistics using posterior distribution of ability estimates (plausible values)
## Zack Williams, 06/05/2020
# By default, draws 2000 PV matrices, calculates mean, median, mode, 95% HDI, SE, and empirical reliability
# PVs drawn with Metropolis-Hastings algorithm, but plausible.type="normal" uses a normal approximation and may be faster
# Summary estimates not exactly the same as those from fscores(), but correlation > 0.99
fscores.pv <- function(obj,plausible.draws=2000,plausible.type="MH",
                       ests="all",ci=0.95,ci.type="HDI",round=3,cores=detectCores()-1,seed=1234,...){
  nf <- extract.mirt(obj,"nfact")
  fnames <- extract.mirt(obj,"factorNames")
  message("Generating plausible values... ")
  pv_list_raw <- fscores(obj,plausible.draws=plausible.draws,
                         plausible.type=plausible.type,technical=list(set.seed=seed),...)
  message("Done! Calculating scores for each factor individually...")
  pvs <- do.call(cbind,pv_list_raw) # collapse PVs into one matrix (1 row per subj.)
  # Meat of the function
  # Loop through factors: for each, pull together full distribution, then find individual summary vals
  result <- lapply(fnames,function(f){
    n <- which(fnames==f)
    cols <- seq(n,ncol(pvs),by = nf) # gather all columns of this factor
    fscore_dists <- pvs[,cols]
    scores <- pbsapply(1:nrow(fscore_dists),function(i){
      theta <- fscore_dists[i,] # All vals in a row (theta scores for one individual)
      quantiles <- c(50*(1-ci),100-50*(1-ci))
      if(ci.type %in% c("HDI","hdi")){
        CI <- unlist(bayestestR::ci(theta,ci=ci,method="HDI"))[-1] # CrI (95% HDI by default)
        names(CI) <- paste0("HDI.",quantiles,"%")
      } else{ # anything else, report ETI (note SI not supported because questionable utility?)
        CI <- unlist(bayestestR::ci(theta,ci=ci,method="ETI"))[-1]
        names(CI) <- paste0("Q.",quantiles,"%") # Name ETI estimates differently so they stick out
      }
      summ <- c(unlist(bayestestR::point_estimate(theta,centrality=ests)),CI,"SE"=sd(theta),"r_xx"=max(1-var(theta),0))
    },cl=cores)
    if(!is.null(round)){ # default: rounds to 3 digits
      scores <- round(t(scores),round)
    } else{ # use round=NULL to return unrounded vals
      scores <- t(scores)
    }
    return(scores)
  })
  names(result) <- fnames
  result$Reliability <- t(sapply(1:nf,function(i){
    X <- result[[i]]
    return(round(summary(X[,ncol(X)]),round))
  }))
  rownames(result$Reliability) <- fnames
  return(result)
}

## marginal_rxx
# Calculates marginal reliability for multidimensional models
# Zack Williams, updated 05/08/21 - now allows for fully marginal reliability
# Updated 01/15/22 to include selection of subset of items
marginal_rxx <- function(mod,density=dnorm,Theta.lim=c(-6,6),dig=3,conditional=F,which.items = 1:extract.mirt(mod, "nitems"),cores=detectCores()-1){
  nfact <- extract.mirt(mod,"nfact") # number of factors
  
  # Function to integrate over (density-weighted reliability values)
  if(conditional==T | nfact==1){ # Use 1-dimensional integral (default; much faster)
    # Information evaluation function
    fn <- function(theta, mod, den, which.factor,...){
      degrees <- rep(90,extract.mirt(mod,"nfact"))
      degrees[which.factor] <- 0
      # Creates thetas in one dimension with all other theta vals set to 0
      th <- do.call(cbind,c(rep(list(0),which.factor-1),
                            list(theta),
                            rep(list(0),extract.mirt(mod,"nfact")-which.factor)))
      th[,which.factor] <- theta
      TI <- testinfo(mod, th,degrees=degrees,which.items = which.items)
      TI / (TI + 1) * den(theta, ...)
    }
  } else { # Fully marginal reliability (integrates over all dimensions)
    # Note that this ignores the density function, as only multivariate normal is supported
    fn <- function(theta_mat,mod, degrees, which.factor){
      degrees <- rep(90,extract.mirt(mod,"nfact"))
      degrees[which.factor] <- 0
      mod_coefs <- coef(mod,simplify=T)
      TI <- testinfo(mod, t(theta_mat),degrees=degrees, which.items = which.items)
      matrix(TI / (TI + 1) * pbapply(theta_mat,2,function(th){mvnfast::dmvn(th,mu=mod_coefs$mean,sigma = mod_coefs$cov)}))
    }
  }
  # Loop through factors and calculate marginal reliability
  ## Conditional or 1 factor (for 1D integration)
  if(conditional==TRUE | nfact==1){
    marg_rxx <- sapply(1:nfact,function(which.factor){
      return(integrate(fn, lower = -Inf, upper=Inf, mod=mod, den=density, which.factor=which.factor)$value)
    })
  } else{ # Full marginal (does multidimensional integration)
    message("Performing multidimensional integration. May take a while...")
    marg_rxx <- pbsapply(1:nfact,function(which.factor){
      return(hcubature(fn, lower = matrix(rep(Theta.lim[1],nfact)), upper = matrix(rep(Theta.lim[2],nfact)),
                       mod=mod, absError = 1e-3, vectorInterface = T, which.factor=which.factor)$integral) # Marginal Reliability 
    },cl=cores)
  }
  
  names(marg_rxx) <- extract.mirt(mod,"factorNames")
  # Now return final product
  if(is.na(dig) | is.null(dig)){ # if dig is NA or NULL, don't round
    return(marg_rxx)
  } else{
    return(round(marg_rxx,dig))
  }
}


## marginal_rxx.boot: percentile bootstrapped CI for Marginal rxx
## Zack Williams, 07/04/2020
# BCa not implemented at this time
# default 2000 bootstrapped replications
marginal_rxx.boot <- function(mod,density=dnorm,B=1000,ci=0.95,dig=3,cores=parallel::detectCores()-1,seed=12345,...){
  set.seed(seed)
  rawdata <- mod@Data$data
  rxx_raw <- marginal_rxx(mod,density=density,...)
  modpars <- mod2values(mod)
  
  ci.low <- (1-ci)/2
  ci.hi <- 1 - ci.low
  
  bs_vals <- pbsapply(1:B,function(b){
    newdata <- rawdata[sample(1:nrow(rawdata),replace = T),]
    # Add in prior if specified
    if(length(mod@Model$parprior)>0){
      pr <- mod@Model$parprior
    } else{
      pr <- NULL
    }
    bmod <- mirt(newdata,model=mod@Model$model,pars=modpars,parprior=pr,
                 itemtype = mod@Model$itemtype,verbose = F)
    return(marginal_rxx(bmod,density = density))
  },cl=cores)
  
  CI <- quantile(bs_vals,probs = c(ci.low,ci.hi))
  names(CI) <- paste0("CI.",names(CI))
  return(round(c("rxx"=rxx_raw,CI),dig))
}

### mi.missForest: Multiple Imputation using MissForest
# Takes a dataframe, uses missForest to impute values, returns list of dataframes
# default 10 multiple imputations
# Automatically takes characters and turns them to factors. 
# ordered = vector of column names/numbers that should be turned into ordered categories 
# ... - input to actual missForest function
mi.missForest <- function(DF,m=10,ordered=NULL,cores=detectCores()-1,seed=12345,file=NULL,...){
  # If given "ordered" as column names, convert to column numbers
  if(class(ordered)=="character"){
    ordered <- sapply(ordered,function(X){grep(paste0("^",X,"$"),names(DF))})
  }
  # Search for character variables, turn to factors
  for(i in 1:ncol(DF)){
    if((i %in% ordered | class(DF[,i])=="character") & length(unique(DF[,i])) <= 53){
      # Since missForest can't handle categorical predictors with >53 categories, skip those columns
      DF[,i] <- as.factor(DF[,i])
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
    saveRDS(result,file = file)
  }
  return(result)
}


#### mirt.robust / multipleGroup.robust / bfactor.robust
### Estimates mirt model using robust (weighted) maximum marginal likelihood (weighted by person fit stats)
### Based on Hong, M. R., & Cheng, Y. (2019). Robust maximum marginal likelihood 
# (RMML) estimation for item response theory models. Behavior Research Methods, 51(2), 573-588.
# https://link.springer.com/article/10.3758/s13428-018-1150-4
### Needs complete data to work - uses missForest to impute values in cases where input has missing data
### Note that this will also make the mod fit less well using C2, etc., as it causes misfit for the "bad" people.
mirt.robust <- function(data,fscore.method="MAP",...){
  if(any(is.na(data))){ # Missing values, need to impute
    message("Missing values detected. Imputing data using missForest...")
    data <- data.frame(apply(quiet(mi.missForest(data,m=1,ordered=colnames(data)))[[1]],2,as.numeric))
    message(" Done!\n")
  }
  mod <- mirt(data,...)
  # Calculating person fit statistics
  per.fit <- personfit(mod, method = fscore.method)$Zh
  # Calculating weight
  weight <- pnorm(per.fit)*nrow(data)/sum(pnorm(per.fit))
  # Robust estimation
  robust.mod <- mirt(data,...,survey.weights = weight)
  # Add person fit statistics to object output
  robust.mod@Fit$PersonFit <- personfit(robust.mod, method = fscore.method)$Zh 
  return(robust.mod)
}

multipleGroup.robust <- function(data,fscore.method="MAP",...){
  if(any(is.na(data))){ # Missing values, need to impute
    message("Missing values detected. Imputing data using missForest...")
    data <- data.frame(apply(quiet(mi.missForest(data,m=1,ordered=colnames(data)))[[1]],2,as.numeric))
    message(" Done!\n")
  }
  mod <- multipleGroup(data,...)
  # Calculating person fit statistics
  per.fit <- personfit(mod, method = fscore.method)$Zh
  # Calculating weight
  weight <- pnorm(per.fit)*nrow(data)/sum(pnorm(per.fit))
  # Robust estimation
  robust.mod <- multipleGroup(data,...,survey.weights = weight)
  # Add person fit statistics to object output
  robust.mod@Fit$PersonFit <- personfit(robust.mod, method = fscore.method)$Zh 
  return(robust.mod)
}

bfactor.robust <- function(data,fscore.method="MAP",...){
  if(any(is.na(data))){ # Missing values, need to impute
    message("Missing values detected. Imputing data using missForest...")
    data <- data.frame(apply(quiet(mi.missForest(data,m=1,ordered=colnames(data)))[[1]],2,as.numeric))
    message(" Done!\n")
  }
  mod <- bfactor(data,...)
  # Calculating person fit statistics
  per.fit <- personfit(mod, method = fscore.method)$Zh
  # Calculating weight
  weight <- pnorm(per.fit)*nrow(data)/sum(pnorm(per.fit))
  # Robust estimation
  robust.mod <- bfactor(data,...,survey.weights = weight)
  # Add person fit statistics to object output
  robust.mod@Fit$PersonFit <- personfit(robust.mod, method = fscore.method)$Zh 
  return(robust.mod)
}

# Probability of obtaining two significantly different test scores for IRT score estimates (PTDS_IRT)
# Zack Williams, 10/04/2020
# Based on a hella-obscure paper by Müller (2006): https://doi.org/10.1177/0013164405284034
# Calculates the percentage of scores 
### Interpretation guidelines:
# ‘very poor’ for <30%
# ‘poor’ between 30% and 45%
# ‘moderately’ for 45% up to 60%
# ‘good’ for 60% up to 75%
# ‘very good’ for 75% up to 90% 
# 'excellent' for >90% 
PTDS.irt <- function(scores,SEs,bootstrap=FALSE,ci=0.95,B=2000,cores=detectCores()-1){
  # Check to make sure input is good
  if(class(scores) != "numeric" | class(SEs) != "numeric"){
    stop("Both 'scores' and 'SEs' arguments should be numeric vectors.")
  }
  if(length(scores) != length(SEs)){
    stop("Input vectors are not the same length!")
  }
  # Find NA values and remove them
  not.na <- !(is.na(scores) | is.nan(scores) | is.na(SEs) | is.nan(SEs))
  
  scores <- scores[not.na]
  SEs <- SEs[not.na]
  combos <- combn(1:length(scores),2)
  # Go through combinations, check for significant differences
  out <- pbapply(combos,2,function(i){
    theta1 <- scores[i[1]]
    theta2 <- scores[i[2]]
    SE1 <- SEs[i[1]]
    SE2 <- SEs[i[2]]
    
    abs((theta1 - theta2)/sqrt(SE1^2 + SE2^2)) > 1.96
  },cl = cores)
  if(bootstrap){
    message("Bootstrapping...")
    bs <- pbsapply(1:B,function(X){mean(sample(out,length(out),replace=T))},cl=cores)
    CI <- quantile(bs,c((1-ci)/2,1-(1-ci)/2))
  } else{
    CI <- NULL
  }
  return(round(c("PTDS.IRT"=mean(out),CI),3))
}


