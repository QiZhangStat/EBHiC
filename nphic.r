## function for creating B-spline
bsplfun <- function(xrange=c(0,1),breaks=seq(xrange[1],xrange[2],length.out=100),k=4,ncores=1){ 
          nbreaks <- length(breaks)
          nbasis <- nbreaks+k-2
          step.break <- mean(abs(diff(breaks)))
          breaks.extend <- c(breaks[1]-step.break*((k-1):1),breaks,breaks[nbreaks]+step.break*(1:(k-1)))
          coef0 <- matrix(1,1,1)
          bspl <- lapply(1:length(breaks.extend),function(ii) list(breaks=breaks.extend[ii:(ii+1)],coef=coef0))
          for(kk in 2:k){
                    bspl.kk <- mclapply(1:(length(bspl)-1), function(ii) bsplfun.updt(ii,bspl),mc.cores=ncores,mc.preschedule=FALSE)
                            bspl <- bspl.kk
                  }
           out <- bspl[1:nbasis]
           for(ii in 1:(k-1)){
              outii <- out[[ii]]
              breaksii <- outii[['breaks']][(k+1-ii):(k+1)]
              coefii <- outii[['coef']][,(k-ii+1):k]
              if(ii==1){coefii <- matrix(coefii,ncol=1)}
              out[[ii]] <- list(breaks=breaksii,coef=coefii)
            }
           for(ii in 1:(k-1)){
              outii <- out[[nbasis+1-ii]]
              breaksii <- outii[['breaks']][1:(ii+1)]
              coefii <- outii[['coef']][,1:ii]
              if(ii==1){coefii <- matrix(coefii,ncol=1)}
              out[[nbasis+1-ii]] <- list(breaks=breaksii,coef=coefii)
            }
           for(ii in 1:nbasis){
              outii <- out[[ii]]
              breaksii <- outii[['breaks']]
              coefii <- outii[['coef']]
              intii <- sum(sapply(1:(length(breaksii)-1), function(jj) sum((breaksii[jj+1]^(1:k)-breaksii[jj]^(1:k))*coefii[,jj]/(1:k))))
              coefii <- coefii/intii
              out[[ii]] <- list(breaks=breaksii,coef=coefii)
            }
           out <- out[-c(1,length(out))] ### remove the boundary basis.
           attributes(out) = list(breaks=breaks,k=k)
           return(out)
        }

bsplfun.updt <- function(ii,bspl){
      bs1 <- bspl[[ii]]
          breaks1 <- bs1[['breaks']]
          coef1 <- bs1[['coef']]
          bs2 <- bspl[[ii+1]]
          breaks2 <- bs2[['breaks']]
          coef2 <- bs2[['coef']]
		  ## we assume that length(breaks1) and length(breaks2) are the same
          kk <- length(breaks1)
		  ## we assume that breaks1[-1] == breaks2[1:(kk-1)]
          breaks <- c(breaks1[1],breaks2) 
          coef <- rbind(0,cbind(coef1,0))/(breaks1[kk]-breaks1[1]) - rbind(cbind(coef1,0),0)*breaks1[1]/(breaks1[kk]-breaks1[1]) -rbind(0,cbind(0,coef2))/(breaks2[kk]-breaks2[1]) + rbind(cbind(0,coef2),0)*breaks2[kk]/(breaks2[kk]-breaks2[1])
          out <- list(breaks=breaks,coef=coef)
          return(out)
    }


## perform integration over one spline basis fnction
iPoisBsFun <- function(bk,x,cf,t=1){
  k <- length(cf)
  ti = 1/t
  out <- sapply(1:k, function(jj) cf[jj]*t*exp(lgamma(x+t*jj)-lgamma(x+1))*(pgamma(bk[2]^ti,shape=t*jj+x)-pgamma(bk[1]^ti,shape=t*jj+x)))
  return(out)
}


## integration of Binomial + B spline
intPoisBspl <- function(bs,x,t=1){ 
  breaks <- bs[['breaks']]
  cf <- bs[['coef']]
  out <- sapply(1:(length(breaks)-1), function(ii) sum(iPoisBsFun(breaks[ii:(ii+1)],x,cf[,ii],t)))
  return(sum(out))
}


## get the design matrix c_ij
getDesignMtx <- function(bs,x,t=1,ncores=1){ 
nbasis <- length(bs)
ndt <- length(x)
out <- mclapply(1:ndt, function(ii) sapply(bs, function(jj) intPoisBspl(jj,x[ii],t)),mc.cores=ncores,mc.preschedule=FALSE)
out <- t(as.matrix(as.data.frame(out)))
dimnames(out) <- NULL
return(out)
}


## calculate posterior non-central moment of lambda, using the design matrix and weights (of spline components) as the input.
getMomt <- function(pi,d0,dk,x,m=1){
    fct = sapply(x, function(xx) exp(sum(log(xx+1:m))))
    out = fct*(dk%*%matrix(pi,length(pi),1))/(d0%*%matrix(pi,length(pi),1))
    return(out)
}


## NP-Pois model. We remove the basis that only covers bin at the boundary
emPoisBspl <- function(x,t=1,bspl=NULL, pi.init=NULL,ncores=1,err.max=0.000001,iter.max=200){
  n <- length(x)
  if(is.null(bspl)){
    breaks=seq(min(x)^t,max(x)^t,length.out=101)
    k = 4
	## create the break points and coefficients of the B spline functions
    bspl <- bsplfun(range(breaks),breaks=breaks,k=k,ncores=ncores) 
  }else{
    k = attr(bspl,'k')
  }
  nb <- length(bspl)
  if(is.null(pi.init)){
    pi.init <- hist(x^t,breaks=seq(min(x)^t,max(x)^t,length.out=nb+1))$counts/n
  }
  ## increase the range by 5 to facilitate calculation of moments
  xall = 0:(max(x)+5)  
  dmtxall = getDesignMtx(bspl,xall,t,ncores)
  ## getDesignMtx(bspl,x,ncores=ncores) ### get the design matrix
  dmtx <- dmtxall[x+1,] 
  err <- Inf
  iter <- 0
  ll.init <- -Inf
  ll.all <- ll.init
  err.all <- err

  while(err>err.max&iter<iter.max){
    dtot <-  dmtx%*%matrix(pi.init,length(pi.init),1)
    post <- dmtx*(matrix(pi.init,n,nb,byrow=TRUE))/matrix(dtot,n,nb,byrow=FALSE)
    pi <- colSums(post)
    pi <- pi/sum(pi)
	## use the variable triangular kernel to smooth pi
    kern = vTriAng(pi)
    pi = kern%*%pi
    ll <- sum(apply(dmtx,1, function(ii) log(sum(ii*pi))))
    err <- ll-ll.init
    err <- max(err, err/abs(ll))
    ll.all <- c(ll.all,ll)
    err.all <- c(err.all,err)
    pi.init <- pi
    ll.init <- ll
    iter <- iter+1
  }
  ### calculate the posterior mean and variances of lambda.
  ## get the design matrix
  d1 <- dmtxall[x+2,]
  d2 <- dmtxall[x+3,]
  lambda.post.mean = getMomt(pi.init,dmtx,d1,x,m=1)
  lambda.post.var = getMomt(pi.init,dmtx,d2,x,m=2)-lambda.post.mean^2
  
  out <- list(lambda.mean=lambda.post.mean,lambda.var=lambda.post.var,x=x,f=as.numeric(dmtx%*%pi.init),xall=xall,fall=as.numeric(dmtxall%*%pi.init),pi=as.numeric(pi.init),post=post,bspl=bspl,kern=kern,dmtx=dmtx,dmtxall=dmtxall,ll.all=ll.all,err.all=err.all, convergence=list(err=err,err.max=err.max,converged=(err<err.max),niter=iter,ll=ll.init),controls=list(k=k,nbasis=nb,breaks=breaks))
  return(out)
}


## triangular kernel on 1:nl with variable range a such that the probability mass to be smoothed over is at least wthr. The "bandwidth" h could be fixed at 1.
vTriAng <- function(w,wthr=1/length(w),h=1){
  nl = length(w)
  sumwindow <- function(w,ii,a){
    wd = max(1,ii-a):min(ii+a,nl)
    sum(w[wd])
  }
  triAng <- function(ii,a){
    kii = rep(0,nl+2*a)
    val = (a+1)^h-c(a:0,1:a)^h
    kii[ii:(ii+2*a)] = val
    kii = kii[(a+1):(nl+a)]
    kii = kii/sum(kii)
    kii
  }
  avec = sapply(1:nl, function(ii) which(sapply(1:(nl-1), sumwindow, w=w,ii=ii)>=wthr)[1])
  out = matrix(sapply(1:nl, function(ii) triAng(ii,avec[ii])),ncol=nl,nrow=nl,byrow=F)
  return(out)  
}


## triangular kernel on 1:nl with fixed range. The "bandwidth" h could be fixed at 1.
fTriAng <- function(nl,a=1,h=1){
  triAng <- function(ii,a){
    kii = rep(0,nl+2*a)
    val = (a+1)^h-c(a:0,1:a)^h
    kii[ii:(ii+2*a)] = val
    kii = kii[(a+1):(nl+a)]
    kii = kii/sum(kii)
    kii
  }
  out = matrix(sapply(1:nl, triAng,a=a),ncol=nl,nrow=nl,byrow=F) ### colSums = 1
  return(out)  
}


## EM algorithm for estimating B-spline coefficent for the smoothed Poisson
emPoisBsplMB <- function(x.list,nbreaks=50,k=4,ncores=1,err.max=0.000001,iter.max=200,t=NULL,tau.range = NULL, thr.pi=NULL,para.kern=c(1,1)){  ### NPBin model. We remove the basis that only covers bin at the boundary
    xmax.margin = 5
    err.max = 1e-6
    t = 1
	## 1 to 100 list names
    d = as.numeric(names(x.list))
	## find the largest number in each list: the largest contact counts given a distance
    xmax = sapply(x.list,max)
	## find the minimum and maximum number from xmax: the largest contact count and smallest contact count over all distances
    rangemax = range(xmax)
  if(is.null(t)){
    if(is.null(tau.range)){
	## for the bin with the smallest range, the effective range of lambda^t is larger than 4, so roughly 3/4 breaks are useful. Note that the breaks between 0 and 1 are less useful than the others.
    tau2=4  
	## comparing with the bin with the largest range, the bin with the smallest range has at least (1/3)*(3/4) proportion of effective breaks.
    tau1=3  
    }else{
    tau2 = tau.range[2]
    tau1 = tau.range[1]
    }
    tmin = log(tau2)/log(rangemax[1])
    tmax = log(tau1)/(log(rangemax[2])-log(rangemax[1]))
    ## Note that it is possible that tmin>tmax, in either case, their average could be a reasonable choice 
    t = min((tmin+tmax)/2,1)
  }
	## The largest contact counts given a distance to the power t
    xmaxt = xmax^t
	## From 0 to largest contact counts divided by 100
    breaks = seq(0,rangemax[2]^t,length.out=nbreaks)
	## create the break points and coefficients of the B spline functions
    bspl <- bsplfun(range(breaks),breaks=breaks,k=k,ncores=ncores) 
	## number of break points
    nb = length(bspl)
	## number of different distances
    nl = length(x.list)
    ## number of contacts at each distance
    nvec = sapply(x.list,length)
	## number of contacts overall
    n = sum(nvec)
	## 1 to n vector
    xall = 0:(rangemax[2]+xmax.margin)
	## for each largest count at each distance, we add the margin to it and get the list of 1 to the largest count + margin for each distance
    xall.list = lapply(xmax, function(mm) 0:(mm+xmax.margin))
	## get the design matrix
    dmtxall = getDesignMtx(bspl,xall,t,ncores)
	## design matrix for each distances
    dmtx.list = lapply(x.list, function(x) dmtxall[x+1,])
	## design matrix for largest count at each distance + margin
    dmtxall.list = lapply(xall.list, function(x) dmtxall[x+1,])
	## 0 to n divided by 101
    brk.init = seq(0,rangemax[2]^t,length.out=nb+1)
	## B-spline coefficients initialization
    pi.init.list = lapply(1:nl, function(ii) hist(x.list[[ii]]^t,breaks=brk.init,plot=F)$counts/nvec[ii])
	## triangular kernel on 1:nl with fixed range. The "bandwidth" h could be fixed at 1.
    kernl = fTriAng(nl,a=para.kern[1],h=para.kern[2])

	## initializing error to infinite
    err <- Inf
	## initializing number of iter
    iter <- 0
    ll.init <- rep(-Inf,nl)
    ll.all <- list()
    ll.all[[1]] <- ll.init
    err.all <- err
	## EM algorithm to update B-spline coefficients until difference between the new and old likelihood converges and the iteration less than a given number of times
    while(err>err.max&iter<iter.max){
      iterB <- function(ii){
		## design matrix for each distance
        dmxii = dmtx.list[[ii]]
		## length of design matrix for each distance
        nii = length(x.list[[ii]])
		## initializing B-spline coefficients for each distance
        piii = pi.init.list[[ii]]
		## find the number of "effective basis" within the range of the data in this bin
        nbii = which(breaks>=xmaxt[ii])[1]+k-4  
		## matrix of number of bin * number of effective basis
        dmxii = dmxii[,1:nbii]
        piii = piii[1:nbii]
        dtot = dmxii%*%matrix(piii,ncol=1)
        post = dmxii*matrix(piii,nvec[[ii]],nbii,byrow=TRUE)/matrix(dtot,nvec[[ii]],nbii,byrow=FALSE)        
        pi = colSums(post)
		## raw estimate of B-spline coefficients
        pi = pi/sum(pi)
		## new variable triangular kernel to smooth the spline weights B-spline coefficients over the range of lambda.
        kernpi = vTriAng(pi,wthr=log(nii)/nii) 
        pi = kernpi%*%pi  
        out = list(post=matrix(0,nvec[[ii]],nb),pi=rep(0,nb))
        out[['post']][,1:nbii] = post
        out[['pi']][1:nbii] = pi    
        out
      }
    ## apply iterB to each distance
    out.iter = mclapply(1:nl, iterB, mc.cores=ncores, mc.preschedule=F)
    post.list = lapply(1:nl, function(ii) out.iter[[ii]][['post']])
	## raw bins 
    pi.mtx = as.matrix(as.data.frame(lapply(1:nl, function(ii) out.iter[[ii]][['pi']])))
	## use fixed bandwidth tri-angular kernel to smooth over bins.
    pi.mtx = pi.mtx%*%kernl  
	## convert pi matrix to pi list according to distances
    pi.list = lapply(1:nl, function(ii) pi.mtx[,ii])     
	## get etsimate of likelihood
    ll = sapply(1:nl, function(ii) sum(apply(dmtx.list[[ii]],1, function(jj) log(sum(jj*pi.list[[ii]])))))
	## get difference between new likelihood and old likelihood
    err.vec <- ll-ll.init
    err <- max(err.vec, err.vec/abs(ll))
    ll.all[[iter+1]] <- ll
    err.all <- c(err.all,err)
	## get the coverged B-spline coefficients list
    pi.init.list <- pi.list
    ll.init <- ll
    iter <- iter+1
    }
	
	## calculate the posterior mean and variances of lambda.
    d1.list <- lapply(x.list,function(xx) dmtxall[xx+2,]) 
    d2.list <- lapply(x.list,function(xx) dmtxall[xx+3,]) 
	## obtain posterior lumbda mean for each distance and each bin pair
	lambda.post.mean.list = mclapply(1:nl, function(ii) getMomt(pi.init.list[[ii]],dmtx.list[[ii]],d1.list[[ii]],x.list[[ii]],m=1),mc.cores=ncores,mc.preschedule=F)
	## obtain posterior lumbda variance for each distance and each bin pair
	lambda.post.var.list  = mclapply(1:nl, function(ii) getMomt(pi.init.list[[ii]],dmtx.list[[ii]],d2.list[[ii]],x.list[[ii]],m=2)-lambda.post.mean.list[[ii]]^2,mc.cores=ncores,mc.preschedule=F)
	f.list    = mclapply(1:nl, function(ii) as.numeric(dmtx.list[[ii]]%*%pi.init.list[[ii]]),   mc.cores=ncores,mc.preschedule=F)
	fall.list = mclapply(1:nl, function(ii) as.numeric(dmtxall.list[[ii]]%*%pi.init.list[[ii]]),mc.cores=ncores,mc.preschedule=F) ##     fall.list=mclapply(1:nl, function(ii) as.numeric(dmtxall%*%pi.init.list[[ii]]),mc.cores=ncores,mc.preschedule=F)
    out <- list(lambda.mean=lambda.post.mean.list, lambda.var=lambda.post.var.list, dist=d, x=x.list, f=f.list, xall=xall.list, fall=fall.list, pi=pi.init.list, post=post.list, bspl=bspl, dmtx=dmtx.list, dmtxall=dmtxall,t=t,kern.bin=kernl,ll.all=ll.all,err.all=err.all, convergence=list(err=err,err.max=err.max,converged=(err<err.max),niter=iter,ll=ll.init),controls=list(k=k,nbasis=nb,breaks=breaks))
    return(out)
}


dzinb <- function(x,size,prob,p0){
    dnp = dnbinom(x,size=size,prob=prob,log=F)*(1-p0)
    ## for dnbinom, mu = size*(1-prob)/prob 
    dnp[x==0] = p0+dnbinom(0,size=size,prob=prob,log=F)*(1-p0)
    return(dnp)
}

logit <- function(p) log(p/(1-p))

expit <- function(x) 1/(1+exp(-x))


## estimate the null for one single bin using zero-inflated-negative-binomial-distribution
estNull.single <- function(x,f,xall,fall,init=NULL,iter.max=200,err.max=1e-6,q=0.5,pnull = 0.5){ 
  qnull = round(quantile(x[x>0],probs=pnull))
  if(is.null(init))
  {
    xt = x[x<=qnull]  
    m0 = mean(xt==0)
    m1 <- mean(xt)
	## initial values of the log- shape and scale paramters of gamma distn
    init <- c(log(m1),0,logit(m0)) 
  }
  nx  = length(xall) #qnull+2 #
  ell = function(ipt)
  {
    r  = exp(ipt[1])
    p  = expit(ipt[2])
    p0 = expit(ipt[3])
    f0all = dzinb(xall,size=r,prob=p,p0=p0)
    v0=f0all/fall
    vec = (abs(v0[2:(nx-1)]-v0[1:(nx-2)])^q+abs(v0[3:nx]-v0[2:(nx-1)])^q)*fall[2:(nx-1)] 
    out = mean(vec)
    out
  }
  cio = log(mean(x>0))-log(mean(x)) 
  optout <- constrOptim(theta = init,f=ell,ui=matrix(c(-1,1,0),nrow=1),ci=cio,method = "Nelder-Mead")      ### use constrained optimization
  ### for ZINB, the mean of the null model m01 = (1-p0)*r*(1-p)/p should be smaller than the mean of the data, which implies that -log(r)+logit(p)>= -log(mean(x)) + log(mean(x>0)), as an upper bound of p0 should be mean(x==0)
  coef.opt = c(exp(optout[['par']][1]),expit(optout[['par']][2]),expit(optout[['par']][3]))
  f0 <- dzinb(x,size=coef.opt[1],prob=coef.opt[2],p0=coef.opt[3])
  f0all <- dzinb(xall,size=coef.opt[1],prob=coef.opt[2],p0=coef.opt[3])
  idnull = (x>0)&(x<=qnull)
  pi0 = min(1,1/quantile(f0[idnull]/f[idnull],probs=0.98))
  coef <- c(r=coef.opt[1],p=coef.opt[2],p0=coef.opt[3],pi0=pi0)
  out <- list(coef.null=coef)
  out[['pi0']] <- pi0
  out[['f0']] <- f0
  out[['locfdr']] <- pi0*f0/f
  idq = (x<=median(x))
  out[['locfdr']][idq] <- 1  #### set locfdr to 1 if the count is less than the median, (or first quartile (or another predefined quantile)
  out[['locfdr']][out[['locfdr']]>1] <- 1
  out[['f0all']] <- f0all
  out[['convergence.null']] <- list(opt.out=optout)
  out[['out.mtx']] <- data.frame(x=x,f=f,f0=f0,locfdr=out[['locfdr']])#,FDR=out[['FDR']])
  return(out)
}


## estimate the null based on the loss function.
estNullMB.loss <- function(mod, ncores=ncores,init=NULL,iter.max=200,err.max=1e-6,q=0.5,pnull = 0.2){
  x = mod[['x']]
  f = mod[['f']]
  nl = length(x)
  d = as.numeric(names(x))
  fall = mod[['fall']]
  xall = mod[['xall']]
  ## estimate the null for one single bin using zero-inflated-negative-binomial-distribution
  output = mclapply(1:nl, function(ii) estNull.single(x=x[[ii]],f=f[[ii]],xall=xall[[ii]],fall=fall[[ii]],init=init,iter.max=iter.max,err.max=err.max,q=q,pnull=pnull),mc.cores=ncores,mc.preschedule=F)
  out = mod
  out[['coef.null']] = lapply(output, function(ii) ii[['coef.null']]) 
  out[['f0']] 		 = lapply(output, function(ii) ii[['f0']])
  out[['f0all']] 	 = lapply(output, function(ii) ii[['f0all']])
  out[['locfdr']] 	 = lapply(output, function(ii) ii[['locfdr']]) 
  out[['out.mtx.list']] = lapply(output, function(ii) ii[['out.mtx']]) 
  return(out)
}


## loci-pair false discovery rate
locfdr2FDR <- function(locfdr,pthr=0.5){
n = length(locfdr)
out = rep(1,n)
id = which(locfdr<=pthr)
locfdr.id = locfdr[id]
out.id = sapply(1:length(id), function(ii) mean(locfdr.id[locfdr.id<=locfdr.id[ii]]))
out[id] = out.id
return(out)
}


rank2nhit <- function(r,id){
out <- sapply(r, function(y) sum((r<=y)&id))
return(out)
}


### a wrapper function for EBHiC using Fit-Hi-C output file as the input.
ebhicWrapper <- function(dt, nbreaks=101,ncores=1,iter.max=100,init=NULL,error.max=1e-6,para.kern=c(1,1),q=2,pnull=0.9){
  ### the wrapper of the empirical bayes model of HiC, input the chr data, and output the sameobject with additional columns representing our output.
  
  ## obtian the distance between different pairs of bins
  distance = abs((dt[,"fragmentMid2"]-dt[,"fragmentMid1"]))
  ## combine poitions of bins, number of contacts, p_value and q_value
  dt = dt[,c('fragmentMid1','fragmentMid2', 'contactCount', 'p_value', 'q_value')]
  ## rename the columns in dt
  colnames(dt) = c('fragmentMid1','fragmentMid2', 'contactCount', 'pval_fithic', 'qval_fithic')
  ## include distances into dt
  dt$distance = distance
  ## grouping infomation in dt according to distance
  dt.list = split(dt, f=distance)
  rm(dt)
  ## obtain number of different distances
  nl = length(dt.list)
  d  = names(dt.list)
  ## grouping contact counts according to distance
  x.list = lapply(dt.list, function(xx) xx[,'contactCount'])
  #  k = 4 # order of spline
  #  nbk = 101
  #      print(182)
  
  modtemp.mb = emPoisBsplMB(x.list, nbreaks=nbreaks, k=4, ncores=ncores, iter.max=iter.max, para.kern=para.kern)## amount of smoothing inscrease as one or two of the two parameters increase.
  #	  print(185)
  rm(x.list)
  ## estimate the null based on the loss function 
  nulltemp.mb = estNullMB.loss(modtemp.mb,ncores=ncores,init=init,iter.max=iter.max,err.max=err.max,q=q,pnull = pnull)
  #	  print(188)
  rm(modtemp.mb)
  ## get the mean, variance and locfdr of lambda
  lambda.mean = nulltemp.mb[['lambda.mean']]
  lambda.var = nulltemp.mb[['lambda.var']]
  locfdr = nulltemp.mb[['locfdr']]
  #      print(191)
  
  out.list = lapply(1:nl, function(ii){
    out.ii 		  = dt.list[[ii]]
    out.ii$locfdr 	  = locfdr[[ii]]
    out.ii$mean.post = lambda.mean[[ii]]
    out.ii$var.post  = lambda.var[[ii]]
    out.ii})
  #      print(199)
  
  output = do.call('rbind',out.list)
  output$qvalue = locfdr2FDR(output$locfdr,pthr=0.5)
  #      print(203)
  
  return(list(mod=nulltemp.mb,output=output))
}

