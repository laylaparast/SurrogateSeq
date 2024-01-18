delta.e.estimate = function(sone=NULL,szero=NULL, szerop, yzerop, extrapolate = TRUE,mat = NULL, n1=NULL, n0=NULL) {
	if(!is.null(mat)){
		sone = mat[1:n1]
		szero = mat[(n1+1):(n1+n0)]
	}
	h.paper4 = bw.nrd(szerop)*(length(szerop)^(-0.11))
	mu.s0 = sapply(szero,pred.smooth.2,kernel.use=szerop, bw=h.paper4, outcome=yzerop)
  	if(sum(is.na(mu.s0))>0 & extrapolate){
  		print(paste("Note: ", sum(is.na(mu.s0)), " values extrapolated."))
    	c.mat = cbind(szero, mu.s0)
    	for(o in 1:length(mu.s0)) {
    		if(is.na(mu.s0[o])){
    			distance = abs(szero - szero[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			mu.s0[o] = new.est[1]   #in case there are multiple matches
    		}
    	}
		}
	mu.s1 = sapply(sone,pred.smooth.2,kernel.use=szerop, bw=h.paper4, outcome=yzerop)
  	if(sum(is.na(mu.s1))>0 & extrapolate){
  		print(paste("Note: ", sum(is.na(mu.s1)), " values extrapolated."))
    	c.mat = cbind(sone, mu.s1)
    	for(o in 1:length(mu.s1)) {
    		if(is.na(mu.s1[o])){
    			distance = abs(sone - sone[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			mu.s1[o] = new.est[1]   #in case there are multiple matches
    		}
    	}
		}
	delta.e = mean(mu.s1) - mean(mu.s0)	
	n.b1 = length(sone); n.b0 = length(szero); n.b = n.b1+n.b0
	var.first.term = (n.b/(n.b1^2))*(sum(mu.s1^2) - (1/n.b1)*(sum(mu.s1)^2))
	var.second.term = (n.b/(n.b0^2))*(sum(mu.s0^2) - (1/n.b0)*(sum(mu.s0)^2)) 
	var.closed = (var.first.term + var.second.term)/n.b
	#note that what this return is sigma_e/sqrt(n_b), then Z = delta.e/THIS THING
	sd.e.closed = sqrt(var.closed)
	delta.e.z = delta.e/sd.e.closed
	delta.e.p=2*(1- pnorm(abs(delta.e.z)))
	return(list("delta.e" = delta.e, "sd.closed" = sd.e.closed, "delta.e.z" = delta.e.z, "delta.e.p" = delta.e.p))
}

pred.smooth.2 <-function(kernel.use,kernel.apply, bw,outcome) { 	
    return(sum(Kern.FUN(kernel.use,kernel.apply,bw=bw)*outcome)/sum(Kern.FUN(kernel.use,kernel.apply,bw=bw)))
 
  }

Kern.FUN <- function(zz,zi,bw) 
  { 
    out = (VTM(zz,length(zi))- zi)/bw
	dnorm(out)/bw
           
  }
  
 VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }

gs.boundaries = function(szerop, sonep, yzerop, n.stg, B.norm=1e6, alpha=0.05, pp=0.4, inf.fraction = (1:n.stg)/n.stg, plot=FALSE){
	bdr.naive = rep(qnorm(1-0.05/2),n.stg)
	bdr.bonf = rep(qnorm(1-0.05/(2*n.stg)),n.stg)
	
	# Compute covariance matrix
	cov.stuffA4=cov.surr.gs(s0.4.est=s0.studya, s1.4.est=s1.studya, 	sa.0=s0.studya, ya.0=StudyA.aids$y0, nb.0=nrow(StudyB.aids$s0), nb.1=nrow(StudyB.aids$s1), full.matrix=TRUE, naive=TRUE)

	# Simulate MVN null paths
	paths.norm4=mvrnorm(n=B.norm, mu=rep(0,n.stg), Sigma=cov.stuffA4$cov.stand.samp)  

	# Compute boundaries
	bdr.Poc=bdr.gs.mc.gen(alpha=alpha, mc.paths=paths.norm4, w.vec=rep.int(1, times=n.stg))  

    bdr.OF=bdr.gs.mc.gen(alpha=alpha, mc.paths=paths.norm4, w.vec=sqrt(1/inf.fraction))  
	bdr.WT=bdr.gs.mc.gen(alpha=alpha, mc.paths=paths.norm4, w.vec=inf.fraction^(pp-.5))  
    
    if(plot){
    	mydata = as.data.frame(cbind(c(1:n.stg), rep(0,n.stg)))
    	names(mydata) = c("Time", "Stat")
	ymax = max(bdr.naive[1],bdr.bonf[1], bdr.OF$bndry.vec, bdr.Poc$bndry.vec, bdr.WT$bndry.vec)
	ymin = -ymax

    p1 =	ggplot(mydata, aes(Time, Stat)) + geom_point(size=0, color="white") + geom_hline(yintercept=bdr.naive[1], linetype="dashed", color = "red") + geom_hline(yintercept=-bdr.naive[1], linetype="dashed", color = "red") + geom_hline(yintercept=bdr.bonf[1], linetype="dashed", color = "blue") + geom_hline(yintercept=-1*bdr.bonf[1], linetype="dashed", color = "blue") + stat_function(fun = function(xx) bdr.OF$cons*(xx/n.stg)^(-1/2),linetype="dashed", color = "purple") + stat_function(fun = function(xx) -bdr.OF$cons*(xx/n.stg)^(-1/2),linetype="dashed", color = "purple") + stat_function(fun = function(xx) bdr.Poc$cons,linetype="dashed", color = "black")+ stat_function(fun = function(xx) -bdr.Poc$cons,linetype="dashed", color = "black") + stat_function(fun = function(xx) bdr.WT$cons*(xx/n.stg)^(.4-1/2),linetype="dashed", color = "grey") + stat_function(fun = function(xx) -bdr.WT$cons*(xx/n.stg)^(.4-1/2),linetype="dashed", color = "grey") + ylim(ymin,ymax) + labs(title = "Boundaries", x = "Time/Stages", y = expression(paste("Test statistic, ",Z[E]))) + annotate(geom = 'text', label = "Naive = red; Bonferroni = blue; O'Brien-Fleming = purple; Pocock = black; Wang-Tsiatis = grey", x = -Inf, y = Inf, hjust = 0, vjust = 1, size=2)
    



    }
    if(!plot) {p1 = NULL}
	return(list("Naive" = bdr.naive, "Bonf" = bdr.bonf, "Pocock" = bdr.Poc$bndry.vec, "OBrien_Fleming" = bdr.OF$bndry.vec, "Wang_Tsiatis" = bdr.WT$bndry.vec, plot = p1))

}

cov.surr.gs = function(s0.4.est, s1.4.est, sa.0, ya.0, nb.0, nb.1, full.matrix=TRUE, naive=FALSE) {
 # Return a vector of the variances and, if full.matrix=TRUE, the standardized covariance matrix for the group sequential statistic for Layla's surrogate marker problem.  Both of these are based on Jay's delta method approximations.  See my notes for formulas.
  #
  # INPUTS
  # s0.4.est, s1.4.est = s0 and s1 matrices for estimating means and covariances of S_0, S_1 in the study B data.  For designing tests (e.g., finding boundaries) these may come from Study A data, but for analyzing tests these may come from study B data. No. of columns is the number of stages, number of rows may differ from rows in sa.0
  # sa.0, ya.0 = study A data
  # nb.0, nb.1 = study B sample sizes, i.e., the number of surrogates observed at each stage for the test statistic we're computing covariance for
  # naive = set to TRUE to compute covariance for "cumulative" test statistic, FALSE for "naive" statistic that only uses study B data from timepoint J at the J-th analysis.
  #
  # OUTPUTS
  # var.vec.del, cov.stand.del = a vector of the variances and, if full.matrix=TRUE (otherwise will be NULL), the covariance matrix of the standardized test statistic computed by the delta method.
  # var.vec.samp, cov.stand.samp = the same but computed by the sample mean and covariance of s0.4.est and s1.4.est.
  #
  #  J.B. 23.May.23; last modified 3.Jul.23
  
   
  n.stg=ncol(sa.0)
  na.0=length(ya.0)

  theta.0=apply(s0.4.est,2,mean)
  theta.1=apply(s1.4.est,2,mean)
  covS.0=cov(s0.4.est)
  covS.1=cov(s1.4.est)
  
  d.0=vector(mode="numeric", length=n.stg)
  d.1=vector(mode="numeric", length=n.stg)
  
  h=(1/na.0)^(1/3)
  
  # make d's
  for (j in 1:n.stg) {
    d.0[j]=(sum((sa.0[,j]-theta.0[j])*dnorm((sa.0[,j]-theta.0[j])/h)*ya.0) * sum(dnorm((sa.0[,j]-theta.0[j])/h)) - sum((sa.0[,j]-theta.0[j])*dnorm((sa.0[,j]-theta.0[j])/h)) * sum(dnorm((sa.0[,j]-theta.0[j])/h)*ya.0)) / (h*sum(dnorm((sa.0[,j]-theta.0[j])/h)))^2
    d.1[j]=(sum((sa.0[,j]-theta.1[j])*dnorm((sa.0[,j]-theta.1[j])/h)*ya.0) * sum(dnorm((sa.0[,j]-theta.1[j])/h)) - sum((sa.0[,j]-theta.1[j])*dnorm((sa.0[,j]-theta.1[j])/h)) * sum(dnorm((sa.0[,j]-theta.1[j])/h)*ya.0)) / (h*sum(dnorm((sa.0[,j]-theta.1[j])/h)))^2
  }

  # Make variance vector
  var.vec.del=vector(mode="numeric", length=n.stg)
  var.vec.samp=vector(mode="numeric", length=n.stg)
  
  mu.hat0=matrix(nrow=nrow(s0.4.est), ncol=n.stg)
  mu.hat1=matrix(nrow=nrow(s1.4.est), ncol=n.stg) # fill these with the kernel estimators of s0.4.est, s1.4.est
  for (j in 1:n.stg) {
    #cat("j=",j,"\n",sep='')
    mu.hat0[,j]=sapply(X=s0.4.est[,j], FUN=kern.estJ, band.h=h, sa.vec=sa.0[,j], ya.vec=ya.0)
    mu.hat1[,j]=sapply(X=s1.4.est[,j], FUN=kern.estJ, band.h=h, sa.vec=sa.0[,j], ya.vec=ya.0)
  }
    
  cov.muS.0=cov(mu.hat0)
  cov.muS.1=cov(mu.hat1)
  
  if (!naive) { # "cumulative" test statistic
    for (j in 1:n.stg) {
      var.vec.del[j]=(1/nb.0)*(t(d.0[1:j]) %*% covS.0[1:j,1:j] %*% d.0[1:j]) + (1/nb.1)*(t(d.1[1:j]) %*% covS.1[1:j,1:j] %*% d.1[1:j])
      var.vec.samp[j]=(1/nb.0)*sum(cov.muS.0[1:j,1:j]) + (1/nb.1)*sum(cov.muS.1[1:j,1:j])
    } 
  } else { # the naive test statistic
    var.vec.del=d.0^2*diag(covS.0)/nb.0 + d.1^2*diag(covS.1)/nb.1
    var.vec.samp=diag(cov.muS.0)/nb.0 + diag(cov.muS.1)/nb.1
  }
  
  # Make covariance matrix?
  if (!full.matrix) {  # don't compute full matrix
    cov.stand.del=NULL
    cov.stand.samp=NULL
  } else {  # compute full matrix
    # un-standardized cov matrices
    cov.unstand.del = diag(var.vec.del) 
    cov.unstand.samp = diag(var.vec.samp)
    if (!naive) { # "cumulative" test statistic
      for (j in 2:n.stg) {
        for (jj in 1:(j-1)) {
          cov.unstand.del[j,jj]=(1/nb.0)*(t(d.0[1:j]) %*% covS.0[1:j,1:jj] %*% d.0[1:jj]) + (1/nb.1)*(t(d.1[1:j]) %*% covS.1[1:j,1:jj] %*% d.1[1:jj])
          cov.unstand.del[jj,j]=cov.unstand.del[j,jj]
          cov.unstand.samp[j,jj]=(1/nb.0)*sum(cov.muS.0[1:j,1:jj]) + (1/nb.1)*sum(cov.muS.1[1:j,1:jj])
          cov.unstand.samp[jj,j]=cov.unstand.samp[j,jj]
        }
      }
    } else { # the naive test statistic
      for (j in 2:n.stg) {
        for (jj in 1:(j-1)) {
          cov.unstand.del[j,jj]=(1/nb.0)*d.0[j]*covS.0[j,jj]*d.0[jj] + (1/nb.1)*d.1[j]*covS.1[j,jj]*d.1[jj]
          cov.unstand.del[jj,j]=cov.unstand.del[j,jj]
        }
      }
      cov.unstand.samp=cov.muS.0/nb.0 + cov.muS.1/nb.1
    }
    
    # standardize
    cov.stand.del=cov2cor(cov.unstand.del)
    cov.stand.samp=cov2cor(cov.unstand.samp)
  }

  return(list("var.vec.del"=var.vec.del, "cov.stand.del"=cov.stand.del, "var.vec.samp"=var.vec.samp, "cov.stand.samp"=cov.stand.samp))
}

kern.estJ=function(sb.arg, band.h, sa.vec, ya.vec) {
  # Takes in scalar sb.arg and computes kernel density estimate using other args.
  #
  # J.B. 29.Jun.23; last modified 3.Jul.23.
  
  #cat("length(sa.vec)=", length(sa.vec), ", length(sb.arg)=", length(sb.arg),"\n", sep='')
  w=dnorm((sa.vec-sb.arg)/band.h) # weights
  
  if (sum(w)>0) { # usual
    #cat("length(ya.vec)=", length(ya.vec),", length(w)=", length(w), "\n", sep='')
    return(weighted.mean(x=ya.vec, w=w))
  } else { # weights have underflowed so extrapolate
    return(ya.vec[which.min(abs(sa.vec-sb.arg))]) # the y going with the s closest to sb.arg
  }
}

bdr.gs.mc.gen = function(alpha=.05, mc.paths, w.vec) {
 # Return the boundaries for any group sequential test of the null vs. 2-sided alternative whose boundaries take the form of a single constant times a known weight vector, which is w.vec.  These include Pocock (w.vec=(1,1,..)), O'Brien-Fleming (w.vec=(sqrt(n.stg/1), sqrt(n.stg/2), ..., 1)), etc.  It does this by returning quantiles of the sample paths of the (null) test statistic paths in mc.paths.   
  #
  # INPUTS
  # alpha = desired rejection probability of the test
  # mc.paths = matrix of sample paths, each row being a sample path, no. of columns is number of stages
  #
  # OUTPUTS
  # cons = the constant in the boundary vector cons*w.vec
  # bndry.vec = the boundary vector cons*w.vec
  #
  #  J.B. 4.Jul.23; last modified 6.Nov.23
  
  d=dim(mc.paths)
  n.stg=d[2]  # no. cols
  B=d[1]  # no. rows
  
  mc.paths.adj=abs(mc.paths) %*% diag(1/w.vec)  # j-th test stat gets multiplied by 1/w.vec[j]
  max.stat=apply(mc.paths.adj,1, function(x) max(x))  # max of each (adjusted) path
  cons=quantile(x=max.stat, probs=1-alpha)
  
  return(list("cons"=cons, "bndry.vec"= cons*w.vec))
}

gs.boundaries.fut = function(szerop, sonep, yzerop, n.stg, B.norm=1e6, alpha=0.05, pp=0.4, inf.fraction = (1:n.stg)/n.stg, j.star=1, alpha0=(j.star/n.stg)*alpha, plot=FALSE){
	bdr.naive = rep(qnorm(1-0.05/2),n.stg)
	bdr.bonf = rep(qnorm(1-0.05/(2*n.stg)),n.stg)
	
	# Compute covariance matrix
	cov.stuffA4=cov.surr.gs(s0.4.est=s0.studya, s1.4.est=s1.studya, 	sa.0=s0.studya, ya.0=StudyA.aids$y0, nb.0=nrow(StudyB.aids$s0), nb.1=nrow(StudyB.aids$s1), full.matrix=TRUE, naive=TRUE)

	# Simulate MVN null paths
	paths.norm4=mvrnorm(n=B.norm, mu=rep(0,n.stg), Sigma=cov.stuffA4$cov.stand.samp)  

	# Compute boundaries
	bdr.Poc=bdr.gs.mc.fut(c1=NULL, c2=NULL,pp=.5,j.star=j.star,alpha=alpha,alpha0=.05*.75, mc.paths=paths.norm4, inf.fraction=inf.fraction, n.stg=n.stg)  

    bdr.OF=bdr.gs.mc.fut(c1=NULL, c2=NULL, pp=0, n.stg=n.stg, j.star=j.star, alpha=.05, alpha0=.05*.1, mc.paths=paths.norm4,inf.fraction=inf.fraction)  
	bdr.WT=bdr.gs.mc.fut(c1=NULL, c2=NULL, pp=pp, n.stg=n.stg, j.star=j.star, alpha=.05, mc.paths=paths.norm4, inf.fraction=inf.fraction)  
    
    Pocock.futility = bdr.Poc$a
     Pocock.nullrejection = bdr.Poc$b
     OBrien_Fleming.futility = bdr.OF$a
     OBrien_Fleming.nullrejection = bdr.OF$b
     Wang_Tsiatis.futility = bdr.WT$a
     Wang_Tsiatis.nullrejection = bdr.WT$b
     
    if(plot){
    	mydata = as.data.frame(cbind(c(1:n.stg), rep(0,n.stg)))
    	names(mydata) = c("Time", "Stat")
		ymax = max(bdr.naive[1],bdr.bonf[1], bdr.OF$bndry.vec, bdr.Poc$bndry.vec, bdr.WT$bndry.vec)
	ymin = -ymax

    p1 =	ggplot(mydata, aes(Time, Stat)) + geom_point(size=0, color="white") + geom_hline(yintercept=bdr.naive[1], linetype="dashed", color = "red") + geom_hline(yintercept=-bdr.naive[1], linetype="dashed", color = "red") + geom_hline(yintercept=bdr.bonf[1], linetype="dashed", color = "blue") + geom_hline(yintercept=-1*bdr.bonf[1], linetype="dashed", color = "blue") + stat_function(fun = function(xx) bdr.OF$c1*(xx/n.stg)^(-1/2),linetype="dashed", color = "purple") + stat_function(fun = function(xx) -bdr.OF$c1*(xx/n.stg)^(-1/2),linetype="dashed", color = "purple") + stat_function(fun = function(xx) ((bdr.OF$c1+bdr.OF$c2)*(xx/n.stg)^(.5)-bdr.OF$c2*(xx/n.stg)^(0-.5)),linetype="dashed", color = "purple") + stat_function(fun = function(xx) -((bdr.OF$c1+bdr.OF$c2)*(xx/n.stg)^(.5)-bdr.OF$c2*(xx/n.stg)^(0-.5)),linetype="dashed", color = "purple") + stat_function(fun = function(xx) bdr.Poc$c1,linetype="dashed", color = "black")+ stat_function(fun = function(xx) -bdr.Poc$c1,linetype="dashed", color = "black") + stat_function(fun = function(xx) ((bdr.Poc$c1+bdr.Poc$c2)*(xx/n.stg)^(.5)-bdr.Poc$c2*(xx/n.stg)^(.5-.5)),linetype="dashed", color = "black")+ stat_function(fun = function(xx) -((bdr.Poc$c1+bdr.Poc$c2)*(xx/n.stg)^(.5)-bdr.Poc$c2*(xx/n.stg)^(.5-.5)),linetype="dashed", color = "black") + stat_function(fun = function(xx) bdr.WT$c1*(xx/n.stg)^(pp-1/2),linetype="dashed", color = "grey") + stat_function(fun = function(xx) -bdr.WT$c1*(xx/n.stg)^(pp-1/2),linetype="dashed", color = "grey") + stat_function(fun = function(xx) ((bdr.WT$c1+bdr.WT$c2)*(xx/n.stg)^(.5)-bdr.WT$c2*(xx/n.stg)^(pp-.5)),linetype="dashed", color = "grey") + stat_function(fun = function(xx) -((bdr.WT$c1+bdr.WT$c2)*(xx/n.stg)^(.5)-bdr.WT$c2*(xx/n.stg)^(pp-.5)),linetype="dashed", color = "grey") + ylim(ymin,ymax) + labs(title = "Boundaries", x = "Time/Stages", y = expression(paste("Test statistic, ",Z[E]))) + annotate(geom = 'text', label = "Naive = red; Bonferroni = blue; O'Brien-Fleming = purple; Pocock = black; Wang-Tsiatis = grey", x = -Inf, y = Inf, hjust = 0, vjust = 1, size=2)

    }
    if(!plot) {p1 = NULL}
	return(list("Naive" = bdr.naive, "Bonf" = bdr.bonf, "Pocock.futility" = Pocock.futility, "Pocock.nullrejection" = Pocock.nullrejection, "OBrien_Fleming.futility" = OBrien_Fleming.futility, "OBrien_Fleming.nullrejection" = OBrien_Fleming.nullrejection, "Wang_Tsiatis.futility" = Wang_Tsiatis.futility, "Wang_Tsiatis.nullrejection" = Wang_Tsiatis.nullrejection, plot = p1))

}

bdr.gs.mc.fut = function(c1=NULL, c2=NULL, pp=.4, n.stg, j.star=1, alpha=.05, alpha0=(j.star/n.stg)*alpha, mc.paths, inf.fraction=(1:n.stg)/n.stg, N.iter.max=100, alpha.tol=.02*alpha) {
  
 # Return the boundaries for a group sequential test of the null vs. 2-sided alternative *with futility stopping* with boundaries given by the Wang-Tsiatis power family as in Jennison & Turnbull, Exp (5.3) on page 118.  It does this by returning quantiles of the sample paths of the (null) test statistic paths in mc.paths.   
  #
  # INPUTS
  # c1, c2 = these are the constants determining the outer boundary b[j] = c1 * (j/J)^{pp-1/2} and futility boundaries a[j] = (c1+c2) * (j/J)^{1/2} - c2 * (j/J)^{pp-1/2} for j >= j.star, where J is the max no of stages (AKA n.stg). If c1 is null, it is found as the upper alpha0 quantile of the max over the first j.star stages.
  # pp = parameter in the power family of boundaries. pp=.5 gives boundaries similar to Pocock, pp=0 gives similar to O'Brien-Fleming.
  # n.stg = maximum number of analyses (or "stages")
  # j.star = earliest stage at which futility stopping is allowed. Should be <= n.stg-1 (there is already "futility stopping" at the n.stg-th stage anyway).
  # alpha = desired type I error probability
  # alpha0 = the part of alpha that c1 is chosen to spend in first j.star stages
  # mc.paths = matrix of sample paths, each row being a sample path, no. of columns should be no. of stages, n.stg
  # inf.fraction = information fraction used in computing boundaries, an n.stg-long vector. Should be positive, non-decreasing, with last element 1.  Default is 1/n.stg, 2/n.stg, ..., 1.
  # N.iter.max = max no. of iterations for finding c2
  # alpha.tol = the tolerance for stopping search for c2
  #
  # OUTPUTS
  # a,b = the futility and null-rejection boundary vectors
  # prej = prob. of rejecting the null (at any stage)
  # EM, se.M
  # c1, c2 = constants used in boundaries a, b
  #
  #  J.B. 4.Jul.23; last modified 11.Jan.24
  
   
  #d=dim(mc.paths)
  #B=d[1]  # no. rows
  
  if (is.null(c1)) {
    if (j.star==1) { # 1x1 matrices cause problems so do this by hand
      paths.c1= abs(mc.paths[,1:j.star]) * (inf.fraction[1:j.star])^(-(pp-.5))
    } else { # j.star >= 2 so do regular matrix multiplication
      paths.c1=apply(abs(mc.paths[,1:j.star]) %*% diag((inf.fraction[1:j.star])^(-(pp-.5))), 1, function(x) max(x))
    }
    c1=quantile(x=paths.c1, probs=1-alpha0)
  }
  
  b=c1*(inf.fraction^(pp-.5))
  if (is.null(c2)) {  # Find of c2
    # def aux function
    p.rej.cand=function(c2.cand) {  # prob of rejecting null, using cs.cand as candidate value for c2
      a.cand=rep.int(x=0, times=n.stg)
      r=inf.fraction[j.star:n.stg]
      a.cand[j.star:n.stg]=pmax(0,(c1+c2.cand)*sqrt(r) - c2.cand*(r^(pp-.5)))
      return(op.char.gs.fut(b.vec=b, a.vec=a.cand, paths=mc.paths)$prej)
    }
    c2.L=-c1; p.L=p.rej.cand(c2.L)
    c2.R=c1/((inf.fraction[n.stg-1])^(pp-1)-1); p.R=p.rej.cand(c2.R)
    if (p.L>alpha) {
      a=rep.int(0,times=n.stg) # no futility stopping
    } else if (p.R<alpha) { # use c2=c2.R for some (conservative) futility stopping
      c2=c2.R
      a=rep.int(x=0, times=n.stg)
      r=inf.fraction[j.star:n.stg]
      a[j.star:n.stg]=pmax(0,(c1+c2)*sqrt(r) - c2*(r^(pp-.5)))
    } else {
      # p.L <= alpha and p.R >= alpha, so start bisection method
      N.iter=1
      c2=.5*(c2.L+c2.R)
      p.c2=p.rej.cand(c2)
      while ((N.iter<N.iter.max) & (abs(p.c2-alpha) > alpha.tol)) {
        N.iter=N.iter+1
        if (p.c2>alpha) c2.R=c2 else c2.L=c2
        c2=.5*(c2.L+c2.R)
        p.c2=p.rej.cand(c2)
      }
      #cat("c2 = ", c2, ", P(rej) = ", p.c2,"\n", sep="")
      a=rep.int(x=0, times=n.stg)
      r=inf.fraction[j.star:n.stg]
      a[j.star:n.stg]=pmax(0,(c1+c2)*sqrt(r) - c2*(r^(pp-.5)))
    }
  } else { # c2 was given 
    a=rep.int(x=0, times=n.stg)
    r=inf.fraction[j.star:n.stg]
    a[j.star:n.stg]=pmax(0,(c1+c2)*sqrt(r) - c2*(r^(pp-.5)))
  }
  
  # call op.char.gs.fut
  op.char=op.char.gs.fut(b.vec=b, a.vec=a, paths=mc.paths)  # OUTPUTS EM, se.M, prej
    
  return(list("EM"=op.char$EM, "se.M"=op.char$se.M, "prej"=op.char$prej, "a"=a, "b"=b, "c1"=c1, "c2"=c2))
}

op.char.gs.fut = function(b.vec,a.vec,paths) {
 # Compute the operating characteristics on the group sequential test with futility stopping statistics in paths: The expected stopping stage no., plus the probability of rejecting the null in favor of the 2-sided alternative. This is for a general GS test which uses the boundaries in bndry.vec.
  #
  # INPUTS
  # a.vec, b.vec = futility and "null-rejection" boundaries, respectively.  Both should be >= 0, and a.vec[n.stg] = b.vec[n.stg]. a.vec[j]=0 means no futility stopping at stage j.
  # paths = matrix of test statistic sample paths, each row being a sample path, no. of columns is max number of stages
  #
  # OUTPUTS
  # EM = expected stopping stage number
  # se.M = standard error of stopping time
  # prej = probability of rejecting null (at any stage)
  #
  #  J.B. 6.Jul.23; last modified 6.Jul.23
  
  d=dim(paths)
  n.stg=d[2]  # no. cols
  B=d[1]  # no. rows
  
  paths.abs=abs(paths)
  
  M.a=apply(paths.abs %*% diag(1/a.vec), 1, function(x) which(x<1)[1]) # where was first crossing of futility boundary (not necessarily pre-stopping time), or NA if no crossing
  M.b=apply(paths.abs %*% diag(1/b.vec), 1, function(x) which(x>=1)[1]) # where was first crossing of null-rejection boundary (not necessarily pre-stopping time), or NA if no crossing
  
  M.a[is.na(M.a)] = n.stg+1  # deal w/ NA's before taking mins and comparing
  M.b[is.na(M.b)] = n.stg+1
  
  M=pmin(M.a, M.b)
  EM=mean(M)
  se.M=sd(M)/sqrt(B)
  
  prej=sum(M.b<M.a)/B
  
  return(list("EM"=EM, "se.M"=se.M, "prej"=prej))
}