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