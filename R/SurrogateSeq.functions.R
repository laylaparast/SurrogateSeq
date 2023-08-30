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
	return(list("delta.e" = delta.e, "sd.e" = sqrt(var(mu.s1)/length(sone)+var(mu.s0)/length(szero)), "sd.closed" = sd.e.closed, "delta.e.z" = delta.e.z, "delta.e.p" = delta.e.p))
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

