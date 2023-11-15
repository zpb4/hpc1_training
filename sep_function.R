#Implementation of SEP from Hutson et al. (2019)
library(gsl)
library(fGarch)
library(stats)

#helper functions
kfun<-function(beta,alpha){
  kinv<-gamma(1+((beta+1)/2))*2^(1+0.5*(1+beta))/(4*alpha*(1-alpha))
  k<-1/kinv
  return(k)
}

zfun<-function(x,theta,sigma){
  z<-(x-theta)/sigma
  return(z)
}

#SEP pdf function
sep_pdf<-function(x,theta,sigma,beta,alpha){
  k<-kfun(beta,alpha)
  z<-zfun(x,theta,sigma)
  y=(k/sigma)*exp(-0.5*(abs(z)+(2*alpha-1)*z)^(2/(1+beta)))
  return(y)
}

#SEP cdf function
sep_cdf<-function(x,theta,sigma,beta,alpha){
  z<-zfun(x,theta,sigma)
  y<-c()
  for(i in 1:length(z)){
    if(z[i]<=0){
      gam_x<-2^((2/(beta+1))-1) * ((alpha-1)*z[i])^(2/(beta+1))
      gam_a<-((beta+1)/2)
      num<-alpha*gamma_inc(gam_a,gam_x)
      den<-gamma((beta+1)/2)
      y[i]<-num/den
    }
    if(z[i]>0){
      gam_x<-2^((2/(beta+1))-1)*(1/(alpha*z[i]))^(-2/(beta+1))
      gam_a<-((beta+1)/2)
      num<-(1-alpha)*gamma_inc(gam_a,gam_x,)
      den<-gamma((beta+1)/2)
      y[i]<-1-(num/den)
    }
  }
  return(y)
}

#SEP quantile function
sep_quant<-function(u,theta,sigma,beta,alpha){
  q<-c()
  for(i in 1:length(u)){
    if(u[i]<=alpha){
      num<-(2^(0.5-(1/(beta+1))) * sqrt(qgamma((1-(u[i]/alpha)),((beta+1)/2))))^(beta+1)
      den<-1-alpha
      q[i]<-theta - sigma*(num/den)
    }
    if(u[i]>alpha){
      num<-(2^((1/(beta+1))-0.5) / sqrt(qgamma((1-((1-u[i])/(1-alpha))),((beta+1)/2))))^(-(beta+1))
      den<-alpha
      q[i]<-theta + sigma*(num/den)
    }
  }
  return(q)
}

sep_quant2<-function(u,theta,sigma,beta,alpha){
  q<-c()
  for(i in 1:length(u)){
    if(u[i]<=alpha){
      num<-(2^(0.5-(1/(beta+1))) * sqrt(gamma((beta+1)/2)/gamma_inc((beta+1)/2,u[i]/alpha)))^(beta+1)
      den<-1-alpha
      q[i]<-theta - sigma*(num/den)
    }
    if(u[i]>alpha){
      num<-(2^((1/(beta+1))-0.5) / sqrt(gamma((beta+1)/2)/gamma_inc((beta+1)/2,(1-u[i])/(1-alpha))))^(-(beta+1))
      den<-alpha
      q[i]<-theta + sigma*(num/den)
    }
  }
  return(q)
}

#SEP sampling function
sep_samp<-function(n,theta,sigma,beta,alpha){
  u<-runif(n)
  samps<-sep_quant(u,theta,sigma,beta,alpha)
  return(samps)
}


#MLE functions

#The functions within ----- are from Hutson et al., but did not work for some reason
#----------------------------------------------------------------------------
#log likelihood function per observation
#sep_ll<-function(x,theta,sigma,beta,alpha){
  #ll<-(-1)*((abs(x-theta)+(2*alpha-1)*(x-theta))^(2/(beta+1)))/(2*sigma) + log(alpha*(1-alpha)) - 
    #(1+0.5*(beta+1))*log(2) - log(gamma((beta+1)/2+1)) - log(sigma)
  #return(ll)
#}

#log likelihood function for MLE with analytic sigma estimation
#sep_mle<-function(x,pars){
  #mu=pars[1]
  #sigma=pars[2]
  #beta=pars[3]
  #alpha=pars[4]
  #ll=sum(sep_ll(x,mu,sigma,beta,alpha))
  #if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  #return(ll)
#}

#log likelihood function for MLE with analytic sigma estimation
#sep_mle_sigest<-function(x,pars){
  #mu=pars[1]
  #beta=pars[2]
  #alpha=pars[3]
  #sigma=sep_sigma_est(x,mu,beta,alpha)
  #ll=sum(sep_ll(x,mu,sigma,beta,alpha))
  #if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  #return(ll)
#}
#-----------------------------------------------------------------------------

#log likelihood function for MLE using canonical form with SEP pdf
sep_mle<-function(x,pars){
  theta=pars[1]
  sigma=pars[2]
  beta=pars[3]
  alpha=pars[4]
  ll=sum(log(sep_pdf(x,theta,sigma,beta,alpha)))
  if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  return(ll)
}

#analytic sigma estimate
sep_sigma_est<-function(x,theta,beta,alpha){
  n<-length(x)
  num<-sum((abs(x-theta)+(2*alpha-1)*(x-theta))^(2/(beta+1)))
  den<-n*(beta+1)
  sigma<-(num/den)^((beta+1)/2)
  return(sigma)
}

#log likelihood function for MLE using canonical form with SEP pdf and analytic sigma estimate
sep_mle_sigest<-function(x,pars){
  theta=pars[1]
  beta=pars[2]
  alpha=pars[3]
  sigma=sep_sigma_est(x,theta,beta,alpha)
  ll=sum(log(sep_pdf(x,theta,sigma,beta,alpha)))
  if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  return(ll)
}

######################################END##################################