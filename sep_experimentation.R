#setwd('h:/Projects/MA SWM collaboration - Phase 2/logratio_experimentation/')
source('sep_function.R')
library(fGarch)

#Script to experiment with manually coded SEP distribution (Hutson et al., 2019) as compared to
#SGED distribution from R fGarch package

#SGED params
#mu: location parameter, 'mean'
#sigma: scale parameter, 'sd'
#nu: shape parameter, 'nu'
#xi: skew parameter, 'xi'
#mu and sigma behave like location and scale parameters in other distributions
#nu is fat-tailed Laplace at 1, Gaussian at 2, and Uniform as nu->inf 
#nu < 1 going towards a min of ~0.1 is progressively more fat-tailed Cauchy distribution
#xi is 1 for symmetric, 0.1 for highly left skewed, and 10 for highly right skewed

#SEP params
#theta: location parameter, 'mode'
#sigma: scale parameter, 'sigma'
#beta: shape parameter, 'beta'
#alpha: skew parameter, 'alpha'
#theta and sigma behave like location and scale parameters in SGED distribution but are not necessarily the same value
#beta is fat-tailed Laplace at 1, Gaussian at 0, and Uniform at -1 
#beta > 1 is progressively more fat-tailed Cauchy distribution
#alpha is O.5 for symmetric, 1 for highly left skewed, and 0 for highly right skewed


#SGED parameterizations
#Gaussian(0,1): beta=0, xi=1
#Laplace(0,1): beta=1, xi=1
#Uniform(0,1): beta=-1, xi=1
#add skew to any above by xi=10 for right skewed or xi=0.1 for left skewed


#1) Generate distributions from R fGarch 'skew generalized error distribution (SGED)' which is alternate parameterization of SEP
mn=0    #(-inf:inf)
sdev=1  #0:inf
nu=0.5  #0.1:inf
xi=1    #0.1:10
samps=1000

x<-rsged(samps,mean=mn,sd=sdev,nu=nu,xi=xi)

#2) Fit SEP model to SGED generated distribution via MLE
#SEP optimization
#theta bounds: -inf:inf, sigma bounds: >0, beta bounds: -1:10, alpha bounds: 0:1
#pars[theta,sigma,beta,alpha]

par_upr<-c(5,10,5,1)
par_est<-c(0,1,0,0.5)
par_lwr<-c(-5,0,-1,0)

sep_fit<-optim(par=par_est,sep_mle,x=x,
               method = 'L-BFGS-B',lower = par_lwr,upper = par_upr,
               control = list(fnscale=-1,maxit=100000))
sep_fit$par
theta<-sep_fit$par[1]
sigma<-sep_fit$par[2]
beta<-sep_fit$par[3]
alpha<-sep_fit$par[4]

#SEP parameterizations
#Gaussian(0,1): beta=0, alpha=0.5
#Laplace(0,1): beta=1, alpha=0.5
#Uniform(0,1): beta=-1, alpha=0.5
#add skew to any above by alpha=0 for right skewed or alpha=1 for left skewed

#3) Compare generated SGED vs fitted SEP distributions

#histogram comparison between R SGED sample (x) and fitted SEP pdf (red line)
png('histogram-comparison.png',width=512,height=256)
par(mfrow=c(1,2))
hist(x,breaks=c(-100,seq(-5,5,0.1),100),freq=FALSE,xlim=c(-5,5), main='SGED v SEP pdf',xlab='')
lines(seq(-5,5,0.01),sep_pdf(seq(-5,5,0.01),theta,sigma,beta,alpha),col='red')

hist(sep_samp(samps,theta,sigma,beta,alpha),breaks=c(-100,seq(-5,5,0.1),100),freq=FALSE,xlim=c(-5,5),col='yellow', main='SEP samples v SEP pdf',xlab='')
lines(seq(-5,5,0.01),sep_pdf(seq(-5,5,0.01),theta,sigma,beta,alpha),col='red')

dev.off()

#check that cdf and quantile functions behave appropriately
F<-sep_cdf(sort(x),theta,sigma,beta,alpha)
Q<-sep_quant2(F,theta,sigma,beta,alpha)

hist(Q,breaks=c(-100,seq(-5,5,0.1),100),freq=FALSE,xlim=c(-5,5), main='SEP cdf-quant v SEP pdf',xlab='')
lines(seq(-5,5,0.01),sep_pdf(seq(-5,5,0.01),theta,sigma,beta,alpha),col='red')
          
#P-P plot
png('P-P_plot.png')
par(mfrow=c(1,1))
emp_cdf<-rank(sort(x))/(length(x)+1)
est_cdf<-sep_cdf(sort(x),theta,sigma,beta,alpha)

plot(est_cdf,emp_cdf,xlim=c(0,1),ylim=c(0,1))
abline(0,1,col='gray')

dev.off()


#####################################END####################################
