import mpmath
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import uniform
from scipy.stats import gamma


def kfun(beta,alpha):
  kinv = math.gamma(1+((beta+1)/2)) * 2**(1+0.5*(1+beta)) / (4*alpha*(1-alpha))
  k = 1/kinv
  return k

def zfun(x,theta,sigma):
  z = (x-theta)/sigma
  return z

#SEP pdf function
def sep_pdf(x,theta,sigma,beta,alpha):
  k = kfun(beta,alpha)
  z = zfun(x,theta,sigma)
  y = (k/sigma) * np.exp(-0.5*(abs(z)+(2*alpha-1)*z)**(2/(1+beta)))
  return y


#SEP cdf function
def sep_cdf(x,theta,sigma,beta,alpha):
  z = zfun(x,theta,sigma)
  y = np.zeros(shape=[len(x)])
  for i in range(0,len(z)):
    if z[i]<=0:
      gam_x = 2**((2/(beta+1))-1) * ((alpha-1)*z[i])**(2/(beta+1))
      gam_a = ((beta+1)/2)
      num = alpha*mpmath.gammainc(gam_a,gam_x)
      den = math.gamma((beta+1)/2)
      y[i] = num/den
    if z[i]>0:
      gam_x = 2**((2/(beta+1))-1)*(1/(alpha*z[i]))**(-2/(beta+1))
      gam_a = ((beta+1)/2)
      num = (1-alpha)*mpmath.gammainc(gam_a,gam_x,)
      den = math.gamma((beta+1)/2)
      y[i] = 1-(num/den)
  return(y)


#SEP quantile function
def sep_quant(u,theta,sigma,beta,alpha):
  q = np.zeros(shape=[len(u)])
  for i in range(0,len(u)):
    if u[i] <= alpha:
      num = (2**(0.5-(1/(beta+1))) * np.sqrt(gamma.isf((1-(u[i]/alpha)),((beta+1)/2))))**(beta+1)
      den = 1-alpha
      q[i] = theta - sigma*(num/den)
    if u[i] > alpha:
      num = (2**((1/(beta+1))-0.5) / np.sqrt(gamma.isf((1-((1-u[i])/(1-alpha))),((beta+1)/2))))**(-(beta+1))
      den = alpha
      q[i] = theta + sigma*(num/den)

  return(q)


#SEP sampling function
def sep_samp(n,theta,sigma,beta,alpha):
  u = uniform.rvs(size=n)
  samps = sep_quant(u,theta,sigma,beta,alpha)
  return(samps)

#analytic sigma estimate
def sep_sigma_est(x,theta,beta,alpha):
  n = len(x)
  num = np.sum((np.abs(x-theta)+(2*alpha-1)*(x-theta))**(2/(beta+1)))
  den = n*(beta+1)
  sigma = (num/den)**((beta+1)/2)
  return(sigma)

#mle function with all params
def sep_mle(pars,x):
  theta = pars[0]
  sigma = pars[1]
  beta = pars[2]
  alpha = pars[3]
  try:
      ll = np.sum(np.log(sep_pdf(x,theta,sigma,beta,alpha)))
  except ZeroDivisionError:
      ll = (-1e9)
  return(-ll)

#mle function with analytic sigma estimate
def sep_mle_sigest(pars,x):
  theta = pars[0]
  beta = pars[1]
  alpha = pars[2]
  sigma = sep_sigma_est(x,theta,beta,alpha)
  try:
      ll = np.sum(np.log(sep_pdf(x,theta,sigma,beta,alpha)))
  except ZeroDivisionError:
      ll = (-1e9)
  return(-ll)

#------------------------------------------------------------------------------
