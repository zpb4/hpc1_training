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
def sep_mle(pars):
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
def sep_mle_sigest(pars):
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
#define SEP parameterization to test here:
theta=0.5
sigma=0.5
beta=1
alpha=0.6

#plotting inputs
bins=np.linspace(-100,100,1000)
x1 = np.linspace(-100,100,10000)
x = sep_samp(10000,theta,sigma,beta,alpha)

#plot empirical histogram of SEP sampled values (x) vs SEP density
"""
emp_hist = plt.hist(x,bins=bins,density=True)
y = sep_pdf(x1,theta,sigma,beta,alpha)
pdf_lin = plt.plot(x1,y,color='red')
plt.xlim(xlm[0],xlm[1])
plt.show()
"""

#run sampled values (x) through SEP cdf and then quantile function to ensure correct operation
F = sep_cdf(x,theta,sigma,beta,alpha)
Q = sep_quant(F,theta,sigma,beta,alpha)

#plot derived quantiles of x from F and Q transforms above and compare to original SEP pdf
"""
q_hist = plt.hist(Q,bins=bins,density=True)
y = sep_pdf(x1,theta,sigma,beta,alpha)
pdf_lin = plt.plot(x1,y,color='red')
plt.xlim(-5,5)
plt.show()
"""

"""
Test MLE with log likelihood function
"""

#scipy imports
from scipy.optimize import Bounds
from scipy.optimize import minimize
from scipy.optimize import BFGS
from scipy.optimize import SR1
from scipy import optimize

"""
#optimizing the full parameter set [theta,sigma,beta,alpha] performed worse than the LL function
#with the analytic sigma estimate below

#MLE with all parameters via 'sep_mle' function
bounds = Bounds([-5, 0.000001, -0.99, 0.01], [5, 10, 5, 0.99])
pars = np.array([0, 1, 0, 0.5])
res = minimize(sep_mle, pars, method='trust-constr',jac="2-point",hess=BFGS(),options={'verbose': 1}, bounds=bounds)
estimated parameters
theta_est = res.x[0]
sigma_est = res.x[1]
beta_est = res.x[2]
alpha_est = res.x[3]
"""

"""
#Global optimization of MLE function
#did not always work well

bounds = [(-5,5),(0.0001,10),(-0.99,5),(0.01,0.99)]
res = optimize.shgo(sep_mle,bounds)
"""

#MLE with analytic sigma estimate via 'sep_mle_sigest' function
bounds = Bounds([-5, -0.99, 0.01], [5, 5, 0.99])
pars = np.array([0, 0, 0.5])
res = minimize(sep_mle_sigest, pars, method='trust-constr',jac="2-point",hess=BFGS(),options={'verbose': 1}, bounds=bounds)
theta_est = res.x[0]
beta_est = res.x[1]
alpha_est = res.x[2]
sigma_est = sep_sigma_est(x,theta_est,beta_est,alpha_est)



print(theta,sigma,beta,alpha)
print(theta_est,sigma_est,beta_est,alpha_est)

#compare empirical distribution (histogram) to fitted distribution (red line)
"""
emp_hist = plt.hist(x,bins=bins,density=True)
y = sep_pdf(x1,theta_est,sigma_est,beta_est,alpha_est)
pdf_lin = plt.plot(x1,y,color='red')
plt.xlim(-5,5)
plt.show()
"""


