# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:12:25 2023

@author: zpb4
"""
import matplotlib.pyplot as plt
import numpy as np
from sep_function import *
#------------------------------------------------------------------------------
#define SEP parameterization to test here:
theta=0.5
sigma=0.5
beta=1
alpha=0.6

plot_out=False

#plotting inputs
bins=np.linspace(-100,100,1000)
x1 = np.linspace(-100,100,10000)
x = sep_samp(10000,theta,sigma,beta,alpha)
xlm=(-5,5)

#plot empirical histogram of SEP sampled values (x) vs SEP density
if plot_out==True:
    emp_hist = plt.hist(x,bins=bins,density=True)
    y = sep_pdf(x1,theta,sigma,beta,alpha)
    pdf_lin = plt.plot(x1,y,color='red')
    plt.xlim(xlm[0],xlm[1])
    #plt.show()
    plt.savefig('plot/py-pdf-check.png')

#run sampled values (x) through SEP cdf and then quantile function to ensure correct operation
F = sep_cdf(x,theta,sigma,beta,alpha)
Q = sep_quant(F,theta,sigma,beta,alpha)

#plot derived quantiles of x from F and Q transforms above and compare to original SEP pdf
if plot_out==True:
    q_hist = plt.hist(Q,bins=bins,density=True)
    y = sep_pdf(x1,theta,sigma,beta,alpha)
    pdf_lin = plt.plot(x1,y,color='red')
    plt.xlim(-5,5)
    #plt.show()
    plt.savefig('plot/py-cdf-quant-check.png')  #save figure


"""
Test MLE with log likelihood function
"""

#scipy imports
from scipy.optimize import Bounds
from scipy.optimize import minimize
from scipy.optimize import BFGS


#MLE with analytic sigma estimate via 'sep_mle_sigest' function
bounds = Bounds([-5, -0.99, 0.01], [5, 5, 0.99])
pars = np.array([0, 0, 0.5])
res = minimize(sep_mle_sigest, pars,args=(x), method='trust-constr',jac="2-point",hess=BFGS(),options={'verbose': 1}, bounds=bounds)
theta_est = res.x[0]
beta_est = res.x[1]
alpha_est = res.x[2]
sigma_est = sep_sigma_est(x,theta_est,beta_est,alpha_est)


#compare original parameters to MLE estimated parameters
print(theta,sigma,beta,alpha)
print(theta_est,sigma_est,beta_est,alpha_est)

np.save('out/py-sep-mle-est.npy',([theta_est,sigma_est,beta_est,alpha_est]))

#compare empirical distribution (histogram) to fitted distribution (red line)
if plot_out==True:
    emp_hist = plt.hist(x,bins=bins,density=True)
    y = sep_pdf(x1,theta_est,sigma_est,beta_est,alpha_est)
    pdf_lin = plt.plot(x1,y,color='red')
    plt.xlim(-5,5)
    #plt.show()
    plt.savefig('plot/py-mle-check.png') #save figure


#-------------------------------------------------------------------------------------------------