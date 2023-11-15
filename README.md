# hpc1_training
This package contains function scripts for the Skew Exponential Power (SEP) distribution as detailed in:  
   
  Hutson, A. D. (2019). An alternative skew exponential power distribution formulation. Communications in Statistics - Theory and Methods, 48(12), 3005–3024. https://doi.org/10.1080/03610926.2018.1473600  
     
These scripts are intended to support initial training in HPC operations  
* There is a R version of the pdf, cdf, quantile, and sampling functions and a Python equivalent ('sep_function.R', 'sep_function.py')
* For the R version, there is a separate 'SEP experimentation' script to compare the functionally equivalent Skew Generalized Error Distribution (SGED) from a prebuilt package in R to the Maximum Likelihood Estimated (MLE) distributional fit for the SEP. This demonstrates the differences in the values of the 4 parameters that are similar, but different, between the SEP and SGED
* The Python version includes comparison plots (commmented out because they don't work on Hopper) and a printout comparison between the parameters of a sample generated from SEP sampling function and the same distribution's parameter estimates via MLE
* There are also bash (*.sh) scripts that run the R and Python scripts ('run_sep-fun_r.sh','run_sep-fun_py.sh')
* As part of demonstrating HPC operations with these simple scripts, recommend highlighting
  * How to ensure package/module compatibility for both R and Python
  * Good coding practices such as use of functions, use of 'source' in R to bring in pre-built functions, and use of R-data structure (.rds) files to save output data
  * Challeng in seeing output data (i.e. plots) when operating on remote machines
  * Error handling in R/Python versus SLURM scheduled jobs
