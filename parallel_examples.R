
library(parallel)
library(doParallel)

source('sep_function.R')

#determine no of available cores
n.cores <- parallel::detectCores(logical=FALSE)
print(n.cores)

#start internal cluster
cl<-makeCluster(n.cores,type='PSOCK')
print(cl)

#register cluster for 'DoParallel' that contains 'foreach' functionality
doParallel::registerDoParallel(cl = cl)
foreach::getDoParRegistered()

#create and nxK vector of standard normal samples
k_samps<-100
samp_size<-10000

sample_vec<-rep(samp_size,k_samps)

norm_samps<-sapply(sample_vec,rnorm) #vectorize apply of rnorm function
norm_samps_lst<-lapply(sample_vec,rnorm) #list-vector apply of rnorm function

#---------------------------------------------------------------
#run the fit sequence with a standard for loop and output time
fit_vec<-vector('list',dim(norm_samps)[2])

sepfit_noparallel<-system.time(
    for(i in 1:dim(norm_samps)[2]){
      fit_vec[[i]]<-SEPfit(norm_samps[,i])
    }
)

print(paste('for-loop',sepfit_noparallel[3],'sec'))

#---------------------------------------------------------------
#run the fit sequence with a standard for loop and output time
sepfit_noparallel_vec<-system.time(
  apply(norm_samps,2,SEPfit)
)

print(paste('apply',sepfit_noparallel_vec[3],'sec'))

#----------------------------------------------------------------
#run with a 'foreach' loop
sepfit_foreach<-system.time(
  foreach_vec<-foreach(i = 1:dim(norm_samps)[2],.combine='c',.inorder=F,.packages=c('stats'),.export=c('.GlobalEnv'))%dopar%{
    out<-SEPfit(norm_samps[,i])
    return(out)
  }
)

print(paste('foreach',sepfit_foreach[3],'sec'))

#----------------------------------------------------------------
#run with parallel 'mclapply'; only works on UNIX (not windows)
sepfit_mclapply<-system.time(
  mclapply_vec<-mclapply(norm_samps_lst,SEPfit,mc.cores=n.cores)
)

print(paste('mclapply',sepfit_mclapply[3],'sec'))

#----------------------------------------------------------------------
#run with parallel 'ClusterApply'
src_in<-clusterEvalQ(cl=cl,source('sep_function.R')) #loads 'sep_function.R' on all worker nodes
fun_in<-clusterExport(cl=cl,varlist = c('SEPfit','sep_pdf','sep_mle','kfun','zfun')) #individually loads functions in the current environment
sepfit_clapply<-system.time(
  clapply_vec<-clusterApply(cl=cl,norm_samps_lst,SEPfit)
)

print(paste('clapply',sepfit_clapply[3],'sec'))

#----------------------------------------------------------------------
#run with parallel 'ClusterApplyLB'; LB means load balancing, a more flexible implementation
src_in<-clusterEvalQ(cl=cl,source('sep_function.R')) #loads 'sep_function.R' on all worker nodes

sepfit_clapplylb<-system.time(
  clapply_vec<-clusterApplyLB(cl=cl,norm_samps_lst,SEPfit)
)

print(paste('clapply-lb',sepfit_clapplylb[3],'sec'))


#----------------------------------------------------------------------
#run with parallel 'parCapply'; LB means load balancing, a more flexible implementation
src_in<-clusterEvalQ(cl=cl,source('sep_function.R')) #loads 'sep_function.R' on all worker nodes

sepfit_parcapply<-system.time(
  clapply_vec<-parCapply(cl=cl,norm_samps,SEPfit)
)

print(paste('parcapply',sepfit_parcapply[3],'sec'))



#stop the cluster at the end
stopCluster(cl)

#####################################END###################################
