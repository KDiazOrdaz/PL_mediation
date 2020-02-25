require(parallel)

####### Run the Simulation for Bias_variance plots ######


N.reps = 1000
abc_list<- list(c(1,1,1),c(0,0,0))
N_vals = as.list(c(100,500,1000))
fns_list = list('gd_bin','gd_bin_M','gd_bin_XY','gd_bin_Y','gd_bin_MX')

## list all deired simulation params
sim_list = expand.grid(fns_list,N_vals,abc_list,1:N.reps)
sim_list = apply(sim_list,1,function(x) list(gdata = as.character(x[[1]]), nabc = c(x[[2]],x[[3]])))

## Run the simulation

system.time ( bias_var_results <- mclapply(sim_list,function(x){
  df = do.call(x$gdata,list(x$nabc))
  output = c(fit_beta(df),boot_beta(df))
  return(list(output,x$nabc,x$gdata))
  
  },mc.cores=12)
) #Estimated 14 hours run time


##format results
l = length(bias_var_results)
test_small = unlist(bias_var_results,recursive = F)
results_full =  as.data.frame(list(do.call(rbind,test_small[(3*(1:l)-2)]),
                                   do.call(rbind,test_small[(3*(1:l)-1)]),
                                   as.character(test_small[(3*(1:l))]))    )

names(results_full) <- c('beta1','beta2','beta3','nide',
                         'var1.hat','var2.hat','var3.hat','var4.hat',
                         'var1.boot','var2.boot','var3.boot','var4.boot',
                         'N','a','b','c','gdata')

results_summary = aggregate(cbind(beta1,beta2,beta3,nide,var1.hat,var2.hat,var3.hat,var4.hat,
                                  var1.boot,var2.boot,var3.boot,var4.boot) ~ N + a + b + c + gdata,
                            data=results_full,FUN=function(set){
                              mean(set)
                            })

results_summary2 = aggregate(cbind(beta1,beta2,beta3,nide,var1.hat,var2.hat,var3.hat,var4.hat,
                                  var1.boot,var2.boot,var3.boot,var4.boot) ~ N + a + b + c + gdata,
                            data=results_full,FUN=function(set){
                              sd(set)/sqrt(N.reps)
                            })


## Save the data
write.table(results_full,file='Sim_data/run_bv_1.txt')
write.table(results_summary ,file='Sim_data/run_bv_1_mean_sum.txt')
write.table(results_summary2,file='Sim_data/run_bv_1_sd_sum.txt')
  




