require(parallel)
# detectCores() #12

N.reps = 10000
large = 1
small = 0.1



#### Hypoithsis: alpha=0 (no-mediaiton) ####
med_prop = 0 #hypothesised mediation proportion

abc_list<- list(c(0,0,1),c(0,large,1),c(large,0,1),c(small,small,1))
N_vals = as.list(seq(50,500,by=50))
fns_list = list('gd_bin','gd_bin_M','gd_bin_XY')
## list all deired simulation params (more highly parellelised)
sim_list = expand.grid(fns_list,N_vals,abc_list,1:N.reps)
sim_list = apply(sim_list,1,function(x) list(gdata = as.character(x[[1]]), nabc = c(x[[2]],x[[3]])))


## Run the simulation
system.time ( test<- mclapply(sim_list,function(x){
  df = do.call(x$gdata,list(x$nabc))
  output = tryCatch({
    gen_score_theta_bin(df,med_prop) #,true.param
  },error = function(e){
    #print('using scaled one')
    tryCatch({
      gen_score_theta_bin(df,med_prop,step.scale=0.01) #,true.param
    },error = function(e){
      print('A problem...')
      c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0)
    })
  })
  
  return(list(output,x$nabc,x$gdata))
},mc.cores=12)

)

##format results
l = length(test)
test_small = unlist(test,recursive = F)
results_full =  as.data.frame(list(do.call(rbind,test_small[(3*(1:l)-2)]),
                                   do.call(rbind,test_small[(3*(1:l)-1)]),
                                   as.character(test_small[(3*(1:l))]))    )

names(results_full) <- c('CUE_score','Tstep_score','Robust_wald','T1','T2','T3',
                         'Sobel','LR','T1_ols','T2_ols','T3_ols',
                         'beta1','beta2','beta3','beta1_cue','beta2_cue','beta3_cue',
                         'beta1_2s','beta2_2s','beta3_2s','N','a','b','c','gdata')

results_summary = aggregate(cbind(CUE_score,Tstep_score,Robust_wald) ~ N + a + b + c + gdata,
                            data=results_full,FUN=function(set)mean(set>qchisq(0.95,df=1)))
## Save the data
write.table(results_full,file='run_10k_full_bin_1.txt')
write.table(results_summary,file='run_10k_sum_bin_1.txt')



#### Repeat for test alpha=1 (no-direct effect) ####
med_prop = 1 

abc_list<- list(c(0,0,0),c(0,large,0),c(large,0,0),c(small,small,small))
N_vals = as.list(seq(50,500,by=50))
fns_list = list('gd_bin','gd_bin_MX','gd_bin_Y')

## list all deired simulation params
sim_list = expand.grid(fns_list,N_vals,abc_list,1:N.reps)
sim_list = apply(sim_list,1,function(x) list(gdata = as.character(x[[1]]), nabc = c(x[[2]],x[[3]])))


## Run the simulation
system.time ( test2 <- mclapply(sim_list,function(x){
  df = do.call(x$gdata,list(x$nabc))
  output = tryCatch({
    gen_score_theta_bin(df,med_prop) #,true.param
  },error = function(e){
    #print('using scaled one')
    tryCatch({
      gen_score_theta_bin(df,med_prop,step.scale=0.01) #,true.param
    },error = function(e){
      print('A problem...')
      c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0)
    })
  })
  
  return(list(output,x$nabc,x$gdata))
},mc.cores=12)

)

##format results
l = length(test2)
test_small_2 = unlist(test2,recursive = F)
results_full_2 =  as.data.frame(list(do.call(rbind,test_small_2[(3*(1:l)-2)]),
                                   do.call(rbind,test_small_2[(3*(1:l)-1)]),
                                   as.character(test_small_2[(3*(1:l))]))    )

names(results_full_2) <- c('CUE_score','Tstep_score','Robust_wald','T1','T2','T3',
                         'Sobel','LR','T1_ols','T2_ols','T3_ols',
                         'beta1','beta2','beta3','beta1_cue','beta2_cue','beta3_cue',
                         'beta1_2s','beta2_2s','beta3_2s','N','a','b','c','gdata')


results_summary_2 = aggregate(cbind(CUE_score,Tstep_score,Robust_wald) ~ N + a + b + c + gdata,
                            data=test3,FUN=function(set)mean(set>qchisq(0.95,df=1)))
## Save the data
write.table(test3,file='run_10k_full_bin_2.txt')
write.table(results_summary_2,file='run_10k_sum_bin_2.txt')
