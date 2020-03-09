require(ggplot2)
require(gridExtra)
require(reshape2)

fns = list("gd_bin","gd_bin_M","gd_bin_MX","gd_bin_XY","gd_bin_Y" )
## Simulated Data sets for No-mediation and No-Natural-direct-effect hypohteses
results_nomed <- read.table(file='Sim_data/run_10k_full_bin_1.txt')
results_nonde <- read.table(file='Sim_data/run_10k_full_bin_3.txt')


results_nomed_sum = aggregate(cbind(CUE_score,Tstep_score,Robust_wald,Sobel,LR) ~ N + a + b + c + gdata,
                            data=results_nomed,FUN=function(set)mean(set>qchisq(0.95,df=1)))

results_nonde_sum = aggregate(cbind(CUE_score,Tstep_score,T3,T3_ols) ~ N + a + b + c + gdata,
                              data=results_nonde,FUN=function(set)mean(set>qchisq(0.95,df=1)))

write.table(results_nomed_sum,file='Sim_data/run_10k_sum_bin_1.txt')
write.table(results_nonde_sum,file='Sim_data/run_10k_sum_bin_3.txt')



## Get summary tables

results_nomed_sum <- read.table(file='Sim_data/run_10k_sum_bin_1.txt')
results_nonde_sum <- read.table(file='Sim_data/run_10k_sum_bin_3.txt')

melt_res <- function(df) melt(df,id.vars = c('N','a','b','c','gdata'),
                              variable.name = 'Test',value.name='Rej')

results_nomed_sum = melt_res(results_nomed_sum)
results_nonde_sum = melt_res(results_nonde_sum)



hyp_plot <- function(df){
  df$Group = paste0('beta = \n(',df$a,',',df$b,',',df$c,')')
  df$Group = factor(df$Group,levels = #Put plots not in alphabetical order
                      levels(factor(df$Group))[c(1,2,4,3)])
  
  df$gdata = factor(df$gdata,levels=fns[c(1,3,4,2,5)],labels=
                      c('MXY','MX','XY','M','Y'))
  
  p = ggplot(data=df) + 
    geom_point(aes(x=N,y=Rej,col=Test),show.legend = TRUE,alpha=0.8) +
    geom_hline(yintercept = 0.05,lty='dashed')+
    facet_grid(gdata~Group)+
    theme_bw()+
    xlab('Sample size, n')+
    ylab('Rejection Proportion')+
    theme(plot.title=element_text(size=16,face="bold"),
          #axis.text.x=element_text(face="bold"),
          #axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=0,face="bold") ,
          strip.text.x = element_text(margin = margin(2, 0, 2, 0))
    )
  
  return(p)
}

pdf('plots/hyp_plt_nomed.pdf',width=8,height=7)
hyp_plot(results_nomed_sum)
dev.off()

pdf('plots/hyp_plt_nonde.pdf',width=8,height=7)
hyp_plot(results_nonde_sum)
dev.off()


