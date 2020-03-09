###Bias and Variance plots

require(ggplot2)
require(gridExtra)
require(reshape2)
require(latex2exp)

#Import simulated data
bv_df <- merge(read.table('Sim_data/run_bv_1_mean_sum.txt'),
               read.table('Sim_data/run_bv_1_sd_sum.txt'),
               by=c('N','a','b','c','gdata'),
               suffixes=c('.MCmean','.MCsd'))


#### Create data frame of biases ####
fns = list("gd_bin","gd_bin_M","gd_bin_MX","gd_bin_XY","gd_bin_Y" )
params = sapply(1:4,function(i)paste0('beta',i))
N.reps=1000


bias_plot_df = bv_df[1:5]

bias_plot_df [,sapply(1:4,function(i)paste0('beta',i))] <- 
  bv_df[,c('beta1.MCmean','beta2.MCmean','beta3.MCmean','nide.MCmean')] -
  cbind(bv_df[,c('a','b','c')],bv_df$a*bv_df$b)

bias_plot_df [,sapply(1:4,function(i)paste0('beta',i,'.sd'))] <- 
  bv_df[,c('beta1.MCsd','beta2.MCsd','beta3.MCsd','nide.MCsd')] #see[1]

biases = melt(bias_plot_df,c('N','a','b','c','gdata'),value.name = 'bias.sd',
              measure.vars = sapply(params,paste0,'.sd') ,variable.name = 'param')

biases$param <- gsub(".sd","",biases$param)

biases = merge(biases,melt(bias_plot_df,c('N','a','b','c','gdata'),value.name = 'bias',
                           measure.vars = params,variable.name = 'param'))
biases$param <- gsub("bias","beta",biases$param)
rm(bias_plot_df)

valids <- function(df){
  incl = !(  (df$gdata=="gd_bin_M" & df$param %in% c('beta3'))               |
               (df$gdata=="gd_bin_Y"  & df$param %in% c('beta1','beta4'))    )
  return(df[incl,])
}

biases = valids(biases)



#### Bias plot ####


bias_plot <- function(df){
  scale = qnorm(0.975) #95% confidence interval
  alph = 0.6
  
  #Rename and reorder factor variables
  df$Group = paste0('N=',df$N,'\n(',df$a,',',df$b,',',df$c,')')
  
  df$Group = factor(df$Group,levels = #Put plots not in alphabetical order
                      levels(factor(df$Group))[c(1,2,5,6,3,4)])
  
  df$gdata = factor(df$gdata,levels=fns[c(1,3,4,2,5)],labels=
                      c('MXY','MX','XY','M','Y'))
  
  df$param = factor(df$param,levels=params[c(4,3,2,1)])
  
  #Make plot
  p=ggplot(data=df,aes(x=param,y=bias,ymin=bias-scale*bias.sd,ymax=bias+scale*bias.sd,color=param))+
    geom_pointrange(aes(col=param),show.legend = TRUE,shape=5,size=0.1)+
    geom_errorbar(aes(ymin=bias-scale*bias.sd,ymax=bias+scale*bias.sd,col=param),width=0.2,cex=0.5,show.legend = TRUE)+
    geom_hline(aes(yintercept =0))+
    xlab('')+
    ylab('Bias')+
    #facet_wrap(~Group,strip.position="left",nrow=6,scales = "free_y") +
    facet_grid(Group~gdata,switch='y')+
    theme_bw()+
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          #axis.text.x=element_text(face="bold"),
          #axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold") ,
          strip.text.x = element_text(margin = margin(2, 0, 2, 0))
    )+coord_flip()
  return(p)
}

pdf('plots/bias_plt.pdf',width=8,height=7)
bias_plot(biases)
dev.off()
  
 
#### Create data frame of vairances ####

#[1] https://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf
#for standard error in the sample variance (and also unbiased standard deviation estimate)  
  
var_plot_df = bv_df[1:5]

var_plot_df [,sapply(1:4,function(i)paste0('MC',i))] <- 
  bv_df[,c('beta1.MCsd','beta2.MCsd','beta3.MCsd','nide.MCsd')]^2*N.reps

var_plot_df [,sapply(1:4,function(i)paste0('MC',i,'.sd'))] <- 
  sqrt(2/(N.reps-1))*bv_df[,c('beta1.MCsd','beta2.MCsd','beta3.MCsd','nide.MCsd')]^2*N.reps  #see[1]


#var_plot_df [,sapply(1:4,function(i)paste0('MC',i,'.sd'))] <- 
#  results_summary3[,c('beta1','beta2','beta3','nide')]


var_plot_df [,sapply(1:4,function(i)paste0('Asym',i))] <- 
  bv_df[,sapply(1:4,function(i)paste0('var',i,'.hat.MCmean'))]

var_plot_df [,sapply(1:4,function(i)paste0('Asym',i,'.sd'))] <- 
  bv_df[,sapply(1:4,function(i)paste0('var',i,'.hat.MCsd'))]

var_plot_df [,sapply(1:4,function(i)paste0('Boot',i))] <- 
  bv_df[,sapply(1:4,function(i)paste0('var',i,'.boot.MCmean'))]

var_plot_df [,sapply(1:4,function(i)paste0('Boot',i,'.sd'))] <- 
  bv_df[,sapply(1:4,function(i)paste0('var',i,'.boot.MCsd'))]


melt_merge <- function(st){
  var1 = melt(var_plot_df,c('N','a','b','c','gdata'),value.name = 'Var.sd',
              measure.vars = sapply(1:4,function(x) paste0(st,x,'.sd')) ,
              variable.name = 'param')
  var1$param <- var1$param <- gsub(".sd","",var1$param)
  var2 = melt(var_plot_df,c('N','a','b','c','gdata'),value.name = 'Var',
              measure.vars = sapply(1:4,function(x) paste0(st,x)) ,
              variable.name = 'param')
  var1 = merge(var1,var2)
  var1$Method <- st
  var1$param <- gsub(st,"beta",var1$param)
  return(var1)
}

variances = do.call(rbind,lapply(list('MC','Asym','Boot'),melt_merge))
variances = valids(variances)
rm(var_plot_df)
variances[,c('Var','Var.sd')] <- variances[,c('Var','Var.sd')] *cbind(variances$N,variances$N)


#### Create plot of vairances ####


variance_plot <- function(df,par){ #take par as an asgument to create different plots for each parameter
  scale = qnorm(0.975) #95% confidence interval
  alph = 0.6
  df = subset(df,param==par)
  
  #Rename and reorder factor variables
  df$Group = paste0('N=',df$N,'\n(',df$a,',',df$b,',',df$c,')')
  
  df$Group = factor(df$Group,levels = #Put plots not in alphabetical order
                      levels(factor(df$Group))[c(1,2,5,6,3,4)])
  
  df$gdata = factor(df$gdata,levels=fns[c(1,3,4,2,5)],labels=
                      c('MXY','MX','XY','M','Y'))
  
  df$Method = factor(df$Method)
  
  lat_name = paste0(gsub('beta','$\\\\hat{\\\\beta_',par),'}$')
  if(par=='beta4') lat_name = '$\\hat{\\beta}_1\\hat{\\beta}_2$'
  
  #Make plot
  p=ggplot(data=df,aes(x=Method,y=Var,ymin=Var-scale*Var.sd,ymax=Var+scale*Var.sd,color=Method))+
    geom_pointrange(aes(col=Method),show.legend = TRUE,shape=5,size=0.1)+
    geom_errorbar(aes(ymin=Var-scale*Var.sd,ymax=Var+scale*Var.sd,color=Method),width=0.2,cex=0.5,show.legend = TRUE)+
    #geom_hline(aes(yintercept =1))+
    xlab('')+
    ylab(TeX(paste0('Variance of ',lat_name)))+
    #facet_wrap(~Group,strip.position="left",nrow=6,scales = "free_y") +
    facet_grid(Group~gdata,switch='y')+
    theme_bw()+
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          #axis.text.x=element_text(face="bold"),
          #axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold") ,
          strip.text.x = element_text(margin = margin(2, 0, 2, 0))
    )+coord_flip()
  return(p)
}

variance_plot(variances,'beta1')


sapply(params,function(par){ #saveplots
  pdf(paste0('plots/var_plt_',par,'.pdf'),width=8,height=7)
  print(variance_plot(variances,par))
  dev.off()
})







