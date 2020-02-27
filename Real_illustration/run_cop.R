
df = cop_treat_complete

Z_formula = "~Gender+age+hads_dep_base+hads_anx_base+PSEQ_overall_base+
HEIQ_base+CPAQ_base+Housing+English+Employment+
Education+employed_ftedu+centre"

Z_formula = "~Gender+poly(age,2)+hads_anx_base+employed_ftedu*centre+hads_dep_base+
    PSEQ_overall_base+
HEIQ_base+CPAQ_base"

gamx = glm(paste0('X',Z_formula),family=binomial,data=df)$coefficients
modm = glm(paste0('M',Z_formula,'+X'),family=gaussian,data=df)
mody = glm(paste0('Y',Z_formula,'+X+M'),family=gaussian,data=df)

T1_ols = (summary(modm)$coefficients['X','t value'])^2
T2_ols = (summary(mody)$coefficients['M','t value'])^2
T3_ols = (summary(mody)$coefficients['X','t value'])^2
Sobel  = T1_ols*T2_ols /(T1_ols +T2_ols )
LR = min(T1_ols,T2_ols)
ols_stats = c(Sobel,LR,T1_ols,T2_ols,T3_ols)

initm = modm$coefficients
inity = mody$coefficients

beta1 = initm['X']
beta23 = inity[c('M','X')]
gamm = initm[!(names(initm)%in% c('X'))]
gamy = inity[!(names(inity)%in% c('M','X'))]

gamma0 = c(gamx,gamx,gamm,gamm,gamy,gamy)  #Intial MLE estimate for gamma
beta0 = c(beta1,beta23)

theta = c(beta0,gamma0) #MLE estimate of theta

#beta0[1]*beta0[2] /(beta0[1]*beta0[2] + beta0[3])




Max.it = 1000
prec = 1e-12
step.scale = 0.01
converged = FALSE
i = 1

while(i<Max.it& !converged){
  step = CUE_vec_J_cop(Z_formula,df,theta)
  theta = theta - solve(step$J_unconstr,step$vec_unconstr)
  if(sum((abs(step$vec)>=prec))==0){converged = TRUE}
  i = i+1
}

#if (!converged){cat('Unconstrained not converged',signif(max(step$vec),3),'\n')}
beta_unc = theta[1:3]
gamma_unc = theta[-c(1:3)]
step = CUE_vec_J_cop(Z_formula,df,theta)
A.mat = step$A.mat
RSTest = step$RSTest
T_stats = step$T_stats




## Do the CUE bit
med_prop = 0

theta = c(beta_unc,gamma_unc,0)
converged = FALSE
i = 1
while(i<Max.it& !converged){
  step = Lagrangian_cue_cop(Z_formula,df,theta,med_prop)
  theta = theta - step.scale* solve(step$J,step$vec) #step size reduced by 0.25 see: https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
  
  if(sum((abs(step$vec)>=prec))==0){converged = TRUE}
  i = i+1
}

beta_cue = theta[1:3]
cue_score_indirect = CUE_vec_J_cop(Z_formula,df,theta[-length(theta)])$score


### Do the 2step bit
med_prop = 0
beta = c(beta_unc,0)
converged = FALSE
i = 1
while(i<2*Max.it& !converged){
  if(sum((abs(step$vec)>=prec))==0){converged = TRUE}
  i = i+1
  step = Lagrangian_2step_cop(Z_formula,df,beta,med_prop,gamma=gamma_unc,A = A.mat)
  beta = beta - step.scale* solve(step$J,step$vec) #step size reduced by 0.25 see: https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
}
if (!converged){cat('2-step not converged',signif(max(step$vec),3),'\n')}
beta_2step = beta[1:3]
Tstep_indirect = step$score


#### Direct effect test
#cue
med_prop = 1
theta = c(beta_unc,gamma_unc,0)
converged = FALSE
i = 1
while(i<Max.it& !converged){
  step = Lagrangian_cue_cop(Z_formula,df,theta,med_prop)
  theta = theta - step.scale* solve(step$J,step$vec) #step size reduced by 0.25 see: https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
  
  if(sum((abs(step$vec)>=prec))==0){converged = TRUE}
  i = i+1
}
cue_score_direct = CUE_vec_J_cop(Z_formula,df,theta[-length(theta)])$score

#2step
med_prop = 1
beta = c(beta_unc,0)
converged = FALSE
i = 1
while(i<Max.it& !converged){
  step = Lagrangian_2step_cop(Z_formula,df,beta,med_prop,gamma=gamma_unc,A = A.mat)
  beta = beta - step.scale* solve(step$J,step$vec) #step size reduced by 0.25 see: https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
  if(sum((abs(step$vec)>=prec))==0){converged = TRUE}
  i = i+1
}
if (!converged){cat('2-step not converged',signif(max(step$vec),3),'\n')}
beta_2step_2 = beta[1:3]
Tstep_direct = derivs_2step_cop(Z_formula,df,beta_2step_2,gamma_unc,A=A.mat)$score









out = c(ols_stats ,
        cue_score_indirect,
          Tstep_indirect,
         RSTest, T_stats,cue_score_direct,Tstep_direct,
         beta_unc,beta_cue,beta_2step)

names(out) <- c( 'Sobel','LR','T1_ols','T2_ols','T3_ols',
                  'CUE_score','Tstep_score','Robust_wald','T1','T2','T3','cue_score_direct',
                 'Tstep_direct',
                         'beta1','beta2','beta3','beta1_cue','beta2_cue','beta3_cue',
                         'beta1_2s','beta2_2s','beta3_2s')



pvals = 1-pchisq(out[1:13],df=1)
p_frame = data.frame(P_value=round(pvals,7),Signif = ifelse(pvals<=0.1,ifelse(pvals<=0.05,ifelse(pvals<=0.01,ifelse(pvals<=0.001,'***','**'),'*'),'.'),' '))

v_frame = round(data.frame(OLS=beta0,OLS_sd = sqrt(beta0^2/ols_stats[3:5]),
                 G=beta_unc,G_sd = sqrt(beta_unc^2/T_stats),
                 CUE=beta_cue,Tstep = beta_2step),3)


nide = c(beta0[1]*beta0[2],
  sqrt((beta0[1]*beta0[2])^2 /out['Sobel']),
  beta_unc[1]*beta_unc[2],
  sqrt((beta_unc[1]*beta_unc[2])^2/out['Robust_wald']),
  0,0)
v_frame = rbind(v_frame,nide)

  

print(list(P_values = p_frame, Parameter_estimates=v_frame) )
