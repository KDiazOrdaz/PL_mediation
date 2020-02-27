#Good CUE code
#gives the CUE vector and J matrix for newton raphson
#gives the unconstrained vector and J for newton raphson

#CUE_fitting function addapted to allow for a formula for confounders

CUE_vec_J_cop <- function(formula,df,theta,brow=F){
  
  #### Read in parameters and construct residuals ####
  beta = theta[1:3]
  gamma = theta[-(1:3)]
  beta1 = beta[1]
  beta2 = beta[2]
  delta  = beta[3]
  
  p = length(gamma)/6
  gamma_x1 = gamma[1:p]
  gamma_x2 = gamma[(p+1)  :(2*p)]
  gamma_m1 = gamma[(2*p+1):(3*p)]
  gamma_m2 = gamma[(3*p+1):(4*p)]
  gamma_y1 = gamma[(4*p+1):(5*p)]
  gamma_y2 = gamma[(5*p+1):(6*p)]
  
  X = df$X
  M = df$M
  Y = df$Y
  N = length(M)
  Z.matrix = model.matrix(as.formula(formula),data = df)
  
  odds_x1  = Z.matrix%*%gamma_x1 #The log-odds
  odds_x2  = Z.matrix%*%gamma_x2
  f1       = Z.matrix%*%gamma_m1
  f2       = Z.matrix%*%gamma_m2
  g1       = Z.matrix%*%gamma_y1 
  g2       = Z.matrix%*%gamma_y2 
  
  
  e1 = exp(odds_x1)
  e2 = exp(odds_x2)
  nt = 1 # Number of trials for X~ Binomial (nt,p(Z))

  x.res1   = X-nt*e1/(1+e1)
  x.res2   = X-nt*e2/(1+e2)
  m.res1   = M - beta1*X - f1 
  m.res2   = M - beta1*X - f2 
  y.res1   = Y - beta2*M - delta*X - g1
  y.res2   = Y - beta2*M - delta*X - g2
  
  
  hdot1     = nt*e1/(1+e1)^2  #derivative of E[X|Z] wrt. log-odds (=odds_x1)
  hdot2     = nt*e2/(1+e2)^2
  
  #### Estimating Equations and first derivatives ####
  
  U1 = sum(x.res1*m.res1)/N
  U2 = sum(m.res2*y.res1)/N
  U3 = sum(x.res2*y.res2)/N
  U = c(U1,U2,U3)
  
  
  V1 = t(Z.matrix) %*% (m.res1*hdot1)/N
  V2 = t(Z.matrix) %*% x.res1/N
  V3 = t(Z.matrix) %*% y.res1/N
  V4 = t(Z.matrix) %*% m.res2/N
  V5 = t(Z.matrix) %*% (y.res2*hdot2)/N
  V6 = t(Z.matrix) %*% x.res2/N
  
  sig11 = sum(x.res1^2*m.res1^2)/N
  sig22 = sum(m.res2^2*y.res1^2)/N
  sig33 = sum(x.res2^2*y.res2^2)/N
  
  sig12 = sum(x.res1*m.res1*m.res2*y.res1)/N
  sig13 = sum(x.res1*m.res1*x.res2*y.res2)/N
  sig23 = sum(m.res2*y.res1*x.res2*y.res2)/N
  
  Sig = matrix(c(sig11,sig12,sig13,
                 sig12,sig22,sig23,
                 sig13,sig23,sig33),ncol=3,byrow = T)#/N
  
  dU1dbeta = -c(sum(X*x.res1),0,0)/N
  dU2dbeta = -c(sum(X*y.res1),sum(M*m.res2),sum(X*m.res2))/N
  dU3dbeta = -c(0,sum(M*x.res2),sum(X*x.res2))/N
  
  p0 = rep.int(0,p)
  dU1dgamma = -c(V1,p0,V2,p0,p0,p0)
  dU2dgamma = -c(p0,p0,p0,V3,V4,p0)
  dU3dgamma = -c(p0,V5,p0,p0,p0,V6)
  

  dUdbeta = matrix(c(dU1dbeta,dU2dbeta,dU3dbeta),nrow=3,byrow = TRUE)
  dUdgamma = matrix(c(dU1dgamma,dU2dgamma,dU3dgamma),nrow=3,byrow = TRUE)
  
  dUdtheta = cbind(dUdbeta,dUdgamma)
  
  dsig11dbeta = -2*c(sum(X*x.res1^2*m.res1),0,0)/N
  dsig22dbeta = -2*c(sum(X*y.res1^2*m.res2),sum(M*m.res2^2*y.res1),sum(X*m.res2^2*y.res1))/N
  dsig33dbeta = -2*c(0,sum(M*x.res2^2*y.res2),sum(X*x.res2^2*y.res2))/N
  
  dsig12dbeta = c(sum(-X*x.res1*y.res1*(m.res1+m.res2)),
                  sum(-M*x.res1*m.res1*m.res2),
                  sum(-X*x.res1*m.res1*m.res2))/N
  dsig13dbeta = c(sum(-X*x.res1*x.res2*y.res2),
                  sum(-M*x.res1*x.res2*m.res1),
                  sum(-X*x.res1*x.res2*m.res1))/N
  dsig23dbeta = c(sum(-X*x.res2*y.res1*y.res2),
                  sum(-M*x.res2*m.res2*(y.res1+y.res2)),
                  sum(-X*x.res2*m.res2*(y.res1+y.res2)))/N
  
  dsig11dgamx1 = t(Z.matrix) %*% (m.res1^2*hdot1*x.res1)/N
  dsig11dgamm1 = t(Z.matrix) %*% (x.res1^2*m.res1)/N
  dsig22dgamm2 = t(Z.matrix) %*% (y.res1^2*m.res2)/N
  dsig22dgamy1 = t(Z.matrix) %*% (m.res2^2*y.res1)/N
  dsig33dgamx2 = t(Z.matrix) %*% (y.res2^2*x.res2*hdot2)/N
  dsig33dgamy2 = t(Z.matrix) %*% (x.res2^2*y.res2)/N  
  
  dsig11dgamma = -2*c(dsig11dgamx1,p0,dsig11dgamm1,p0,p0,p0)
  dsig22dgamma= -2*c(p0,p0,p0,dsig22dgamm2,dsig22dgamy1,p0)
  dsig33dgamma = -2*c(p0,dsig33dgamx2,p0,p0,p0,dsig33dgamy2)
  
  
  dsig12dgamma = -c(t(Z.matrix) %*% (hdot1*m.res1*m.res2*y.res1),
                    p0,
                    t(Z.matrix) %*% (x.res1*m.res2*y.res1),
                    t(Z.matrix) %*% (x.res1*m.res1*y.res1),
                    t(Z.matrix) %*% (x.res1*m.res1*m.res2),
                    p0)/N
  dsig13dgamma = -c(t(Z.matrix) %*% (hdot1*m.res1*y.res2*x.res2),
                    t(Z.matrix) %*% (hdot2*m.res1*y.res2*x.res1),
                    t(Z.matrix) %*% (x.res1*y.res2*x.res2),
                    p0,
                    p0,
                    t(Z.matrix) %*% (m.res1*x.res2*x.res1))/N
  dsig23dgamma = -c(p0,
                    t(Z.matrix) %*% (hdot2*y.res1*m.res2*y.res2),
                    p0,
                    t(Z.matrix) %*% (x.res2*y.res1*y.res2),
                    t(Z.matrix) %*% (m.res2*x.res2*y.res2),
                    t(Z.matrix) %*% (m.res2*x.res2*y.res1))/N
  
  dsig11dtheta = c(dsig11dbeta,dsig11dgamma)
  dsig22dtheta = c(dsig22dbeta,dsig22dgamma)
  dsig33dtheta = c(dsig33dbeta,dsig33dgamma)
  dsig12dtheta = c(dsig12dbeta,dsig12dgamma)
  dsig13dtheta = c(dsig13dbeta,dsig13dgamma)
  dsig23dtheta = c(dsig23dbeta,dsig23dgamma)
  
  A = solve(Sig)
  AU = A%*%U
  B = dUdtheta - 
    AU[1]*rbind(dsig11dtheta,dsig12dtheta,dsig13dtheta) -
    AU[2]*rbind(dsig12dtheta,dsig22dtheta,dsig23dtheta) -
    AU[3]*rbind(dsig13dtheta,dsig23dtheta,dsig33dtheta)  
  
  f_out = N*(t(AU)%*%U)
  dfdbeta = N*(t(AU)%*%(dUdbeta + B[,1:3]))
  
  
  #### Lots of second derivatives ####
  MX = sum(M*X)
  XX = sum(X*X)
  XZ = t(Z.matrix) %*% (X)
  MZ = t(Z.matrix) %*% (M)
  
  p000 = c(0,0,0,p0,p0,p0,p0,p0,p0)
  
  dU1dbetadtheta = matrix(c(0,0,0,t(Z.matrix) %*% (hdot1*X),p0,p0,p0,p0,p0,
                            p000,p000),byrow=T,nrow=3)/N
  dU2dbetadtheta = matrix(c(0,MX,XX,p0,p0,p0,p0,XZ ,p0,
                            MX,0,0,p0,p0,p0,MZ ,p0,p0,
                            XX,0,0,p0,p0,p0,XZ ,p0,p0),byrow=T,nrow=3)/N
  dU3dbetadtheta = matrix(c( p000,
                             0,0,0,p0,t(Z.matrix) %*% (hdot2*M) ,p0,p0,p0,p0,
                             0,0,0,p0,t(Z.matrix) %*% (hdot2*X) ,p0,p0,p0,p0),byrow=T,nrow=3)/N
  
  
  dsig11db1db1 = 2*sum(X^2* x.res1^2)
  dsig22db1db1 = 2*sum(X^2* y.res1^2)
  dsig22db1db2 = 4*sum(X*M* y.res1*m.res2)
  dsig22db1db3 = 4*sum(X^2* y.res1*m.res2)
  dsig22db2db2 = 2*sum(M^2* m.res2^2)
  dsig22db2db3 = 2*sum(X*M* m.res2^2)
  dsig22db3db3 = 2*sum(X^2* m.res2^2)
  dsig33db2db2 = 2*sum(M^2* x.res2^2)
  dsig33db2db3 = 2*sum(X*M* x.res2^2)
  dsig33db3db3 = 2*sum(X^2* x.res2^2)
  
  dsig11db1dgamx1 = 4*t(Z.matrix) %*% (hdot1*X*x.res1*m.res1)
  dsig11db1dgamm1 = 2*t(Z.matrix) %*% (X*x.res1^2)
  dsig22db1dgamm2 = 2*t(Z.matrix) %*% (X*y.res1^2) 
  dsig22db2dgamm2 = 4*t(Z.matrix) %*% (M*y.res1*m.res2)
  dsig22db3dgamm2 = 4*t(Z.matrix) %*% (X*y.res1*m.res2)
  dsig22db1dgamy1 = 4*t(Z.matrix) %*% (X*y.res1*m.res2)
  dsig22db2dgamy1 = 2*t(Z.matrix) %*% (M*m.res2^2)
  dsig22db3dgamy1 = 2*t(Z.matrix) %*% (X*m.res2^2)
  dsig33db2dgamx2 = 4*t(Z.matrix) %*% (M*y.res2*x.res2*hdot2)
  dsig33db3dgamx2 = 4*t(Z.matrix) %*% (X*y.res2*x.res2*hdot2)
  dsig33db2dgamy2 = 2*t(Z.matrix) %*% (M*x.res2^2)
  dsig33db3dgamy2 = 2*t(Z.matrix) %*% (X*x.res2^2)
  
  dsig11dbetadtheta = matrix(c( dsig11db1db1,0,0,dsig11db1dgamx1, p0,dsig11db1dgamm1,p0,p0,p0,
                                p000,p000),byrow=T,nrow=3)/N
  
  dsig22dbetadtheta = matrix(c(dsig22db1db1,dsig22db1db2,dsig22db1db3,p0,p0,p0,dsig22db1dgamm2,dsig22db1dgamy1,p0,
                               dsig22db1db2,dsig22db2db2,dsig22db2db3,p0,p0,p0,dsig22db2dgamm2,dsig22db2dgamy1,p0,
                               dsig22db1db3,dsig22db2db3,dsig22db3db3,p0,p0,p0,dsig22db3dgamm2,dsig22db3dgamy1,p0),byrow=T,nrow=3)/N
  
  dsig33dbetadtheta = matrix(c(p000,
                               0,dsig33db2db2,dsig33db2db3,p0,dsig33db2dgamx2,p0,p0,p0,dsig33db2dgamy2,
                               0,dsig33db2db3,dsig33db3db3,p0,dsig33db3dgamx2,p0,p0,p0,dsig33db3dgamy2),byrow=T,nrow=3)/N
  
  dsig12db1db1 = 2*sum(X^2* x.res1*y.res1)
  dsig12db1db2 = sum(X*M*   x.res1*(m.res1+m.res2))
  dsig12db1db3 = sum(X^2*   x.res1*(m.res1+m.res2))
  dsig13db1db2 = sum(X*M*   x.res1*x.res2)
  dsig13db1db3 = sum(X^2*   x.res1*x.res2)
  dsig23db1db2 = sum(X*M*     x.res2*(y.res1+y.res2))
  dsig23db1db3 = sum(X^2*     x.res2*(y.res1+y.res2))
  dsig23db2db2 = 2*sum(M^2*   x.res2*m.res2)
  dsig23db2db3 = 2*sum(X*M*   x.res2*m.res2)
  dsig23db3db3 = 2*sum(X^2*   x.res2*m.res2)
  
  dsig12db1dgamx1 = t(Z.matrix) %*% (hdot1*X*y.res1*(m.res1+m.res2))
  dsig12db1dgamy1 = t(Z.matrix) %*% (X*x.res1*(m.res1+m.res2))
  dsig12db1dgamm1 = t(Z.matrix) %*% (X*x.res1*y.res1)
  dsig12db1dgamm2 = dsig12db1dgamm1
  dsig12db2dgamx1 = t(Z.matrix) %*% (hdot1*M*m.res1*m.res2)
  dsig12db2dgamm1 = t(Z.matrix) %*% (M*x.res1*m.res2)
  dsig12db2dgamm2 = t(Z.matrix) %*% (M*x.res1*m.res1)
  dsig12db3dgamx1 = t(Z.matrix) %*% (hdot1*X*m.res1*m.res2)
  dsig12db3dgamm1 = t(Z.matrix) %*% (X*x.res1*m.res2)
  dsig12db3dgamm2 = t(Z.matrix) %*% (X*x.res1*m.res1)
  
  dsig13db1dgamx1 = t(Z.matrix) %*% (hdot1*X*y.res2*x.res2)
  dsig13db1dgamx2 = t(Z.matrix) %*% (X*hdot2*y.res2*x.res1)
  dsig13db1dgamy2 = t(Z.matrix) %*% (X*x.res1*x.res2)
  dsig13db2dgamx1 = t(Z.matrix) %*% (hdot1*M*x.res2*m.res1)
  dsig13db2dgamx2 = t(Z.matrix) %*% (M*hdot2*m.res1*x.res1)
  dsig13db2dgamm1 = t(Z.matrix) %*% (M*x.res2*x.res1)
  dsig13db3dgamx1 = t(Z.matrix) %*% (hdot1*X*x.res2*m.res1)
  dsig13db3dgamx2 = t(Z.matrix) %*% (X*hdot2*m.res1*x.res1)
  dsig13db3dgamm1 = t(Z.matrix) %*% (X*x.res2*x.res1)
  
  dsig23db1dgamx2 = t(Z.matrix) %*% (X*hdot2*y.res1*y.res2)
  dsig23db1dgamy1 = t(Z.matrix) %*% (X*x.res2*y.res2)
  dsig23db1dgamy2 = t(Z.matrix) %*% (X*x.res2*y.res1)
  dsig23db2dgamx2 = t(Z.matrix) %*% (M*hdot2*m.res2*(y.res1+y.res2))
  dsig23db2dgamm2 = t(Z.matrix) %*% (M*x.res2*(y.res1+y.res2))
  dsig23db2dgamy1 = t(Z.matrix) %*% (M*x.res2*m.res2)
  dsig23db2dgamy2 = dsig23db2dgamy1
  dsig23db3dgamx2 = t(Z.matrix) %*% (X*hdot2*m.res2*(y.res1+y.res2))
  dsig23db3dgamm2 = t(Z.matrix) %*% (X*x.res2*(y.res1+y.res2))
  dsig23db3dgamy1 = t(Z.matrix) %*% (X*x.res2*m.res2)
  dsig23db3dgamy2 = dsig23db3dgamy1
  
  
  
  dsig12dbetadtheta = matrix(c(dsig12db1db1,dsig12db1db2,dsig12db1db3,dsig12db1dgamx1,p0,dsig12db1dgamm1,dsig12db1dgamm2,dsig12db1dgamy1,p0,
                               dsig12db1db2,0,0                      ,dsig12db2dgamx1,p0,dsig12db2dgamm1,dsig12db2dgamm2,p0,p0,
                               dsig12db1db3,0,0                      ,dsig12db3dgamx1,p0,dsig12db3dgamm1,dsig12db3dgamm2,p0,p0),byrow=T,nrow=3)/N
  
  dsig13dbetadtheta = matrix(c(0,dsig13db1db2,dsig13db1db3,dsig13db1dgamx1 ,dsig13db1dgamx2,p0              ,p0,p0,dsig13db1dgamy2 ,
                               dsig13db1db2,0,0           ,dsig13db2dgamx1 ,dsig13db2dgamx2,dsig13db2dgamm1 ,p0,p0,p0,
                               dsig13db1db3,0,0           ,dsig13db3dgamx1 ,dsig13db3dgamx2,dsig13db3dgamm1 ,p0,p0,p0),byrow=T,nrow=3)/N
  
  dsig23dbetadtheta = matrix(c(0           ,dsig23db1db2,dsig23db1db3,p0,dsig23db1dgamx2,p0,p0              ,dsig23db1dgamy1,dsig23db1dgamy2,
                               dsig23db1db2,dsig23db2db2,dsig23db2db3,p0,dsig23db2dgamx2,p0,dsig23db2dgamm2 ,dsig23db2dgamy1,dsig23db2dgamy2,
                               dsig23db1db3,dsig23db2db3,dsig23db3db3,p0,dsig23db3dgamx2,p0,dsig23db3dgamm2 ,dsig23db3dgamy1,dsig23db3dgamy2),byrow=T,nrow=3)/N
  
  
  
  d2fdbetadgamma = N*(
    2*t(B[,1:3])%*%A%*%B + 
    2* AU[1]*dU1dbetadtheta + 
    2* AU[2]*dU2dbetadtheta + 
    2* AU[3]*dU3dbetadtheta -
    AU[1]^2 *dsig11dbetadtheta - 
    AU[2]^2 *dsig22dbetadtheta - 
    AU[3]^2 *dsig33dbetadtheta - 
    2*AU[1]*AU[2]*dsig12dbetadtheta - 
    2*AU[1]*AU[3]*dsig13dbetadtheta - 
    2*AU[2]*AU[3]*dsig23dbetadtheta)


  #### nuisance parameters ####
  dV1dbeta = - t(Z.matrix) %*%cbind((X*hdot1),0,0)
  dV2dbeta = -                cbind(p0       ,0,0)
  dV3dbeta = - t(Z.matrix) %*%cbind(0        ,M,X)
  dV4dbeta = - t(Z.matrix) %*%cbind(X        ,0,0)
  dV5dbeta = - t(Z.matrix) %*%cbind(0        ,(M*hdot2),(X*hdot2))
  dV6dbeta = -                cbind(p0       ,0,0)
  
  dVdbeta = matrix(rbind(dV1dbeta,dV2dbeta,dV3dbeta,dV4dbeta,dV5dbeta,dV6dbeta),ncol=3)/N
  
  pp0 = matrix(0,nrow=p,ncol=p)
  
  Z.sq = t(Z.matrix) %*%Z.matrix
  Z.hdot1 = t(Z.matrix) %*%(Z.matrix*c(hdot1))
  Z.hdot2 = t(Z.matrix) %*%(Z.matrix*c(hdot1))
  
  Z.hdotdot1 = -t(Z.matrix) %*% (Z.matrix*c(-((hdot1*(e1-1)/(e1+1))*m.res1)))
  Z.hdotdot2 = -t(Z.matrix) %*% (Z.matrix*c(-((hdot2*(e2-1)/(e2+1))*y.res2)))
  
  dV1dgamma = -cbind(Z.hdotdot1,pp0,Z.hdot1,pp0,pp0,pp0)
  dV2dgamma = -cbind(Z.hdot1,pp0,pp0,pp0,pp0,pp0)
  dV3dgamma = -cbind(pp0,pp0,pp0,pp0,Z.sq,pp0)
  dV4dgamma = -cbind(pp0,pp0,pp0,Z.sq,pp0,pp0)
  dV5dgamma = -cbind(pp0,Z.hdotdot2,pp0,pp0,pp0,Z.hdot2)
  dV6dgamma = -cbind(pp0,Z.hdot2,pp0,pp0,pp0,pp0)

  dVdgamma = matrix(rbind(dV1dgamma,dV2dgamma,dV3dgamma,dV4dgamma,dV5dgamma,dV6dgamma),ncol=6*p)/N
  
  
  #### Return output ####
  vec = c(dfdbeta,V1,V2,V3,V4,V5,V6)
  J = rbind(d2fdbetadgamma,cbind(dVdbeta,dVdgamma))
  
  vec_unconstr = c(U,V1,V2,V3,V4,V5,V6)
  J_unconstr = rbind(cbind(dUdbeta,dUdgamma),cbind(dVdbeta,dVdgamma))
  

  inv_dUdbeta = solve(dUdbeta)
  #PhiPhiT = t(inv_dUdbeta)%*%Sig%*%inv_dUdbeta
  PhiPhiT = inv_dUdbeta%*%Sig%*%t(inv_dUdbeta)
  
  T1 = N*(beta1^2)/PhiPhiT[1,1]
  T2 = N*(beta2^2)/PhiPhiT[2,2]
  T3 = N*(delta^2)/PhiPhiT[3,3]
  rho = PhiPhiT[1,2]/sqrt(PhiPhiT[1,1]*PhiPhiT[2,2])
  
  #RSTest = N*(beta1*beta2)^2/(beta1^2*sig22 + beta2^2*sig11 + 2*beta1*beta2*sig12)
  #RSTest = T1*T2/(T1+T2)
  RSTest = T1*T2/(T1+T2+2*rho*sqrt(T1*T2))
  
  if(brow) {browser()}
  return(list(vec=vec,J=J,score=f_out,
              vec_unconstr=vec_unconstr,J_unconstr=J_unconstr,
              A.mat=A, var=diag(Sig)/N ,
              RSTest = RSTest,
              T_stats = c(T1,T2,T3) ))
  
  #v1 = c(sig11,sig22,sig33,sig12,sig13,sig23)
  #J1 = rbind(dsig11dtheta,dsig22dtheta,dsig33dtheta,dsig12dtheta,dsig13dtheta,dsig23dtheta)
  #v1 = vec
  #J1 = J
  #return(list(vec=v1,J=J1))
  
}

Lagrangian_cue_cop <-function(formula,df,theta,med_prop=NULL){
  if(!is.null(med_prop)){#if there is a constraint do this:
    d = length(theta)
    lambda = theta[d]
    step = CUE_vec_J_cop(formula,df,theta[-d])
    vec = c(step$vec,0)
    J = cbind(rbind(step$J,0),0)
    
    #constraint g=0
    g = theta[1]*theta[2]*(1-med_prop)-med_prop*theta[3]
    dgdb = c(theta[2]*(1-med_prop),theta[1]*(1-med_prop),-med_prop)
    d2gdb2 = matrix(c(0,(1-med_prop),0,
                      (1-med_prop),0,0,
                      0,0,0),nrow=3)
    
    vec[1:3] <- vec[1:3] +lambda*dgdb
    vec[d] <- g
    
    J[1:3,1:3] <- J[1:3,1:3] + lambda*d2gdb2
    J[1:3,d] <- dgdb
    J[d,1:3] <- dgdb
    
    
  }else{
    step = CUE_vec_J_cop(formula,df,theta)
    vec = step$vec
    J = step$J
  }
  
  return(list(vec=vec,J=J))
}




### 2step code for cop


derivs_2step_cop <- function(formula,df,beta,gamma,A){
  beta1 = beta[1]
  beta2 = beta[2]
  delta  = beta[3]
  
  p = length(gamma)/6
  gamma_x1 = gamma[1:p]
  gamma_x2 = gamma[(p+1)  :(2*p)]
  gamma_m1 = gamma[(2*p+1):(3*p)]
  gamma_m2 = gamma[(3*p+1):(4*p)]
  gamma_y1 = gamma[(4*p+1):(5*p)]
  gamma_y2 = gamma[(5*p+1):(6*p)]
  
  X = df$X
  M = df$M
  Y = df$Y
  N = length(M)
  Z.matrix = model.matrix(as.formula(formula),data = df)
  
  odds_x1  = Z.matrix%*%gamma_x1 #The log-odds
  odds_x2  = Z.matrix%*%gamma_x2
  f1       = Z.matrix%*%gamma_m1
  f2       = Z.matrix%*%gamma_m2
  g1       = Z.matrix%*%gamma_y1 
  g2       = Z.matrix%*%gamma_y2 
  
  
  e1 = exp(odds_x1)
  e2 = exp(odds_x2)
  nt = 1 # Number of trials for X~ Binomial (nt,p(Z))
  
  x.res1   = X-nt*e1/(1+e1)
  x.res2   = X-nt*e2/(1+e2)
  m.res1   = M - beta1*X - f1 
  m.res2   = M - beta1*X - f2 
  y.res1   = Y - beta2*M - delta*X - g1
  y.res2   = Y - beta2*M - delta*X - g2
  
  
  hdot1     = nt*e1/(1+e1)^2  #derivative of E[X|Z] wrt. log-odds (=odds_x1)
  hdot2     = nt*e2/(1+e2)^2
  
  U = c(sum(x.res1*m.res1)/N,sum(m.res2*y.res1)/N,sum(x.res2*y.res2)/N)
  
  dU1dbeta = -c(sum(X*x.res1),0,0)/N
  dU2dbeta = -c(sum(X*y.res1),sum(M*m.res2),sum(X*m.res2))/N
  dU3dbeta = -c(0,sum(M*x.res2),sum(X*x.res2))/N
  dUdbeta = matrix(c(dU1dbeta,dU2dbeta,dU3dbeta),nrow=3,byrow = TRUE)
  
  MX = sum(M*X)
  XX = sum(X*X)
  dU2dbetadbeta = matrix(c(0,MX,XX,MX,0,0,XX,0,0),nrow=3)/N
  AU = A%*%U
  
  return(list(score = N*(t(AU)%*%U),
              vec = N*2*(t(AU)%*%dUdbeta),
              J   = N*(2*t(dUdbeta)%*%A%*%dUdbeta + 2*AU[2]*dU2dbetadbeta)))
}


Lagrangian_2step_cop <- function(formula,df,beta,med_prop=NULL,...){
  lambda = beta[4]
  step = derivs_2step_cop(formula,df,beta[1:3],...)
  vec = c(step$vec,0)
  J = cbind(rbind(step$J,0),0)
  #constraint g=0
  g = beta[1]*beta[2]*(1-med_prop)-med_prop*beta[3]
  dgdb = c(beta[2]*(1-med_prop),beta[1]*(1-med_prop),-med_prop)
  d2gdb2 = matrix(c(0,(1-med_prop),0,
                    (1-med_prop),0,0,
                    0,0,0),nrow=3)
  vec[1:3] <- vec[1:3] +lambda*dgdb
  vec[4] <- g
  J[1:3,1:3] <- J[1:3,1:3] + lambda*d2gdb2
  J[1:3,4] <- dgdb
  J[4,1:3] <- dgdb
  return(list(vec = vec, J = J,score=step$score))
}


