rm(list=ls())

library(rootSolve)

LVpeypred2=function(t,y,parameters)
{
  dy1=parameters[1]*y[1]*(1-(y[1]/parameters[5]))-parameters[2]*y[1]*y[2] #Proie=N=x1 et prédateurs=P=x2 ET parms1=r ,parms2=a, parms3=e et parms4=m
  dy2=parameters[2]*parameters[3]*y[1]*y[2]-parameters[4]*y[2]
  list(c(dy1,dy2))
}

Equi<-function(r,a,e,m,K){
  p1<-c(0,0)
  p2<-c(K,0) 
  p3<-c((m/(e*a)),((r*K-(r*m)/(a*e))/(a*K)))
  rbind(p1,p2,p3)
}


r=0.2
a=0.2
e=1
m=0.4
K=10

Equi(r,a,e,m,K)

stability<-function(r,a,e,m,K){
  Matrix=Equi(r,a,e,m,K)
  Eigen<-vector(length=nrow(Matrix))
  for (i in 1:nrow(Matrix)){
    Jac<-jacobian.full(Matrix[i,], func=LVpeypred2, dy = NULL, time = 0, parms = c(r,a,e,m,K), 
                  pert = 1e-8)
    Tr<-Jac[1,1]+Jac[2,2]
    Det<-Jac[1,1]*Jac[2,2]-Jac[1,2]*Jac[2,1]
    delta<-Tr^(2)-4*Det
    if (Det<0){stab<-"Point selle"}
      else{if(Tr<0){stab<- -1}
      else{stab<- 1}}
    Eigen[i]<-stab
  }
  return(cbind(Matrix,Eigen))
}

stability(r,a,e,m,K)
