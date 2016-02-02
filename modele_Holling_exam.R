rm(list=ls())

library(rootSolve)
library(deSolve)
library(phaseR)

Holling1=function(t,x,parms)
{
  dx1=x[1]*(2-x[1])-(x[1]*x[2])/(1+x[1]) #parms1=mu 
  dx2=(x[1]*x[2])/(1+x[1])-parms[1]*x[2]
  list(c(dx1,dx2))
}

Holling2=function(t,y,parameters)
{
  dy1=y[1]*(2-y[1])-(y[1]*y[2])/(1+y[1]) #parms1=mu 
  dy2=(y[1]*y[2])/(1+y[1])-parameters[1]*y[2]
  list(c(dy1,dy2))
}

ptEqHolling<-function(mu){
  p1<-c(0,0)
  p2<-c(2,0) 
  p3<-c(mu/(1-mu),(2-3*mu)/(1-mu)^2)
  if(mu>0 && mu<(2/3)){eq<-rbind(p1,p2,p3)}
  return(eq)
}

mu=2.89/6

Eq<-ptEqHolling(mu)

HollingStability<-function(mu){
  Eigen<-vector(length=nrow(Eq))
  Delta<-vector(length=nrow(Eq))
  for (i in 1:nrow(Eq)){
    Jac<-jacobian.full(Eq[i,], func=Holling1, time = 0, parms=c(mu))
    Tr<-Jac[1,1]+Jac[2,2]
    Det<-Jac[1,1]*Jac[2,2]-Jac[1,2]*Jac[2,1]
    delta<-Tr^(2)-4*Det
    if (Det<0 && Tr!=0){stab<-"PS"}  #Point Selle
    if (Tr==0){stab<-"C"}   #Centres
    if(delta>0 && Det>0 && Tr<0){stab<- "NAS"}  #Noeud Asymptotiquement Stable
    if(delta<0 && Det>0 && Tr<0){stab<- "FAS"}  #Foyer Asymptotiquement Stable
    if (delta>0 && Det>0 && Tr>0){stab<- "NI"}  #Noeud Instable
    if (delta<0 && Det>0 && Tr>0){stab<- "FI"}  #Foyer Instable
    Eigen[i]<-stab
    Delta[i]<-delta
    print(Jac)
    #print(Det)
    #print(Tr)
    #print(delta)
  }
  return(cbind(Eq,Eigen,Delta))#& Delta ?
}

HollingStability(mu)



temps=seq(0,100,by=0.1)
init1=c(1,10)
init2=c(20,5)

solution1=lsoda(y=init1,times=temps,func=Holling1,parms=c(mu))
solution2=lsoda(y=init2,times=temps,func=Holling1,parms=c(mu))
plot(solution1[,2],solution1[,3],col=4,type='l',lwd=1,xlim=c(0,20),ylim=c(0,10),
     ylab="Prédateurs",
     xlab='Proies',
     main="Portrait de phase du modèle Holling")
lines(solution2[,2],solution2[,3],lty=2,col=2,lwd=2)
legend('topright',legend=c("1prey,10pred", "20prey,5pred"),lwd=2,col=c(4,2),lty=c(1,2) )

flowField(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30),parameters=c(mu),points=25)
"nullclines(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30), parameters =c(mu), points = 25,
system =\"two.dim\", colour = c(\"red\", \"blue\"), add = TRUE)"
trajectory(deriv=Holling1, y0 =c(30,15), n = NULL, t.start = 0, t.end=200, t.step = 0.01,
           parameters = c(mu), system = "two.dim", colour = "purple",lwd=1,lty=2,
           add = TRUE)

points(c(init1[1],init2[1]),c(init1[2],init2[2]),col=c(4,2),lwd=3)

for (i in 1:nrow(Eq)){
  points(Eq[i,1],Eq[i,2],col=1,lwd=3)
}



x11()
plot(temps,solution1[,2],type='l',col="blue",ylab="N et P",
     xlab='temps',ylim=c(0,20),xlim=c(0,100),
     main="Chroniques du modèle Holling")
lines(temps,solution1[,3],lty=2,col="blue")
lines(temps,solution2[,2],lty=1,col="red")
lines(temps,solution2[,3],lty=2,col="red")
legend('topright',legend=c("Prey-1prey,10pred",
                           "Pred-1prey,10pred", "Prey-20prey,5pred", "Pred-20prey,5pred"),
       lwd=2,col=c(4,4,2,2),lty=c(1,2,1,2) )

#------------------------------------------------For various mu--------------------------------------------------
mu2=1/6

Eq<-ptEqHolling(mu2)

HollingStability(mu2)

temps=seq(0,100,by=0.1)
init1=c(1,10)
init2=c(20,5)

x11()
solution1=lsoda(y=init1,times=temps,func=Holling1,parms=c(mu2))
solution2=lsoda(y=init2,times=temps,func=Holling1,parms=c(mu2))
plot(solution1[,2],solution1[,3],col=4,type='l',lwd=1,xlim=c(0,20),ylim=c(0,10),
     ylab="Prédateurs",
     xlab='Proies',
     main="Portrait de phase du modèle Holling")
lines(solution2[,2],solution2[,3],lty=2,col=2,lwd=2)
legend('topright',legend=c("1prey,10pred", "20prey,5pred"),lwd=2,col=c(4,2),lty=c(1,2) )

flowField(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30),parameters=c(mu2),points=25)
"nullclines(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30), parameters =c(mu2), points = 25,
system =\"two.dim\", colour = c(\"red\", \"blue\"), add = TRUE)"
trajectory(deriv=Holling1, y0 =c(30,15), n = NULL, t.start = 0, t.end=200, t.step = 0.01,
           parameters = c(mu2), system = "two.dim", colour = "purple",lwd=1,lty=2,
           add = TRUE)

points(c(init1[1],init2[1]),c(init1[2],init2[2]),col=c(4,2),lwd=3)

for (i in 1:nrow(Eq)){
  points(Eq[i,1],Eq[i,2],col=1,lwd=3)
}



x11()
plot(temps,solution1[,2],type='l',col="blue",ylab="N et P",
     xlab='temps',ylim=c(0,20),xlim=c(0,100),
     main="Chroniques du modèle Holling")
lines(temps,solution1[,3],lty=2,col="blue")
lines(temps,solution2[,2],lty=1,col="red")
lines(temps,solution2[,3],lty=2,col="red")
legend('topright',legend=c("Prey-1prey,10pred",
                           "Pred-1prey,10pred", "Prey-20prey,5pred", "Pred-20prey,5pred"),
       lwd=2,col=c(4,4,2,2),lty=c(1,2,1,2) )


#------------------------------------------------For various mu--------------------------------------------------
mu3=3.5/6

Eq<-ptEqHolling(mu3)

HollingStability(mu3)

temps=seq(0,200,by=0.1)
init1=c(1,10)
init2=c(20,5)

x11()
solution1=lsoda(y=init1,times=temps,func=Holling1,parms=c(mu3))
solution2=lsoda(y=init2,times=temps,func=Holling1,parms=c(mu3))
plot(solution1[,2],solution1[,3],col=4,type='l',lwd=1,xlim=c(0,20),ylim=c(0,10),
     ylab="Prédateurs",
     xlab='Proies',
     main="Portrait de phase du modèle Holling")
lines(solution2[,2],solution2[,3],lty=2,col=2,lwd=2)
legend('topright',legend=c("1prey,10pred", "20prey,5pred"),lwd=2,col=c(4,2),lty=c(1,2) )

flowField(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30),parameters=c(mu3),points=25)
"nullclines(deriv=Holling2, x.lim=c(0,25),y.lim=c(0,30), parameters =c(mu3), points = 25,
system =\"two.dim\", colour = c(\"red\", \"blue\"), add = TRUE)"
trajectory(deriv=Holling1, y0 =c(30,15), n = NULL, t.start = 0, t.end=200, t.step = 0.01,
           parameters = c(mu3), system = "two.dim", colour = "purple",lwd=1,lty=2,
           add = TRUE)

points(c(init1[1],init2[1]),c(init1[2],init2[2]),col=c(4,2),lwd=3)

for (i in 1:nrow(Eq)){
  points(Eq[i,1],Eq[i,2],col=1,lwd=3)
}



x11()
plot(temps,solution1[,2],type='l',col="blue",ylab="N et P",
     xlab='temps',ylim=c(0,20),xlim=c(0,100),
     main="Chroniques du modèle Holling")
lines(temps,solution1[,3],lty=2,col="blue")
lines(temps,solution2[,2],lty=1,col="red")
lines(temps,solution2[,3],lty=2,col="red")
legend('topright',legend=c("Prey-1prey,10pred",
                           "Pred-1prey,10pred", "Prey-20prey,5pred", "Pred-20prey,5pred"),
       lwd=2,col=c(4,4,2,2),lty=c(1,2,1,2) )


#------------------------------------------------For various mu--------------------------------------------------
mu4=0.0000001/6

Eq<-ptEqHolling(mu4)

HollingStability(mu4)
