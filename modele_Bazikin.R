rm(list=ls())

library(rootSolve)
library(deSolve)
library(phaseR)

Bazikin1=function(t,x,parms)
{
  dx1=x[1]-(x[1]*x[2]/(1+parms[2]*x[1]))-parms[1]*x[1]*x[1] #parms1=E ,parms2=a, parms3=gamma 
  dx2=-parms[3]*x[2]*(1-x[1]/(1+parms[2]*x[1]))
  list(c(dx1,dx2))
}

Bazikin2=function(t,y,parameters)
{
  dy1=y[1]-(y[1]*y[2]/(1+parameters[2]*y[1]))-parameters[1]*y[1]*y[1] #parms1=E ,parms2=a, parms3=gamma 
  dy2=-parameters[3]*y[2]*(1-y[1]/(1+parameters[2]*y[1]))
  list(c(dy1,dy2))
}

Bazikinpt<-function(E){
  p1<-c(0,0)
  p2<-c(1/E,0) 
  p3<-c(2,2-4*E)
  if(E==0){eq<-rbind(p1,p3)}
  if(E>0 && E<0.5){eq<-rbind(p1,p2,p3)}
  if(E>=0.5){eq<-rbind(p1,p2)}
  return(eq)
}

E=1/6
a=0.5
gamma=1

Eq<-Bazikinpt(E)

BazikinStability<-function(E,a,gamma){
  Eigen<-vector(length=nrow(Eq))
  Delta<-vector(length=nrow(Eq))
  for (i in 1:nrow(Eq)){
    Jac<-jacobian.full(Eq[i,], func=Bazikin1, time = 0, parms=c(E,a,gamma))
    Tr<-Jac[1,1]+Jac[2,2]
    Det<-Jac[1,1]*Jac[2,2]-Jac[1,2]*Jac[2,1]
    delta<-Tr^(2)-4*Det
    if (Det<0 && Tr!=0){stab<-"PS"}
    if (Tr==0){stab<-"C"}
    if(Det>0 && Tr<0){stab<- "AS"}
    if (Det>0 && Tr>0){stab<- "I"}
    Eigen[i]<-stab
    Delta[i]<-delta
    print(Jac)
    print(Det)
    print(delta)
  }
  return(cbind(Eq,Eigen,Delta))
}

BazikinStability(E,a,gamma)



temps=seq(0,100,by=0.1)
init1=c(1,10)
init2=c(20,5)

solution1=lsoda(y=init1,times=temps,func=Bazikin1,parms=c(E,a,gamma))
solution2=lsoda(y=init2,times=temps,func=Bazikin1,parms=c(E,a,gamma))
plot(solution1[,2],solution1[,3],col=4,type='l',lwd=1,xlim=c(0,20),ylim=c(0,10),
     ylab="Herbivores",
     xlab='Biomasse',
     main="Portrait de phase du modèle Bazikin")
lines(solution2[,2],solution2[,3],lty=2,col=2,lwd=2)
legend('topright',legend=c("1prey,10Herbi", "20prey,5Herbi"),lwd=2,col=c(4,2),lty=c(1,2) )

flowField(deriv=Bazikin2, x.lim=c(0,25),y.lim=c(0,30),parameters=c(E,a,gamma),points=25)
"nullclines(deriv=Bazikin2, x.lim=c(0,25),y.lim=c(0,30), parameters =c(E,a,gamma), points = 25,
           system =\"two.dim\", colour = c(\"red\", \"blue\"), add = TRUE)"
trajectory(deriv=Bazikin1, y0 =c(30,15), n = NULL, t.start = 0, t.end=200, t.step = 0.01,
           parameters = c(E,a,gamma), system = "two.dim", colour = "purple",lwd=1,lty=2,
           add = TRUE)
points(c(1,20),c(10,5),col=c(4,2),lwd=3)

for (i in 1:nrow(Eq)){
  points(Eq[i,1],Eq[i,2],col=1,lwd=3)
}



x11()
plot(temps,solution1[,2],type='l',col="blue",ylab="N et P",
     xlab='temps',ylim=c(0,20),xlim=c(0,100),
     main="Chroniques du modèle Bazikin")
lines(temps,solution1[,3],lty=2,col="blue")
lines(temps,solution2[,2],lty=1,col="red")
lines(temps,solution2[,3],lty=2,col="red")
legend('topright',legend=c("Prey-1prey,10Herbi",
                           "Herbi-1prey,10Herbi", "20prey,5Herbi", "20prey,5Herbi"),
       lwd=2,col=c(4,4,2,2),lty=c(1,2,1,2) )

