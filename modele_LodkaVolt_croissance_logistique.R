rm(list=ls())

library(deSolve)
library(phaseR)

LVpeypred1=function(t,x,parms)
{
  dx1=parms[1]*x[1]*(1-(x[1]/parms[5]))-parms[2]*x[1]*x[2] #Proie=N=x1 et prédateurs=P=x2 ET parms1=r ,parms2=a, parms3=e et parms4=m ,5=k
  dx2=parms[2]*parms[3]*x[1]*x[2]-parms[4]*x[2]
  list(c(dx1,dx2))
}

LVpeypred2=function(t,y,parameters)
{
  dy1=parameters[1]*y[1]*(1-(y[1]/parameters[5]))-parameters[2]*y[1]*y[2] #Proie=N=x1 et prédateurs=P=x2 ET parms1=r ,parms2=a, parms3=e et parms4=m
  dy2=parameters[2]*parameters[3]*y[1]*y[2]-parameters[4]*y[2]
  list(c(dy1,dy2))
}

temps=seq(0,200,by=0.1)

init1=c(10,2)
init2=c(10,10)#ou 1,2

r1=0.2
a1=0.2
e1=1
m1=0.4
K1=10# ou 1
solution1=lsoda(y=init1,times=temps,func=LVpeypred1,parms=c(r1,a1,e1,m1,K1))
solution2=lsoda(y=init2,times=temps,func=LVpeypred1,parms=c(r1,a1,e1,m1,K1))

plot(solution1[,2],solution1[,3],col=4,type='l',lwd=2,xlim=c(0,20),ylim=c(0,10),
     ylab="Densité de population P",
     xlab='Densité de pop de N',
     main="Portrait de phase du modèle LV_logistique")
lines(solution2[,2],solution2[,3],lty=1,col=2,lwd=2)
legend('topright',legend=c("10prey,2pred", "10prey,10pred"),lwd=2,col=c(4,2),lty=c(1,1) )

flowField(deriv=LVpeypred2, x.lim=c(0,200),y.lim=c(0,30),parameters=c(r1,a1,e1,m1,K1),points=100)
nullclines(deriv=LVpeypred2, x.lim=c(0,200),y.lim=c(0,30), parameters =c(r1,a1,e1,m1,K1), points = 200,
           system = "two.dim", colour = c("red", "blue"), add = TRUE)
trajectory(deriv=LVpeypred2, y0 =c(60,15), n = NULL, t.start = 0, t.end=200, t.step = 0.01,
           parameters = c(r1,a1,e1,m1,K1), system = "two.dim", colour = "purple",lwd=3,
           add = TRUE)


plot(temps,solution1[,2],type='l',col="blue",ylab="N et P",
     xlab='temps',ylim=c(0,20),xlim=c(0,200),
     main="Chroniques du modèle LV_logistique")
lines(temps,solution1[,3],lty=2,col="blue")
lines(temps,solution2[,2],lty=1,col="red")
lines(temps,solution2[,3],lty=2,col="red")
legend('topright',legend=c("Prey-10prey,2pred",
                           "Pred-10prey,2pred", "Prey-10prey,10pred", "Pred-10prey,10pred"),
       lwd=2,col=c(4,4,2,2),lty=c(1,2,1,2) )
