rm(list=ls())

library(deSolve)
library(phaseR)


LVpeypred1=function(t,x,parms)
{
  dx1=parms[1]*x[1]-parms[2]*x[1]*x[2] #Proie=N=x1 et prédateurs=P=x2 ET parms1=r ,parms2=a, parms3=e et parms4=m
  dx2=parms[2]*parms[3]*x[1]*x[2]-parms[4]*x[2]
  list(c(dx1,dx2))
}

LVpeypred2=function(t,y,parameters)
{
  dy1=parameters[1]*y[1]-parameters[2]*y[1]*y[2] #Proie=N=x1 et prédateurs=P=x2 ET parms1=r ,parms2=a, parms3=e et parms4=m
  dy2=parameters[2]*parameters[3]*y[1]*y[2]-parameters[4]*y[2]
  list(c(dy1,dy2))
}

temps=seq(0,100,by=0.1)

init1=c(10,20)
init2=c(20,10)

r1=0.8
a1=0.2
e1=0.2
m1=0.4
solution1=lsoda(y=init1,times=temps,func=LVpeypred1,parms=c(r1,a1,e1,m1))

solution2=lsoda(y=init2,times=temps,func=LVpeypred1,parms=c(r1,a1,e1,m1))

plot(solution1[,2],solution1[,3],col=4,type='l',lwd=2,ylim=c(0,30),xlim=c(0,100),
     ylab="Densité de population P",
     xlab='Densité de pop de N',
     main="Portrait de phase du modèle LV")
lines(solution2[,2],solution2[,3],lty=1,col=2,lwd=2)
abline(v=m1/(a1*e1),lty=2,col=1,lwd=1)
abline(h=r1/a1,lty=1,col=1,lwd=1)
points(c(10,20),c(20,10),col=c(4,2))
legend('topright',legend=c("start: 10prey,20pred", "start: 20prey,10pred"),lwd=2,col=c(4,2),lty=c(1,1) )

flowField(deriv=LVpeypred2, x.lim=c(0,100),y.lim=c(0,30),parameters=c(r1,a1,e1,m1),points=20)
nullclines(deriv=LVpeypred2, x.lim=c(0,100),y.lim=c(0,30), parameters =c(r1,a1,e1,m1), points = 200,
           system = "two.dim", colour = c("red", "blue"), add = TRUE)
trajectory(deriv=LVpeypred2, y0 =c(60,15), n = NULL, t.start = 0, t.end=100, t.step = 0.01,
           parameters = c(r1,a1,e1,m1), system = "two.dim", colour = "purple",lwd=3,
           add = TRUE)

plot(temps,solution1[,2],type='l',col="blue",ylab="N et P",
     xlab='temps',
     main="Chroniques du modèle LV")
lines(temps,solution1[,3],lty=2,col="blue")
lines(temps,solution2[,2],lty=1,col="red")
lines(temps,solution2[,3],lty=2,col="red")
legend('topright',legend=c("Prey-10prey,20pred",
                           "Pred-10prey,20pred", "Prey-20prey,10pred", "Pred-20prey,10pred"),
       lwd=2,col=c(4,4,2,2),lty=c(1,2,1,2) )
