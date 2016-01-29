rm(list=ls())

Eq3=function(t,x,parms)
{
  dx1=parms[1]*x[1]*(1-(x[1]/parms[2])^parms[3])
  list(dx1)
}

temps=seq(0,100,by=0.1)

init1=c(50)

r1=0.5
K1=300
theta1=0.5
solution1=lsoda(y=init1,times=temps,func=Eq3,parms=c(r1,K1,theta1))

theta2=1
solution2=lsoda(y=init1,times=temps,func=Eq3,parms=c(r1,K1,theta2))

theta3=2
solution3=lsoda(y=init1,times=temps,func=Eq3,parms=c(r1,K1,theta3))

theta4=10
solution4=lsoda(y=init1,times=temps,func=Eq3,parms=c(r1,K1,theta4))

r5=-0.5
theta5=0.5
solution5=lsoda(y=init1,times=temps,func=Eq3,parms=c(r5,K1,theta5))

plot(temps,solution1[,2],ylim=c(0,500),type="l",
     ylab="Densité de population (N)",
     main="Modèle de Verhuslt")
lines(temps,solution2[,2],lwd=2)
lines(temps,solution3[,2],lwd=2,lty=2)
lines(temps,solution4[,2],lwd=1,lty=3)
lines(temps,solution5[,2],lwd=2,lty=3)
legend(15,270,legend=c("NO=50 ; r=0.5 ; K=300 ; theta=0.5","NO=50 ; r=0.5 ; K=300 ; theta=1",
                       "NO=50 ; r=0.5 ; K=300 ; theta=2","NO=50 ; r=0.5 ; K=300 ; theta=10",
                       "NO=50 ; r=-0.5 ; K=300 ; theta=0.5"),lty=c(1,1,2,3,3),lwd=c(1,2,2,1,2))