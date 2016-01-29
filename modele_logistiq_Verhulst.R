rm(list=ls())

Eq2=function(t,x,parms)
{
  dx1=parms[1]*x[1]*(1-(x[1]/parms[2]))
  list(dx1)
}

temps=seq(0,100,by=0.1)

init1=c(50)
init2=c(400)

r1=0.5
K1=300
solution1=lsoda(y=init1,times=temps,func=Eq2,parms=c(r1,K1))

r2=1
K2=300
solution2=lsoda(y=init1,times=temps,func=Eq2,parms=c(r2,K2))

r3=2
K3=300
solution3=lsoda(y=init1,times=temps,func=Eq2,parms=c(r3,K3))

r4=0.5
K4=300
solution4=lsoda(y=init2,times=temps,func=Eq2,parms=c(r4,K4))

r5=0.5
K5=400
solution5=lsoda(y=init2,times=temps,func=Eq2,parms=c(r5,K5))

r6=0.5
K6=400
solution6=lsoda(y=init1,times=temps,func=Eq2,parms=c(r6,K6))

r7=-0.5
K7=300
solution7=lsoda(y=init1,times=temps,func=Eq2,parms=c(r7,K7))

plot(temps,solution1[,2],ylim=c(0,500),type="l",
     ylab="Densité de population (N)",
     main="Modèle de Verhuslt")
lines(temps,solution2[,2],lwd=2)
lines(temps,solution3[,2],lwd=2,lty=2)
lines(temps,solution4[,2],lwd=1,lty=3)
lines(temps,solution5[,2],lwd=2,lty=3)
lines(temps,solution6[,2],lwd=1,lty=5)
lines(temps,solution7[,2],lty=2)
legend(15,270,legend=c("NO=50 ; r=0.5 ; K=300","NO=50 ; r=1 ; K=300",
                       "NO=50 ; r=2 ; K=300","NO=400 ; r=0.5 ; K=300",
                       "NO=400 ; r=0.5 ; K=400","NO=50 ; r=0.5 ; K=400",
                       "NO=50 ; r=-0.5 ; K=300"),lty=c(1,1,2,3,3,5,2),lwd=c(1,2,2,1,2,1,1))