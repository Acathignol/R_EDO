rm(list=ls())

Eq1=function(t,x,parms)
{
  dx1=parms[1]*x[1]
  list(dx1)
}

temps=seq(0,20,by=0.1)

init1=c(50)

r1=0.1
solution1=lsoda(y=init1,times=temps,func=Eq1,parms=c(r1))

r2=0
solution2=lsoda(y=init1,times=temps,func=Eq1,parms=c(r2))

r3=-0.1
solution3=lsoda(y=init1,times=temps,func=Eq1,parms=c(r3))

plot(temps,solution1[,2],ylim=c(0,120),type="l",
     ylab="Densité de population (N)",
     main="Modèle de Malthus pour différentes valeurs de r")
lines(temps,solution2[,2],type="l",col=2)
lines(temps,solution3[,2],type="l",col=3)
legend("topleft",legend=c("r=0.1","r=0","r=-0.1"),col=c(1,2,3),lty=1)