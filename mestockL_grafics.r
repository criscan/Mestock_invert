rm(list=ls())
library(areaplot)

system('mod_xxx -ind xxx.dat -nox')  # for model running

source('read.admb.R')
source('read.admbFit.R')
source('por_recluta_r.R')


data <-read.rep('for_R.rep')
attach(data)


target=0.4
#Desembarques--------------------------------------------------------------------------
plot(Yrs, Desembarques[1,], cex.lab=1.5, type="h", xlab="Año",ylab="Miles t", col = "gray",
     ylim = c(0,max(Desembarques[1,])*1.1),lwd=6,
     main="Desembarques", cex.main = 1.5)
lines(Yrs, Desembarques[2,],col="red")

#Indices--------------

par(mfrow = c(2, 2))
sim=20

ubi=which(Bcrucero[1,]>0)
maxi=max(Bcrucero[1:2,ubi])*1.1;
mini=0;


plot(Yrs[ubi],Bcrucero[2,ubi],ylab="Indice",xlab="Año",main="Campaña",pch=sim,type="l",ylim=c(mini,maxi),
     col="red",lwd=2)

lines(Yrs[ubi],Bcrucero[1,ubi],pch=sim,type="b")

#CPUE------------------------------------------------------
ubi=which(CPUE[1,]>0)
maxi=max(CPUE[1:2,ubi])*1.1;
mini=0;


plot(Yrs[ubi],CPUE[2,ubi],ylab="Indice",xlab="Año",main="CPUE",pch=sim,type="l",ylim=c(mini,maxi),
     col="red",lwd=2)

lines(Yrs[ubi],CPUE[1,ubi],pch=sim,type="b")


#Talla promedio------------------------------------------------------

maxi=max(Lmed_flo[2:3,])*1.05
mini=min(Lmed_flo[2:3,])*.95

y=Lmed_flo[2,]
x=Lmed_flo[1,]
plot(x,y,ylab="Talla promedio",xlab="Año",main="Talla promedio flota",
     pch=sim,type="b",ylim=c(mini,maxi))
lines(x,Lmed_flo[3,],lwd=2,col="red")


maxi=max(Lmed_srv[2:3,])*1.05
mini=min(Lmed_srv[2:3,])*.95

y=Lmed_srv[2,]
x=Lmed_srv[1,]
plot(x,y,ylab="Talla promedio",xlab="Año",main="Talla promedio campaña",
     pch=sim,type="b",ylim=c(mini,maxi))
lines(x,Lmed_srv[3,],lwd=2,col="red")



#Comps_tallas_f----------------------------------------------------------------------
n=length(Lmed_flo[1,])
dl=0.5*(Tallas[2]-Tallas[1])

par(mfcol = c(3, 3))
for (i in 1:n)
{
  areaplot(Tallas,Frecs_capt_obs[i,],main=paste(Lmed_flo[1,i]),lwd=0.5, col="lightgray",ylab="", xlab="Talla",
           ylim=c(0,max(c(Frecs_capt_obs[i,],Frecs_capt_pred[i,]))))
  lines(Tallas,Frecs_capt_pred[i,],col="red",lwd=2)
}


#Comps_tallas_s----------------------------------------------------------------------
n=length(Lmed_srv[1,])

par(mfcol = c(4, 3))

for (i in 1:n)
{
  areaplot(Tallas,Frecs_srv_obs[i,],main=paste(Lmed_srv[1,i]),lwd=0.5, col="lightgray",ylab="", xlab="Talla",
           ylim=c(0,max(c(Frecs_srv_obs[i,],Frecs_srv_pred[i,]))))
  lines(Tallas,Frecs_srv_pred[i,],col="red",lwd=2)
}


#Comps_tallas_marginal----------------------------------------------------------------------

par(mfrow = c(1, 2))
areaplot(Tallas, Frec_marg_flo[1,],col="gray",lwd=0.5,xlab="Talla",ylab="Proporcion",main="Frec marginal capturas")
lines(Tallas,Frec_marg_flo[2,],col="red",lwd=2)

areaplot(Tallas, Frec_marg_srv[1,],col="gray",lwd=0.5,xlab="Talla",ylab="Proporcion",main="Frec marginal campañas")
lines(Tallas,Frec_marg_srv[2,],col="red",lwd=2)

#qqplot_tallas------------

par(mfrow=c(2,2))
Residuals <- as.vector(Frecs_capt_obs)-as.vector(Frecs_capt_pred)
Residuals=Residuals/sd(Residuals)
xfit<-seq(min(Residuals),max(Residuals),length=40) 
yfit<-dnorm(xfit) 
hist(Residuals, freq=FALSE, xlab="Residuales std (CompsL capturas)", ylab="Frequencia",30)
lines(xfit, yfit,col="red",lwd=2)
qqnorm(Residuals,xlab="Cuantiles teóricos", ylab="Residuales stds (CompsL capturas)",col="gray")
lines(xfit, xfit,col="red",lwd=2)

Residuals <- as.vector(Frecs_srv_obs)-as.vector(Frecs_srv_pred)
Residuals=Residuals/sd(Residuals)
xfit<-seq(min(Residuals),max(Residuals),length=40) 
yfit<-dnorm(xfit) 
hist(Residuals, freq=FALSE, xlab="Residuales std (CompsL campaña)", ylab="Frequencia",30)
lines(xfit, yfit,col="red",lwd=2)
qqnorm(Residuals,xlab="Cuantiles teóricos", ylab="Residuales stds (CompsL campaña)",col="gray")
lines(xfit, xfit,col="red",lwd=2)



#Selectividad----------------------------------------------------------------------

par(mfrow=c(2,2))

persp(Yrs,Tallas,Sel_f,theta=60,phi=40,expand=0.6, ticktype = "detailed",main="Selectividad capturas",
      xlab="Años",ylab="Talla",zlab="Proporcion",col = "lightblue", border=NA, shade=0.75)
matplot(Tallas,t(Sel_f),type="l",lty=1,xlab="Talla",ylab="Proporción",lwd = 1,col="black")
lines(Tallas,msex,col="red",lwd = 2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


persp(Yrs,Tallas,Sel_srv,theta=60,phi=40,expand=0.6, ticktype = "detailed",main="Selectividad campañas",
      xlab="Años",ylab="Talla",zlab="Proporcion",col = "lightblue", border=NA, shade=0.75)
matplot(Tallas,t(Sel_srv),type="l",lty=1,xlab="Talla",ylab="Proporción",lwd = 1,col="black")
lines(Tallas,msex,col="red",lwd = 2)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)



#Reclutamientos-------------

par(mfrow = c(1, 2))
sdR=Reclus[3,]
li=Reclus[1,]-1.96*sdR
ls=Reclus[1,]+1.96*sdR
nyrs=length(Bio_reprod[1,])
plot(Yrs,Reclus[2,],ylim = c(0,max(ls)),type="l",
     ylab="Reclutamientos",xlab="Año",pch = 16,lty=3,col="red")
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])
polygon(x,y,col="#DCDCDC",border="#DCDCDC")
lines(Yrs,Reclus[2,],lwd="2",col="red")
lines(Yrs,Reclus[1,],lwd="2")

barplot(Dev_log_R~Yrs,xlab="Años",ylab="Anomalía(log)")
abline(h=0,lwd = 1)
box()


par(mfrow = c(1, 1))

alfa=4*h*R0/(5*h-1)
beta=(1-h)*B0/(5*h-1)

ssb=seq(0,max(c(B0,max(Bio_reprod))),B0/50)
rec=alfa*ssb/(beta+ssb)

plot(Bio_reprod[1,seq(1,length(Yrs)-lag)],Reclus[1,seq(lag+1,length(Yrs))],
     ylim=c(0,max(Reclus[1,])),type="b",
     xlim=c(0,max(c(B0,max(Bio_reprod)))),
     xlab="Biomasa reproductiva",
     ylab="Reclutamiento")
lines(ssb,rec,col="blue",lwd=2)
text(Bio_reprod[1,seq(1,length(Yrs)-lag)],Reclus[1,seq(lag+1,length(Yrs))],
     paste(Yrs[seq(lag+1,length(Yrs))]),cex=0.8)

text(Bio_reprod[1,length(Yrs)-lag],Reclus[2,length(Yrs)],
     paste(Yrs[length(Yrs)]),cex=0.8,col="red")
abline(v=0.4*B0,col="red",lty=2)

#Analisis por recluta BIOMASA REPRODUCTIVA------------
par(mfrow = c(1, 1))

Sel=Sel_f[length(Yrs),]
M=Linf_k_M_h[3]
ypr_out<-por_recluta_r(Pre_r,Ptrans,Sel,msex,PesoL,target,h,M,dt,Tallas)
attach(ypr_out)
Ftar=Ftar
BPRtar=BPRtar
YPRtar=YPRtar

plot(Fcr,Y/max(Y),type="l", col="green", lwd=2, main="Analisis por recluta", xlab="Mortalidad por pesca", ylab="BPR, YPR relativos",
     ylim = c(0,1))
lines(Fcr,B/max(B), col="magenta", lwd=2)
text(Ftar,0.05,paste("Fref=",round(Ftar,2)))
abline(h = target, lty = 2,lwd=1)
abline(v = Ftar, lty = 2,lwd=1)
abline(v = Mort_F[1,length(Yrs)], lty = 2,lwd=1,col="red")
text(Mort_F[1,length(Yrs)],0.1,paste("Fcr=",round(Mort_F[1,length(Yrs)],2)),
     col="red")



#Biomasa, SPR y F con IC REPRODUCTIVA-----------------------------------------------------------
par(mfrow = c(2, 1))

li=Bio_reprod[1,]-1.96*Bio_reprod[2,]
ls=Bio_reprod[1,]+1.96*Bio_reprod[2,]

nyrs=length(Bio_reprod)

plot(Yrs,Bio_reprod[1,],ylim = c(0,max(Bio_reprod[1,]+1.96*Bio_reprod[2,])*1.01),type="l",
     ylab="Biomasa",xlab="Año",pch = 16,cex=1,lwd=2,main="Biomasa")

x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Bio_reprod[1,],lwd=2)
abline(h = target*B0, col = "green",lwd=2,lty=2)
abline(h = 0.5*target*B0, col = "red",lwd=2,lty=2)

# lines(Yrs,Bio_tot,col="blue",lty=2,lwd=2)
# lines(Yrs,Bio_explot[1,],,col="tomato",lty=2,lwd=2)
# lines(Yrs,Desembarques[1,],type="h",lwd=5,col="black")
# lines(Yrs,Desembarques[1,],type="h",lwd=5,col="black")
# lines(Yrs,Bcrucero[1,],type="p",lwd=5,col="black")



li=Mort_F[1,]-1.96*Mort_F[2,]
ls=Mort_F[1,]+1.96*Mort_F[2,]
nyrs=length(Mort_F)

plot(Yrs,Mort_F[1,],ylim = c(0,max(Mort_F[1,]+1.96*Mort_F[2,])*1.01),type="l",
     ylab="F",xlab="Año",pch = 16,cex=1,lwd=2,main="Mortalidad por pesca")
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Mort_F[1,],lwd=2)
abline(h = Ftar, col = "red",lty = 2,lwd=2)



#Kobe REPRODUCTIVA---------------------------------------------------------------------
par(mfrow = c(1, 1))
SPR=RPRLP
BRMS=B0*target
nysim=length(BRep_proy[1,])

nyrs=length(Yrs)
  
plot(0,0,pch = 16,ylab="F/Frms",xlab="B/Brms", xlim = c(0,max(SPR[1,]/target)), ylim = c(0,max(Mort_F[1,]/Ftar,1)*1.5), 
     type="o",col="white",lty="dashed",main=paste("B/Brms=",round(SPR[1,][length(Yrs)]/target,2),
                                                  " F/Frms=",round(Mort_F[1,][length(Yrs)]/Ftar,2)))

polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1") #amarillo
polygon(c(1,1.1*max(SPR[1,]/target),1.1*max(SPR[1,]/target),1),c(0,0,1,1),col="green") #verde
polygon(c(1,1.1*max(SPR[1,]/target),1.1*max(SPR[1,]/target),1),c(1,1,1.5*max(1,Mort_F[1,]/Ftar),1.5*max(1,Mort_F[1,]/Ftar)),col="orange") #amarillo
polygon(c(0,1,1,0),c(1,1,1.5*max(1,Mort_F[1,]/Ftar),1.5*max(1,Mort_F[1,]/Ftar)),col="red") #rojo

lines(SPR[1,]/target,Mort_F[1,]/Ftar,pch = 16, type="o",col="black",lty="dashed")
lines(SPR[1,nyrs]/target,Mort_F[1,][length(Yrs)]/Ftar,type="p",col="blue",pch = 16)
text(SPR[1,]/target*.95,Mort_F[1,]/Ftar,paste(Yrs),cex=0.8)


X0=Bio_reprod[1,nyrs]/BRMS
Y0=Mort_F[1,nyrs]/Ftar
cvF=Mort_F[2,nyrs]/Mort_F[1,nyrs];
cvB=Bio_reprod[2,nyrs]/Bio_reprod[1,nyrs];

arrows(X0,Y0-1.96*cvF*Y0,X0,Y0+1.96*cvF*Y0,
       length = 0.05, code = 3, angle = 90, lwd=2, col="blue")

arrows(X0-1.96*cvB*X0,Y0,X0+1.96*cvB*X0,Y0,
       length = 0.05, code = 3, angle = 90,lwd=2, col="blue")

box()




#N0--------------------------------------------------------------------------------------
par(mfrow = c(1, 1))

Lm=rep(0,4)

n=N0[1,]/max(N0[1,])
Lm[1]=Lr_Sr_beta[1]

for (i in 2:15)
{
  n=(n*exp(-Linf_k_M_h[3]))%*%Ptrans
  Lm[i]=sum(n*Tallas)/sum(n)}

matplot(Tallas,t(N0),type="l",lty=1,lwd=2,col="black",
        main=paste("Lr=",round(Lr_Sr_beta[1],2)," Sr=",round(Lr_Sr_beta[2],2),
                   "Loo=",round(Linf_k_M_h[1],2)," k=",round(Linf_k_M_h[2],2)))  
abline(v=Lm,col="red",lty=2)


#Proyecciones---------------------------------
par(mfcol = c(2, 1))

nanos=length(Yrs)
nysim=dim(BRep_proy)[2]
yproy=seq(Yrs[length(Yrs)]+1,Yrs[length(Yrs)]+nysim)
vecto=seq(1,nanos)

plot(Yrs[vecto],Bio_reprod[1,vecto]/B0,type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
     ylim=c(0,max(c(Bio_reprod[1,vecto],BRep_proy)/B0)),lwd=3,lty=1,xlab="Año", ylab="B/B0",
     main="Biomasa",col="darkgreen")
  abline(h=target,  col = "green",lty = 2)
  abline(h=0.5*target,  col = "red",lty = 2)
  abline(v=Yrs[length(Yrs)]+1,lty = 2,lwd=1)
  matlines(yproy,t(BRep_proy/B0),lwd="2",lty=1)
  
  RMS=YPRtar*R0
  
  plot(Yrs[vecto],Desembarques[1,vecto],type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
       ylim=c(0,max(max(Capt_proy),Desembarques)),lwd=3,lty=1,xlab="Año", ylab="Capturas",
       main="Capturas",col="darkgreen")
  matlines(yproy,t(Capt_proy),lwd=2,lty=1,type="l")
  abline(h=RMS,  col = "green",lty = 2)
  abline(v=Yrs[length(Yrs)]+1,lty = 2,lwd=1)
  
  legend("topright",paste(round(Fmult,2)),lty=1,bty="n",col=seq(1,length(Fmult)),
         cex=1.0,lwd=2)
  


#Genera excel--------------------------------------------------------------------------------------

cv_bio=Bio_reprod[2,length(Yrs)]/Bio_reprod[1,length(Yrs)]
risk.cp=1-pnorm(BRep_proy[,1]/B0,target,cv_bio*BRep_proy[,1]/B0)
risk.lp=1-pnorm(BRep_proy[,nysim]/B0,target,cv_bio*BRep_proy[,nysim]/B0)

Fproy=Mort_F[1,nyrs]*Fmult


Variables=data.frame(Mult_Eff=Fmult,Fproy=Fproy,B_B0.cp=BRep_proy[,1]/B0,B_B0.lp=BRep_proy[,nysim]/B0, 
                     Captura.cp=Capt_proy[,1],Captura.lp=Capt_proy[,nysim],p_Brms.cp=risk.cp,p_Brms.lp=risk.lp)

write.csv(Variables, 'Decision.csv',row.names = F)


ICi=Bio_reprod[1,]-1.96*Bio_reprod[2,]
ICs=Bio_reprod[1,]+1.96*Bio_reprod[2,]
p_low=1-pnorm((Bio_reprod[1,]),BRMS,(Bio_reprod[2,]))
p_crush=1-pnorm((Bio_reprod[1,]),0.5*BRMS,(Bio_reprod[2,]))
p_high=pnorm((Mort_F[1,]),Ftar,(Mort_F[2,]))

Variables2=data.frame(years=Yrs,Biom_rep=Bio_reprod[1,],ICi=ICi,ICs=ICs,Biom_explot=Bio_explot[1,],
                      p_Brms=p_low,
                      p_Blim=p_crush,
                      p_Frms=p_high,
                      Reclutamiento=Reclus[2,],
                      Fcr=Mort_F[1,],
                      F_Fmrs=Mort_F[1,]/Ftar,
                      B_Brms=Bio_reprod[1,]/BRMS)

write.csv(Variables2, 'Var_Pobl.csv',row.names = F)


# riesgo=1-pnorm(Red_stock[1,],target,Red_stock[2,])
# riesgo_crash=1-pnorm(Red_stock[1,],0.5*target,Red_stock[2,])

# Variables2=data.frame(Mult_Eff=Fmult,Fref=Fmult*Mort_F[1,nyrs],B_B0=BRep_proy[,length(BRep_proy[1,])]/,Riesgo=riesgo,Riesgo_colapso=riesgo_crash,Captura_cp=Capt_proy[1,], 
#                       Captura_lp=colMeans(Capt_proy[seq(nysim,nysim-4,-1),]))
# 
# 
# write.csv(Variables2, paste('Manejo_',name,'.csv'), row.names = F)



