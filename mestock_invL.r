rm(list=ls())
library(areaplot)

system('mod_centolla -ind centolla2.dat -nox')  # for model running

source('read.admb.R')
source('read.admbFit.R')
source('por_recluta_r.R')


data <-read.rep('for_R.rep')
attach(data)

par(family = "serif") # Cambia la fuente a "serif"


target=0.4
#Desembarques--------------------------------------------------------------------------
plot(Yrs, Desembarques[1,]/1000, cex.lab=1.5, type="h", xlab="Año",ylab="Miles t", col = "gray",
     ylim = c(0,max(Desembarques[1,]/1000)*1.1),lwd=6,
     main="Desembarques", cex.main = 1.5)
lines(Yrs, Desembarques[2,]/1000,col="red")

#Indices--------------

par(mfrow = c(2, 1))
sim=20

ubi=which(CPUE[1,]>0)
maxi=max(CPUE[1:2,ubi])*1.1;
mini=0;


plot(Yrs[ubi],CPUE[2,ubi],ylab="Indice",xlab="Año",main="CPUE",pch=sim,type="l",ylim=c(mini,maxi),
     col="red",lwd=2)

lines(Yrs[ubi],CPUE[1,ubi],pch=sim,type="b")

y=CPUE[1,]
cv=cv_CPUE

ubi=which(Bcrucero[1,]>0)
maxi=max(Bcrucero[1:2,ubi])*1.1;
mini=0;


plot(Yrs[ubi],Bcrucero[2,ubi],ylab="Indice",xlab="Año",main="Campaña",pch=sim,type="l",ylim=c(mini,maxi),
     col="red",lwd=2)

lines(Yrs[ubi],Bcrucero[1,ubi],pch=sim,type="b")


#Talla promedio------------------------------------------------------

par(mfrow = c(2, 1))

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

par(mfcol = c(4, 3))
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
par(mfrow = c(2, 1))

if (length(Sel_f)>length(PesoL)){
matplot(Tallas,t(Sel_f),col="gray",type="l",lwd=2,lty = 1,xlab="Talla",ylab="Proporcion",main="Selectividad Flota")
lines(Tallas,msex,col="black",lwd = 2)
 } else {
plot(Tallas,Sel_f,col="gray",type="l",lwd=2,lty = 1,xlab="Talla",ylab="Proporcion",main="Selectividad Flota")
lines(Tallas,msex,col="black",lwd = 2)
}
abline(v=90,col="green",lwd=2)
abline(v=110,col="red",lwd=2)


if (length(Sel_srv)>length(PesoL)){
  matplot(Tallas,t(Sel_srv),col="gray",type="l",lwd=2,lty = 1,xlab="Talla",ylab="Proporcion",main="Selectividad Campaña")
  lines(Tallas,msex,col="black",lwd = 2)
} else {
  plot(Tallas,Sel_srv,col="gray",type="l",lwd=2,lty = 1,xlab="Talla",ylab="Proporcion",main="Selectividad Flota")
  lines(Tallas,msex,col="black",lwd = 2)
}
abline(v=90,col="green",lwd=2)
abline(v=110,col="red",lwd=2)

#Reclutamientos-------------
par(family = "serif") # Cambia la fuente a "serif"
par(mfrow = c(1, 1))
sdR=Reclus[3,]
li=Reclus[1,]-1.96*sdR
ls=Reclus[1,]+1.96*sdR
nyrs=length(Bio_reprod)
plot(Yrs,Reclus[2,],ylim = c(0,max(ls)),type="l",
     ylab="Reclutamientos (escala relativa)",main= "Reclutamientos", xlab="Año",cex.lab= 1.2,pch = 16,lty=3,col="red")
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])
polygon(x,y,col="#DCDCDC",border="#DCDCDC")
lines(Yrs,Reclus[2,],lwd="2",col="red")
lines(Yrs,Reclus[1,],lwd="2")

abline(h = c(0, 10, 20, 30,40,50), col = "black",lty = 2,lwd=1)
abline(v = c (2000, 2005, 2010, 2015, 2020),col = "black",lty = 2,lwd=1)
box()

barplot(Dev_log_R~Yrs,xlab="Año",ylab="Anomalía (log)", cex.lab= 1.2, main= "Anomalías del Reclutamiento", ylim = c(-1.2, 1.2))
abline(h=0,lwd = 1)



box()
abline(h = c(0.0, -0.5,0.5,1.0, -1.0), col = "black",lty = 2,lwd=1)
abline(v = c (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),col = "gray" ,lty = 2,lwd=1)

par(mfrow = c(1, 1))

alfa=4*h*R0/(5*h-1)
beta=(1-h)*B0/(5*h-1)

ssb=seq(0,max(c(B0,max(Bio_reprod))),B0/50)
rec=alfa*ssb/(beta+ssb)

plot(Bio_reprod[1,seq(1,length(Yrs)-lag)],main= "Relación Stock-recluta", cex.main= 1.5, Reclus[1,seq(lag+1,length(Yrs))],
     ylim=c(0,max(Reclus[1,])),type="b",
     xlim=c(0,max(c(B0,max(Bio_reprod)))),
     xlab="Biomasa reproductiva (t)", cex.lab=1.25,
     ylab="Reclutamiento (escala relativa)")
lines(ssb,rec,col="blue",lwd=2)
text(Bio_reprod[1,seq(1,length(Yrs)-lag)],Reclus[1,seq(lag+1,length(Yrs))],
     paste(Yrs[seq(lag+1,length(Yrs))]),cex=0.8)

text(Bio_reprod[1,length(Yrs)-lag],Reclus[2,length(Yrs)],
     paste(Yrs[length(Yrs)]),cex=0.8,col="red")
abline(v=0.4*B0,col="red",lty=3, lwd=2.5)

abline(h = c(0, 5,10,15, 20), col = "gray",lty = 2,lwd=1)
abline(v = c (0, 5000, 10000, 15000,  20000,25000),col = "gray",lty = 2,lwd=1)

#Analisis por recluta  BIOMASA REPRODUCTIVA------------
par(mfrow = c(1, 1))

Sel=Sel_f[length(Yrs),]
M=Linf_k_M_h[3]
ypr_out<-por_recluta_r(Pre_r,Ptrans,Sel,msex,PesoL,target,h,M,dt,Tallas)
attach(ypr_out)
Ftar=Ftar
BPRtar=BPRtar
YPRtar=YPRtar
par(family = "serif") # Cambia la fuente a "serif"http://127.0.0.1:33695/graphics/plot_zoom_png?width=706&height=503
plot(Fcr,Y/max(Y),type="l", col="green", lwd=2, main="Análisis por recluta", cex.main=1.5, 
     xlab="Mortalidad por pesca (F)", ylab="BPR, YPR (escala relativa)", cex.lab=1.2,
     ylim = c(0,1))
lines(Fcr,B/max(B), col="magenta", lwd=2)
text(Ftar,0.05,paste("Fref=",round(Ftar,2)), font=2)
abline(h = target, lty = 2,lwd=1)
abline(v = Ftar, lty = 2,lwd=1)
#abline(v = Mort_F[1,length(Yrs)], lty = 2,lwd=1,col="red")

abline(h = c(0.0,0.2,0.4,0.6,0.8,1.0), col = "gray",lty = 3,lwd=0.8)
abline(v = c (0.0,0.5,1.0,1.5),col = "gray",lty = 3,lwd=0.8)
points(Ftar, target, col = "black", pch = 19, cex = 1.5) 

#Biomasa, SPR y F con IC REPRODUCTIVA-----------------------------------------------------------
par(mfrow = c(1, 1))

li=Bio_reprod[1,]-1.96*Bio_reprod[2,]
ls=Bio_reprod[1,]+1.96*Bio_reprod[2,]

nyrs=length(Bio_reprod)

plot(Yrs,Bio_reprod[1,],ylim = c(0,max(Bio_reprod[1,]+1.96*Bio_reprod[2,])*1.01),type="l",
     ylab="Biomasa (t)",xlab="Año",cex.lab=1.25, pch = 16,cex=1,lwd=2,main="Biomasa")



x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Bio_reprod[1,],lwd=2)
abline(h = target*B0, col = "green",lty = 2,lwd=2)
abline(h = 0.5*target*B0, col = "red",lty = 2,lwd=2)

li=Mort_F[1,]-1.96*Mort_F[2,]
ls=Mort_F[1,]+1.96*Mort_F[2,]
nyrs=length(Mort_F)

abline(h = c(0, 5000, 10000, 15000, 20000, 25000), col = "black",lty = 2,lwd=1)
abline(v = c (2000, 2005, 2010, 2015, 2020),col = "black",lty = 2,lwd=1)

plot(Yrs,Mort_F[1,],ylim = c(0,max(Mort_F[1,]+1.96*Mort_F[2,])*1.01),type="l",
     ylab="F",xlab="Año",cex.lab=1.25,pch = 16,cex=1,lwd=2,main="Mortalidad por pesca")
x=c(Yrs,Yrs[seq(length(Yrs),1,-1)])
y=c(li,ls[seq(length(ls),1,-1)])

polygon(x,y,col="#DCDCDC",border="#DCDCDC")  

lines(Yrs,Mort_F[1,],lwd=2)
abline(h = Ftar, col = "red",lty = 2,lwd=2)

abline(h = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), col = "black",lty = 2,lwd=1)
abline(v = c (2000, 2005, 2010, 2015, 2020),col = "black",lty = 2,lwd=1)

#Kobe REPRODUCTIVA---------------------------------------------------------------------
par(mfrow = c(1, 1))
SPR=RPRLP
BRMS=B0*target
nysim=length(BRep_proy[1,])

nyrs=length(Yrs)
  
plot(SPR[1,]/target,Mort_F[1,]/Ftar,pch = 16,ylab="F/Frms",xlab="B/Brms", xlim = c(0,max(SPR[1,]/target)), ylim = c(0,max(Mort_F[1,]/Ftar)*1.5), cex.lab=1.2,
     type="o",col="black",lty="dashed",main=paste("B/Brms=",round(SPR[1,][length(Yrs)]/target,2),
                                                  " F/Frms=",round(Mort_F[1,][length(Yrs)]/Ftar,2)))

polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1") #amarillo
polygon(c(1,1.1*max(SPR[1,]/target),1.1*max(SPR[1,]/target),1),c(0,0,1,1),col="green") #verde
polygon(c(1,1.1*max(SPR[1,]/target),1.1*max(SPR[1,]/target),1),c(1,1,1.5*max(Mort_F[1,]/Ftar),1.5*max(Mort_F[1,]/Ftar)),col="orange") #amarillo
polygon(c(0,1,1,0),c(1,1,1.5*max(Mort_F[1,]/Ftar),1.5*max(Mort_F[1,]/Ftar)),col="red") #rojo

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
Lm[1]=sum(n*Tallas)/sum(n)

plot(Tallas,n,type="l",lty=1,col="black",xlab="Talla(mm)",ylab="Densidad",
        main=paste("Lr=",round(Lr_Sr_beta[1],2)," Sr=",round(Lr_Sr_beta[2],2)),lwd=2)

for (i in 2:22)
{
  n=(n*exp(-Linf_k_M_h[3]))%*%Ptrans
  lines(Tallas,n,lwd=2)
  
  Lm[i]=sum(n*Tallas)/sum(n)}
  
#abline(v=Lm,col="red",lty=2)

#Biomasa RMS vs Capturas en equilibrio

plot(B/max(B),Y*R0,type="l",lwd=2,xlab="Biomasa relativa (B/B0)",ylab="Capturas en equilibrio (t)",cex.lab=1.2,
     main=paste("Brms/B0=",round(B[which(Y==max(Y))]/max(B),2),"  RMS=",round(max(Y*R0),0)))
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(h = max(Y*R0), lty = 2,lwd=1.5, col="red")
abline(v = B[which(Y==max(Y))]/max(B), lty = 2,lwd=1.5, col="red")

#Proyecciones---------------------------------

par(family = "serif") # Cambia la fuente a "serif"
par(mfcol = c(1, 1))

nyrs=10
nanos=length(Yrs)
nysim=dim(BRep_proy)[2]
yproy=seq(Yrs[length(Yrs)]+1,Yrs[length(Yrs)]+nysim)
vecto=seq(1,25)

plot(Yrs[vecto],Bio_reprod[1,vecto]/B0,type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
     ylim=c(0,max(c(Bio_reprod[1,vecto],BRep_proy)/B0)),lwd="2",lty=1,xlab="Año", ylab="B/B0", cex.lab=1.2,
     main="Biomasa", cex.main=1.5)
abline(h=target,  col = "green",lty = 2)
abline(h=0.5*target,  col = "red",lty = 2)
abline(v=Yrs[length(Yrs)]+1,lty = 2,lwd=1)
matlines(yproy,t(BRep_proy/B0),lwd="2",lty=1)

legend("topleft",paste(round(Fmult,2)),lty=1,bty="n",col=seq(1,length(Fmult)),
       cex=0.5,lwd=2)

RMS=YPRtar*R0

plot(Yrs[vecto],Desembarques[1,vecto],type="l",xlim=c(min(Yrs[vecto]),max(yproy)),
     ylim=c(0,max(Desembarques[1,vecto])),lwd="2",lty=1,xlab="Año", ylab="Capturas (t)", cex.lab=1.2,
     main="Capturas", , cex.main=1.5)
matlines(yproy,t(Capt_proy),lwd=2,lty=1)
abline(h=RMS,  col = "green",lty = 2)

legend("topleft",paste(round(Fmult,2)),lty=1,bty="n",col=seq(1,length(Fmult)),
       cex=0.5,lwd=2)


# Probavilidades
par(mfrow = c(2, 1))

p_low=pnorm(SPR[1,nyrs],target,cvB*SPR[1,nyrs])
p_high=pnorm(Mort_F[1,nyrs],Ftar,Mort_F[2,nyrs])

eje=seq(0,1,by=0.005)
d1=dnorm(eje,SPR[1,nyrs],cvB*SPR[1,nyrs])
p1=pnorm(eje,SPR[1,nyrs],cvB*SPR[1,nyrs])

eje2=seq(0,2*Mort_F[1,nyrs],by=0.01)
d2=dnorm(eje2,Mort_F[1,nyrs],Mort_F[2,nyrs])
p2=pnorm(eje2,Mort_F[1,nyrs],Mort_F[2,nyrs])


areaplot(eje/target,d1,col="cyan",lwd=2,ylab="Riesgo", xlab="B/Brms",
         main=paste("p(B<Brms)=",round(p_low,3)))
abline(v=SPR[1,nyrs]/target)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(v=1,col="blue",lwd=2)

areaplot(eje2/Ftar,d2,col="pink",lwd=2,ylab="Riesgo", xlab="F/Frms",
         main=paste("p(F>Frms)=",round(p_high,3)))
abline(v=Mort_F[1,nyrs]/Ftar)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
abline(v=1,col="red",lwd=2)


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




