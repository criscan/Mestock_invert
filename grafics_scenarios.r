rm(list=ls()) # erasure all objects

source('read.admb.R')
source('read.admbFit.R')
source('por_recluta_r.R')
library(matlib)

nsc=9 #defino el numero de escenarios
nyrs=19 #numero de a?os

Bio=matrix(0, ncol = nyrs, nrow = nsc)
Deple=matrix(0, ncol = nyrs, nrow = nsc)
Fmort=matrix(0, ncol = nyrs, nrow = nsc)
tabla=matrix(0, ncol = 18, nrow = nsc)
tabla2=matrix(0, ncol = 5, nrow = nsc)


# for (i in 8:nsc)
# {
#   system(paste('modela_v1 -ind lango_S',i-1,'.dat -nox',sep=""))
#   shell(paste("copy for_R.rep for_RS",i-1,".rep",sep=""))
#   shell(paste("copy modela_v1.std lango_S",i-1,".std",sep=""))
#   shell(paste("copy modela_v1.par lango_S",i-1,".par",sep=""))
# }


for (i in 1:nsc)
{
  
  name=paste('for_RS',i-1,'.rep',sep="")
  data <-read.rep(name)
  
  attach(data)
  
  
  Sel=Sel_f[length(Yrs),]
  M=Linf_k_M_h[3]
  ypr_out<-por_recluta_r(Pre_r,Ptrans,Sel,msex,PesoL,0.4,h,M,dt,Tallas)
  attach(ypr_out)
  F40=0.87*M #Ftar[1]
  
  
  Bio[i,1:length(Yrs)]=Bio_reprod[1,]
  Deple[i,1:length(Yrs)]=RPRLP[1,]
  
  
  Fmort[i,1:length(Yrs)]=Mort_F[1,]/F40
  tabla[i,]=Likeval
  tabla2[i,1:5]=c(B0, R0, RPRLP[1,nyrs-1],Mort_F[1,nyrs-1]/F40, Bio_reprod[1,nyrs-1])
  
  
}



write.csv(tabla, 'Likeval.csv', append = FALSE, sep = " ", dec = ".",
          row.names = F, col.names = T)

write.csv(tabla2, 'Vars_casos.csv', append = FALSE, sep = " ", dec = ".",
          row.names = F, col.names = T)


par(mfrow = c(1, 2))
matplot(Yrs,t(Deple),type="l",lwd=2,ylim=c(0,max(Deple)),lty=1,col=c(1:nsc),
        ylab="B/B0",xlab="Año")
abline(h=0.4,col="red")
#grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)

matplot(Yrs,t(Fmort),type="l",lwd=2,ylim=c(0,max(Fmort)),col=c(1:nsc),lty=1,
        ylab="F/Frms",xlab="Año")
legend("topleft",c("S0","S1","S2","S3","S4","S5","S6","S7","S8"),col=c(1:nsc),
       lty=1,lwd=2,cex=0.8,bty="n")
abline(h=1.0,col="red")
#grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)



