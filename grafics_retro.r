rm(list=ls()) # erasure all objects

source('read.admb.R')
source('read.admbFit.R')
source('por_recluta_r.R')
library(matlib)

nsc=6 #defino el numero de escenarios
nyrs=23 #numero de años

Bio=matrix(NA, ncol = nyrs, nrow = nsc)
Deple=matrix(NA, ncol = nyrs, nrow = nsc)
Fmort=matrix(NA, ncol = nyrs, nrow = nsc)
tabla=matrix(NA, ncol = 15, nrow = 7)
tabla2=matrix(NA, ncol = 5, nrow = 7)



# 
#  for (i in 1:nsc)
#  {
#    system(paste('mod_centolla -ind S2_',i-1,'.dat -nox',sep=""))
#    shell(paste("copy for_R.rep for_S2_",i-1,".rep -nox",sep=""))
#    }


t=2000:2022

for (i in 1:nsc)
{
  
  name=paste('for_S2_',i-1,'.rep',sep="")
  data <-read.rep(name)

  attach(data)
  

  Bio[i,1:length(Yrs)]=Bio_reprod[1,]
  Deple[i,1:length(Yrs)]=RPRLP[1,]

}

# 
# 
# write.csv(tabla, 'Likeval.csv', append = FALSE, sep = " ", dec = ".",
#           row.names = F, col.names = T)
# 
# write.csv(tabla2, 'Vars_casos.csv', append = FALSE, sep = " ", dec = ".",
#           row.names = F, col.names = T)
# 

par(mfrow = c(1, 2))
matplot(t,t(Bio),type="l",lwd=2,ylim=c(0,30000),lty=1,col=c(1:6),
        ylab="Biomasa",xlab="Año")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)


matplot(t,t(Deple),type="l",lwd=2,ylim=c(0,1.5),lty=1,col=c(1:6),
        ylab="B/B0",xlab="Año")
abline(h=0.4,col="red")
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)





