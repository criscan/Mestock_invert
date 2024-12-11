por_recluta_r<-function(Pre_r,Ptrans,Sel,msex,PesoL,tar,h,M,dt,Tallas){


  N=seq(1,length(Tallas))
  N0=N
  Ntar=N
  Ctar=N
  
  Fcr=seq(0.001,5*M,0.01)
  B=Fcr
  BE=Fcr
  Y=Fcr
  R=Fcr
  
  # Calcula B0
  
  R0=1
  N0=(R0*Pre_r)%*%solve(diag(length(Tallas))-Ptrans%*%(diag(length(Tallas))*exp(-M)))
  B0=sum(N0*msex*PesoL*exp(-dt*M))
  alfa=4*h*R0/(5*h-1);
  beta=(1-h)/(5*h-1)*B0;
  aux=0
  
  for (j in 1:length(Fcr)) {
    
    F=Fcr[j]*Sel
    Z=F+M
    N=R0*Pre_r     
    
    for (t in 1:15) {
    N=(N*exp(-Z))%*%Ptrans+R0*Pre_r
    }
    
    N=(R0*Pre_r)%*%solve(diag(length(Tallas))-Ptrans%*%(diag(exp(-Z))))
    
    N=as.vector(N)
    C=N*F/Z*(1-exp(-Z))
    
    BE[j]=sum(N*Sel*PesoL*exp(-0.5*Z))
    B[j]=alfa*sum(N*msex*PesoL*exp(-dt*Z))-beta
    R[j]=alfa*B[j]/(beta+B[j])
    Y[j]=R[j]*sum(C*PesoL)

    
    if (B[j]/B[1]>0.99*tar){Ftar=Fcr[j]
    BPRtar=B[j]
    YPRtar=Y[j]
    }
    
    if (Y[j]>aux){
      BRMS=B[j]
      YRMS=Y[j]
      FRMS=Fcr[j]
      }

    aux=Y[j]
    
  }
  

  outputs=list(Fcr=Fcr, B=B, BE=BE, Y=Y, Ftar=Ftar, BPRtar=BPRtar, YPRtar=YPRtar,
               FRMS=FRMS, BRMS=BRMS, YRMS=YRMS)

  return(outputs)  

  
}

