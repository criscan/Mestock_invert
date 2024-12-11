GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");

TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 


DATA_SECTION
 init_int ymin
 init_int ymax

 number ntime
 !!ntime=ymax-ymin+1;

 init_int ntallas

 init_matrix mdatos(1,ntime,1,11)
 init_vector Tallas(1,ntallas)

 init_int N_ftc
 init_ivector nanos_ftc(1,N_ftc)
 init_matrix Ctot(1,N_ftc,1,ntallas)

 init_int N_fts1
 init_ivector nanos_fts1(1,N_fts1)
 init_matrix Ncru(1,N_fts1,1,ntallas)

 init_number sigmaR
 init_number L50m
 init_number L95m
 init_number dts
 init_number log_aw
 init_number bw
 init_int    minedad
 init_int    nedades
 
 init_vector Par_bio(1,7)
 init_vector cv_Par_bio(1,7)
 init_int    fase_Loo
 init_int    fase_k
 init_int    fase_Lr
 init_int    fase_sr
 init_int    fase_beta
 init_int    fase_M
 init_int    fase_h

 
  number log_Linf_prior
  number log_k_prior
  number log_Lr_prior
  number log_sr_prior
  number log_beta_prior
  number log_M_prior
  number log_h_prior

  
  !! log_Linf_prior = log(Par_bio(1));
  !! log_k_prior = log(Par_bio(2));
  !! log_Lr_prior = log(Par_bio(3));
  !! log_sr_prior = log(Par_bio(4));
  !! log_beta_prior= log(Par_bio(5));
  !! log_M_prior= log(Par_bio(6));
  !! log_h_prior= log(Par_bio(7));


 init_number L50prior
 init_number s1prior
 init_number s2prior
 init_number cvf1
 init_number cvf2
 init_number cvf3
 init_int fase_f1
 init_int fase_f2
 init_int fase_f3

 number log_L50prior
 number log_s1prior
 number log_s2prior


 !! log_L50prior = log(L50prior);
 !! log_s1prior = log(s1prior);
 !! log_s2prior = log(s2prior);


 init_number L50priorC
 init_number s1priorC
 init_number s2priorC
 init_number cvc1
 init_number cvc2
 init_number cvc3
 init_int fase_c1
 init_int fase_c2
 init_int fase_c3

 number log_L50priorC
 number log_s1priorC
 number log_s2priorC


 !! log_L50priorC = log(L50priorC);
 !! log_s1priorC = log(s1priorC);
 !! log_s2priorC = log(s2priorC);
 
 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)

 init_int    nbloques2
 init_vector ybloques2(1,nbloques2)

 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)

 init_int    nqbloquesc
 init_vector yqbloquesc(1,nqbloquesc)

 init_number qprior
 init_number cv_log_qprior

 number log_qprior
 !! log_qprior = log(qprior);


 // FASES DE ESTIMACION 
 init_int    opt_qf
 init_int    opt_qc

 init_int    opt_F
 init_int    opt_devRt
 init_int    fasepropR0

 init_int    npbr
 init_vector pbr(1,npbr)
 init_int ntime_sim



INITIALIZATION_SECTION

  log_Linf       log_Linf_prior
  log_k          log_k_prior
  log_Lr         log_Lr_prior
  log_sr         log_sr_prior
  log_beta       log_beta_prior
  log_M          log_M_prior
  log_h          log_h_prior
  log_L50        log_L50prior 
  log_L50c       log_L50priorC 
  log_sigma1     log_s1prior 
  log_sigma2     log_s2prior
  log_sigma1c    log_s1priorC 
  log_sigma2c    log_s2priorC
  log_beta       log_beta_prior 
  log_qcru       log_qprior

PARAMETER_SECTION


 init_vector log_L50(1,nbloques1,fase_f1)  
 init_vector log_sigma1(1,nbloques1,fase_f2)
 init_vector log_sigma2(1,nbloques1,fase_f3)


 init_vector log_L50c(1,nbloques2,fase_c1)  
 init_vector log_sigma1c(1,nbloques2,fase_c2)
 init_vector log_sigma2c(1,nbloques2,fase_c3)

// parametros reclutamientos y mortalidades)
 init_number log_Rmed(1)
 init_bounded_dev_vector log_desv_Rt(1,ntime,-10,10,opt_devRt)
 init_bounded_vector log_F(1,ntime,-20,0.7,opt_F) // log  mortalidad por pesca por flota
 init_number log_propR0(fasepropR0)


// capturabilidades
 init_vector log_qflo(1,nqbloques,opt_qf)
 init_vector log_qcru(1,nqbloquesc,opt_qc)

// Crecimiento
 init_number log_Linf(fase_Loo)
 init_number log_k(fase_k)
 init_number log_Lr(fase_Lr)
 init_number log_sr(fase_sr)
 init_number log_beta(fase_beta)
 init_number log_M(fase_M)
 init_number log_h(fase_h)

//---------------------------------------------------------------------------------
//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector BMcru(1,ntime)
 vector Brec(1,ntime)
 vector pred_CPUE(1,ntime);
 vector pred_Bcru(1,ntime);
 vector pred_Desemb(1,ntime);
 vector likeval(1,18);
 vector Neq(1,ntallas);

 vector Rpred(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_year(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector delta(1,ntallas)
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)
 vector pre(1,ntallas)
 vector msex(1,ntallas)
 vector Wmed(1,ntallas)

 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector BDo(1,ntime);
 vector prior(1,7)
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector Bcru(1,ntime);

 vector Lmed_obs(1,N_ftc);
 vector Lmed_pred(1,N_ftc);
 vector Lmed_obsc(1,N_fts1);
 vector Lmed_predc(1,N_fts1);

 vector edades(1,nedades)
 vector nm(1,ntime)
 vector nmc(1,ntime)
 vector dt_cpue(1,ntime)
 vector dt_camp(1,ntime)

 matrix cv_index(1,4,1,ntime)

 matrix S1(1,nbloques1,1,ntallas)
 matrix S2(1,nbloques2,1,ntallas)

 matrix Sel(1,ntime,1,ntallas)
 matrix Selc(1,ntime,1,ntallas)
 
 matrix F(1,ntime,1,ntallas)
 matrix Z(1,ntime,1,ntallas)
 matrix S(1,ntime,1,ntallas)


 matrix N(1,ntime,1,ntallas)

 matrix NM(1,ntime,1,ntallas)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)
 matrix NVcru(1,ntime,1,ntallas)
 matrix No(1,2*nedades,1,ntallas)


 matrix pred_Ctot(1,ntime,1,ntallas)

 matrix pobs(1,N_ftc,1,ntallas)
 matrix ppred(1,N_ftc,1,ntallas)
 matrix pobsc(1,N_fts1,1,ntallas)
 matrix ppredc(1,N_fts1,1,ntallas)

 matrix T(1,ntallas,1,ntallas)

 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,ntallas)

 number suma1
 number suma2
 number suma3
 number suma4

 number penalty

 number So
 number alfa
 number beta

 number Linf
 number k
 number Linfh
 number M
 number Lr
 number sr
 number h

 number BDp
 number Npplus
 number Bp_anch 

 number nm1;
 number cuenta1;
 number alfa_sr;
 number beta_sr;
 number Rp;
 number rango;

 vector Np(1,ntallas)
 vector Zpbr(1,ntallas)
 vector Fpbr(1,ntallas)
 vector Sp(1,ntallas)

 matrix Bp(1,npbr,1,ntime_sim)
 vector CTPp(1,ntallas)
 matrix Yp(1,npbr,1,ntime_sim)

 
 objective_function_value f
  
 sdreport_vector Bexplot(1,ntime)
 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 sdreport_vector RPRlp(1,ntime) // 
 sdreport_vector Reclutas(1,ntime)
 sdreport_vector Fmort(1,ntime)
 
 sdreport_number SSBo
 sdreport_vector CBAtmas1(1,npbr)
 sdreport_vector BDBD0LP(1,npbr)
 sdreport_vector BDBD0tmas2(1,npbr)


PRELIMINARY_CALCS_SECTION

 yrs=column(mdatos,1);
 Desemb=column(mdatos,2);
 CPUE=column(mdatos,4);
 Bcru=column(mdatos,6);
 nm=column(mdatos,8);
 nmc=column(mdatos,9);
 dt_cpue=column(mdatos,10);
 dt_camp=column(mdatos,11);
 

 edades.fill_seqadd(minedad,1);

 cv_index(1)=column(mdatos,3);
 cv_index(2)=column(mdatos,5);
 cv_index(3)=column(mdatos,7);

 Wmed=exp(log_aw)*pow(Tallas,bw);
 rango=(L95m-L50m);
 msex=1/(1+exp(-log(19)*(Tallas-L50m)/rango));

 Unos_tallas=1;// lo uso en operaciones matriciales con tallas
 Unos_year=1;// lo uso en operaciones matriciales con el año


RUNTIME_SECTION
  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-3,1e-5
  maximum_function_evaluations 100,100,200,3000,3500

PROCEDURE_SECTION
// se listan las funciones que contienen los calculos
 Eval_Trans_talla_talla();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();

 if(last_phase()){Eval_CTP();}


//-----------------------------------------------------------------
FUNCTION Eval_Trans_talla_talla

   Linf=exp(log_Linf);
   k=exp(log_k);
   Lr=exp(log_Lr);
   sr=exp(log_sr);
   beta=exp(log_beta);
   M=exp(log_M);
   h=exp(log_h);

 int i, j;
 
// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas

  for (i=1;i<=ntallas;i++){

  if(delta(i)<0){delta(i)=0.01;}

  }

  Lesp=Tallas+delta; // talla esperada luego del crecimiento
  sigmaL=delta*beta;  

  for (i=1;i<=ntallas;i++){
    for (j=1;j<=ntallas;j++){
      if(i==j){
         T(i,j)=1.0;}}
   }


  for (i=1;i<=ntallas;i++){

    for (j=1;j<=ntallas;j++){
     if(sigmaL(i)>0){
     T(i,j)=mfexp(-0.5*square((Lesp(i)-Tallas(j))/sigmaL(i)));}}
   }


  for (j=1;j<=ntallas;j++){
  T(j)/=sum(T(j));
  } 


//----------------------------------------------------------------------

FUNCTION Eval_selectividad
 int i,j;

 // FLOTA...................

 for (j=1;j<=nbloques1;j++){

 S1(j)=exp(-0.5*square(Tallas-exp(log_L50(j)))/square(exp(log_sigma1(j))));
 
 // S1(j)=1./(1+exp(-log(19)*(Tallas-exp(log_L50(j)))/exp(log_sigma1(j))));

    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }}



 }

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }
   }

 // CRUCERO...................

 for (j=1;j<=nbloques2;j++){

  S2(j)=exp(-0.5*square(Tallas-exp(log_L50c(j)))/square(exp(log_sigma1c(j))));
 //  S2(j)=1./(1+exp(-log(19)*(Tallas-exp(log_L50c(j)))/exp(log_sigma1c(j))));


    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50c(j))){
      S2(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50c(j)))/square(exp(log_sigma2c(j))));
      }}

  }



   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques2;j++){
              if (yrs(i)>=ybloques2(j)){
                Selc(i)=S2(j);}
       }
   }


FUNCTION Eval_mortalidades

 Fmort=exp(log_F);

 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_tallas));

 Z=F+M;

 S=mfexp(-1.0*Z);


FUNCTION Eval_abundancia
 int i, j;

  Lr=Par_bio(3);
  sr=Par_bio(4);

  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);


// genero la poblacion virginal de LP;
  No(1)=pre*exp(log_Rmed);//-0.5*square(sigmaR));

  for (int j=2;j<=2*nedades;j++){
  No(j)=(No(j-1)*exp(-1.*M))*T;
  }


  SSBo=sum(elem_prod(colsum(No)*mfexp(-dts*M),elem_prod(Wmed,msex)));
  alfa_sr=4*h*exp(log_Rmed)/(5*h-1);//
  beta_sr=(1-h)*SSBo/(5*h-1);// Reclutamiento


// genero la poblacion equilibrio inicial;
  No(1)=exp(log_propR0)*pre*exp(log_Rmed);//-0.5*square(sigmaR));

  for (int j=2;j<=2*nedades;j++){
  No(j)=elem_prod(No(j-1),exp(-Z(1)))*T;
  }



// -----------------primer año
  Reclutas(1)=mfexp(log_Rmed);//-0.5*square(sigmaR));
  Rpred(1)=Reclutas(1);

  N(1)=colsum(No);
  NMD(1)=elem_prod(elem_prod(N(1),mfexp(-dts*Z(1))),msex);
  BD(1)=sum(elem_prod(Wmed,NMD(1)));


// --------------------dinamica anual
  for (i=2;i<=ntime;i++){

  Reclutas(i)=mfexp(log_Rmed+log_desv_Rt(i));//-0.5*square(sigmaR));
  Rpred(i)=Reclutas(i);
  
  if(i>minedad){

  Rpred(i)=(alfa_sr*BD(i-minedad)/(beta_sr+BD(i-minedad)));
  Reclutas(i)=Rpred(i)*mfexp(log_desv_Rt(i));}//-0.5*square(sigmaR)); }

  N(i)=elem_prod(N(i-1),S(i-1))*T+pre*Reclutas(i);
  NMD(i)=elem_prod(elem_prod(N(i),mfexp(-dts*Z(i))),msex);
  BD(i)=sum(elem_prod(Wmed,NMD(i)));

  } //


FUNCTION Eval_biomasas

 for (int i=1;i<=ntime;i++){ 
 NVflo(i)=elem_prod(elem_prod(N(i),mfexp(-1.*dt_cpue(i)*Z(i))),Sel(i));
 NVcru(i)=elem_prod(elem_prod(N(i),mfexp(-1.*dt_camp(i)*Z(i))),Selc(i));
 }

// vectores de biomasas derivadas

 BMflo=Wmed*trans(NVflo);
 Bexplot=Wmed*trans(elem_prod(N,Sel));
 BMcru=Wmed*trans(NVcru);
 BT=Wmed*trans(N);

 RPRlp=BD/SSBo;


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y año
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// vectores de desembarques predichos por año
 pred_Desemb=Wmed*trans(pred_Ctot);

// matrices de proporcion de capturas por talla y año
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));

 for (int i=1;i<=N_ftc;i++){
 ppred(i)=pred_Ctot(nanos_ftc(i)-ymin+1)/sum(pred_Ctot(nanos_ftc(i)-ymin+1));
 }


 pobsc=elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 for (int j=1;j<=N_fts1;j++){
 ppredc(j)=NVcru(nanos_fts1(j)-ymin+1)/sum(NVcru(nanos_fts1(j)-ymin+1));
 }

 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);

 Lmed_predc=Tallas*trans(ppredc);
 Lmed_obsc=Tallas*trans(pobsc);


FUNCTION Eval_indices
 
   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*BMflo(i);}
       }
   }


   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesc;j++){
              if (yrs(i)>=yqbloquesc(j)){
                 pred_Bcru(i)=exp(log_qcru(j))*BMcru(i);}
       }
   }


FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {
  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}

  if (Bcru(i)>0){
    suma2+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(3,i));}
 }



FUNCTION Eval_funcion_objetivo

 suma3=0; suma4=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;//Bcru

 likeval(3)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb

 for (int i=1;i<=N_ftc;i++){
 suma3+=-nm(nanos_ftc(i)-ymin+1)*sum(elem_prod(pobs(i),log(ppred(i)+1e-10)));
 }

 for (int i=1;i<=N_fts1;i++){
 suma4+=-nmc(nanos_fts1(i)-ymin+1)*sum(elem_prod(pobsc(i),log(ppredc(i)+1e-10)));
 }

 likeval(4)=suma3;//
 likeval(5)=suma4;//

// lognormal Ninicial y Reclutas
 if(active(log_desv_Rt)){
 likeval(6)=1./(2*square(sigmaR))*norm2(log_desv_Rt);}

 if (active(log_L50)){
 likeval(7)=1./(2*square(cvf1))*norm2(log_L50-log_L50prior);}

 if (active(log_sigma2)){
 likeval(8)=1./(2*square(cvf3))*norm2(log_sigma2-log_s2prior);}

 if (active(log_sigma2c)){
 likeval(9)=1./(2*square(cvc3))*norm2(log_sigma2c-log_s2priorC);}

 if (active(log_L50c)){
 likeval(10)=1./(2*square(cvc1))*norm2(log_L50c-log_L50priorC);}

 if (active(log_Linf)){
 likeval(11)=1./(2*square(cv_Par_bio(1)))*square(log_Linf-log_Linf_prior);}

 if (active(log_k)){
 likeval(12)=1./(2*square(cv_Par_bio(2)))*square(log_k-log_k_prior);}

 if(active(log_Lr)){
 likeval(13)=1./(2*square(cv_Par_bio(3)))*square(log_Lr-log_Lr_prior);}

 if(active(log_sr)){
 likeval(14)=1./(2*square(cv_Par_bio(4)))*square(log_sr-log_sr_prior);}

 if(active(log_beta)){
 likeval(15)=1./(2*square(cv_Par_bio(5)))*square(log_beta-log_beta_prior);}

 if(active(log_M)){
 likeval(16)=1./(2*square(cv_Par_bio(6)))*square(log_M-log_M_prior);}

 if(active(log_h)){
 likeval(17)=1./(2*square(cv_Par_bio(7)))*square(log_h-log_h_prior);}

 if (active(log_qcru)){
 likeval(18)=1./(2*square(cv_log_qprior))*norm2(log_qcru-log_qprior);}

 if (active(log_F)){
 penalty+=1000*norm2(log_F-mean(log_F));}


 f=(sum(likeval)+penalty);

 if(last_phase){
 f=sum(likeval)+200*square(log_desv_Rt(ntime));}


FUNCTION  Eval_CTP

//-----------------------------------------------------------------
  for (int i=1;i<=npbr;i++){ // ciclo de PBR

  Np=N(ntime);
  Sp=S(ntime);
  Fpbr=F(ntime)*pbr(i);//
  Zpbr=Fpbr+M;

  for (int j=1;j<=ntime_sim;j++){ // ciclo de años

  if(j<=minedad){
  Rp=(alfa_sr*BD(ntime-minedad+j)/(beta_sr+BD(ntime-minedad+j)));
  }
  if(j>minedad){
  Rp=(alfa_sr*Bp(i,j-minedad)/(beta_sr+Bp(i,j-minedad)));
  }
  Np=elem_prod(Np,Sp)*T+pre*Rp; //
  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dts*Zpbr)),elem_prod(msex,Wmed)));

  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  Yp(i,j)=sum(elem_prod(CTPp,Wmed));
  Sp=exp(-1.*Zpbr);
  }
 
 CBAtmas1(i)=Yp(i,1);// recomendacion CBA(t+1)
 BDBD0LP(i)=Bp(i,ntime_sim)/SSBo;// %B0LP
 BDBD0tmas2(i)=Bp(i,2)/SSBo;// %B0LP

 }



REPORT_SECTION

 report << "Yrs" << endl;
 report << yrs << endl;
 report << "Bcrucero" << endl;
 report << Bcru << endl;
 report << pred_Bcru << endl;
 report << "CPUE" << endl;
 report << CPUE << endl;
 report << pred_CPUE << endl;

FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;



 ofstream print_R("for_R.rep");

 print_R << "Yrs" << endl;
 print_R << yrs << endl;
 print_R << "Bcrucero" << endl;
 print_R << Bcru << endl;
 print_R << pred_Bcru << endl;
 print_R << "CPUE" << endl;
 print_R << CPUE << endl;
 print_R << pred_CPUE << endl;
 print_R << "Desembarques" << endl;
 print_R << Desemb << endl;
 print_R << pred_Desemb << endl;
 print_R << "Lmed_flo" << endl;
 print_R << nanos_ftc << endl;
 print_R << Lmed_obs << endl;
 print_R << Lmed_pred << endl;
 print_R << "Lmed_srv" << endl;
 print_R << nanos_fts1 << endl;
 print_R << Lmed_obsc << endl;
 print_R << Lmed_predc << endl;
 print_R << "Bio_reprod" << endl;
 print_R << BD << endl;
 print_R << BD.sd << endl;
 print_R << "Bio_tot" << endl;
 print_R << BT << endl;
 print_R << "Bio_explot" << endl;
 print_R << Bexplot << endl;
 print_R << Bexplot.sd << endl;
 print_R << "Reclus" << endl;
 print_R << Reclutas<< endl;
 print_R << Rpred<< endl;
 print_R << Reclutas.sd<< endl;

 print_R << "Dev_log_R" << endl;
 print_R << log_desv_Rt << endl;

 print_R << "Pre_r" << endl;
 print_R << pre << endl;


 print_R << "Mort_F " << endl;
 print_R <<Fmort << endl;
 print_R <<Fmort.sd << endl;


 print_R<<"Tallas"<<endl;
 print_R<<Tallas<<endl;

 print_R<<"Abundancia_talla"<<endl;
 print_R<<N<<endl;

 print_R<<"Sel_f"<<endl;
 print_R<<Sel<<endl;

 print_R<<"Sel_srv"<<endl;
 print_R<<Selc<<endl;

 print_R<<"msex"<<endl;
 print_R<<msex<<endl;

 print_R<<"PesoL"<<endl;
 print_R<<Wmed<<endl;

 print_R<<"h"<<endl;
 print_R<<h<<endl;

 print_R<<"dt"<<endl;
 print_R<<dts<<endl;

 print_R<<"q_survey"<<endl;
 print_R<<exp(log_qcru)<<endl;

 print_R<<"lag"<<endl;
 print_R<<minedad<<endl;

 print_R<<"Lmuda"<<endl;
 print_R<<Lesp<<endl;

 print_R<<"sLmuda"<<endl;
 print_R<<sigmaL<<endl;

 print_R << "Frecs_capt_obs" << endl;
 print_R << pobs<< endl;
 print_R << "Frecs_capt_pred" << endl;
 print_R << ppred<< endl;
 print_R << "Frecs_srv_obs" << endl;
 print_R << pobsc<< endl;
 print_R << "Frecs_srv_pred" << endl;
 print_R << ppredc<< endl;

 print_R << "Frec_marg_flo" << endl;
 print_R << colsum(pobs)<< endl;
 print_R << colsum(ppred)<< endl;

 print_R << "Frec_marg_srv" << endl;
 print_R << colsum(pobsc)<< endl;
 print_R << colsum(ppredc)<< endl;
 print_R << "B0" << endl;
 print_R << SSBo << endl;
 print_R << "R0" << endl;
 print_R << exp(log_Rmed) << endl;
 print_R << "RPRLP " << endl;
 print_R << RPRlp << endl;
 print_R << RPRlp.sd << endl;
 print_R << "Lmed" << endl;
 print_R << Lesp << endl;
 print_R << sigmaL << endl;
 print_R << "Tallas_pre" << endl;
 print_R << Tallas<< endl;
 print_R << pre<< endl;
 print_R << "Ptrans" << endl;
 print_R << T << endl;
 print_R << "Lr_Sr_beta" << endl;
 print_R <<Lr<<" "<<sr<<" "<<beta<< endl;
 print_R << "Linf_k_M_h" << endl;
 print_R <<Linf<<" "<<k<<" "<<M<<" "<<h<< endl;
 print_R << "alfaSR_betaSR" << endl;
 print_R <<alfa_sr<<" "<<beta_sr<< endl;


// ESTIMA nm y CV

  suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=N_ftc;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 print_R << "nm_flota_crucero" <<endl;
 print_R <<pow(nm1,1/cuenta1)<< endl;


 suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=N_fts1;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppredc(i),1-ppredc(i)));
      suma2=norm2(pobsc(i)-ppredc(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 print_R <<pow(nm1,1/cuenta1)<< endl;

 print_R << "cv_capt" << endl;
 print_R << cv_index(1) << endl;
 print_R << "cv_CPUE" << endl;
 print_R << cv_index(2) << endl;
 print_R << "cv_cru" << endl;
 print_R << cv_index(3) << endl;
 print_R << "Fmult" << endl;
 print_R << pbr << endl;
 print_R << "BRep_proy" << endl;
 print_R << Bp << endl;
 print_R << "Capt_proy" << endl;
 print_R << Yp << endl;
 print_R << "Likeval" << endl;
 print_R <<  likeval << endl;
 print_R << "MaxGrad" << endl;
 print_R <<objective_function_value::pobjfun->gmax<<endl;
 print_R << "FunObj" << endl;
 print_R <<  f << endl;
 print_R << "N0" << endl;
 print_R << No << endl;

