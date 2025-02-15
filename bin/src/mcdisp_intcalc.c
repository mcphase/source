//calculate intensities for given energy

//***********************************************************************//
// initialize intensity calculation for a specific Q vector (calls d*1calc)
//***********************************************************************
// returns 1 on success and zero on failure
//***********************************************************************
void intcalc_ini(inimcdis & ini,par & inputpars,mdcf & md,int do_Erefine,double & epsilon,int do_verbose,int do_gobeyond,int calc_rixs,int do_phonon,Vector & hkl,int qcounter)
{int i,j,k,l,m,jmin,i1,j1,tn; Vector qijk(1,3);double QQ;
  hkl2ijk(qijk,hkl, inputpars.cs.abc);
    QQ=Norm(qijk);

// polarisation factor
    ComplexMatrix pol(1,3,1,3);pol=0;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qijk(i)*qijk(j)/QQ/QQ;//(qijk*qijk);
    }}

 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
  char filename[MAXNOFCHARINLINE];
 float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
 if(do_verbose==1) printf("#initializing intensity calculation\n");
  double Gamman,gamma,gammaP; 
  complex<double> imaginary(0,1),dP1;
  FILE * fin; 
  Vector mf(1,ini.nofcomponents);
      int mqdim=3;if(calc_rixs)mqdim=9;
      ComplexVector L1(1,3),S1(1,3),P1(1,3),mq1_dip(1,mqdim),mq1(1,mqdim);

  // Precalculate values of Debye-Waller and Form Factors for this Q-vector to save calls to (*inputpars.jjj[ion]).* functions
  double DBWF[inputpars.cs.nofatoms+1];
  complex <double> * SL;SL=new complex<double>[inputpars.cs.nofatoms+1];
  for(l=1;l<=inputpars.cs.nofatoms;++l) {
     DBWF[l] = (*inputpars.jjj[l]).debyewallerfactor(QQ);
     SL[l] = complex <double> ((*inputpars.jjj[l]).SLR,(*inputpars.jjj[l]).SLI);
  }

  snprintf(filename,MAXNOFCHARINLINE,"./results/%smcdisp.trs",ini.prefix);  
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.cs.nofatoms;++l){
  if((fin = fopen(filename,"rb"))==NULL){snprintf(filename,MAXNOFCHARINLINE,"./results/mcdisp.trs");
                  fin = fopen_errchk(filename,"rb");}

      // do calculation for atom s=(ijkl)
      for(int ll=1;ll<=ini.nofcomponents;++ll){mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} 
                                               //mf ... exchange field vector of atom s

  if(do_Erefine) // clear chi0 matrices (put sign of qcounter negative
  {(*inputpars.jjj[l]).chi0(md.chi0pointer(i,j,k,l),ini.emin, fabs(epsilon/2),md.nofEstps,epsilon,
                  qijk,-qcounter,nn[6],ini.T,mf,ini.Hext, md.est(i,j,k,l),i,j,k,l);}

  jmin=0;
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3]&&l==(int)nn[4])
    {tn=(int)nn[5];++jmin;  
    // calculate delta(single ion excitation energy), 
    // Malphabeta(transition matrix elements)

//      fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
//      if(nn[6]<SMALL_QUASIELASTIC_ENERGY){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");}
      
        j1=(*inputpars.jjj[l]).transitionnumber; // try calculation for transition  j
       
      (*inputpars.jjj[l]).transitionnumber=tn; // try calculation for transition  tn
      if(do_verbose==1)(*inputpars.jjj[l]).transitionnumber=-tn;

  if(do_Erefine) // fill chi0 matrices 
  {(*inputpars.jjj[l]).chi0(md.chi0pointer(i,j,k,l),ini.emin, fabs(epsilon/2), md.nofEstps,epsilon,
                  qijk,qcounter,nn[6],ini.T,mf,ini.Hext, md.est(i,j,k,l),i,j,k,l);}

  if(calc_rixs)// RIXS
      {if((*inputpars.jjj[l]).drixs1calc(qijk,ini.T,mq1_dip,md.est(i,j,k,l))!=0)
       // calculate <-|Rijomega|+> see haverkort paper: transition operator for RIXS
       {if(do_verbose)printf("#calculating RIXS intensity for ion %s\n",(*inputpars.jjj[l]).sipffilename);
       (*inputpars.jjj[l]).FF_type=4; // put FFTYPE to RIXS to indicate that this is implemented
       }
       else{ if(do_verbose)printf("#warning mcdisp - function drixs1 not implemented for single ion module of ion %s, no intensity from this ion\n",(*inputpars.jjj[l]).sipffilename);
           mq1_dip=0;mq1_dip(1)= complex <double> (1e-10,0.0);(*inputpars.jjj[l]).FF_type=1;
           }
      }else
      {// neutron intensities
      // try dipole approximation for this ion
      // dipole approx: <M(Q)>=<L>*FL(Q)+2*<S>*FS(Q) with FL(Q)=(j0+j2) and FS(Q)=j0
      if((*inputpars.jjj[l]).dL1calc(ini.T,mf,ini.Hext,L1,md.est(i,j,k,l))!=0&&
         (*inputpars.jjj[l]).dS1calc(ini.T,mf,ini.Hext,S1,md.est(i,j,k,l))!=0)
         {(*inputpars.jjj[l]).FF_type=3;mq1_dip=(*inputpars.jjj[l]).F(-QQ)*L1+2.0*(*inputpars.jjj[l]).F(QQ)*S1;
         if(do_verbose)printf("#L1=%g FL(Q)=%g S1=%g FS(Q)=%g\n",Norm2(L1),(*inputpars.jjj[l]).F(-QQ),Norm2(S1),(*inputpars.jjj[l]).F(QQ));}
      else  {ComplexVector m1(1,3); // try dipole approximation using dm1calc
                                    // gJ=0 dipole approx: <M(Q)>=<M>*F(Q) with F(Q)=j0  
                                    // gJ>0 dipole approx: <M(Q)>=<M>*F(Q) with F(Q)=j0-(1-2/gJ)j2                
             if((*inputpars.jjj[l]).dm1calc(ini.T,mf,ini.Hext,m1,md.est(i,j,k,l))!=0)
             {mq1_dip=m1*(*inputpars.jjj[l]).F(QQ);(*inputpars.jjj[l]).FF_type=2; 
             }
             else {if(do_verbose)printf("#warning mcdisp - functions dmq1,dm1calc,dL1calcd,S1calc not implemented for single ion module of ion %s, no magnetic neutron intensity from this ion\n",(*inputpars.jjj[l]).sipffilename);
                   mq1_dip=0;mq1_dip(1)= complex <double> (1e-10,0.0);(*inputpars.jjj[l]).FF_type=1;
                  }
             }
      // try if going beyond dip approximation for formfactor works
     if(do_gobeyond){
        if((*inputpars.jjj[l]).dMQ1calc(qijk,ini.T,mq1,nn[6],md.est(i,j,k,l))!=0)// calculate <-|M(Q)|+>
        {if(do_verbose)printf("#going beyond dipole approx for ion %s Nu=%g Nu_dip=%g \n",(*inputpars.jjj[l]).sipffilename,Norm2(mq1),Norm2(mq1_dip));
         (*inputpars.jjj[l]).FF_type*=-1; // put FFTYPE negative to indicate that going beyond works
        }
        else{ if(do_verbose)printf("#warning mcdisp - function dmq1 not implemented for single ion module of ion %s, only doing dipolar intensity\n",(*inputpars.jjj[l]).sipffilename);
             mq1=mq1_dip;
            }
         mq1=pol*mq1;// take into account polarisation factor ... only Mper scatters
                  } // do_gobeyond
         if(ini.outS<5)mq1_dip=pol*mq1_dip;// take into account polarisation factor ... only Mper scatters

      // try to calculate phonon intensity
      if(do_phonon){
                  if((*inputpars.jjj[l]).dP1calc(ini.T,mf,ini.Hext,P1,md.est(i,j,k,l))==0)
                  {dP1=0;gammaP=SMALL_GAMMA;}
                  else
                  {dP1=0;for(int pp1=1;pp1<=3;++pp1)dP1+=DBWF[l]*SL[l]*(P1(pp1)*qijk(pp1));
                   gammaP=real(dP1*conj(dP1));if(gammaP>SMALL_NORM)dP1/=sqrt(gammaP);
                  }
                   }
     } // if calc_rixs

      if(do_gobeyond){mq1*=DBWF[l]; // multiply Debye Waller factor
                      Gamman=Norm2(mq1);if(Gamman>SMALL_NORM)mq1/=sqrt(Gamman);
                      }
      mq1_dip*=DBWF[l]; // multiply Debye Waller factor
      gamma=Norm2(mq1_dip);if(gamma>SMALL_NORM)mq1_dip/=sqrt(gamma);

       (*inputpars.jjj[l]).transitionnumber=j1; // put back transition number for 1st transition
       j1=md.baseindex(i,j,k,l,jmin); 
      
         // treat correctly case for  energy loss
	 if (nn[6]<0){if(do_gobeyond)mq1=mq1.Conjugate();
                      if(do_phonon)dP1=conj(dP1);
                      mq1_dip=mq1_dip.Conjugate();}
                            // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
        if (nn[6]>SMALL_QUASIELASTIC_ENERGY)
	    {if(do_gobeyond)md.sqrt_Gamma(i,j,k)(j1)=sqrt(Gamman);
             if(do_phonon)md.sqrt_GammaP(i,j,k)(j1)=sqrt(gammaP);
             md.sqrt_Gamma_dip(i,j,k)(j1)=sqrt(gamma);}
	    else if (nn[6]<-SMALL_QUASIELASTIC_ENERGY)
            {if(do_gobeyond)md.sqrt_Gamma(i,j,k)(j1)=imaginary*sqrt(Gamman);
             if(do_phonon)md.sqrt_GammaP(i,j,k)(j1)=imaginary*sqrt(gammaP);
             md.sqrt_Gamma_dip(i,j,k)(j1)=imaginary*sqrt(gamma);}
	    else
	    { //quasielastic line needs gamma=SMALL_QUASIELASTIC_ENERGY .... because Mijkl and therefore gamma have been set to 
              // wn/kT instead of wn-wn'=SMALL_QUASIELASTIC_ENERGY*wn/kT (in jjjpar.cpp -mdcalc routines)
	      //set fix delta but keep sign
	      if (nn[6]<0){
  			     if(do_gobeyond)md.sqrt_Gamma(i,j,k)(j1)=imaginary*sqrt(SMALL_QUASIELASTIC_ENERGY*Gamman);
  			     if(do_phonon)md.sqrt_GammaP(i,j,k)(j1)=imaginary*sqrt(SMALL_QUASIELASTIC_ENERGY*gammaP);
  			     md.sqrt_Gamma_dip(i,j,k)(j1)=imaginary*sqrt(SMALL_QUASIELASTIC_ENERGY*gamma);
                           }
	      else        {  if(do_gobeyond)md.sqrt_Gamma(i,j,k)(j1)=sqrt(SMALL_QUASIELASTIC_ENERGY*Gamman);
  			     if(do_phonon)md.sqrt_GammaP(i,j,k)(j1)=sqrt(SMALL_QUASIELASTIC_ENERGY*gammaP);
  			     md.sqrt_Gamma_dip(i,j,k)(j1)=sqrt(SMALL_QUASIELASTIC_ENERGY*gamma);
                          }
	    }	                   
        if(do_gobeyond)for(m=1;m<=mqdim;++m){md.dMQs(i,j,k)(mqdim*(j1-1)+m)=mq1(m);}
        if(do_phonon)md.dPs(i,j,k)(1*(j1-1)+1)=dP1;
    
        for(m=1;m<=mqdim;++m){md.dMQ_dips(i,j,k)(mqdim*(j1-1)+m)=mq1_dip(m);}       
    }}}
    fclose(fin);
  }}}}
 delete []SL;
}

#include <time.h>

//**************************************************************************/
#ifdef _THREADS
#if defined  (__linux__) || defined (__APPLE__)
void *intcalc_approx(void *input)
#else
DWORD WINAPI intcalc_approx(void *input)
#endif
#else
double intcalc_approx(ComplexMatrix & chi,ComplexMatrix & chibey,ComplexMatrix & chiPhon,
        Matrix & pol, double & intensitybey,double & intensityP,
        mfcf & qee_real,mfcf & qee_imag,ComplexMatrix & Echargedensity,
        mfcf & qsd_real,mfcf & qsd_imag,ComplexMatrix & Espindensity,
        mfcf & qod_real,mfcf & qod_imag,ComplexMatrix & Eorbmomdensity,
        mfcf & qep_real,mfcf & qep_imag,ComplexMatrix & Ephonon,
        mfcf & qem_real,mfcf & qem_imag,ComplexMatrix & Emagmom,
        mfcf & qes_real,mfcf & qes_imag,ComplexMatrix & Espin,
        mfcf & qel_real,mfcf & qel_imag,ComplexMatrix & Eorbmom,
        int dimA, const ComplexMatrix &Tau, int level,double en, const inimcdis & ini,
       const par & inputpars,Vector & hkl, mdcf & md,int do_verbose,int calc_rixs,int do_phonon,double & QQ)
#endif
{//calculates approximate intensity for energylevel i - according to chapter 8.2 mcphas manual

#ifdef _THREADS
   intcalcapr_input *myinput; myinput = (intcalcapr_input *)input;
   int thread_id = myinput->thread_id;
   double intensitybey = myinput->intensitybey;
   double intensityP = myinput->intensityP;
   #define chi (*thrdat.chi[thread_id])
   #define chibey (*thrdat.chibey[thread_id])
   #define chiPhon (*thrdat.chiPhon[thread_id])
   #define pol (*thrdat.pol[thread_id])
   #define qee_real (*thrdat.qee_real[thread_id])
   #define qee_imag (*thrdat.qee_imag[thread_id])
   #define qsd_real (*thrdat.qsd_real[thread_id])
   #define qsd_imag (*thrdat.qsd_imag[thread_id])
   #define qod_real (*thrdat.qod_real[thread_id])
   #define qod_imag (*thrdat.qod_imag[thread_id])
   #define qep_real (*thrdat.qep_real[thread_id])
   #define qep_imag (*thrdat.qep_imag[thread_id])
   #define qem_real (*thrdat.qem_real[thread_id])
   #define qem_imag (*thrdat.qem_imag[thread_id])
   #define qes_real (*thrdat.qes_real[thread_id])
   #define qes_imag (*thrdat.qes_imag[thread_id])
   #define qel_real (*thrdat.qel_real[thread_id])
   #define qel_imag (*thrdat.qel_imag[thread_id])
   #define Echargedensity (*thrdat.Echargedensity[thread_id])
   #define Espindensity (*thrdat.Espindensity[thread_id])
   #define Eorbmomdensity (*thrdat.Eorbmomdensity[thread_id])
   #define Ephonon (*thrdat.Ephonon[thread_id])
   #define Emagmom (*thrdat.Emagmom[thread_id])
   #define Espin (*thrdat.Espin[thread_id])
   #define Eorbmom (*thrdat.Eorbmom[thread_id])
   int level =  myinput->level;//, dimA = myinput->dimA, do_verbose = myinput->do_verbose;
   #define Tau (*thrdat.Tau[thread_id])
   double en = myinput->En; 
   int calc_rixs = myinput->calc_rixs;
   int do_phonon = myinput->do_phonon;
   #define ini (*thrdat.ini[thread_id])
   #define inputpars (*thrdat.inputpars[thread_id])
   #define hkl thrdat.hkl
   #define md (*thrdat.md[thread_id])
   double QQ;
#endif
 int mqdim=3;if(calc_rixs)mqdim=9;
 int i,j,i1,j1,k1,l1,t1,i2,j2,k2,l2,t2,s,ss,b,bb;
 double intensity=1.2; 
 double ki,kf;
 complex <double> sumS;
 //complex <double> chileft;
 //complex <double> chileftbey;
 Vector qijk(1,3);
  hkl2ijk(qijk,hkl, (Vector&)inputpars.cs.abc);
 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
//   qijk(1)=hkl(1)*2*PI/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)*2*PI/inputpars.b;
//    qijk(3)=hkl(3)*2*PI/inputpars.c;
    QQ=Norm(qijk);
 // init eigenvector to zero
  if(ini.calculate_chargedensity_oscillation){qee_real.clear();qee_imag.clear();}
  if(ini.calculate_spindensity_oscillation){qsd_real.clear();qsd_imag.clear();}
  if(ini.calculate_orbmomdensity_oscillation){qod_real.clear();qod_imag.clear();}
  if(ini.calculate_phonon_oscillation){qep_real.clear();qep_imag.clear();}
  if(ini.calculate_magmoment_oscillation){qem_real.clear();qem_imag.clear();} 
  if(ini.calculate_spinmoment_oscillation){qes_real.clear();qes_imag.clear();}
  if(ini.calculate_orbmoment_oscillation){qel_real.clear();qel_imag.clear();}


 // Added code to re-use previously calculated values of sqrt(gamma)*U and conj(U)*conj(sqrt(gamma)). mdl 110705
 int maxb=-1,bval,/*ncel=-1,*/nval; complex<double> defval(-0.1,0.); md.ncel=-1;// int tval;
 for(i2=1;i2<=ini.mf.na();++i2) for(j2=1;j2<=ini.mf.nb();++j2) for(k2=1;k2<=ini.mf.nc();++k2) { 
   bval=md.baseindex_max(i2,j2,k2); if(bval>maxb) maxb=bval; 
   nval=md.in(i2,j2,k2); if(nval>md.ncel) md.ncel=nval; } md.ncel++;
   if(maxb==md.nofcomponents) maxb++; // To ensure matrix is not square so overloaded operator ComplexMatrix=(complex<double>) sets all elements to 1e-16 not just diagonal.
   if(md.Ug==0) { md.Ug = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.Ug[i]=0; }
   if(md.gU==0) { md.gU = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.gU[i]=0; }
   if(md.bUg==0) { md.bUg = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.bUg[i]=0; }
   if(md.bgU==0) { md.bgU = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.bgU[i]=0; }
   if(md.PUg==0) { md.PUg = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.PUg[i]=0; }
   if(md.PgU==0) { md.PgU = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.PgU[i]=0; }
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){ int in1=md.in(i1,j1,k1);
   if(md.gU[in1]==0) { md.gU[in1] = new ComplexMatrix(1,mqdim,1,maxb); *md.gU[in1]=defval; }
   if(md.Ug[in1]==0) { md.Ug[in1] = new ComplexMatrix(1,mqdim,1,maxb); *md.Ug[in1]=defval; }
   if(md.bgU[in1]==0) { md.bgU[in1] = new ComplexMatrix(1,mqdim,1,maxb);*md.bgU[in1]=defval; }
   if(md.bUg[in1]==0) { md.bUg[in1] = new ComplexMatrix(1,mqdim,1,maxb);*md.bUg[in1]=defval; }
   if(md.PgU[in1]==0) { md.PgU[in1] = new ComplexMatrix(1,1,1,maxb);*md.PgU[in1]=defval; }
   if(md.PUg[in1]==0) { md.PUg[in1] = new ComplexMatrix(1,1,1,maxb);*md.PUg[in1]=defval; } }}}

chi=0;chibey=0;chiPhon=0;
int ssm1,in1,in2;

// determine chi
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
   in1=md.in(i1,j1,k1);      
 for(l1=1;l1<=md.nofatoms;++l1){
 for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      s=index_s(i1,j1,k1,l1,t1,md,ini);// sm1=s-1;//s3=sm1*mqdim;
      b=md.baseindex(i1,j1,k1,l1,t1);
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
    in2=md.in(i2,j2,k2);
  for(l2=1;l2<=md.nofatoms;++l2){
  for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      ss=index_s(i2,j2,k2,l2,t2,md,ini);ssm1=ss-1;
      bb=md.baseindex(i2,j2,k2,l2,t2);
   
      if(do_phonon)
      { if((*md.PgU[in1])(1,b)==defval)  (*md.PgU[in1])(1,b)  = conj(md.sqrt_GammaP(i1,j1,k1)(b))
                                                                 * md.dPs(i1,j1,k1)(b);
        if((*md.PUg[in2])(1,bb)==defval) (*md.PUg[in2])(1,bb) = conj(md.dPs(i2,j2,k2)(bb))
                                                                 * md.sqrt_GammaP(i2,j2,k2)(bb);
         chiPhon(1,1)+= PI * (*md.PgU[in1])(1,b) * Tau(s,level) * en * conj(Tau(ss,level)) * (*md.PUg[in2])(1,bb);          
              } // i,j,do_phonon

      if(intensitybey>0)
      {for(j=1;j<=mqdim;++j){
         if((*md.bUg[in2])(j,bb)==defval) (*md.bUg[in2])(j,bb) = conj(md.dMQs(i2,j2,k2)((bb-1)*mqdim+j))
                                                                 * md.sqrt_Gamma(i2,j2,k2)(bb);
         for(i=1;i<=mqdim;++i){
        if((*md.bgU[in1])(i,b)==defval)  (*md.bgU[in1])(i,b)  = conj(md.sqrt_Gamma(i1,j1,k1)(b))
                                                                 * md.dMQs(i1,j1,k1)((b-1)*mqdim+i);
                       
         chibey(i,j)+= PI * (*md.bgU[in1])(i,b) * Tau(s,level) * en * conj(Tau(ss,level)) * (*md.bUg[in2])(j,bb);         
      }}} // i,j,intensitybey

      for(j=1;j<=mqdim;++j){
        if((*md.Ug[in2])(j,bb)==defval) (*md.Ug[in2])(j,bb) = conj(md.dMQ_dips(i2,j2,k2)((bb-1)*mqdim+j))
                                                                 * md.sqrt_Gamma_dip(i2,j2,k2)(bb);
         for(i=1;i<=mqdim;++i){
        if((*md.gU[in1])(i,b)==defval)  (*md.gU[in1])(i,b)  = conj(md.sqrt_Gamma_dip(i1,j1,k1)(b))
                                                                 * md.dMQ_dips(i1,j1,k1)((b-1)*mqdim+i);
                     
         chi(i,j)+= PI * (*md.gU[in1])(i,b) * Tau(s,level) * en * conj(Tau(ss,level)) * (*md.Ug[in2])(j,bb);         
      }}
                     
    for(j=1;j<=md.nofcomponents;++j){ 
     if(ssm1*md.nofcomponents+j==1){if(ini.calculate_chargedensity_oscillation)for(i=1;i<=CHARGEDENS_EV_DIM;++i)
                                        {qee_real.mf(i1,j1,k1)(CHARGEDENS_EV_DIM*(l1-1)+i)+=real(Echargedensity(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qee_imag.mf(i1,j1,k1)(CHARGEDENS_EV_DIM*(l1-1)+i)+=imag(Echargedensity(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_spindensity_oscillation)for(i=1;i<=3*SPINDENS_EV_DIM;++i)
                                        {qsd_real.mf(i1,j1,k1)(3*SPINDENS_EV_DIM*(l1-1)+i)+=real(Espindensity(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qsd_imag.mf(i1,j1,k1)(3*SPINDENS_EV_DIM*(l1-1)+i)+=imag(Espindensity(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_orbmomdensity_oscillation)for(i=1;i<=3*ORBMOMDENS_EV_DIM;++i)
                                        {qod_real.mf(i1,j1,k1)(3*ORBMOMDENS_EV_DIM*(l1-1)+i)+=real(Eorbmomdensity(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qod_imag.mf(i1,j1,k1)(3*ORBMOMDENS_EV_DIM*(l1-1)+i)+=imag(Eorbmomdensity(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_phonon_oscillation)for(i=1;i<=PHONON_EV_DIM;++i)
                                        {qep_real.mf(i1,j1,k1)(PHONON_EV_DIM*(l1-1)+i)+=real(Ephonon(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qep_imag.mf(i1,j1,k1)(PHONON_EV_DIM*(l1-1)+i)+=imag(Ephonon(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_magmoment_oscillation)for(i=1;i<=MAGMOM_EV_DIM;++i)
                                        {qem_real.mf(i1,j1,k1)(MAGMOM_EV_DIM*(l1-1)+i)+=real(Emagmom(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qem_imag.mf(i1,j1,k1)(MAGMOM_EV_DIM*(l1-1)+i)+=imag(Emagmom(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        // check consistency:
                                        // printf("Emagmom(s=%i,alpha=%i)= %g %+g i \n",s,i,
                                        //          real(Emagmom(s,i)),imag(Emagmom(s,i)));
                                         }//printf("at E=%g\n",en);  
                                     if(ini.calculate_spinmoment_oscillation)for(i=1;i<=SPIN_EV_DIM;++i)
                                        {qes_real.mf(i1,j1,k1)(SPIN_EV_DIM*(l1-1)+i)+=real(Espin(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qes_imag.mf(i1,j1,k1)(SPIN_EV_DIM*(l1-1)+i)+=imag(Espin(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_orbmoment_oscillation)for(i=1;i<=ORBMOM_EV_DIM;++i)
                                        {qel_real.mf(i1,j1,k1)(ORBMOM_EV_DIM*(l1-1)+i)+=real(Eorbmom(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qel_imag.mf(i1,j1,k1)(ORBMOM_EV_DIM*(l1-1)+i)+=imag(Eorbmom(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     }
                                    }
   }}}}
  }}}
 }}}


  //complex<double> im(0,1.0);

 // determine dsigma in barns per cryst unit cell !
 //divide by number of crystallographic unit cells  (ini.mf.n()) in magnetic unit cell


   double bose;
   if (fabs(en/KB/ini.T)>SMALL_QUASIELASTIC_ENERGY*0.1)
   {bose=1.0/(1.0-exp(-en*(1.0/KB/ini.T)));  
   }else{//quasielastic needs special treatment
         //if(fabs(en/KB/ini.T)>1e-30){bose=ini.T*KB/en;}
         //else{bose=1e30;} .... this is no good
         //(problem: quasielastic intensity depends on value of SMALL_QUASIELASTIC_ENERGY !!)
	 // in principle this SMALL_QUASIELASTIC_ENERGY in denominator has to cancel with epsilon
	 // in population matrix Mijkl(i.e. gamma) ... therefore we skip it:
	 // (for small energies delta_s the md.sqrt_gamma has been set = sqr(SMALL_QUASIELASTIC_ENERGY*gamma) and this is
	 // inserted into the calculation of chi above)
//         bose=ini.T*KB/(SMALL_QUASIELASTIC_ENERGY*0.1); // changed by MR 9.3.11 because factor omegar introduced in DMD formulas
//   bose=ini.T*KB;
     bose=ini.T*KB/en;  // replaced by the following line MR 9.3.11
   }
//  bose=fabs(bose);// this is to correctly consider omegar/|omegar| in the formula for chi'' ... introduced 2.4.10 MR
                    // deleted by MR 9.3.11 because new derivation of DMD gives factor omegar
                    //                      only (introduced above in bose) and
                    //                      we remove normalisation of tau in eigenvalueGeneral
                    //                      so that Tau is in units of sqrt of meV^-1


if (calc_rixs){// use 1-9 components of chi to store result !!! (other components do not count          
               chi*=bose/(double)ini.mf.n();
              } 
       else{double factor=0.5*bose*3.65/4.0/PI/(double)ini.mf.n()/PI/2.0;
            chi*=factor;
            if(ini.outS<5){
            sumS=Trace(chi);intensity=fabs(real(sumS)); // if polarisation factor is already in Mperp we need Trace instead of of sum
                      if (real(sumS)<-0.1){fprintf(stderr,"ERROR mcdisp: dipolar approx intensity %g negative,E=%g, bose=%g\n",real(sumS),en,bose);exit(1);}
                      if (fabs(imag(sumS))>0.1){fprintf(stderr,"ERROR mcdisp: dipolar approx intensity %g %+g iimaginary\n",real(sumS),imag(sumS));exit(1);}
                          }else
                          {intensity=-1;} // no dipole intensity output if scattering function is not including pol factor !
            if(intensitybey>0){chibey*=factor;
                               sumS=Trace(chibey);intensitybey=fabs(real(sumS));// if polarisation factor is already in Mperp we need Trace instead of of sum
                               intensitybey=fabs(real(sumS)); if (real(sumS)<-0.1){fprintf(stderr,"ERROR mcdisp: intensity in beyond dipolar approx formalism %g negative,E=%g, bose=%g\n\n",real(sumS),en,bose);exit(1);}
                                                   if (fabs(imag(sumS))>0.1){fprintf(stderr,"ERROR mcdisp: intensity  in beyond dipolar approx formalism %g %+g iimaginary\n",real(sumS),imag(sumS));exit(1);}
                              }
            if(do_phonon){sumS=chiPhon(1,1)/PI/2.0/(double)ini.mf.n();sumS*=bose;// printf("chiPhon= %g i %g ",real(chiPhon(1,1)),imag(chiPhon(1,1)));
                          intensityP=fabs(real(sumS)); if (real(sumS)<-0.1){fprintf(stderr,"ERROR mcdisp: phonon intensity in approx formalism %g negative,E=%g, bose=%g\n\n",real(sumS),en,bose);exit(1);}
                                                   if (fabs(imag(sumS))>0.1){fprintf(stderr,"ERROR mcdisp: phonon  intensity  in  approx formalism %g %+g iimaginary\n",real(sumS),imag(sumS));exit(1);}
                         }                              
// here should be entered factor  k/k' 
if (ini.ki==0)
{if (ini.kf*ini.kf+0.4811*en<0)
 {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with kf=const=%g/A at (%g,%g,%g)\n",en,ini.kf,hkl(1),hkl(2),hkl(3));
  intensity=0;
  intensitybey=0;
  intensityP=0;
 }
 else
 { 
 ki=sqrt(ini.kf*ini.kf+0.4811*en);
 if(intensity>0)intensity*=ini.kf/ki;
 if(intensitybey>0)intensitybey*=ini.kf/ki;
 if(intensityP>0)intensityP*=ini.kf/ki;
 }
}
else
{if (ini.ki*ini.ki-0.4811*en<0)
 {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)\n",en,ini.ki,hkl(1),hkl(2),hkl(3));
    intensity=0;
    intensitybey=0;
    intensityP=0;
 }
 else
 {kf=sqrt(ini.ki*ini.ki-0.4811*en);
  if(intensity>0)intensity*=kf/ini.ki;
  if(intensitybey>0)intensitybey*=kf/ini.ki;
  if(intensityP>0)intensityP*=kf/ini.ki;
 }
}
               } // if calc_rixs
#ifdef _THREADS
//printf("intcalc: %i %g :",thread_id,intensitybey);
myinput->intensity=intensity;
myinput->intensitybey=intensitybey;
myinput->intensityP=intensityP;
myinput->QQ=QQ;
#undef md
#undef hkl
#undef Tau
#undef ini
#undef inputpars
#undef chi
#undef chibey
#undef chiPhon
#undef pol
#undef qee_real
#undef qee_imag
#undef qsd_real
#undef qsd_imag
#undef qod_real
#undef qod_imag
#undef qep_real
#undef qep_imag
#undef qem_real
#undef qem_imag
#undef qes_real
#undef qes_imag
#undef qel_real
#undef qel_imag
#undef Echargedensity
#undef Espindensity
#undef Eorbmomdensity
#undef Ephonon
#undef Emagmom
#undef Espin
#undef Eorbmom
MUTEX_LOCK(&mutex_loop);
thrdat.thread_id = thread_id;
EVENT_SIG(checkfinish);
MUTEX_UNLOCK(&mutex_loop);
#if defined  (__linux__) || defined (__APPLE__)
pthread_exit(NULL);
#else
return 0;
#endif
#else
return intensity;	
#endif
}




//**************************************************************************/
#ifdef _THREADSREFINE
#if defined  (__linux__) || defined (__APPLE__)
void *intcalc_Erefine(void *input)
#else
DWORD WINAPI intcalc_Erefine(void *input)
#endif
#else
double intcalc_Erefine(ComplexMatrix & ch, int Estp,inimcdis & ini,par & inputpars,jq & J,Vector & q,Vector & hkl,mdcf & md,int do_verbose,double epsilon)
#endif
{int i,j,i1,j1,k1,l1,i2,j2,k2,l2,s,ss,b,bb;
 double intensity=1.2;
 double QQ,ki,kf;

#ifdef _THREADSREFINE
   intcalcapr_input *myinput; myinput = (intcalcapr_input *)input;
   int thread_id = myinput->thread_id;
   //int dimA = myinput->dimA;
   int Estp  = myinput->Estp; 
   #define ini (*thrdat.ini[thread_id])
   #define inputpars (*thrdat.inputpars[thread_id])
   #define J (*thrdat.J[thread_id])
   #define q thrdat.q
   #define hkl thrdat.hkl
   #define md (*thrdat.md[thread_id])
   #define ch (*thrdat.ch[thread_id])
   int do_verbose = myinput->do_verbose;
   double epsilon = myinput->epsilon; 
#endif
 int dim=md.nofcomponents*md.nofatoms*md.n();
// determine chi
   ComplexMatrix chi(1,dim,1,dim);
   ComplexMatrix Ac(1,dim,1,dim);
   ComplexMatrix Acinv(1,dim,1,dim);
   ComplexMatrix Bc(1,dim,1,dim);
   ComplexMatrix cc1(1,md.nofcomponents,1,md.nofcomponents);
   Ac=0;Bc=0;
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
     for(l1=1;l1<=ini.nofatoms;++l1){
     s=md.inM(i1,j1,k1,l1)*md.nofcomponents;
    for(i=1;i<=md.nofcomponents;++i){
     Ac(s+i,s+i)=1; // set diagonal elements 1 (make Ac a unit matrix)
    for(j=1;j<=md.nofcomponents;++j){
     Bc(s+i,s+j)=(*md.chi0pointer(i1,j1,k1,l1)[Estp])(i,j);
     }}
      b=(md.baseindex(i1,j1,k1,l1,1)-1)*md.nofcomponents;

  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
       for(l2=1;l2<=ini.nofatoms;++l2){
       ss=md.inM(i2,j2,k2,l2)*md.nofcomponents;
      bb=(md.baseindex(i2,j2,k2,l2,1)-1)*md.nofcomponents;
     for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
      cc1(i,j)=J.mati(J.in(i1,j1,k1),J.in(i2,j2,k2))(b+i,bb+j);}}
      cc1=(*md.chi0pointer(i1,j1,k1,l1)[Estp])*cc1;
     for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
      Ac(s+i,ss+j)-=cc1(i,j);
    }}
   }}}}   
 }}}}

if(do_verbose){printf("diagonalising matrix A for Estep %i\n",Estp);//myPrintComplexMatrix(stdout,Ac); 
              }
 chi=Ac.Inverse()*Bc;

 // determine chi'' and S (bose factor)
   complex<double> im(0,1.0);
   double en=ini.emin+Estp*epsilon/2;
   complex<double> z(en,epsilon);
   ComplexMatrix S(1,dim,1,dim);
   complex<double> bose;
   bose=1.0/(1.0-exp(-z*(1.0/KB/ini.T)));
   S=bose/(im)*(chi-chi.Transpose().Conjugate());
   ch=0;
 // polarization factor
// neutrons only sense first 3x3 part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
 Matrix pol(1,md.nofcomponents,1,md.nofcomponents);
 Vector qijk(1,3);
  hkl2ijk(qijk,hkl, inputpars.cs.abc);
 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
 pol=0; double qsqr=qijk*qijk;
//    qijk(1)=hkl(1)/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)/inputpars.b;
//    qijk(3)=hkl(3)/inputpars.c;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qijk(i)*qijk(j)/qsqr;//(qijk*qijk);
    }}
    QQ=Norm(qijk);
// yes and for intermediate coupling we need another polarization factor
// because neutrons sense the first 6x6 part of S
    Matrix polICIC(1,md.nofcomponents,1,md.nofcomponents);
    Matrix polICn(1,md.nofcomponents,1,md.nofcomponents);
    Matrix polnIC(1,md.nofcomponents,1,md.nofcomponents);
    polICIC=0;polICn=0;polnIC=0;
    for(i=1;i<=6&&i<=md.nofcomponents;++i){
    for(j=1;j<=6&&j<=md.nofcomponents;++j){polICIC(i,j)=pol((i+1)/2,(j+1)/2);
                      if(i==1||i==3||i==5){polICIC(i,j)*=2.0;} // this accounts for the 
                      if(j==1||j==3||j==5){polICIC(i,j)*=2.0;} // fact that gs=2 and gl=1
    }}
    for(i=1;i<=3&&i<=md.nofcomponents;++i){
    for(j=1;j<=6&&j<=md.nofcomponents;++j){polnIC(i,j)=pol(i,(j+1)/2);
                      if(j==1||j==3||j==5){polnIC(i,j)*=2.0;} // fact that gs=2 and gl=1
    }}
    for(i=1;i<=6&&i<=md.nofcomponents;++i){
    for(j=1;j<=3&&j<=md.nofcomponents;++j){polICn(i,j)=pol((i+1)/2,j);
                      if(i==1||i==3||i==5){polICn(i,j)*=2.0;} // this accounts for the 
    }}

 //multiply polarization factor, formfactor and debeywallerfactor
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 for(l1=1;l1<=md.nofatoms;++l1){
       s=md.inM(i1,j1,k1,l1)*md.nofcomponents;
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  for(l2=1;l2<=md.nofatoms;++l2){
          ss=md.inM(i2,j2,k2,l2)*md.nofcomponents;
    for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
      ch(i,j)+=chi(s+i,ss+j);
      if((*inputpars.jjj[l1]).gJ==0&&(*inputpars.jjj[l2]).gJ==0)
      {S(s+i,ss+j)*=polICIC(i,j); 
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l1]).debyewallerfactor(QQ); //  debey waller factor
       if(i==2||i==4||i==6){S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l2]).debyewallerfactor(QQ); // debey waller factor
       if(j==2||j==4||j==6){S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
      }
      if((*inputpars.jjj[l1]).gJ==0&&(*inputpars.jjj[l2]).gJ!=0)
      {S(s+i,ss+j)*=polICn(i,j); 
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l1]).debyewallerfactor(QQ); //  debey waller factor
       if(i==2||i==4||i==6){S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
       S(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
      }
      if((*inputpars.jjj[l1]).gJ!=0&&(*inputpars.jjj[l2]).gJ==0)
      {S(s+i,ss+j)*=polnIC(i,j); 
       S(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debyewallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // debey waller factor
       if(j==2||j==4||j==6){S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
      }
      if((*inputpars.jjj[l1]).gJ!=0&&(*inputpars.jjj[l2]).gJ!=0)
      {S(s+i,ss+j)*=pol(i,j);
       S(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debyewallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
       S(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
      }
    }}   
  }
  }}}
 }
 }}}

 // determine dsigma in barns / cryst unit
 //divide by number of crytallographic unit cells  (ini.mf.n()) in magnetic unit cell
intensity=abs(Sum(S))/ini.mf.n()/PI/2.0*3.65/4.0/PI; 

// here should be entered factor  k/k' + absolute scale factor
if (ini.ki==0)
{if (ini.kf*ini.kf+0.4811*en<0)
   {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with kf=const=%g/A at (%g,%g,%g)\n",en,ini.kf,hkl(1),hkl(2),hkl(3));
    intensity=0;}
 else
 { ki=sqrt(ini.kf*ini.kf+0.4811*en);
   intensity*=ini.kf/ki;
 }
}
else
{if (ini.ki*ini.ki-0.4811*en<0)
   {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)\n",en,ini.ki,hkl(1),hkl(2),hkl(3));
    intensity=0;}
 else
 { 
  kf=sqrt(ini.ki*ini.ki-0.4811*en);
  intensity*=kf/ini.ki;
 }
}


#ifdef _THREADSREFINE
#undef ini
#undef inputpars
#undef J
#undef q
#undef hkl
#undef md
#undef ch
myinput->intensity=intensity;
MUTEX_LOCK(&mutex_loop);
thrdat.thread_id = thread_id;
EVENT_SIG(checkfinish);
MUTEX_UNLOCK(&mutex_loop);
#if defined  (__linux__) || defined (__APPLE__)
pthread_exit(NULL);
#else
return 0;
#endif
#else
return intensity;	
#endif
}
// *******************************************************************************************
// for rixs convert a vector from u123 coordinate system to ijk for a given azimuth
// eis,p and eos,p are polarisation vectors for sigma/pi plarisation in terms of
                                // the ijk coordinate system ijk form an euclidian righthanded 
                                //coordinate system j||b, k||(a x b) and i normal to j and k,
                                // where abc denote the crystal lattice vectors as defined in mcphas.j
/*RIXS intensity components:azimuth is defined as in Longfield et al. PRB 66 054417 (2002)
coordinate system u1,u2,u3: the scattering plane, defined by the
direction of the incident and final wave vectors k and k', contains u1 lying
perpendicular to Q and in the sense of k, and u3 antiparallel to the scattering
vector Q= k'-k.
angles for azimuth=0: alpha_i=angle(ai,u3), delta_i=angle(ai_perp,u1)
(where a1,a2,a3=a,b,c and ai_perp is the projection of ai onto the plane
perpendicular to Q. In the chosen experimental geometry
azimuth=0, when a1=a points to the x-ray source.
*/
void u123_to_ijk(ComplexVector & e,
              double azimuth,  // azimuth in radians
              Vector & qijk,   // Q in terms of coordinate system ijk
              Vector & hkl,   // miller indices of Q
              Vector & abc,   // lattice constants a b c and angles alpha beta gamma
              double &QQ     // absolute value of Q
             )
{// 1. write for azimuth=0 the vectors u1 u2 u3 in terms of ijk coordinate system
Vector u1(1,3),u2(1,3),u3(1,3),a(1,3),e1(1,3); // to store the vector components 
e1=0;e1(1)=1;
// u3 || Q
u3=qijk/QQ;
// u2 || Q x a, a from ddbdc2ijk
dadbdc2ijk(a,e1,abc); // transforms (100) given in terms of abc to ijk coordinate system
xproduct(u2,qijk,a); u2=u2/Norm(u2); // warning 
                                 //!! might be zero !! then choose e.g. b to point into source drection
if(Norm(u2)<SMALL_QUASIELASTIC_ENERGY){e1(1)=0;e1(2)=1; dadbdc2ijk(a,e1,abc);xproduct(u2,qijk,a);  }
// u1 = u2x u3
xproduct(u1,u2,u3);
// 2. for azimuth 0 write now down the components of vector e in ijk coordinate system
Vector eijkr(1,3),eijki(1,3);
eijkr=real(e(1))*u1+real(e(2))*u2+real(e(3))*u3;
eijki=imag(e(1))*u1+imag(e(2))*u2+imag(e(3))*u3;

// 3. now rotate vector eijk around qijk by the angle azimuth (radians)
Vector x(1,3),p(1,3),y(1,3),pp(1,3);double PP;
pp=(eijkr*u3)*u3;p=eijkr-pp;PP=Norm(p);x=p/PP;xproduct(y,u3,x);
    eijkr=pp+PP*cos(azimuth)*x+PP*sin(azimuth)*y;
pp=(eijki*u3)*u3;p=eijki-pp;PP=Norm(p);x=p/PP;xproduct(y,u3,x);
    eijki=pp+PP*cos(azimuth)*x+PP*sin(azimuth)*y;
for(int i=1;i<=3;++i){e(i)=complex<double>(eijkr(i),eijki(i));}

}

// *******************************************************************************************
// rixs function to calculate the polarisation vectors ei eo for sigma pi polarisation in terms of ijk
void calc_eps(ComplexVector & eis,
              ComplexVector & eip,
              ComplexVector & eir,
              ComplexVector & eil,
              ComplexVector & eos,
              ComplexVector & eop,
              ComplexVector & eor,
              ComplexVector & eol,
              inimcdis & ini,
              double azimuth,  // azimuth in radians
              Vector & qijk,   // Q in terms of coordinate system ijk
              Vector & hkl,   // miller indices of Q
              Vector & abc,   // lattice constants a b c and angles alpha beta gamma
              double &QQ,     // absolute value of Q
              double &en) // energy transfer in meV
                                // eis,p and eos,p are polarisation vectors for sigma/pi plarisation in terms of
                                // the ijk coordinate system ijk form an euclidian righthanded 
                                //coordinate system j||b, k||(a x b) and i normal to j and k,
                                // where abc denote the crystal lattice vectors as defined in mcphas.j
{/*RIXS intensity components: defined as in Longfield et al. PRB 66 054417 (2002)
coordinate system u1,u2,u3: the scattering plane, defined by the
direction of the incident and final wave vectors k and k', contains u1 lying
perpendicular to Q and in the sense of k, and u3 parallel to the scattering
vector Q= k-k'.
*/
// STEP 1 write down eips eip eos eop in terms of coordinate system u1 u2 u3
eis(1)=0; // sigma is parallel to u2
eis(2)=1;
eis(3)=0;
eos=eis;    
double ki=0,kf=0;
// determine ki kf from energy
if (ini.ki==0)
{if (ini.kf+5.0679e-7*en<0)
 {fprintf(stderr,"warning mcdisp - calculation of rixs intensity: energy transfer %g meV cannot be reached with kf=const=%g/A at (%g,%g,%g)\n",en,ini.kf,hkl(1),hkl(2),hkl(3));
 }
 else
 {kf=ini.kf;ki=ini.kf+5.0679e-7*en; // this holds for photons
 }
}
else
{if (ini.ki-5.0679e-7*en<0)
 {fprintf(stderr,"warning mcdisp - calculation of rixs intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)\n",en,ini.ki,hkl(1),hkl(2),hkl(3));
 }
 else
 {ki=ini.ki;kf=ini.ki-5.0679e-7*en;
 }
}
double ca=(ki*ki+QQ*QQ-kf*kf)/(2*ki*QQ);// c = cos s = sin
double cb=(kf*kf+QQ*QQ-ki*ki)/(2*kf*QQ);
double sb=sqrt(1-cb*cb);
double sa=sqrt(1-ca*ca);

eip(1)=ca;
eip(2)=0;
eip(3)=-sa;
eop(1)=-cb;
eop(2)=0;
eop(3)=-sb;
                                
 complex<double> imaginary(0,1);
double sq2m1=1/sqrt(2);
eor=sq2m1*(eos+imaginary*eop);
eol=sq2m1*(eos-imaginary*eop);
eir=sq2m1*(eis+imaginary*eip);
eil=sq2m1*(eis-imaginary*eip);

// STEP 2 transform the vectors to coordinate system i j k
u123_to_ijk(eis,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eip,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eil,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eir,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eos,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eop,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eor,azimuth,qijk,hkl,abc,QQ);
u123_to_ijk(eol,azimuth,qijk,hkl,abc,QQ);
}


// *******************************************************************************************
// rixs function to calculate rxs intensity from rixs susceptibility
double calc_irix(ComplexVector &ei,ComplexVector & eo,ComplexMatrix & chi)
{double Intensity=0;
// chi(1-9,1-9) ... Rij components 1,2,3,...9 are equivalent
//     1     2     3    4       5   6      7     8    9
// to eiei',eiej',eiek',ejei',ejej',ejek',ekei',ekej',ekek' 
 ComplexVector ee(1,9),dummy(1,9);int k=1;
 for(int i=1;i<=3;++i)for(int j=1;j<=3;++j){ee(k)=ei(i)*conj(eo(j));++k;}
for(int i=1;i<=9;++i){dummy(i)=0;for(int j=1;j<=9;++j){dummy(i)+=chi(i,j)*ee(j);
}} //matvec mult
 //( programmed here by indices because only the indices 1...9 of chi contain valid data 
 Intensity=real(ee*dummy); // note a*b of complex vectors conjugates the second vector
                              // so we do not need to conjugate the vector here
 return Intensity;
}










//*************** OLD OLD OLD .... removed june 2013 !!! . intcalc() is not used any more !!!

//**************************************************************************/
#ifdef _THREADSREFINE
#if defined  (__linux__) || defined (__APPLE__)
void *intcalc(void *input)
#else
DWORD WINAPI intcalc(void *input)
#endif
#else
double intcalc(ComplexMatrix & ch,int dimA, double en,inimcdis & ini,par & inputpars,jq & J,Vector & q,Vector & hkl,mdcf & md,int do_verbose,double epsilon)
#endif
{int i,j,i1,j1,k1,l1,t1,i2,j2,k2,l2,t2,s,ss,bmax,bbmax,b;
 double intensity=1.2;
 double QQ,ki,kf;

#ifdef _THREADSREFINE
   intcalcapr_input *myinput; myinput = (intcalcapr_input *)input;
   int thread_id = myinput->thread_id;
   int dimA = myinput->dimA;
   double en = myinput->En; 
   #define ini (*thrdat.ini[thread_id])
   #define inputpars (*thrdat.inputpars[thread_id])
   #define J (*thrdat.J[thread_id])
   #define q thrdat.q
   #define hkl thrdat.hkl
   #define md (*thrdat.md[thread_id])
   #define ch (*thrdat.ch[thread_id])
   // int do_verbose = myinput->do_verbose;
   double epsilon = myinput->epsilon; 
#endif
 
 complex<double> z(en,epsilon);
 complex<double> eps(epsilon*3,0);
 // determine chi
   ComplexMatrix chi(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   ComplexMatrix Ac(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   ComplexMatrix Acinv(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   ComplexMatrix Bc(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   Ac=0;Bc=0;
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
   bmax=md.baseindex_max(i1,j1,k1);
   ComplexMatrix chi0c(1,md.nofcomponents*bmax,1,md.nofcomponents*bmax);
   ComplexMatrix dd(1,md.nofcomponents*bmax,1,md.nofcomponents*bmax);
   ComplexMatrix cc(1,md.nofcomponents*bmax,1,md.nofcomponents*bmax);
   cc=0; dd=0;
   s=(index_s(i1,j1,k1,1,1,md,ini)-1)*md.nofcomponents;

 // to be removed - put to jjjpar !!! 
   for(l1=1;l1<=md.nofatoms;++l1){
   for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      b=md.baseindex(i1,j1,k1,l1,t1);   
   for(i=1;i<=md.nofcomponents;++i)
   {
     if (md.delta(i1,j1,k1)(b)>SMALL_QUASIELASTIC_ENERGY)
     { //normal inelastic intensity
      cc(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=1.0/(md.delta(i1,j1,k1)(b)-z);
      dd(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.0;
     }
    else if (md.delta(i1,j1,k1)(b)<-SMALL_QUASIELASTIC_ENERGY)
     {cc(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.0;
      dd(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=1.0/(-md.delta(i1,j1,k1)(b)+z);
     }
    else
     { 
     //quasielastic intensity ...  artificially we introduce a splitting epsilon !!! compare Jensen 91 p 158
     // factor 0.5 because every transition is counted as half positive and half negative energy...
     cc(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.5*eps/(eps-z);
     dd(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=-0.5*eps/(eps+z);
     }
    }}}
    chi0c=md.M(i1,j1,k1)*cc+md.M(i1,j1,k1).Transpose()*dd; 
   // myPrintComplexMatrix(stdout,chi0c); 
  // end of code to be romved !!!

   // for i1 j1 k1 site in magnetic unit cell:
   // assemble chi0c from the chi0's of the different atoms in the crystallographic basis
   //chi0c=0;
   //for(l1=1;l1<=md.nofatoms;++l1){
   //chi0c= from series of pointers to matrices with single ion chi for specified energy ...
   // mind we should change all this code here to work with lower dimension: t1 is not needed and can be omitted 
   // and matrices should get smaller !!!
   // }
   // myPrintComplexMatrix(stdout,chi0c); 
     

//myPrintComplexMatrix(stdout,cc); 
    for(i=1;i<=md.nofcomponents*bmax;++i){
     Ac(s+i,s+i)=1; // set diagonal elements 1 (make Ac a unit matrix)
    for(j=1;j<=md.nofcomponents*bmax;++j){
     Bc(s+i,s+j)=chi0c(i,j);
     }}

  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
     ss=(index_s(i2,j2,k2,1,1,md,ini)-1)*md.nofcomponents;
     bbmax=md.baseindex_max(i2,j2,k2);
     ComplexMatrix cc1(1,md.nofcomponents*bmax,1,md.nofcomponents*bbmax);
  
     cc1=chi0c*J.mati(J.in(i1,j1,k1),J.in(i2,j2,k2));
    for(i=1;i<=md.nofcomponents*bmax;++i){for(j=1;j<=md.nofcomponents*bbmax;++j){
      Ac(s+i,ss+j)-=cc1(i,j);
    }}
   }}}   
 }}}


//myPrintComplexMatrix(stdout,Ac); 
 chi=Ac.Inverse()*Bc;

 // determine chi'' and S (bose factor)
   complex<double> im(0,1.0);
   ComplexMatrix S(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   complex<double> bose;
   bose=1.0/(1.0-exp(-z*(1.0/KB/ini.T)));
   S=bose/(im)*(chi-chi.Transpose().Conjugate());
   ch=0;
 // polarization factor
// neutrons only sense first 3x3 part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
 Matrix pol(1,md.nofcomponents,1,md.nofcomponents);
 Vector qijk(1,3);
  hkl2ijk(qijk,hkl, inputpars.cs.abc);
 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
 pol=0; double qsqr=qijk*qijk;
//    qijk(1)=hkl(1)/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)/inputpars.b;
//    qijk(3)=hkl(3)/inputpars.c;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qijk(i)*qijk(j)/qsqr;//(qijk*qijk);
    }}
    QQ=Norm(qijk);
// yes and for intermediate coupling we need another polarization factor
// because neutrons sense the first 6x6 part of S
    Matrix polICIC(1,md.nofcomponents,1,md.nofcomponents);
    Matrix polICn(1,md.nofcomponents,1,md.nofcomponents);
    Matrix polnIC(1,md.nofcomponents,1,md.nofcomponents);
    polICIC=0;polICn=0;polnIC=0;
    for(i=1;i<=6&&i<=md.nofcomponents;++i){
    for(j=1;j<=6&&j<=md.nofcomponents;++j){polICIC(i,j)=pol((i+1)/2,(j+1)/2);
                      if(i==1||i==3||i==5){polICIC(i,j)*=2.0;} // this accounts for the 
                      if(j==1||j==3||j==5){polICIC(i,j)*=2.0;} // fact that gs=2 and gl=1
    }}
    for(i=1;i<=3&&i<=md.nofcomponents;++i){
    for(j=1;j<=6&&j<=md.nofcomponents;++j){polnIC(i,j)=pol(i,(j+1)/2);
                      if(j==1||j==3||j==5){polnIC(i,j)*=2.0;} // fact that gs=2 and gl=1
    }}
    for(i=1;i<=6&&i<=md.nofcomponents;++i){
    for(j=1;j<=3&&j<=md.nofcomponents;++j){polICn(i,j)=pol((i+1)/2,j);
                      if(i==1||i==3||i==5){polICn(i,j)*=2.0;} // this accounts for the 
    }}

 //multiply polarization factor, formfactor and debeywallerfactor
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 for(l1=1;l1<=md.nofatoms;++l1){
   for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      s=(index_s(i1,j1,k1,l1,t1,md,ini)-1)*md.nofcomponents;
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  for(l2=1;l2<=md.nofatoms;++l2){
   for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      ss=(index_s(i2,j2,k2,l2,t2,md,ini)-1)*md.nofcomponents;
    for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
      ch(i,j)+=chi(s+i,ss+j);
      if((*inputpars.jjj[l1]).gJ==0&&(*inputpars.jjj[l2]).gJ==0)
      {S(s+i,ss+j)*=polICIC(i,j); 
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l1]).debyewallerfactor(QQ); //  debey waller factor
       if(i==2||i==4||i==6){S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l2]).debyewallerfactor(QQ); // debey waller factor
       if(j==2||j==4||j==6){S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
      }
      if((*inputpars.jjj[l1]).gJ==0&&(*inputpars.jjj[l2]).gJ!=0)
      {S(s+i,ss+j)*=polICn(i,j); 
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l1]).debyewallerfactor(QQ); //  debey waller factor
       if(i==2||i==4||i==6){S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l1]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
       S(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
      }
      if((*inputpars.jjj[l1]).gJ!=0&&(*inputpars.jjj[l2]).gJ==0)
      {S(s+i,ss+j)*=polnIC(i,j); 
       S(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debyewallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
       S(s+i,ss+j)*=0.5*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // debey waller factor
       if(j==2||j==4||j==6){S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(-QQ);}else{S(s+i,ss+j)*=(*inputpars.jjj[l2]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
      }
      if((*inputpars.jjj[l1]).gJ!=0&&(*inputpars.jjj[l2]).gJ!=0)
      {S(s+i,ss+j)*=pol(i,j);
       S(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debyewallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
       S(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
      }
    }}   
  }}
  }}}
 }}
 }}}

 // determine dsigma in barns / cryst unit
 //divide by number of crytallographic unit cells  (ini.mf.n()) in magnetic unit cell
intensity=abs(Sum(S))/ini.mf.n()/PI/2.0*3.65/4.0/PI; 

// here should be entered factor  k/k' + absolute scale factor
if (ini.ki==0)
{if (ini.kf*ini.kf+0.4811*en<0)
   {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with kf=const=%g/A at (%g,%g,%g)\n",en,ini.kf,hkl(1),hkl(2),hkl(3));
    intensity=0;}
 else
 { ki=sqrt(ini.kf*ini.kf+0.4811*en);
   intensity*=ini.kf/ki;
 }
}
else
{if (ini.ki*ini.ki-0.4811*en<0)
   {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)\n",en,ini.ki,hkl(1),hkl(2),hkl(3));
    intensity=0;}
 else
 { 
  kf=sqrt(ini.ki*ini.ki-0.4811*en);
  intensity*=kf/ini.ki;
 }
}


#ifdef _THREADSREFINE
#undef ini
#undef inputpars
#undef J
#undef q
#undef hkl
#undef md
#undef ch
myinput->intensity=intensity;
MUTEX_LOCK(&mutex_loop);
thrdat.thread_id = thread_id;
EVENT_SIG(checkfinish);
MUTEX_UNLOCK(&mutex_loop);
#if defined  (__linux__) || defined (__APPLE__)
pthread_exit(NULL);
#else
return 0;
#endif
#else
return intensity;	
#endif
}








