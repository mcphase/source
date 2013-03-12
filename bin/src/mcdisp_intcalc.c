//calculate intensities for given energy




//***********************************************************************//
// initialize intensity calculation for going beyond dipole approximation
//***********************************************************************
// returns 1 on success and zero on failure
//***********************************************************************
void intcalc_ini(inimcdis & ini,par & inputpars,mdcf & md,int do_verbose,int do_gobeyond,int calc_rixs,Vector & hkl)
{int i,j,k,l,m,jmin,i1,j1,tn; Vector qijk(1,3);double QQ;
 Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                  abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
 hkl2ijk(qijk,hkl, abc);
    QQ=Norm(qijk);
 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
 
//    qijk(1)=hkl(1)*2*PI/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)*2*PI/inputpars.b;
//    qijk(3)=hkl(3)*2*PI/inputpars.c;

 float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
 if(do_verbose==1) printf("#initializing intensity calculation\n");
  double Gamman,gamma; 
  complex<double> imaginary(0,1);
  FILE * fin; 
  Vector mf(1,ini.nofcomponents);
      int mqdim=3;if(calc_rixs)mqdim=9;
      ComplexVector L1(1,3),S1(1,3),mq1_dip(1,mqdim),mq1(1,mqdim);

 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.nofatoms;++l){
  fin = fopen_errchk ("./results/mcdisp.trs","rb");
  jmin=0;
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3]&&l==(int)nn[4])
    {tn=(int)nn[5];++jmin;  
    // calculate delta(single ion excitation energy), 
    // Malphabeta(transition matrix elements)

//      fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
//      if(nn[6]<SMALL){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");}
      
        j1=(*inputpars.jjj[l]).transitionnumber; // try calculation for transition  j
       
      // do calculation for atom s=(ijkl)
      for(int ll=1;ll<=ini.nofcomponents;++ll){mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} 
                                               //mf ... exchange field vector of atom s
      (*inputpars.jjj[l]).transitionnumber=tn; // try calculation for transition  tn
      if(do_verbose==1)(*inputpars.jjj[l]).transitionnumber=-tn;

  if(calc_rixs)
      {if((*inputpars.jjj[l]).drixs1calc(qijk,ini.T,mq1,md.est(i,j,k,l))!=0)
       // calculate <-|Rijomega|+> see haverkort paper: transition operator for RIXS
      {if(do_verbose)printf("#calculating RIXS intensity for ion %s\n",(*inputpars.jjj[l]).sipffilename);
       (*inputpars.jjj[l]).FF_type*=4; // put FFTYPE to RIXS to indicate that this is implemented
      }
     else{ if(do_verbose)printf("#warning mcdisp - function drixs1 not implemented for single ion module of ion %s, no intensity from this ion\n",(*inputpars.jjj[l]).sipffilename);
           mq1=0;mq1(1)= complex <double> (1e-10,0.0);(*inputpars.jjj[l]).FF_type=1;mq1=mq1_dip;
          }
      }else
      {
      // try dipole approximation for this ion
      // dipole approx: <M(Q)>=<L>*FL(Q)+2*<S>*FS(Q) with FL(Q)=(j0+j2) and FS(Q)=j0
      if((*inputpars.jjj[l]).dL1calc(ini.T,mf,ini.Hext,L1,md.est(i,j,k,l))!=0&&
         (*inputpars.jjj[l]).dS1calc(ini.T,mf,ini.Hext,S1,md.est(i,j,k,l))!=0)
         {mq1=(*inputpars.jjj[l]).F(-QQ)*L1+2.0*(*inputpars.jjj[l]).F(QQ)*S1;(*inputpars.jjj[l]).FF_type=3;}
      else  {ComplexVector m1(1,3); // try dipole approximation using dm1calc
                                    // gJ=0 dipole approx: <M(Q)>=<M>*F(Q) with F(Q)=j0  
                                    // gJ>0 dipole approx: <M(Q)>=<M>*F(Q) with F(Q)=j0-(1-2/gJ)j2                
             if((*inputpars.jjj[l]).dm1calc(ini.T,mf,ini.Hext,m1,md.est(i,j,k,l))!=0)
             {mq1=m1*(*inputpars.jjj[l]).F(QQ);(*inputpars.jjj[l]).FF_type=2;}
             else {printf("#warning mcdisp - functions dmq1,dm1calc,dL1calcd,S1calc not implemented for single ion module of ion %s, no magnetic neutron intensity from this ion\n",(*inputpars.jjj[l]).sipffilename);
                   mq1=0;mq1(1)= complex <double> (1e-10,0.0);(*inputpars.jjj[l]).FF_type=1;
                  }
             }
      mq1_dip=mq1; 
      // try if going beyond dip approximation for formfactor works
     if(do_gobeyond){
      if((*inputpars.jjj[l]).dMQ1calc(qijk,ini.T,mq1,md.est(i,j,k,l))!=0)// calculate <-|M(Q)|+>
      {if(do_verbose)printf("#going beyond dipole approx for ion %s\n",(*inputpars.jjj[l]).sipffilename);
       (*inputpars.jjj[l]).FF_type*=-1; // put FFTYPE negative to indicate that going beyond works
      }
     else{ if(do_verbose)printf("#warning mcdisp - function dmq1 not implemented for single ion module of ion %s, only doing dipolar intensity\n",(*inputpars.jjj[l]).sipffilename);
           mq1=mq1_dip;
          }
                    }
       } // if calc_rixs
   
      Gamman=Norm2(mq1);mq1/=sqrt(Gamman);
      gamma=Norm2(mq1_dip);mq1_dip/=sqrt(gamma);
       (*inputpars.jjj[l]).transitionnumber=j1; // put back transition number for 1st transition
       j1=md.baseindex(i,j,k,l,jmin); 
      
         // treat correctly case for neutron energy loss
	 if (nn[6]<0){mq1=mq1.Conjugate();mq1_dip=mq1_dip.Conjugate();}
                            // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
        if (nn[6]>SMALL)
	    {md.sqrt_Gamma(i,j,k)(mqdim*j1)=sqrt(Gamman);md.sqrt_Gamma_dip(i,j,k)(mqdim*j1)=sqrt(gamma);}
	    else if (nn[6]<-SMALL)
            {md.sqrt_Gamma(i,j,k)(mqdim*j1)=imaginary*sqrt(Gamman);md.sqrt_Gamma_dip(i,j,k)(mqdim*j1)=imaginary*sqrt(gamma);}
	    else
	    { //quasielastic line needs gamma=SMALL .... because Mijkl and therefore gamma have been set to 
              // wn/kT instead of wn-wn'=SMALL*wn/kT (in jjjpar.cpp -mdcalc routines)
	      //set fix delta but keep sign
	      if (nn[6]>0){
  			     md.sqrt_Gamma(i,j,k)(mqdim*j1)=sqrt(SMALL*Gamman);
  			     md.sqrt_Gamma_dip(i,j,k)(mqdim*j1)=sqrt(SMALL*gamma);
                           }
	      else        {
                             md.sqrt_Gamma(i,j,k)(mqdim*j1)=imaginary*sqrt(SMALL*Gamman);
                             md.sqrt_Gamma_dip(i,j,k)(mqdim*j1)=imaginary*sqrt(SMALL*gamma);
	                     }
	    }

        for(m=1;m<=mqdim;++m){md.dMQs(i,j,k)(mqdim*(j1-1)+m)=mq1(m);}
        for(m=1;m<=mqdim;++m){md.dMQ_dips(i,j,k)(mqdim*(j1-1)+m)=mq1_dip(m);}
       
    }}}
    fclose(fin);
  }}}}
  
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
double intcalc_approx(ComplexMatrix & chi,ComplexMatrix & chibey,Matrix & pol, double & intensitybey,
        mfcf & qee_real,mfcf & qee_imag,ComplexMatrix & Echargedensity,
        mfcf & qsd_real,mfcf & qsd_imag,ComplexMatrix & Espindensity,
        mfcf & qod_real,mfcf & qod_imag,ComplexMatrix & Eorbmomdensity,
        mfcf & qep_real,mfcf & qep_imag,ComplexMatrix & Ephonon,
        mfcf & qem_real,mfcf & qem_imag,ComplexMatrix & Emagmom,
        mfcf & qes_real,mfcf & qes_imag,ComplexMatrix & Espin,
        mfcf & qel_real,mfcf & qel_imag,ComplexMatrix & Eorbmom,
        int dimA, const ComplexMatrix &Tau, int level,double en, const inimcdis & ini,const par & inputpars,Vector & hkl,/*const*/ mdcf & md,int do_verbose,int calc_rixs,double & QQ)
#endif
{//calculates approximate intensity for energylevel i - according to chapter 8.2 mcphas manual
int mqdim=3;if(calc_rixs)mqdim=9;

#ifdef _THREADS
   intcalcapr_input *myinput; myinput = (intcalcapr_input *)input;
   int thread_id = myinput->thread_id;
   double intensitybey = myinput->intensitybey;
   #define chi (*thrdat.chi[thread_id])
   #define chibey (*thrdat.chibey[thread_id])
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
   #define ini (*thrdat.ini[thread_id])
   #define inputpars (*thrdat.inputpars[thread_id])
   #define hkl thrdat.hkl
   #define md (*thrdat.md[thread_id])
   double QQ;
#endif

 int i,j,i1,j1,k1,l1,t1,i2,j2,k2,l2,t2,s,ss,s3,ss3,b,bb;
 double intensity=1.2; 
 double ki,kf;
 complex <double> sumS;
 complex <double> chileft;
 complex <double> chileftbey;
 Vector qijk(1,3);
 Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                  abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
 hkl2ijk(qijk,hkl, abc);
 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
//   qijk(1)=hkl(1)*2*PI/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)*2*PI/inputpars.b;
//    qijk(3)=hkl(3)*2*PI/inputpars.c;
    QQ=Norm(qijk);
 
 // init eigenvector to zero
  qee_real.clear();qee_imag.clear();
  qsd_real.clear();qsd_imag.clear();
  qod_real.clear();qod_imag.clear();
  qep_real.clear();qep_imag.clear();
  qem_real.clear();qem_imag.clear(); 
  qes_real.clear();qes_imag.clear();
  qel_real.clear();qel_imag.clear();


 // Added code to re-use previously calculated values of sqrt(gamma)*U and conj(U)*conj(sqrt(gamma)). mdl 110705
 int maxb=-1,bval,/*ncel=-1,*/nval; complex<double> defval(-0.1,0.), tval; md.ncel=-1;
 for(i2=1;i2<=ini.mf.na();++i2) for(j2=1;j2<=ini.mf.nb();++j2) for(k2=1;k2<=ini.mf.nc();++k2) { 
   bval=md.baseindex_max(i2,j2,k2); if(bval>maxb) maxb=bval; 
   nval=md.in(i2,j2,k2); if(nval>md.ncel) md.ncel=nval; } md.ncel++;
   if(maxb==md.nofcomponents) maxb++; // To ensure matrix is not square so overloaded operator ComplexMatrix=(complex<double>) sets all elements to 1e-16 not just diagonal.
   if(md.Ug==0) { md.Ug = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.Ug[i]=0; }
   if(md.gU==0) { md.gU = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.gU[i]=0; }
   if(md.bUg==0) { md.bUg = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.bUg[i]=0; }
   if(md.bgU==0) { md.bgU = new ComplexMatrix *[md.ncel+1]; for(i=1;i<=md.ncel;i++) md.bgU[i]=0; }
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){ int in1=md.in(i1,j1,k1);
   if(md.gU[in1]==0) { md.gU[in1] = new ComplexMatrix(1,mqdim,1,maxb); *md.gU[in1]=defval; }
   if(md.Ug[in1]==0) { md.Ug[in1] = new ComplexMatrix(1,mqdim,1,maxb); *md.Ug[in1]=defval; }
   if(md.bgU[in1]==0) { md.bgU[in1] = new ComplexMatrix(1,mqdim,1,maxb);*md.bgU[in1]=defval; }
   if(md.bUg[in1]==0) { md.bUg[in1] = new ComplexMatrix(1,mqdim,1,maxb);*md.bUg[in1]=defval; } }}}

// determine chi
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
    
 for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
    for(l1=1;l1<=md.nofatoms;++l1){
    for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
    for(l2=1;l2<=md.nofatoms;++l2){
    for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      ss=index_s(i2,j2,k2,l2,t2,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      bb=md.baseindex(i2,j2,k2,l2,t2);
      int in1=md.in(i1,j1,k1), in2=md.in(i2,j2,k2);

      if(intensitybey>0)
      {for(j=1;j<=mqdim;++j){for(i=1;i<=mqdim;++i){
        if((*md.bgU[in1])(i,b)==defval)  (*md.bgU[in1])(i,b)  = conj(md.sqrt_Gamma(i1,j1,k1)(mqdim*b))
                                                                 * md.dMQs(i1,j1,k1)((b-1)*mqdim+i);
        if((*md.bUg[in2])(j,bb)==defval) (*md.bUg[in2])(j,bb) = conj(md.dMQs(i2,j2,k2)((bb-1)*mqdim+j))
                                                                 * md.sqrt_Gamma(i2,j2,k2)(mqdim*bb);
                        
        //chileftbey=conj(md.sqrt_Gamma(i1,j1,k1)(mqdim*b))*md.dMQs(i1,j1,k1)((b-1)*mqdim+i)*Tau(s,level);
        //chibey((s-1)*mqdim+i,(ss-1)*mqdim+j)=
        //     PI*chileftbey*en*conj(Tau(ss,level))*conj(md.dMQs(i2,j2,k2)((bb-1)*mqdim+j))*md.sqrt_Gamma(i2,j2,k2)(mqdim*bb);
         chibey((s-1)*mqdim+i,(ss-1)*mqdim+j) = PI * (*md.bgU[in1])(i,b) * Tau(s,level) * en * conj(Tau(ss,level)) * (*md.bUg[in2])(j,bb);
        // en inserted  MR 9.3.11
      }}} // i,j,intensitybey

      if(!calc_rixs){
      for(j=1;j<=mqdim;++j){for(i=1;i<=mqdim;++i){
        if((*md.gU[in1])(i,b)==defval)  (*md.gU[in1])(i,b)  = conj(md.sqrt_Gamma_dip(i1,j1,k1)(mqdim*b))
                                                                 * md.dMQ_dips(i1,j1,k1)((b-1)*mqdim+i);
        if((*md.Ug[in2])(j,bb)==defval) (*md.Ug[in2])(j,bb) = conj(md.dMQ_dips(i2,j2,k2)((bb-1)*mqdim+j))
                                                                 * md.sqrt_Gamma_dip(i2,j2,k2)(mqdim*bb);
                        
        //chileft=conj(md.sqrt_gamma(i1,j1,k1)(mqdim*b))*md.dMQ_dips(i1,j1,k1)((b-1)*mqdim+i)*Tau(s,level);
        //chi((s-1)*mqdim+i,(ss-1)*mqdim+j)=
        //     PI*chileft*en*conj(Tau(ss,level))*conj(md.dMQ_dips(i2,j2,k2)((bb-1)*mqdim+j))*md.sqrt_gamma(i2,j2,k2)(mqdim*bb);
         chi((s-1)*mqdim+i,(ss-1)*mqdim+j) = PI * (*md.gU[in1])(i,b) * Tau(s,level) * en * conj(Tau(ss,level)) * (*md.Ug[in2])(j,bb);
        // en inserted  MR 9.3.11
      }}
                     }
    for(j=1;j<=md.nofcomponents;++j){ 
     if((ss-1)*md.nofcomponents+j==1){if(ini.calculate_chargedensity_oscillation)for(i=1;i<=CHARGEDENS_EV_DIM;++i)
                                        {qee_real.mf(i1,j1,k1)(CHARGEDENS_EV_DIM*(l1-1)+i)+=real(Echargedensity(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qee_imag.mf(i1,j1,k1)(CHARGEDENS_EV_DIM*(l1-1)+i)+=imag(Echargedensity(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_spindensity_oscillation)for(i=1;i<=SPINDENS_EV_DIM;++i)
                                        {qsd_real.mf(i1,j1,k1)(SPINDENS_EV_DIM*(l1-1)+i)+=real(Espindensity(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qsd_imag.mf(i1,j1,k1)(SPINDENS_EV_DIM*(l1-1)+i)+=imag(Espindensity(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_orbmomdensity_oscillation)for(i=1;i<=ORBMOMDENS_EV_DIM;++i)
                                        {qod_real.mf(i1,j1,k1)(ORBMOMDENS_EV_DIM*(l1-1)+i)+=real(Eorbmomdensity(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qod_imag.mf(i1,j1,k1)(ORBMOMDENS_EV_DIM*(l1-1)+i)+=imag(Eorbmomdensity(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_phonon_oscillation)for(i=1;i<=PHONON_EV_DIM;++i)
                                        {qep_real.mf(i1,j1,k1)(PHONON_EV_DIM*(l1-1)+i)+=real(Ephonon(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qep_imag.mf(i1,j1,k1)(PHONON_EV_DIM*(l1-1)+i)+=imag(Ephonon(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     if(ini.calculate_magmoment_oscillation)for(i=1;i<=MAGMOM_EV_DIM;++i)
                                        {qem_real.mf(i1,j1,k1)(MAGMOM_EV_DIM*(l1-1)+i)+=real(Emagmom(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         qem_imag.mf(i1,j1,k1)(MAGMOM_EV_DIM*(l1-1)+i)+=imag(Emagmom(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
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

  complex<double> im(0,1.0);

 // polarization factor
// neutrons only sense first mqdimxmqdim part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
    pol=0; double qsqr=qijk*qijk;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qijk(i)*qijk(j)/qsqr;//(qijk*qijk);
    }}
  
  // Precalculate values of Debye-Waller and Form Factors for this Q-vector to save calls to (*inputpars.jjj[ion]).* functions
  double DBWF[md.nofatoms+1];//, Fqm[md.nofatoms+1], Fqp[md.nofatoms+1];
  for(l1=1;l1<=md.nofatoms;++l1) {
     DBWF[l1] = (*inputpars.jjj[l1]).debyewallerfactor(QQ);
  }

 //multiply polarization factor, formfactor and debyewallerfactor
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 for(l1=1;l1<=md.nofatoms;++l1){
 for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
//   s=((((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*md.nofatoms+(l1-1));
      s=(index_s(i1,j1,k1,l1,t1,md,ini)-1);s3=s*mqdim;
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  for(l2=1;l2<=md.nofatoms;++l2){
  for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
//   ss=((((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*md.nofatoms+(l2-1));
      ss=(index_s(i2,j2,k2,l2,t2,md,ini)-1);ss3=ss*mqdim;

if(calc_rixs){
             for(i=1;i<=mqdim;++i){for(j=1;j<=mqdim;++j){
             chibey(s3+i,ss3+j) *= ( DBWF[l1] * DBWF[l2] ); 
                         }} // i,j,intesitybey            
      }else{
     if(intensitybey>0) {for(i=1;i<=mqdim;++i){for(j=1;j<=mqdim;++j){
           chibey(s3+i,ss3+j) *= ( pol(i,j) * DBWF[l1] * DBWF[l2] ); 
                         }}} // i,j,intesitybey
            for(i=1;i<=mqdim;++i){for(j=1;j<=mqdim;++j){
              chi(s3+i,ss3+j) *= ( pol(i,j) * DBWF[l1] * DBWF[l2] ); 
                         }}
            }
  }}
  }}}
 }}
 }}}


 // determine dsigma in barns per cryst unit cell !
 //divide by number of crystallographic unit cells  (ini.mf.n()) in magnetic unit cell


   double bose;
   if (fabs(en/KB/ini.T)>SMALL*0.1)
   {bose=1.0/(1.0-exp(-en*(1.0/KB/ini.T)));  
   }else{//quasielastic needs special treatment
         //if(fabs(en/KB/ini.T)>1e-30){bose=ini.T*KB/en;}
         //else{bose=1e30;} .... this is no good
         //(problem: quasielastic intensity depends on value of SMALL !!)
	 // in principle this SMALL in denominator has to cancel with epsilon
	 // in population matrix Mijkl(i.e. gamma) ... therefore we skip it:
	 // (for small energies delta_s the md.sqrt_gamma has been set = sqr(SMALL*gamma) and this is
	 // inserted into the calculation of chi above)
//         bose=ini.T*KB/(SMALL*0.1); // changed by MR 9.3.11 because factor omegar introduced in DMD formulas
//   bose=ini.T*KB;
     bose=ini.T*KB/en;  // replaced by the following line MR 9.3.11
   }
//  bose=fabs(bose);// this is to correctly consider omegar/|omegar| in the formula for chi'' ... introduced 2.4.10 MR
                    // deleted by MR 9.3.11 because new derivation of DMD gives factor omegar
                    //                      only (introduced above in bose) and
                    //                      we remove normalisation of tau in eigenvalueGeneral
                    //                      so that Tau is in units of sqrt of meV^-1

if (calc_rixs){sumS=Sum(chibey)/(double)ini.mf.n();sumS*=bose;
             intensity=fabs(real(sumS)); // do not know if this is right
             chibey=(bose/(double)ini.mf.n())*chibey;
             // determine chi(6x6 matrix in lagrange notation)  from chibeyond(9x9 matrix)
             chi(1,1)=chibey(1,1);// Rij components 1,2,3,...9 are equivalent to eiei',eiej',eiek',ejei',ejej',ejek',ekei',ekej',ekek' 
             chi(1,2)=chibey(2,2);//  Int components Iuv u=1,...,6=eiei,  ejej,  ekek,  eiej+ejei,    eiek+ekei,   ejek+ekej
             chi(1,3)=chibey(3,3);//                     v=1,...,6=ei'ei',ej'ej',ek'ek',ei'ej'+ej'ei',ei'ek'+ek'ei',ej'ek'+ek'ej'
             chi(1,4)=chibey(1,2)+chibey(2,1);
             chi(1,5)=chibey(1,3)+chibey(3,1);
             chi(1,6)=chibey(2,3)+chibey(3,2);
             chi(2,1)=chibey(4,4);
             chi(2,2)=chibey(5,5);
             chi(2,3)=chibey(6,6);
             chi(2,4)=chibey(4,5)+chibey(5,4);
             chi(2,5)=chibey(4,6)+chibey(6,4);
             chi(2,6)=chibey(5,6)+chibey(6,5);
             chi(3,1)=chibey(7,7);
             chi(3,2)=chibey(8,8);
             chi(3,3)=chibey(9,9);
             chi(3,4)=chibey(7,8)+chibey(8,7);
             chi(3,5)=chibey(7,9)+chibey(9,7);
             chi(3,6)=chibey(8,9)+chibey(9,8);
             chi(4,1)=chibey(1,4)+chibey(4,1);
             chi(4,2)=chibey(2,5)+chibey(5,2);
             chi(4,3)=chibey(3,6)+chibey(6,3);
             chi(4,4)=chibey(1,5)+chibey(5,1)+chibey(4,2)+chibey(2,4);
             chi(4,5)=chibey(1,6)+chibey(6,1)+chibey(3,4)+chibey(4,3);
             chi(4,6)=chibey(2,6)+chibey(6,2)+chibey(3,5)+chibey(5,3);
             chi(5,1)=chibey(1,7)+chibey(7,1);
             chi(5,2)=chibey(2,8)+chibey(8,2);
             chi(5,3)=chibey(3,9)+chibey(9,3);
             chi(5,4)=chibey(1,8)+chibey(8,1)+chibey(2,7)+chibey(7,2);
             chi(5,5)=chibey(1,9)+chibey(9,1)+chibey(3,7)+chibey(7,3);
             chi(5,6)=chibey(2,9)+chibey(9,2)+chibey(3,8)+chibey(8,3);
             chi(6,1)=chibey(4,7)+chibey(7,4);
             chi(6,2)=chibey(5,8)+chibey(8,5);
             chi(6,3)=chibey(6,9)+chibey(9,6);
             chi(6,4)=chibey(4,8)+chibey(8,4)+chibey(5,7)+chibey(7,5);
             chi(6,5)=chibey(4,9)+chibey(9,4)+chibey(6,7)+chibey(7,6);
             chi(6,6)=chibey(5,9)+chibey(9,5)+chibey(6,8)+chibey(8,6);
            }
       else{sumS=Sum(chi)/PI/2.0*3.65/4.0/PI/(double)ini.mf.n();sumS*=0.5*bose;
intensity=fabs(real(sumS));
                      if (real(sumS)<-0.1){fprintf(stderr,"ERROR mcdisp: dipolar approx intensity %g negative,E=%g, bose=%g\n",real(sumS),en,bose);exit(1);}
                      if (fabs(imag(sumS))>0.1){fprintf(stderr,"ERROR mcdisp: dipolar approx intensity %g %+g iimaginary\n",real(sumS),imag(sumS));exit(1);}
if(intensitybey>0){sumS=Sum(chibey)/PI/2.0*3.65/4.0/PI/(double)ini.mf.n();sumS*=0.5*bose;
                   intensitybey=fabs(real(sumS)); if (real(sumS)<-0.1){fprintf(stderr,"ERROR mcdisp: intensity in beyond dipolar approx formalism %g negative,E=%g, bose=%g\n\n",real(sumS),en,bose);exit(1);}
                                                   if (fabs(imag(sumS))>0.1){fprintf(stderr,"ERROR mcdisp: intensity  in beyond dipolar approx formalism %g %+g iimaginary\n",real(sumS),imag(sumS));exit(1);}
                  }
// here should be entered factor  k/k' + absolute scale factor
if (ini.ki==0)
{if (ini.kf*ini.kf+0.4811*en<0)
 {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with kf=const=%g/A at (%g,%g,%g)\n",en,ini.kf,hkl(1),hkl(2),hkl(3));
  intensity=0;
  intensitybey=0;
 }
 else
 { 
 ki=sqrt(ini.kf*ini.kf+0.4811*en);
 intensity*=ini.kf/ki;
 intensitybey*=ini.kf/ki;
 }
}
else
{if (ini.ki*ini.ki-0.4811*en<0)
 {fprintf(stderr,"warning mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)\n",en,ini.ki,hkl(1),hkl(2),hkl(3));
    intensity=0;
    intensitybey=0;
 }
 else
 {kf=sqrt(ini.ki*ini.ki-0.4811*en);
  intensity*=kf/ini.ki;
  intensitybey*=kf/ini.ki;
 }
}
} // if calc_rixs
#ifdef _THREADS
myinput->intensity=intensity;
myinput->intensitybey=intensitybey;
myinput->QQ=QQ;
#undef md
#undef hkl
#undef Tau
#undef ini
#undef inputpars
#undef chi
#undef chibey
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
#ifdef _THREADS
#if defined  (__linux__) || defined (__APPLE__)
void *intcalc(void *input)
#else
DWORD WINAPI intcalc(void *input)
#endif
#else
double intcalc(int dimA, double en,inimcdis & ini,par & inputpars,jq & J,Vector & q,Vector & hkl,mdcf & md,int do_verbose,double epsilon)
#endif
{int i,j,i1,j1,k1,l1,t1,i2,j2,k2,l2,t2,s,ss,bmax,bbmax,b;
 double intensity=1.2;
 double QQ,ki,kf;

#ifdef _THREADS
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
   int do_verbose = myinput->do_verbose;
   double epsilon = myinput->epsilon; 
#endif

 complex<double> z(en,epsilon);
 complex<double> eps(epsilon/4,0);
 // determine chi
   ComplexMatrix chi(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   ComplexMatrix Ac(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   ComplexMatrix Acinv(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   ComplexMatrix Bc(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
   Ac=0;Bc=0;
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
//   s=(ini.mf.nb()*ini.mf.nc()*(i1-1)+ini.mf.nc()*(j1-1)+k1-1)*md.nofcomponents*md.nofatoms;
   bmax=md.baseindex_max(i1,j1,k1);
   ComplexMatrix chi0c(1,md.nofcomponents*bmax,1,md.nofcomponents*bmax);
   ComplexMatrix dd(1,md.nofcomponents*bmax,1,md.nofcomponents*bmax);
   ComplexMatrix cc(1,md.nofcomponents*bmax,1,md.nofcomponents*bmax);
   cc=0; dd=0;
   s=(index_s(i1,j1,k1,1,1,md,ini)-1)*md.nofcomponents;

   for(l1=1;l1<=md.nofatoms;++l1){
   for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      b=md.baseindex(i1,j1,k1,l1,t1);   
   for(i=1;i<=md.nofcomponents;++i)
   {
     if (md.delta(i1,j1,k1)(b)>SMALL)
     { //normal inelastic intensity
      cc(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=1.0/(md.delta(i1,j1,k1)(b)-z);
      dd(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.0;
     }
    else if (md.delta(i1,j1,k1)(b)<-SMALL)
     {cc(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.0;
      dd(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=1.0/(-md.delta(i1,j1,k1)(b)+z);
     }
    else
     { 
     //quasielastic intensity ...  artificially we introduce a splitting epsilon !!! compare Jensen 91 p 158
     // factor 0.5 because every transition is counted as half positive and half negative energy...
     cc(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.5*eps/(eps-z);
     dd(md.nofcomponents*(b-1)+i,md.nofcomponents*(b-1)+i)=0.5*eps/(eps+z);
     }
    }}}

    chi0c=md.M(i1,j1,k1)*cc+md.M(i1,j1,k1).Transpose()*dd; 
//myPrintComplexMatrix(stdout,cc); 
    for(i=1;i<=md.nofcomponents*bmax;++i){
     Ac(s+i,s+i)=1; // set diagonal elements 1 (make Ac a unit matrix)
    for(j=1;j<=md.nofcomponents*bmax;++j){
     Bc(s+i,s+j)=chi0c(i,j);
     }}

  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
//   ss=(ini.mf.nb()*ini.mf.nc()*(i2-1)+ini.mf.nc()*(j2-1)+k2-1)*md.nofcomponents*md.nofatoms;
     ss=(index_s(i2,j2,k2,1,1,md,ini)-1)*md.nofcomponents;
     bbmax=md.baseindex_max(i2,j2,k2);
     ComplexMatrix cc1(1,md.nofcomponents*bmax,1,md.nofcomponents*bbmax);
  
  // for(l2=1;l2<=md.nofatoms;++l2)
  // for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2)

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
 //  bose=1.0;
   S=bose/(im)*(chi-chi.Transpose().Conjugate());
   
 // polarization factor
// neutrons only sense first 3x3 part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
 Matrix pol(1,md.nofcomponents,1,md.nofcomponents);
 Vector qijk(1,3);
 Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                  abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
 hkl2ijk(qijk,hkl, abc);
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
//   s=((((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*md.nofatoms+(l1-1))*md.nofcomponents;
      s=(index_s(i1,j1,k1,l1,t1,md,ini)-1)*md.nofcomponents;
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  for(l2=1;l2<=md.nofatoms;++l2){
   for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
//   ss=((((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*md.nofatoms+(l2-1))*md.nofcomponents;
      ss=(index_s(i2,j2,k2,l2,t2,md,ini)-1)*md.nofcomponents;
    for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
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


#ifdef _THREADS
#undef ini
#undef inputpars
#undef J
#undef q
#undef hkl
#undef md
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
