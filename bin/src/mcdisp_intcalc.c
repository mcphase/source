//calculate intensities for given energy




//***********************************************************************//
// initialize intensity calculation for going beyond dipole approximation
//***********************************************************************
// returns 1 on success and zero on failure
//***********************************************************************
int intcalc_beyond_ini(inimcdis & ini,par & inputpars,mdcf & md,int do_verbose,Vector & hkl)
{int i,j,k,l,m,n,jmin,i1,j1,tn; Vector qijk(1,3);
 Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                  abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
 hkl2ijk(qijk,hkl, abc);
 // transforms Miller indices (in terms of reciprocal lattice abc*)
 // to Q vector in ijk coordinate system
 
//    qijk(1)=hkl(1)*2*PI/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)*2*PI/inputpars.b;
//    qijk(3)=hkl(3)*2*PI/inputpars.c;

 float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
 if(do_verbose==1) printf("#calculating intensity beyond dipole approximation\n");
// determine unitary transformation Matrix V (q)  Gamma and N for going beyond dip interaction
  Vector Gamma(1,ini.nofcomponents);
  complex<double> imaginary(0,1);
  // transition matrix Nij
  ComplexMatrix Nijkl(1,ini.nofcomponents,1,ini.nofcomponents);
  // transformation matrix Vij
  ComplexMatrix Vijkl(1,ini.nofcomponents,1,ini.nofcomponents);
  FILE * fin; //
  Vector mf(1,ini.nofcomponents);
   
  int sort=0;int maxiter=1000000;



 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  md.V(i,j,k)=0;  md.N(i,j,k)=0;
  for(l=1;l<=inputpars.nofatoms;++l){
  fin = fopen_errchk ("./results/mcdisp.trs","rb");
  jmin=0;
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3]&&l==(int)nn[4])
    {tn=(int)nn[5];++jmin;  
    // calculate delta(single ion excitation energy), 
    // Malphabeta(transition matrix elements)

      // do calculation for atom s=(ijkl)
//      for(ll=1;ll<=ini.nofcomponents;++ll)
//       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s

//      fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
//      if(nn[6]<SMALL){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");}
      
        j1=(*inputpars.jjj[l]).transitionnumber; // try calculation for transition  j
        (*inputpars.jjj[l]).transitionnumber=-tn; // try calculation for transition  j
        if(do_verbose==1)(*inputpars.jjj[l]).transitionnumber=tn;
      int nnt;
      nnt=(*inputpars.jjj[l]).dncalc(qijk,ini.T,Nijkl,md.est(i,j,k,l));
//       myPrintComplexMatrix(stdout,Nijkl); 

      (*inputpars.jjj[l]).transitionnumber=j1; // put back transition number for 1st transition
      if(nnt==0)
      {if(do_verbose)printf("#warning mcdisp - function dncalc not implemented for single ion module, only doing dipolar intensity\n");
       fclose(fin);return 0;}
      else
      {
       j1=md.baseindex(i,j,k,l,jmin); 
      
//       if(fabs(fabs(d)-fabs(nn[6]))>SMALLEDIF)
//        {fprintf(stderr,"ERROR mcdisp: reading mcdisp.trs with transition energy delta %g meV differnt from internal calculation %g meV %g\n",nn[6],d);	 
//         exit(EXIT_FAILURE);}
//       md.delta(i,j,k)(j1)=nn[6]; // set delta
     // diagonalizeMs to get unitary transformation matrix Us
     myEigenSystemHermitean (Nijkl,Gamma,Vijkl,sort=1,maxiter); 
	// conjugate:note the eigensystemhermitean returns eigenvectors as column vectors, but
	// the components need to be complex conjugated 

         // treat correctly case for neutron energy loss
	 if (nn[6]>=0){Vijkl=Vijkl.Conjugate();}

     if (Gamma(ini.nofcomponents)>=0&&fabs(Gamma(ini.nofcomponents-1))<SMALL) 
                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
                           {if (nn[6]>SMALL)
			    {md.sqrt_Gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=sqrt(Gamma(ini.nofcomponents));// gamma(ini.nofcomponents)=sqr(gamma^s)
                            }
			    else if (nn[6]<-SMALL)
                            {md.sqrt_Gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=imaginary*sqrt(Gamma(ini.nofcomponents));// gamma(ini.nofcomponents)=sqr(gamma^s)
                            }
 			    else
			    { //quasielastic line needs gamma=SMALL .... because Mijkl and therefore gamma have been set to 
			      // wn/kT instead of wn-wn'=SMALL*wn/kT (in jjjpar.cpp -mdcalc routines)
			      //set fix delta but keep sign
			          if (nn[6]>0){//md.delta(i,j,k)(j1)=SMALL;
  			     md.sqrt_Gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=sqrt(SMALL*Gamma(ini.nofcomponents));
                                              }
				  else        {//md.delta(i,j,k)(j1)=-SMALL;
                             md.sqrt_Gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=imaginary*sqrt(SMALL*Gamma(ini.nofcomponents));
			                      }
			    }
			   }else 
                           {fprintf(stderr,"ERROR eigenvalue of single ion matrix <0: ev1=%g ev2=%g ev3=%g ... evn=%g\n",Gamma(1),Gamma(2),Gamma(3),Gamma(ini.nofcomponents));
                            exit(EXIT_FAILURE);}
        for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){
        md.V(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Vijkl(m,n);
        md.N(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Nijkl(m,n);
        }}    
       }
    }}}
    fclose(fin);

  }}}}
  return 1;
}

//**************************************************************************/
#ifdef _THREADS
#ifdef __linux__
void *intcalc_approx(void *input)
#else
DWORD WINAPI intcalc_approx(void *input)
#endif
#else
double intcalc_approx(ComplexMatrix & chi,ComplexMatrix & chibey,Matrix & pol,Matrix & polICIC,Matrix & polICn,Matrix & polnIC, double & intensitybey,mfcf & ev_real,mfcf & ev_imag,mfcf & eev_real,mfcf & eev_imag,ComplexMatrix & Ec,int dimA, const ComplexMatrix &Tau, int level,double en, const inimcdis & ini,const par & inputpars,Vector & hkl,const mdcf & md,int do_verbose,double & QQ)
#endif
{//calculates approximate intensity for energylevel i - according to chapter 8.2 mcphas manual

#ifdef _THREADS
   intcalcapr_input *myinput; myinput = (intcalcapr_input *)input;
   int thread_id = myinput->thread_id;
   double intensitybey = myinput->intensitybey;
   #define chi (*thrdat.chi[thread_id])
   #define chibey (*thrdat.chibey[thread_id])
   #define pol (*thrdat.pol[thread_id])
   #define polICIC (*thrdat.polICIC[thread_id])
   #define polICn (*thrdat.polICn[thread_id])
   #define polnIC (*thrdat.polnIC[thread_id])
   #define ev_real (*thrdat.ev_real[thread_id])
   #define ev_imag (*thrdat.ev_imag[thread_id])
   #define eev_real (*thrdat.eev_real[thread_id])
   #define eev_imag (*thrdat.eev_imag[thread_id])
   #define Ec (*thrdat.Ec[thread_id])
   int level =  myinput->level;//, dimA = myinput->dimA, do_verbose = myinput->do_verbose;
   #define Tau (*thrdat.Tau[thread_id])
   double en = myinput->En; 
   #define ini (*thrdat.ini[thread_id])
   #define inputpars (*thrdat.inputpars[thread_id])
   #define hkl thrdat.hkl
   #define md (*thrdat.md[thread_id])
   double QQ;
#endif

 int i,j,i1,j1,k1,l1,t1,i2,j2,k2,l2,t2,s,ss,b,bb;
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
  ev_real.clear();ev_imag.clear();
  eev_real.clear();eev_imag.clear();

// determine chi
  //  chi=0;chibey=0;
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){


//     stau=(ini.mf.nb()*ini.mf.nc()*(i1-1)+ini.mf.nc()*(j1-1)+k1-1)*md.nofatoms;
//     s=stau*md.nofcomponents;
     
 for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
//     sstau=(ini.mf.nb()*ini.mf.nc()*(i2-1)+ini.mf.nc()*(j2-1)+k2-1)*md.nofatoms;
//     ss=sstau*md.nofcomponents;
    for(l1=1;l1<=md.nofatoms;++l1){
    for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
    for(l2=1;l2<=md.nofatoms;++l2){
    for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      ss=index_s(i2,j2,k2,l2,t2,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      bb=md.baseindex(i2,j2,k2,l2,t2);
        
    
    for(j=1;j<=md.nofcomponents;++j){
     if((ss-1)*md.nofcomponents+j==1){for(i=1;i<=ini.extended_eigenvector_dimension;++i)
                                        {eev_real.mf(i1,j1,k1)(ini.extended_eigenvector_dimension*(l1-1)+i)+=real(Ec(s,i)*Tau(s,level))*sqrt(fabs(en));// add this transition
                                         eev_imag.mf(i1,j1,k1)(ini.extended_eigenvector_dimension*(l1-1)+i)+=imag(Ec(s,i)*Tau(s,level))*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                        }
                                     }
    for(i=1;i<=md.nofcomponents;++i){
     chileft=conj(md.sqrt_gamma(i1,j1,k1)(md.nofcomponents*b,md.nofcomponents*b))*md.U(i1,j1,k1)((b-1)*md.nofcomponents+i,(b-1)*md.nofcomponents+md.nofcomponents)*Tau(s,level);
if(intensitybey>0)chileftbey=conj(md.sqrt_Gamma(i1,j1,k1)(md.nofcomponents*b,md.nofcomponents*b))*md.V(i1,j1,k1)((b-1)*md.nofcomponents+i,(b-1)*md.nofcomponents+md.nofcomponents)*Tau(s,level);

     chi((s-1)*md.nofcomponents+i,(ss-1)*md.nofcomponents+j)=
     PI*chileft*en*conj(Tau(ss,level))*conj(md.U(i2,j2,k2)((bb-1)*md.nofcomponents+j,(bb-1)*md.nofcomponents+md.nofcomponents))*md.sqrt_gamma(i2,j2,k2)(md.nofcomponents*bb,md.nofcomponents*bb);
  // en inserted  MR 9.3.11

     // here we fill the eigenvector mf with the information from chi
     if((ss-1)*md.nofcomponents+j==1){ev_real.mf(i1,j1,k1)(md.nofcomponents*(l1-1)+i)+=real(chileft)*sqrt(fabs(en));// add this transition
                                      ev_imag.mf(i1,j1,k1)(md.nofcomponents*(l1-1)+i)+=imag(chileft)*sqrt(fabs(en));// *sqrt(fabs(en)) inserted 13.3.2011 MR
                                      }

if(intensitybey>0){  chibey((s-1)*md.nofcomponents+i,(ss-1)*md.nofcomponents+j)=
     PI*chileftbey*en*conj(Tau(ss,level))*conj(md.V(i2,j2,k2)((bb-1)*md.nofcomponents+j,(bb-1)*md.nofcomponents+md.nofcomponents))*md.sqrt_Gamma(i2,j2,k2)(md.nofcomponents*bb,md.nofcomponents*bb);}
  // en inserted  MR 9.3.11
    }}
   }}}}
  }}}
 }}}

  complex<double> im(0,1.0);
 //  chi'' to  S (bose factor) ... fluctuation dissipation theorem
//myPrintComplexMatrix(stdout,chi); 
//myPrintComplexMatrix(stdout,Tau); 
//   S=bose*2*chi;                          replaced by putting bose factor to final sumS
//if(intensitybey>0)  Sbey=bose*2*chibey;   MR 2.4.10
//S=chi; if(intensitybey>0)  Sbey=chibey; substituted chi for S to safe computation time

 // polarization factor
// neutrons only sense first 3x3 part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
    pol=0;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qijk(i)*qijk(j)/(qijk*qijk);
    }}
// yes and for intermediate coupling we need another polarization factor
// because neutrons sense the first 6x6 part of S
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



 //multiply polarization factor, formfactor and debyewallerfactor
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

      //--------------------------------------------------------------------------------------------------
      if((*inputpars.jjj[l1]).gJ==0&&(*inputpars.jjj[l2]).gJ==0)
      {chi(s+i,ss+j)*=polICIC(i,j);
       chi(s+i,ss+j)*=0.5*(*inputpars.jjj[l1]).debyewallerfactor(QQ); // multiply (2S+L) with factor 1/2 to be conformant
                                                                    // to gj/2F(Q)<J>=M/2F(Q) in case of gj>0 (see below),debye waller factor
if(intensitybey>0){chibey(s+i,ss+j)*=polICIC(i,j);
                   chibey(s+i,ss+j)*=(*inputpars.jjj[l1]).debyewallerfactor(QQ);} //  debey waller factor
       if(i==2||i==4||i==6){chi(s+i,ss+j)*=(*inputpars.jjj[l1]).F(-QQ);}else{chi(s+i,ss+j)*=(*inputpars.jjj[l1]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
       chi(s+i,ss+j)*=0.5*(*inputpars.jjj[l2]).debyewallerfactor(QQ); // multiply (2S+L) with factor 1/2 to be conformant
                                                                    // to gj/2F(Q)<J>=M/2F(Q) in case of gj>0 (see below),debye waller factor
if(intensitybey>0) chibey(s+i,ss+j)*=(*inputpars.jjj[l2]).debyewallerfactor(QQ); // debey waller factor
       if(j==2||j==4||j==6){chi(s+i,ss+j)*=(*inputpars.jjj[l2]).F(-QQ);}else{chi(s+i,ss+j)*=(*inputpars.jjj[l2]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
      }
      //--------------------------------------------------------------------------------------------------
      if((*inputpars.jjj[l1]).gJ==0&&(*inputpars.jjj[l2]).gJ!=0)
      {chi(s+i,ss+j)*=polICn(i,j);
       chi(s+i,ss+j)*=0.5*(*inputpars.jjj[l1]).debyewallerfactor(QQ); // multiply (2S+L) with factor 1/2 to be conformant
                                                                    // to gj/2F(Q)<J>=M/2F(Q) in case of gj>0 (see below),debye waller factor
if(intensitybey>0){       chibey(s+i,ss+j)*=polICn(i,j);
       chibey(s+i,ss+j)*=(*inputpars.jjj[l1]).debyewallerfactor(QQ); }//  debey waller factor
       if(i==2||i==4||i==6){chi(s+i,ss+j)*=(*inputpars.jjj[l1]).F(-QQ);}else{chi(s+i,ss+j)*=(*inputpars.jjj[l1]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
       chi(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
if(intensitybey>0)  chibey(s+i,ss+j)*=(*inputpars.jjj[l2]).debyewallerfactor(QQ); // and debey waller factor
      }
      //--------------------------------------------------------------------------------------------------
      if((*inputpars.jjj[l1]).gJ!=0&&(*inputpars.jjj[l2]).gJ==0)
      {chi(s+i,ss+j)*=polnIC(i,j);
       chi(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debyewallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
       chi(s+i,ss+j)*=0.5*(*inputpars.jjj[l2]).debyewallerfactor(QQ);// multiply (2S+L) with factor 1/2 to be conformant
                                                                    // to gj/2F(Q)<J>=M/2F(Q) in case of gj>0 (see below),debye waller factor
if(intensitybey>0){       chibey(s+i,ss+j)*=polnIC(i,j);
       chibey(s+i,ss+j)*=(*inputpars.jjj[l1]).debyewallerfactor(QQ); // and  + debey waller factor
       chibey(s+i,ss+j)*=(*inputpars.jjj[l2]).debyewallerfactor(QQ); }// debey waller factor
       if(j==2||j==4||j==6){chi(s+i,ss+j)*=(*inputpars.jjj[l2]).F(-QQ);}else{chi(s+i,ss+j)*=(*inputpars.jjj[l2]).F(QQ);}
                               // mind here we should use different formfactors for spin and orbital components !!!
                               // formfactor +QQ..spin formfactor (j0), -QQ .. orbital formfactor (j0+j2)
      }
      //--------------------------------------------------------------------------------------------------
      if((*inputpars.jjj[l1]).gJ!=0&&(*inputpars.jjj[l2]).gJ!=0)
      {chi(s+i,ss+j)*=pol(i,j);
       chi(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debyewallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
       chi(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debyewallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
if(intensitybey>0){       chibey(s+i,ss+j)*=pol(i,j);
       chibey(s+i,ss+j)*=(*inputpars.jjj[l1]).debyewallerfactor(QQ); // and + debey waller factor
       chibey(s+i,ss+j)*=(*inputpars.jjj[l2]).debyewallerfactor(QQ); }// and  + debey waller factor
      }
      //--------------------------------------------------------------------------------------------------
    }}   
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

sumS=Sum(chi)/PI/2.0*3.65/4.0/PI/(double)ini.mf.n();sumS*=2.0*bose;
intensity=fabs(real(sumS));
                      if (real(sumS)<-0.1){fprintf(stderr,"ERROR mcdisp: dipolar approx intensity %g negative,E=%g, bose=%g\n",real(sumS),en,bose);exit(1);}
                      if (fabs(imag(sumS))>0.1){fprintf(stderr,"ERROR mcdisp: dipolar approx intensity %g %+g iimaginary\n",real(sumS),imag(sumS));exit(1);}
if(intensitybey>0){sumS=Sum(chibey)/PI/2.0*3.65/4.0/PI/(double)ini.mf.n();sumS*=2.0*bose;
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
#undef polICIC
#undef polICn
#undef polnIC
#undef ev_real
#undef ev_imag
#undef eev_real
#undef eev_imag
#undef Ec
MUTEX_LOCK(&mutex_loop);
thrdat.thread_id = thread_id;
EVENT_SIG(checkfinish);
MUTEX_UNLOCK(&mutex_loop);
#ifdef __linux__
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
#ifdef __linux__
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
 pol=0;
//    qijk(1)=hkl(1)/inputpars.a; // only correct for ortholattices !!!!
//    qijk(2)=hkl(2)/inputpars.b;
//    qijk(3)=hkl(3)/inputpars.c;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qijk(i)*qijk(j)/(qijk*qijk);
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
#ifdef __linux__
pthread_exit(NULL);
#else
return 0;
#endif
#else
return intensity;	
#endif
}
