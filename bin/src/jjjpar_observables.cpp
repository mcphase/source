/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/************************* OBSERVABLES **************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


/****************************************************************************/
/****************************************************************************/
// 0. phonon displacement p in A
/****************************************************************************/
/****************************************************************************/
int jjjpar::pcalc (Vector &mom, double & T, Vector &  Hxc,Vector & Hext ,ComplexMatrix & parstorage)
{ switch (module_type)
  {case 1: 
   case 2:
   case 4: 
   case 3: 
   case 5: fprintf(stderr,"Warning: phonons in internal modules not implemented, continuing ... \n");
          return 0;break;
   default: if (p==NULL) {mom=0;return false;} 
            else{(*p)(&mom,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);return true;}
  }
}


int  jjjpar::dp1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & p1,ComplexMatrix & ests)
{float delta=maxE;p1(1)=complex <double> (ninit,pinit); 
 switch (module_type)
  {case 0: if(dp1==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: phonons  not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*dp1)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&p1,&delta,&ests);}
           break;
   case 1:
   case 2:
   case 3:
   case 4:
   case 5: 
   default:if(transitionnumber<0)fprintf(stderr,"Warning: phonons in internal modules not implemented, continuing ... \n");
          return 0;break;
   }
}
/****************************************************************************/
/****************************************************************************/
// 1. MAGNETIC MOMENT in units  muB
/****************************************************************************/
/****************************************************************************/
int jjjpar::mcalc (Vector &mom, double & T, Vector &  Hxc,Vector & Hext ,ComplexMatrix & parstorage)
{double lnZ,U;
 switch (module_type)
  {case 1: kramer(mom,T,Hxc,Hext,lnZ,U);mom*=gJ;return true;break;
   case 2:
   case 4: (*iops).Jcalc(mom,T,Hxc,Hext,parstorage);mom*=gJ;return true;break;
   case 3: brillouin(mom,T,Hxc,Hext,lnZ,U);mom*=gJ;return true;break;
   case 5: cluster_Icalc(mom,T,Hxc,Hext,lnZ,U);mom*=gJ;return true;break; 
                                       // currently only 3 operators in cluster
                                       // implemented which are Ma Mb Mc
                                       // of (coupled) kramers doublet sipfs
                                       // ... in future this has to be substituted by cluster _mcalc
                                       // and operators of magnetic moments have to be handled specially
   default: if (m==NULL) {mom=0;return false;} 
            else{(*m)(&mom,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);return true;}
  }
}

int  jjjpar::dm1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & m1,ComplexMatrix & ests)
{float delta=maxE;m1(1)=complex <double> (ninit,pinit);
 ComplexVector uu1(1,m1.Hi());int nnt,i;
 switch (module_type)
  {case 0: if(dm1==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: dm1 calc  is not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*dm1)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&m1,&delta,&ests);}
           break;
   case 1: nnt=kramerdm(transitionnumber,T,Hxc,Hext,m1,delta);m1*=gJ;return nnt;break;
   case 2:
   case 4: uu1(1)=m1(1);
           nnt=(*iops).dJ1calc(transitionnumber,T,Hxc,Hext,uu1,delta,ests);
           for (i=1;i<=m1.Hi();++i)m1(i)=gJ*uu1(i);return nnt;break;
   case 3: nnt=brillouindm(transitionnumber,T,Hxc,Hext,m1,delta);m1*=gJ;return nnt;break;
   case 5:if(transitionnumber<0)fprintf(stderr,"Problem: dm1 calc in internal module cluster not implemented, continuing ... \n");break;
  default:if(transitionnumber<0)fprintf(stderr,"Problem: dm1 calc in internal module ... not implemented, continuing ... \n");
          break;
   }
return 0;
}

int jjjpar::Lcalc (Vector &Lmom, double & T, Vector &  Hxc,Vector & Hext ,ComplexMatrix & parstorage)
{double lnZ,U;
 switch (module_type)
  {case 1: kramer(Lmom,T,Hxc,Hext,lnZ,U);Lmom*=(2.0-gJ);return true;break;
   case 2:
   case 4: (*iops).Jcalc(Lmom,T,Hxc,Hext,parstorage);Lmom*=(2.0-gJ);return true;break;
   case 3: brillouin(Lmom,T,Hxc,Hext,lnZ,U);Lmom*=(2.0-gJ);return true;break;
   case 5: cluster_Icalc(Lmom,T,Hxc,Hext,lnZ,U);Lmom*=(2.0-gJ);return true;break; 
                                       // currently only 3 operators in cluster
                                       // implemented which are Ma Mb Mc
                                       // of (coupled) kramers doublet sipfs
                                       // ... in future this has to be substituted by cluster _mcalc
                                       // and operators of magnetic moments have to be handled specially
   default: if (L==NULL) {Lmom=0;return false;} 
            else{(*L)(&Lmom,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);return true;}
  }
}

int  jjjpar::dL1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & L1,ComplexMatrix & ests)
{float delta=maxE;L1(1)=complex <double> (ninit,pinit);
  switch (module_type)
  {static int washere=0;
   case 0: if(dL1==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: dL1 calc  is not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*dL1)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&L1,&delta,&ests);}
           break;
   case 1:
   case 2:
   case 3:
   case 4:
   case 5: 
   default: if(transitionnumber<0&& washere==0){washere=1;fprintf(stderr,"Problem: dL1calc in internal modules not implemented, continuing ... \n");}
          return 0;break;
   }
}


int jjjpar::Scalc (Vector &Smom, double & T, Vector &  Hxc,Vector & Hext ,ComplexMatrix & parstorage)
{double lnZ,U;
 switch (module_type)
  {case 1: kramer(Smom,T,Hxc,Hext,lnZ,U);Smom*=(gJ-1.0);return true;break;
   case 2:
   case 4: (*iops).Jcalc(Smom,T,Hxc,Hext,parstorage);Smom*=(gJ-1.0);return true;break;
   case 3: brillouin(Smom,T,Hxc,Hext,lnZ,U);Smom*=(gJ-1.0);return true;break;
   case 5: cluster_Icalc(Smom,T,Hxc,Hext,lnZ,U);Smom*=(gJ-1.0);return true;break; 
                                       // currently only 3 operators in cluster
                                       // implemented which are Ma Mb Mc
                                       // of (coupled) kramers doublet sipfs
                                       // ... in future this has to be substituted by cluster _mcalc
                                       // and operators of magnetic moments have to be handled specially
   default: if (S==NULL) {Smom=0;return false;} 
            else{(*S)(&Smom,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);return true;}
  }
}

int  jjjpar::dS1calc (double & T,Vector &  Hxc,Vector & Hext, ComplexVector & S1,ComplexMatrix & ests)
{float delta=maxE;S1(1)=complex <double> (ninit,pinit);
 switch (module_type)
  {case 0: if(dS1==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: dS1 calc  is not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*dS1)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&S1,&delta,&ests);}
           break;
   case 1:
   case 2:
   case 3:
   case 4:
   case 5: 
   default:if(transitionnumber<0)fprintf(stderr,"Problem: dS1calc in internal modules not implemented, continuing ... \n");
          return 0;break;
   }
}

/****************************************************************************/
/****************************************************************************/
//2 . NEUTRON SCATTERING OPERATOR  --------------------------------------
/****************************************************************************/
/****************************************************************************/


/****************************************************************************/
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
// according to stored eigenstate matrix est
// input: Qvec ..... Q Vector components 123=xyz=cab
/****************************************************************************/
int jjjpar::MQ(ComplexVector & Mq, Vector & Qvec)
{double J0,J2,J4,J6;
 double Q,th,ph;
            Q = Norm(Qvec); //dspacing
//            d = 2.0 * PI / Q; s=0.5 / d;
      J0=j0(Q);
      J2=j2(Q);
      J4=j4(Q);
      J6=j6(Q);
            complex<double>dummy;
switch (module_type)
  {case 0:  if (mq==NULL) {mom=0;return false;} 
            else{getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // for external module we must provide th and ph with respect
                                                       // to abc coordinate system
            (*mq)(&Mq,&th,&ph,&J0,&J2,&J4,&J6,&est);
            if(Norm(Zc)==0){fprintf(stderr,"WARNING mcdiff: Z(K) coefficients not found or zero in file %s\n",sipffilename);}

            return true;}break;
   case 2:  getpolar(Qvec(3),Qvec(1),Qvec(2),Q,th,ph); // internal module cfield needs transformation because
                                                       // of its convention ijk||yzx where  j||b k||axb and i||jxk
            Mq=(*iops).MQ(th,ph,J0,J2,J4,J6,Zc,est);
             // cfield module provides Mq(123)=Mq(xyz)
             // we must transform this to mcdiff internal ijk||yzx coordinate system
            dummy=Mq(3);Mq(3)=Mq(1);Mq(1)=Mq(2);Mq(2)=dummy;
            if(Norm(Zc)==0){fprintf(stderr,"WARNING mcdiff: Z(K) coefficients not found or zero in file %s\n",sipffilename);}
            return true;break;
   case 4:  getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // for so1ion we must th and ph with respect to abc coordinate system
             if(Norm(Zc)==0){fprintf(stderr,"WARNING mcdiff: Z(K) coefficients not found or zero in file %s\n",sipffilename);}
            Mq=(*iops).MQ(th,ph,J0,J2,J4,J6,Zc,est);return true;break;
   default: return false; // all other internal modules do not currently provide mq
  }
}

/****************************************************************************/
// returns transition element  dMQ(1..3) in order to be able to go beyond
//
// dipolar approximation in mcdisp - it requires a call to eigenstates first
//
//on input
//    transitionnumber has to be set correctly to that one which is to be computed
//    sign(transitionnumber)... 1... without printout, -1 with extensive printout
//    est		matrix with eigenstates, eigenvalues [meV], population numbers
//    T                 temperature
//     Q                 components of Qvector in euclidian coordinates 123=abc
//  on output
//    int   	total number of transitions
//    dMQ	<-|M(Q)|+> sqrt(n- - n+),  n+,n- population numbers
//    with Q the scattering operator according to Lovesey 11.4, p 222, eq 6.87b
//     (note that  <M(Q)>=-2x<Q>_TH in units of mb)
//    .... occupation number of states (- to + transition chosen according to transitionnumber)
//
/****************************************************************************/
int jjjpar::dMQ1calc(Vector & Qvec,double & T, ComplexVector & dMQ,ComplexMatrix & ests)
{double delta=maxE;dMQ(1)=complex <double> (ninit,pinit);
 double J0,J2,J4,J6;
 double Q,th,ph;
 int i;     complex<double>dummy; // introduced 3.4.10 MR
            Q = Norm(Qvec); //dspacing
      //      d = 2.0 * PI / Q; s=0.5 / d;
      J0=j0(Q);
      J2=j2(Q);
      J4=j4(Q);
      J6=j6(Q);
	 // calculate th and ph (polar angles of Q with respect to xyz of CEF)
 switch (module_type)
  {static int washere=0;

   case 0:if (ddnn!=NULL){getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph);
                          return (*ddnn)(&transitionnumber,&th,&ph,&J0,&J2,&J4,&J6,&ests,&T,&dMQ,&delta);break;}
          else {return 0;}
   case 2:  getpolar(Qvec(3),Qvec(1),Qvec(2),Q,th,ph); // for internal module cfield xyz||cba and we have to give dMQ1 polar angles with respect to xyz
            i=(*iops).dMQ1(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,dMQ);
            // and we have to switch indices in matrix nat(1..3,1..3) to conform with xyz||cba changed MR 3.4.10
            dummy=dMQ(1);dMQ(1)=dMQ(2);dMQ(2)=dMQ(3);dMQ(3)=dummy; // changed MR 3.4.10
            return i;break;
   case 4:  getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // for internal module so1ion xyz||abc and we have to give dMQ1 polar angles with respect to xyz
            return (*iops).dMQ1(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,dMQ);break;
   default: if(washere==0){fprintf(stderr,"Warning in scattering operator function dMQcalc - for ion %s \ngoing beyond dipolar approximation is not implemented\n",sipffilename);
                           washere=1;}
            return 0;
  }

}



/************************************************************************************/
//  RETURN TOTAL FORMFACTOR,
//    however if gJ=0 and Q>0 return spin form factor FS(Q)=<j0(Q)>
//            if gJ=0 and Q<0 return angular  form factor FL(Q)=<j0(Q)>+<j2(Q)>
//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
/************************************************************************************/
   double jjjpar::F(double Q)
// {if(gJ==0&&Q>0){return j0(Q);} // in case of intermediate coupling return spin form factor
//  if(gJ==0&&Q<0){return j0(-Q)+j2(-Q);} // in case of intermediate coupling return angular form factor
// return (j0(Q) + j2(Q) * (2 / gJ - 1)); // formfactor F(Q) for rare earth
   // Rewrote to use saved value if Q same as previous call.
   {
     unsigned int iq=nsaved; //iqmax = (unsigned int)(Qsaved[MAXSAVEQ]==-1e16?nsaved:6), iq=nsaved;
     for(unsigned int ic=MAXSAVEQ; ic--;) {
//      if(fabs(Q-Qsaved[iq])<1e-6) return Fsaved[iq];
        if(Q==Qsaved[iq]) return Fsaved[iq];  // Starts search at last saved value, indexed by nsaved.
        if(!(iq--)) iq=MAXSAVEQ-1;
     }
     double Fval;
     if(gJ) {                                       // gJ!=0
        Fval = (j0(Q) + j2(Q) * (2 / gJ - 1)); }
     else {
        if(Q>0) Fval = j0(Q);
        else Fval = j0(-Q)+j2(-Q); }
     Qsaved[nsaved] = Q; Fsaved[nsaved--] = Fval;   // nsaved starts at MAXSAVEQ
     if(!nsaved) nsaved=MAXSAVEQ-1;
     return Fval;
   }
   double jjjpar::j0(double Q)
  {double value=0,s; if(fabs(Q)<0.1)return 1.0;
   if(Np(1)!=0){value=jl(0,Q);}// here enter calculation from radial wave function parameters
   else
   {s=Q/4/PI;    value=magFFj0(1)*exp(-magFFj0(2)*s*s)+magFFj0(3)*exp(-magFFj0(4)*s*s)+magFFj0(5)*exp(-magFFj0(6)*s*s)+magFFj0(7);
   }return value;
  }
   double jjjpar::j1(double Q)
  {double value=0;   if(fabs(Q)<0.1)return 0.0;
   if(Np(1)!=0){value=jl(1,Q);}// here enter calculation from radial wave function parameters
   return value;
  }
   double jjjpar::j2(double Q)
  {double value=0,s;  if(fabs(Q)<0.1)return 0.0;
   s=Q/4/PI;
   if(Np(1)!=0){value=jl(2,Q);}// here enter calculation from radial wave function parameters
    else
   { value=magFFj2(1)*exp(-magFFj2(2)*s*s)+magFFj2(3)*exp(-magFFj2(4)*s*s)+magFFj2(5)*exp(-magFFj2(6)*s*s)+magFFj2(7);
    value*=s*s;
   }return value;
  }
   double jjjpar::j3(double Q)
  {double value=0; if(fabs(Q)<0.1)return 0.0;
   if(Np(1)!=0){value=jl(3,Q);}// here enter calculation from radial wave function parameters
   return value;
  }
   double jjjpar::j4(double Q)
  {double value=0,s;  if(fabs(Q)<0.1)return 0.0;
     s=Q/4/PI;
         if(Np(1)!=0){value=jl(4,Q);}// here enter calculation from radial wave function parameters
     else
   {  value=magFFj4(1)*exp(-magFFj4(2)*s*s)+magFFj4(3)*exp(-magFFj4(4)*s*s)+magFFj4(5)*exp(-magFFj4(6)*s*s)+magFFj4(7);
      value*=s*s;
   }return value;
  }
   double jjjpar::j5(double Q)
  {double value=0; if(fabs(Q)<0.1)return 0.0;
      if(Np(1)!=0){value=jl(5,Q);}// here enter calculation from radial wave function parameters
      return value;
  }
   double jjjpar::j6(double Q)
  {double value=0,s;   if(fabs(Q)<0.1)return 0.0;
     s=Q/4/PI;
         if(Np(1)!=0){value=jl(6,Q);}// here enter calculation from radial wave function parameters
     else
   {  value=magFFj6(1)*exp(-magFFj6(2)*s*s)+magFFj6(3)*exp(-magFFj6(4)*s*s)+magFFj6(5)*exp(-magFFj6(6)*s*s)+magFFj6(7);
      value*=s*s;
    } return value;
  }

   double jjjpar::jl(int l,double QA){
    int p,q, pmax=0;
    double Q=QA*0.5292;// convert Q from 1/A into 1/a0
    Vector coeff(1,9);
    for(p=1;p<=9;++p){if(Np(p)!=0){pmax=p;
                                   coeff(p)=Cp(p)*pow(2.0*Xip(p)/Q,Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of <j%i(Q)> failed due to Xi%i<=0 - continuing with <j%i(Q)>=0\n",l,p,l);return 0;}
                     }            }
    if(pmax==0){fprintf (stderr,"Warning: calculation of <j%i(Q)> failed - continuing with <j%i(Q)>=0\n",l,l);return 0;}

    long double value=0;
    for(p=1;p<=pmax;++p){
    for(q=1;q<=pmax;++q){
   if((int)Np(p)+(int)Np(q)<l+2){jl_lmax=l-1;return 0.0;}
  // condition: Np(p)+Np(q)>l+1- otherwise jl(Q) not possible to calculate
  // with this wave function
    value+=coeff(p)*coeff(q)*tl(l,(int)Np(p)+(int)Np(q),(Xip(p)+Xip(q))/Q);
                          }}
     return (double)value;
   }



long double jjjpar::tl(int l,int N,long double x)
     {double value=0.0;
      switch (l)
       { case 0: value=sn(1,N,x);break;
         case 1: value=sn(2,N,x)-cn(1,N,x);break;
         case 2: value=3*sn(3,N,x)-3*cn(2,N,x)-sn(1,N,x);break;
         case 3: value=cn(1,N,x)-15*cn(3,N,x)-6*sn(2,N,x)+15*sn(4,N,x);break;
         case 4: value=10*cn(2,N,x)-105*cn(4,N,x)+sn(1,N,x)-45*sn(3,N,x)+105*sn(5,N,x);break;
         case 5: value=-cn(1,N,x)+105*cn(3,N,x)-945*cn(5,N,x)+15*sn(2,N,x)-420*sn(4,N,x)+945*sn(6,N,x);break;
         case 6: value=-21*cn(2,N,x)+1260*cn(4,N,x)-10395*cn(6,N,x)-sn(1,N,x)+210*sn(3,N,x)-4725*sn(5,N,x)+10395*sn(7,N,x);break;
        default: fprintf(stderr,"Error function jjjpar:tl - value l=%i not implemented\n",l);exit(1);
       }
     return value;
     }
/*
long double jjjpar::sn(int n,int N,long double x)
   {complex <double> c(x,-1.0);
    long double value;
    value=(double)factorial(N-n)*imag(pow(c,-N+n-1));
    return value;
   }
long double jjjpar::cn(int n,int N,long double x)
   {complex <double> c(x,-1.0);
    long double value;
    value=(double)factorial(N-n)*real(pow(c,-N+n-1));
    return value;
   }
*/
// ------------------------------------------------------------------------- //
// Rewrite ::sn() and ::cn() to avoid using complex numbers, to run faster
// ------------------------------------------------------------------------- //
/* Complex powers, from Mathematica: z=x+I y; Do[Print[ComplexExpand[z^ex]], {ex, 1, 8}] (type I as <ESC>ii<Esc>)
 * x+I y
 * x^2+2 I x y-y^2
 * x^3-3 x y^2+I (3 x^2 y-y^3)
 * x^4-6 x^2 y^2+y^4+I (4 x^3 y-4 x y^3)
 * x^5-10 x^3 y^2+5 x y^4+I (5 x^4 y-10 x^2 y^3+y^5)
 * x^6-15 x^4 y^2+15 x^2 y^4-y^6+I (6 x^5 y-20 x^3 y^3+6 x y^5)
 * x^7-21 x^5 y^2+35 x^3 y^4-7 x y^6+I (-7 x^6 y+35 x^4 y^3-21 x^2 y^5+y^7)
 * x^8-28 x^6 y^2+70 x^4 y^4-28 x^2 y^6+y^8+I (-8 x^7 y+56 x^5 y^3-56 x^3 y^5+8 x y^7)
 * when y=-1:                          For powers z^{-p}, take -ve of Im part and div by (1+x^2)^p.
 * x                           +I( -1 )
 * -1+x^2                      +I( -2 x )
 *  -3 x+x^3                   +I(  1-3 x^2 )
 *  1-6 x^2+x^4                +I(  4 x-4 x^3 )
 *  5 x-10 x^3+x^5             +I( -1+10 x^2-5 x^4 )
 *  -1+15 x^2-15 x^4+x^6       +I( -6 x+20 x^3-6 x^5 )
 * -7 x+35 x^3-21 x^5+x^7      +I( 1-21 x^2+35 x^4-7 x^6)
 * 1-28 x^2+70 x^4-28 x^6+x^8  +I( 8 x-56 x^3+56 x^5-8 x^7)
 */

long double jjjpar::sn(int n,int N,long double x)    // Need imaginary part
 {
    long double denom=1.; if((-N+n-1)<0) denom=pow(1+x*x,-(-N+n-1));
    switch(-N+n-1) {
      case  0: return 0.; break;
      case  1: return (long double)factorial(N-n) *  -1.; break;
      case -1: return (long double)factorial(N-n) * (1. / denom); break;
      case  2: return (long double)factorial(N-n) *  -2*x; break;
      case -2: return (long double)factorial(N-n) * ( 2*x / denom); break;
      case  3: return (long double)factorial(N-n) *   (1-3*x*x); break;
      case -3: return (long double)factorial(N-n) * (-(1-3*x*x) / denom); break;
      case  4: return (long double)factorial(N-n) *   4*x * (1-x*x); break;
      case -4: return (long double)factorial(N-n) * (-4*x * (1-x*x) / denom); break;
      case  5: return (long double)factorial(N-n) *   (-1 + x*x * (10 - 5*x*x)); break;
      case -5: return (long double)factorial(N-n) * (-(-1 + x*x * (10 - 5*x*x)) / denom); break;
      case  6: return (long double)factorial(N-n) *   x * (-6 + x*x * (20 - 6*x*x)); break;
      case -6: return (long double)factorial(N-n) * (-x * (-6 + x*x * (20 - 6*x*x)) / denom);  break;
      case  7: return (long double)factorial(N-n) *   ( 1 + x*x * (-21 + x*x * (35 - 7*x*x))); break;
      case -7: return (long double)factorial(N-n) * (-( 1 + x*x * (-21 + x*x * (35 - 7*x*x))) / denom); break;
      case  8: return (long double)factorial(N-n) *   x * (8 + x*x * (-56 + x*x * (56 - 8*x*x))); break;
      case -8: return (long double)factorial(N-n) * (-x * (8 + x*x * (-56 + x*x * (56 - 8*x*x))) / denom); break;
      default: //fprintf(stderr,"jjjpar::sn() Bad power %i\n",-N+n-1); exit(-1);
         complex <double> c(x,-1.0); return (long double)(factorial((double)(N-n))*imag(pow(c,-N+n-1.)));
    }
 }
long double jjjpar::cn(int n,int N,long double x)    // Need real part
 {
    long double denom=1.; if((-N+n-1)<0) denom=pow(1+x*x,-(-N+n-1));
    switch(-N+n-1) {
      case  0: return 1.; break;
      case  1: return (long double)factorial(N-n) *  x; break;
      case -1: return (long double)factorial(N-n) * (x / denom); break;
      case  2: return (long double)factorial(N-n) *  (-1+x*x); break;
      case -2: return (long double)factorial(N-n) * ((-1+x*x) / denom); break;
      case  3: return (long double)factorial(N-n) *  x * (-3+x*x); break;
      case -3: return (long double)factorial(N-n) * (x * (-3+x*x) / denom); break;
      case  4: return (long double)factorial(N-n) *  (1 + x*x * (-6+x*x)); break;
      case -4: return (long double)factorial(N-n) * ((1 + x*x * (-6+x*x)) / denom); break;
      case  5: return (long double)factorial(N-n) *  x * (5 + x*x *(-10+x*x)); break;
      case -5: return (long double)factorial(N-n) * (x * (5 + x*x *(-10+x*x)) / denom); break;
      case  6: return (long double)factorial(N-n) *  (-1 + x*x * (15 + x*x * (-15+x*x))); break;
      case -6: return (long double)factorial(N-n) * ((-1 + x*x * (15 + x*x * (-15+x*x))) / denom); break;
      case  7: return (long double)factorial(N-n) *  x * (-7 + x*x *(35 + x*x * (-21 + x*x))); break;
      case -7: return (long double)factorial(N-n) * (x * (-7 + x*x *(35 + x*x * (-21 + x*x))) / denom); break;
      case  8: return (long double)factorial(N-n) *  ( 1 + x*x * (-28 + x*x * (70 + x*x * (-28+x*x)))); break;
      case -8: return (long double)factorial(N-n) * (( 1 + x*x * (-28 + x*x * (70 + x*x * (-28+x*x)))) / denom); break;
      default: //fprintf(stderr,"jjjpar::cn() Bad power %i\n",-N+n-1); exit(-1);
         complex <double> c(x,-1.0); return (long double)(factorial((double)(N-n))*imag(pow(c,-N+n-1.)));
    }
 }

// formfactor information print to fout
void jjjpar::FFinfo(FILE * fout) // has to be consistent with mcdisp_intcalc settings FF_type
                                 // and mcdiff settings of FF_type
{ if(FF_type!=0){
  switch(FF_type)
   {case 1: fprintf(fout,"no contribution to magnetic neutron/xray intensity Imag");break;
    case -1: fprintf(fout,"beyond dip.appr. for magnetic n-intensity Imag, no contribution to Imag_dip for this ion, FF coefficients:");break;
    case 3: fprintf(fout,"contribution to magnetic n-intensity Imag calc. in dip approx:  <M(Q)>=<L>*FL(Q)+2*<S>*FS(Q) with FL(Q)=(j0+j2) and FS(Q)=j0, FF coefficients:");break;
    case -3: fprintf(fout,"beyond dip.appr. for magnetic n-intensity Imag, Imag_dip calc. with:  <M(Q)>=<L>*FL(Q)+2*<S>*FS(Q) with FL(Q)=(j0+j2) and FS(Q)=j0, FF coefficients:");break;
    case 2: if(gJ!=0){fprintf(fout,"contribution to magnetic n-intensity Imag calc. in dip approx: <M(Q)>=<M>*F(Q) with F(Q)=j0-(1-2/gJ)j2  FF coefficients:");}
            else {fprintf(fout,"contribution to magnetic n-intensity Imag calc. in dip approx: <M(Q)>=<M>*F(Q) with F(Q)=j0  FF coefficients:");}
            break;
    case -2: if(gJ!=0){fprintf(fout,"beyond dip.appr. for magnetic n-intensity Imag, Imag_dip calc. with: <M(Q)>=<M>*F(Q) with F(Q)=j0-(1-2/gJ)j2  FF coefficients:");}
            else {fprintf(fout,"beyond dip.appr. for magnetic n-intensity Imag, Imag_dip calc. with: <M(Q)>=<M>*F(Q) with F(Q)=j0  FF coefficients:");}
            break;
    case 4:fprintf(fout,"ion included in rixs intensity calculation");break;
    default: fprintf(stderr,"Error: inconsistency in formfactor calculation\n"); exit(1);
   }

  if(FF_type!=1&&FF_type!=4){  
  if(Np(1)!=0){
  fprintf(fout," - formfactor calc. from radial wave function parameters in %s: <jl(Q)> considered up to l=%i",sipffilename,jl_lmax);
  }
  else
  {
  for (int j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",magFFj0(j));}
  for (int j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",magFFj2(j));}
  }                }
    }
  fprintf(fout, "\n");
}


#define DIV4PI 0.07957747154594766788444  // Value of 1/4/PI to save some computation time avoiding division

/************************************************************************************/
//   debyewallerfactor = exp(-2 * DWF *s*s)      (sf ~ exp(-2 DWF sin^2(theta) / lambda^2)=EXP (-W),  (2*DWF=B=8 pi^2 <u^2>)
/************************************************************************************/
   double jjjpar::debyewallerfactor(double & Q)
   {
     if(DWF==0) return 1.;  // Quick exit
     double s; unsigned int iq=DBWnsaved;        // Starts search at last saved value, indexed by DBWnsaved.
     for(unsigned int ic=MAXSAVEQ; ic--;) {
        if(Q==DBWQsaved[iq]) return DBWsaved[iq]; if(!(iq--)) iq=MAXSAVEQ-1; }
//  s=Q/4/PI;
    s=Q*DIV4PI;
    double DBW=exp(-2*DWF*s*s);
     DBWQsaved[DBWnsaved] = Q; DBWsaved[DBWnsaved--] = DBW;   // nsaved starts at MAXSAVEQ
     if(!DBWnsaved) DBWnsaved=MAXSAVEQ;
     return DBW;
//  return exp(-2*DWF*s*s);
   }

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
// DENSITIES  ----------------------------------------------------------
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/************************************************************************************/
// evaluate radial wave function
/************************************************************************************/
double jjjpar::radial_wavefunction(double rr) // rr given in Angstroems, returns R(r) in units of 1/A^1.5
   {//printf("%g ",rr);

    double R=0;int p;double a0=0.5292;
    int ok=0;
    double r=rr/a0;// r is the distance in units of a0
    for(p=1;p<=9;++p){if(Np(p)!=0){ok=1;
                                   R+=exp(-Xip(p)*r)*pow(r,Np(p)-1)*Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"\n\nWarning: calculation of radial wave function R(r=%g) failed due to Xi%i<=0 - continuing with R(r=%g)=0\n\n\n",r,p,r);return 0;}
                     }            }
    // now we have R in units of 1/a0^1.5
    R/=sqrt(a0*a0*a0);
    // now we have R in units of 1/A^1.5
//printf("%g ",R);
    if(R==0){R=1e-30;}
    if (ok==1) return R;


//  we have to find the 4f wavefunction R4f(r) for each single ion and the Zlm, cfield has nothing: so we have
//     to take this from chrgplt.bas - a little problem: how do we get the correct R4f(r) ? for a first attempt
//     we could just take the same for all RE.

    static int washere=0;
    if(washere==0){washere=1;fprintf (stderr,"\n\n!! Warning !!: radial wave function parameters not found, will use 4f hydrogen radial wave function\n\n\n");}
double rs;
//k^2 = 11 / 10 * 11 / 9 * 11 / 8 * 11 / 7 * 11 / 6 * 11 / 5 * 11 / 4 * 11 / 3 * 11 / 2 * 11 / 1 * 11
rs = rr * exp(-rr);
R = 280.4 * rs * rs * rs * rs  * exp(-1.5 * rr);
//printf("R4f(%g)=%g\n ",rr,R);
    if(R==0){R=1e-30;}
return R;
   }
/************************************************************************************/
   //functions to calculate radial matrix elements <r^n> from radial wave function in units of a0=0.5292 A
/************************************************************************************/
   double jjjpar::rk_from_radial_wavefunction(int k)
   {int p,q, pmax=0;
    Vector coeff(1,9);

    for(p=1;p<=9;++p){if(Np(p)!=0){pmax=p;
                                   coeff(p)=Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of <r^%i> failed due to Xi%i<=0 - continuing with <r^%i>=0\n",k,p,k);return 0;}
                     }            }
    if(pmax==0){fprintf (stderr,"Warning: calculation of <r^%i> failed - continuing with <r^%i>=0\n",k,k);return 0;}
    double rk=0;
    for(p=1;p<=pmax;++p){
    for(q=1;q<=pmax;++q){// the following does not work because of large factorials and powers
                        // rk+=coeff(p)*coeff(q)*factorial((int)Np(p)+(int)Np(q)+k)/pow(Xip(p)+Xip(q),Np(p)+Np(q)+k+1);
                        // therefore we have to substitute it by  calculating the product "by hand"
                        double product=1.0;
                        int nnk=(int)Np(p)+(int)Np(q)+k;
                        int i;
                        for(i=1;i<=nnk;++i){product*=i/(Xip(p)+Xip(q));}
                       rk+=coeff(p)*coeff(q)*product/(Xip(p)+Xip(q));
                        }}
     return rk;
   }

   int jjjpar::r2_from_radial_wavefunction() {r2=rk_from_radial_wavefunction(2);if(module_type==2||module_type==4){(*iops).r2=r2;}return true;}
   int jjjpar::r4_from_radial_wavefunction() {r4=rk_from_radial_wavefunction(4);if(module_type==2||module_type==4){(*iops).r4=r4;}return true;}
   int jjjpar::r6_from_radial_wavefunction() {r6=rk_from_radial_wavefunction(6);if(module_type==2||module_type==4){(*iops).r6=r6;}return true;}

void jjjpar::save_radial_wavefunction(const char * filename)
   {double r=0.1;
    FILE * fout;
    if (radial_wavefunction(r)==0){fprintf(stderr,"Warning: save_radial_wavefunction not possible\n");return;}
    fout=fopen_errchk(filename,"w");
    fprintf(fout,"# radial wave function for %s\n",sipffilename);
    fprintf(fout,"# the radial wave function is expanded as \n");
    fprintf(fout,"# R(r)=sum_p C_p R_Np,XIp(r)\n");
    fprintf(fout,"# R_Np,XIp(r)=r^(Np-1).exp(-XIp * r).(2 * XIp)^(Np+0.5)/sqrt(2Np!)\n");
    fprintf(fout,"# radial wave function parameters Np XIp Cp values are\n");
    fprintf(fout,"# tabulated in clementi & roetti Atomic data and \n");
    fprintf(fout,"# nuclear data tables 14 (1974) 177-478 for the transition metals\n");
    fprintf(fout,"# for rare earth parameters can be found in Freeman and Watson PR 127 (1962) 2058\n");
    fprintf(fout,"# and O. Sovers, J. Phys. Chem. Solids Vol 28 (1966) 1073\n");
    fprintf(fout,"# the parameters used are (a0=0.5292 A): \n");
    int p;
    for(p=1;p<=9;++p){if(Np(p)!=0){fprintf(fout,"#! N%i=%g XI%i=%g /a0 C%i=%g\n",p,Np(p),p,Xip(p),p,Cp(p));}}
    fprintf(fout,"# r[A]  vs R(r)[1/A^1.5]\n");
    for(r=0.01;r<=10;r*=1.05){fprintf(fout,"%8.8g  %8.8g\n",r,radial_wavefunction(r));}
    fclose(fout);
   }

// sum over different Zlm using the coefficients a(l,m)
double jjjpar::zlmsum(Matrix & a, double & teta, double & fi)
{double ro,ct,ct2,st,st2,sfi,cfi,cf2,sf2;
 ct = cos(teta);                      // z/r
 ct2 = ct * ct;
 st = sin(teta);   // y/r=st sfi
 st2 = st * st;
 sfi = sin(fi);    // x/r=st cfi
 sf2=sfi*sfi;
 cfi = cos(fi);
 cf2=cfi*cfi;
 ro = a(0, 0) * cnst(0,0);

 ro = ro + a(1,-1) * cnst(1,-1) * st *sfi;
 ro = ro + a(1,0) * cnst(1,0) * ct;
 ro = ro + a(1,1) * cnst(1,1) * st *cfi;

 ro = ro + a(2, -2)* cnst(2,-2)  * 2 * st2 * sfi * cfi;
 ro = ro + a(2, -1)* cnst(2,-1)  * st * sfi * ct;
 ro = ro + a(2, 0)* cnst(2,0)  * (3 * ct2 - 1);
 ro = ro + a(2, 1)* cnst(2,1)  * st * cfi * ct;
 ro = ro + a(2, 2)* cnst(2,2)  * st2 * (cfi * cfi - sfi * sfi);

 ro = ro + a(3,-3) * cnst(3,-3) * (3* st2*cf2 - st2*sf2)*st*sfi ;
 ro = ro + a(3,-2) * cnst(3,-2) * 2*sfi*cfi*st2*ct;
 ro = ro + a(3,-1) * cnst(3,-1) * st*sfi*(5*ct2-1);
 ro = ro + a(3,0) * cnst(3,0) * ct *(5*ct2-3);
 ro = ro + a(3,1) * cnst(3,1) * st*cfi*(5*ct2-1);
 ro = ro + a(3,2) * cnst(3,2) * (st2*cf2-st2*sf2)*ct;
 ro = ro + a(3,3) * cnst(3,3) * st*cfi*(st2*cf2-3*st2*sf2);


 ro = ro + a(4, -4)* cnst(4,-4) * st2 * st2 * 4 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 ro = ro + a(4, -3)* cnst(4,-3) * ct * st * st2 * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 ro = ro + a(4, -2)* cnst(4,-2) * (7 * ct2 - 1) * 2 * st2 * cfi * sfi;
 ro = ro + a(4, -1)* cnst(4,-1) * st * sfi * ct * (7 * ct2 - 3);
 ro = ro + a(4, 0)* cnst(4,0) * (35 * ct2 * ct2 - 30 * ct2 + 3);
 ro = ro + a(4, 1) * cnst(4,1) * st * cfi * ct * (7 * ct2 - 3);
 ro = ro + a(4, 2)* cnst(4,2)  * (7 * ct2 - 1) * st2 * (cfi * cfi - sfi * sfi);
 ro = ro + a(4, 3)* cnst(4,3)  * ct * st * st2 * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 ro = ro + a(4, 4) * cnst(4,4) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);

// x/r=st cfi     y/r=st sfi    ct= z/r

 ro = ro + a(5,-5) * cnst(5,-5) *st*sfi*(5*st2*st2*cf2*cf2-10*st2*cf2*st2*sf2+st2*st2*sf2*sf2);
ro = ro + a(5,-4) * cnst(5,-4) *4*st*sfi*(st2*cf2-st2*sf2)*ct;
ro = ro + a(5,-3) * cnst(5,-3) *st*sfi*(3*st2*cf2-st2*sf2)*(9*ct2-1);
ro = ro + a(5,-2) * cnst(5,-2) *2*st*cfi*st*sfi*ct*(3*ct2-1);
ro = ro + a(5,-1) * cnst(5,-1) *st*sfi*(21*ct2*ct2-14*ct2+1);
ro = ro + a(5,0) * cnst(5,0) *ct*(63*ct2*ct2-70*ct2+15);
ro = ro + a(5,1) * cnst(5,1) *st*cfi*(21*ct2*ct2-14*ct2+1);
ro = ro + a(5,2) * cnst(5,2) *(st2*cf2-st2*sf2)*ct*(3*ct2-1);
ro = ro + a(5,3) * cnst(5,3) *st*cfi*(st2*cf2-3*st2*sf2)*(9*ct2-1);
ro = ro + a(5,4) * cnst(5,4) *(st2*st2*cf2*cf2-6*st2*cf2*st2*sf2+st2*st2*sf2*sf2)*ct;
ro = ro + a(5,5) * cnst(5,5) *st*cfi*(st2*st2*cf2*cf2-10*st2*cf2*st2*sf2+5*st2*st2*sf2*sf2);


 ro = ro + a(6, -6)* cnst(6,-6) * st2 * st2 * st2 * (6 * cfi * cfi * cfi * cfi * cfi * sfi - 20 * cfi * cfi * cfi * sfi * sfi * sfi + 6 * cfi * sfi * sfi * sfi * sfi * sfi);
 ro = ro + a(6, -5)* cnst(6,-5) * ct * st * st2 * st2 * (5 * cfi * cfi * cfi * cfi * sfi - 10 * cfi * cfi * sfi * sfi * sfi + sfi * sfi * sfi * sfi * sfi);
 ro = ro + a(6, -4)* cnst(6,-4) * (11 * ct2 - 1) * 4 * st2 * st2 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
 ro = ro + a(6, -3) * cnst(6,-3)* (11 * ct * ct2 - 3 * ct) * st2 * st * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
 ro = ro + a(6, -2)* cnst(6,-2) * 2 * st2 * sfi * cfi * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2);
 ro = ro + a(6, -1)* cnst(6,-1) * ct * st * sfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 ro = ro + a(6, 0)* cnst(6,0) * (231 * ct2 * ct2 * ct2 - 315 * ct2 * ct2 + 105 * ct2 - 5);
 ro = ro + a(6, 1)* cnst(6,1)  * ct * st * cfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
 ro = ro + a(6, 2)* cnst(6,2)  * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2) * st2 * (cfi * cfi - sfi * sfi);
 ro = ro + a(6, 3)* cnst(6,3)  * (11 * ct * ct2 - 3 * ct) * st2 * st * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
 ro = ro + a(6, 4)* cnst(6,4)  * (11 * ct2 - 1) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);
 ro = ro + a(6, 5)* cnst(6,5)  * ct * st * st2 * st2 * (cfi * cfi * cfi * cfi * cfi - 10 * cfi * cfi * cfi * sfi * sfi + 5 * cfi * sfi * sfi * sfi * sfi);
 ro = ro + a(6, 6)* cnst(6,6) * st2 * st2 * st2 * (cfi * cfi * cfi * cfi * cfi * cfi - 15 * cfi * cfi * cfi * cfi * sfi * sfi + 15 * cfi * cfi * sfi * sfi * sfi * sfi - sfi * sfi * sfi * sfi * sfi * sfi);
 return ro;
 }



/****************************************************************************/
/****************************************************************************/
// 3.CHARGE DENSITIES  ----------------------------------------------------------
/****************************************************************************/
/****************************************************************************/

//***********************************************************************
// function to calculate coefficients of expansion of chargedensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
//***********************************************************************
int jjjpar::chargedensity_coeff (Vector &mom, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage)
{mom=0;
 switch (module_type)
  {case 1: fprintf(stderr,"Problem: chargedensity  in module kramer is not possible, continuing ... \n");
           return false;break;
   case 2:
   case 4: (*iops).chargedensity_coeffcalc(mom,T,Hxc,Hext,parstorage);
           break;
   case 3: fprintf(stderr,"Problem: chargedensity  in module brillouin is not possible, continuing ... \n");
           return false;break;
   case 0: if(cd_m==NULL){fprintf(stderr,"Problem: chargedensity  is not possible in module %s, continuing ... \n",modulefilename);
           return false;} else {(*cd_m)(&mom,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);}
           break;
   case 5:fprintf(stderr,"Problem: chargedensity is not possible in module cluster, continuing ... \n");
           return false;break;
   default:fprintf(stderr,"Problem: chargedensity is not possible in module, continuing ... \n");
           return false;break;
  }

// Indices for chargedensity
//            0 not used
//          0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
printf("# coefficients for density calculation\n");
printf("#chargedensity is expanded in tesseral harmonics Zlm\n\
#   ro(r)= sum_lm a(l,m) R^2(r) Zlm(Omega)\n#\n ");
for(int i=1;i<=CHARGEDENS_EV_DIM;++i)printf("#! a(%i,%i) =%12.6f\n",k[i],q[i],myround(mom(i)));
printf("\n");
return true;
}

int jjjpar::dchargedensity_coeff1(double & T,Vector &  Hxc,Vector & Hext, ComplexVector & chargedensity_coeff1,ComplexMatrix & ests)
{float delta=maxE;chargedensity_coeff1(1)=complex <double> (ninit,pinit);
 switch (module_type)
  {case 1: if(transitionnumber<0)fprintf(stderr,"Problem: chargedensity  in module kramer is not possible, continuing ... \n");
           return 0;break;
   case 2:
   case 4: return(*iops).dchargedensity_coeff1calc(transitionnumber,T,Hxc,Hext,chargedensity_coeff1,delta,ests);
           break;
   case 3: if(transitionnumber<0)fprintf(stderr,"Problem: chargedensity  in module brillouin is not possible, continuing ... \n");
           return 0;break;
   case 0: if(cd_dm==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: chargedensity  is not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*cd_dm)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&chargedensity_coeff1,&delta,&ests);}
           break;
   case 5:if(transitionnumber<0)fprintf(stderr,"Problem: chargedensity is not possible in module cluster, continuing ... \n");
           return 0;break;
   default:fprintf(stderr,"Problem: chargedensity is not possible in module, continuing ... \n");
           return 0;break;
  }
}

//***********************************************************************
// sub for calculation of charge density given a radiu R and polar angles teta,
// fi and expansion coefficients alm
//***********************************************************************
double jjjpar::chargedensity_calc (double & teta,double & fi,double & R, Vector & moments)
{double ro,rr;

if((module_type==0)&&(ro_calc!=NULL)){(*ro_calc)(&ro,&teta,&fi,&R,&moments,&gJ,&ABC,&sipffilename);return ro;}
if (R>4.0||R<0){ro=0;}else{
 int l,m;
 Matrix a(0,6,-6,6); 
 for(l=1;l<=5;l+=2){for(m=-l;m<=l;++m){a(l,m)=0;}}
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for spindensity
//          0 not used
//          0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// R given in Angstroems, returns R(r) in units of 1/A^1.5
   rr=radial_wavefunction(R);
   rr=rr*rr;// then the chargedensity will be in units of 1/A^3
int i;
 if(moments.Hi()==28)
 {
 for (i=1;i<=28;++i){a(k[i],q[i])=moments(i);}
  ro=-rr*zlmsum(a,teta,fi); // minus, because electrons are negative
}
 else
 {fprintf(stderr,"Error jjjpar.chargedensity_calc: dimension of moments=%i must be 28\n",moments.Hi());
  exit(EXIT_FAILURE);
 }
 }
return ro;
}


/****************************************************************************/
/****************************************************************************/
// 4.SPIN DENSITIES  ----------------------------------------------------------
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
// function to calculate coefficients of expansion of spindensity in terms
// of Zlm R^2(r) at a given temperature T and  effective field H
/****************************************************************************/
int jjjpar::spindensity_coeff (Vector &mom,int xyz, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage)
{mom=0;
 switch (module_type)
  {case 1: fprintf(stderr,"Problem: spindensity  in module kramer is not possible, continuing ... \n");
           return false;break;
   case 2:
   case 4: fprintf(stderr,"Problem: spindensity  in module so1ion and cfield do not work, continuing ... \n");
           return false;break;
// comment on module so1ion/cfield:
//fprintf(stderr,"Problem: calcmagdensity>0 in %s, spindensity  in module so1ion and cfield do not work correctly yet, quitting... \n",sipffilename);
//     exit(EXIT_FAILURE);}  // here I quit because it is yet unclear if the formulas programmed in are correct
                           // I assume that there is proportionality between
                           //         sum_i(2si+li) Zlm(Omega_i) to
                           // and
                           //               1/2 (J Olm(J) + Olm(J) J)   (see ionpars.cpp, function cfield() )
                           // with coefficients being gJ*tetan(l)
                           // ... this is probably not correct: 1. Wigner eckhardt holds only for spherical tensor operators and not for products of such
                           //                                   2. even if it holds, then the coefficients have to be calculated for spin and orbital contributions and added
                           // todo: check Wigner Eckhardt theorem and if it holds calculate coefficients
                           //
   case 3: fprintf(stderr,"Problem: spindensity  in module brillouin is not possible, continuing ... \n");
           return false;break;
   case 0: if(sd_m==NULL){fprintf(stderr,"Problem: spindensity  is not possible in module %s, continuing ... \n",modulefilename);
           return false;} else {(*sd_m)(&mom,&xyz,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);}
           break;
   case 5:fprintf(stderr,"Problem: spindensity is not possible in module cluster, continuing ... \n");
           return false;break;
   default:fprintf(stderr,"Problem: spindensity is not possible in module, continuing ... \n");
           return false;break;
  }
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
for(int i=1;i<=SPINDENS_EV_DIM;++i){printf("#! aS%i(%i,%i) =%12.6f\n",xyz,k[i],q[i],myround(mom(i)));}
return true;
}

int jjjpar::dspindensity_coeff1(double & T,Vector &  Hxc,Vector & Hext, ComplexVector & spindensity_coeff1,ComplexMatrix & ests)
{float delta=maxE;spindensity_coeff1(1)=complex <double> (ninit,pinit);
 switch (module_type)
  {case 1: if(transitionnumber<0)fprintf(stderr,"Problem: spindensity  in module kramer is not possible, continuing ... \n");
           return 0;break;
   case 2:
   case 4: if(transitionnumber<0)fprintf(stderr,"Problem: spindensity  in module  so1ion/cfeld is not possible, continuing ... \n");
           return 0;break;
   case 3: if(transitionnumber<0)fprintf(stderr,"Problem: spindensity  in module brillouin is not possible, continuing ... \n");
           return 0;break;
   case 0: if(sd_dm==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: spindensity  is not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*sd_dm)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&spindensity_coeff1,&delta,&ests);}
           break;
   case 5:if(transitionnumber<0)fprintf(stderr,"Problem: spindensity is not possible in module cluster, continuing ... \n");
           return 0;break;
   default:fprintf(stderr,"Problem: spindensity is not possible in module, continuing ... \n");
           return 0;break;
  }
}

//***********************************************************************
// sub for calculation of spin density component given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
double jjjpar::spindensity_calc (double & teta,double & fi,double & R, Vector & moments)
{double ro,rr;

if (R>4.0||R<0){ro=0;}else{

 Matrix a(0,6,-6,6); 
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849 
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// R given in Angstroems, returns R(r) in units of 1/A^1.5
   rr=radial_wavefunction(R);
   rr=rr*rr;// then the spindensity will be in units of 1/A^3
int i;
 if(moments.Hi()==49)
 {
 for (i=1;i<=49;++i){a(k[i],q[i])=moments(i);}
   ro=rr*zlmsum(a,teta,fi);
 }
 else
 {fprintf(stderr,"Error jjjpar.spindensitycalc: dimension of moments=%i must be 49\n",moments.Hi());
  exit(EXIT_FAILURE);
 }
 }
return ro;
}

//***********************************************************************
// sub for calculation of spin density vector given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
Vector  jjjpar::spindensity_calc (double & teta,double & fi,double & R, Vector & momentsx,Vector & momentsy,Vector & momentsz)
{double rr;
 static Vector mm(1,3);
if (R>4.0||R<0){mm=0;}else{

 Matrix a(0,6,-6,6);
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
// R given in Angstroems, returns R(r) in units of 1/A^1.5
   rr=radial_wavefunction(R);
   rr=rr*rr;// then the spindensity will be in units of 1/A^3
int i;
  for (i=1;i<=49;++i){a(k[i],q[i])=momentsx(i);}
   mm(1)=rr*zlmsum(a,teta,fi);
  for (i=1;i<=49;++i){a(k[i],q[i])=momentsy(i);}
   mm(2)=rr*zlmsum(a,teta,fi);
  for (i=1;i<=49;++i){a(k[i],q[i])=momentsz(i);}
   mm(3)=rr*zlmsum(a,teta,fi);
 }
return mm;
}

//***********************************************************************
// subs for calculation gradient of spin  density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix jjjpar::gradspindensity_calc(double & teta,double & fi,double & R, Vector & momentx, Vector & momenty, Vector & momentz)
{static Matrix grad(1,3,1,3);
 if (R>3.9||R<0){grad=0;}else{
 Vector m0(1,3);Vector m1(1,3);Vector m2(1,3);Vector m3(1,3);
 double d=0.01; // differential in Angstroem
 double teta1,teta2,teta3,fi1,fi2,fi3,R1,R2,R3;
 double ct,st,sf,cf;
 ct = cos(teta); st = sin(teta);   // y/r=st sfi
 sf = sin(fi);   cf = cos(fi);
 // now we use Jacobi Matrix (dr,dth,dfi)=J (dx,dy,dz)
  // dx                      dy             dz
 R1=R+d*st*cf;          R2=R+d*st*sf;       R3=R+d*ct;
 teta1=teta+d*ct*cf/R;  teta2=teta+d*ct*sf/R; teta3=teta-d*st/R;
 if(st>0){fi1=fi-d*sf/st/R;fi2=fi+d*cf/st/R;}else{fi1=fi;fi2=fi;} fi3=fi;

 m0=spindensity_calc(teta,fi,R,momentx,momenty,momentz);
 m1=spindensity_calc(teta1,fi1,R1,momentx,momenty,momentz);
 m2=spindensity_calc(teta2,fi2,R2,momentx,momenty,momentz);
 m3=spindensity_calc(teta3,fi3,R3,momentx,momenty,momentz);
 //            d/dx     d/dy      d/dz
 m1=(m1-m0)/d; m2=(m2-m0)/d; m3=(m3-m0)/d;grad=MatrixfromVectors(m1,m2,m3);
 }
 return grad;
}

/****************************************************************************/
/****************************************************************************/
// 5.ORBITAL MOMENT DENSITIES  ----------------------------------------------------------
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
// function to calculate coefficients of expansion of orbital moment density in terms
// of Zlm F(r) at a given temperature T and  effective field H
/****************************************************************************/
int jjjpar::orbmomdensity_coeff (Vector &mom,int xyz, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage)
{mom=0;
 switch (module_type)
  {case 1: fprintf(stderr,"Problem: orbmomdensity  in module kramer is not possible, continuing ... \n");
           return false;break;
   case 2:
   case 4: fprintf(stderr,"Problem: orbmomdensity  in module so1ion and cfield do not work, continuing ... \n");
           return false;break;
// comment on module so1ion/cfield:
//fprintf(stderr,"Problem: calcmagdensity>0 in %s, orbmomdensity  in module so1ion and cfield do not work correctly yet, quitting... \n",sipffilename);
//     exit(EXIT_FAILURE);}  // here I quit because it is yet unclear if the formulas programmed in are correct
                           // I assume that there is proportionality between
                           //         sum_i(2si+li) Zlm(Omega_i) to
                           // and
                           //               1/2 (J Olm(J) + Olm(J) J)   (see ionpars.cpp, function cfield() )
                           // with coefficients being gJ*tetan(l)
                           // ... this is probably not correct: 1. Wigner eckhardt holds only for spherical tensor operators and not for products of such
                           //                                   2. even if it holds, then the coefficients have to be calculated for spin and orbital contributions and added
                           // todo: check Wigner Eckhardt theorem and if it holds calculate coefficients
                           //
   case 3: fprintf(stderr,"Problem: orbmomdensity  in module brillouin is not possible, continuing ... \n");
           return false;break;
   case 0: if(od_m==NULL){fprintf(stderr,"Problem: orbmomdensity  is not possible in module %s, continuing ... \n",modulefilename);
           return false;} else {(*od_m)(&mom,&xyz,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&parstorage);}
           break;
   case 5:fprintf(stderr,"Problem: orbmomdensity is not possible in module cluster, continuing ... \n");
           return false;break;
   default:fprintf(stderr,"Problem: orbmomdensity is not possible in module, continuing ... \n");
           return false;break;
  }
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){printf("#! aL%i(%i,%i) =%12.6f\n",xyz,k[i],q[i],myround(mom(i)));}
return true;
}
int jjjpar::dorbmomdensity_coeff1(double & T,Vector &  Hxc,Vector & Hext, ComplexVector & orbmomdensity_coeff1,ComplexMatrix & ests)
{float delta=maxE;orbmomdensity_coeff1(1)=complex <double> (ninit,pinit);
 switch (module_type)
  {case 1: if(transitionnumber<0)fprintf(stderr,"Problem: orbmomdensity  in module kramer is not possible, continuing ... \n");
           return 0;break;
   case 2:
   case 4: if(transitionnumber<0)fprintf(stderr,"Problem: orbmomdensity  in module  so1ion/cfeld is not possible, continuing ... \n");
           return 0;break;
   case 3: if(transitionnumber<0)fprintf(stderr,"Problem: orbmomdensity  in module brillouin is not possible, continuing ... \n");
           return 0;break;
   case 0: if(od_dm==NULL){if(transitionnumber<0)fprintf(stderr,"Problem: orbmomdensity  is not possible in module %s, continuing ... \n",modulefilename);
           return 0;} else {return (*od_dm)(&transitionnumber,&T,&Hxc,&Hext,&gJ,&ABC,&sipffilename,&orbmomdensity_coeff1,&delta,&ests);}
           break;
   case 5:if(transitionnumber<0)fprintf(stderr,"Problem: orbmomdensity is not possible in module cluster, continuing ... \n");
           return 0;break;
   default:fprintf(stderr,"Problem: orbmomdensity is not possible in module, continuing ... \n");
           return 0;break;
  }
}

//***********************************************************************
// subs for calculation gradient of orbital moment density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix jjjpar::gradorbmomdensity_calc(double & teta,double & fi,double & R, Vector & momentx, Vector & momenty, Vector & momentz)
{static Matrix grad(1,3,1,3);
 if (R>3.9||R<0){grad=0;}else{
 Vector m0(1,3);Vector m1(1,3);Vector m2(1,3);Vector m3(1,3);
 double d=0.01; // differential in Angstroem
 double teta1,teta2,teta3,fi1,fi2,fi3,R1,R2,R3;
 double ct,st,sf,cf;
 ct = cos(teta);  st = sin(teta);   // y/r=st sfi
 sf = sin(fi);    cf = cos(fi);
 // now we use Jacobi Matrix (dr,dth,dfi)=J (dx,dy,dz)
  // dx                      dy             dz
 R1=R+d*st*cf;          R2=R+d*st*sf;       R3=R+d*ct;
 teta1=teta+d*ct*cf/R;  teta2=teta+d*ct*sf/R; teta3=teta-d*st/R;
 if(st>0){fi1=fi-d*sf/st/R;fi2=fi+d*cf/st/R;}else{fi1=fi;fi2=fi;} fi3=fi;

 m0=orbmomdensity_calc(teta,fi,R,momentx,momenty,momentz);
 m1=orbmomdensity_calc(teta1,fi1,R1,momentx,momenty,momentz);
 m2=orbmomdensity_calc(teta2,fi2,R2,momentx,momenty,momentz);
 m3=orbmomdensity_calc(teta3,fi3,R3,momentx,momenty,momentz);
 //            d/dx     d/dy      d/dz
 m1=(m1-m0)/d; m2=(m2-m0)/d; m3=(m3-m0)/d;grad=MatrixfromVectors(m1,m2,m3);
 }
 return grad;
}

//***********************************************************************
// sub for calculation of orbital moment density component given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
double jjjpar::orbmomdensity_calc (double & teta,double & fi,double & R, Vector & moments)
{double ro,rr;

if (R>4.0||R<0){ro=0;}else{

Matrix a(0,6,-6,6);
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

int i;
for (i=1;i<=49;++i){a(k[i],q[i])=moments(i);}

// R given in Angstroems, returns R(r) in units of 1/A^1.5
  rr=Fr(R);
  // then the orbital moment density will be in units of 1/A^3
  ro=rr*zlmsum(a,teta,fi);
 }
return ro;
}

//***********************************************************************
// sub for calculation of orbital moment density vector given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
Vector jjjpar::orbmomdensity_calc (double & teta,double & fi,double & R, Vector & momentsx, Vector & momentsy, Vector & momentsz)
{double rr;
 static Vector mm(1,3);
if (R>4.0||R<0){mm=0;}else{

Matrix a(0,6,-6,6);
//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

int i;
// R given in Angstroems, returns R(r) in units of 1/A^1.5
  rr=Fr(R);
  // then the orbital moment density will be in units of 1/A^3

 for (i=1;i<=49;++i){a(k[i],q[i])=momentsx(i);}
  mm(1)=rr*zlmsum(a,teta,fi);
 for (i=1;i<=49;++i){a(k[i],q[i])=momentsy(i);}
  mm(2)=rr*zlmsum(a,teta,fi);
 for (i=1;i<=49;++i){a(k[i],q[i])=momentsz(i);}
  mm(3)=rr*zlmsum(a,teta,fi);
 }
return mm;
}

/************************************************************************************/
// evaluate F(r) for orbital momentum density (see mynotes, balcar 1975)
/************************************************************************************/
double jjjpar::Fr(double rr) // evaluate F(r)=1/r integral_r^inf dx R^2(x)
                      // r in units of Angstroems, F(r) in units of 1/A^3
   {//printf("%g ",rr);

    double F_R=0;int p,pp,i;double a0=0.5292;
    double fp,fpp,sumi;
    int ok=0;
    double r=rr/a0;// r is the distance in units of a0
    for(p=1;p<=9;++p){if(Np(p)!=0){ok=1;
                                   if(Xip(p)<=0){fprintf (stderr,"\n\nWarning: calculation of radial integral F(r=%g) failed due to Xi%i<=0 - continuing with F(r=%g)=0\n\n\n",r,p,r);return 0;}
    fp=exp(-Xip(p)*r)*pow(r,Np(p)-1)*Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));

    for(pp=1;pp<=9;++pp){if(Np(pp)!=0){
     fpp=exp(-Xip(pp)*r)*pow(r,Np(pp)-1)*Cp(pp)*pow(2.0*Xip(pp),Np(pp)+0.5)/sqrt((double)factorial(2*(int)Np(pp)));
     sumi=0;for(i=1;i<=Np(p)+Np(pp)-2;++i){sumi+=pow(r,-i)/factorial((int)Np(p)+(int)Np(pp)-2-i)/pow(Xip(p)+Xip(pp),i+1);}
     sumi*=factorial(Np(p)+Np(pp)-2)/r;
                   F_R+=fp*fpp*sumi;
                          }             }
                     }            }
    // now we have F_R in units of 1/a0^3
    F_R/=a0*a0*a0;
    // now we have F_R in units of 1/A^3
    //printf("%g ",F_R);
    if (ok==1) return F_R;


//  we have to find the 4f wavefunction R4f(r) for each single ion and the Zlm, cfield has nothing: so we have
//     to take this from chrgplt.bas - a little problem: how do we get the correct R4f(r) ? for a first attempt
//     we could just take the same for all RE.

    static int washere=0;
    if(washere==0){washere=1;fprintf (stderr,"\n\n!! Warning !!: radial wave function parameters not found, will use 4f hydrogen radial wave function\n\n\n");}
double rs;
//k^2 = 11 / 10 * 11 / 9 * 11 / 8 * 11 / 7 * 11 / 6 * 11 / 5 * 11 / 4 * 11 / 3 * 11 / 2 * 11 / 1 * 11
rs = rr * exp(-rr);
fp = 280.4 * rs * rs * rs * rs  * exp(-1.5 * rr);

sumi=0;for(i=1;i<=8;++i){sumi+=pow(r,-i)/factorial(8-i)/pow(11.0,i+1);}
sumi*=factorial(8)/r;
F_R+=fp*fp*sumi;

    // now we have F_R in units of 1/a0^3
    F_R/=a0*a0*a0;

return F_R;
   }

/****************************************************************************/
/****************************************************************************/
// 4.CURRENT DENSITIES  ----------------------------------------------------------
/****************************************************************************/
/****************************************************************************/

//***********************************************************************
// sub for calculation of orbital current density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
Vector jjjpar::currdensity_calc (double & teta,double & fi,double & R, Vector & momentlx, Vector & momently, Vector & momentlz)
{double ro,rr;
 static Vector mm(1,3);
if (R>4.0||R<0){mm=0;}else{
 
 Matrix ax(0,6,-6,6),ay(0,6,-6,6),az(0,6,-6,6);
 Matrix bx(0,6,-6,6),by(0,6,-6,6),bz(0,6,-6,6);
 Matrix dx(0,6,-6,6),dy(0,6,-6,6),dz(0,6,-6,6);
 double ct,st,sf,cf,fp;
 ct = cos(teta);                      // z/r
 st = sin(teta);   // y/r=st sfi
 sf = sin(fi);    // x/r=st cfi
 cf = cos(fi);
 Vector Jr(1,3),Jth(1,3),Jfi(1,3); // Jacobi matrix
 Jr(1)=st*cf; Jr(2)=st*sf; Jr(3)=ct;
 Jth(1)=ct*cf;Jth(2)=ct*sf;Jth(3)=-st;
 if(st>0){Jfi(1)=-sf/st;Jfi(2)=cf/st;Jfi(3)=0;}else{Jfi=0;}
 Vector res(1,3),a(1,3);

//a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

int i;
  // R given in Angstroems, returns R(r) in units of 1/A^1.5
   ro=radial_wavefunction(R);
   ro=0.9274e-1*ro*ro/R;// then the spindensity will be in units of 1/A^3  ro=R(R);
   rr=0.9274e-1*Fr(R)/R;
  // then the orbital moment density will be in units of 1/A^3
  // then orbital current density will be in mb/A^4
  // transform to milliAmpere/A^2 by multiplying with 0.9274e-4

for (i=1;i<=49;++i){
ax(k[i],q[i])=momentlx(i);
ay(k[i],q[i])=momently(i);
az(k[i],q[i])=momentlz(i);
}

for (i=1;i<=49;++i){
a(1)=momentlx(i);
a(2)=momently(i);
a(3)=momentlz(i);
xproduct(res,a,Jr);
bx(k[i],q[i])=res(1);
by(k[i],q[i])=res(2);
bz(k[i],q[i])=res(3);
// now blm coefficients are calculated

// we need now the dlm
dx(k[i],q[i])=res(1);
dy(k[i],q[i])=res(2);
dz(k[i],q[i])=res(3);

     // add second term m* Jfi x al-m
if(q[i]!=0){
     a(1)=q[i]*ax(k[i],-q[i]);
     a(2)=q[i]*ay(k[i],-q[i]);
     a(3)=q[i]*az(k[i],-q[i]);
     xproduct(res,Jfi,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
     // add third term cottheta |m| Jth x alm
     if(st>0){
     a(1)=fabs(q[i])*ct/st*ax(k[i],q[i]);
     a(2)=fabs(q[i])*ct/st*ay(k[i],q[i]);
     a(3)=fabs(q[i])*ct/st*az(k[i],q[i]);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
             }
           }
// insert here addition of flm term !!
if(q[i]==-1){
     a(1)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*sf*ax(k[i],0);
     a(2)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*sf*ay(k[i],0);
     a(3)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*sf*az(k[i],0);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }
if(q[i]==+1){
     a(1)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*cf*ax(k[i],0);
     a(2)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*cf*ay(k[i],0);
     a(3)=-sqrt(k[i]*(double)(k[i]+1)/2.0)*cf*az(k[i],0);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }
if(q[i]>+1){
     fp=-sqrt(k[i]*(double)(k[i]+1)-q[i]*(q[i]-1));
     a(1)=fp*cf*ax(k[i],q[i]-1)-fp*sf*ax(k[i],-q[i]+1);
     a(2)=fp*cf*ay(k[i],q[i]-1)-fp*sf*ay(k[i],-q[i]+1);
     a(3)=fp*cf*az(k[i],q[i]-1)-fp*sf*az(k[i],-q[i]+1);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }
if(q[i]<-1){
     fp=-sqrt(k[i]*(double)(k[i]+1)-q[i]*(q[i]+1));
     a(1)=fp*cf*ax(k[i],q[i]+1)+fp*sf*ax(k[i],-q[i]-1);
     a(2)=fp*cf*ay(k[i],q[i]+1)+fp*sf*ay(k[i],-q[i]-1);
     a(3)=fp*cf*az(k[i],q[i]+1)+fp*sf*az(k[i],-q[i]-1);
     xproduct(res,Jth,a);
      dx(k[i],q[i])+=res(1);
      dy(k[i],q[i])+=res(2);
      dz(k[i],q[i])+=res(3);
            }

}

  mm(1)=ro*zlmsum(bx,teta,fi)+rr*zlmsum(dx,teta,fi);
  mm(2)=ro*zlmsum(by,teta,fi)+rr*zlmsum(dy,teta,fi);
  mm(3)=ro*zlmsum(bz,teta,fi)+rr*zlmsum(dz,teta,fi);

 }
return mm;
}

//***********************************************************************
// subs for calculation gradient of current density given a radiu R and polar angles teta,
// fi and expansion coeff. of Zlm R^2(r)
//***********************************************************************
 Matrix jjjpar::gradcurrdensity_calc(double & teta,double & fi,double & R, Vector & momentx, Vector & momenty, Vector & momentz)
{static Matrix grad(1,3,1,3);
 if (R>3.9||R<0){grad=0;}else{
 Vector m0(1,3);Vector m1(1,3);Vector m2(1,3);Vector m3(1,3);
 double d=0.01; // differential in Angstroem
 double teta1,teta2,teta3,fi1,fi2,fi3,R1,R2,R3;
 double ct,st,sf,cf;
 ct = cos(teta);  st = sin(teta);   // y/r=st sfi
 sf = sin(fi);    cf = cos(fi);
 // now we use Jacobi Matrix (dr,dth,dfi)=J (dx,dy,dz)
  // dx                      dy             dz
 R1=R+d*st*cf;          R2=R+d*st*sf;       R3=R+d*ct;
 teta1=teta+d*ct*cf/R;  teta2=teta+d*ct*sf/R; teta3=teta-d*st/R;
 if(st>0){fi1=fi-d*sf/st/R;fi2=fi+d*cf/st/R;}else{fi1=fi;fi2=fi;} fi3=fi;

 m0=currdensity_calc(teta,fi,R,momentx,momenty,momentz);
 m1=currdensity_calc(teta1,fi1,R1,momentx,momenty,momentz);
 m2=currdensity_calc(teta2,fi2,R2,momentx,momenty,momentz);
 m3=currdensity_calc(teta3,fi3,R3,momentx,momenty,momentz);
 //            d/dx     d/dy      d/dz
 m1=(m1-m0)/d; m2=(m2-m0)/d; m3=(m3-m0)/d;grad=MatrixfromVectors(m1,m2,m3);
 }
 return grad;
}


/****************************************************************************/
// returns transition element  drixs(1..9) in order to calculate the RIXS intensity
//
//  - it requires a call to eigenstates first
//
//on input
//    transitionnumber has to be set correctly to that one which is to be computed
//    sign(transitionnumber)... 1... without printout, -1 with extensive printout
//    est		matrix with eigenstates, eigenvalues [meV], population numbers
//    T                 temperature
//     Q                 components of Qvector in euclidian coordinates 123=abc
//  on output
//    int   	total number of transitions
//    drixs	<-|Rij|+> sqrt(n- - n+),  n+,n- population numbers
//               Rij transition operator according to Haverkort RIXS PRL
//    .... occupation number of states (- to + transition chosen according to transitionnumber)
//
/****************************************************************************/
int jjjpar::drixs1calc(Vector & Qvec,double & T, ComplexVector & drixs,ComplexMatrix & ests)
{double delta=maxE;drixs(1)=complex <double> (ninit,pinit);
 double J0,J2,J4,J6;
 double Q,th,ph;
        
            Q = Norm(Qvec); //dspacing
      //      d = 2.0 * PI / Q; s=0.5 / d;
      J0=j0(Q); // formfactor coefficients left here in case these are needed in future ??
      J2=j2(Q);
      J4=j4(Q);
      J6=j6(Q);
	 // calculate th and ph (polar angles of Q with respect to xyz of CEF)
 switch (module_type)
  {static int washere=0;

   case 0:if (rixs!=NULL){getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph);
                          return (*rixs)(&transitionnumber,&th,&ph,&J0,&J2,&J4,&J6,&ests,&T,&drixs,&delta);break;}
          else {return 0;}
   case 4:  getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // for internal module so1ion xyz||abc and we have to give dMQ1 polar angles with respect to xyz
            return (*iops).cfielddrixs1(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,drixs);break;
   case 2:  // for cfield because of coordinate rotation (complicated because of tensor) not implemented, not necessary I believe !
   default: drixs=0;if(washere==0){fprintf(stderr,"Warning in scattering operator function drixs1calc - for ion %s \ndoing RIXS  is not implemented\n",sipffilename);
                           washere=1;}
            return 0;
  }

}
