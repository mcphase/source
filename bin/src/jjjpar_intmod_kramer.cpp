//------------------------------------------------------------------------------------------------
//routine Icalc for kramers doublet
//------------------------------------------------------------------------------------------------
void jjjpar::kramer (Vector & Jret,double & T, Vector &  gjmbHxc,Vector & Hext, double & lnZ, double & U)
{ /*on input
    ABC(1...3)  A,M,Ci....saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/
  double alpha, betar, betai, lambdap,lambdap_KBT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, jap, jam, jbp, jbm, jcp, jcm,Z;
  double alpha_lambdap,alphaplambdap,alphaxlambdap;
  Vector gjmbH(1,3);
  gjmbH=gjmbHxc+gJ*MU_B*Hext;

  alpha = ABC[2] * gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
  betai = -ABC[3] * gjmbH[3];
  lambdap2 = alpha * alpha + betar * betar + betai * betai;
  lambdap = sqrt (lambdap2);
  lambdap_KBT=lambdap/KB/T;
  if (lambdap_KBT>700){lambdap_KBT=700;}
  if (lambdap_KBT<-700){lambdap_KBT=-700;}
  expm = exp (lambdap_KBT);
  expp = 1/expm; //=exp (-lambdap_KBT);
  Z = expp + expm;
  lnZ=log(Z);
  np = expp / Z;
  nm = expm / Z;
  U=lambdap*(np-nm); // energy

//  nennerp = (alpha - lambdap) * (alpha - lambdap) + betar * betar + betai * betai;
//  nennerm = (alpha + lambdap) * (alpha + lambdap) + betar * betar + betai * betai;
    alphaxlambdap=alpha*lambdap;
    alpha_lambdap=alpha-lambdap;
    alphaplambdap=alpha+lambdap;
    nennerp=  2.0*(-alphaxlambdap+lambdap2);    
    nennerm=  2.0*(alphaxlambdap+lambdap2);    

  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
      jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
      jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
      if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = ABC[2];
	}
      jcp = 0;
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
      jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
      jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
      if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = -ABC[2];
	}
      jcm = 0;
    }

  Jret[1] = np * jap + nm * jam;
  Jret[2] = np * jbp + nm * jbm;
  Jret[3] = np * jcp + nm * jcm;
//  printf ("Ha=%g Hb=%g Hc=%g Ja=%g Jb=%g Jc=%g \n", 
//     gjmbH[1]/MU_B/gjJ, gjmbH[2]/MU_B/gjJ, gjmbH[3]/MU_B/gjJ, J[1], J[2], J[3]);
}

int jjjpar::kramerdm(int & transitionnumber,double & T,Vector &  gjmbHxc,Vector & Hext,ComplexVector & u1,float & delta)
{ 
  /*on input
    transitionnumber ... number of transition to be computed - meaningless for kramers doublet, because there is only 1 transition
    ABC[i]	saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    u1(i)	<-|(Ji-<Ji>)|+> sqrt(tanh(delta/2kT))
*/
  double alpha, betar, betai, lambdap,lambdap_KBT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, nenner;
  complex<double> ja,jb,jc,i(0,1), jap, jam, jbp, jbm, jcp, jcm;
  double alpha_lambdap,alphaplambdap,alphaxlambdap;
  double Z;
  double lnz,u;
  Vector gjmbH(1,3);
  gjmbH=gjmbHxc+gJ*MU_B*Hext;

  static Vector Jret(1,3);
  Jret=0;
  // clalculate thermal expectation values (needed for quasielastic scattering)
  if(T>0){kramer(Jret,T,gjmbHxc,Hext,lnz,u);}else{T=-T;}
  int pr;
  pr=1;
  if (transitionnumber<0) {pr=0;transitionnumber*=-1;}

  alpha = ABC[2]* gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
  betai = -ABC[3] * gjmbH[3];
  lambdap2 = alpha * alpha + betar * betar + betai * betai;
  lambdap = sqrt (lambdap2);


  
  lambdap_KBT=lambdap/KB/T;
  if (lambdap_KBT>700){lambdap_KBT=700;}
  if (lambdap_KBT<-700){lambdap_KBT=-700;}
  expm = exp (lambdap_KBT);
  expp = 1/expm; //=exp (-lambdap_KBT);
  Z = expp + expm;
  np = expp / Z;
  nm = expm / Z;

    alphaxlambdap=alpha*lambdap;
    alpha_lambdap=alpha-lambdap;
    alphaplambdap=alpha+lambdap;
    nennerp=  2.0*(-alphaxlambdap+lambdap2);    
    nennerm=  2.0*(alphaxlambdap+lambdap2);    


if (transitionnumber==2)
{ delta=2*lambdap; //set delta !!!


    nenner=sqrt(nennerp*nennerm);

  if (nenner > SMALL)
    {
      ja = -ABC[1] * 2.0*(alpha * betar+i * betai * lambdap) / nenner;
      jb = -ABC[2] * 2.0 * (betar*betar+betai*betai) / nenner;
      jc = -ABC[3] * 2.0*(alpha*betai -i *betar*lambdap) / nenner;
    }
  else
    {
      if (alpha > SMALL)
	{ja = ABC[1];  // <-| is the ground state
  	 jb = 0;
         jc = -i*ABC[3];
	}
      else
	{ja = ABC[1];  // <+| is the ground state
  	 jb = 0;
         jc = i*ABC[3]; 	
	}
    }
 if (delta>SMALL)
 {// now lets calculate mat
 u1(1)=ja*sqrt(nm-np);
 u1(2)=jb*sqrt(nm-np);
 u1(3)=jc*sqrt(nm-np);
 }else
 {// quasielastic scattering needs epsilon * nm / KT ....
 u1(1)=ja*sqrt(nm/KB/T);
 u1(2)=jb*sqrt(nm/KB/T);
 u1(3)=jc*sqrt(nm/KB/T);
 }
}
else
{ delta=-SMALL; // transition within the same level
  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
      jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
      jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
      if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = -ABC[2];
	}
      jcp = 0;
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
      jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
      jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
      if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = ABC[2];
	}
      jcm = 0;
    }
 if (transitionnumber==1)
 {// now lets calculate mat
 u1(1)=(jam-Jret(1))*sqrt(nm/KB/T);
 u1(2)=(jbm-Jret(2))*sqrt(nm/KB/T);
 u1(3)=(jcm-Jret(3))*sqrt(nm/KB/T);
 }else{
 // now lets calculate mat
 u1(1)=(jap-Jret(1))*sqrt(np/KB/T);
 u1(2)=(jbp-Jret(2))*sqrt(np/KB/T);
 u1(3)=(jcp-Jret(3))*sqrt(np/KB/T);
 }
}
if (pr==1) printf ("delta=%4.6g meV\n",delta);

return 3; // kramers doublet has always exactly one transition + 2 levels (quasielastic scattering)!
}
/**************************************************************************/

Matrix jjjpar::krameropmat (int & n ,Vector &  gjmbHxc,Vector & Hext)
{
 /* on input
    ABC(1...3)  A,M,Ci....saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    n		which operator 0=Hamiltonian, 1,2,3=J1,J2,J3
    gjmbH	vector of effective field [meV]
  on output    
    operator matrix of Hamiltonian, J1, J2, J3 depending on n
*/
  // put matrix to format needed for library diagonalize function
//   for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=i1;++j1){
//    mat1(j1,i1)=imag(M(i1,j1));
//    mat1(i1,j1)=real(M(i1,j1));
//   }}
 // setup matrix and diagonalize
 //  Driver routine to compute the  eigenvalues and normalized eigenvectors
 //  of a complex Hermitian matrix z.The real parts of the elements must be
 //  stored in the lower triangle of z,the imaginary parts (of the elements
 //  corresponding to the lower triangle) in the positions
 //  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
 //  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
 //  True, otherwise  not ordered for sort = False. The real  and imaginary
 //  parts of the eigenvectors are  returned in  the columns of  zr and zi.
 //  The storage requirement is 3*n*n + 4*n complex numbers.
 //  All matrices and vectors have to be allocated and removed by the user.
 //  They are checked for conformance !
 // void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi,
 // 			   int sort, int maxiter)
 Vector gjmbH(1,3);
  gjmbH=gjmbHxc+gJ*MU_B*Hext;

Matrix opmat(1,2,1,2);
switch(n)
{case 0: opmat(1,1)= ABC[3]*gjmbH[3];      opmat(1,2)=ABC[2]*gjmbH[2];
         opmat(2,1)= -ABC[1]*gjmbH[1];     opmat(2,2)=-ABC[3]*gjmbH[3];break;
 case 1: opmat(1,1)= 0      ;opmat(1,2)=0;
         opmat(2,1)= ABC[1];opmat(2,2)=0;break;
 case 2: opmat(1,1)= 0      ;opmat(1,2)=-ABC[2];
         opmat(2,1)= 0      ;opmat(2,2)=0;break;
 case 3: opmat(1,1)= -ABC[3];opmat(1,2)=0;
         opmat(2,1)= 0     ;opmat(2,2)=ABC[3];break;
 default: fprintf(stderr,"ERROR operator calculation in module kramer - functio krameropmat: n=%i\n",n);exit(EXIT_FAILURE);
}
return opmat;
}
