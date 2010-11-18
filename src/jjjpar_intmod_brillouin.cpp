//------------------------------------------------------------------------------------------------
//routine mcalc for brillouin function
//------------------------------------------------------------------------------------------------
void jjjpar::brillouin (Vector & Jret, double & T, Vector & gjmbH, double & lnZ, double & U)
{ /*on input
    ABC(1)  J=S....Spin quantum number
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/

// check dimensions of vector
if(Jret.Hi()!=3||gjmbH.Hi()!=3||ABC.Hi()!=1)
   {fprintf(stderr,"Error loadable module brillouin.so: wrong number of dimensions - check number of columns in file mcphas.j or number of parameters in single ion property file\n");
    exit(EXIT_FAILURE);}
    
double JJ,KBT,XJ,gmhkt,Jav,gmh,Z,X;

// program brillouin function for S=J=ABC(1)
JJ=ABC[1];
KBT=T*KB;
gmh=Norm(gjmbH);
gmhkt=gmh/KBT;

if(JJ*gmhkt>100||gmhkt>100){Jav=JJ;lnZ=JJ*gmhkt;}
else
{X=exp(gmhkt);
 XJ=exp(JJ*gmhkt);

// printf("1-X=%g gmhkt=%g\n",1-X,gmhkt);

 if (X<=1.000001){Z=2*JJ+1;Jav=0;}
 else
 {Z=(XJ*X-1/XJ)/(X-1.0);
  Jav=JJ*(XJ*X*X-1/XJ)+(JJ+1)*X*(1.0/XJ-XJ);
  Jav/=(X-1);
  Jav/=(XJ*X-1/XJ);
 }
//for (i=-JJ*2;i<=+0.000001;++i)
//{dd=i*gmhkt;
// if (dd<-700){expp=0;}else{expp=exp(dd);}
// Z += expp; //this is not yet Z, a factor exp(J gJ Heff/kT) is missing
//}


//Jav=0;
//for (i=-JJ*2;i<=+0.000001;++i)
//{dd=i*gmhkt;
// if (dd<-700){expp=0;}else{expp=exp(dd);}
// Jav+=(JJ+i)*expp/Z;
//}
//Z*=exp(JJ*gmhkt); //this is now the correct Z
lnZ=log(Z);
}

U=-gmh*Jav;



if (gmh>0)
{ Jret[1] = Jav*gjmbH(1)/gmh;
  Jret[2] = Jav*gjmbH(2)/gmh;
  Jret[3] = Jav*gjmbH(3)/gmh;
 }
 else
 {Jret=0;}
//  printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", H[1], H[2], H[3], m[1], m[2], m[3]);
}

/**************************************************************************/
// for mcdisp this routine is needed
/**************************************************************************/

int jjjpar::brillouindm(int & tn,double & T,Vector & gjmbH,ComplexMatrix & mat,float & delta)
{ 
  /*on input
    tn          transition-number
    ABC(1)      S=J spin quantum number
    g_J		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splittings [meV] 
    mat(i,j)	transition matrix elements ...
*/

static Vector Jret(1,3);
int pr;

  pr=1;
  if (tn<0) {pr=0;tn*=-1;}

  double JJ,KBT,XJ,gmhkt,gmh,Z,R,X,sinth,hxxyy,jjkt;
  complex <double> i(0,1),bx,by,bz;

// program brillouin function for S=J=ABC(1)
  JJ=ABC[1];
  KBT=T*KB;
  gmh=Norm(gjmbH);
  gmhkt=gmh/KBT;
  X=exp(gmhkt);
  XJ=exp(JJ*gmhkt);
// calculate Z and R
if (X==1.0){Z=2*JJ+1;R=0;}
else
{if(X>1e50){Z=XJ;R=-2.0*JJ*XJ;}
 else
 {Z=(XJ*X-1/XJ)/(X-1.0);
  R=JJ*(1/XJ-XJ*X*X)+(JJ+1)*X*(XJ-1.0/XJ);
  R/=0.5*(X-1)*(X-1);
 }
}
//printf("module brillouin -dmcalc: Z=%g R=%g X=%g XJ=%g |gjmbH|=%g\n",Z,R,X,XJ,gmh);

// calculate coefficients bx,by,bz
 hxxyy=gjmbH(1)*gjmbH(1)+gjmbH(2)*gjmbH(2);
 if (hxxyy/gjmbH(3)/gjmbH(3)>SMALL*SMALL)
 {sinth=sqrt(hxxyy)/gmh;
  bx=-gjmbH(2)+i*gjmbH(1)*gjmbH(3)/gmh;
  bx/=2*gmh*sinth;
  by=gjmbH(1)+i*gjmbH(2)*gjmbH(3)/gmh;
  by/=2*gmh*sinth;
  }
 else
 {sinth=0;by=0.5;
  if(gjmbH(3)>0)
  {bx=0.5*i;}
  else
  {bx=-0.5*i;}
 }
  bz=-i*sinth*0.5;
// -----------------------------------------
//printf("module brillouin -dmcalc: bx=%g by=%g bz=%g\n",bx,by,bz);
if (tn==2) // transition to finite energy
 {delta=gmh; //set delta !!!

 if (delta>SMALL)
  {// now lets calculate mat
  mat(1,1)=bx*conj(bx)*(-R/Z);
  mat(1,2)=bx*conj(by)*(-R/Z);
  mat(1,3)=bx*conj(bz)*(-R/Z);
  mat(2,1)=by*conj(bx)*(-R/Z);
  mat(2,2)=by*conj(by)*(-R/Z);
  mat(2,3)=by*conj(bz)*(-R/Z);
  mat(3,1)=bz*conj(bx)*(-R/Z);
  mat(3,2)=bz*conj(by)*(-R/Z);
  mat(3,3)=bz*conj(bz)*(-R/Z);
  } else
  {// quasielastic scattering needs epsilon * nm / KT ....
  jjkt=0.6666667*JJ*(JJ+1)/KBT;
  mat(1,1)=bx*conj(bx)*jjkt;
  mat(1,2)=bx*conj(by)*jjkt;
  mat(1,3)=bx*conj(bz)*jjkt;
  mat(2,1)=by*conj(bx)*jjkt;
  mat(2,2)=by*conj(by)*jjkt;
  mat(2,3)=by*conj(bz)*jjkt;
  mat(3,1)=bz*conj(bx)*jjkt;
  mat(3,2)=bz*conj(by)*jjkt;
  mat(3,3)=bz*conj(bz)*jjkt;
  }
 }
 else
 { delta=-SMALL; // tn=1 ... transition within the same level
   if(X==1.0){jjkt=JJ*(2*JJ*JJ+3*JJ+1)/3/KBT/(2*JJ+1);}
   else {if(X>1e50)
         {jjkt=JJ*JJ*KBT;}
         else 
         {jjkt=(1-2*JJ-2*JJ*JJ)/XJ;
         jjkt+=JJ*JJ/X/XJ;
	 jjkt+=(JJ*JJ+2*JJ+1)*X/XJ;
	 jjkt-=(JJ+1)*(JJ+1)*XJ;
	 jjkt+=(2*JJ*JJ+2*JJ-1)*XJ*X;
	 jjkt-=JJ*JJ*XJ*X*X;
	 jjkt*=X/(1-X)/(1-X);
	 jjkt/=(1/XJ-X*XJ)*KBT;
         }
        }
 // now lets calculate mat
 mat(1,1)=gjmbH(1)*gjmbH(1)*jjkt;
 mat(1,2)=gjmbH(1)*gjmbH(2)*jjkt;
 mat(1,3)=gjmbH(1)*gjmbH(3)*jjkt;
 mat(2,1)=gjmbH(2)*gjmbH(1)*jjkt;
 mat(2,2)=gjmbH(2)*gjmbH(2)*jjkt;
 mat(2,3)=gjmbH(2)*gjmbH(3)*jjkt;
 mat(3,1)=gjmbH(3)*gjmbH(1)*jjkt;
 mat(3,2)=gjmbH(3)*gjmbH(2)*jjkt;
 mat(3,3)=gjmbH(3)*gjmbH(3)*jjkt;
 }
if (pr==1) printf ("delta=%4.6g meV\n",delta);

return 2;
// brillouin function has 2 effective transitions
}
