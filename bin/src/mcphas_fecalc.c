 /*****************************************************************************/
// here the free energy is calculated for a given (initial) spinconfiguration
// using the meanfield algorithm / Monte Carlo 

double evalfe(double & U,double & Eelastic,spincf & sps,mfcf & mf,inipar & ini, par & inputpars,double & T, Vector * lnzi, Vector * ui)
// calculate free energy fe and energy u
{ Vector d1(1,inputpars.cs.nofcomponents),meanfield(1,inputpars.cs.nofcomponents);
int i,j,k,l,m1,s;double fe;
 fe=0;U=0; // initialize fe and u
for (i=1;i<=sps.na();++i){for (j=1;j<=sps.nb();++j){for (k=1;k<=sps.nc();++k)
{s=sps.in(i,j,k);
 for(l=1;l<=inputpars.cs.nofatoms;++l)
 {fe-=KB*T*lnzi[s][l];// sum up contributions from each ion
  if(T==0)fe+=lnzi[s][l]; // for MC just sum lnzi
  U+=ui[s][l];//fprintf(stdout,"ui(%i,%i)=%g ",s,l,ui[s][l]);
// correction term
  for(m1=1;m1<=inputpars.cs.nofcomponents;++m1)
   {d1[m1]=sps.m(i,j,k)[inputpars.cs.nofcomponents*(l-1)+m1];
  meanfield[m1]=mf.mf(i,j,k)[inputpars.cs.nofcomponents*(l-1)+m1];}
 
// add correction term
  fe+=0.5*(meanfield*d1);
  U+=0.5*(meanfield*d1);

 // printf ("Hi=%g Hj=%g Hk=%g ma=%g mb=%g mc=%g \n", meanfield[1], meanfield[2], meanfield[3], d1[1], d1[2], d1[3]);
 }
}}}
fe/=(double)sps.n(); //normalise to primitiv crystal unit cell
U/=(double)sps.n();//fprintf(stdout,"fe=%g\n",fe);

if(ini.doeps){Eelastic=sps.epsilon*inputpars.Cel*sps.epsilon;
              fe+=Eelastic;U+=Eelastic;              
              fe-=sps.epsilon*mf.epsmf;
              U-=sps.epsilon*mf.epsmf;
} // add elastic energy and magnetoelastic energy

fe/=sps.nofatoms; //normalise to meV/ion (ion=subsystem)
U/=sps.nofatoms;
Eelastic/=sps.nofatoms;
return fe;
 }


void corrfunc(Matrix & jj,int n, int l,par & inputpars,spincf & sps)
// jj a nofcomponents x nofcomponents Matrix
// l ... sublattic (atom) index
// n ... neighbour number in neighbour list
{int i,j,k,i1,j1,k1,l1,i2,j2; Vector xyz(1,3);Vector d(1,3);Vector d_rint(1,3);div_t result;
      //calculate spincorrelation function of neighbour n of sublattice l
     l1=(*inputpars.jjj[l]).sublattice[n];  // ... yes, this is the sublattice of the neighbour
xyz=(*inputpars.jjj[l]).dn[n]+(*inputpars.jjj[l]).xyz-(*inputpars.jjj[l1]).xyz; 
							// ---> xyz is then integer denoting
							// difference vector of abc lattice unit cells
                                                        // changed 22.9.2024 MR !!!
    // 1. transform xyz to primitive lattice
     d=inputpars.rez*(const Vector&)xyz;
jj=0;
     for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); //round relative position to integer numbers to correct small numerical errors

     // go through magnetic unit cell and sum up the contribution of every atom
      for(i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k){
         // now we are at atom ijk of sublattice l: the spin is  sps.m(i,j,k)(1..3+3*(l-1))
	 // which i1 j1 k1 l1 is the neighbour in distance d ?
	
	i1=i+(int)(d_rint(1));
	j1=j+(int)(d_rint(2));
	k1=k+(int)(d_rint(3));
        while (i1<=0) {i1+=sps.na();}result=div(i1,sps.na());i1=result.rem;
        while (j1<=0) {j1+=sps.nb();}result=div(j1,sps.nb());j1=result.rem;
        while (k1<=0) {k1+=sps.nc();}result=div(k1,sps.nc());k1=result.rem;
        result=div(i1,sps.na());i1=result.rem; if(i1==0)i1=sps.na();
        result=div(j1,sps.nb());j1=result.rem; if(j1==0)j1=sps.nb();
        result=div(k1,sps.nc());k1=result.rem; if(k1==0)k1=sps.nc();

	// sum up correlation function
           for(i2=1;i2<=inputpars.cs.nofcomponents;++i2)
                {
//jj(i2+inputpars.cs.nofcomponents*inputpars.cs.nofcomponents*(l-1))+=
jj(i2,i2)+=
	        (sps.m(i,j,k)(i2+inputpars.cs.nofcomponents*(l-1)))*
		(sps.m(i1,j1,k1)(i2+inputpars.cs.nofcomponents*(l1-1))); //<JaJa>,<JbJb> ...
               }
           //    k2=inputpars.cs.nofcomponents;
           for(i2=1;i2<=inputpars.cs.nofcomponents-1;++i2)
              {for(j2=i2+1;j2<=inputpars.cs.nofcomponents;++j2)
                        {//++k2;
//jj(k2+inputpars.cs.nofcomponents*inputpars.cs.nofcomponents*(l-1))+=
jj(i2,j2)+=
			 (sps.m(i,j,k)(i2+inputpars.cs.nofcomponents*(l-1)))*
			 (sps.m(i1,j1,k1)(j2+inputpars.cs.nofcomponents*(l1-1))); //<JaJb>
			//++k2;
//jj(k2+inputpars.cs.nofcomponents*inputpars.cs.nofcomponents*(l-1))+=
jj(j2,i2)+=
			 (sps.m(i,j,k)(j2+inputpars.cs.nofcomponents*(l-1)))*
			 (sps.m(i1,j1,k1)(i2+inputpars.cs.nofcomponents*(l1-1))); //<JbJa>
			}
	      }
      }}} // divide by number of basis sets in magnetic unit cell
jj/=sps.n();
 }




 

void calc_mfijk(Vector & m,spincf & sps,int & i,int & j,int& k,int & exstr,inipar & ini,int sdim,Matrix & GG, Matrix * jj,
         par & inputpars,int & diagonalexchange,int ll=0)
{int di,dj,dk,l;if(ll>0){ // only clear the specific meanfield for atom ll
                          int lm1m3=inputpars.cs.nofcomponents*(ll-1);
                        for(int m1=1;m1<=inputpars.cs.nofcomponents;++m1){m(lm1m3+m1)=0;}
                        } 
   else{m=0;}
   for (int i1=1;i1<=sps.na();++i1){if (i1>=i){di=i1-i;}else{di=sps.na()-i+i1;}
                               for (int j1=1;j1<=sps.nb();++j1){if (j1>=j){dj=j1-j;}else{dj=sps.nb()-j+j1;}
			                                    for (int k1=1;k1<=sps.nc();++k1){if (k1>=k){dk=k1-k;}else{dk=sps.nc()-k+k1;}

    l=sps.in(di,dj,dk);//di dj dk range from 0 to to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!
                       // and index a difference between crystal unit cell positions in the
                       // magnetic supercell

////      if(r==1){fprintf(stdout,"l=%i di=%i dj=%i dk=%i\n",l,di,dj,dk);    myPrintMatrix(stdout,jj[l]);getchar();}

if(ll>0){int lm1m3=inputpars.cs.nofcomponents*(ll-1);
     // only do one atom and not all
    // here the contribution of the crystal unit cell i1 j1 k1 (i1,j1,k1 indicate the
     // position of the crystal unit cell in the magnetic supercell) to the mean field
     // of the crystal unit cell i j k is calculated by one matrix multiplication
     if (diagonalexchange==0||inputpars.cs.nofatoms>1)
     {for(int m0=1;m0<=inputpars.cs.nofcomponents;++m0)
          for(int m1=1;m1<=inputpars.cs.nofatoms*inputpars.cs.nofcomponents;++m1)
             {m(lm1m3+m0)+=jj[l](lm1m3+m0,m1)*sps.m(i1,j1,k1)(m1);
            if(exstr>0&&ini.linepsjj==0){for(int bb=1;bb<=6;++bb)
              m(lm1m3+m0)+=jj[l+(sdim+2)*bb](lm1m3+m0,m1)*sps.epsilon(bb)*sps.m(i1,j1,k1)(m1);
            }}
     }else
     {//do the diagonal elements separately to accellerate the sum
      for(int m1=1;m1<=inputpars.cs.nofcomponents;++m1)
         {m(m1)+=sps.m(i1,j1,k1)(m1)*jj[l](m1,m1);
           if(exstr>0&&ini.linepsjj==0){for(int bb=1;bb<=6;++bb)
            m(m1)+=jj[l+(sdim+2)*bb](m1,m1)*sps.epsilon(bb)*sps.m(i1,j1,k1)(m1);
            }
         }
     }
}else{ //do all atoms
     // here the contribution of the crystal unit cell i1 j1 k1 (i1,j1,k1 indicate the
     // position of the crystal unit cell in the magnetic supercell) to the mean field
     // of the crystal unit cell i j k is calculated by one matrix multiplication
     if (diagonalexchange==0||inputpars.cs.nofatoms>1)
     {m+=jj[l]*(const Vector&)sps.m(i1,j1,k1);
if(exstr>0&&ini.linepsjj==0){for(int bb=1;bb<=6;++bb)
            m+=jj[l+(sdim+2)*bb]*sps.epsilon(bb)*(const Vector&)sps.m(i1,j1,k1);
            }
     }else
     {//do the diagonal elements separately to accellerate the sum
      for(int m1=1;m1<=inputpars.cs.nofcomponents;++m1)
         {m(m1)+=sps.m(i1,j1,k1)(m1)*jj[l](m1,m1);
if(exstr>0&&ini.linepsjj==0){for(int bb=1;bb<=6;++bb)
            m(m1)+=jj[l+(sdim+2)*bb](m1,m1)*sps.epsilon(bb)*sps.m(i1,j1,k1)(m1);
            }
         }
     }
 }
    }}}
  if(ini.doeps&&ini.linepscf==0){m+=sps.epsilon*GG;}
}

void calc_spsijk(Vector & m,mfcf & mf,int & i,int & j,int & k,par & inputpars,double & T,Vector & Hex,
 Vector * lnzi,Vector * ui,int  s,int  ss,ComplexMatrix **Icalcpars,ComplexVector **stat=NULL)
{  Vector d1(1,inputpars.cs.nofcomponents);
    Vector moment(1,inputpars.cs.nofcomponents);
  for(int l=1;l<=inputpars.cs.nofatoms;++l)
  {int lm1m3;
   lm1m3=inputpars.cs.nofcomponents*(l-1);
   for(int m1=1;m1<=inputpars.cs.nofcomponents;++m1)
   {d1[m1]=mf.mf(i,j,k)[lm1m3+m1];}
if(stat==0)
 (*inputpars.jjj[l]).Icalc(moment,T,d1,Hex,lnzi[s][l],ui[s][l],(*Icalcpars[inputpars.cs.nofatoms*ss+l-1]));
else
 (*inputpars.jjj[l]).Icalc(moment,T,d1,Hex,lnzi[s][l],ui[s][l],(*Icalcpars[inputpars.cs.nofatoms*ss+l-1]),stat[l]);

   if(isnan(lnzi[s][l])){fprintf (stderr, "Icalc returns lnzi=nan for s=%i l=%i\n",s,l);exit (EXIT_FAILURE);}
   if(isnan(ui[s][l])){fprintf (stderr, "Icalc returns ui=nan for s=%i l=%i\n",s,l);exit (EXIT_FAILURE);}
   for(int m1=1;m1<=inputpars.cs.nofcomponents;++m1)
   {m(lm1m3+m1)=moment[m1];}
  }
}


double fecalc(double & U, double & Eelastic, int & r,double & spinchange,Vector Hex,double T,inipar & ini,par & inputpars,
             spincf & sps,mfcf & mf,testspincf & testspins, qvectors & testqs,physproperties * physprops)
{/*on input:
    T		Temperature[K]
    Hex		Vector of external magnetic field [T] in ijk coordinates
    inputpars	exchange and other parameters
    sps		initial spinconfiguration
    testspins	all other testspinconfigurations
  on return:
    returns free energy[meV]
    sps		selfconsistently stabilized spinconfiguration (may be different
		from initial spinconfiguration)
    u		mangetic energy[meV]

 */

 double fe,dE; // free energy
 Matrix GG(1,6,1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);
 Vector sigma(1,6); // external stress tensor in Voigt notation and units meV/pVol
 Vector diff(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms),d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
 Vector meanfield(1,inputpars.cs.nofcomponents);
                 Matrix II(1,inputpars.cs.nofcomponents,1,inputpars.cs.nofcomponents);
//Matrix III(1,3,1,3);Vector dn(1,3);int sl;

 char text[MAXNOFCHARINLINE];char outfilename [MAXNOFCHARINLINE]; // some text variable
 int i,j,k,l,s,sdim,m,n,m1;
 r=0;
 div_t result; // some modulo variable
 float    sta=1000000; // initial value of standard deviation
 float staold=2000000;
 float bigstep;
 float smallstep;
 int slowct=10;
 float stepratio=1.0;
 ++ini.nofcalls;
 spinchange=0; // initial value of spinchange
 sdim=sps.in(sps.na(),sps.nb(),sps.nc()); // dimension of spinconfigurations
 Vector  * lnzi; lnzi=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){lnzi[i]=Vector(1,inputpars.cs.nofatoms);} // partition sum for every atom
 Vector  * ui; ui=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){ui[i]=Vector(1,inputpars.cs.nofatoms);ui[i]=0;} // magnetic energy for every atom
 ComplexMatrix ** Icalcpars;Icalcpars=new ComplexMatrix*[inputpars.cs.nofatoms*sdim+2];

// for each ion in the supercell make a copy of the parstorage matrix 
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {for (l=1;l<=inputpars.cs.nofatoms;++l){
  Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1]=
    new ComplexMatrix((*inputpars.jjj[l]).Icalc_parstorage.Rlo(),
                      (*inputpars.jjj[l]).Icalc_parstorage.Rhi(),
                      (*inputpars.jjj[l]).Icalc_parstorage.Clo(),
                      (*inputpars.jjj[l]).Icalc_parstorage.Chi());
  (*Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1])=(*inputpars.jjj[l]).Icalc_parstorage;

//if((*Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1])!=(*inputpars.jjj[l]).Icalc_parstorage)
// {printf("error in matrix copy\n");exit(1);}

  }}}}
 int diagonalexchange=1;
 FILE * fout;
 time_t time_of_last_output=0;
         
 spincf  spsold(sps.na(),sps.nb(),sps.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents); // spinconf variable to store old sps
 mfcf  mfold(mf.na(),mf.nb(),mf.nc(),inputpars.cs.nofatoms,inputpars.cs.nofcomponents); // spinconf variable to store old mf
// initialize mfold with large number
   for(s=0;s<=mfold.in(mfold.na(),mfold.nb(),mfold.nc());++s){ mfold.mi(s)=1000;}
   
 spsold=sps;

if(ini.doeps){ // set coupling matrix 
for(i=1;i<=6;++i)
 for(l=1;l<=inputpars.cs.nofatoms;++l)
  for(m=1;m<=inputpars.cs.nofcomponents;++m)
   GG(i,(l-1)*inputpars.cs.nofcomponents+m)=(*(*inputpars.jjj[l]).G)(i,m);
// invert elastic constants 
 //printf("#Inverting Elastic Constants Matrix\n");
  inputpars.CelInv=inputpars.Cel.Inverse();
            
// initialize epsilon mean field to zero
  mf.epsmf=0;

// set stress tensor in units of meV/primitive unit Cell, i.e. convert input strain s1-s6 from GPa to meV/pVol
// 1GPa=1e+9Pa=1e+9J/m^3
// 1meV= 1.60218e-22 J
// 1 A= 1e-10 m
// 1meV/pVol=1.60218e-22 J/A^3 x  A^3/pVol = 1.60218e+8  J/m^3 x  A^3/pVol = 1.60218e-1 GPa x  A^3/pVol
for(i=1;i<=6;++i)sigma(i)=Hex(6+i);
sigma*=inputpars.cs.pVol()/1.60218e-1;

}

// coupling coefficients jj[](a-c) berechnen
// for (r=0;r<=sdim;++r)
int exstr=0;if(ini.ipx!=NULL){exstr=6;}
 Matrix * jj; jj= new Matrix [(sdim+2)*(1+exstr)];
 for(i=0;i<=(sdim+2)*(1+exstr)-1;++i){jj[i]=Matrix(1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms,1,inputpars.cs.nofcomponents*inputpars.cs.nofatoms);} // coupling coeff.variable
   if (jj == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

   for(s=0;s<=(sdim+2)*(1+exstr)-1;++s){jj[s]=0;} //clear jj(j,...)

   for(m=1;m<=inputpars.cs.nofatoms;++m)
   {if ((*inputpars.jjj[m]).diagonalexchange==0){diagonalexchange=0;} // if any ion has anisotropic exchange - calculate anisotropic
    if(exstr>0){if ((*(*ini.ipx).jjj[m]).diagonalexchange==0){diagonalexchange=0;}
                if ((*(*ini.ipy).jjj[m]).diagonalexchange==0){diagonalexchange=0;}
                if ((*(*ini.ipz).jjj[m]).diagonalexchange==0){diagonalexchange=0;}
               }
    for(l=1;l<=(*inputpars.jjj[m]).paranz;++l)
    {//sum up l.th neighbour interaction of atom m
                                             // atom m = sublattice m
	n=(*inputpars.jjj[m]).sublattice[l]; // n set to sublattice of neighbor l

    // determine s (index of difference between crystal unit cells in the magnetic supercell)
    // start with calculating the difference vector xyz of origins of crystal unit cells
                   // bugfix GdVO3: sign of 2nd term changed and last term added 12.12.07
     xyz=(*inputpars.jjj[m]).dn[l]+(*inputpars.jjj[m]).xyz-(*inputpars.jjj[n]).xyz;
         // distance of neighbour l
                                    // xyz of sublattice m
                                                            // xyz of sublattice n

    // transform distance vector xyz to primitive lattice
     d=inputpars.rez*(const Vector&)xyz;

     for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); //round relative position to integer numbers (to do
                                             // something sensible if not integer, i.e. if sublattice
					     // of neighbour has not been identified by par.cpp)

        i=(int)(d_rint(1));
	j=(int)(d_rint(2));
	k=(int)(d_rint(3));
        // here we have the difference between crystal unitc cells ijk in the magnetic
        // supercell given by the indices i j k: if they point out of the magnetic supercell
        // they are folded back into it in the next 3 lines: this is allowed  because it is
        // irrelevant for the mean field summation
        // where the neighbor actually sits, but only on which sublattice it sits...
        while (i<=0) {i+=sps.na();}result=div(i,sps.na());i=result.rem; // only distance is important ...
        while (j<=0) {j+=sps.nb();}result=div(j,sps.nb());j=result.rem;
        while (k<=0) {k+=sps.nc();}result=div(k,sps.nc());k=result.rem;
      // s is determined from a vector ijk connecting the different crystal unit cells
	s=sps.in(i,j,k); //ijk range here from 0 to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!

        //     myPrintMatrix(stdout,(*inputpars.jjj[m]).jij[l]);

	// sum up the contribution of the interaction parameter to the interaction matrix jj[s] to be
        // used in the meanfield calculation below
	for(i=1;i<=inputpars.cs.nofcomponents;++i){for(j=1;j<=inputpars.cs.nofcomponents;++j){
	  jj[s](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=(*inputpars.jjj[m]).jij[l](i,j);

//remark: function par:jij(l) returns exchange constants (*inputpars.jjj[1]).jij[l](1-9)
        }}

    }

if(exstr>0){

for(l=1;l<=(*(*ini.ipx).jjj[m]).paranz;++l)
    {//sum up l.th neighbour interaction of atom m
                                             // atom m = sublattice m
	n=(*(*ini.ipx).jjj[m]).sublattice[l]; // n set to sublattice of neighbor l

    // determine s (index of difference between crystal unit cells in the magnetic supercell)
    // start with calculating the difference vector xyz of origins of crystal unit cells
                   // bugfix GdVO3: sign of 2nd term changed and last term added 12.12.07
     xyz=(*(*ini.ipx).jjj[m]).dn[l]+(*(*ini.ipx).jjj[m]).xyz-(*(*ini.ipx).jjj[n]).xyz;
         // distance of neighbour l
                                    // xyz of sublattice m
                                                            // xyz of sublattice n

    // transform distance vector xyz to primitive lattice
     d=inputpars.rez*(const Vector&)xyz;

     for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); //round relative position to integer numbers (to do
                                             // something sensible if not integer, i.e. if sublattice
					     // of neighbour has not been identified by par.cpp)

        i=(int)(d_rint(1));
	j=(int)(d_rint(2));
	k=(int)(d_rint(3));
        // here we have the difference between crystal unitc cells ijk in the magnetic
        // supercell given by the indices i j k: if they point out of the magnetic supercell
        // they are folded back into it in the next 3 lines: this is allowed  because it is
        // irrelevant for the mean field summation
        // where the neighbor actually sits, but only on which sublattice it sits...
        while (i<=0) {i+=sps.na();}result=div(i,sps.na());i=result.rem; // only distance is important ...
        while (j<=0) {j+=sps.nb();}result=div(j,sps.nb());j=result.rem;
        while (k<=0) {k+=sps.nc();}result=div(k,sps.nc());k=result.rem;
      // s is determined from a vector ijk connecting the different crystal unit cells
	s=sps.in(i,j,k); //ijk range here from 0 to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!

        //     myPrintMatrix(stdout,(*inputpars.jjj[m]).jij[l]);

	// sum up the contribution of the interaction parameter to the interaction matrix jj[s] to be
        // used in the meanfield calculation below
	for(i=1;i<=inputpars.cs.nofcomponents;++i){for(j=1;j<=inputpars.cs.nofcomponents;++j){

// sum up exchange striction part dJalphabeta(ij)/dRalpha depsalphagamma Rijgamma/depsbeta
// distance vector xyz (fractional lattice coordinates) transform to ijk system
Vector Rij(1,3);dadbdc2ijk(Rij,xyz,inputpars.cs.abc);

// beta = 1  (xx)
jj[s+(sdim+2)*1](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=(*(*ini.ipx).jjj[m]).jij[l](i,j)*Rij(1);
// beta =2 (yy)
jj[s+(sdim+2)*2](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=(*(*ini.ipy).jjj[m]).jij[l](i,j)*Rij(2);
// beta =3 (zz)
jj[s+(sdim+2)*3](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=(*(*ini.ipz).jjj[m]).jij[l](i,j)*Rij(3);
// beta =4 (2yz=2zy)
jj[s+(sdim+2)*4](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=0.5*(*(*ini.ipy).jjj[m]).jij[l](i,j)*Rij(3);
jj[s+(sdim+2)*4](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=0.5*(*(*ini.ipz).jjj[m]).jij[l](i,j)*Rij(2);
// beta =5 (2xz=2zx)
jj[s+(sdim+2)*5](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=0.5*(*(*ini.ipx).jjj[m]).jij[l](i,j)*Rij(3);
jj[s+(sdim+2)*5](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=0.5*(*(*ini.ipz).jjj[m]).jij[l](i,j)*Rij(1);
// beta =6 (2xy=2yx)
jj[s+(sdim+2)*6](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=0.5*(*(*ini.ipx).jjj[m]).jij[l](i,j)*Rij(2);
jj[s+(sdim+2)*6](inputpars.cs.nofcomponents*(m-1)+i,inputpars.cs.nofcomponents*(n-1)+j)+=0.5*(*(*ini.ipy).jjj[m]).jij[l](i,j)*Rij(1);

       
//remark: function par:jij(l) returns exchange constants (*inputpars.jjj[1]).jij[l](1-9)
        }}

    }  

  }          

   }


if (ini.displayall==1)   // display spincf if button is pressed
 {   strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
     fout = fopen_errchk (outfilename, "w");
     snprintf(text,MAXNOFCHARINLINE,"fecalc:%i spins, iteration %i initial values sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
     sps.eps(fout,text);
     fclose (fout);
      fprintf(stdout,"%s\n",text);
      sps.print(stdout);
     snprintf(text,MAXNOFCHARINLINE,"fecalc:%i meanfields, iteration %i initial values sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
      fprintf(stdout,"%s\n",text);
      mf.print(stdout);
  
     sleep(200);
 }

// mf loop for selfconsistency **********************************************************
// mf loop for selfconsistency **********************************************************
// mf loop for selfconsistency **********************************************************

for (r=1;sta>ini.maxstamf;++r)
{if (spinchange>ini.maxspinchange)
    {delete []jj;delete []lnzi;delete []ui;
          for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.cs.nofatoms;++l){
      delete Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []Icalcpars;
     if (verbose==1) {fprintf(stderr,"feDIV!MAXspinchangE");}++ini.nofmaxspinchangeDIV;
     return 2*FEMIN_INI+1;}

 //1. calculate mf from sps (and calculate sta) |||||||||||||||||||||||||||||||||||||||||||
 sta=0;dE=0; if(ini.doeps)mf.epsmf=0;
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {  if(ini.doeps)mf.epsmf+=GG*sps.m(i,j,k);
  
  calc_mfijk(mf.mf(i,j,k),sps,i,j,k,exstr,ini,sdim,GG,jj,inputpars,diagonalexchange);


  diff=mf.mf(i,j,k)-mfold.mf(i,j,k);sta+=diff*diff;
 // dE-=0.5*diff*(const Vector&)sps.m(i,j,k); // here we tried to calculate dE - energy difference for the step
  diff*=stepratio;mf.mf(i,j,k)=mfold.mf(i,j,k)+diff;//step gently ... i.e. scale change of MF with stepratio
  
  }}}
  // normalize mf.epsmf to crystallographic primitive unit cell (in accordance with elastic constants!)
  mf.epsmf/=sps.n();
// if(r==1){mf.print_human_readable(stdout);}

  
  sta=sqrt(sta/sps.n()/inputpars.cs.nofatoms);
  bigstep=fmodf(ini.bigstep-0.0001,1.0);
  if (ini.bigstep>1.0){smallstep=bigstep/(ini.bigstep-bigstep);}else{smallstep=bigstep/5;}
  if (r==1) {stepratio=smallstep;dE=0;} //in first loop initialize stepratio to smallstep
  if (staold<sta&&stepratio==bigstep){stepratio=smallstep;slowct=10;}//if sta increases then set stepratio to bigstep
  if (staold>sta&&stepratio<bigstep){--slowct;if (slowct<=0)stepratio=bigstep;} // at least for 10 cycles
  staold=sta;
// dE/=(double)sps.n()*sps.nofatoms; //normalise to formula unit the energy difference for this step
// if (dE>KB*T&&r>10&&stepratio<bigstep)printf("sta=%g dE=%g r=%i stepratio=%g spinschange=%g\n",sta,dE,r,stepratio,spinchange);
// ---> printing this dE  yields the result, that dE is > KB*T always when the strucuture
// is oscillating and finally diverges because of MAXSPINCHANGE reached.

 if ((ini.maxnofmfloops<=1&&r==1)||(r==2&&ini.maxnofmfloops==2)){sta=0;} // end loop on first calculation of MF from sps if no MF looping required
else
{mfold=mf;
//2. calculate sps from mf --------------------------------------------------------------
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {diff=sps.m(i,j,k);

   calc_spsijk(sps.m(i,j,k),mf,i,j,k,inputpars,T,Hex,lnzi,ui,sps.in(i,j,k),sps.in(i-1,j-1,k-1),Icalcpars);

  diff-=sps.m(i,j,k);
  spinchange+=sqrt(diff*diff)/sps.n();
  }}}

  if(ini.doeps){// here should come the exchange striction: calculate correlation function
                // and multiply with  corresponding derivative of two ion interaction
                // --> and add to mf.epsmf
                // corrfunc(Matrix & jj,int n, int l,par & inputpars,spincf & sps)
                // jj a nofcomponents x nofcomponents Matrix
                // l ... sublattic (atom) index
                // n ... neighbour number in neighbour list
if(ini.ipx!=NULL){
                for(l=1;l<=inputpars.cs.nofatoms;++l)for(n=1;n<=(*(*ini.ipx).jjj[l]).paranz;++n)
                {//calculate spincorrelation function of neighbour n of sublattice l
                   corrfunc(II,n,l,(*ini.ipx),sps);double dldlssumx=0,dldlssumy=0,dldlssumz=0;
                   for(int dl=1;dl<=inputpars.cs.nofcomponents;++dl)
                    for(int dls=1;dls<=inputpars.cs.nofcomponents;++dls){
                        dldlssumx+=(*(*ini.ipx).jjj[l]).jij[n](dl,dls)*II(dl,dls);
                        dldlssumy+=(*(*ini.ipy).jjj[l]).jij[n](dl,dls)*II(dl,dls);
                        dldlssumz+=(*(*ini.ipz).jjj[l]).jij[n](dl,dls)*II(dl,dls);}
mf.epsmf(1)+=0.5*(*(*ini.ipx).jjj[l]).dr[n](1)*dldlssumx;  // we can take dr from ipx because ipy ipz have all the same dr !
mf.epsmf(2)+=0.5*(*(*ini.ipy).jjj[l]).dr[n](2)*dldlssumy;
mf.epsmf(3)+=0.5*(*(*ini.ipz).jjj[l]).dr[n](3)*dldlssumz;
//if((*(*ini.ipx).jjj[l]).dn[n](3)==0){III=II(1,3,1,3);dn=(*(*ini.ipx).jjj[l]).dn[n];sl=(*inputpars.jjj[l]).sublattice[n];}
mf.epsmf(4)+=0.25*(*(*ini.ipz).jjj[l]).dr[n](2)*dldlssumz;
mf.epsmf(4)+=0.25*(*(*ini.ipy).jjj[l]).dr[n](3)*dldlssumy;
mf.epsmf(5)+=0.25*(*(*ini.ipz).jjj[l]).dr[n](1)*dldlssumz;
mf.epsmf(5)+=0.25*(*(*ini.ipx).jjj[l]).dr[n](3)*dldlssumx;
mf.epsmf(6)+=0.25*(*(*ini.ipy).jjj[l]).dr[n](1)*dldlssumy;
mf.epsmf(6)+=0.25*(*(*ini.ipx).jjj[l]).dr[n](2)*dldlssumx;                    
               }}

                sps.epsilon=inputpars.CelInv*(mf.epsmf+sigma);
               }

  //treat program interrupts
  #ifdef _THREADS
  MUTEX_LOCK (&mutex_tests);
  #endif
  checkini(testspins,testqs,ini);
  #ifdef _THREADS
  MUTEX_UNLOCK (&mutex_tests);
  #endif
if (ini.displayall==1)  // if all should be displayed - write sps picture to file .spins.eps
 {
     strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
     fout = fopen_errchk (outfilename, "w");
     snprintf(text,MAXNOFCHARINLINE,"fecalc:%i spins, iteration %i sta=%g ini.maxstamf=%g spinchange=%g",sps.n(),r,sta,ini.maxstamf,spinchange);
     sps.eps(fout,text);
     fclose (fout);

      fprintf(stdout,"%s\n",text);
      sps.print(stdout);
   snprintf(text,MAXNOFCHARINLINE,"... as calculated from %i meanfields, iteration %i sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
      fprintf(stdout,"%s\n",text);
      mf.print(stdout);
  
     sleep(200);
 }
   //for verbose mode do some outputs
 if (verbose==1)
 {if (time(0)-time_of_last_output>2)
  {time_of_last_output=time(0);
   strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"fe_status.dat");
     fout = fopen_errchk (outfilename, "a");
     fe=evalfe(U,Eelastic,sps,mf,ini,inputpars, T,lnzi,ui);

  #ifndef _THREADS
   fprintf(fout,"%i %g %g %g %g %g %g\n",(int)time(0),log((double)r)/log(10.0),log(sta)/log(10.0),log(spinchange+1e-10)/log(10),stepratio,100*(double)ini.successrate/ini.nofcalls,fe);
   #else
   htcalc_input *tin; int thrid;
   if ((tin=(htcalc_input*)THRLC_GET(threadSpecificKey))==THRLC_GET_FAIL) thrid = 0; else thrid = tin->thread_id+1;
   fprintf(fout,"%i %g %g %g %g %g %g %i\n",(int)time(0),log((double)r)/log(10.0),log(sta)/log(10.0),log(spinchange+1e-10)/log(10),stepratio,100*(double)ini.successrate/ini.nofcalls,fe,thrid);
   #endif
   fclose(fout);
  }
 }
if (r>ini.maxnofmfloops){if(ini.nofrndtries>=0) 
    {delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.cs.nofatoms;++l){
      delete Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []Icalcpars;

     if (verbose==1) {fprintf(stderr,"feDIV!MAXlooP");

                     }
     ++ini.nofmaxloopDIV;
     return 2*FEMIN_INI;}
                         else{sta=0;}} // continue with Monte Carlo
}
}
// end mf loop for selfconsistency **********************************************************
// end mf loop for selfconsistency **********************************************************
// end mf loop for selfconsistency **********************************************************
//printf ("hello end of selfconsistency loop after %i iterations\n",r);
//for(int ec=1;ec<=6;++ec)printf("mf.eps(%i)=%g ",ec,mf.epsmf(ec));printf("\nsl=%i\n",sl);
//myPrintMatrix(stdout,III);
// myPrintVector(stdout,dn);



// initialize ui[s][l] for maxnofmfloops==1  (without changing sps)
if(ini.maxnofmfloops==1)for (i=1;i<=sps.na();++i)for(j=1;j<=sps.nb();++j)for(k=1;k<=sps.nc();++k)
  calc_spsijk(diff,mf,i,j,k,inputpars,T,Hex,lnzi,ui,sps.in(i,j,k),sps.in(i-1,j-1,k-1),Icalcpars);

// ***************************************************************************
// do real Monte Carlo simulation - only for negative ini.nofrndtries !!!
// ***************************************************************************
if(ini.nofrndtries<0)
{ 
// Real Monte Carlo Remarks
// - one has to take energy eigenstates and not arbitrary states of the subsystems(ions),
//   if one takes arbitrary randum quantum states, then a schottky crystal field 
//   peak in specific heat cannot be desribed, because energies are continously and not
//   discrete
// - free energy calculation is not trivial, the metropolis random walk delivers
//   not the boltzmann probabilities, only physical observables
// - careful choose the initial ui[l][s] in this loop depending on maxnofmfloops, if
//   we are to calculate physical properties or just try a monte carlo: probably
//   best approach: take spins and mean fields and find the most nearby energy eigenstate
//   for each ion to represent best the spin (coming from a stabilised mf or from a previous
//   monte carlo run)
double TT=0;
ComplexVector *** states;states= new ComplexVector ** [sdim+2];
for(i=0;i<=sdim+1;++i){states[i]=new ComplexVector * [inputpars.cs.nofatoms+1];
                       for(l=1;l<=inputpars.cs.nofatoms;++l)states[i][l]=NULL;
                      }
for (i=1;i<=sps.na();++i)for(j=1;j<=sps.nb();++j)for(k=1;k<=sps.nc();++k)
 { calc_spsijk(sps.m(i,j,k),mf,i,j,k,inputpars,TT,Hex,lnzi,ui,sps.in(i,j,k),sps.in(i-1,j-1,k-1),Icalcpars,states[sps.in(i,j,k)]);
  // subtract from ui the exchange energy term, because it refers to old exchange fields
  for(l=1;l<=inputpars.cs.nofatoms;++l){
      int lm1m3=inputpars.cs.nofcomponents*(l-1);
      for(m1=1;m1<=inputpars.cs.nofcomponents;++m1)
       ui[sps.in(i,j,k)][l]+=mf.mf(i,j,k)(lm1m3+m1)*sps.m(i,j,k)(lm1m3+m1);
                                       }
 }

// recalculate all mean fields for the initial spin configuration to obtain correct initial energy E0
for (i=1;i<=sps.na();++i)for(j=1;j<=sps.nb();++j)for(k=1;k<=sps.nc();++k)
  {calc_mfijk(mf.mf(i,j,k),sps,i,j,k,exstr,ini,sdim,GG,jj,inputpars,diagonalexchange);
  // add to ui the new - correct exchange energy term
  for(l=1;l<=inputpars.cs.nofatoms;++l){
      int lm1m3=inputpars.cs.nofcomponents*(l-1);
      for(m1=1;m1<=inputpars.cs.nofcomponents;++m1)
       ui[sps.in(i,j,k)][l]-=mf.mf(i,j,k)(lm1m3+m1)*sps.m(i,j,k)(lm1m3+m1);
                                       }
 }

// now call evalfe to get E0 correctly ..
double E0;
evalfe(E0,Eelastic,sps,mf,ini,inputpars, TT,lnzi,ui);





// energy histogramm
// number of points in histogram
#define EHIST_NOFPOINTS 200
// width in energy 
#define EHIST_WIDTH      0.02*KB*T
int hi[EHIST_NOFPOINTS];int iE;for(iE=0;iE<EHIST_NOFPOINTS;++iE)hi[iE]=0;

// this is a ring storage for histogram: hi[0] stores number of runs for E0( initial
// state energy), larger energies are in hi[1,2,...], smaller energies in hi[99,98,...]


double E=0;  
Vector d1(1,inputpars.cs.nofcomponents);
Vector Imom(1,inputpars.cs.nofcomponents);
Vector mn(1,inputpars.cs.nofatoms*inputpars.cs.nofcomponents);
Vector mo(1,inputpars.cs.nofatoms*inputpars.cs.nofcomponents);
Vector En(1,inputpars.cs.nofatoms);
U=0;


Vector totalJ(1,inputpars.cs.nofcomponents);
Vector dtotalJ(1,inputpars.cs.nofcomponents);
Vector magmom(1,3); 
 Vector dmagmom(1,3),mom(1,3);
Vector tJ(1,inputpars.cs.nofcomponents);

if(physprops!=NULL){(*physprops).totalJ=0; // total operator moment <I>
    totalJ=sps.totalJ();
                    (*physprops).m=0; // magnetic moment
    magmom=0;
    for (l=1;l<=inputpars.cs.nofatoms;++l){
    // go through magnetic unit cell and sum up the contribution of every atom
    for(i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k){
      (*inputpars.jjj[l]).mcalc(dmagmom,TT, d1,Hex,(*Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1]),states[sps.in(i,j,k)][l]);
     magmom+=dmagmom;
    }}}}
    magmom/=(double)sps.n()*(double)sps.nofatoms;
                   }
// remember to introduce into ICalc T=0 beahviour (random state energy and moment!)
int nofMC=(-ini.nofrndtries-1)*sps.n();bool keep=false;
// Monte Carlo random loop ---------------------------------------------------------
for(r=1;r<nofMC;++r) 
{// choose a random spin (i.e. primitive unit cell with nofatoms spins:)
	                 i=rndint(sps.na());
		         j=rndint(sps.nb());
		         k=rndint(sps.nc());
                         l=rndint(inputpars.cs.nofatoms);
		         
// randomize  atoms in unit cell ijk with ICalc T=0 and calculate single ion energy change
 double expEKT;s=sps.in(i,j,k);dE=0;
 int lm1m3=inputpars.cs.nofcomponents*(l-1);
 
// here we update single ion energy ui[s][l] because it has been (long ago) calculated with very ancient exchange fields
 // mf. This should be updated to account for the spin changes on other sites
 // - subtract from ui the old exchange term:
 for(m1=1;m1<=inputpars.cs.nofcomponents;++m1) ui[s][l]+=mf.mf(i,j,k)(lm1m3+m1)*sps.m(i,j,k)(lm1m3+m1);
 // update mean field for spin l (only)
 calc_mfijk(mf.mf(i,j,k),sps,i,j,k,exstr,ini,sdim,GG,jj,inputpars,diagonalexchange,l);
 // update ui[s][l] with current exchange field
  for(m1=1;m1<=inputpars.cs.nofcomponents;++m1){d1(m1)=mf.mf(i,j,k)(lm1m3+m1);
                                                ui[s][l]-=d1(m1)*sps.m(i,j,k)(lm1m3+m1);}

 if(physprops!=NULL){(*inputpars.jjj[l]).mcalc(mom,TT, d1,Hex,(*Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1]),states[s][l]);
 }
 // now create randomly a new moment Imom and En(l) to be taken as new energy ui[s][l]
  mn=sps.m(i,j,k); // initialize mn to set all nofatoms spins
 ComplexVector savs(1, (*states[s][l]).Hi());savs=(*states[s][l]);
 (*inputpars.jjj[l]).Icalc(Imom,TT,d1,Hex,lnzi[s][l],En(l),(*Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1]),states[s][l]);
  // returns the moment Imom and the energy En(l) of a random chosen Energy eigenstate, 
  // which is returned in states[s][l] (to be used for physical property calculations)
   if(isnan(lnzi[s][l])){fprintf (stderr, "Icalc returns lnzi=nan for s=%i l=%i\n",s,l);exit (EXIT_FAILURE);}
   if(isnan(En(l))){fprintf (stderr, "Icalc returns ui=nan for s=%i l=%i\n",s,l);exit (EXIT_FAILURE);}
    for(m1=1;m1<=inputpars.cs.nofcomponents;++m1){
      mn(lm1m3+m1)=Imom[m1];
     if(physprops!=NULL){dtotalJ(m1)=(Imom[m1]-sps.m(i,j,k)(lm1m3+m1))/(sps.n()*sps.nofatoms);}
        }
if(physprops!=NULL){ dmagmom=-mom;
                        (*inputpars.jjj[l]).mcalc(mom,TT, d1,Hex,(*Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1]),states[s][l]);
                         dmagmom+=mom;dmagmom*=1.0/(sps.n()*sps.nofatoms);

                        }
   

 // update energy of the ensemble
  // single ion term and  two ion terms ( two ion term is included in ui !!
 // the usual correction term with prefactor 0.5 is not necessary: we count only one spin
 // in the mean field of all others, yet want to calculate the energy change of the
 // whole supercell  )
// But: there is one problem - the "self energy" term in the interaction energy -
//      -0.5*sps.m(ijk)(l)*jj[0](l)*sps.m(ijk)(l) which changes in the new configuration
//      to -0.5*Imom*jj[0](l)*Imom ... the contribution of ion l to this term is
//      not included correctly in dE=En(l)-ui[s][l], because En(l) calculated above 
//      estimates this term as -Imom*mf=-Imom*jj[0](l)*sps.m(ijk)(l), i.e. the
//      En(l) has to be modified to include this term correctly.
//      and  ui[s][l] calculated previously
//      estimates this term as -sps.m(ijk)(l)*mf=-sps.m(ijk)(l)*jj[0](l)*sps.m(ijk)(l), i.e. the
//      ui[s][l] used for dE calculation  has also to be modified to include this term correctly.
double selfenergycorrection_En=0,selfenergycorrection_ui=0;
for(m1=1;m1<=inputpars.cs.nofcomponents;++m1)
  {double d=0,dd=0,ddd,dds; for(int m2=1;m2<=inputpars.cs.nofcomponents;++m2){
               //d+=jj[0](lm1m3+m1,lm1m3+m2)*(sps.m(i,j,k)(lm1m3+m2)-0.5*Imom(m2));
               //dd+=jj[0](lm1m3+m1,lm1m3+m2)*(sps.m(i,j,k)(lm1m3+m2)-0.5*sps.m(i,j,k)(lm1m3+m2));
               ddd=jj[0](lm1m3+m1,lm1m3+m2)*sps.m(i,j,k)(lm1m3+m2);
               dds=jj[0](lm1m3+m1,lm1m3+m2)*Imom(m2);
               d+=ddd-0.5*dds;
              // ds+=ddd-dds;
               dd+=0.5*ddd;
                                                                     }
   selfenergycorrection_En+=Imom(m1)*d;
   // Encorr+=Imom(m1)*ds;
   selfenergycorrection_ui+=sps.m(i,j,k)(lm1m3+m1)*dd;
}

    dE=En(l)+selfenergycorrection_En-ui[s][l]-selfenergycorrection_ui;
// En(l)=En(l)+Encorr; // Encorr not necessary, because new ui[s][l] should refer to the old mf !

 // monte carlo stepping according to metropolis 1953 - check if we should do the step (i.e. keep=true)
 keep=true;double xi;
 if(dE>0){if(dE/(KB*T)<HUGE_EXP)expEKT=exp(-dE/(KB*T)); else  expEKT=0;
         xi=rnd(1);if(xi>expEKT){keep=false;(*states[s][l])=savs;}
         }  
 
if(keep==true)
{
/*
{// ************* 
// alternative way of calculating dE - for check
 // subtract from ui the exchange energy term, because it refers to old exchange fields
for (int i=1;i<=sps.na();++i)for(int j=1;j<=sps.nb();++j)for(int k=1;k<=sps.nc();++k)
 {   for(int l=1;l<=inputpars.cs.nofatoms;++l){
      int lm1m3=inputpars.cs.nofcomponents*(l-1);
      for(m1=1;m1<=inputpars.cs.nofcomponents;++m1)
       ui[sps.in(i,j,k)][l]+=mf.mf(i,j,k)(lm1m3+m1)*sps.m(i,j,k)(lm1m3+m1);
                                       }
 }

 sps.m(i,j,k)=mn;
 
// recalculate all mean fields for the new spin configuration to obtain correct initial energy EE
for (int i=1;i<=sps.na();++i)for(int j=1;j<=sps.nb();++j)for(int k=1;k<=sps.nc();++k)
  {calc_mfijk(mf.mf(i,j,k),sps,i,j,k,exstr,ini,sdim,GG,jj,inputpars,diagonalexchange);
  // add to ui the new - correct exchange energy term
  for(int l=1;l<=inputpars.cs.nofatoms;++l){
      int lm1m3=inputpars.cs.nofcomponents*(l-1);
      for(m1=1;m1<=inputpars.cs.nofcomponents;++m1)
       ui[sps.in(i,j,k)][l]-=mf.mf(i,j,k)(lm1m3+m1)*sps.m(i,j,k)(lm1m3+m1);
                                       }
 }
 ui[s][l]=En(l);

double EE;evalfe(EE,Eelastic,sps,mf,ini,inputpars, TT,lnzi,ui);
double dEE=(EE-E0)*sps.n()*sps.nofatoms-E;
if(fabs(dEE-dE)>0.01){printf("r=%i dE=%g dEE=%g\n",r,dE,dEE);exit(0);}
    
// *************
}*/
ui[s][l]=En(l);sps.m(i,j,k)=mn;
E+=dE;
if(physprops!=NULL){totalJ+=dtotalJ;
                    magmom+=dmagmom;
                   }
//if(E<(-0.4-E0)*sps.n()*sps.nofatoms){printf("r=%i E=%g \n",r,E/(sps.n()*sps.nofatoms)+E0);exit(0);}

}


// update energy histogram
// determine point iE
if(physprops!=NULL)
{iE=(int)rint(E/(EHIST_WIDTH*sps.n()*sps.nofatoms));
while(iE<0)iE+=EHIST_NOFPOINTS;
while(iE>EHIST_NOFPOINTS-1)iE-=EHIST_NOFPOINTS;
    ++hi[iE];
}
 // here make averages of observables
 // energy
 U+=E;//printf("%i %g %g %g\n",r,E,xi,En(1));
if(physprops!=NULL){
 // <I>
                   (*physprops).totalJ+=totalJ;
// <M> magnetic moment
                    (*physprops).m+=magmom;
     }
 } // next Monte Carlo r ---------------------------------------------------------
 // now determine absolute Energy


 
// Energy
 if(nofMC>0)U/=nofMC;   U/=sps.n()*sps.nofatoms; U+=E0;

// Operators <I>
 if(physprops!=NULL)if(nofMC>0){
 (*physprops).totalJ=(*physprops).totalJ*(1.0/nofMC);
 (*physprops).m=(*physprops).m*(1.0/nofMC);
                    }

 // Z and fe ... DIFFICULT - I do not know how to calculate !!!??????? 
// see histogram method in Bi paper 2015
if(physprops!=NULL){
strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"mcphas_mc.hst");
fout = fopen_errchk (outfilename, "w");fprintf(fout,"# Monte Carlo Energy Histogram\n# E  p(E)\n");
    for(iE=0;iE<EHIST_NOFPOINTS;++iE)
{if(iE<EHIST_NOFPOINTS/2)E=E0+iE*EHIST_WIDTH;
 else     E=E0+(iE-EHIST_NOFPOINTS)*EHIST_WIDTH;
fprintf(fout ,"%g %g\n",E,(double)hi[iE]/nofMC);
//Z=gE*exp(-E/(KB*T))*nofMC/hi[iE];  // does not work because I do not know degeneracy gE
}
fclose (fout);
}
// if(verbose)printf(" log Z=%g ",log(Z));
// for the moment just take U
fe=U;

//printf("stopping MC run with U=%g meV /ion fe = %g meV /ion",U,fe);
for(i=0;i<=sdim+1;++i){for(l=1;l<=inputpars.cs.nofatoms;++l)if(states[i][l]!=NULL) delete states[i][l];
                       delete [] states[i];
                      }
delete [] states;
}
// ***************************************************************************
// END Monte Carlo simulation 
// ***************************************************************************
else
{// calculate free energy 
fe=evalfe(U,Eelastic,sps,mf,ini,inputpars, T,lnzi,ui);
}



if (ini.displayall==1)
 {
     strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
      fout = fopen_errchk (outfilename, "w");
        snprintf(text,MAXNOFCHARINLINE,"fecalc:%i spins, iteration %i sta=%g spinchange=%g fe=%g",sps.n(),r,sta,spinchange,fe);
      sps.eps(fout,text);
      fclose (fout);
      fprintf(stdout,"%s\n",text);
      sps.print(stdout);
       snprintf(text,MAXNOFCHARINLINE,"fecalc:%i meanfields, iteration %i sta=%g spinchange=%g fe=%g",sps.n(),r,sta,spinchange,fe);
      fprintf(stdout,"%s\n",text);
      mf.print(stdout);
      sleep(200);
  }
 delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.cs.nofatoms;++l){
 delete Icalcpars[inputpars.cs.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []Icalcpars;
++ini.successrate;
return fe;
}

