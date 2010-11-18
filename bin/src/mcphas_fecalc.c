/*****************************************************************************/
// here the free energy is calculated for a given (initial) spinconfiguration
// using the meanfield algorithm
double fecalc(Vector  Hex,double T,par & inputpars,
             spincf & sps,mfcf & mf,double & u,testspincf & testspins, qvectors & testqs)
{/*on input:
    T		Temperature[K]
    Hex		Vector of external magnetic field [T]
    inputpars	exchange and other parameters
    sps		initial spinconfiguration
    testspins	all other testspinconfigurations
  on return:
    returns free energy[meV]
    sps		selfconsistently stabilized spinconfiguration (may be different
		from initial spinconfiguration)
    u		mangetic energy[meV]

 */
 double fe; // free energy
 Vector diff(1,inputpars.nofcomponents*inputpars.nofatoms),d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
 Vector meanfield(1,inputpars.nofcomponents),moment(1,inputpars.nofcomponents),d1(1,inputpars.nofcomponents);
 char text[60]; // some text variable
 int i,j,k,i1,j1,k1,di,dj,dk,l,r,s,sdim,m,n,m1;
 div_t result; // some modulo variable
 float    sta=1000000; // initial value of standard deviation
 float staold=2000000;
 float bigstep;
 float smallstep;
 int slowct=10;
 float stepratio=1.0;
 static int successrate=0;
 static int nofcalls=0;
 ++nofcalls;
 float spinchange=0; // initial value of spinchange
 sdim=sps.in(sps.na(),sps.nb(),sps.nc()); // dimension of spinconfigurations
 Vector  * lnzi; lnzi=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){lnzi[i]=Vector(1,inputpars.nofatoms);} // partition sum for every atom
 Vector  * ui; ui=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){ui[i]=Vector(1,inputpars.nofatoms);} // magnetic energy for every atom
 ComplexMatrix ** mcalcpars;mcalcpars=new ComplexMatrix*[inputpars.nofatoms*sdim+2];
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {for (l=1;l<=inputpars.nofatoms;++l){
  mcalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1]=new ComplexMatrix((*inputpars.jjj[l]).mcalc_parstorage.Rlo(),(*inputpars.jjj[l]).mcalc_parstorage.Rhi(),(*inputpars.jjj[l]).mcalc_parstorage.Clo(),(*inputpars.jjj[l]).mcalc_parstorage.Chi());
  (*mcalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1])=(*inputpars.jjj[l]).mcalc_parstorage;
  }}}}
 int diagonalexchange=1;
 FILE * fin_coq;
 time_t time_of_last_output=0;

 spincf  spsold(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents); // spinconf variable to store old sps
 mfcf  mfold(mf.na(),mf.nb(),mf.nc(),inputpars.nofatoms,inputpars.nofcomponents); // spinconf variable to store old mf
 spsold=sps;

// coupling coefficients jj[](a-c) berechnen
// for (r=0;r<=sdim;++r)

 Matrix * jj; //if (inputpars.diagonalexchange()==0){i=9;}else{i=3;}
  jj= new Matrix [(sdim+1)+1];for(i=0;i<=sdim+1;++i){jj[i]=Matrix(1,inputpars.nofcomponents*inputpars.nofatoms,1,inputpars.nofcomponents*inputpars.nofatoms);} // coupling coeff.variable
   if (jj == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

   // initialize mfold with zeros
   for(s=0;s<=mfold.in(mfold.na(),mfold.nb(),mfold.nc());++s){mfold.mi(s)=0;}
   for(s=0;s<=sdim;++s){jj[s]=0;} //clear jj(j,...)

   for(m=1;m<=inputpars.nofatoms;++m)
   {if ((*inputpars.jjj[m]).diagonalexchange==0){diagonalexchange=0;} // if any ion has anisotropic exchange - calculate anisotropic
    for(l=1;l<=(*inputpars.jjj[m]).paranz;++l)
    {//sum up l.th neighbour interaction of atom m
                                             // atom m = sublattice m
	n=(*inputpars.jjj[m]).sublattice[l]; // n set to sublattice of neighbor l

    // determine s (index of difference between crystal unit cells in the magnetic supercell)
    // start with calculating the difference vector of origins of crystal unit cells
                   // bugfix GdVO3: sign of 2nd term changed and last term added 12.12.07
     xyz=(*inputpars.jjj[m]).dn[l]+(*inputpars.jjj[m]).xyz-(*inputpars.jjj[n]).xyz;

    // transform distance vector to primitive lattice
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
        while (i<=0) i+=sps.na();result=div(i,sps.na());i=result.rem; // only distance is important ...
        while (j<=0) j+=sps.nb();result=div(j,sps.nb());j=result.rem;
        while (k<=0) k+=sps.nc();result=div(k,sps.nc());k=result.rem;
      // s is determined from a vector ijk connecting the different crystal unit cells
	s=sps.in(i,j,k); //ijk range here from 0 to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!

        //     myPrintMatrix(stdout,(*inputpars.jjj[m]).jij[l]);

	// sum up the contribution of the interaction parameter to the interaction matrix jj[s] to be
        // used in the meanfield calculation below
	for(i=1;i<=inputpars.nofcomponents;++i){for(j=1;j<=inputpars.nofcomponents;++j){
	  jj[s](inputpars.nofcomponents*(m-1)+i,inputpars.nofcomponents*(n-1)+j)+=(*inputpars.jjj[m]).jij[l](i,j);

	//remark: function par:jij(l) returns exchange constants (*inputpars.jjj[1]).jij[l](1-9)
        }}

    }
   }



if (ini.displayall==1)   // display spincf if button is pressed
 {
     fin_coq = fopen_errchk ("./results/.spins.eps", "w");
     sprintf(text,"fecalc:%i spins, iteration 0",sps.n());
     sps.eps(fin_coq,text);
     fclose (fin_coq);
     sleep(2);

// sps.display(text);
 }


// loop for selfconsistency
for (r=1;sta>ini.maxstamf;++r)
{if (r>ini.maxnofmfloops)
    {delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
      delete mcalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []mcalcpars;

     if (verbose==1) fprintf(stderr,"feDIV!MAXlooP");
     return 20000;}
 if (spinchange>ini.maxspinchange)
    {delete []jj;delete []lnzi;delete []ui;
          for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
      delete mcalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []mcalcpars;
     if (verbose==1) fprintf(stderr,"feDIV!MAXspinchangE");
     return 20001;}

 //1. calculate mf from sps (and calculate sta)
 sta=0;
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {mf.mf(i,j,k)=0;
  for (l=1;l<=inputpars.nofatoms;++l){
    if(inputpars.gJ(l)==0)              {
     for (i1=1;i1<=6&&i1<=inputpars.nofcomponents;++i1){
            if(i1==2||i1==4||i1==6){ mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+i1]=Hex(i1/2)*MU_B;}
	    else                   { mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+i1]=2*Hex((i1+1)/2)*MU_B;}
  				                       }
            				}
    else                                {
     for (i1=1;i1<=3&&i1<=inputpars.nofcomponents;++i1){
               mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+i1]=Hex(i1)*inputpars.gJ(l)*MU_B;
  				                       }
					}
				     }
  for (i1=1;i1<=sps.na();++i1){if (i<i1){di=i-i1+sps.na();}else{di=i-i1;}
                               for (j1=1;j1<=sps.nb();++j1){if (j<j1){dj=j-j1+sps.nb();}else{dj=j-j1;}
			                                    for (k1=1;k1<=sps.nc();++k1){if (k<k1){dk=k-k1+sps.nc();}else{dk=k-k1;}
    l=sps.in(di,dj,dk);//di dj dk range from 0 to to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!
                       // and index a difference between crystal unit cell positions in the
                       // magnetic supercell

     // here the contribution of the crystal unit cell i1 j1 k1 (i1,j1,k1 indicate the
     // position of the crystal unit cell in the magnetic supercell) to the mean field
     // of the crystal unit cell i j k is calculated by one matrix multiplication
     if (diagonalexchange==0||inputpars.nofatoms>1)
     {mf.mf(i,j,k)+=jj[l]*(const Vector&)sps.m(i1,j1,k1);
     }else
     {//do the diagonal elements separately to accellerate the sum
      for(m1=1;m1<=inputpars.nofatoms*inputpars.nofcomponents;++m1)
         {mf.mf(i,j,k)(m1)+=sps.m(i1,j1,k1)(m1)*jj[l](m1,m1);}
     }
    }}}
  diff=mf.mf(i,j,k)-mfold.mf(i,j,k);sta+=diff*diff;
  diff*=stepratio;mf.mf(i,j,k)=mfold.mf(i,j,k)+diff;//step gently ... i.e. scale change of MF with stepratio
  }}}
  mfold=mf;
  sta=sqrt(sta/sps.n()/inputpars.nofatoms);
  //printf("sta=%g\n",sta);
  bigstep=fmodf(ini.bigstep-0.0001,1.0);
  if (ini.bigstep>1.0){smallstep=bigstep/(ini.bigstep-bigstep);}else{smallstep=bigstep/5;}
  if (r==1) {stepratio=smallstep;} //in first loop initialize stepratio to smallstep
  if (staold<sta&&stepratio==bigstep){stepratio=smallstep;slowct=10;}//if sta increases then set stepratio to bigstep
  if (staold>sta&&stepratio<bigstep){--slowct;if (slowct<=0)stepratio=bigstep;} // at least for 10 cycles
  staold=sta;

//2. calculate sps from mf
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {diff=sps.m(i,j,k);s=sps.in(i,j,k);
  for(l=1;l<=inputpars.nofatoms;++l)
  {int lm1m3;
   lm1m3=inputpars.nofcomponents*(l-1);
   for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {d1[m1]=mf.mf(i,j,k)[lm1m3+m1];}
   (*inputpars.jjj[l]).mcalc(moment,T,d1,lnzi[s][l],ui[s][l],(*mcalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1]));
   for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {sps.m(i,j,k)(lm1m3+m1)=moment[m1];}
  }
  diff-=sps.m(i,j,k);
  spinchange+=sqrt(diff*diff)/sps.n();
  }}}

  //treat program interrupts
  #ifdef _THREADS
  MUTEX_LOCK (&mutex_tests);
  #endif
  checkini(testspins,testqs);
  #ifdef _THREADS
  MUTEX_UNLOCK (&mutex_tests);
  #endif
if (ini.displayall==1)  // if all should be displayed - write sps picture to file .spins.eps
 {
      fin_coq = fopen_errchk ("./results/.spins.eps", "w");
     sprintf(text,"fecalc:%i spins, iteration %i sta=%g",sps.n(),r,sta);
     sps.eps(fin_coq,text);
     fclose (fin_coq);
     sleep(2);
 }
   //for verbose mode do some outputs
 if (verbose==1)
 {if (time(0)-time_of_last_output>2)
  {time_of_last_output=time(0);
   fin_coq= fopen_errchk ("./results/.fe_status.dat","a");
   #ifndef _THREADS
   fprintf(fin_coq,"%i %g %g %g %g %g\n",(int)time(0),log((double)r)/log(10.0),log(sta)/log(10.0),spinchange,stepratio,100*(double)successrate/nofcalls);
   #else
   htcalc_input *tin; int thrid;
   if ((tin=(htcalc_input*)THRLC_GET(threadSpecificKey))==THRLC_GET_FAIL) thrid = 0; else thrid = tin->thread_id+1;
   fprintf(fin_coq,"%i %g %g %g %g %g %i\n",(int)time(0),log((double)r)/log(10.0),log(sta)/log(10.0),spinchange,stepratio,100*(double)successrate/nofcalls,thrid);
   #endif
   fclose(fin_coq);
  }
 }

}

//printf ("hello end of mf procedure");

// calculate free energy fe and energy u
fe=0;u=0; // initialize fe and u
for (i=1;i<=sps.na();++i){for (j=1;j<=sps.nb();++j){for (k=1;k<=sps.nc();++k)
{s=sps.in(i,j,k);
 for(l=1;l<=inputpars.nofatoms;++l)
 {fe-=KB*T*lnzi[s][l];// sum up contributions from each ion
  u+=ui[s][l];
// correction term
  for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {d1[m1]=sps.m(i,j,k)[inputpars.nofcomponents*(l-1)+m1];
  meanfield[m1]=mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+m1];}
  // subtract external field (only necessary for magnetic field, not for quadrupolar fields,
  // because the Cf parameters are treated separately in mcalc and not as part of the quadrupolar
  // field)
  if(inputpars.gJ(l)==0)
  {for(m1=1;m1<=6&&m1<=inputpars.nofcomponents;++m1)

                           {if(m1==2||m1==4||m1==6) {meanfield[m1]-=Hex[m1/2]*MU_B;}

                            else                    {meanfield[m1]-=2*Hex[(m1+1)/2]*MU_B;}

                           }

                          }
 else

                          {
for(m1=1;m1<=3&&m1<=inputpars.nofcomponents;++m1)
                           {meanfield[m1]-=Hex[m1]*inputpars.gJ(l)*MU_B;}

                          }

  // add correction term
  fe+=0.5*(meanfield*d1);
  u+=0.5*(meanfield*d1);
//  printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", H[1], H[2], H[3], m[1], m[2], m[3]);
 }
}}}
fe/=(double)sps.n()*sps.nofatoms; //normalise to formula unit
u/=(double)sps.n()*sps.nofatoms;

if (ini.displayall==1)
 {
      fin_coq = fopen_errchk ("./results/.spins.eps", "w");
       sprintf(text,"fecalc:%i spins, iteration %i, fe=%gmeV",sps.n(),i,fe);
       sps.eps(fin_coq,text);
      fclose (fin_coq);
      sleep(2);
// sps.display(text);
  }

 delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
 delete mcalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []mcalcpars;
++successrate;
return fe;
}

