/*****************************************************************************/
// here the free energy is calculated for a given (initial) spinconfiguration
// using the meanfield algorithm

double evalfe(physproperties & physprop,spincf & sps,mfcf & mf,inipar & ini, par & inputpars,double & T, Vector * lnzi, Vector * ui)
// calculate free energy fe and energy u
{ Vector d1(1,inputpars.nofcomponents),meanfield(1,inputpars.nofcomponents);
  int i,j,k,l,m1,s;
 physprop.fe=0;physprop.u=0; // initialize fe and u
for (i=1;i<=sps.na();++i){for (j=1;j<=sps.nb();++j){for (k=1;k<=sps.nc();++k)
{s=sps.in(i,j,k);
 for(l=1;l<=inputpars.nofatoms;++l)
 {physprop.fe-=KB*T*lnzi[s][l];// sum up contributions from each ion
  physprop.u+=ui[s][l];
// correction term
  for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {d1[m1]=sps.m(i,j,k)[inputpars.nofcomponents*(l-1)+m1];
  meanfield[m1]=mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+m1];}
  // add correction term
  physprop.fe+=0.5*(meanfield*d1);
  physprop.u+=0.5*(meanfield*d1);
 // printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", meanfield[1], meanfield[2], meanfield[3], d1[1], d1[2], d1[3]);
 }
}}}
physprop.fe/=(double)sps.n(); //normalise to primitiv crystal unit cell
physprop.u/=(double)sps.n();

if(ini.doeps){physprop.Eel=sps.epsilon*inputpars.Cel*sps.epsilon;
              physprop.fe+=physprop.Eel;physprop.u+=physprop.Eel;
} // add elastic energy

physprop.fe/=sps.nofatoms; //normalise to formula unit
physprop.u/=sps.nofatoms;
physprop.Eel/=sps.nofatoms;
return physprop.fe;
 }

double fecalc(Vector  Hex,double T,inipar & ini,par & inputpars,
             spincf & sps,mfcf & mf,physproperties & physprop,testspincf & testspins, qvectors & testqs)
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
 double fe,dE; // free energy
 Matrix GG(1,6,1,inputpars.nofcomponents*inputpars.nofatoms);
 Vector diff(1,inputpars.nofcomponents*inputpars.nofatoms),d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
 Vector moment(1,inputpars.nofcomponents), d1(1,inputpars.nofcomponents),meanfield(1,inputpars.nofcomponents);
 char text[MAXNOFCHARINLINE];char outfilename [MAXNOFCHARINLINE]; // some text variable
 int i,j,k,i1,j1,k1,di,dj,dk,l,r=0,s,sdim,m,n,m1;
 div_t result; // some modulo variable
 float    sta=1000000; // initial value of standard deviation
 float staold=2000000;
 float bigstep;
 float smallstep;
 int slowct=10;
 float stepratio=1.0;
 ++nofcalls;
 float spinchange=0; // initial value of spinchange
 sdim=sps.in(sps.na(),sps.nb(),sps.nc()); // dimension of spinconfigurations
 Vector  * lnzi; lnzi=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){lnzi[i]=Vector(1,inputpars.nofatoms);} // partition sum for every atom
 Vector  * ui; ui=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){ui[i]=Vector(1,inputpars.nofatoms);} // magnetic energy for every atom
 ComplexMatrix ** Icalcpars;Icalcpars=new ComplexMatrix*[inputpars.nofatoms*sdim+2];

 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {for (l=1;l<=inputpars.nofatoms;++l){
  Icalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1]=new ComplexMatrix((*inputpars.jjj[l]).Icalc_parstorage.Rlo(),(*inputpars.jjj[l]).Icalc_parstorage.Rhi(),(*inputpars.jjj[l]).Icalc_parstorage.Clo(),(*inputpars.jjj[l]).Icalc_parstorage.Chi());
  (*Icalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1])=(*inputpars.jjj[l]).Icalc_parstorage;
  }}}}
 int diagonalexchange=1;
 FILE * fin_coq;
 time_t time_of_last_output=0;
         
 spincf  spsold(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents); // spinconf variable to store old sps
 mfcf  mfold(mf.na(),mf.nb(),mf.nc(),inputpars.nofatoms,inputpars.nofcomponents); // spinconf variable to store old mf
 spsold=sps;

if(ini.doeps){ // set coupling matrix 
for(i=1;i<=6;++i)
 for(l=1;l<=inputpars.nofatoms;++l)
  for(m=1;m<=inputpars.nofcomponents;++m)
   GG(i,(l-1)*inputpars.nofcomponents+m)=(*(*inputpars.jjj[l]).G)(i,m);
// invert elastic constants multiplied by number of atoms in magnetic unitcell - exit if not possible
 // printf("#Inverting Elastic Constants Matrix\n");
  inputpars.CelInv=(1.0/(sps.na()*sps.nb()*sps.nc()))*inputpars.Cel.Inverse();
            //
// initialize epsilon to zero
  sps.epsilon=0;mf.epsmf=0;
}

// coupling coefficients jj[](a-c) berechnen
// for (r=0;r<=sdim;++r)
 
 Matrix * jj; jj= new Matrix [(sdim+1)+1];
 for(i=0;i<=sdim+1;++i){jj[i]=Matrix(1,inputpars.nofcomponents*inputpars.nofatoms,1,inputpars.nofcomponents*inputpars.nofatoms);} // coupling coeff.variable
   if (jj == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

   // initialize mfold with zeros
   for(s=0;s<=mfold.in(mfold.na(),mfold.nb(),mfold.nc());++s){mfold.mi(s)=1000;}
   for(s=0;s<=sdim;++s){jj[s]=0;} //clear jj(j,...)

   for(m=1;m<=inputpars.nofatoms;++m)
   {if ((*inputpars.jjj[m]).diagonalexchange==0){diagonalexchange=0;} // if any ion has anisotropic exchange - calculate anisotropic
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
 {   strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
     fin_coq = fopen_errchk (outfilename, "w");
     sprintf(text,"fecalc:%i spins, iteration %i initial values sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
     sps.eps(fin_coq,text);
     fclose (fin_coq);
      fprintf(stdout,"%s\n",text);
      sps.print(stdout);
     sprintf(text,"fecalc:%i meanfields, iteration %i initial values sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
      fprintf(stdout,"%s\n",text);
      mf.print(stdout);
  
     sleep(200);
 }

// loop for selfconsistency
for (r=1;sta>ini.maxstamf;++r)
{if (r>ini.maxnofmfloops)
    {delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
      delete Icalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []Icalcpars;

     if (verbose==1) fprintf(stderr,"feDIV!MAXlooP");++nofmaxloopDIV;
     return 2*FEMIN_INI;}
 if (spinchange>ini.maxspinchange)
    {delete []jj;delete []lnzi;delete []ui;
          for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
      delete Icalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []Icalcpars;
     if (verbose==1) fprintf(stderr,"feDIV!MAXspinchangE");++nofmaxspinchangeDIV;
     return 2*FEMIN_INI+1;}

 //1. calculate mf from sps (and calculate sta)
 sta=0;dE=0; if(ini.doeps)mf.epsmf=0;
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {mf.mf(i,j,k)=0;  if(ini.doeps)mf.epsmf+=GG*sps.m(i,j,k);
  for (i1=1;i1<=sps.na();++i1){if (i1>=i){di=i1-i;}else{di=sps.na()-i+i1;}
                               for (j1=1;j1<=sps.nb();++j1){if (j1>=j){dj=j1-j;}else{dj=sps.nb()-j+j1;}
			                                    for (k1=1;k1<=sps.nc();++k1){if (k1>=k){dk=k1-k;}else{dk=sps.nc()-k+k1;}

    l=sps.in(di,dj,dk);//di dj dk range from 0 to to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!
                       // and index a difference between crystal unit cell positions in the
                       // magnetic supercell

////      if(r==1){fprintf(stdout,"l=%i di=%i dj=%i dk=%i\n",l,di,dj,dk);    myPrintMatrix(stdout,jj[l]);getchar();}

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
  if(ini.doeps){mf.mf(i,j,k)+=sps.epsilon*GG;}
  diff=mf.mf(i,j,k)-mfold.mf(i,j,k);sta+=diff*diff;
 // dE-=0.5*diff*(const Vector&)sps.m(i,j,k); // here we tried to calculate dE - energy difference for the step
  diff*=stepratio;mf.mf(i,j,k)=mfold.mf(i,j,k)+diff;//step gently ... i.e. scale change of MF with stepratio
  }}}
  mfold=mf;
  sta=sqrt(sta/sps.n()/inputpars.nofatoms);
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

//2. calculate sps from mf
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {diff=sps.m(i,j,k);s=sps.in(i,j,k);
  for(l=1;l<=inputpars.nofatoms;++l)
  {int lm1m3;
   lm1m3=inputpars.nofcomponents*(l-1);
   for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {d1[m1]=mf.mf(i,j,k)[lm1m3+m1];}
   (*inputpars.jjj[l]).Icalc(moment,T,d1,Hex,lnzi[s][l],ui[s][l],(*Icalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1]));
   for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {sps.m(i,j,k)(lm1m3+m1)=moment[m1];}
  }
  diff-=sps.m(i,j,k);
  spinchange+=sqrt(diff*diff)/sps.n();
  }}}
  if(ini.doeps){sps.epsilon=inputpars.CelInv*mf.epsmf;
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
     fin_coq = fopen_errchk (outfilename, "w");
     sprintf(text,"fecalc:%i spins, iteration %i sta=%g ini.maxstamf=%g spinchange=%g",sps.n(),r,sta,ini.maxstamf,spinchange);
     sps.eps(fin_coq,text);
     fclose (fin_coq);

      fprintf(stdout,"%s\n",text);
      sps.print(stdout);
   sprintf(text,"... as calculated from %i meanfields, iteration %i sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
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
     fin_coq = fopen_errchk (outfilename, "a");
     fe=evalfe(physprop,sps,mf,ini,inputpars, T,lnzi,ui);

  #ifndef _THREADS
   fprintf(fin_coq,"%i %g %g %g %g %g %g\n",(int)time(0),log((double)r)/log(10.0),log(sta)/log(10.0),log(spinchange+1e-10)/log(10),stepratio,100*(double)successrate/nofcalls,fe);
   #else
   htcalc_input *tin; int thrid;
   if ((tin=(htcalc_input*)THRLC_GET(threadSpecificKey))==THRLC_GET_FAIL) thrid = 0; else thrid = tin->thread_id+1;
   fprintf(fin_coq,"%i %g %g %g %g %g %g %i\n",(int)time(0),log((double)r)/log(10.0),log(sta)/log(10.0),log(spinchange+1e-10)/log(10),stepratio,100*(double)successrate/nofcalls,fe,thrid);
   #endif
   fclose(fin_coq);
  }
 }

}

//printf ("hello end of selfconsistency loop after %i iterations\n",r);
fe=evalfe(physprop,sps,mf,ini,inputpars, T,lnzi,ui);


if (ini.displayall==1)
 {
     strcpy(outfilename,"./results/.");strcpy(outfilename+11,ini.prefix);
     strcpy(outfilename+11+strlen(ini.prefix),"spins.eps");
      fin_coq = fopen_errchk (outfilename, "w");
        sprintf(text,"fecalc:%i spins, iteration %i sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
      sps.eps(fin_coq,text);
      fclose (fin_coq);
      fprintf(stdout,"%s\n",text);
      sps.print(stdout);
       sprintf(text,"fecalc:%i meanfields, iteration %i sta=%g spinchange=%g",sps.n(),r,sta,spinchange);
      fprintf(stdout,"%s\n",text);
      mf.print(stdout);
      sleep(200);
  }

 delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
 delete Icalcpars[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []Icalcpars;
++successrate;
return fe;
}

