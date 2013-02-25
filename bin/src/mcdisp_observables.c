//******************************************************************************
// FUNCTIONS FOR OBSERVABLES
//******************************************************************************

// some function to write fileheader efficiently
void writeheader(par & inputpars,FILE * fout)
{  time_t curtime;
  struct tm *loctime;
   fprintf(fout, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
  
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
  fprintf(fout, "# List of atomic positions dr1 dr2 dr3, moments m \n");
  fprintf(fout, "# Debye Waller factor (sqr(Intensity)~|sf| ~sum_i ()i exp(-2 DWFi sin^2(theta) / lambda^2)=EXP (-Wi),\n# units DWF [A^2], relation to other notations 2*DWF=B=8 pi^2 <u^2>)\n");
  fprintf(fout, "#  and  Lande factors total angular momentum J (=0 if dipole approximation is used) <j0> and <j2> formfactor\n# coefficients\n");
  fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3] DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  fprintf(fout, "#                         ...with j||b, k||(a x b) and i normal to k and j\n");
 
 for (int i = 1;i<=inputpars.nofatoms;++i)
 {if((double)(i)/50==(double)(i/50))
  {
  fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3] DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  
  }
  fprintf(fout, "# %6.3f %6.3f %6.3f %6.3f %6.3f ",(*inputpars.jjj[i]).xyz(1),(*inputpars.jjj[i]).xyz(2),(*inputpars.jjj[i]).xyz(3),(*inputpars.jjj[i]).DWF,(*inputpars.jjj[i]).gJ);
  if((*inputpars.jjj[i]).Np(1)!=0){
  fprintf(fout," - formfactor calc. from radial wave function parameters in %s: <jl(Q)> considered up to l=%i",(*inputpars.jjj[i]).sipffilename,(*inputpars.jjj[i]).jl_lmax);
  }
  else
  {
  for (int j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",(*inputpars.jjj[i]).magFFj0(j));}
  for (int j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",(*inputpars.jjj[i]).magFFj2(j));}
  }
  fprintf(fout, "\n");
 }
}


void fillE(int jmin,int i,int j, int k, int l,int ev_dim,ComplexVector & observable_coeff1,
    par & inputpars,ComplexMatrix & observable_Mijkl ,mdcf & md,
    Vector & observable_gamma,double & observable_gamman,ComplexMatrix &observable_Uijkl,int maxiter,float * nn,
    inimcdis & ini,Vector & gamma,ComplexMatrix & Eobservable)
{complex<double> imaginary(0,1);  int sort;
    // here we if required calculate the higher dimension matrix used to do the
       // extension of chi to higher value of (uncoupled) nofcomponents in intcalc_approx ... needed for observablefluctuations, extended eigenvectors ...
          {observable_coeff1=0;observable_coeff1(1)=1e-10;}
          observable_Mijkl = observable_coeff1^observable_coeff1;observable_gamman=Norm2(observable_coeff1);observable_coeff1/=sqrt(observable_gamman);
          myEigenSystemHermitean (observable_Mijkl,observable_gamma,observable_Uijkl,sort=1,maxiter);
          if (fabs(observable_gamman-observable_gamma(ev_dim))>SMALL){fprintf(stderr,"ERROR eigenvalue of extended single ion matrix extM inconsistent: analytic value gamma= %g numerical diagonalisation of extM gives gamma= %g\n",observable_gamman,observable_gamma(ev_dim));
                           exit(EXIT_FAILURE);}
          for(int ii=observable_Uijkl.Rlo(); ii<=observable_Uijkl.Rhi(); ii++){if (fabs(abs(observable_coeff1(ii))-abs(observable_Uijkl(ii,ev_dim)))>SMALL)
                                                {fprintf(stderr,"ERROR eigenvector of extended single ion matrix M inconsistent\n");
                                                 myPrintComplexVector(stderr,observable_coeff1);observable_coeff1=observable_Uijkl.Column(ev_dim);myPrintComplexVector(stderr,observable_coeff1);exit(EXIT_FAILURE);}
                                                 observable_Uijkl(ii,ev_dim)=observable_coeff1(ii);}
         if (nn[6]<0) // if transition energy is less than zero do a conjugation of the matrix
	 {   for(int ii=observable_Uijkl.Rlo(); ii<=observable_Uijkl.Rhi(); ii++) for(int jj=observable_Uijkl.Clo(); jj<=observable_Uijkl.Chi(); jj++) observable_Uijkl[ii][jj]=conj(observable_Uijkl[ii][jj]);
	 }
        // here we should also go for complex conjugate for the vector
         complex <double> observable_gammas;
         if (gamma(ini.nofcomponents)>=0&&fabs(gamma(ini.nofcomponents-1))<SMALL) 
                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
                           {if (nn[6]>SMALL)
			    {observable_gammas=sqrt(observable_gamma(ev_dim));}
			    else if (nn[6]<-SMALL)
                            {observable_gammas=imaginary*sqrt(observable_gamma(ev_dim));}
 			    else
			    { //quasielastic line needs gamma=SMALL .... because Mijkl and therefore gamma have been set to 
			      // wn/kT instead of wn-wn'=SMALL*wn/kT (in jjjpar.cpp -mdcalc routines)
			      //set fix delta but keep sign
			          if (nn[6]>0){observable_gammas=sqrt(observable_gamma(ev_dim));}
				  else        {observable_gammas=imaginary*sqrt(observable_gamma(ev_dim));}
			    }
			   }else 
                           {fprintf(stderr,"ERROR eigenvalue of single ion matrix <0: ev1=%g ev2=%g ev3=%g ... evn=%g\n",gamma(1),gamma(2),gamma(3),gamma(ini.nofcomponents));
                            exit(EXIT_FAILURE);}
         //printf("observable_gamma %g %+g i\n",real(observable_gammas),imag(observable_gammas));
        for(int k1=1;k1<=ev_dim;++k1){Eobservable(index_s(i,j,k,l,jmin,md,ini),k1)=observable_gammas*observable_Uijkl(k1,ev_dim);}
}

FILE * evfileinit(const char * filename,const char*filemode,par & inputpars,const char * filetype,int ev_dim)
  {FILE * fout = fopen_errchk (filename,filemode);
   writeheader(inputpars,fout);  fprintf(fout,"#!<--mcphas.mcdisp.%s-->\n",filetype);
          fprintf (fout, "#!spins_wave_amplitude=1.0\n");
          fprintf (fout, "#!spins_show_ellipses=1.0\n");
          fprintf (fout, "#!spins_show_static_moment_direction=1.0\n");
          fprintf (fout, "#!extended_eigenvector_dimension=%i\n",ev_dim); 
          fprintf (fout, "#!dispersion displayytext=E(meV)\n#Ha[T] Hb[T] Hc[T] T[K] h k l Q[A^-1] energy[meV] int_dipapprFF) [barn/sr/f.u.] int_beyonddipappr [barn/sr/f.u.]  f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
  return fout;
  }

void print_ev(FILE * fout,int i,inimcdis & ini,Vector & hkl,double QQ,Vector & En,Vector & ints,Vector & intsbey,mfcf & qee_real,mfcf & qee_imag)
                               {
                     fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g  %4.4g  %4.4g\n",myround(ini.Hext(1)),myround(ini.Hext(2)),myround(ini.Hext(3)),myround(ini.T),myround(hkl(1)),myround(hkl(2)),myround(hkl(3)),
                              myround(QQ),myround(En(i)),myround(1e-8,ints(i)),myround(1e-8,intsbey(i)));
                     fprintf (fout, "#eigenvector real part\n");
                     qee_real.print(fout); // here we printout the eigenvector of the excitation
                     fprintf (fout, "#eigenvector imaginary part\n");
                     qee_imag.print(fout); // 
                     fprintf (fout, "#\n");
                                }

