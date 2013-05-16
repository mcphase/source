//**************************************************************************************
// populates Matrix Eobservable  ... maybe put this to mcdisp.c ....
//**************************************************************************************
void fillE(int jmin,int i,int j, int k, int l,int ev_dim,ComplexVector & observable_coeff1,
    par & inputpars,ComplexMatrix & observable_Mijkl ,mdcf & md,
    Vector & observable_gamma,double & observable_gamman,ComplexMatrix &observable_Uijkl,int maxiter,float * nn,
    inimcdis & ini,Vector & gamma,ComplexMatrix & Eobservable)
{complex<double> imaginary(0,1);  int sort;
    // here we if required calculate the higher dimension matrix used to do the
       // extension of chi to higher value of (uncoupled) nofcomponents in intcalc_approx ... needed for observablefluctuations, extended eigenvectors ...
       //   {observable_coeff1=0;observable_coeff1(1)=1e-10;}
          observable_Mijkl = observable_coeff1^observable_coeff1;observable_gamman=Norm2(observable_coeff1);observable_coeff1/=sqrt(observable_gamman);
          myEigenSystemHermitean (observable_Mijkl,observable_gamma,observable_Uijkl,sort=1,maxiter);
          if (fabs(observable_gamman-observable_gamma(ev_dim))>SMALL_QUASIELASTIC_ENERGY){fprintf(stderr,"WARNING eigenvalue of extended single ion matrix observable_Mijkl inconsistent: analytic value gamma= %g numerical diagonalisation of extM gives gamma= %g\n",observable_gamman,observable_gamma(ev_dim));}
          for(int ii=observable_Uijkl.Rlo(); ii<=observable_Uijkl.Rhi(); ii++){if (fabs(abs(observable_coeff1(ii))-abs(observable_Uijkl(ii,ev_dim)))>SMALL_QUASIELASTIC_ENERGY)
                                                {fprintf(stderr,"ERROR eigenvector of extended single ion matrix M inconsistent\n");
                                                 myPrintComplexVector(stderr,observable_coeff1);observable_coeff1=observable_Uijkl.Column(ev_dim);myPrintComplexVector(stderr,observable_coeff1);exit(EXIT_FAILURE);}
                                                 observable_Uijkl(ii,ev_dim)=observable_coeff1(ii);}
         if (nn[6]<0) // if transition energy is less than zero do a conjugation of the matrix
	 {   for(int ii=observable_Uijkl.Rlo(); ii<=observable_Uijkl.Rhi(); ii++) for(int jj=observable_Uijkl.Clo(); jj<=observable_Uijkl.Chi(); jj++) observable_Uijkl[ii][jj]=conj(observable_Uijkl[ii][jj]);
	 }
        // here we should also go for complex conjugate for the vector
         complex <double> observable_gammas;
         if (gamma(ini.nofcomponents)>=0&&fabs(gamma(ini.nofcomponents-1))<SMALL_QUASIELASTIC_ENERGY) 
                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystemHermitean
			   // sorts the eigenvalues according to ascending order !!!
                           {if (nn[6]>SMALL_QUASIELASTIC_ENERGY)
			    {observable_gammas=sqrt(observable_gamma(ev_dim));}
			    else if (nn[6]<-SMALL_QUASIELASTIC_ENERGY)
                            {observable_gammas=imaginary*sqrt(observable_gamma(ev_dim));}
 			    else
			    { //quasielastic line needs gamma=SMALL_QUASIELASTIC_ENERGY .... because Mijkl and therefore gamma have been set to 
			      // wn/kT instead of wn-wn'=SMALL_QUASIELASTIC_ENERGY*wn/kT (in jjjpar.cpp -mdcalc routines)
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

//******************************************************************************
// FUNCTIONS FOR OUTPUT of CALCULATION RESULTS
//******************************************************************************

// some function to write fileheader efficiently
void writeheader(par & inputpars,FILE * fout)
{  time_t curtime;
  struct tm *loctime;
   fprintf(fout, "#output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
  
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
  fprintf(fout, "# List of atomic positions dr1 dr2 dr3\n");
  fprintf(fout, "# Debye Waller factor (sqr(Intensity)~|sf| ~sum_i ()i exp(-2 DWFi sin^2(theta) / lambda^2)=EXP (-Wi),\n# units DWF [A^2], relation to other notations 2*DWF=B=8 pi^2 <u^2>)\n");
  fprintf(fout, "#  and  Lande factors total angular momentum J, <j0> and <j2> formfactor\n# coefficients\n");
  fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3] DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  
 for (int i = 1;i<=inputpars.nofatoms;++i)
 {if((double)(i)/50==(double)(i/50))
  {
  fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3] DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  
  }
  fprintf(fout, "# %6.3f %6.3f %6.3f %6.3f %6.3f ",(*inputpars.jjj[i]).xyz(1),(*inputpars.jjj[i]).xyz(2),(*inputpars.jjj[i]).xyz(3),(*inputpars.jjj[i]).DWF,(*inputpars.jjj[i]).gJ);
  (*inputpars.jjj[i]).FFinfo(fout);
 }
}

void writeheaders(FILE * foutqom,FILE * foutqei,FILE * foutdstot,FILE * foutds,par & inputpars,inimcdis & ini,int & calc_rixs,int & do_Erefine)                   
{//open output files and write headers for output files mcdisp.qei qex qom dsigma.tot dsigma
  fprintf(foutqom,"#!<--mcphas.mcdisp.qom-->\n");

 if(calc_rixs){fprintf(foutqei,"#!<--mcphas.mcdisp.qex-->\n");
               writeheader(inputpars,foutqei);
               fprintf(foutqei,"#RIXS intensity components:azimuth is defined as in Longfield et al. PRB 66 054417 (2002)\n"
                               "#coordinate system u1,u2,u3: the scattering plane, defined by the\n"
                               "#direction of the incident and final wave vectors k and k', contains u1 lying\n"
                               "#perpendicular to Q and in the sense of k, and u3 parallel to the scattering\n"
                               "#vector Q= k-k'.\n"
                               "#angles for azimuth=0: alpha_i=angle(ai,u3), delta_i=angle(ai_perp,u1)\n"
                               "#(where a1,a2,a3=a,b,c and ai_perp is the projection of ai onto the plane\n"
                               "#perpendicular to Q. In the chosen experimental geometry\n"
                               "#azimuth=0, when a1=a points to the x-ray source.\n");
             fprintf (foutqei, "#dispersion displayytext=E(meV)\n#displaylines=false \n");
             ini.print_usrdefcolhead(foutqei);
             fprintf (foutqei,"h   k   l Q[A^-1] energy[meV] " 
                              "Irix: Isigmasigma azimuthsigmasigma   Isigmapi azsigmapi  "
                              " Ipisigma azpisigma   Ipipi azpipi          Irightright azrightright"
                              " Irightleft azreightleft Ileftright azleftright Ileftleft azlefteft(deg) "
                              "  [a.u./sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n");

            fprintf (foutqom, "#dispersion \n");ini.print_usrdefcolhead(foutqom);fprintf(foutqom,"h k l  energies[meV]\n");

       }else{
             writeheader(inputpars,foutqom);
            fprintf (foutqom, "#dispersion \n");ini.print_usrdefcolhead(foutqom);fprintf(foutqom,"h k l  energies[meV] > intensities Imag (full calc) [barn/sr/f.u.]   f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n");

            fprintf(foutqei,"#!<--mcphas.mcdisp.qei-->\n");
            writeheader(inputpars,foutqei);
            fprintf (foutqei, "#dispersion displayytext=E(meV)\n#displaylines=false \n");
            ini.print_usrdefcolhead(foutqei);fprintf(foutqei,"h k l Q[A^-1] energy[meV] Imag_dip [barn/sr/f.u.]"
                                                          " Imag [barn/sr/f.u.] Inuc [barn/sr/f.u.]");
            switch(ini.outS)
            {case 0:break;
             case 1:fprintf(foutqei,"Smag_dip: Sxxreal(Q,omega) Sxximag Sxyreal Sxyimag Sxzreal Sxzimag ...Szzimag [barn/sr/f.u.]");break;
             case 2:fprintf(foutqei,"Smag: Sxxreal(Q,omega) Sxximag Sxyreal Sxyimag Sxzreal Sxzimag ...Szzimag [barn/sr/f.u.]");break;
             case 3:fprintf(foutqei,"Smag_dip: Suureal(Q,omega) Suuimag Suvreal Suvimag Suwreal Suwimag ...Swwimag [barn/sr/f.u.]");break;
             case 4:fprintf(foutqei,"Smag: Suureal(Q,omega) Suuimag Suvreal Suvimag Suwreal Suwimag ...Swwimag [barn/sr/f.u.]");break;
            }
            fprintf(foutqei,"f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n");

           fprintf(foutdstot,"#!<--mcphas.mcdisp.dsigma.tot-->\n");writeheader(inputpars,foutdstot);
           fprintf (foutdstot, "#!Total Scattering Cross Section (obtained by DMD method) in energy range [emin=%g ; emax=%g]\n",ini.emin,ini.emax);
           ini.print_usrdefcolhead(foutdstot);
           fprintf(foutdstot,"h k l  dsigma_mag_dip/dOmeg dsigma_mag/dOmeg[barn/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n");

            if (do_Erefine==1){
                       fprintf(foutds,"#!<--mcphas.mcdisp.dsigma-->\n");
                         writeheader(inputpars,foutds);
                         fprintf (foutds, "#Scattering Cross Section \n");
                         ini.print_usrdefcolhead(foutds);
                         fprintf(foutds,"h k l  energy[meV] dsigma/dOmegadE' [barn/mev/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n"); 
                             }  
          }
}


//****************************************************************************************************************************
void writehklblocknumber(FILE * foutqom,FILE * foutqei,FILE * foutdstot,FILE * foutds,
                    FILE * foutqee,FILE * foutqsd,FILE * foutqod,FILE * foutqep,FILE * foutqem,FILE * foutqes,FILE * foutqel,
                    inimcdis & ini,int & calc_rixs,int & do_Erefine,int & counter)
 {  if(ini.hklfile_start_index[0]>0)for(int is=1;is<=ini.hklfile_start_index[0];++is)if(ini.hklfile_start_index[is]==counter)
                       {fprintf(foutqei,"#!hklblock_number=%i\n",is);fprintf(foutqom,"#!hklblock_number=%i\n",is);                                       
                        if(!calc_rixs){fprintf(foutdstot,"#!hklblock_number=%i\n",is);
                                       if (do_Erefine==1){fprintf(foutds,"#!hklblock_number=%i\n",is);}
                                      }
                        if(ini.calculate_chargedensity_oscillation)fprintf(foutqee,"#!hklblock_number=%i\n",is);
                        if(ini.calculate_spindensity_oscillation)fprintf(foutqsd,"#!hklblock_number=%i\n",is);
                        if(ini.calculate_orbmomdensity_oscillation)fprintf(foutqod,"#!hklblock_number=%i\n",is);
                        if(ini.calculate_phonon_oscillation)fprintf(foutqep,"#!hklblock_number=%i\n",is);
                        if(ini.calculate_magmoment_oscillation)fprintf(foutqem,"#!hklblock_number=%i\n",is);
                        if(ini.calculate_spinmoment_oscillation)fprintf(foutqes,"#!hklblock_number=%i\n",is);
                        if(ini.calculate_orbmoment_oscillation)fprintf(foutqel,"#!hklblock_number=%i\n",is);
                        }
 }                    

//****************************************************************************************************************************
// some function to efficiently output standard deviations
void staout(FILE*fout,double & sta,double & sta_int,double & sta_without_antipeaks,double & sta_int_without_antipeaks,double & sta_without_weights,double & sta_int_without_weights,double & sta_without_antipeaks_weights,double & sta_int_without_antipeaks_weights)
    {
    fprintf(fout,"#definitions:\n#\n");
    fprintf(fout,"#sta                               = sum_i |weight(i)|*[Eexp(i) - nearestEcalc]^[2*sign(weight(i))]\n");
    fprintf(fout,"#sta_int                           = sum_i |weight(i)|*[Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^[2*sign(weight(i))]\n",SMALLINT);
    fprintf(fout,"#sta_without_antipeaks             = sum_i_with_weight(i)>0  weight(i)*[Eexp(i) - nearestEcalc]^2\n");
    fprintf(fout,"#sta_int_without_antipeaks         = sum_i_with_weight(i)>0  weight(i)*[Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^2\n",SMALLINT);
    fprintf(fout,"#sta_without_weights               = sum_i [Eexp(i) - nearestEcalc]^[2*sign(weight(i))]\n");
    fprintf(fout,"#sta_int_without_weights           = sum_i [Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^[2*sign(weight(i))]\n",SMALLINT);
    fprintf(fout,"#sta_without_antipeaks_weights     = sum_i_with_weight(i)>0 [Eexp(i) - nearestEcalc]^2\n");
    fprintf(fout,"#sta_int_without_antipeaks_weights = sum_i_with_weight(i)>0 [Eexp(i) - nearestEcalc_with_Int>%gb/srf.u.]^2\n#\n",SMALLINT);
    fprintf(fout,"#!sta= %8.6g \n",sta);
    fprintf(fout,"#!sta_int= %8.6g \n",sta_int);
    fprintf(fout,"#!sta_without_antipeaks= %8.6g \n",sta_without_antipeaks);
    fprintf(fout,"#!sta_int_without_antipeaks= %8.6g \n",sta_int_without_antipeaks);
    fprintf(fout,"#!sta_without_weights= %8.6g \n",sta_without_weights);
    fprintf(fout,"#!sta_int_without_weights= %8.6g \n",sta_int_without_weights);
    fprintf(fout,"#!sta_without_antipeaks_weights= %8.6g \n",sta_without_antipeaks_weights);
    fprintf(fout,"#!sta_int_without_antipeaks_weights= %8.6g \n",sta_int_without_antipeaks_weights);
    }

//****************************************************************************************************
// initialize eigenvector file for observable
//****************************************************************************************************
FILE * evfileinit(const char * filemode,const char*filename,par & inputpars,const char * filetype,int ev_dim)
  {FILE * fout = fopen_errchk (filename,filemode);
          fprintf(fout,"#!<--mcphas.mcdisp.%s-->\n",filetype);
          fprintf (fout, "#!spins_wave_amplitude=1.0\n");
          fprintf (fout, "#!spins_show_ellipses=1.0\n");
          fprintf (fout, "#!spins_show_static_moment_direction=1.0\n");
          fprintf (fout, "#!extended_eigenvector_dimension=%i\n",ev_dim); 
          writeheader(inputpars,fout);
          fprintf (fout, "#!dispersion displayytext=E(meV)\n#Ha[T] Hb[T] Hc[T] T[K] h k l Q[A^-1] energy[meV] int_dipapprFF) [barn/sr/f.u.] int_beyonddipappr [barn/sr/f.u.]  f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
  return fout;
  }

//****************************************************************************************************
// print eigenvector of observable to file
//****************************************************************************************************
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

