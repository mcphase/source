// methods for class inimcdis 
#include "inimcdis.hpp"
#include <martin.h>
#include "../../version"

#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__FreeBSD__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <windows.h>
#endif


int usrdefcols[]={8, 1,2,3,4,5,6,7,8}; // user defined output columns (first number is number of usr def output columns)
                                             // in files mcdisp.qei,qex,qom,dsigma,dsigma.tot
int colcod[]=    {-1,5,6,7,4,12,13,14,8}; // field to store code for assigning type of data to columns of output,
                                           // set default values here (see list below for different types)
                                           // using the out5 out6  ... commands in mcdisp.par these codes can be modified
#define COLHEADDIM 29	
// different output data for columns 10 and 11
const char * colhead []= {  "Qinc[1/A] ", //  0
                            "Qx[1/A]   ",  //   1
                            "Qy[1/A]   ", //    2
                            "Qz[1/A]   ", //    3                                                  
                            "T[K]      ", //    4                                                  
                            "Ha[T]     ", //    5                                                  
                            "Hb[T]     ", //    6                                                  
                            "Hc[T]     ", //    7 
                            "|Q|[1/A]  ",  //    8                                                                  
                            "hprim  ",  //    9                                                                  
                            "kprim  ",  //    10                                                                  
                            "lprim  ",  //    11                                                                 
                            "h ",  //    12                                                                 
                            "k ",  //    13                                                                 
                            "l ",  //    14                                                                 
                            "Hi[T]   ", //    15                                                  
                            "Hj[T]   ", //    16                                                  
                            "Hk[T]   ", //    17 
                            "Ea[kV/mm] ",  //    18      
                            "Eb[kV/mm] ",  //    19      
                            "Ec[kV/mm] ",  //    20      
                            "Ei[kV/mm] ",  //    21      
                            "Ej[kV/mm] ",  //    22      
                            "Ek[kV/mm] ",  //    23      
                            "s1[GPa] ",  //    24     
                            "s2[GPa] ",  //    25      
                            "s3[GPa] ",  //    26      
                            "s4[GPa] ",  //    27      
                            "s5[GPa] ",  //    28      
                            "s6[GPa] "  //    29     
                               };

// different output data for user defined columns ...
double inimcdis::setcolvalue(int & i,Vector & Qvec, double & Qincr,Vector & qprim,Vector & hkl)
{      switch (i) {
case 0:  return Qincr;break;
case 1:  return Qvec(1);break;
case 2:  return Qvec(2);break;
case 3:  return Qvec(3);break;
case 4:  return T;break;
case 5:  return Habc(1);break;
case 6:  return Habc(2);break;
case 7:  return Habc(3);break;
case 8:  return Norm(Qvec);break;
case 9:  return qprim(1);break;
case 10:  return qprim(2);break;
case 11:  return qprim(3);break;
case 12:  return hkl(1);break;
case 13:  return hkl(2);break;
case 14:  return hkl(3);break;
case 15:  return Hext(1);break;
case 16:  return Hext(2);break;
case 17:  return Hext(3);break;
case 18:  return Eabc(1);break;
case 19:  return Eabc(2);break;
case 20:  return Eabc(3);break;
case 21:  return Hext(4);break;
case 22:  return Hext(5);break;
case 23:  return Hext(6);break;
case 24:  return Hext(7);break;
case 25:  return Hext(8);break;
case 26:  return Hext(9);break;
case 27:  return Hext(10);break;
case 28:  return Hext(11);break;
case 29:  return Hext(12);break;
default: fprintf(stderr,"Error mcdisp: unknown column code\n");exit(EXIT_FAILURE);
                    }

return 0;
}


 // *************************************************************************
 // ************************ inimcdis *************************************
 // *************************************************************************
 // class of initial parameters for program mcdisp
 // *************************************************************************


void inimcdis::helpexit() // type info and error exit 
{     printf (" \n %s \n",MCDISPVERSION);
    printf ("use as: mcdisp\n"); 
    printf (" or as: mcdisp [options] [file]\n");
    printf ("  [file] ... input file with mean field set (default mcdisp.mf)\n");
    printf ("Options:\n");
    printf (" -jq           ... calculate J(Q) (Fourier transform of 2ion coupling) store in mcdisp.jq largest evalue and eigenvector\n");
    printf ("                   if energies are given for hkls in mcdisp.par, output file mcdisp_scaled.jq contains scaled parameters\n");
    printf ("                   such that energy of first hkl set corresponds to highest eigenvalue of J(Q)\n");
    printf (" -jqe          ... calculate J(Q) (Fourier transform of 2ion coupling) store in mcdisp.jq all eigenvalues \n");
    printf (" -jqm          ... calculate J(Q) (Fourier transform of 2ion coupling) store in mcdisp.jq all components \n");
    printf (" -cd           ... for calculateing J(Q) (Fourier transform of 2ion coupling) add classical dipole \n");
    printf ("                   interaction using Ewald summation, Bowden J.Phys.C:solid state phys. 14(1981) L827 \n");
    printf ("                   only the first three interaction operators I1.I2.I3 are considered and it \n");
    printf ("                   is assumed that gJ*I1,gJ*I2,gJ*I3 are the components of the magnetic moment (gJ given in sipf file)\n");
    printf (" -max n        ... restrict single ion susceptibility to n lowest\n");
    printf ("                   lying transitions starting from the ground state\n");
    printf (" -minE E       ... an energy range may be given by minE and maxE: only\n");
    printf (" -maxE E           single ion transitions within this energy range will \n");
    printf ("                   be considered\n");
    printf (" -r            ... refine energies\n");
    printf (" -x            ... calculate resonant inelastic x-ray intensities (maximized with respect to azimuth) instead of neutron intensities\n");
    printf (" -xa   stp     ... calculate resonant inelastic x-ray intensities with complete azimuth dependence for each reflection (stp in deg)\n");
    printf (" -xaf  az      ... calculate resonant inelastic x-ray intensities at specified azimuth (deg) for each reflection\n"
            " -X[observable]... to calculate omega and Q dependent  susceptibility tensor X''(Q,omega) for an observable \n"
            "                   observable can be: M (magnetic moment) pel (electricalpolarisation)\n"
            "                   (e.g. use -Xpel for RAMAN and inelastic X-ray scattering intensity on phonons). \n"
            "                      -XM for magnetic susceptibility e.g. to calculate EPR spectra \n"
            );
    printf (" -d            ... calculate intensities in dipole approximation only\n");
    printf (" -v            ... verbose\n");
    printf (" -a            ... do not overwrite output files in results - append results\n");
    printf (" -A            ... do not overwrite output files - compare hkl's to be calculated with existing list in mcdisp.qom and\n");
    printf ("                   continue calculation at last matching q vector\n");
    printf (" -c            ... only create single ion transition file ./results/mcdisp.trs and exit\n");
    printf (" -t            ... read single ion transition file ./results/mcdisp.trs (do not create it)\n");
    printf (" -ninit n      ... maximum number n of (low energy) initial states (single ion transitions)\n");
    printf ("                   (not functional with all single ion modules)\n");
    printf (" -pinit p      ... minimum populationnumber p of initial state (single ion transitions)\n");
    printf ("                   in order to be considered (not functional with all single ion modules)\n");
    printf (" -prefix 001   ... prefix for parameters to be read from mcdisp.par and used for creation of output files\n"
            "                   (useful for running in parallel calculations for different zones: e.g. put in\n"
            "                   mcdisp.par instead of #!hklline= several statements #!001hklline= ... #!002hklline=\n"
            "                   and start several jobs of mcdisp with -prefix 001, -prefix 002 simultaneously, afterwards merge\n"
            "                   output files, e.g. *mcdisp.qei  with appendfile)\n");
    printf (" -ignore_non_hermitian_matrix_error   ... ignores error when energies get complex due to unphysical mf groundstate\n");
    printf ("\n");
    printf ("Note: files which must be in current directory -\n");
    printf ("      ./mcdisp.par, ./mcphas.j, directory ./results\n");
      exit (EXIT_FAILURE);
} 

// print user defined column headers
void inimcdis::print_usrdefcolhead(FILE *fout)
{fprintf(fout,"#");
 for(int i=1;i<=usrdefcols[0];++i)fprintf(fout,"%i%*s",i,(int)strlen(colhead[colcod[usrdefcols[i]]])-1,"");
 for(int i=usrdefcols[0]+1;i<=usrdefcols[0]+2;++i)fprintf(fout,"    %i    ",i);
 for(int i=usrdefcols[0]+3;i<=usrdefcols[0]+5;++i)fprintf(fout,"               %i   ",i);
fprintf(fout,"\n#");
 for(int i=1;i<=usrdefcols[0];++i)fprintf(fout,"%s",colhead[colcod[usrdefcols[i]]]);
}

// print characteristic external parameter string
void inimcdis::mfstring(char *str,size_t t)
{
 snprintf(str,t,"T=%4.4g Hi=%4.4g Hj=%4.4g Hk=%4.4g",
              T,Hext(1),Hext(2),Hext(3));
 if(Norm(Eabc)>SMALL_FIELD)snprintf(str+strlen(str),t-strlen(str)," Ei=%4.4g Ej=%4.4g Ek=%4.4g",
              Hext(4),Hext(5),Hext(6));
  if(Norm(Hext(7,12))>SMALL_FIELD)snprintf(str+strlen(str),t-strlen(str)," s1=%4.4g s2=%4.4g s3=%4.4g s4=%4.4g s5=%4.4g s6=%4.4g",
              Hext(7),Hext(8),Hext(9),Hext(10),Hext(11),Hext(12));
}

// print user defined columns
void inimcdis::print_usrdefcols(FILE *fout,Vector &Qvec, double & Qincr, Vector & qprim,Vector & hkl,bool withtext)
{
 for(int i=1;i<=usrdefcols[0];++i)
 if(withtext)fprintf(fout,"%s=%4.4g ",colhead[colcod[usrdefcols[i]]],myround(setcolvalue(colcod[usrdefcols[i]],Qvec,Qincr,qprim,hkl)));
 else fprintf(fout,"%4.4g ",myround(setcolvalue(colcod[usrdefcols[i]],Qvec,Qincr,qprim,hkl)));
}
// save parameters (which were read from mcdisp.par and mcdisp.mf)
void inimcdis::save()
{char savfilename[MAXNOFCHARINLINE];
  snprintf(savfilename,MAXNOFCHARINLINE,"results/_%s%s",prefix,parfile);
  save(savfilename);
}
// save parameters (which were read from mcdisp.par)
void inimcdis::save(const char * filename)
{ printf("Saving %s\n",filename);
  FILE * fout;int i,j;
  fout=fopen(filename,"w");if (fout==NULL) {fprintf(stderr,"ERROR - file %s cannot be opened \n",filename);exit(EXIT_FAILURE);} 
  fprintf(fout,"# Parameter file  mcdisp.par - read by %s\n",MCDISPVERSION);
  fprintf(fout,"#<!--mcdisp.mcdisp.par>\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
  fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# mcdisp calculates the neutron scattering cross section dsigma/dOmegadE' [barn/sr/meV/f.u.]\n");
  fprintf(fout,"#           f.u.=crystallogrpaphic unit cell (r1xr2xr3) for inelastic and diffuse scattering\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# depending on what is kept constant it follows either kf or ki (1/A)\n");
  fprintf(fout,"# for neutrons: E=(hbar k)^2/2m_n=81.8meV/lambda(A)^2=2.072meV (k(1/A))^2 ...  k(1/A)=sqrt(0.483*E(meV))\n");
  fprintf(fout,"# for X-rays:   E=c hbar k   =1.24e06 meV/lambda(A)=1973202 meV k(1/A)    ...  k(1/A)=5.0679e-7*E(meV)\n");
  if(kf!=0){fprintf(fout,"#!kf=%g\n",kf);}else{fprintf(fout,"#!ki=%g\n",ki);}
  fprintf(fout,"# \n");
  fprintf(fout,"# emin and emax define the energy range in which neutron intensities are calculated\n");
  fprintf(fout,"# for full calculation of the dynamical susceptibility (option \"-r\", inversion of the MF-RPA equation \n");
  fprintf(fout,"# for each point in Q-omega space) the minimum and maximum energy has to be given (energy stepwidth is \n");
  fprintf(fout,"# equal to the parameter epsilon given in the command line after \"-r\")\n");

  fprintf(fout,"#\n");
  fprintf(fout,"#!emin=%g\n",emin);

  fprintf(fout,"#!emax=%g\n",emax);

  fprintf(fout,"#\n");
  fprintf(fout,"# optional switches which can be 0 or 1 are\n");
  fprintf(fout,"#!calculate_magmoment_oscillation=%i  creates mcdisp.qem\n",calculate_magmoment_oscillation);
  fprintf(fout,"#!calculate_spinmoment_oscillation=%i  creates mcdisp.qes\n",calculate_spinmoment_oscillation);
  fprintf(fout,"#!calculate_orbmoment_oscillation=%i  creates mcdisp.qeo\n",calculate_orbmoment_oscillation);
  fprintf(fout,"#!calculate_chargedensity_oscillation=%i  creates mcdisp.qee\n",calculate_chargedensity_oscillation);
  fprintf(fout,"#!calculate_spindensity_oscillation=%i  creates mcdisp.qsd\n",calculate_spindensity_oscillation);
  fprintf(fout,"#!calculate_orbmomdensity_oscillation=%i  creates mcdisp.qod\n",calculate_orbmomdensity_oscillation);
  fprintf(fout,"#!calculate_phonon_oscillation=%i  creates mcdisp.qep\n",calculate_phonon_oscillation);
  fprintf(fout,"#!calculate_pel_oscillation=%i  creates mcdisp.qpe\n",calculate_pel_oscillation);
  fprintf(fout,"#\n" 
               "#     out* controls the type of output in user defined columns in files mcdisp.qei,qex,qom,dsigma,dsigma.tot\n");
  for(int i=1;i<=usrdefcols[0];++i)fprintf(fout,"#!out%i=%i \n",usrdefcols[i],colcod[usrdefcols[i]]);
  fprintf(fout,"#     ... in out*=n the numbers n have the following meaning:\n");
  for(i=0;i<=COLHEADDIM;++i){
  fprintf(fout,"#            %i....%s\n",i,colhead[i]);
                   }
  fprintf(fout,"#\n");
  fprintf(fout,"# optional switch outS for control of the output of the magnetic scattering function in results/mcdisp.qei\n");
  fprintf(fout,"#! outS=%i\n",outS);
  fprintf(fout,"# .. valid values are\n"
               "# 0: not output of Salphabeta\n"
               "# 1: output Sperpalphabeta(Q,omega) in dipole approximation, with alpha,beta=x,y,z\n"
               "# 2: output Sperpalphabeta(Q,omega) going beyond dipole approximation (if possible), with alpha,beta=x,y,z\n"
               "# 3: output Sperpalphabeta(Q,omega) in dipole approximation, with alpha,beta=u,v,w\n"
               "# 4: output Sperpalphabeta(Q,omega) going beyond dipole approximation (if possible), with alpha,beta=u,v,w\n"
               "# 5: output Salphabeta(Q,omega) in dipole approximation, with alpha,beta=x,y,z (no output of dip intensity)\n"
               "# 6: output Salphabeta(Q,omega) in dipole approximation, with alpha,beta=u,v,w (no output of dip intensity)\n"
               "# xyz coordinate refer to y||b, z||(a x b) and x normal to y and z\n"
               "# uvw coordinates refer to u||Q=k-k', w perpendicular to the scattering plane\n"
               "#     (as determined by the cross product of subsequent vectors in the input\n"
               "#     q-vector list) and v perpendicular to u and w, such that uvw form a righthanded system\n#\n#\n");
  fprintf(fout,"# Commands such as the following have been read and used to generate the hkl list below:\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# - a Q vector mesh to be mapped in the calculation it can be in Miller indices\n");
  fprintf(fout,"#hmin=0 hmax=1 deltah=0.1\n");
  fprintf(fout,"#kmin=0 kmax=1 deltak=0.1\n");
  fprintf(fout,"#lmin=0 lmax=1 deltal=0.1\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# or in Qx Qy Qz (1/A) where y||b z||axb and x perpendicular to y and z\n");
  fprintf(fout,"#Qxmin=0 Qxmax=1 deltaQx=0.1\n");
  fprintf(fout,"#Qymin=0 Qymax=1 deltaQy=0.1\n");
  fprintf(fout,"#Qzmin=0 Qzmax=1 deltaQz=0.1\n");
  fprintf(fout,"# - file(s) containing list of Q vectors with (optional) energies of observed excitations to be fitted\n");
  fprintf(fout,"# h k l [E(meV) [statistical_weight  [intensity [fwhm ]]]]\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# hklfile=file1\n");
  fprintf(fout,"# hklfile=file2\n");
  fprintf(fout,"# ...\n#\n");
  fprintf(fout,"# or\n");
  fprintf(fout,"# Qx Qy Qz(1/A) [E(meV) [statistical_weight  [intensity [fwhm ]]]]\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# QxQyQzfile=file1\n");
  fprintf(fout,"# QxQyQzfile=file2\n");
  fprintf(fout,"# ...\n#\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# - some lines in reciprocal space \n");
  fprintf(fout,"#\n");
  fprintf(fout,"#hklline=h1=0 k1=1 l1=0 to hN=1 kN=1 lN=0 Nstp=21\n"); 
  fprintf(fout,"#hklline=h1=0 k1=2 l1=0 to hN=1 kN=1 lN=0 Nstp=21\n"); 
  fprintf(fout,"#\n");
  fprintf(fout,"#QxQyQzline=Qx1=0 Qy1=1 Qz1=0 to QxN=1 QyN=1 QzN=0 Nstp=21\n"); 
  fprintf(fout,"#QxQyQzline=Qx1=0 Qy1=2 Qz1=0 to QxN=1 QyN=1 QzN=0 Nstp=21\n"); 
  fprintf(fout,"#\n");
  fprintf(fout,"# - some planes in reciprocal space \n");
  fprintf(fout,"#\n");
  fprintf(fout,"#hklplane=h0=0 k0=1 l0=0 to hN=1 kN=1 lN=0 Nstp=21 to hM=1 kM=0 lM=3 Mstp=21  \n"); 
  fprintf(fout,"# or\n");
  fprintf(fout,"#QxQyQzplane=Qx0=0 Qy0=1 Qz0=0 to QxN=1 QyN=1 QzN=0 Nstp=21 to QxM=1 QyM=0 QzM=3 Mstp=21  \n"); 
  fprintf(fout,"# or\n");
  fprintf(fout,"# - a list of Q vectors with (optional) energies of observed excitations to be fitted\n");
  fprintf(fout,"# h k l [E(meV) [statistical_weight  [intensity [fwhm ]]]]\n");
     for (j=1;j<=nofhkls;++j) 
  	      {if(hkls[j][0]<=3){for(i=1;i<=hkls[j][0];++i)fprintf(fout,"%g ",hkls[j][i]);fprintf(fout,"\n");} // print hkl
               else             { int k;
               for(k=NOFHKLCOLUMNS;k<=hkls[j][0];k+=NOFHKLCOLUMNS-3){
               for(i=1;i<=3;++i){fprintf(fout,"%g ",hkls[j][i]);} // print hkl
               if(do_jqf){for(i=4;i<=3+hkls[j][k-NOFHKLCOLUMNS+7];++i)fprintf(fout,"%g ",hkls[j][k-NOFHKLCOLUMNS+i]);
                          fprintf(fout,"\n"); 
                         }else{
                                 fprintf(fout,"%g ",hkls[j][k-NOFHKLCOLUMNS+4]);// print E
                                 fprintf(fout,"%g ",hkls[j][k-NOFHKLCOLUMNS+5]);// print weight
                                 for(i=6;i<=NOFHKLCOLUMNS&&hkls[j][k-NOFHKLCOLUMNS+i]>0;++i)fprintf(fout,"%g ",hkls[j][k-NOFHKLCOLUMNS+i]);
               fprintf(fout,"\n");
                                 }
                                                                    }
                                }
	      }
  fprintf(fout,"\n");
  fclose (fout);

  fout=fopen(mf_file,"w");
  fprintf(fout,"# Parameter file  mcdisp.mf - read by %s\n",MCDISPVERSION);
  fprintf(fout,"#<!--mcdisp.mcdisp.mf>\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
  fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"#'T'             temperature T(K)\n");
  fprintf(fout,"#'Ha' 'Hb' 'Hc'  magnetic field(T)\n");
  fprintf(fout,"#'n'             number of atoms in magnetic unit cell\n");
  fprintf(fout,"#'nofatoms'      number of atoms in primitive crystal unit cell\n");
  fprintf(fout,"#'nofcomponents' dimension of moment vector of a magnetic atoms\n");
  fprintf(fout,"#! T=%g  Ha=%g  Hb=%g Hc=%g  n=%i nofatoms=%i nofcomponents=%i\n",T,Hext(1),Hext(2),Hext(3),mf.n(),nofatoms,nofcomponents);
  mf.print(fout);
  fclose (fout);

}


void inimcdis::read_hkl_list(FILE * finhkl,double ** hkls,int readqxqyqz,int & do_jqfile,Vector & abc)
{int i,j;
 float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;  
                 while (feof(finhkl)==0)
                     {if ((i=inputline(finhkl,nn))>=3)
	              {if(readqxqyqz){Vector qijk(1,3),hkl(1,3); // transform to hkl
                                      qijk(1)=nn[1];qijk(2)=nn[2];qijk(3)=nn[3];
                                      ijk2hkl(hkl,qijk,abc);
                                      nn[1]=hkl(1);nn[2]=hkl(2);nn[3]=hkl(3);
                                     }
                       // here check if hkl already in list and if yes, extend its energies
                     if(nofhkls>0&&fabs(hkls[nofhkls][1]-nn[1])+fabs(hkls[nofhkls][2]-nn[2])+fabs(hkls[nofhkls][3]-nn[3])<1e-9)
                       {if(i>3) // only do the energy field extension if energies are given ...
                        {int nold=hkls[nofhkls][0]; // remember old number of hkl
                         hkls[nofhkls+1]=new double [nold+1]; // remember old hkl and energies in appended hkls
                         for(j=0;j<=nold;++j){hkls[nofhkls+1][j]=hkls[nofhkls][j];}
                         delete []hkls[nofhkls]; // free old hkls and claim memory with more storage 
                         hkls[nofhkls]=new double [nold+NOFHKLCOLUMNS-3+1];hkls[nofhkls][0]=nold+NOFHKLCOLUMNS-3;
                         for(j=1;j<=nold;++j){hkls[nofhkls][j]=hkls[nofhkls+1][j];} // fill storage with old values
                         for(j=4;j<=i&&j<=NOFHKLCOLUMNS;++j){hkls[nofhkls][nold+j-3]=nn[j];} // add new values
                         for(j=i+1;j<=NOFHKLCOLUMNS;++j){hkls[nofhkls][nold+j-3]=0.0;} // put weight,int,fwhm ... to zero unless entered
                         if(i==4){hkls[nofhkls][nold+2]=1.0;} // put weight to 1 if not entered
                         if(do_jqfile){if(i>6){hkls[nofhkls][nold+4]=3;}
                                       else{hkls[nofhkls][nold+4]=i-3;}}//store number of J(Q) eigenvalues in fwhm column
                         delete []hkls[nofhkls+1]; // free memory for remembering old hkl and energies
                        }
                       }
                       else 
                       {// a new set of hkl starts
                       ++nofhkls;
	               hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                       hkls[nofhkls][0]=NOFHKLCOLUMNS;if (i==3)hkls[nofhkls][0]=3;
                       for(j=1;j<=i&&j<=NOFHKLCOLUMNS;++j){hkls[nofhkls][j]=nn[j];}
                       for(j=i+1;j<=NOFHKLCOLUMNS;++j){hkls[nofhkls][j]=0.0;} // put weight,int,fwhm ... to zero unless entered
                       if(i==4){hkls[nofhkls][5]=1.0;} // put weight to 1 if not entered
                       if(do_jqfile){if(i>6){hkls[nofhkls][7]=3;}
                                     else{hkls[nofhkls][7]=i-3;}}//store number of J(Q) eigenvalues in fwhm column
                       }
	              }
                     }
}

// findnewmatch = true: check if new match (which is not listed in the n prefixes of lofpref) can be found in instr using
//                      char prefix (which may contain a wildcard '*'), if it can be found increase n
//                      and put new match prefix in lofpref and put findnewmatch to false and return
//                      extract_with_prefix for the new match
//                      if it cannot be found return extract(instr,parameter)
// findnewmatch = false: return extract_with_prefix for the prefix=lofpref[n], i.e. prefix+var has to be found to return 0
int inimcdis::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,int & var,parser & ob)
{int ret; size_t s=MAXNOFCHARINLINE;char val[MAXNOFCHARINLINE];
ret=extract_match(findnewmatch,  n,lofpref ,instr,pref,  parameter, val,s,1);
if(ret==0){var=(int)ob.eval_exp(val);}
return ret;
}
int inimcdis::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,float & var,parser & ob)
{int ret; size_t s=MAXNOFCHARINLINE;char val[MAXNOFCHARINLINE];
ret=extract_match(findnewmatch,  n,lofpref ,instr,pref,  parameter, val,s,1);
if(ret==0)var=(float)ob.eval_exp(val);
return ret;
}
int inimcdis::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,double & var,parser & ob)
{int ret; size_t s=MAXNOFCHARINLINE;char val[MAXNOFCHARINLINE];
ret=extract_match(findnewmatch,  n,lofpref ,instr,pref,  parameter, val,s,1);
if(ret==0)var=ob.eval_exp(val);
return ret;
}
int inimcdis::extract_match(bool & findnewmatch, int & n,char**lofpref ,char * instr,char * pref, const char * parameter,char * var,size_t ns,int  m)
{
 if(findnewmatch==true)
 {if(pref[0]!='\0')
  {char prefvar[MAXNOFCHARINLINE];char * s,*p;
  snprintf(prefvar,MAXNOFCHARINLINE,"%s%s",pref,parameter);
  s=wstrstr(instr,prefvar);
  if(s!=NULL)
  { // ok there seems to be a new match - see if it is already in the list
   snprintf(prefvar,MAXNOFCHARINLINE,"%s",parameter);
   p=wstrstr(instr,prefvar);
   int i=0;for(char * t=s;t<p;++t){prefvar[i]=*t;++i;}prefvar[i]='\0'; // put the prefix to prefvar
    i=0;bool mm=false;while(i<n&&mm==false){mm=match(prefvar,lofpref[i]);
    ++i;}
    if(i==n&&mm==false){lofpref[i]=new char [strlen(prefvar)+2];snprintf(lofpref[i],strlen(prefvar)+1,"%s",prefvar);
                          snprintf(prefix,MAXNOFCHARINLINE,"%s",prefvar);
                     ++n; findnewmatch=false;printf("new matching prefix found: %s\n",lofpref[i]);
            }
   }// no new match -> extract without prefix
   else
   {
    return extract(instr,parameter,var,ns, m);
   }
  }
 else
 {if(n==0){int i=0;lofpref[i]=new char [1];lofpref[i][0]='\0';
                          prefix[0]='\0';
                     ++n; findnewmatch=false;}
 }
 }
 if(n>0)
 {
 if(findnewmatch==true) return extract(instr,parameter,var,ns,m);  // if we still have to find new match: return without prefix 
  else return extract_with_prefix(instr,lofpref[n-1],parameter,var,ns,m); // if a new match has been found: return parameters with this new prefix
 }
 else
 {return extract_with_prefix(instr,pref,parameter,var,ns, m);
 }
}

// *************************************************************************
int inimcdis::load (int & nofinis,char**lofpref,char * mffile,char * pref,int & do_jqfile,Vector & abc,int & nofcomp,int & nofat,int  verbose)
{ bool findnewmatch=true;if(nofinis==-1){findnewmatch=false;}
   errno=1;do_jqf=do_jqfile;nofthreads=0;outcolset=false;
  char instr[MAXNOFCHARINLINE],hklfile[MAXNOFCHARINLINE],hklline[MAXNOFCHARINLINE],somestring[MAXNOFCHARINLINE];
  int nofhkllists=1; parser ob;
  Hext=0; Habc=0;Eabc=0;
  FILE *fin,*finhkl;float N,M,h0,k0,l0,h1,k1,l1,hN,kN,lN,hM,kM,lM;
   strcpy(prefix,pref); // set prefix
 // ****************************** read mf configuration from mffile *****************************************  
 fin=fopen(mffile,"rb");
   if (fin==NULL) {fprintf(stderr,"#Warning - file %s not found - trying to read mcdisp.mf\n",mffile);
   snprintf(mffile,MAXNOFCHARINLINE,"mcdisp.mf");fin=fopen(mffile,"rb");
    if (fin==NULL) {fprintf(stderr,"Warning - file %s not found - doing calculation at T=300K assuming zero mean and external fields / stress\n",mffile);
    }
   }
 snprintf(mf_file,MAXNOFCHARINLINE,"results/_%s",mffile); // store filename of mffile used
 if(fin==NULL){snprintf(instr,MAXNOFCHARINLINE,"#!T=300 Ha=0 Hb=0 Hc=0 n=1 spins nofatoms=1 in primitive basis nofcomponents=%i - configuration",nofcomponents);
 T=300; Hext=0;nofatoms=nofat;nofcomponents=nofcomp;
}else {instr[0]='#';  instr[1]='\0';
  while(instr[strspn(instr," \t")]=='#'&&instr[strspn(instr," \t#")]!='!'){fgets(instr,MAXNOFCHARINLINE,fin);}
  parseline(instr,ob);
  extract(instr,"T",T,ob); 
  extract(instr,"Ha",Habc[1],ob); 
  extract(instr,"Hb",Habc[2],ob);
  extract(instr,"Hc",Habc[3],ob); 
  extract(instr,"Ea",Eabc[1],ob); 
  extract(instr,"Eb",Eabc[2],ob);
  extract(instr,"Ec",Eabc[3],ob); 
  extract(instr,"Hi",Hext[1],ob); 
  extract(instr,"Hj",Hext[2],ob);
  extract(instr,"Hk",Hext[3],ob); 
  extract(instr,"Ei",Hext[4],ob); 
  extract(instr,"Ej",Hext[5],ob);
  extract(instr,"Ek",Hext[6],ob); 
  extract(instr,"s1",Hext[7],ob); 
  extract(instr,"s2",Hext[8],ob);
  extract(instr,"s3",Hext[9],ob); 
  extract(instr,"s4",Hext[10],ob); 
  extract(instr,"s5",Hext[11],ob);
  extract(instr,"s6",Hext[12],ob); 
  
  crosscheck_H_E(Hext,Habc,Eabc,abc); 

  if(verbose!=0)printf("# reading mean field configuration from file %s mf=gj muB heff [meV]\n#%s \n",mffile,instr);
  
  
  
  extract(instr,"nofatoms",nofatoms,ob); 
  extract(instr,"nofcomponents",nofcomponents,ob); 
  if(nofcomponents!=nofcomp){fprintf(stderr,"ERROR loading mean field configuration nofcomponents from mcphas.j (%i) different from file %s (%i)\n",nofcomp,mffile,nofcomponents);exit(EXIT_FAILURE);}
  if(nofatoms!=nofat){fprintf(stderr,"ERROR loading mean field configuration nofatoms from mcphas.j (%i) different from file %s (%i)\n",nofat,mffile,nofatoms);exit(EXIT_FAILURE);}
  if(mf.load(fin)==0)
   {fprintf(stderr,"ERROR loading mean field configuration\n");exit(EXIT_FAILURE);}
  fclose(fin);
 }

 strcpy(info,instr);

 //********************************  
    errno = 0; 
  
  // **************** initialize parameters to default values ********************************************
  emin=-DBL_MAX;emax=DBL_MAX;
  ki=0;kf=0;
  calculate_magmoment_oscillation=0;
  calculate_spinmoment_oscillation=0;
  calculate_orbmoment_oscillation=0;
  calculate_chargedensity_oscillation=0;
  calculate_spindensity_oscillation=0;
  calculate_orbmomdensity_oscillation=0;
  calculate_phonon_oscillation=0;
  calculate_pel_oscillation=0;
  outS=0;
  qmin=0;qmax=0;deltaq=0;
 // ******************************** reading parameters  from mcdisp.par run 1 to determine i=nofhkls-estimate****************************************************
  int i=0,hklblock=0,QxQyQzblock=0,j;
  if(verbose)printf("reading file %s\n",parfile);
  fin = fopen(parfile, "rb"); 
if (fin==NULL) { return 1;
}else
{   
  while (fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
  {++i; // i is used to estimate an upper boundary for the number of hkls in the hkl list 
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"emin",emin,ob); 
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"emax",emax,ob); 
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"ki",ki,ob); 
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kf",kf,ob); 
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_magmoment_oscillation",calculate_magmoment_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_spinmoment_oscillation",calculate_spinmoment_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_orbmoment_oscillation",calculate_orbmoment_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_chargedensity_oscillation",calculate_chargedensity_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_spindensity_oscillation",calculate_spindensity_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_orbmomdensity_oscillation",calculate_orbmomdensity_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_phonon_oscillation",calculate_phonon_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"calculate_pel_oscillation",calculate_pel_oscillation,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"outS",outS,ob);
     extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"nofthreads",nofthreads,ob);
     for(int j=1;j<=usrdefcols[0];++j) // extract user defined output columns
     {snprintf(somestring,MAXNOFCHARINLINE,"out%i",usrdefcols[j]);
      if(0==extract_match( findnewmatch,nofinis,lofpref,instr,prefix, somestring,colcod[usrdefcols[j]],ob))outcolset=true;
    }

     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hmin",qmin[1],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kmin",qmin[2],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"lmin",qmin[3],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hmax",qmax[1],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kmax",qmax[2],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"lmax",qmax[3],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltah",deltaq[1],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltak",deltaq[2],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltal",deltaq[3],ob); 
     if(hklblock==9){++nofhkllists;hklblock=0;i+=(int)ceil(fabs((qmax(1)-qmin(1))/deltaq(1)+1)*fabs((qmax(2)-qmin(2))/deltaq(2)+1)*fabs((qmax(3)-qmin(3))/deltaq(3)+1));}

     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qxmin",qmin[1],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qymin",qmin[2],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qzmin",qmin[3],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qxmax",qmax[1],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qymax",qmax[2],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qzmax",qmax[3],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltaQx",deltaq[1],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltaQy",deltaq[2],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltaQz",deltaq[3],ob); 
     if(QxQyQzblock==9){++nofhkllists;QxQyQzblock=0;i+=(int)ceil(fabs((qmax(1)-qmin(1))/deltaq(1)+1)*fabs((qmax(2)-qmin(2))/deltaq(2)+1)*fabs((qmax(3)-qmin(3))/deltaq(3)+1));}

     if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hklfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");while (fgets(hklfile,MAXNOFCHARINLINE,finhkl)!=NULL)++i;
                  fclose(finhkl);++nofhkllists;
                 }
     if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"QxQyQzfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");while (fgets(hklfile,MAXNOFCHARINLINE,finhkl)!=NULL)++i;
                  fclose(finhkl);++nofhkllists;
                 }
     if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hklline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) Nsteps=21
                 {if(!extract(instr,"Nstp",N))i+=N;++nofhkllists;
                 }
     if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"QxQyQzline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) Nsteps=21
                 {if(!extract(instr,"Nstp",N))i+=N;++nofhkllists;
                 }
     if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hklplane",hklline,MAXNOFCHARINLINE-1,1))  
                 {if(!extract(instr,"Nstp",N)&&!extract(instr,"Mstp",M))i+=N*M;++nofhkllists;
                 }
     if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"QxQyQzplane",hklline,MAXNOFCHARINLINE-1,1))  
                 {if(!extract(instr,"Nstp",N)&&!extract(instr,"Mstp",M))i+=N*M;++nofhkllists;
                 }
  }
  fclose (fin);
 }
 // ************************************ end reading parameters first run to determine i (estimate for nofhkls) *****************************************
 // check parameters
  if (ki==0) {if (kf==0) kf=100;
              if(verbose)fprintf(stdout,"#Calculating intensities for  kf=const=%4.4g/A\n",kf);
	     }
	     else
	     {kf=0;
	      if(verbose)fprintf(stdout,"#Calculating intensities for ki=const=%4.4g/A\n",ki);
	     }
  // Checks if nofthreads set in mcdisp.par, if not check environment or use system calls
  if(nofthreads<1) {
    char* c_nofthreads=getenv("MCPHASE_NOFTHREADS");  // Check if system environment variable set from dos.bat/lin.bat
    if (c_nofthreads)
       nofthreads = atoi(c_nofthreads);
    else {  
#if defined(__linux__)                               // System-dependent calls to find number of processors (from GotoBLAS)
       nofthreads = get_nprocs();
#elif defined(__FreeBSD__) || defined(__APPLE__)
       int m[2]; size_t len;
       m[0] = CTL_HW; m[1] = HW_NCPU; len = sizeof(int);
       sysctl(m, 2, &nofthreads, &len, NULL, 0);
#else
       SYSTEM_INFO sysinfo; GetSystemInfo(&sysinfo);
       nofthreads = sysinfo.dwNumberOfProcessors;
#endif
    }
    if(nofthreads<1||nofthreads>255) nofthreads=1;             // All else fails: use only 1 thread
  }
  if(verbose)printf("# nofthreads=%i\n",nofthreads);
  // ***************************************** reread mcdisp.par creating the hkl list ******************************************************************
  if(hkls!=NULL){ for (int ii=1;ii<=nofhkls;++ii) { delete []hkls[ii];}delete []hkls;}
  hkls=new double *[i+10]; // dimension the list
  //printf("nofhkl estimate i=%i\n",i);
  nofhkls=0;hklblock=0;QxQyQzblock=0;
  if(hklfile_start_index!=NULL) delete []hklfile_start_index;
  hklfile_start_index= new int [nofhkllists+1];
  hklfile_start_index[0]=nofhkllists;
  nofhkllists=0;Vector hkl(1,3),qijk(1,3);

       // now read the hkls in mcdisp.par
      ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
      fin = fopen(parfile, "rb");read_hkl_list(fin,hkls,0,do_jqfile,abc); fclose(fin); 
      if(nofhkls==0){--nofhkllists;}

  fin = fopen(parfile, "rb"); // if in mcdisp.par we find a hklfile= ... insert hkl from this file into list
            while (fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
               {// treat hklblocks
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hmin",qmin[1],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kmin",qmin[2],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"lmin",qmin[3],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hmax",qmax[1],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"kmax",qmax[2],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"lmax",qmax[3],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltah",deltaq[1],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltak",deltaq[2],ob); 
     hklblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltal",deltaq[3],ob); 
     if(hklblock==9){++nofhkllists;hklblock=0;hklfile_start_index[nofhkllists]=nofhkls+1;
                    printf("# ... hklblock hklmin(%g %g %g) to hklmax(%g %g %g) with hklstepsize (%g %g %g)\n",qmin(1),qmin(2),qmin(3),qmax(1),qmax(2),qmax(3),deltaq(1),deltaq(2),deltaq(3));
                   for(h1=qmin(1);h1<=qmax(1);h1+=deltaq(1))
                   for(k1=qmin(2);k1<=qmax(2);k1+=deltaq(2))
                   for(l1=qmin(3);l1<=qmax(3);l1+=deltaq(3))
                                    {++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=h1;
                                     hkls[nofhkls][2]=k1;
                                     hkls[nofhkls][3]=l1;
                                    }
                    }

     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qxmin",qmin[1],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qymin",qmin[2],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qzmin",qmin[3],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qxmax",qmax[1],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qymax",qmax[2],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"Qzmax",qmax[3],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltaQx",deltaq[1],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltaQy",deltaq[2],ob); 
     QxQyQzblock+=1-extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"deltaQz",deltaq[3],ob); 
     if(QxQyQzblock==9){++nofhkllists;hklblock=0;hklfile_start_index[nofhkllists]=nofhkls+1;
                    printf("# ... hklblock hklmin(%g %g %g) to hklmax(%g %g %g) with hklstepsize (%g %g %g)\n",qmin(1),qmin(2),qmin(3),qmax(1),qmax(2),qmax(3),deltaq(1),deltaq(2),deltaq(3));
                   for(qijk(1)=qmin(1);qijk(1)<=qmax(1);qijk(1)+=deltaq(1))
                   for(qijk(2)=qmin(2);qijk(2)<=qmax(2);qijk(2)+=deltaq(2))
                   for(qijk(3)=qmin(3);qijk(3)<=qmax(3);qijk(3)+=deltaq(3))
                                    {++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     ijk2hkl(hkl, qijk,abc);
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=hkl(1);
                                     hkls[nofhkls][2]=hkl(2);
                                     hkls[nofhkls][3]=hkl(3);
                                    }
                    }
                // treat hklplane statements
                if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hklplane",hklline,MAXNOFCHARINLINE-1,1)) //#hklplane=h0=0 k0=1 l0=0 to hN=1 kN=1 lN=0 Nstp=21 to hM=1 kM=0 lM=3 Mstp=21  
                 { if(extract(instr,"h0",h0)){printf("error mcdisp reading %s: in hklplane - h0 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"k0",k0)){printf("error mcdisp reading %s: in hklplane - k0 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"l0",l0)){printf("error mcdisp reading %s: in hklplane - l0 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"hN",hN)){printf("error mcdisp reading %s: in hklplane - hN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"kN",kN)){printf("error mcdisp reading %s: in hklplane - kN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"lN",lN)){printf("error mcdisp reading %s: in hklplane - lN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"hM",hM)){printf("error mcdisp reading %s: in hklplane - hM not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"kM",kM)){printf("error mcdisp reading %s: in hklplane - kM not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"lM",lM)){printf("error mcdisp reading %s: in hklplane - lM not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading %s: in hklplane - Nstp  not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Mstp",M)){printf("error mcdisp reading %s: in hklplane - Mstp  not found",parfile);exit (EXIT_FAILURE);}
                   printf("# ... hklplane (%g %g %g) to (%g %g %g) with %g points to (%g %g %g) with %g points\n",h0,k0,l0,hN,kN,lN,N,hM,kM,lM,M);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(int ii=1;ii<=N;++ii)for(j=1;j<=M;++j){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=h0+(hN-h0)*(ii-1)/(N-1)+(hM-h0)*(j-1)/(M-1);
                                     hkls[nofhkls][2]=k0+(kN-k0)*(ii-1)/(N-1)+(kM-k0)*(j-1)/(M-1);
                                     hkls[nofhkls][3]=l0+(lN-l0)*(ii-1)/(N-1)+(lM-l0)*(j-1)/(M-1);
                                    }
                 }

                // treat QxQyQzplane statements
                if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"QxQyQzplane",hklline,MAXNOFCHARINLINE-1,1)) //#QxyQzplane=Qx0=0 Qy0=1 Qz0=0 to QxN=1 QyN=1 QzN=0 Nstp=21 to QxM=1 QyM=0 QzM=3 Mstp=21
                 { if(extract(instr,"Qx0",h0)){printf("error mcdisp reading %s: in QxyQzplane - Qx0 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Qy0",k0)){printf("error mcdisp reading %s: in QxyQzplane - Qy0 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Qz0",l0)){printf("error mcdisp reading %s: in QxyQzplane - Qz0 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QxN",hN)){printf("error mcdisp reading %s: in QxyQzplane - QxN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QyN",kN)){printf("error mcdisp reading %s: in QxyQzplane - QyN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QzN",lN)){printf("error mcdisp reading %s: in QxyQzplane - QzN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QxM",hM)){printf("error mcdisp reading %s: in QxyQzplane - QxM not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QyM",kM)){printf("error mcdisp reading %s: in QxyQzplane - QyM not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QzM",lM)){printf("error mcdisp reading %s: in QxyQzplane - QzM not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading %s: in QxyQzplane - Nstp  not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Mstp",M)){printf("error mcdisp reading %s: in QxyQzplane - Mstp  not found",parfile);exit (EXIT_FAILURE);}
                   printf("# ... QxyQzplane (%g %g %g) to (%g %g %g) with %g points to (%g %g %g) with %g points\n",h0,k0,l0,hN,kN,lN,N,hM,kM,lM,M);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(int ii=1;ii<=N;++ii)for(j=1;j<=M;++j){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     qijk(1)=h0+(hN-h0)*(ii-1)/(N-1)+(hM-h0)*(j-1)/(M-1);
                                     qijk(2)=k0+(kN-k0)*(ii-1)/(N-1)+(kM-k0)*(j-1)/(M-1);
                                     qijk(3)=l0+(lN-l0)*(ii-1)/(N-1)+(lM-l0)*(j-1)/(M-1);
                                     ijk2hkl(hkl, qijk,abc);
                                     hkls[nofhkls][1]=hkl(1);
                                     hkls[nofhkls][2]=hkl(2);
                                     hkls[nofhkls][3]=hkl(3);
                                    }
                 }


                // treat hklline statements
                if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hklline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) N=21
                 { if(extract(instr,"h1",h1)){printf("error mcdisp reading %s: in hklline - h1 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"k1",k1)){printf("error mcdisp reading %s: in hklline - k1 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"l1",l1)){printf("error mcdisp reading %s: in hklline - l1 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"hN",hN)){printf("error mcdisp reading %s: in hklline - hN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"kN",kN)){printf("error mcdisp reading %s: in hklline - kN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"lN",lN)){printf("error mcdisp reading %s: in hklline - lN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading %s: in hklline - Nstp  not found",parfile);exit (EXIT_FAILURE);}
                    printf("# ... hklline (%g %g %g) to (%g %g %g) with %g points\n",h1,k1,l1,hN,kN,lN,N);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(int ii=1;ii<=N;++ii){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=h1+(hN-h1)*(ii-1)/(N-1);
                                     hkls[nofhkls][2]=k1+(kN-k1)*(ii-1)/(N-1);
                                     hkls[nofhkls][3]=l1+(lN-l1)*(ii-1)/(N-1);
                                    }
                 }
                // treat QxQyQzline statements
                if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"QxQyQzline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) N=21
                 { if(extract(instr,"Qx1",h1)){printf("error mcdisp reading %s: in hklline - h1 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Qy1",k1)){printf("error mcdisp reading %s: in hklline - k1 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Qz1",l1)){printf("error mcdisp reading %s: in hklline - l1 not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QxN",hN)){printf("error mcdisp reading %s: in hklline - hN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QyN",kN)){printf("error mcdisp reading %s: in hklline - kN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"QzN",lN)){printf("error mcdisp reading %s: in hklline - lN not found",parfile);exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading %s: in hklline - N  not found",parfile);exit (EXIT_FAILURE);}
                    printf("# ... QxQyQzline (%g %g %g)/A to (%g %g %g)/A with %g points\n",h1,k1,l1,hN,kN,lN,N);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(int ii=1;ii<=N;++ii){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     qijk(1)=h1+(hN-h1)*(ii-1)/(N-1);
                                     qijk(2)=k1+(kN-k1)*(ii-1)/(N-1);
                                     qijk(3)=l1+(lN-l1)*(ii-1)/(N-1);
                                     ijk2hkl(hkl, qijk,abc);
                                     hkls[nofhkls][1]=hkl(1);
                                     hkls[nofhkls][2]=hkl(2);
                                     hkls[nofhkls][3]=hkl(3);
                                    }
                 }
                 // treat hklfile statements
                if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"hklfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                  read_hkl_list(finhkl,hkls,0,do_jqfile,abc);
                  fclose(finhkl);
                 }
                 // treat QxQyQzfile statements
                if(!extract_match( findnewmatch,nofinis,lofpref,instr,prefix,"QxQyQzfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                  read_hkl_list(finhkl,hkls,1,do_jqfile,abc);
                  fclose(finhkl);
                 }

              }
       fclose (fin);
       // now read also the hkls in mcdisp.par
   //   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
   //   fin = fopen(parfile, "rb");read_hkl_list(fin,hkls,0,do_jqfile,abc); fclose(fin); 
      if(nofhkls==0){++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                nofhkls=1;hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=1;
                                     hkls[nofhkls][2]=0;
                                     hkls[nofhkls][3]=0;} 
  if(i+9<nofhkls){fprintf(stderr,"Error mcdisp: nofhkl =%i > estimate %i\n",nofhkls,i);exit(EXIT_FAILURE);}
     
return 0;
}


// *************************************************************************
//load parameters from file
int inimcdis::load (char * mffile,char * pref,int & do_jqfile,Vector & abc,int & nofcomp,int & nofat)
{int n=-1;char ** lp; lp=NULL;
return load(n,lp,mffile,pref,do_jqfile, abc,nofcomp,nofat,0);
}


//constructor ... load initial parameters from file
inimcdis::inimcdis(const char * file,char * pref,char * mffile,
                   int & do_jqfile,bool inc_cd,Vector & abc,
                   int & nofcomp,int & nofat)
{hkls=NULL;hklfile_start_index=NULL;Hext=Vector(1,HEXT_DIMENSION);Habc=Vector(1,3);Eabc=Vector(1,3);
 qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);mf=mfcf(1,1,1,nofat,nofcomp);
  include_cd=inc_cd;
  parfile= new char [MAXNOFCHARINLINE]; 
  mf_file= new char [MAXNOFCHARINLINE]; 
  snprintf(parfile,MAXNOFCHARINLINE,"%s%s",pref,file);
 info= new char [MAXNOFCHARINLINE];
  prefix = new char[MAXNOFCHARINLINE];
  strcpy(prefix,pref);
  if(load(mffile,pref,do_jqfile, abc,nofcomp,nofat)!=0){if(pref[0]!='\0'){fprintf(stderr,"File %s not found - trying %s\n",parfile,file);
                strcpy(parfile,file);}
                if(load(mffile,pref,do_jqfile, abc,nofcomp,nofat)!=0){
    fprintf(stderr,"# Warning: Cannot load file %s - using default values ! \n",parfile); 
// insert  here default values and save into mcdisp.par
emin=-100;emax=100;ki=0;kf=100;colcod[1]=5;colcod[2]=6;colcod[3]=7;colcod[4]=4;
nofhkls=0;
info[0]='\0';prefix[0]='\0';
save(parfile); 
  }
 }
}

//kopier-konstruktor 
inimcdis::inimcdis (const inimcdis & p)
{do_jqf=p.do_jqf;
 include_cd=p.include_cd;
 parfile= new char [MAXNOFCHARINLINE];
  strcpy(parfile,p.parfile);
 mf_file= new char [MAXNOFCHARINLINE];
  strcpy(mf_file,p.mf_file);
 info= new char [MAXNOFCHARINLINE];strcpy(info,p.info);
 prefix= new char [MAXNOFCHARINLINE]; strcpy(prefix,p.prefix);  
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  Eabc=Vector(1,3);Habc=Vector(1,3);
  Hext=Vector(1,HEXT_DIMENSION); Hext=p.Hext;
  Eabc=p.Eabc;Habc=p.Habc;
  qmin=p.qmin;
  qmax=p.qmax;
  emin=p.emin;
  emax=p.emax;
  kf=p.kf;
  ki=p.ki;
  calculate_magmoment_oscillation=p.calculate_magmoment_oscillation;
  calculate_spinmoment_oscillation=p.calculate_spinmoment_oscillation;
  calculate_orbmoment_oscillation=p.calculate_orbmoment_oscillation;
  calculate_chargedensity_oscillation=p.calculate_chargedensity_oscillation;
  calculate_spindensity_oscillation=p.calculate_spindensity_oscillation;
  calculate_orbmomdensity_oscillation=p.calculate_orbmomdensity_oscillation;
  calculate_phonon_oscillation=p.calculate_phonon_oscillation;
  calculate_pel_oscillation=p.calculate_pel_oscillation;
  outS=p.outS;
    deltaq=p.deltaq;  
  nofthreads=p.nofthreads;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  nofhkls=p.nofhkls;
  int i,j;
   if(p.hkls!=NULL)
    {  hkls=new double *[nofhkls+10];
      for (j=1;j<=nofhkls;++j) 
  	      {
   if ((int)p.hkls[j][0]==3){hkls[j]=new double [NOFHKLCOLUMNS+1];}
               else {hkls[j]=new double [(int)p.hkls[j][0]+1];}
               for(i=0;i<=p.hkls[j][0];++i)
         	    {hkls[j][i]=p.hkls[j][i];}
	      }
     }else hkls=NULL;
    if(p.hklfile_start_index!=NULL)
    {  int nofhkllists=p.hklfile_start_index[0];
       hklfile_start_index= new int [nofhkllists+1];hklfile_start_index[0]=nofhkllists;
      for (j=1;j<=nofhkllists;++j) hklfile_start_index[j]=p.hklfile_start_index[j]; 
    } else hklfile_start_index=NULL;
  mf=mfcf(1,1,1,nofatoms,nofcomponents);mf=p.mf;T=p.T;
}

//destruktor
inimcdis::~inimcdis ()
{if (parfile!=NULL)delete []parfile;
 if (mf_file!=NULL)delete []mf_file;
 if (info!=NULL)delete []info;
 if(prefix!=NULL) delete []prefix;
 int i;
//if (nofhkls==1)
if(hkls!=NULL)
 { for (i=1;i<=nofhkls;++i) 
   { delete []hkls[i];}
   delete []hkls;
   delete  []hklfile_start_index;
 }
}


//***************************************************************
//constructor ... load initial parameters from file
inimdpars::inimdpars (const char * file,char * pref,char * mffile,
             int & do_jqfile,bool inc_cd,Vector & abc,
             int & nofcomponents,int & nofatoms, int & verbose)
{ inis=new inimcdis*[MAXNOFINIS];
  char * lofprefixes[MAXNOFINIS];
  int nofinisold=-1;nofinis=0;
// here we have to load inis[1...nofinis] with different prefixes matching pref - until no new matching
// prefix is found ...
  while(nofinisold<nofinis&&nofinis<MAXNOFINIS)
  {inis[nofinis]=new inimcdis(file,pref,mffile,do_jqfile,inc_cd,abc,nofcomponents,nofatoms);
   nofinisold=nofinis;(*inis[nofinis]).load(nofinis,lofprefixes,mffile,pref,do_jqfile,abc,nofcomponents,nofatoms,verbose);
  }
 
// remove last inis, because it does not contain a new prefix
if(nofinis>0&&nofinis<MAXNOFINIS){delete inis[nofinis];} 
if(nofinis==0)nofinis=1;

}


//kopier-konstruktor  inimcdiss
inimdpars::inimdpars (const inimdpars & p)
{ nofinis=p.nofinis;
  inis=new inimcdis*[MAXNOFINIS];
  for(int i=0;i<nofinis;++i)inis[i]=new inimcdis((*p.inis[i]));
}

//destruktor inimcdiss
inimdpars::~inimdpars ()
{//printf("hello destruktor inimcdiss %i\n",nofinis);  
 for(int i=0;i<nofinis;++i)delete  inis[i];
delete []inis;
//printf("hello destruktor inimcdis\n");  
 }
