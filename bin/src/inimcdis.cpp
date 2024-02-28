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


int usrdefcols[]={4, 1,2,3,4}; // user defined output columns (first number is number of usr def output columns)
                                             // in files mcdisp.qei,qex,qom,dsigma,dsigma.tot
int colcod[]=    {-1,5,6,7,4}; // field to store code for assigning type of data to columns of output,
                                           // set default values here (see list below for different types)
                                           // using the out11 command in mcdiff.in these codes can be modified
#define COLHEADDIM 11
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
                            "lprim  "  //    11                                                                 
                           };

// different output data for user defined columns ...
double inimcdis::setcolvalue(int i,Vector & Qvec, double & Qincr,Vector & qprim)
{
         switch (i) {
case 0:  return Qincr;break;
case 1:  return Qvec(1);break;
case 2:  return Qvec(2);break;
case 3:  return Qvec(3);break;
case 4:  return T;break;
case 5:  return Hext(1);break;
case 6:  return Hext(2);break;
case 7:  return Hext(3);break;
case 8:  return Norm(Qvec);break;
case 9:  return qprim(1);break;
case 10:  return qprim(2);break;
case 11:  return qprim(3);break;
default: fprintf(stderr,"Error mcdisp: unknown column code\n");exit(EXIT_FAILURE);
                    }

return 0;
}


 // *************************************************************************
 // ************************ inipar *************************************
 // *************************************************************************
 // class of initial parameters for program mcphas
 // *************************************************************************

// print user defined column headers
void inimcdis::print_usrdefcolhead(FILE *fout)
{fprintf(fout,"#");
 for(int i=1;i<=usrdefcols[0];++i)fprintf(fout,"%s",colhead[colcod[usrdefcols[i]]]);
}

// print user defined column headers
void inimcdis::print_usrdefcols(FILE *fout,Vector &Qvec, double & Qincr, Vector & qprim)
{
 for(int i=1;i<=usrdefcols[0];++i)fprintf(fout,"%4.4g ",myround(setcolvalue(colcod[usrdefcols[i]],Qvec,Qincr,qprim)));
}
// save parameters (which were read from mcdisp.par)
void inimcdis::save()
{save(savfilename);
}
// save parameters (which were read from mcdisp.par)
void inimcdis::save(const char * filename)
{  FILE * fout;int i,j;
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
  fout=fopen("results/_mcdisp.mf","w");
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


void inimcdis::read_hkl_list(FILE * finhkl,double ** hkls,int readqxqyqz,int do_jqfile,Vector & abc)
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
// *************************************************************************
//constructor ... load initial parameters from file
inimcdis::inimcdis (const char * file,char * spinfile,char * pref,int do_jqfile,Vector & abc)
{ errno=1;do_jqf=do_jqfile;
  char instr[MAXNOFCHARINLINE],hklfile[MAXNOFCHARINLINE],hklline[MAXNOFCHARINLINE],somestring[MAXNOFCHARINLINE];
  int nofhkllists=1;Hext=Vector(1,3);
  FILE *fin,*finhkl;float N,M,h0,k0,l0,h1,k1,l1,hN,kN,lN,hM,kM,lM;
  prefix= new char [strlen(pref)+1]; strcpy(prefix,pref); // set prefix
 // ****************************** read mf configuration from spinfile *****************************************  
  fin=fopen(spinfile,"rb");
   if (fin==NULL) {fprintf(stderr,"#Warning - file %s not found - trying to read mcdisp.mf\n",spinfile);
   snprintf(spinfile,MAXNOFCHARINLINE,"mcdisp.mf");fin=fopen(spinfile,"rb");
    if (fin==NULL) {fprintf(stderr,"ERROR - file %s not found \n",spinfile);exit(EXIT_FAILURE);}
   }
  instr[0]='#';  
  while(instr[strspn(instr," \t")]=='#'&&instr[strspn(instr," \t#")]!='!'){fgets(instr,MAXNOFCHARINLINE,fin);}
  extract(instr,"T",T); 
  extract(instr,"Ha",Hext[1]);
  extract(instr,"Hb",Hext[2]);
  extract(instr,"Hc",Hext[3]);
  info= new char [strlen(instr)+1];strcpy(info,instr);
  printf("#%s \n# reading mean field configuration mf=gj muB heff [meV]\n",instr);
  nofatoms=1;nofcomponents=3;nofthreads=0;
  extract(instr,"nofatoms",nofatoms); 
  extract(instr,"nofcomponents",nofcomponents); 
  mf=mfcf(1,1,1,nofatoms,nofcomponents); 
  if(mf.load(fin)==0)
   {fprintf(stderr,"ERROR loading mean field configuration\n");exit(EXIT_FAILURE);}
  fclose(fin);
 //********************************  
  savfilename= new char [strlen(file)+strlen(prefix)+11];
  snprintf(savfilename,MAXNOFCHARINLINE,"results/_%s%s",prefix,file);
  errno = 0;
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
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
  outS=0;
  qmin=0;qmax=0;deltaq=0;
 // ******************************** reading parameters  from mcdisp.par ****************************************************
  int i=0,hklblock=0,QxQyQzblock=0,j;
  printf("reading file %s\n",file);
  fin = fopen(file, "rb"); 
if (fin==NULL) {fprintf(stderr,"# Warning - file %s not found - using default values ! \n",file);
 emin=-100;emax=100;ki=0;kf=100;colcod[1]=5;colcod[2]=6;colcod[3]=7;colcod[4]=4;
nofhkls=0;save(file); 
}else
{   
  while (fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
  {++i; // i is used to estimate an upper boundary for the number of hkls in the hkl list 
     extract_with_prefix(instr,prefix,"emin",emin); 
     extract_with_prefix(instr,prefix,"emax",emax); 
     extract_with_prefix(instr,prefix,"ki",ki); 
     extract_with_prefix(instr,prefix,"kf",kf); 
     extract_with_prefix(instr,prefix,"calculate_magmoment_oscillation",calculate_magmoment_oscillation);
     extract_with_prefix(instr,prefix,"calculate_spinmoment_oscillation",calculate_spinmoment_oscillation);
     extract_with_prefix(instr,prefix,"calculate_orbmoment_oscillation",calculate_orbmoment_oscillation);
     extract_with_prefix(instr,prefix,"calculate_chargedensity_oscillation",calculate_chargedensity_oscillation);
     extract_with_prefix(instr,prefix,"calculate_spindensity_oscillation",calculate_spindensity_oscillation);
     extract_with_prefix(instr,prefix,"calculate_orbmomdensity_oscillation",calculate_orbmomdensity_oscillation);
     extract_with_prefix(instr,prefix,"calculate_phonon_oscillation",calculate_phonon_oscillation);
     extract_with_prefix(instr,prefix,"outS",outS);
     extract_with_prefix(instr,prefix,"nofthreads",nofthreads);
     for(int j=1;j<=usrdefcols[0];++j) // extract user defined output columns
     {snprintf(somestring,MAXNOFCHARINLINE,"out%i",usrdefcols[j]);
      extract(instr, somestring,colcod[usrdefcols[j]]);
     }

     hklblock+=1-extract_with_prefix(instr,prefix,"hmin",qmin[1]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"kmin",qmin[2]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"lmin",qmin[3]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"hmax",qmax[1]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"kmax",qmax[2]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"lmax",qmax[3]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"deltah",deltaq[1]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"deltak",deltaq[2]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"deltal",deltaq[3]); 
     if(hklblock==9){++nofhkllists;hklblock=0;i+=(int)ceil(fabs((qmax(1)-qmin(1))/deltaq(1)+1)*fabs((qmax(2)-qmin(2))/deltaq(2)+1)*fabs((qmax(3)-qmin(3))/deltaq(3)+1));}

     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qxmin",qmin[1]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qymin",qmin[2]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qzmin",qmin[3]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qxmax",qmax[1]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qymax",qmax[2]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qzmax",qmax[3]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"deltaQx",deltaq[1]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"deltaQy",deltaq[2]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"deltaQz",deltaq[3]); 
     if(QxQyQzblock==9){++nofhkllists;QxQyQzblock=0;i+=(int)ceil(fabs((qmax(1)-qmin(1))/deltaq(1)+1)*fabs((qmax(2)-qmin(2))/deltaq(2)+1)*fabs((qmax(3)-qmin(3))/deltaq(3)+1));}

     if(!extract_with_prefix(instr,prefix,"hklfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");while (fgets(hklfile,MAXNOFCHARINLINE,finhkl)!=NULL)++i;
                  fclose(finhkl);++nofhkllists;
                 }
     if(!extract_with_prefix(instr,prefix,"QxQyQzfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");while (fgets(hklfile,MAXNOFCHARINLINE,finhkl)!=NULL)++i;
                  fclose(finhkl);++nofhkllists;
                 }
     if(!extract_with_prefix(instr,prefix,"hklline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) Nsteps=21
                 {if(!extract(instr,"Nstp",N))i+=N;++nofhkllists;
                 }
     if(!extract_with_prefix(instr,prefix,"QxQyQzline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) Nsteps=21
                 {if(!extract(instr,"Nstp",N))i+=N;++nofhkllists;
                 }
     if(!extract_with_prefix(instr,prefix,"hklplane",hklline,MAXNOFCHARINLINE-1,1))  
                 {if(!extract(instr,"Nstp",N)&&!extract(instr,"Mstp",M))i+=N*M;++nofhkllists;
                 }
     if(!extract_with_prefix(instr,prefix,"QxQyQzplane",hklline,MAXNOFCHARINLINE-1,1))  
                 {if(!extract(instr,"Nstp",N)&&!extract(instr,"Mstp",M))i+=N*M;++nofhkllists;
                 }
  }
  fclose (fin);
 }
 // ************************************ end reading parameters *****************************************
 // check parameters
  if (ki==0) {if (kf==0) kf=100;
              fprintf(stdout,"#Calculating intensities for  kf=const=%4.4g/A\n",kf);
	     }
	     else
	     {kf=0;
	      fprintf(stdout,"#Calculating intensities for ki=const=%4.4g/A\n",ki);
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
  printf("# nofthreads=%i\n",nofthreads);
  // reread mcdisp.par creating the hkl list ******************************************************************
  hkls=new double *[i+10]; // dimension the list
  nofhkls=0;hklblock=0;QxQyQzblock=0;
  hklfile_start_index= new int [nofhkllists+1];hklfile_start_index[0]=nofhkllists;
  nofhkllists=0;Vector hkl(1,3),qijk(1,3);
  fin = fopen(file, "rb"); // if in mcdisp.par we find a hklfile= ... insert hkl from this file into list
            while (fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
               {// treat hklblocks
     hklblock+=1-extract_with_prefix(instr,prefix,"hmin",qmin[1]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"kmin",qmin[2]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"lmin",qmin[3]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"hmax",qmax[1]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"kmax",qmax[2]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"lmax",qmax[3]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"deltah",deltaq[1]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"deltak",deltaq[2]); 
     hklblock+=1-extract_with_prefix(instr,prefix,"deltal",deltaq[3]); 
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

     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qxmin",qmin[1]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qymin",qmin[2]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qzmin",qmin[3]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qxmax",qmax[1]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qymax",qmax[2]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"Qzmax",qmax[3]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"deltaQx",deltaq[1]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"deltaQy",deltaq[2]); 
     QxQyQzblock+=1-extract_with_prefix(instr,prefix,"deltaQz",deltaq[3]); 
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
                if(!extract_with_prefix(instr,prefix,"hklplane",hklline,MAXNOFCHARINLINE-1,1)) //#hklplane=h0=0 k0=1 l0=0 to hN=1 kN=1 lN=0 Nstp=21 to hM=1 kM=0 lM=3 Mstp=21  
                 { if(extract(instr,"h0",h0)){printf("error mcdisp reading mcdisp.par: in hklplane - h0 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"k0",k0)){printf("error mcdisp reading mcdisp.par: in hklplane - k0 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"l0",l0)){printf("error mcdisp reading mcdisp.par: in hklplane - l0 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"hN",hN)){printf("error mcdisp reading mcdisp.par: in hklplane - hN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"kN",kN)){printf("error mcdisp reading mcdisp.par: in hklplane - kN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"lN",lN)){printf("error mcdisp reading mcdisp.par: in hklplane - lN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"hM",hM)){printf("error mcdisp reading mcdisp.par: in hklplane - hM not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"kM",kM)){printf("error mcdisp reading mcdisp.par: in hklplane - kM not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"lM",lM)){printf("error mcdisp reading mcdisp.par: in hklplane - lM not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading mcdisp.par: in hklplane - Nstp  not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Mstp",M)){printf("error mcdisp reading mcdisp.par: in hklplane - Mstp  not found");exit (EXIT_FAILURE);}
                   printf("# ... hklplane (%g %g %g) to (%g %g %g) with %g points to (%g %g %g) with %g points\n",h0,k0,l0,hN,kN,lN,N,hM,kM,lM,M);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(i=1;i<=N;++i)for(j=1;j<=M;++j){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=h0+(hN-h0)*(i-1)/(N-1)+(hM-h0)*(j-1)/(M-1);
                                     hkls[nofhkls][2]=k0+(kN-k0)*(i-1)/(N-1)+(kM-k0)*(j-1)/(M-1);
                                     hkls[nofhkls][3]=l0+(lN-l0)*(i-1)/(N-1)+(lM-l0)*(j-1)/(M-1);
                                    }
                 }

                // treat QxQyQzplane statements
                if(!extract_with_prefix(instr,prefix,"QxQyQzplane",hklline,MAXNOFCHARINLINE-1,1)) //#QxyQzplane=Qx0=0 Qy0=1 Qz0=0 to QxN=1 QyN=1 QzN=0 Nstp=21 to QxM=1 QyM=0 QzM=3 Mstp=21
                 { if(extract(instr,"Qx0",h0)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - Qx0 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Qy0",k0)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - Qy0 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Qz0",l0)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - Qz0 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QxN",hN)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - QxN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QyN",kN)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - QyN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QzN",lN)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - QzN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QxM",hM)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - QxM not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QyM",kM)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - QyM not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QzM",lM)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - QzM not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - Nstp  not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Mstp",M)){printf("error mcdisp reading mcdisp.par: in QxyQzplane - Mstp  not found");exit (EXIT_FAILURE);}
                   printf("# ... QxyQzplane (%g %g %g) to (%g %g %g) with %g points to (%g %g %g) with %g points\n",h0,k0,l0,hN,kN,lN,N,hM,kM,lM,M);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(i=1;i<=N;++i)for(j=1;j<=M;++j){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     qijk(1)=h0+(hN-h0)*(i-1)/(N-1)+(hM-h0)*(j-1)/(M-1);
                                     qijk(2)=k0+(kN-k0)*(i-1)/(N-1)+(kM-k0)*(j-1)/(M-1);
                                     qijk(3)=l0+(lN-l0)*(i-1)/(N-1)+(lM-l0)*(j-1)/(M-1);
                                     ijk2hkl(hkl, qijk,abc);
                                     hkls[nofhkls][1]=hkl(1);
                                     hkls[nofhkls][2]=hkl(2);
                                     hkls[nofhkls][3]=hkl(3);
                                    }
                 }


                // treat hklline statements
                if(!extract_with_prefix(instr,prefix,"hklline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) N=21
                 { if(extract(instr,"h1",h1)){printf("error mcdisp reading mcdisp.par: in hklline - h1 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"k1",k1)){printf("error mcdisp reading mcdisp.par: in hklline - k1 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"l1",l1)){printf("error mcdisp reading mcdisp.par: in hklline - l1 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"hN",hN)){printf("error mcdisp reading mcdisp.par: in hklline - hN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"kN",kN)){printf("error mcdisp reading mcdisp.par: in hklline - kN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"lN",lN)){printf("error mcdisp reading mcdisp.par: in hklline - lN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading mcdisp.par: in hklline - Nstp  not found");exit (EXIT_FAILURE);}
                    printf("# ... hklline (%g %g %g) to (%g %g %g) with %g points\n",h1,k1,l1,hN,kN,lN,N);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(i=1;i<=N;++i){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     hkls[nofhkls][1]=h1+(hN-h1)*(i-1)/(N-1);
                                     hkls[nofhkls][2]=k1+(kN-k1)*(i-1)/(N-1);
                                     hkls[nofhkls][3]=l1+(lN-l1)*(i-1)/(N-1);
                                    }
                 }
                // treat QxQyQzline statements
                if(!extract_with_prefix(instr,prefix,"QxQyQzline",hklline,MAXNOFCHARINLINE-1,1))  // #!hklline=(h1=0 k1=0 l1=1) to (hN=0 kN=0 lN=2) N=21
                 { if(extract(instr,"Qx1",h1)){printf("error mcdisp reading mcdisp.par: in hklline - h1 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Qy1",k1)){printf("error mcdisp reading mcdisp.par: in hklline - k1 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Qz1",l1)){printf("error mcdisp reading mcdisp.par: in hklline - l1 not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QxN",hN)){printf("error mcdisp reading mcdisp.par: in hklline - hN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QyN",kN)){printf("error mcdisp reading mcdisp.par: in hklline - kN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"QzN",lN)){printf("error mcdisp reading mcdisp.par: in hklline - lN not found");exit (EXIT_FAILURE);}
                   if(extract(instr,"Nstp",N)){printf("error mcdisp reading mcdisp.par: in hklline - N  not found");exit (EXIT_FAILURE);}
                    printf("# ... QxQyQzline (%g %g %g)/A to (%g %g %g)/A with %g points\n",h1,k1,l1,hN,kN,lN,N);
                   ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                   for(i=1;i<=N;++i){++nofhkls; hkls[nofhkls]=new double [NOFHKLCOLUMNS+1];
                                     hkls[nofhkls][0]=3;
                                     qijk(1)=h1+(hN-h1)*(i-1)/(N-1);
                                     qijk(2)=k1+(kN-k1)*(i-1)/(N-1);
                                     qijk(3)=l1+(lN-l1)*(i-1)/(N-1);
                                     ijk2hkl(hkl, qijk,abc);
                                     hkls[nofhkls][1]=hkl(1);
                                     hkls[nofhkls][2]=hkl(2);
                                     hkls[nofhkls][3]=hkl(3);
                                    }
                 }
                 // treat hklfile statements
                if(!extract_with_prefix(instr,prefix,"hklfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                  read_hkl_list(finhkl,hkls,0,do_jqfile,abc);
                  fclose(finhkl);
                 }
                 // treat QxQyQzfile statements
                if(!extract_with_prefix(instr,prefix,"QxQyQzfile",hklfile,MAXNOFCHARINLINE-1,1))
                 {finhkl=fopen_errchk(hklfile,"rb");++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
                  read_hkl_list(finhkl,hkls,1,do_jqfile,abc);
                  fclose(finhkl);
                 }

              }
       fclose (fin);
       // now read also the hkls in mcdisp.par
      ++nofhkllists;hklfile_start_index[nofhkllists]=nofhkls+1;
      fin = fopen(file, "rb");read_hkl_list(fin,hkls,0,do_jqfile,abc); fclose(fin); 
  save();
      if(nofhkls==0){fprintf(stderr,"ERROR mcdisp: no hkl's found in mcdisp.par - please edit and insert\n");exit(EXIT_FAILURE);}      
}

//kopier-konstruktor 
inimcdis::inimcdis (const inimcdis & p)
{ do_jqf=p.do_jqf;
  savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);
 info= new char [strlen(p.info)+1];strcpy(info,p.info);
 prefix= new char [strlen(p.prefix)+1]; strcpy(prefix,p.prefix);  
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
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
  outS=p.outS;
    deltaq=p.deltaq;  
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  mf=mfcf(1,1,1,nofatoms,nofcomponents);mf=p.mf;T=p.T;
  nofhkls=p.nofhkls;
  int i,j;
      hkls=new double *[nofhkls+10];
      for (j=1;j<=nofhkls;++j) 
  	      {if ((int)p.hkls[j][0]==3){hkls[j]=new double [NOFHKLCOLUMNS+1];}
               else {hkls[j]=new double [(int)p.hkls[j][0]+1];}
               for(i=0;i<=p.hkls[j][0];++i)
         	    {hkls[j][i]=p.hkls[j][i];}
	      }
       int nofhkllists=p.hklfile_start_index[0];
       hklfile_start_index= new int [nofhkllists+1];hklfile_start_index[0]=nofhkllists;
      for (j=1;j<=nofhkllists;++j) hklfile_start_index[j]=p.hklfile_start_index[j]; 
   
}

//destruktor
inimcdis::~inimcdis ()
{delete []savfilename;delete []info;
 delete []prefix;
 int i;
//if (nofhkls==1)
 { for (i=1;i<=nofhkls;++i) 
   { delete []hkls[i];}
   delete []hkls;
   delete  []hklfile_start_index;
 }
}
