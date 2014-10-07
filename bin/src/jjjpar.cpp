// *************************************************************************
// ************************ class jjjpar     *******************************
// *************************************************************************
// jjjpar is a class to store all parameters associated
// with a specific ion, for example CEF parameters
// and exchange parameters 
// it is used by many programs in the package
// moreover, it loads also the user defined single ion module functions (linux only)
#include "jjjpar.hpp"
#include "perlparse.h"
#include "../../version"
#include<par.hpp>

#define MAXNOFNUMBERSINLINE 2500
#define MAXNOFCHARINLINE 7024

#define SMALL 1e-6   
#define SMALL_QUASIELASTIC_ENERGY 1e-6   //!!! must match SMALL_QUASIELASTIC_ENERGY in mcdisp.h and ionpars.hpp !!!
                     // because it is used to decide whether for small transition
		     // energy the matrix Mijkl contains wn-wn' or wn/kT

#include "jjjpar_basmodfunc.cpp" // basic sipf module functions
#include "jjjpar_observables.cpp" // function for physical observables
#include "jjjpar_intmod_kramer.cpp"   // some functions for module_type=1
#include "jjjpar_intmod_brillouin.cpp"// some functions for module_type=3
#include "jjjpar_intmod_cluster.cpp"// some functions for module_type=5





/************************************************************************************/
// returns total angular momentum quantum number J
/************************************************************************************/
double jjjpar::J()
{
 switch (module_type)
  {
   case 2:
   case 4: return (*iops).J;break;
   case 3:  return ABC[1];break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to calculate quantum number J \n");
            exit (EXIT_FAILURE);
  }
}

/************************************************************************************/
// returns stevens parameters of ion
/************************************************************************************/
Vector & jjjpar::tetan ()
{static Vector tt(1,6);
 tt=0;
 switch (module_type)
  {
   case 2:
   case 4: tt(2)=(*iops).alpha;tt(4)=(*iops).beta;tt(6)=(*iops).gamma;
            return tt;break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to calculate stevens parameters alpha beta gamma \n"); 
            exit (EXIT_FAILURE);
  }
}


// additional
// function to look if string s2 lies in string s1, checking the first n characters of s2
int strncomp(const char * s1,const char * s2, size_t n)
{size_t i;
 if (strlen(s1)>=n)
 {for (i=0;i<=strlen(s1)-n;++i)
  {if (strncmp(&s1[i],s2,n)==0){return 0;}
  }
 }
 return strncmp(s1,s2,n);
}




void jjjpar::increase_nofcomponents(int n) // increase nofcomponents by n
{int i,j,k,nold;
  nold=nofcomponents;
  nofcomponents+=n;
  mom.Resize(1,nofcomponents); 

  Matrix * jijstore;
  jijstore = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jijstore[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jijstore == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  for (i=1;i<=paranz;++i)
   {jijstore[i]=0;
    for (j=1;j<=nold;++j)
    {for (k=1;k<=nold;++k)
     {jijstore[i](j,k)=jij[i](j,k);
   }}}
 

 delete []jij;
  jij = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jij[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  for (i=1;i<=paranz;++i)
  {jij[i]=jijstore[i];}

  delete[] jijstore;
  fprintf(stderr,"Warning: increasing nofcomponents not tested  yet ... addition of parameter sets may be erroneous\n");
//  exit(0);

}


void jjjpar::add(jjjpar & b,Vector & abc) // add set b to this (abc: lattice constants)
{int i,j; 
 if(diagonalexchange==1&&b.diagonalexchange==0)diagonalexchange=0;
  if (nofcomponents!=b.nofcomponents)
  { fprintf (stderr, "class jjjpar: function add - nofcomponents does not match (check number of columns)\n"); 
    exit (EXIT_FAILURE);}

 for(i=1;i<=b.paranz;++i)
 {int found=0;
  for(j=1;j<=paranz&&found==0;++j)
  {if(Norm(dn[j]-b.dn[i])<SMALL)
    {//parameter found in list
     jij[j]+=b.jij[i];found=1;
    }
  }
  if (found==0){ // parameter not found in list !!!
                 jjjpar c(1,diagonalexchange,nofcomponents);
                 for(j=1;j<=paranz&&abc(1)*abc(1)*dn[j](1)*dn[j](1)+
		                    abc(2)*abc(2)*dn[j](2)*dn[j](2)+
				    abc(3)*abc(3)*dn[j](3)*dn[j](3)
				    <
				    abc(1)*abc(1)*b.dn[i](1)*b.dn[i](1)+
		                    abc(2)*abc(2)*b.dn[i](2)*b.dn[i](2)+
				    abc(3)*abc(3)*b.dn[i](3)*b.dn[i](3)
				    ;++j); //look for matching distance
		c.dn[1]=b.dn[i];
		c.jij[1]=b.jij[i];
		addpars(j,c);
               }
 }
}

// enlarge the set of parameters 
// inserting a set of exchange parameters
// into field at position number
void jjjpar::addpars (int number, jjjpar & addjjj)
{ Matrix * jijn;
  Vector * dnn;
  int i;
  jijn = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jijn[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  dnn = new Vector[paranz+1];for(i=0;i<=paranz;++i){dnn[i]=Vector(1,3);}
  
  for (i=1;i<=paranz;++i)
  {jijn[i]=jij[i];
   dnn[i]=dn[i];
  }
  
  if (diagonalexchange!=addjjj.diagonalexchange)
  { fprintf (stderr, "class jjjpar: function addpar - diagonalexchange does not match\n"); 
    exit (EXIT_FAILURE);}
  if (nofcomponents!=addjjj.nofcomponents)
  { fprintf (stderr, "class jjjpar: function addpar - nofcomponents does not match (check number of columns)\n"); 
    exit (EXIT_FAILURE);}
      
    paranz+=addjjj.paranz;  // increase parameters   
  
  delete []jij;
  delete []dn;
  delete []sublattice;
  dn = new Vector[paranz+1];for(i=0;i<=paranz;++i){dn[i]=Vector(1,3);}
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jij[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

// setup new field jij, dn
  for (i=1;i<number;++i)
  {jij[i]=jijn[i];dn[i]=dnn[i];}
  
  for (i=number;i<number+addjjj.paranz;++i)
  {jij[i]=addjjj.jij[i-number+1];dn[i]=addjjj.dn[i-number+1];}
  
  for (i=number+addjjj.paranz;i<=paranz;++i)
  {jij[i]=jijn[i-addjjj.paranz];dn[i]=dnn[i-addjjj.paranz];}
  delete []jijn;
  delete []dnn;
}

// remove neighbour from list
void jjjpar::delpar (int number)
{ Matrix * jijn;
  Vector * dnn;
  int * sublatticen;
  int i;
  --paranz;
  jijn = new Matrix[paranz+1];
  dnn = new Vector[paranz+1];
  sublatticen = new int[paranz+1];
  int offset=0;
  for (i=1;i<=paranz;++i)
  {if(i==number)offset=1;
   jijn[i]=jij[i+offset];
   dnn[i]=dn[i+offset];
   sublatticen[i]=sublattice[i+offset];
  }

  delete []jij;
  delete []dn;
  delete []sublattice;
  dn=dnn;jij=jijn;sublattice=sublatticen;
}


/************************************************************************************/
// save/get parameters 
/************************************************************************************/

//saving parameters to file
void jjjpar::save(FILE * file) 
{ int i,i1,j1;
  saveatom(file);
  fprintf(file,"#da[a]   db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
// save the exchange parameters to file (exactly paranz parameters!)
  for  (i=1;i<=paranz;++i)
  {fprintf(file,"%-+8.6g %-+8.6g %-+8.6g  ",myround(dn[i](1)),myround(dn[i](2)),myround(dn[i](3)));
    // format of matrix 
  // 11 22 33 12 21 13 31 23 32 (3x3 matrix)
  // 11 22 33 44 12 21 13 31 14 41 23 32 24 42 34 43 (4x4 matrix)
  // 11 22 33 44 55 12 21 13 31 14 41 15 51 23 32 24 42 25 52 34 43 35 53 45 54 (5x5 matrix)
  // etc ...
  //save diagonal components of exchange matrix
  for(i1=1;i1<=nofcomponents;++i1){fprintf(file,"%-+8.6e ",jij[i](i1,i1));}
  //save off-diagonal components of exchange matrix (if required)
  if (diagonalexchange==0){for(i1=1;i1<=nofcomponents-1;++i1)
                              {for(j1=i1+1;j1<=nofcomponents;++j1)
                               {fprintf(file,"%-+8.6e %-+8.6e ",jij[i](i1,j1),jij[i](j1,i1));
			       }
			      }
                          }
   fprintf(file,"\n");  
  }
}

void jjjpar::saveatom(FILE * file) 
{   fprintf(file,"#! da=%4.6g [a] db=%4.6g [b] dc=%4.6g [c] nofneighbours=%i diagonalexchange=%i sipffilename=%s\n",xyz(1),xyz(2),xyz(3),paranz,diagonalexchange,sipffilename);
}

//save single ion parameter file filename to path*
void jjjpar::save_sipf(const char * path)
{char * savfilename;
 savfilename= new char[strlen(sipffilename)+strlen(path)+2];
 strcpy(savfilename,path);
 strcpy(savfilename+strlen(path),sipffilename);
 FILE * fout; 
 fout = fopen_errchk (savfilename, "w");
 save_sipf(fout);
 fclose(fout);
 delete []savfilename;

}

void jjjpar::save_sipf(FILE * fout)
{ char  instr[MAXNOFCHARINLINE];
 int i;FILE * cfin;
 switch (module_type)
  {case 1: fprintf(fout,"#!MODULE=kramer\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module Kramer for  \n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#***************************************************************\n#\n");
           fprintf(fout,"#\n# this is a crystal field ground state doublet\n");
           fprintf(fout,"# module, parameters are the following 3 matrix\n# elements\n#\n");
           fprintf(fout,"# A=|<+-|Ja|-+>| B=|<+-|Jb|-+>| C=|<+-|Jc|+->|\n");
           fprintf(fout,"A=%10f \n B=%10f \n C=%10f\n\n",ABC(1),ABC(2),ABC(3));
            
          break;
   case 2: fprintf(fout,"#!MODULE=cfield\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module Cfield for  \n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#***************************************************************\n#\n");
           fprintf(fout,"#\n# crystal field paramerized in Stevens formalism\n#\n");
           (*iops).save(fout);
          break;
   case 3: fprintf(fout,"#!MODULE=brillouin\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module Brillouin for\n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#****************************************************************\n#\n");
           fprintf(fout,"#\n# single ion parameterized by Brillouin function\n");
           fprintf(fout,"# BJ(x) with angular momentum number J=S,\n# no crystal field\n#\n");
           fprintf(fout,"J = %g\n\n",ABC(1));
          break;
   case 4: fprintf(fout,"#!MODULE=so1ion\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module So1ion for  \n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#***************************************************************\n#\n");
           fprintf(fout,"#\n# crystal field paramerized in Stevens formalism\n#\n");
           (*iops).save(fout);
          break;
   case 5: //fprintf(fout,"#!MODULE=cluster\n#<!--mcphase.sipf-->\n");
           //fprintf(fout,"#***************************************************************\n");
           //fprintf(fout,"# Single Ion Parameter File for Module Cluster for\n");
           //fprintf(fout,"# %s\n",MCPHASVERSION);
           //fprintf(fout,"# - program to calculate static magnetic properties\n");
           //fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           //fprintf(fout,"# %s\n",MCDISPVERSION);
           //fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           //fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           //fprintf(fout,"# %s\n",MCDIFFVERSION);
           //fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           //fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           //fprintf(fout,"#****************************************************************\n#\n");
           //fprintf(fout,"#\n# single ion subsystem consists of cluster of ions\n");
           //fprintf(fout,"# cluster structure is desribed in file\n#\n");
           //fprintf(fout,"#!structurefile = %s\n\n",clusterfile);
           //break;
           {char clustsavfile[MAXNOFCHARINLINE];sprintf(clustsavfile,"results/_%s",clusterfilename);
           (*clusterpars).save(clustsavfile);
           (*clusterpars).save_sipfs("results/_");}
   default: // in case of external single ion module just save a copy of the input file 
           char *token;
           cfin=fopen_errchk(sipffilename,"rb");
           while(feof(cfin)==false){fgets(instr, MAXNOFCHARINLINE, cfin);
                      // strip /r (dos line feed) from line if necessary
                      while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
                      fprintf(fout,"%s",instr);
                                    }
           fclose(cfin);
   }

  if(module_type>0&&module_type!=5) // in case of internal modules save common information
   {fprintf(fout,"#----------------\n# number of electrons in unfilled shell gJ\n#----------------\nnof_electrons=%i\n\n",nof_electrons);
    fprintf(fout,"#----------------\n# Lande factor gJ\n#----------------\nGJ=%g\n\n",gJ);
    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"# Neutron Scattering Length (10^-12 cm) (can be complex)\n");
    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"SCATTERINGLENGTHREAL=%g\nSCATTERINGLENGTHIMAG=%g\n",SLR,SLI);
    fprintf(fout,"#  ... note: - if an occupancy other than 1.0 is needed, just reduce \n");
    fprintf(fout,"#              the scattering length linear accordingly\n\n");

    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)\n");
    fprintf(fout,"#                      with s=sin(theta)/lambda=Q/4pi\n");
    fprintf(fout,"# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>\n");
    fprintf(fout,"# unit of DWF is [A^2]\n");
    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"DWF=%g\n",DWF);

    fprintf(fout,"#--------------------------------------------------------------------------------------\n");
    fprintf(fout,"# Neutron Magnetic Form Factor coefficients - thanks to J Brown\n");
    fprintf(fout,"#   d = 2*pi/Q      \n");
    fprintf(fout,"#   s = 1/2/d = Q/4/pi   \n");   
    fprintf(fout,"#   sin(theta) = lambda * s\n");
    fprintf(fout,"#    s2= s*s = Q*Q/16/pi/pi\n");
    fprintf(fout,"#\n");
    fprintf(fout,"#   <j0(Q)>=   FFj0A*EXP(-FFj0a*s2) + FFj0B*EXP(-FFj0b*s2) + FFj0C*EXP(-FFj0c*s2) + FFj0D\n");
    fprintf(fout,"#   <j2(Q)>=s2*(FFj2A*EXP(-FFj2a*s2) + FFj2B*EXP(-FFj2b*s2) + FFj2C*EXP(-FFj2c*s2) + FFj2D\n");
    fprintf(fout,"#   <j4(Q)>=s2*(FFj4A*EXP(-FFj4a*s2) + FFj4B*EXP(-FFj4b*s2) + FFj4C*EXP(-FFj4c*s2) + FFj4D\n");
    fprintf(fout,"#   <j6(Q)>=s2*(FFj6A*EXP(-FFj6a*s2) + FFj6B*EXP(-FFj6b*s2) + FFj6C*EXP(-FFj6c*s2) + FFj6D\n");
    fprintf(fout,"#\n");
    fprintf(fout,"#   Dipole Approximation for Neutron Magnetic Formfactor:\n");
    fprintf(fout,"#        -Spin Form Factor       FS(Q)=<j0(Q)>\n");
    fprintf(fout,"#        -Angular Form Factor    FL(Q)=<j0(Q)>+<j2(Q)>\n");
    fprintf(fout,"#        -Rare Earth Form Factor F(Q) =<j0(Q)>+<j2(Q)>*(2/gJ-1)\n\n");
    fprintf(fout,"#--------------------------------------------------------------------------------------\n");
    fprintf(fout,"FFj0A=%+7.4f FFj0a=%+7.4f FFj0B=%+7.4f FFj0b=%+7.4f FFj0C=%+7.4f FFj0c=%+7.4f FFj0D=%+7.4f\n",magFFj0[1],magFFj0[2],magFFj0[3],magFFj0[4],magFFj0[5],magFFj0[6],magFFj0[7]);
    fprintf(fout,"FFj2A=%+7.4f FFj2a=%+7.4f FFj2B=%+7.4f FFj2b=%+7.4f FFj2C=%+7.4f FFj2c=%+7.4f FFj2D=%+7.4f\n",magFFj2[1],magFFj2[2],magFFj2[3],magFFj2[4],magFFj2[5],magFFj2[6],magFFj2[7]);
    fprintf(fout,"FFj4A=%+7.4f FFj4a=%+7.4f FFj4B=%+7.4f FFj4b=%+7.4f FFj4C=%+7.4f FFj4c=%+7.4f FFj4D=%+7.4f\n",magFFj4[1],magFFj4[2],magFFj4[3],magFFj4[4],magFFj4[5],magFFj4[6],magFFj4[7]);
    fprintf(fout,"FFj6A=%+7.4f FFj6a=%+7.4f FFj6B=%+7.4f FFj6b=%+7.4f FFj6C=%+7.4f FFj6c=%+7.4f FFj6D=%+7.4f\n",magFFj6[1],magFFj6[2],magFFj6[3],magFFj6[4],magFFj6[5],magFFj6[6],magFFj6[7]);
    fprintf(fout,"\n\n");

  if(abs(Zc)>1e-10){
    fprintf(fout,"#----------------------------------------------------------------------\n");
    fprintf(fout,"# coefficients of Z(K') according to Lovesey (Neutron Scattering) vol.2\n");
    fprintf(fout,"# chapter 11.6.1 page 233: Z(K)= ZKcK-1 * <jK-1(Q)> + ZKcK+1 * <jK+1(Q)>\n");
    fprintf(fout,"#  ... these coefficients are needed to go beyond dipolar approx.\n");
    fprintf(fout,"#      for the neutron magnetic formfactor in rare earth ions\n");
    fprintf(fout,"#----------------------------------------------------------------------\n");
    fprintf(fout,"Z1c0=%+10.8f  Z1c2=%+10.8f\n",Zc(1),Zc(2));
    fprintf(fout,"		  Z3c2=%+10.8f  Z3c4=%+10.8f\n",Zc(3),Zc(4)); 
    fprintf(fout,"				    Z5c4=%+10.8f  Z5c6=%+10.8f\n",Zc(5),Zc(6)); 
    fprintf(fout,"						      Z7c6=%+10.8f\n\n",Zc(7));
                    }

  if(abs(Np)>1e-10){fprintf(fout,"#---------------------------------------------------------------------------------------------------\n");
                    fprintf(fout,"# radial wave function parameters, for transition metal ions the the values are tabulated in\n");
                    fprintf(fout,"# Clementi & Roetti Atomic data and nuclear data tables 14 (1974) 177-478, the radial wave\n");
                    fprintf(fout,"# function is expanded as R(r)=sum_p Cp r^(Np-1) . exp(-XIp r) . (2 XIp)^(Np+0.5) / sqrt((2Np)!)\n");
                    fprintf(fout,"# for rare earth ions see Freeman & Watson PR 127(1962)2058, Sovers J. Phys. Chem. Sol. 28(1966)1073\n");
                    fprintf(fout,"#---------------------------------------------------------------------------------------------------\n");
                    for(i=Np.Lo();i<=Np.Hi();++i){if(Np(i)!=0){fprintf(fout,"N%i=%i XI%i=%g C%i=%g\n",i,(int)Np(i),i,Xip(i),i,Cp(i));}
                                                 }
                   fprintf(fout,"\n");
                   }

   }

}


/*****************************************************************************************/
//constructor with file handle of mcphas.j
jjjpar::jjjpar(FILE * file,int nofcomps) 
{ jl_lmax=6;
  char instr[MAXNOFCHARINLINE],exchangeindicesstr[MAXNOFCHARINLINE];
  sipffilename= new char [MAXNOFCHARINLINE];
  clusterfilename=new char [MAXNOFCHARINLINE];
  int i,j,i1,j1,k1;
  int symmetricexchange=0,indexexchangenum=0;
  Matrix exchangeindices;
  float nn[MAXNOFNUMBERSINLINE];
  nn[0]=MAXNOFNUMBERSINLINE;
  xyz=Vector(1,3);
  cnst= Matrix(0,6,-6,6);set_zlm_constants(cnst);
  i=6; FF_type=0;
  while(i>0){fgets_errchk (instr, MAXNOFCHARINLINE, file);
             if(instr[strspn(instr," \t")]!='#'){fprintf (stderr, "Error reading mcphas.j - exchangeparameters start before all variables (da,db,dc,nofneighbours,diagonalexchange and sipffilename) have been given\n");
                                                 exit (EXIT_FAILURE);}
             i+=extract(instr,"x",xyz[1])-1;
             i+=extract(instr,"y",xyz[2])-1;
             i+=extract(instr,"z",xyz[3])-1;
             i+=extract(instr,"da",xyz[1])-1;
             i+=extract(instr,"db",xyz[2])-1;
             i+=extract(instr,"dc",xyz[3])-1;
             i+=extract(instr,"nofneighbours",paranz)-1;
             i+=extract(instr,"diagonalexchange",diagonalexchange)-1;
             i+=extract(instr,"cffilename",sipffilename,(size_t)MAXNOFCHARINLINE)-1;
             i+=extract(instr,"sipffilename",sipffilename,(size_t)MAXNOFCHARINLINE)-1;
            }

 // MDL 29.08.10 - Added to check if exchange parameters are indexed.
 //   - Used: set diagonalexchange=2. The add a line:
 //        #! symmetricexchange=1 nofcomponents=3 indexexchange= JaJb JbJc
 //     Or
 //        #! symmetricexchange=0 nofcomponents=3 indexexchange= 1,2 2,1 2,3 3,2
 //     Where symmetricexchange is 1 or 0 depending on whether you want the non-diagonal elements of the exchange
 //     matrix to be symmetric or not. The column index can be either a string or a list of coordinates (indices)
 //     Note that you if you don't specify symmetricexchange, it is assumed to be 0 (nonsymmetric)
 //     Note also that the JaJb syntax is case sensitive - J must be upper case so that Jj is allowed. Index must be lower case
 if(diagonalexchange==2) { 
    i=1; while(i>0) { 
      if(instr[strspn(instr," \t")]!='#') { 
         fprintf (stderr, "Error reading mcphas.j - diagonalexchange==2, but not indexexchange parameters found\n"); exit (EXIT_FAILURE);}
      fgets_errchk (instr, MAXNOFCHARINLINE, file); 
      extract(instr,"symmetricexchange",symmetricexchange);
      if(extract(instr,"indexexchange",exchangeindicesstr,MAXNOFCHARINLINE)==0) { strcpy(exchangeindicesstr,instr); break; }
    }
    indexexchangenum=get_exchange_indices(exchangeindicesstr,&exchangeindices);
 }

 // fgets_errchk (instr, MAXNOFCHARINLINE, file); //removed by MR 30.4.2010

// read single ion parameter file and see which type it is (internal module or loadable)
  transitionnumber=1;
  
  // must be set before getting pars from sipffile (cluster module needs it)
  nofcomponents=nofcomps; // default value for nofcomponents - (important in case nofparameters=0)

  //start reading again at the beginning of the file to get formfactors, debye waller factor
  get_parameters_from_sipfile(sipffilename);

             dn = new Vector[paranz+1];if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);} // 4 lines moved here to make destructor work MR 30.3.10
             for(i1=0;i1<=paranz;++i1){dn[i1]=Vector(1,3);}
             sublattice = new int[paranz+1];if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
             jij = new Matrix[paranz+1];if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
// read the exchange parameters from file (exactly paranz parameters!)
  for  (i=1;i<=paranz;++i)
  {while((j=inputline(file, nn))==0&&feof(file)==0){}; // returns 0 if comment line or eof, exits with error, if input string too long
   // Additional check to see if we are on the last neighbour, as McPhaseExplorer generates bad files without EOL at the end unless you add an empty line
   if(feof(file)!=0 && (i<paranz||(j-3)<nofcomponents)) { 
                      fprintf (stderr, "Error in jjjpar.cpp: input jjj parameters - \n");
                      fprintf(stderr," end of file reached while reading exchange parameter %i(%i)\n",i,paranz);
                      exit (EXIT_FAILURE);
                    }
    if(i==1){// determine nofcomponents from number of parameters read in first line of mcphas.j
             if(diagonalexchange==1){nofcomponents=j-3;}else if(diagonalexchange==0){nofcomponents=(int)sqrt((double)(j-3));}
             if(module_type==1)
	     {// check dimensions of vector if internal kramers is used
              if(nofcomponents!=3)
              {fprintf(stderr,"Error reading mcphas.j: number of dimensions (not equal 3) not compatible with internal single ion module kramer - check number of columns in file mcphas.j\n");
               exit(EXIT_FAILURE);}
             }
             // dimension arrays
             for(i1=0;i1<=paranz;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
            }
   //check if correct number of columns has been read	        
    if((diagonalexchange==1&&nofcomponents!=j-3)||(diagonalexchange==0&&nofcomponents!=(int)sqrt((double)(j-3)))||(diagonalexchange==2&&indexexchangenum!=(j-3)))
              {fprintf(stderr,"Error reading mcphas.j line %i: check number of columns\n",i);
               exit(EXIT_FAILURE);}

  mom=Vector(1,nofcomponents); 

   //(1-3) give the absolute coordinates of the neighbour and are transformed here to relative
   // coordinates !!
   dn[i](1) = nn[1];dn[i](2) = nn[2];dn[i](3) = nn[3];
   jij[i]=0;

  // format of matrix 
  // 11 22 33 12 21 13 31 23 32 (3x3 matrix)
  // 11 22 33 44 12 21 13 31 14 41 23 32 24 42 34 43 (4x4 matrix)
  // 11 22 33 44 55 12 21 13 31 14 41 15 51 23 32 24 42 25 52 34 43 35 53 45 54 (5x5 matrix)
  // etc ...
  //read diagonal components of exchange matrix
  if (diagonalexchange==1||diagonalexchange==0)
  for(i1=1;i1<=nofcomponents;++i1){jij[i](i1,i1)= nn[i1+3];}
  //read off-diagonal components of exchange matrix (if required)
  if (diagonalexchange==0){k1=3+nofcomponents;
                           for(i1=1;i1<=nofcomponents-1;++i1)
                              {for(j1=i1+1;j1<=nofcomponents;++j1)
                               {++k1;jij[i](i1,j1)= nn[k1];
			        ++k1;jij[i](j1,i1)= nn[k1];
			       }
			      }
                          }
  if (diagonalexchange==2 && indexexchangenum>0)  {
     double dt=0.; int ii,jj;
     for(i1=1; i1<=indexexchangenum; i1++) {
        ii=(int)exchangeindices(i1,1); jj=(int)exchangeindices(i1,2); jij[i](ii,jj) = nn[i1+3]; 
        if(i==paranz) { if(ii!=jj) diagonalexchange=0; }} if(i==paranz&&diagonalexchange!=0) diagonalexchange=1;
     if (symmetricexchange==1) {
        for(i1=1; i1<=nofcomponents; i1++) for(j1=i1+1; j1<=nofcomponents; j1++) if(i1!=j1) {
           dt=0.; if(jij[i](j1,i1)!=0) dt=jij[i](j1,i1); else if(jij[i](i1,j1)!=0) dt=jij[i](i1,j1);
           if(dt!=0) { jij[i](j1,i1)=dt; jij[i](i1,j1)=dt; }
     }}}
  }
  for(unsigned int ui=MAXSAVEQ; ui--; ) { Qsaved[ui]=DBWQsaved[ui]=1e16; Fsaved[ui]=DBWsaved[ui]=0; } nsaved=DBWnsaved=MAXSAVEQ-1;
  for(int ii=0; ii<52; ii++) opmatM[ii] = 0;
}

// constructor with filename of singleion parameter  used by mcdiff and charges-chargeplot
jjjpar::jjjpar(double x,double y,double z, char * sipffile, int n)
{xyz=Vector(1,3);xyz(1)=x;xyz(2)=y;xyz(3)=z;jl_lmax=6;
  jij=0; dn=0; sublattice=0;paranz=0;diagonalexchange=1;
  mom=Vector(1,9); mom=0; nofcomponents=n;
  sipffilename= new char [MAXNOFCHARINLINE];
  clusterfilename=new char [MAXNOFCHARINLINE];
  strcpy(sipffilename,sipffile);
  get_parameters_from_sipfile(sipffilename);
   cnst= Matrix(0,6,-6,6);set_zlm_constants(cnst);
  for(unsigned int ui=MAXSAVEQ; ui--; ) { Qsaved[ui]=DBWQsaved[ui]=1e16; Fsaved[ui]=DBWsaved[ui]=0; } nsaved=DBWnsaved=MAXSAVEQ-1;
  for(int ii=0; ii<52; ii++) opmatM[ii] = 0;

}

// constructor with positions scattering length dwf
jjjpar::jjjpar(double x,double y,double z, double slr,double sli, double dwf)
{xyz=Vector(1,3);xyz(1)=x;xyz(2)=y;xyz(3)=z;jl_lmax=6;
 mom=Vector(1,9); mom=0; FF_type=0;
 DWF=dwf;SLR=slr;SLI=sli;
  magFFj0=Vector(1,7);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,7);magFFj2=0;
  magFFj4=Vector(1,7);magFFj4=0;
  magFFj6=Vector(1,7);magFFj6=0;
  Zc=Vector(1,7);Zc=0; 
   cnst= Matrix(0,6,-6,6);set_zlm_constants(cnst);
   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;
   r2=0;r4=0;r6=0;
   modulefilename=new char[MAXNOFCHARINLINE];
  nof_electrons=0; // no electorns by default
  paranz=0;
  sipffilename= new char [MAXNOFCHARINLINE];
  clusterfilename=new char [MAXNOFCHARINLINE];
  module_type=1;
  for(unsigned int ui=MAXSAVEQ; ui--; ) { Qsaved[ui]=DBWQsaved[ui]=1e16; Fsaved[ui]=DBWsaved[ui]=0; } nsaved=DBWnsaved=MAXSAVEQ-1;
  for(int ii=0; ii<52; ii++) opmatM[ii] = 0;
}

//constructor without file
jjjpar::jjjpar(int n,int diag,int nofmom) 
{ sipffilename= new char [MAXNOFCHARINLINE];jl_lmax=6;
  clusterfilename=new char [MAXNOFCHARINLINE];
  diagonalexchange=diag;
  paranz=n;xyz=Vector(1,3);xyz=0; 
   cnst= Matrix(0,6,-6,6);set_zlm_constants(cnst);
  int i1;r2=0;r4=0;r6=0;
  module_type=1;ABC=Vector(1,3);ABC=0;
  transitionnumber=1;
  nofcomponents=nofmom;
  mom=Vector(1,nofcomponents);
  mom=0; FF_type=0;
  nof_electrons=0;// no electorns by default
  modulefilename=new char[MAXNOFCHARINLINE];
  dn = new Vector[n+1];for(i1=0;i1<=n;++i1){dn[i1]=Vector(1,3);}
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[n+1];for(i1=0;i1<=n;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  magFFj0=Vector(1,7);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,7);magFFj2=0;
  magFFj4=Vector(1,7);magFFj4=0;
  magFFj6=Vector(1,7);magFFj6=0;
  Zc=Vector(1,7);Zc=0;
   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;
  DWF=0;gJ=0;maxE=1e10;pinit=0;ninit=1e10;
  for(unsigned int ui=MAXSAVEQ; ui--; ) { Qsaved[ui]=DBWQsaved[ui]=1e16; Fsaved[ui]=DBWsaved[ui]=0; } nsaved=DBWnsaved=MAXSAVEQ-1;
  for(int ii=0; ii<52; ii++) opmatM[ii] = 0;

}

//copy constructor
jjjpar::jjjpar (const jjjpar & pp)
{ int i;jl_lmax=pp.jl_lmax;
  xyz=Vector(1,3);
  nofcomponents=pp.nofcomponents;
  mom=Vector(1,nofcomponents); 
  xyz=pp.xyz;paranz=pp.paranz;
   cnst= Matrix(0,6,-6,6);set_zlm_constants(cnst);
  SLR=pp.SLR;SLI=pp.SLI;
  FF_type=pp.FF_type;
  nof_electrons=pp.nof_electrons;
  modulefilename=new char[MAXNOFCHARINLINE];
  strncpy (modulefilename,pp.modulefilename, MAXNOFCHARINLINE-1);
  clusterfilename=new char [MAXNOFCHARINLINE];
  strncpy(clusterfilename,pp.clusterfilename, MAXNOFCHARINLINE-1);
  diagonalexchange=pp.diagonalexchange;
  gJ=pp.gJ;module_type=pp.module_type;
  ninit=pp.ninit;maxE=pp.maxE;pinit=pp.pinit;
  Np=pp.Np; Xip=pp.Xip;Cp=pp.Cp;
  r2=pp.r2;r4=pp.r4;r6=pp.r6;

  transitionnumber=pp.transitionnumber;
  sipffilename= new char [strlen(pp.sipffilename)+1];
  strcpy(sipffilename,pp.sipffilename);
  if (pp.module_type==3||pp.module_type==1||pp.module_type==0)  ABC=pp.ABC;
  if ((pp.module_type==5||pp.module_type==3||pp.module_type==1||pp.module_type==0) &&
      (pp.Icalc_parstorage.Cols()>0) && (pp.Icalc_parstorage.Rows()>0))
  {  Icalc_parstorage = ComplexMatrix(pp.Icalc_parstorage.Rlo(),pp.Icalc_parstorage.Rhi(),pp.Icalc_parstorage.Clo(),pp.Icalc_parstorage.Chi());
     Icalc_parstorage = pp.Icalc_parstorage;
  }
  if (pp.module_type==2||pp.module_type==4)  {iops=new ionpars(*pp.iops);//((int)(2*(*pp.iops).J+1));iops=pp.iops;
                           int dj;dj=(int)(2*J()+1);
                           est=ComplexMatrix(0,dj,1,dj);est=pp.est;
                           Icalc_parstorage=ComplexMatrix(0,dj,1,dj);Icalc_parstorage=pp.Icalc_parstorage;
                           }
   if (pp.module_type==5) {clusterpars=new par(*pp.clusterpars);dim=pp.dim;
                          // est=ComplexMatrix(pp.est.Rlo(),pp.est.Rhi(),pp.est.Clo(),pp.est.Chi());est=pp.est;
                          // Ia
                          // cluster_ini_Imat(); 
                          dim=1; Vector Hxc(1,(*clusterpars).nofcomponents);Vector Hext(1,3);
                          dnn= new int [(*clusterpars).nofatoms+1];
                          // determine dimension of H matrix
                          for (int n=1;n<=(*clusterpars).nofatoms;++n)
                          {dnn[n]=(*(*clusterpars).jjj[n]).opmat(1,Hxc,Hext).Rhi();
                           dim*=dnn[n];
                          }
                          Ia= new zsMat<double> * [nofcomponents+1];
                          for(int n = 1;n<=nofcomponents;++n){Ia[n]=new zsMat<double>(dim,dim);(*Ia[n])=(*pp.Ia[n]);}
                          // cluster_M
                          cluster_M= new zsMat<double> * [3+3*(*clusterpars).nofatoms+1];
                          for(int a=0;a<=(*clusterpars).nofatoms;++a)for(int n=1;n<=3;++n)
                           {int index_M=a*3+n; // a .... atom index  n ... xyz components of magnetic moment
                            cluster_M[index_M]=new zsMat<double>(dim,dim);
                            (*cluster_M[index_M])=(*pp.cluster_M[index_M]);
                           }                        
                          clusterH = new Matrix(1,dim,1,dim); *clusterH = *pp.clusterH; 
                          oldHext = new Vector(1,3); *oldHext = *pp.oldHext;
                          }
  
//#ifdef __linux__
/*  if (module_type==0)
  {char * error;
   handle=dlopen (sipffilename,RTLD_NOW | RTLD_GLOBAL);
   if (!handle){fprintf (stderr, "Could not load dynamic library\n");
               if ((error=dlerror())!=NULL) 
	         {fprintf (stderr,"%s\n",error);}
	       exit (EXIT_FAILURE);
	      }
*/
   I=pp.I;   du=pp.du;
   estates=pp.estates;    Icalc_parameter_storage=pp.Icalc_parameter_storage;
   mq=pp.mq;    ddnn=pp.ddnn;
   p=pp.p;dP1=pp.dP1;
   m=pp.m;dm1=pp.dm1;
   L=pp.L;dL1=pp.dL1;
   S=pp.S;dS1=pp.dS1;
   cd_m=pp.cd_m;cd_dm=pp.cd_dm;
   sd_m=pp.sd_m;sd_dm=pp.sd_dm;
   od_m=pp.od_m;od_dm=pp.od_dm;
/*  }*/
//#endif
  magFFj0=Vector(1,7);magFFj0=pp.magFFj0;
  magFFj2=Vector(1,7);magFFj2=pp.magFFj2;
  magFFj4=Vector(1,7);magFFj4=pp.magFFj4;
  magFFj6=Vector(1,7);magFFj6=pp.magFFj6;
  Zc=Vector(1,7);Zc=pp.Zc;
  DWF=pp.DWF;  
int i1;
//dimension arrays
  jij = new Matrix[paranz+1];for(i1=0;i1<=paranz;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  dn = new Vector[paranz+1];for(i1=0;i1<=paranz;++i1){dn[i1]=Vector(1,3);}
  if (dn == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  for (i=1;i<=paranz;++i)
  {jij[i]=pp.jij[i];dn[i]=pp.dn[i];sublattice[i]=pp.sublattice[i];}
  for(unsigned int ui=MAXSAVEQ; ui--; ) { Qsaved[ui]=DBWQsaved[ui]=1e16; Fsaved[ui]=DBWsaved[ui]=0; } nsaved=DBWnsaved=MAXSAVEQ-1;
  
  for(i1=0; i1<52; i1++) { 
     if(pp.opmatM[i1]==0) opmatM[i1]=0; 
     else { 
        opmatM[i1] = new Matrix((*pp.opmatM[i1]).Clo(),(*pp.opmatM[i1]).Chi(),(*pp.opmatM[i1]).Rlo(),(*pp.opmatM[i1]).Rhi()); 
        (*opmatM[i1]) = (*pp.opmatM[i1]); } }
  if(pp.dyn_opmat!=NULL && module_type==0) {
    // For the cluster module, during the MF loop, we need to be able to call the external module to recalculate the Hamiltonian for different Hxc, Hext.
    #ifdef __MINGW32__
      handle=LoadLibrary(pp.modulefilename);
      if ((intptr_t)handle<= HINSTANCE_ERROR){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n"); exit (EXIT_FAILURE); } 
      dyn_opmat=(int(*)(int*,char**,Vector*,Vector*,Matrix*))GetProcAddress(handle,"opmat");
      if (dyn_opmat==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function opmat not possible - continuing\n",(int)GetLastError(),modulefilename);}
    #else
      char * error;
      handle=dlopen (pp.modulefilename,RTLD_NOW | RTLD_GLOBAL);
      if (!handle){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n"); if ((error=dlerror())!=NULL) {fprintf (stderr,"%s\n",error);} exit (EXIT_FAILURE); } 
      *(void **)(&dyn_opmat)=dlsym(handle,"opmat");
      if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dyn_opmat=NULL;}
    #endif
  }
}


//destruktor
jjjpar::~jjjpar ()
{// printf("hello destruktor jjjpar\n");  
   if(jij!=0)        delete []jij; //will not work in linux
   if(dn!=0)         delete []dn;  // will not work in linux
   if(sublattice!=0) delete []sublattice;
   delete []sipffilename;// will not work in linux
  delete []clusterfilename;
  delete []modulefilename;// will not work in linux
   if (module_type==5) {for(int i=1;i<=nofcomponents;++i)
                            { delete Ia[i];}
                       for(int a=0;a<=(*clusterpars).nofatoms;++a)for(int n=1;n<=3;++n)
                        {int index_M=a*3+n; // a .... atom index  n ... xyz components of magnetic moment
                         delete cluster_M[index_M];}
                        delete Ia;delete cluster_M; delete []dnn;
                        delete clusterpars;                         
                        delete clusterH; delete oldHext;
                       }
   if (module_type==2||module_type==4) delete iops;
  for(int ii=0; ii<52; ii++) 
     if(opmatM[ii]!=0) { delete opmatM[ii]; opmatM[ii]=0; }
 
//#ifdef __linux__
// i#ifdef __linux__f (module_type==0)dlclose(handle);
//#endif
// printf("hello end destruktor jjjpar\n");  
 
 }


// for test
/*int main(int argc, char **argv)
{
  jjjpar *d=new jjjpar(1,2,3, argv[1]);
}*/
