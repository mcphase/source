/***********************************************************************
 *
 * reduce_unitcell.c - program to reduce unit cell by removing 
 *
 ***********************************************************************/


#include "par.hpp"
#include "martin.h"
#include "mcphas.h"

int verbose=0;
const char * filemode="w";

#ifdef _THREADS
#undef _THREADS
#include "mcphas_fecalc.c"
#define _THREADS
#else
#include "mcphas_fecalc.c"
#endif

// ***************************************************************************
// routine for the case of removing phononic degrees of freedom and
// - renormalising elastic constants
// - calculating (phonon induced) magnetoelastic and quadrupolar interactions
// ***************************************************************************
void delphonons(par & a)
{fprintf(stderr,"# deleting phonons");
// 1. renormalise elastic constants using only phonon degrees of freedom
//   and calling fecalc with doeps for various applied external stresses ...
double U,Eelastic,spinchange; int r; 
 par phon(a); // for phonon --> elastic constants
par pw(a);  // pw working parameterset ---> multipolar and magnetoelastic interactions

 phon.set_nofcomponents(3);
// remove nonphononic degrees of freedom from phon
Vector u0(1,3),Hxc(1,3),Hext(1,6);double T=1;
Vector h1(1,phon.cs.nofcomponents);
Vector Happ(1,HEXT_DIMENSION);Happ=0;
Matrix dis(1,1,1,1);dis=0;dis(1,1)=1;
for(int n=phon.cs.nofatoms;n>0;--n)
{h1=0;(*phon.jjj[n]).Icalc_parameter_storage_init(h1,Happ,T); // initialize eigenstate matrix
 if(!((*phon.jjj[n]).module_type==external_class&&(*phon.jjj[n]).pcalc(u0,  T,  Hxc,Hext,(*phon.jjj[n]).Icalc_parstorage)))
  {phon.delatom(-n,dis);}
}
double U0,Ua,Ub,Uc; // reference energy

char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
 inipar ini("mcphas.ini",prefix,"reduce_unitcell");ini.doeps=true;
// ---------------------------------
// 1. compute elastic constants ....
// ---------------------------------
// only keep phononic degrees of freedom in phon and
// apply various stresses 
 {         spincf sps(1,1,1,phon.cs.nofatoms,3);
          mfcf mf(1,1,1,phon.cs.nofatoms,3);
          Matrix s(1,6,1,6); //myPrintVector(sps.epsilon,"Initial Strain");
if(fecalc(U0,Eelastic,r,spinchange,Happ,T,ini,phon,sps,mf)>FEMIN_INI)
 {fprintf(stderr,"Error reduce_unitcell option delatoms phon: free energy not stable for elastic constants estimate - modify mcphas.ini and restart\n");exit(EXIT_FAILURE);}
           double cel=0.1;  // fixed stress to apply in GPa
          for(int n=1;n<=6;++n){Happ(6+n)+=cel;
 if(fecalc(U,Eelastic,r,spinchange,Happ,T,ini,phon,sps,mf)>FEMIN_INI)
 {fprintf(stderr,"Error reduce_unitcell option delatoms phon: free energy not stable for elastic constants estimate - modify mcphas.ini and restart\n");exit(EXIT_FAILURE);}
                                Happ(6+n)-=cel;//myPrintVector(sps.epsilon,"Strain x=1-6");
 //                               printf("U+sigma.eps-U0=%12.12g U0=%12.12g -sigma.eps=%12.12g\n",(U-U0)*phon.cs.nofatoms+cel*sps.epsilon(n)*phon.cs.pVol()/1.60218e-1,U0*phon.cs.nofatoms,-cel*sps.epsilon(n)*phon.cs.pVol()/1.60218e-1);
//if(n==4)sps.print(stdout);
                                for(int m=1;m<=6;++m)s(n,m)=sps.epsilon(m)/cel;
                               }  // sigma=cel* eps    eps=s*sigma
          // s and cel must be symmetric, thus if s is not - symmetrize it by averaging off diagonal elements
          s=0.5*(s+s.Transpose());
 //# 1GPa=1e+9Pa=1e+9J/m^3
 //# 1meV= 1.60218e-22 J
 //# 1 A= 1e-10 m
 //# 1meV/pVol=1.60218e-22 J/A^3 x  A^3/pVol = 1.60218e+8  J/m^3 x  A^3/pVol = 1.60218e-1 GPa x  A^3/pVol
 //cf meV/Primitive Unit Cell Volume = 1  GPa
 double cf=phon.cs.pVol()*10.0/1.60218;
          a.Cel=cf*s.Inverse();  
          for(int h=1;h<=6;++h)   
          for(int k=1;k<=6;++k){if(a.Cel(h,k)<-1){fprintf(stderr,"elastic constants negative - rerun with stricter limits in mcphas.ini\n");exit(EXIT_FAILURE);}
                               if(fabs(a.Cel(h,k))<1){a.Cel(h,k)=0;}
                               }
 }
// generate magnetoelastic interactions 
// strategy: - take full unreduced interactions and set for all magnetic ions Gcfph=0
//           - determine phonon induced Gcfph^alphagamma(i)
//           - and multipolar Jgammagamma'(ij)by calculating the energy
//             for ui=0 eps=0 and comparing it with result obtained with 
//  a) zero eps_alpha=0 Ogamma(i)=fixed relax selfconsistently ui ---> Jgammagamma(ii) 
//  b)  nonzero eps_alpha=fixed Ogamma(i)=fixed relax ui and compare to energy of a) ----> Gcfph^alphagamma(i)        
//  c) zero eps nonzero Ogamma(i)=fixed Ogamma'(j)=fixed relax ui --->Jgammagamma'(ij)

// How to fix Ogamma(i) in a mean field loop ?
// ... use jjjpar - module set module = fix
 
for(int n=pw.cs.nofatoms;n>0;--n)
{h1=0;(*pw.jjj[n]).Icalc_parameter_storage_init(h1,Happ,T); // initialize eigenstate matrix
 if(!((*pw.jjj[n]).module_type==external_class&&
      (*pw.jjj[n]).pcalc(u0,  T,  Hxc,Hext,(*phon.jjj[n]).Icalc_parstorage)))
  {(*pw.jjj[n]).module_type=fixmom;(*pw.jjj[n]).MF=0; // for magnetic atoms use module fixmom
   for(int al=1;al<=6;++al)for(int g=1;g<=pw.cs.nofcomponents;++g)(*(*pw.jjj[n]).G)(al,g)=0;   // and remove Gcfph not to disturb the following a) b) c)
  }
}

          spincf sps(1,1,1,pw.cs.nofatoms,pw.cs.nofcomponents);
          mfcf mf(1,1,1,pw.cs.nofatoms,pw.cs.nofcomponents);
         
     

     
ini.doeps=0; // set epsilon=0


//******************************** TESTESTESTESTEST *********************************
/*
// test elastic constants for some strains
spincf spsp(1,1,1,phon.cs.nofatoms,phon.cs.nofcomponents);
mfcf mfp(1,1,1,phon.cs.nofatoms,phon.cs.nofcomponents);
 spsp.epsilon=0;ini.doeps=-1;
if(fecalc(U0,Eelastic,r,spinchange,Happ,T,ini,phon,spsp,mfp)>FEMIN_INI)
{fprintf(stderr,"Error reduce_unitcell option delatoms phon: reference energy not stable - modify mcphas.ini and restart\n");exit(EXIT_FAILURE);}
if(verbose){fprintf(stderr,"doeps=-1 U0=%12.12g",U0*phon.cs.nofatoms); } 
         
   for(int al=1;al<=6;++al){ini.doeps=-1; // ini.doeps=-1 will preserve the strain 
                            spsp.epsilon(al)=0.00350592;double Eel=0.5*spsp.epsilon(al)*a.Cel(al,al)*spsp.epsilon(al);
                            if(fecalc(Ub,Eelastic,r,spinchange,Happ,T,ini,phon,spsp,mfp)>FEMIN_INI)
        {fprintf(stderr,"Error reduce_unitcell option -delatoms phon: free energy not stable for %i  - modify mcphas.ini and restart\n",al);exit(EXIT_FAILURE);}
                            printf("Ub-U0=%12.12g  Eel=%12.12g \n",(Ub-U0)*phon.cs.nofatoms,Eel);      
if(al==4){spsp.print(stdout); myPrintVector(spsp.epsilon,"Strain x=1-6");}                       
                            spsp.epsilon(al)=0;
                            ini.doeps=0;
                           }
exit(0);
*/
//******************************** TESTESTESTESTEST *********************************


 sps.epsilon=0;  


if(fecalc(U0,Eelastic,r,spinchange,Happ,T,ini,pw,sps,mf)>FEMIN_INI)
{fprintf(stderr,"Error reduce_unitcell option delatoms phon: reference energy not stable - modify mcphas.ini and restart\n");exit(EXIT_FAILURE);}
if(verbose){fprintf(stderr,"U0=%12.12g",U0);}


// ---------------------------------
// 2 a) zero eps_alpha=0 Ogamma(i)=fixed relax selfconsistently ui ---> Jgammagamma(ii)
// ---------------------------------
Vector dnull(1,3);dnull=0; int nd;
for(int n=pw.cs.nofatoms;n>0;--n)// go through all magnetic ions
 if((*pw.jjj[n]).module_type==fixmom)
  for(int g=1;g<=pw.cs.nofcomponents;++g)
 {(*pw.jjj[n]).MF(g)=1;
       if(fecalc(Ua,Eelastic,r,spinchange,Happ,T,ini,pw,sps,mf)>FEMIN_INI)
        {fprintf(stderr,"Error reduce_unitcell option delatoms phon: free energy not stable for %i - modify mcphas.ini and restart\n",g);exit(EXIT_FAILURE);}
  if((nd=(*a.jjj[n]).index(dnull))==0)nd=(*a.jjj[n]).addpar(dnull,dnull,n);
   (*a.jjj[n]).jij[nd](g,g)+=-2.0*(Ua-U0)*pw.cs.nofatoms;  // nofatoms multiplied because U and f are normalised to meV/atom
//if(g==4){fprintf(stderr,"n=%i nd=%i g=%i Ua=%12.12g U0=%12.12g --> Jii44=%g\n",n,nd,g,Ua,U0,-2.0*(Ua-U0)*pw.cs.nofatoms);sps.print(stderr);}
// ---------------------------------
// 2 b)  nonzero eps_alpha=fixed Ogamma(i)=fixed relax ui and compare to
// ---------------------------------
//      energy of a) ----> Gcfph^alphagamma(i)   
   sps.epsilon=0;     
   for(int al=1;al<=6;++al){ini.doeps=-1; // ini.doeps=-1 will preserve the strain 
                            sps.epsilon(al)=1e-6;double Eelrenormdiveps=0.5*sps.epsilon(al)*a.Cel(al,al);
                            if(fecalc(Ub,Eelastic,r,spinchange,Happ,T,ini,pw,sps,mf)>FEMIN_INI)
        {fprintf(stderr,"Error reduce_unitcell option delatoms phon: free energy not stable for %i %i - modify mcphas.ini and restart\n",g,al);exit(EXIT_FAILURE);}
                              (*(*a.jjj[n]).G)(al,g)+=-(Ub-Ua)*pw.cs.nofatoms/sps.epsilon(al)+Eelrenormdiveps;
                            sps.epsilon(al)=0;
                            ini.doeps=0;
                           }

  (*pw.jjj[n]).MF(g)=0;
 }

// ---------------------------------
// 2 c) zero eps nonzero Ogamma(i)=fixed Ogamma'(j)=fixed relax ui --->Jgammagamma'(ij)
// ---------------------------------
Vector dabc(1,3),drijk(1,3); 
ini.doeps=0;
for(int n=pw.cs.nofatoms;n>0;--n)// go through all magnetic ions
 if((*pw.jjj[n]).module_type==fixmom)
  for(int n1=pw.cs.nofatoms;n1>0;--n1)// go through all magnetic ions
   if((*pw.jjj[n1]).module_type==fixmom)
    for(int g=1;g<=pw.cs.nofcomponents;++g)
     for(int g1=1;g1<=pw.cs.nofcomponents;++g1)
     if(n1!=n||g1!=g)
     {(*pw.jjj[n]).MF(g)=1;
      (*pw.jjj[n1]).MF(g1)=1;
   if(fecalc(Uc,Eelastic,r,spinchange,Happ,T,ini,pw,sps,mf)>FEMIN_INI)
        {fprintf(stderr,"Error reduce_unitcell option delatoms phon: free energy not stable for %i %i %i %i- modify mcphas.ini and restart\n",g,g1,n,n1);exit(EXIT_FAILURE);}
        dabc=(*a.jjj[n1]).xyz-(*a.jjj[n]).xyz;
       if((nd=(*a.jjj[n]).index(dabc))==0){dadbdc2ijk(drijk,dabc,pw.cs.abc);
                                          nd=(*a.jjj[n]).addpar(dabc,drijk,n1);}
(*a.jjj[n]).jij[nd](g,g1)+=-(Uc-U0)*pw.cs.nofatoms
                           -(*a.jjj[n]).jij[(*a.jjj[n]).index(dnull)](g,g)/2
                           -(*a.jjj[n1]).jij[(*a.jjj[n1]).index(dnull)](g1,g1)/2;
if (fabs((*a.jjj[n]).jij[nd](g,g1))<SMALL){(*a.jjj[n]).jij[nd](g,g1)=0;}
      (*pw.jjj[n]).MF(g)=0;
      (*pw.jjj[n1]).MF(g1)=0;
    }
    

// remove all phonons from a before outputting it ...
for(int n=a.cs.nofatoms;n>0;--n)
{h1=0;(*a.jjj[n]).Icalc_parameter_storage_init(h1,Happ,T); // initialize eigenstate matrix
 if(((*a.jjj[n]).module_type==external_class&&(*a.jjj[n]).pcalc(u0,  T,  Hxc,Hext,(*phon.jjj[n]).Icalc_parstorage)))
  {(*a.jjj[n]).module_type=fixmom;// trick to avoid that Cel, interactions are changed when deleting phonon degree of freedom
    a.delatom(-n,dis);}
}

}

/**********************************************************************/
// hauptprogramm
/**********************************************************************/
int main (int argc, char **argv)
{ 
// check command line
  if (argc <= 1)
    { printf (" program reduce_unitcell, output is written to stdout\n \
                use as: reduce_unitcell [option] mcphas.j\n\n \
                This program checks every atom in the unit cell in file mcphas.j and removes\n \
                any atom, which is connected to another by a lattice vector.\n \
                a list of superfluous sipf file is stored in reduce_unitcell_sipf.del \n \
                Options: -nofcomponents 23 fixes the nofcomponents to 23 by \n \
                        reducing (removing entries) or increasing (by filling with zeroes) \n \
                        the exchange parameter tables\n \
                        -i  forces output with indexchange \n \
                        -ni  forces output without indexchange \n \
			-delatoms 1,2,5,7    instead of removing atoms connected by a lattice vector \n \
			remove atoms number 1,2,5 and 7 from the list and also all interactions with those. \n \
                        In case the atoms to be removed have the phonon module, elastic constants are renormalized \n \
                         effective interactions based on the Einstein model are introduced between the remaining atoms,  \n \
                        crystal field phonon interactions with 4f shells can be treated. Mind: the Einstein model is of limited use \n \
			-delatoms 1:3+-0.75:4+-0.25,2,5,7  removes atoms number 1,2,5,7 for atom 1 the interactions\n \
			 are transferred to 75%% to atom 3 and 25%% to atom 4 and interactions to atoms 3 and 4 are\n \
			 removed,   for 2 5 7 all the interactions with other atoms are removed\n \
                         a +-- separator  triggers in the -delatoms mode that charges are transferred resetting\n \
                          CHARGE variablein sipf files (sipf files are rewritten with modified charge)\n   \
                        -delatoms phon\n   \
                         remove all atoms which have a single ion module with the \n   \
                         pcalc function (phonon degrees of freedom), renormalize elastic\n   \
                         constants by appying different components of the stress tensor and computing\n   \
                         self consistent strain, furthermore compute phonon induced magnetoelastic and \n   \
                         quadrupolar interactions by self consistent mean field calculations for various \n   \
                         strains and quadrupolar moments - this approach is to be preferred over the Einstein model \n   \                        
                        -mcdiff   create also mcdiff.in file with reduced unit cell\n \
                        -v  verbose mode\n \
                \n");
      exit (1);
    } else { fprintf (stderr,"#* reduce_unitcell 250123 *\n"); }

int ow=1,i=0; int n=0,noindexchange=0,mcdiff=0;bool delphon=false;
char * token;
char *substr[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];
float ns[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];ns[0]=MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL;
float nscoeff[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL+1];nscoeff[0]=MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL;

 while(argv[ow][0]=='-'){
 if(strcmp(argv[ow],"-nofcomponents")==0){ow+=1;
 // option setting nofcomponents
 n=(int)strtod(argv[ow],NULL);
 if(n<1){fprintf(stderr,"Error program add option nofcomponents=%i is less than 1\n",n);exit(1);}
                                        }
if(strcmp(argv[ow],"-delatoms")==0){ow+=1;
 if(strcmp(argv[ow],"phon")==0){delphon=true;}
 else{
 // substitute all commma with spaces
  while ((token=strrchr(argv[ow],','))!=NULL){++i;substr[i]=token+1;*token='\0';}
  ++i;substr[i]=argv[ow];
     }
}

 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-i")==0){noindexchange=-1;}
 if(strcmp(argv[ow],"-v")==0){verbose=1;}
 if(strcmp(argv[ow],"-mcdiff")==0){mcdiff=1;}
 ++ow;}

 par a(argv[ow]);

if(n>0){a.set_nofcomponents(n);if(verbose){fprintf(stderr,"Setting nofcomponents=%i\n",n);}}

if(delphon){
delphonons(a);
}
else
if(i==0){
 a.reduce_unitcell(verbose);  
}else{
for(n=1;n<=i;++n)
  {// substitute all : with ' ' in substr[n]
 while ((token=strchr(substr[n],':'))!=NULL){*token=' ';}
//fprintf(stderr,"string is %s\n",substr[n]);
// use splitstring to read atom numbers and as "error+-" the nscoeff
 int nn=splitstring(substr[n],ns,nscoeff); // nn-1= number of neighbours onto which to be distributed
 int dim=1;if(nn>1){dim=nn-1;}
   Matrix dis(1,dim,1,2);int an=ns[1];
 for(int ii=2;ii<=nn;++ii){
// do not use next line, but take care of already deleted atoms from previous loops (renumbering of atoms) ...
//dis(ii-1,1)=ns[ii];
int nnew=ns[ii];
for(int ic=1;ic<n;++ic)
 {int ndone=(int)strtod(substr[ic],NULL);if(ndone<ns[ii])--nnew;
 }

dis(ii-1,1)=nnew;
dis(ii-1,2)=nscoeff[ii];
                          }
   if(n<i){int anp1=strtod (substr[n+1], NULL);if((int)an<=(int)anp1){fprintf(stderr,"Error program reduce_unitcell option -delatoms atom numbers have to be given in ascending order\n");exit(1);}}
if(nn>1){ int ian=(int)an;  a.delatom(ian,dis,verbose); // distribute onto neighbours
}
else  
{// do not distribute onto neighbours.
 //In case of phonon module put effective multipolar interaction between other magnetic ions
 // and renormalise elastic constants
 // (to trigger this make ian negative in call to delatoms)
int ian=-(int)an; a.delatom(ian,dis,verbose);}  
           
  }
}

 a.save(stdout,noindexchange);
 if(mcdiff==1){a.save_mcdiff_in("reduce_unitcell");fprintf(stderr,"# created mcdiff.in\n");}

fprintf(stderr,"# end of reduce_unitcell - list of redundant sipf files\n");
fprintf(stderr,"# in file reduce_unitcell_sipf.del, to delete these files use:\n");
fprintf(stderr,"# perl -l -n -e \"unlink\" reduce_unitcell_sipf.del\n");
}


