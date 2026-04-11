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
// for statistics 
int isfull=0;

#include "myev.h"
#include "mcphas_htcalc.c"
#include "mcphas_fecalc.c"
#include "mcphas_physpropcalc.c"


void getU(double & U,Vector & Happ,double & T,
           inipar & ini,par & p,qvectors & testqs,testspincf & testspins,physproperties & physprop,const char * info)
{int j=1;int nofreppoints=ini.nofreppoints,nofconvrep=ini.nofconvrep;
 float rr=fmodf(ini.repeat-0.00001,1.0);
          double maxstamf=ini.maxstamf;int rep;
          int maxnofmfloops=ini.maxnofmfloops;
          double maxspinchange=ini.maxspinchange;
          for(rep=0;rep<=floor(ini.repeat)&&j>0;++rep)
          {j=htcalc(Happ,T,ini,p,testqs,testspins,physprop,0);
           // returns j=0 if successfull
 //  --> if no spinconfiguration has been found at ht point
 // returns j=1 if recalculation of fe yields too different value
 // returns j=2 if prevailing problem is maxnofmfloops reached
 // returns j=3 if prevailing problem is maxspinchange is reached  
          if(rep+1<=floor(ini.repeat)){
           switch (j)
           {case 1: ini.maxstamf*=rr;printf("repeating with maxstamf=%g\n",ini.maxstamf);
                    break;
            case 2: ini.maxnofmfloops/=rr;printf("repeating with maxnofmfloops=%i\n",ini.maxnofmfloops);
                    break;
            case 3: ini.maxspinchange/=rr;printf("repeating with maxspinchange=%g \n",ini.maxspinchange);
                    break;
            default:  ;
           }
                                       }
          } if(rep>1){++nofreppoints;ini.nofreppoints=nofreppoints;
                      if(j==0){++nofconvrep;ini.nofconvrep=nofconvrep;}
                     }
         ini.maxspinchange=maxspinchange;
         ini.maxnofmfloops=maxnofmfloops;
         ini.maxstamf=maxstamf;
         
  if(j>0){fprintf(stderr,"Error reduce_unitcell: self consistent MF calculation not converged for %s  - modify reduce_unitcell.ini and restart\n",info);exit(EXIT_FAILURE);}
 U=physprop.u;
}
// ***************************************************************************
// routine for the case of removing phononic degrees of freedom and
// - renormalising elastic constants
// - calculating (phonon induced) magnetoelastic and quadrupolar interactions
// ... interactions will be calculated for the first nprim atoms in unit celll par a 
// ***************************************************************************
void delphonons(par & a,bool symmetrize, int noindexchange,int nprim)
{fprintf(stderr,"# deleting phonons");
// 1. renormalise elastic constants using only phonon degrees of freedom
//   and calling htcalc with doeps for various applied external stresses ...
double U; int r; 
if (nprim>a.cs.nofatoms){fprintf(stderr,"Error reduce_unitcell - delphonons: nprim=%i>nofatoms=%i\n",nprim,a.cs.nofatoms);exit(1); }
 par phon(a); // for phonon --> elastic constants
par pw(a);  // pw working parameterset ---> multipolar and magnetoelastic interactions
time_t curtime;
  struct tm *loctime; 

 phon.set_nofcomponents(3);
// remove nonphononic degrees of freedom from phon
Vector u0(1,3),Hxc(1,3),Hext(1,6);double T=1;
Vector Happ(1,HEXT_DIMENSION);Happ=0;
Matrix dis(1,1,1,1);dis=0;dis(1,1)=1;

Vector h1_phon(1,phon.cs.nofcomponents);
for(int n=phon.cs.nofatoms;n>0;--n)
{h1_phon=0;(*phon.jjj[n]).Icalc_parameter_storage_init(h1_phon,Happ,T); // initialize eigenstate matrix
 if(!((*phon.jjj[n]).module_type==external_class&&(*phon.jjj[n]).pcalc(u0,  T,  Hxc,Hext,(*phon.jjj[n]).Icalc_parstorage)))
  {phon.delatom(-n,dis);}
}

if(verbose){fprintf(stderr,"Creating file results/reduce_unitcell_phon_mcphas.j\n");
         phon.save("results/reduce_unitcell_phon_mcphas.j",noindexchange);
            fprintf(stderr,"Creating file results/reduce_unitcell_me_mcphas.j\n");
         pw.save("results/reduce_unitcell_me_mcphas.j",noindexchange);
            }

           
double U0,Ua,Ub,Uc; // reference energy
char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
  char outfilename[MAXNOFCHARINLINE];
 inipar ini("reduce_unitcell.ini",prefix,"reduce_unitcell");ini.doeps=true;
if(ini.nofrndtries<0){fprintf(stderr,"# Error reduce_unitcell - nofrndtries<0 - Monte Carlo calculations not possible\n");exit(1); }
 double z,u;
  Vector h1ext(1,HEXT_DIMENSION);h1ext=0;


// initialize output file of mf configurations
FILE * fout;
if(verbose){fprintf(stderr,"Creating file results/reduce_unitcell_phon_mcphas.mf\n");
           fout=fopen_errchk("results/reduce_unitcell_phon_mcphas.mf","w");
 ini.outcolset=false;
ini.defaultcolcode(1,0); // make T the column 1 in reduce_unitcell_phon.mf
for(int n=1;n<=6;++n)ini.defaultcolcode(n+1,12+n);// put in col 2-7 the stress tensor

fprintf(fout, "#output file of program reduce_unitcell ");
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.mf-->\n");
   phon.savelattice(fout);phon.saveatoms(fout);
   ini.savedemagtensor(fout);
   fprintf (fout, "#!show_abc_unitcell=1.0\n");
   fprintf (fout, "#!show_primitive_crystal_unitcell=1.0\n");
   fprintf (fout, "#!show_magnetic_unitcell=1.0\n");
   fprintf (fout, "#!show_atoms=1.0\n");
   fprintf (fout, "#!show_chargedensity=1.0\n");
   fprintf (fout, "#!spins_scale_moment=1.0\n");
   fprintf (fout, "#!scale_view_1=1.0 scale_view_2=1.0 scale_view_3=1.0\n");
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   ini.print_usrdefcolcodes(fout);
   fprintf (fout, " nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsii eps2=epsjj eps3=epskk eps4=2epsjk eps5=2epsik eps6=2epsij\n");
   fprintf (fout, "    #mf1(atom 1) mf1(atom 2) .... selfconsistent Mean field configuration \n"); 
   fprintf (fout, "    #mf2(atom 1) mf2(atom 2) .... UNITS: mf(atom i)=gJ*mu_B*hxc(atom i)[meV] \n"); 
   fprintf (fout, "    #mf3(atom 1) mf3(atom 2) ....         (i.e. divide by gJ and mu_B=0.05788meV/Tesla to get exchange field hxc[Tesla]\n");
   
            }
// ---------------------------------
fprintf(stderr,"# 1. computing elastic constants ....\n");
// ---------------------------------
// only keep phononic degrees of freedom in phon and
// apply various stresses 

float x=1,y=1;Vector M(1,3);Vector P(1,3); M=0;P=0;
         
 {
#ifdef _THREADS
if(NUM_THREADS>256){fprintf(stderr,"Error mcphas: too many threads required - change hardcode limit 256 in mcphas_htcalc.c line 69 and recompile\n");exit(EXIT_FAILURE);}
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) 
                    tin[ithread] = new htcalc_input(0,ithread,&phon);
#endif
 physproperties physprop_phon(ini.nofspincorrs,ini.maxnofhkls,phon.cs);
   testspincf testspins_phon (ini.maxnoftestspincf,"reduce_unitcell_phon.tst",outfilename,phon.cs.nofatoms,phon.cs.nofcomponents);
   Vector Imax_phon(1,phon.cs.nofatoms*phon.cs.nofcomponents);
   Vector Imom_phon(1,phon.cs.nofcomponents);
//determine saturation momentum (used for generation of qvectors testqs_phon)
if(verbose==1){fprintf(stderr,"# determine saturation momentum running singleion calculations for different fields\n");}
T=1.0;for(int l=1;l<=phon.cs.nofatoms;++l){//h1_phon=0;(*phon.jjj[l]).Icalc_parameter_storage_init(h1_phon,h1ext,T); // initialize eigenstate matrix
      for (int im=1;im<=phon.cs.nofcomponents;++im){h1ext=0;h1_phon=0;h1_phon(im)=20*MU_B; //just put some high field
                            (*phon.jjj[l]).Icalc(Imom_phon,T,h1_phon,h1ext,z,u,(*phon.jjj[l]).Icalc_parstorage);
                            Imax_phon(phon.cs.nofcomponents*(l-1)+im)=Imom_phon(im);
                           }
      }
qvectors testqs_phon (ini.qmin,ini.qmax,ini.deltaq,ini.maxqperiod,ini.maxnofspins,phon,Imax_phon,outfilename,verbose);
  ini.testspins=&testspins_phon;  ini.testqs=&testqs_phon;


    Matrix s(1,6,1,6); //myPrintVector(sps.epsilon,"Initial Strain");
          if(verbose){printf("Calculating Reference Energy U0\n");}
    getU(U0,Happ,T,ini,phon,testqs_phon,testspins_phon,physprop_phon,"elastic constants U0");
    if(verbose)fprintf(stderr,"U0=%g meV ",U0);
           double cel=0.01;  // fixed stress to apply in GPa
          for(int n=1;n<=6;++n){Happ(6+n)+=cel;
 getU(U,Happ,T,ini,phon,testqs_phon,testspins_phon,physprop_phon,"elastic constants");
 if(verbose){ini.print_usrdefcols(fout,x,y,T,Happ,phon.cs.abc,M,P,false);
            fprintf (fout, " %i %i %i ",
            physprop_phon.mf.n()*physprop_phon.mf.nofatoms,physprop_phon.mf.nofatoms,physprop_phon.mf.nofcomponents);
            fprintf(fout,"0 %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(physprop_phon.sps.epsilon(1)),myround(physprop_phon.sps.epsilon(2)),myround(physprop_phon.sps.epsilon(3)),myround(physprop_phon.sps.epsilon(4)),myround(physprop_phon.sps.epsilon(5)),myround(physprop_phon.sps.epsilon(6)));
            physprop_phon.mf.print(fout);fprintf(fout,"\n");}

                                Happ(6+n)-=cel;if(verbose)myPrintVector(stderr,physprop_phon.sps.epsilon,"Strain x=1-6");
 //                               printf("U+sigma.eps-U0=%12.12g U0=%12.12g -sigma.eps=%12.12g\n",(U-U0)*phon.cs.nofatoms+cel*sps.epsilon(n)*phon.cs.pVol()/1.60218e-1,U0*phon.cs.nofatoms,-cel*sps.epsilon(n)*phon.cs.pVol()/1.60218e-1);
//if(n==4)sps.print(stdout);
                                for(int m=1;m<=6;++m)s(n,m)=physprop_phon.sps.epsilon(m)/cel;
print_time_estimate_until_end((6-n)/n);
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
          for(int k=1;k<=6;++k){if(h==k&&h<4&&a.Cel(h,k)<-1){fprintf(stderr,"elastic constants negative - rerun with stricter limits in mcphas.ini\n");exit(EXIT_FAILURE);}
                                if(fabs(a.Cel(h,k))<1){a.Cel(h,k)=0;}
                               }
#ifdef _THREADS
for (int ithread=0; ithread<ini.nofthreads; ithread++) delete tin[ithread];
#endif

 }
if(verbose)fprintf(stderr,"\n#elastic constants generated\n");
if(verbose){fclose(fout);
            fprintf(stderr,"creating file results/reduce_unitcell_me_mcphas.mf\n");
            fout=fopen_errchk("results/reduce_unitcell_me_mcphas.mf","w");
           ini.defaultcolcode(1,19); // make x the column 1 in reduce_unitcell_phon.mf
           ini.defaultcolcode(2,20); // make y the column 2 in reduce_unitcell_phon.mf
           ini.defaultcolcode(3,0); // make T the column 3 in reduce_unitcell_phon.mf
           //for(int n=1;n<=6;++n)ini.defaultcolcode(n+1,12+n);// put in col 2-7 the stress tensor

              fprintf(fout, "#output file of program reduce_unitcell ");
               curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
               fprintf(fout,"#!<--mcphas.mcphas.mf-->\n");
               pw.savelattice(fout);pw.saveatoms(fout);
               ini.savedemagtensor(fout);
   fprintf (fout, "#!show_abc_unitcell=1.0\n");
   fprintf (fout, "#!show_primitive_crystal_unitcell=1.0\n");
   fprintf (fout, "#!show_magnetic_unitcell=1.0\n");
   fprintf (fout, "#!show_atoms=1.0\n");
   fprintf (fout, "#!show_chargedensity=1.0\n");
   fprintf (fout, "#!spins_scale_moment=1.0\n");
   fprintf (fout, "#!scale_view_1=1.0 scale_view_2=1.0 scale_view_3=1.0\n");
   fprintf (fout, "# x>0, y ... indicate which interaction operator Ialpha at which ion n is nonzero\n");
   fprintf (fout, "#          x,y=(n-1)*nofcomponents+alpha  (two ion q-interaction)\n");
   fprintf (fout, "#          with alpha=1,...,nofcomponents=%i and n=1,...,nofatoms=%i\n",pw.cs.nofcomponents,pw.cs.nofatoms);
   fprintf (fout, "# x=0 ... only one interaction operator Ialpha ion n is nonzero (selfenergy)\n");
   fprintf (fout, "# x<0, y ... strain epsilon_beta (bet=1...6) nonzero and interaction operator Ialpha at  ion n is nonzero\n");
   fprintf (fout, "#          x=-beta y=(n-1)*nofcomponents+alpha  (crystal field phonon interation Gcfph(epsilon)\n");
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   ini.print_usrdefcolcodes(fout);
   fprintf (fout, " nofspins nofatoms(in primitive basis) nofmeanfield-components errorcode(0=ok,1=failed) eps1=epsii eps2=epsjj eps3=epskk eps4=2epsjk eps5=2epsik eps6=2epsij\n");
   fprintf (fout, "    #mf1(atom 1) mf1(atom 2) .... selfconsistent Mean field configuration \n"); 
   fprintf (fout, "    #mf2(atom 1) mf2(atom 2) .... UNITS: mf(atom i)=gJ*mu_B*hxc(atom i)[meV] \n"); 
   fprintf (fout, "    #mf3(atom 1) mf3(atom 2) ....         (i.e. divide by gJ and mu_B=0.05788meV/Tesla to get exchange field hxc[Tesla]\n");
   
            
       }

// ---------------------------------
// strategy: - take full unreduced interactions and set for all magnetic ions Gcfph=0
//           - determine phonon induced Gcfph^alphagamma(i)
//           - and multipolar Jgammagamma'(ij)by calculating the energy
//             for ui=0 eps=0 and comparing it with result obtained with 
//  a) zero eps_alpha=0 Ogamma(i)=fixed relax selfconsistently ui ---> Jgammagamma(ii) 
//  b)  nonzero eps_alpha=fixed Ogamma(i)=fixed relax ui and compare to energy of a) ----> Gcfph^alphagamma(i)        
//  c) zero eps nonzero Ogamma(i)=fixed Ogamma'(j)=fixed relax ui --->Jgammagamma'(ij)

// How to fix Ogamma(i) in a mean field loop ?
// ... use jjjpar - module set module = fix
physproperties physprop(ini.nofspincorrs,ini.maxnofhkls,pw.cs);
// load testspinconfigurations (nooftstspinconfigurations,init-file,sav-file)
 testspincf testspins (ini.maxnoftestspincf,"reduce_unitcell.tst",outfilename,pw.cs.nofatoms,pw.cs.nofcomponents);
  Vector Imax(1,pw.cs.nofatoms*pw.cs.nofcomponents);
  Vector Imom(1,pw.cs.nofcomponents);
  Vector h1(1,pw.cs.nofcomponents);


for(int n=pw.cs.nofatoms;n>0;--n)
{h1=0;(*pw.jjj[n]).Icalc_parameter_storage_init(h1,Happ,T); // initialize eigenstate matrix
 for (int im=1;im<=pw.cs.nofcomponents;++im){h1ext=0;h1=0;h1(im)=20*MU_B; //just put some high field
                            (*pw.jjj[n]).Icalc(Imom,T,h1,h1ext,z,u,(*pw.jjj[n]).Icalc_parstorage);
                            Imax(pw.cs.nofcomponents*(n-1)+im)=Imom(im);
                           //if(verbose==1)printf("Imax(%i)=%g\n",pw.cs.nofcomponents*(l-1)+im,Imax(pw.cs.nofcomponents*(l-1)+im));
			   }
 if(!((*pw.jjj[n]).module_type==external_class&&
      (*pw.jjj[n]).pcalc(u0,  T,  Hxc,Hext,(*pw.jjj[n]).Icalc_parstorage)))
  {(*pw.jjj[n]).module_type=fixmom;(*pw.jjj[n]).MF=0; // for magnetic atoms use module fixmom
   // and remove Gcfph not to disturb the following a) b) c)
   for(int al=1;al<=6;++al)for(int g=1;g<=pw.cs.nofcomponents;++g)(*(*pw.jjj[n]).G)(al,g)=0;  
  // and remove Jij between magnetic atoms for the same reason
   for(int m=1;m<=(*pw.jjj[n]).paranz;++m){int sl=(*pw.jjj[n]).sublattice[m];
       if(!((*pw.jjj[sl]).module_type==external_class&&
      (*pw.jjj[sl]).pcalc(u0,  T,  Hxc,Hext,(*pw.jjj[sl]).Icalc_parstorage)))(*pw.jjj[sl]).delpar(m);
                                   }

  }
}
qvectors testqs (ini.qmin,ini.qmax,ini.deltaq,ini.maxqperiod,ini.maxnofspins,pw,Imax,outfilename,verbose);
#ifdef _THREADS
if(NUM_THREADS>256){fprintf(stderr,"Error mcphas: too many threads required - change hardcode limit 256 in mcphas_htcalc.c line 69 and recompile\n");exit(EXIT_FAILURE);}
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) 
                    tin[ithread] = new htcalc_input(0,ithread,&pw);
#endif

 ini.testspins=&testspins;  ini.testqs=&testqs;

     
ini.doeps=0; // set epsilon=0

 physprop.sps.epsilon=0;  

getU(U0,Happ,T,ini,pw,testqs,testspins,physprop,"U0-reference");
if(verbose){fprintf(stderr,"\n U0=%12.12g\n",U0);}
 
// ---------------------------------
fprintf(stderr,"# 2. computing multipolar self interaction J(0 0 0)\ and magnetoelastic interaction Gcfph\n");
fprintf(stderr,"#a) zero eps_alpha=0 Ogamma(i)=fixed relax selfconsistently ui ---> Jgammagamma(ii)\n");
fprintf(stderr,"#b) nonzero eps_alpha=fixed Ogamma(i)=fixed relax ui and compare to\n");
fprintf(stderr,"#      energy of a) ----> Gcfph^alphagamma(i) \n");

// ---------------------------------
Vector dnull(1,3);dnull=0; int nd;int nofptstodo=0,nofptsdone=0;
for(int n=pw.cs.nofatoms;n>0;--n)// go through all magnetic ions
 if((*pw.jjj[n]).module_type==fixmom)
  for(int g=1;g<=pw.cs.nofcomponents;++g){++nofptstodo;}

for(int n=pw.cs.nofatoms;n>0;--n)// go through all magnetic ions
 if((*pw.jjj[n]).module_type==fixmom)
  for(int g=1;g<=pw.cs.nofcomponents;++g)
 {(*pw.jjj[n]).MF(g)=1;
          if(verbose){printf("Calculating self interaction J%i(0) for atom %i(%i), i.e. for component %i \n",g,n,pw.cs.nofatoms,g);}
  
       getU(Ua,Happ,T,ini,pw,testqs,testspins,physprop,"selfenergy");
      if(verbose){x=0;y=(n-1)*pw.cs.nofcomponents+g;
             ini.print_usrdefcols(fout,x,y,T,Happ,phon.cs.abc,M,P,false);
            fprintf (fout, " %i %i %i ",
            physprop.mf.n()*physprop.mf.nofatoms,physprop.mf.nofatoms,physprop.mf.nofcomponents);
            fprintf(fout,"0 %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(physprop.sps.epsilon(1)),myround(physprop.sps.epsilon(2)),myround(physprop.sps.epsilon(3)),myround(physprop.sps.epsilon(4)),myround(physprop.sps.epsilon(5)),myround(physprop.sps.epsilon(6)));
            physprop.mf.print(fout);fprintf(fout,"\n");}

  if((nd=(*a.jjj[n]).index(dnull))==0)nd=(*a.jjj[n]).addpar(dnull,dnull,n);
 (*a.jjj[n]).jij[nd](g,g)+=-2.0*(Ua-U0)*pw.cs.nofatoms;  // nofatoms multiplied because U and f are normalised to meV/atom
if (fabs((*a.jjj[n]).jij[nd](g,g))<SMALL){(*a.jjj[n]).jij[nd](g,g)=0;} 
else if (verbose){fprintf(stderr,"\n atom %i I_%i  <--> x=%g y=%g ",n,g,x,y);} 

//if(g==4){fprintf(stderr,"n=%i nd=%i g=%i Ua=%12.12g U0=%12.12g --> Jii44=%g\n",n,nd,g,Ua,U0,-2.0*(Ua-U0)*pw.cs.nofatoms);sps.print(stderr);}
// ---------------------------------
//nonzero eps_alpha=fixed Ogamma(i)=fixed relax ui and compare to
//      energy of a) ----> Gcfph^alphagamma(i)   
// ---------------------------------
   physprop.sps.epsilon=0;     
   for(int al=1;al<=6;++al){
 if(verbose){printf("Calculating magnetoelastic interaction  for atom %i(%i), term Gcfph epsilon_%i I%i \n",n,pw.cs.nofatoms,al,g);}
  
                 ini.doeps=-1; // ini.doeps=-1 will preserve the strain 
                 physprop.sps.epsilon(al)=1e-6;
                 double Eelrenormdiveps=0.5*physprop.sps.epsilon(al)*a.Cel(al,al);
                getU(Ub,Happ,T,ini,pw,testqs,testspins,physprop,"Gcfph");
      if(verbose){x=-al;y=(n-1)*pw.cs.nofcomponents+g;
             ini.print_usrdefcols(fout,x,y,T,Happ,phon.cs.abc,M,P,false);
            fprintf (fout, " %i %i %i ",
            physprop.mf.n()*physprop.mf.nofatoms,physprop.mf.nofatoms,physprop.mf.nofcomponents);
            fprintf(fout,"0 %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(physprop.sps.epsilon(1)),myround(physprop.sps.epsilon(2)),myround(physprop.sps.epsilon(3)),myround(physprop.sps.epsilon(4)),myround(physprop.sps.epsilon(5)),myround(physprop.sps.epsilon(6)));
            physprop.mf.print(fout);fprintf(fout,"\n");}
 (*(*a.jjj[n]).G)(al,g)+=-(Ub-Ua)*pw.cs.nofatoms/physprop.sps.epsilon(al)+Eelrenormdiveps;
if (fabs((*(*a.jjj[n]).G)(al,g))<SMALL){(*(*a.jjj[n]).G)(al,g)=0;} 
else if (verbose){fprintf(stderr,"\nepsilon_%i - atom %i I_%i  <--> x=%g y=%g ",al,n,g,x,y);} 

/*if(al==1&&g==8){
fprintf(stderr,"n=%i  Ub=%12.12g  Ua=%12.12g U0=%12.12g --> G11=%g\n",
  n,Ub*pw.cs.nofatoms,Ua*pw.cs.nofatoms,U0*pw.cs.nofatoms,(*(*a.jjj[n]).G)(al,g));
physprop.sps.print(stderr);}
*/
                  physprop.sps.epsilon(al)=0;
                  ini.doeps=0;
                           }

  (*pw.jjj[n]).MF(g)=0;
 ++nofptsdone;--nofptstodo;print_time_estimate_until_end(nofptstodo/nofptsdone);
 }

// ---------------------------------
fprintf(stderr,"# 3. computing multipolar two ion interactions\n");
fprintf(stderr,"# zero eps nonzero Ogamma(i)=fixed Ogamma'(j)=fixed relax ui --->Jgammagamma'(ij)\n");
// ---------------------------------
Vector dabc(1,3),drijk(1,3); 
Matrix prim_unitcell_ijk(1,3,1,3);
prim_unitcell_ijk=pw.cs.prim_unitcell_ijk();
ini.doeps=0;nofptstodo=0;nofptsdone=0;
//for(int n=pw.cs.nofatoms;n>0;--n)// go through all magnetic ions
for(int n=nprim;n>0;--n)// go through all magnetic ions
 if((*pw.jjj[n]).module_type==fixmom)
  for(int n1=pw.cs.nofatoms;n1>0;--n1)// go through all magnetic ions
   if((*pw.jjj[n1]).module_type==fixmom)
    for(int g=1;g<=pw.cs.nofcomponents;++g)
     for(int g1=1;g1<=pw.cs.nofcomponents;++g1)
     if(n1!=n||g1!=g){++nofptstodo;}

for(int n=nprim;n>0;--n)// go through all magnetic ions
 if((*pw.jjj[n]).module_type==fixmom)
  for(int n1=pw.cs.nofatoms;n1>0;--n1)// go through all magnetic ions
   if((*pw.jjj[n1]).module_type==fixmom)
    for(int g=1;g<=pw.cs.nofcomponents;++g)
     for(int g1=1;g1<=pw.cs.nofcomponents;++g1)
     if(n1!=n||g1!=g)
     {if(verbose){printf("Calculating two ion interaction J%i%i(%i%i), i.e. for atom %i I%i - atom %i I%i \n",g,g1,n,n1,n,g,n1,g1);}
  
      (*pw.jjj[n]).MF(g)=1;
      (*pw.jjj[n1]).MF(g1)=1;
  getU(Uc,Happ,T,ini,pw,testqs,testspins,physprop,"bilinear interaction");
     if(verbose){x=(n-1)*pw.cs.nofcomponents+g;y=(n1-1)*pw.cs.nofcomponents+g1;
             ini.print_usrdefcols(fout,x,y,T,Happ,phon.cs.abc,M,P,false);
            fprintf (fout, " %i %i %i ",
            physprop.mf.n()*physprop.mf.nofatoms,physprop.mf.nofatoms,physprop.mf.nofcomponents);
            fprintf(fout,"0 %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g\n",myround(physprop.sps.epsilon(1)),myround(physprop.sps.epsilon(2)),myround(physprop.sps.epsilon(3)),myround(physprop.sps.epsilon(4)),myround(physprop.sps.epsilon(5)),myround(physprop.sps.epsilon(6)));
            physprop.mf.print(fout);fprintf(fout,"\n");}

   double djij=-(Uc-U0)*pw.cs.nofatoms
                        -(*a.jjj[n]).jij[(*a.jjj[n]).index(dnull)](g,g)/2
                      -(*a.jjj[n1]).jij[(*a.jjj[n1]).index(dnull)](g1,g1)/2;
if(!symmetrize)
{  dabc=(*a.jjj[n1]).xyz-(*a.jjj[n]).xyz;

// either  a) only one representative interaction 
  dadbdc2ijk(drijk,dabc,pw.cs.abc);
  if((nd=(*a.jjj[n]).index(dabc))==0){nd=(*a.jjj[n]).addpar(dabc,drijk,n1);}
(*a.jjj[n]).jij[nd](g,g1)+=djij;

if (fabs((*a.jjj[n]).jij[nd](g,g1))<SMALL){(*a.jjj[n]).jij[nd](g,g1)=0;} 
else if (verbose){fprintf(stderr,"\natom %i I_%i - atom %i I_%i  <--> x=%g y=%g ",n,g,n1,g1,x,y);} 
}
else
{
//  b) 
// here we should generate all nearest neighbours on sublattice n1 and distribute
// the effective interaction on them - this should then be 
// in line with the symmetry of the system
// probe +- one supercell ... calculate distance to neighbours and 
// multiplicity
Vector s(1,3);double rm=1e10;int mult=0;
for(int si=-1;si<=1;++si)
for(int sj=-1;sj<=1;++sj)
for(int sk=-1;sk<=1;++sk)
 {s(1)=si;s(2)=sj;s(3)=sk;
  dabc=(*a.jjj[n1]).xyz+s-(*a.jjj[n]).xyz;
  dadbdc2ijk(drijk,dabc,pw.cs.abc);
  r=Norm(drijk);
  if(fabs(rm-r)<SMALL){++mult;}
  else
  if(r<rm-SMALL)if(r>SMALL){rm=r;mult=1;}
  
 }
for(int si=-1;si<=1;++si)
for(int sj=-1;sj<=1;++sj)
for(int sk=-1;sk<=1;++sk)
 {s(1)=si;s(2)=sj;s(3)=sk;
  dabc=(*a.jjj[n1]).xyz+s-(*a.jjj[n]).xyz;
  dadbdc2ijk(drijk,dabc,pw.cs.abc);
  r=Norm(drijk);
  if(fabs(rm-r)<SMALL){
//-----------
       if((nd=(*a.jjj[n]).index(dabc))==0){nd=(*a.jjj[n]).addpar(dabc,drijk,n1);}
(*a.jjj[n]).jij[nd](g,g1)+=djij/mult;

if (fabs((*a.jjj[n]).jij[nd](g,g1))<SMALL){(*a.jjj[n]).jij[nd](g,g1)=0;}


//-----------
                       }

  }
}


      (*pw.jjj[n]).MF(g)=0;
      (*pw.jjj[n1]).MF(g1)=0;
  ++nofptsdone;--nofptstodo;print_time_estimate_until_end(nofptstodo/nofptsdone); 
 }
fclose(fout);    

// remove all phonons from a before outputting it ...
for(int n=a.cs.nofatoms;n>0;--n)
{h1=0;(*a.jjj[n]).Icalc_parameter_storage_init(h1,Happ,T); // initialize eigenstate matrix
 if(((*a.jjj[n]).module_type==external_class&&(*a.jjj[n]).pcalc(u0,  T,  Hxc,Hext,(*a.jjj[n]).Icalc_parstorage)))
  {(*a.jjj[n]).module_type=fixmom;// trick to avoid that Cel, interactions are changed when deleting phonon degree of freedom
    a.delatom(-n,dis);}
}

if(symmetrize)
{if(verbose){fprintf(stderr,"sort interaction parameters according to ascending distance ...\n");}

// sort exchange parameters ...
a.sort();
if(verbose){fprintf(stderr,"averaging equidistant interaction parameters ...\n");}

// average equidistant neighbour interactions to obey symmetry relations
for(int n=a.cs.nofatoms;n>0;--n)
{for(int i=1;i<=(*a.jjj[n]).paranz;++i)
 {Matrix jav(1,a.cs.nofcomponents,1,a.cs.nofcomponents);
  jav=(*a.jjj[n]).jij[i];int mult=1;
  for(int j=1;i+j<=(*a.jjj[n]).paranz&&
               fabs(Norm((*a.jjj[n]).dr[i+j])-Norm((*a.jjj[n]).dr[i+j-1]))<SMALL;++j)
   {jav+=(*a.jjj[n]).jij[i+j];++mult;}
  jav=(1.0/mult)*jav;
  for (int j=0;j<mult;++j)(*a.jjj[n]).jij[i+j]=jav;
 }
}
}
#ifdef _THREADS
for (int ithread=0; ithread<ini.nofthreads; ithread++) delete tin[ithread];
#endif

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
                         constants by applying different components of the stress tensor and computing\n   \
                         self consistent strain, furthermore compute phonon induced magnetoelastic and \n   \
                         quadrupolar interactions by self consistent mean field calculations for various \n   \
                         strains and quadrupolar moments - this approach is to be preferred over the Einstein model \n   \
                         this approach needs a file reduce_unitcell.ini and optional reduce_unitcell.tst with  \n   \
                         the files correspond in the format exactly to mcphas.ini and mcphas.tst \n \
                        -delatoms phons\n   \
                         same as -delatoms phon but symmetrize interactions by creating equivalent nearest\n \
                         neighbours and averaging calculated two ion interactions \n \
                        -delatoms phone 3 3 3 2 2 2 \n   \
                         same as -delatoms phon but as inital step extend  primitive unit cell \n \
                         to a 3x3x3 supercell - then perform the calculation and \n \
                         ouput the results only for the primitive subcell 2 2 2 \n \
                        -mcdiff   create also mcdiff.in file with reduced unit cell\n \
                        -v  verbose mode\n \
                        -vv very verbose mode\n \
                \n");
      exit (1);
    } else { fprintf (stderr,"#* reduce_unitcell 250123 *\n"); }

int ow=1,i=0; int n=0,noindexchange=0,mcdiff=0,n1=1,n2=1,n3=1,s1=1,s2=1,s3=1;bool delphon=false,symmetrize,extend=false;
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
 if(strcmp(argv[ow],"phon")==0){delphon=true;symmetrize=false;}
 else
 if(strcmp(argv[ow],"phons")==0){delphon=true;symmetrize=true;}
 else
 if(strcmp(argv[ow],"phone")==0){delphon=true;symmetrize=false;extend=true;
                                 ++ow; n1=atoi(argv[ow]);
                                 ++ow; n2=atoi(argv[ow]);
                                 ++ow; n3=atoi(argv[ow]);
                                 ++ow; s1=atoi(argv[ow]);
                                 ++ow; s2=atoi(argv[ow]);
                                 ++ow; s3=atoi(argv[ow]);
                                }
 else{
 // substitute all commma with spaces
  while ((token=strrchr(argv[ow],','))!=NULL){++i;substr[i]=token+1;*token='\0';}
  ++i;substr[i]=argv[ow];
     }
}

 if(strcmp(argv[ow],"-ni")==0){noindexchange=1;}
 if(strcmp(argv[ow],"-i")==0){noindexchange=-1;}
 if(strcmp(argv[ow],"-v")==0){verbose=1;}
 if(strcmp(argv[ow],"-vv")==0){verbose=2;}
 if(strcmp(argv[ow],"-mcdiff")==0){mcdiff=1;}
 ++ow;}

 par a(argv[ow]);

if(n>0){a.set_nofcomponents(n);if(verbose){fprintf(stderr,"Setting nofcomponents=%i\n",n);}}

if(delphon){Matrix p(1,3,1,3);int nprim=a.cs.nofatoms;p=a.cs.r;// remember primitive lattice
  if(extend){a.extend_unitcell(n1,n2,n3,s1,s2,s3);}
 delphonons(a,symmetrize,noindexchange,nprim);
  if(extend){a.cs.r=p;a.reduce_unitcell(verbose);} // go back to original primitive lattice 
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

 a.sort(); // sort output parameters according to ascending distance
 a.save(stdout,noindexchange);
 if(mcdiff==1){a.save_mcdiff_in("reduce_unitcell");fprintf(stderr,"# created mcdiff.in\n");}

fprintf(stderr,"# end of reduce_unitcell - list of redundant sipf files\n");
fprintf(stderr,"# in file reduce_unitcell_sipf.del, to delete these files use:\n");
fprintf(stderr,"# perl -l -n -e \"unlink\" reduce_unitcell_sipf.del\n");
}


