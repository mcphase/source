/****************************************************
 * spins - abalyse and display spinconfiguration at given htpoint
 * Author: Martin Rotter
 ****************************************************/
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <cerrno>
#include <cstdlib>
#include <string>
#include "../../version"
#include <par.hpp>
#include "spincf.hpp"
#include "martin.h"
#include "graphic_parameters.hpp"
#include "cryststruct.hpp"
#include "elements.c"

#include "spins_parse.c"
#include "densities_func.c"
#include "spins_q.c"

// function for bond length elongation int table mode -tL2
void L2_calc(double & tL2,int i0, int j0, int k0, int ii0,FILE * fout, char * outstr, char * tstr, 
              Vector & nmin, Vector & nmax, spincf & savmf,Matrix & pm_unitcell,
             double & T,Vector & h0, Vector & Hext,par & inputpars)
{Vector r0(1,3),r(1,3),p0(1,3),p(1,3),h(1,inputpars.cs.nofcomponents),hh(1,savmf.nofcomponents*savmf.nofatoms);
 double d; 
 // if module allows to calculate position shift of an atom - only then use this for output 
    if(true==(*inputpars.jjj[ii0]).pcalc(p0,T,h0,Hext,(*inputpars.jjj[ii0]).Icalc_parstorage))
     {

 r0=savmf.pos(i0,j0,k0,ii0, inputpars.cs); // position of ion
  p0+=dr(savmf.epsilon,r0); // add the strain dr using position Vector r and given the strain tensor epsilon in Voigt notation
 for(int i1=nmin(1)-1;i1<=nmax(1)+1;++i1) // go through different prim magnetic unit cells
 for(int j1=nmin(1)-1;j1<=nmax(1)+1;++j1)
 for(int k1=nmin(1)-1;k1<=nmax(1)+1;++k1)
    for (int i=1;i<=savmf.na();++i)
    for(int j=1;j<=savmf.nb();++j)
    for(int k=1;k<=savmf.nc();++k){hh=savmf.m(i,j,k);//look at each atom
      for(int ii=1;ii<=inputpars.cs.nofatoms;++ii)
   {r=savmf.pos(i,j,k,ii, inputpars.cs); 
    r+=(double)i1*pm_unitcell.Column(1)+(double)j1*pm_unitcell.Column(2)+(double)k1*pm_unitcell.Column(3);
    d=Norm(r-r0);if(d>0.1&&d<tL2){ h=0; 
   for(int nt=1;nt<=inputpars.cs.nofcomponents;++nt){h(nt)=hh(nt+inputpars.cs.nofcomponents*(ii-1));}
     if(true==(*inputpars.jjj[ii]).pcalc(p,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
     {p+=dr(savmf.epsilon,r);// add displacement due to strain
       // now p-p0 is the displacement and r-r0 the unstrained bond length
       // ... calculate the projection of p-p0 on r-r0 and square it
     double L2=(p-p0)*(r-r0); L2/=d;
       L2*=L2; // square this projection - this is the longitudinal bond elongation
                                  fprintf(fout,"%s %s %g %g\n",outstr,tstr,d,L2);
     }
                                 }
   }}
  }
fprintf(fout,"#\n"); // comment line after each atom
}
                  

         

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
bool load_structure(
    int argc,
    char** argv,
    LoadedInput& s)
{
  if (s.mode != OperationMode::Table&&s.fileMode != FileMode::FC){
  s.fout = fopen_errchk ("./results/spins.out", "w"); // unless it is table -t* or -fc options
   s.cs.print_mcdiff_in_header(s.fout,"spins",0);
fprintf(s.fout,"\
#   0.73250   0.00000                       -1.43200 -1.43200  0.71600  0.00000 -0.71600  1.43200  0.00000\n\
#   0.73250   0.00000                       -2.56800 -2.56800  1.28400  0.00000 -1.28400  2.56800  0.00000\n\
#\n\
#\n\
# %%SECTION 3%% DESCRIPTION OF THE LATTICE\n\
#\n\
# -----------------------------------------------------------------------------\n");

// input file header and conf------------------------------------------------------------------
   headerinput(s.fin,s.fout,s.gp,s.cs,s.out.data(),s.N);}
else
  {
// input file header and conf------------------------------------------------------------------
   headerinput(s.fin,stderr,s.gp,s.cs,s.out.data(),s.N);
  }

  
   
// if(cs.nofatoms<1){fclose (fin);fprintf(stderr,"#!!! Error program spins reading nofatoms=%i - must be >0 !!!\n",cs.nofatoms);exit(1);}
    if (s.cs.nofatoms < 1) {
        fprintf(stderr, "Error spins: no atoms found in input structure\n");
        return false;
    }
 //if(cs.nofcomponents<1){fclose (fin);fprintf(stderr,"#!!! Error program spins reading nofcomponents=%i - must be >0 !!!\n",cs.nofcomponents);exit(1);}

    if (s.cs.nofcomponents < 1) {
        fprintf(stderr, "Error spins: no components found in input structure\n");
        return false;
    }

   //spincf savmf(1,1,1,cs.nofatoms,cs.nofcomponents);
s.savmf.reset(new spincf(1, 1, 1, s.cs.nofatoms, s.cs.nofcomponents));

    // Move the remaining crystal-structure setup here, including
    // s.cs4 initialization and atom-position adjustments.
// ------------------------load spinsconfigurations and check which one is nearest -------------------------------   
for(int i=1;i<=NOF_USERDEF_MCPHAS_COLS;++i)s.aa[i]=1e100;
int apf=0;
for(int i=s.os+1;i<argc;++i){
if(strcmp(argv[i],"-prefix")==0){apf=2;}
                       }
double x,y; 
int i;
if (s.mode==OperationMode::File&&argc<3+s.os){
                 // in this case we have exactly one argument left,
                 // i.e.  s.aa[1] becomes a number of a spinconfig in a file
                 if(argc<2+s.os){fprintf(stderr,"Error spins - option -f / -fc requires at least one argument to identify spin configuration in file\n");exit(EXIT_FAILURE); }
                 s.aa[0]=strtod(argv[1+s.os],NULL);printf("# the configuration number %g\n",s.aa[0]);
                 if(mystrnstr(s.infilename, ".sps\0", MAXNOFCHARINLINE)!=NULL||
                    mystrnstr(s.infilename, ".mf\0", MAXNOFCHARINLINE)!=NULL)s.aa[0]=-s.aa[0];
                 }
                           
else{if(argc-1==2+s.os||argc-1==6+s.os+apf){
               // in this case  x and y in the phasediagram are given in the command line
               // arguments. We have to find out, which columns ix and iy
               //  correspond to x and y and put into corresponding s.aa[ix] and s.aa[iy] the
               // remaining arguments
                s.aa[0]=0;
                i=find_usrdef_out("x",s.out.data());if(i>0)s.aa[i]=strtod(argv[1+s.os],NULL);
                i=find_usrdef_out("y",s.out.data());if(i>0)s.aa[i]=strtod(argv[2+s.os],NULL);
               }
     else
     if(argc-1==4+s.os||argc-1==8+s.os+apf)
     {// in this case  T, Ha, Hb, Hc in the phasediagram are given in the command line
               // arguments. We have to find out, which columns iT iHa iHb  and iHc
               //  correspond to T Ha Hb and Hc and put into corresponding s.aa[iT],... and s.aa[iHc] the
               // remaining arguments
              s.aa[0]=0;
             i=find_usrdef_out("T",s.out.data());if(i>0)s.aa[i]=strtod(argv[1+s.os],NULL);
             i=find_usrdef_out("Ha",s.out.data());if(i>0)s.aa[i]=strtod(argv[2+s.os],NULL);
             i=find_usrdef_out("Hb",s.out.data());if(i>0)s.aa[i]=strtod(argv[3+s.os],NULL);
             i=find_usrdef_out("Hc",s.out.data());if(i>0)s.aa[i]=strtod(argv[4+s.os],NULL);
     }
     else
     if(argc-1==NOF_USERDEF_MCPHAS_COLS+s.os||argc-1==NOF_USERDEF_MCPHAS_COLS+4+s.os+apf)
     {s.aa[0]=0;for(int i=1;i<=NOF_USERDEF_MCPHAS_COLS;++i)s.aa[i]=strtod(argv[i+s.os],NULL);
     } else
     if(s.mode==OperationMode::Table&&argc==s.os+1+apf) // this is true for -t* options and no parameters - should print all
      {s.longtable=true;s.aa[0]=-1;s.nomorethanthis=1;}
     else
     {
      fprintf(stderr,"#Error program spins - wrong number of arguments %i (optional %i)!\n",argc-1,s.os);
      exit(1);
     }
    } 

if(check_for_best(s.fin,s.aa,*s.savmf,x,y,s.T,s.Hext,s.outstr,s.out.data(),s.inputpars->cs.abc,s.nomorethanthis))
  {fclose (s.fin);fprintf(stderr,"#!!! Error program spins - no stable structure found !!!\n");exit(1);}
if(!s.longtable){fclose (s.fin);} // for longtable, i.e. more HT point we leave s.fin open



if(s.nofcomp>s.savmf->nofcomponents){s.nofcomp=s.savmf->nofcomponents;}
if(s.nofcomp<s.savmf->nofcomponents){printf("#! printing only nofcomponents=%i components\n",s.nofcomp);}
  

if(s.fileMode == FileMode::FC){return true;}// loaded enough for -fc option

  
int ii,nt,k,j;
// determine primitive magnetic unit cell
Vector xyz(1,3),dd0(1,3),dd3(1,3),dd(1,3);
if (s.inputpars->cs.r!=s.cs.r){s.cs.r=s.inputpars->cs.r;}
if (s.inputpars->cs.abc!=s.cs.abc){s.cs.abc=s.inputpars->cs.abc;}

// check sipffilenames -> if these are not present in cs (from headerinput), put also atomic positions from s.inputpars
for(ii=1;ii<=s.inputpars->cs.nofatoms;++ii)
{if(s.cs.sipffilenames[ii]==NULL){  
 s.cs.sipffilenames[ii]=new char[MAXNOFCHARINLINE];
 }
 if (strcmp((*s.inputpars->jjj[ii]).sipffilename,s.cs.sipffilenames[ii])!=0)
             {strcpy(s.cs.sipffilenames[ii],(*s.inputpars->jjj[ii]).sipffilename);
              s.cs.x[ii]=s.inputpars->cs.x[ii];
              s.cs.y[ii]=s.inputpars->cs.y[ii];
              s.cs.z[ii]=s.inputpars->cs.z[ii];
// HERE take care about atoms sitting at nearly the same position (a nucleus and a mangetic shell
// from makenn -cfph ... and put the positions exactly at the same correct value 
// (makenn had shifted the magnetic charge cloud by 0.01 A along c in order to enable 
// the correct evaluation of cf-phonon interactions by mcphas and mcdisp. here we correct for
// this shift in order to get the right output of charge densities and in spins.out the
// right positions for doing mcdiff !!
                   for(int i=1;i<ii;++i){//loop all atoms which have been read
                     if ( fabs((s.cs.x[ii]-s.cs.x[i])*s.cs.abc[1])<0.2 &&
                          fabs((s.cs.y[ii]-s.cs.y[i])*s.cs.abc[2])<0.2 &&
                          fabs((s.cs.z[ii]-s.cs.z[i])*s.cs.abc[3])<0.2 )
                         { // atom i and ii are the same atom therefore check if they
                           // are displaced along c and move magnetic atom back to nuclear position
                          if(fabs((s.cs.x[ii]-s.cs.x[i])*s.cs.abc[1])>0.001 ||
                             fabs((s.cs.y[ii]-s.cs.y[i])*s.cs.abc[2])>0.001)
                              {fprintf(stderr,"Error spins.c: atoms %i (%g %g %g) and %i (%g %g %g) too close\n",i,s.cs.x[i],s.cs.y[i],s.cs.z[i],ii,s.cs.x[ii],s.cs.y[ii],s.cs.z[ii]);
                                      return false;}
                            if(s.cs.z[ii]>s.cs.z[i]){s.cs.z[ii]=s.cs.z[i];}
                            else{s.cs.z[i]=s.cs.z[ii];}
                                            }
                                       }
             }
}
s.cs4.abc=s.cs.abc;s.cs4.r=s.cs.r;s.cs4.nofatoms=s.cs.nofatoms;s.cs4.nofcomponents=s.cs.nofcomponents;
s.savmf->calc_prim_mag_unitcell(s.p,s.cs.abc,s.cs.r);
  
 if (s.fileMode == FileMode::F) // option -f all loading done which is needed
 { return true;}
   

// FROM HERE ON IT IS ONLY EXECUTED IF GRAPHICS ARE DESIRED OR -t OPTION IS PRESENT... 

s.gp.read();

  ii=0; 
  // if individual ions of the cluster are to be shown make nofatoms in spincf larger !
  for (j=1;j<=s.inputpars->cs.nofatoms;++j){
         if(s.arrow==4&&(*s.inputpars->jjj[j]).module_clust==true){for(k=1;k<=(*(*s.inputpars->jjj[j]).clusterpars).cs.nofatoms;++k)
                                                 {++ii;par inputpars4((*(*s.inputpars->jjj[j]).clusterpars));
                                                 s.cs4.x[ii]=s.cs.x[j]+(*inputpars4.jjj[k]).xyz(1);
                                                 s.cs4.y[ii]=s.cs.y[j]+(*inputpars4.jjj[k]).xyz(2);
                                                 s.cs4.z[ii]=s.cs.z[j]+(*inputpars4.jjj[k]).xyz(3);s.cs4.sipffilenames[ii]=new char[MAXNOFCHARINLINE];
                                                 strcpy(s.cs4.sipffilenames[ii],(*inputpars4.jjj[k]).sipffilename);
                                                 //check if cluster abc  are the same as s.inputpars
           if(fabs(s.cs4.abc(1)-inputpars4.cs.abc(1))>1e-5)
              {fprintf(stderr,"Error program spins - a=%g in cluster %s not the same as a=%g in mcphas.j\n",s.cs4.abc(1),(*s.inputpars->jjj[j]).sipffilename,inputpars4.cs.abc(1));exit(1);}
           if(fabs(s.cs4.abc(2)-inputpars4.cs.abc(2))>1e-5)
              {fprintf(stderr,"Error program spins - b=%g in cluster %s not the same as b=%g in mcphas.j\n",s.cs4.abc(2),(*s.inputpars->jjj[j]).sipffilename,inputpars4.cs.abc(2));exit(1);}
           if(fabs(s.cs4.abc(3)-inputpars4.cs.abc(3))>1e-5)
              {fprintf(stderr,"Error program spins - c=%g in cluster %s not the same as c=%g in mcphas.j\n",s.cs4.abc(3),(*s.inputpars->jjj[j]).sipffilename,inputpars4.cs.abc(3));exit(1);}
           if(fabs(s.cs4.abc(4)-inputpars4.cs.abc(4))>1e-5)
              {fprintf(stderr,"Error program spins - alpha=%g in cluster %s not the same as alpha=%g in mcphas.j\n",s.cs4.abc(4),(*s.inputpars->jjj[j]).sipffilename,inputpars4.cs.abc(4));exit(1);}
           if(fabs(s.cs4.abc(5)-inputpars4.cs.abc(5))>1e-5)
              {fprintf(stderr,"Error program spins - beta=%g in cluster %s not the same as beta=%g in mcphas.j\n",s.cs4.abc(5),(*s.inputpars->jjj[j]).sipffilename,inputpars4.cs.abc(5));exit(1);}
           if(fabs(s.cs4.abc(6)-inputpars4.cs.abc(6))>1e-5)
              {fprintf(stderr,"Error program spins - gamma=%g in cluster %s not the same as gamma=%g in mcphas.j\n",s.cs4.abc(6),(*s.inputpars->jjj[j]).sipffilename,inputpars4.cs.abc(6));exit(1);}                                                 
                                                 }
                                               }
                                               else
                                              {++ii;s.cs4.x[ii]=s.cs.x[j];s.cs4.y[ii]=s.cs.y[j];s.cs4.z[ii]=s.cs.z[j];s.cs4.sipffilenames[ii]=s.cs.sipffilenames[j];}
              } s.cs4.nofatoms=ii;

    return true;
}

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
/**********************************************************************/
// hauptprogramm
/**********************************************************************/
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
int main (int argc, char **argv)
{ double x,y;
  char  tstr [MAXNOFCHARINLINE];

LoadedInput LI;

    if (!parse_options(argc, argv, LI)) {
        if (LI.fin) fclose(LI.fin);
        if (LI.fout && LI.fout != stdout) fclose(LI.fout);
        return EXIT_FAILURE;
    }

switch (LI.mode) {
case OperationMode::Eigenvector:
             return spins_q(argc,argv,LI);
              break;
default:
  if (!load_structure(argc, argv, LI)) {
        if (LI.fin) fclose(LI.fin);
        if (LI.fout && LI.fout != stdout) fclose(LI.fout);
        return EXIT_FAILURE;
    }
 }
   // ----------------------------output configuration ----------------------------------------------------------------
    fprintf(stdout,"#! %s - configuration\n",LI.outstr);
    fprintf(stdout,"# Components of the Demagnetisation Tensor in SI Units\n");
    fprintf(stdout,"# refering to ijk coordinate system\n");
    fprintf(stdout,"# defined by  j||b, k||(a x b) and i normal to k and j\n");
    fprintf(stdout,"#!Nii=%g Nij=%g Nik=%g\n",LI.N(1,1),LI.N(1,2),LI.N(1,3));
    fprintf(stdout,"#!Njj=%g Njk=%g\n",LI.N(2,2),LI.N(2,3));
    fprintf(stdout,"#!Nkk=%g\n",LI.N(3,3));

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
//    Print   Header
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
switch (LI.mode) {
case OperationMode::File:
 switch (LI.fileMode) {
  case FileMode::FC:   // for option -fc limit nofat (atoms to be printed)
   {if(LI.nofat>LI.savmf->nofatoms){LI.nofat=LI.savmf->nofatoms;}
   if(LI.nofat<LI.savmf->nofatoms){printf("#! printing only nofatoms=%i atoms\n",LI.nofat);}
   if(strcmp(argv[LI.os]+strlen(argv[LI.os])-3,".mf")==0)
   {LI.savmf->print_commented(stdout,"Hex",LI.minl,LI.maxl,LI.maxn,LI.limit,1,LI.nofat);}
   else{LI.savmf->print_commented(stdout,"I",LI.minl,LI.maxl,LI.maxn,LI.limit,1,LI.nofat);}
        break; }
 // .............................................................................                                
 case FileMode::F:
   { LI.inputpars->savelattice(LI.fout);
   fprintf (LI.fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(LI.fout,"#! strain tensor: eps1=%4.4g (epsii) eps2=%4.4g (epsjj) eps3=%4.4g (epskk) eps4=%4.4g (2epsjk) eps5=%4.4g (2epsik) eps6=%4.4g (2epsij)\n",
    myround(LI.savmf->epsilon(1)),myround(LI.savmf->epsilon(2)),myround(LI.savmf->epsilon(3)),myround(LI.savmf->epsilon(4)),myround(LI.savmf->epsilon(5)),myround(LI.savmf->epsilon(6)));

  fprintf(LI.fout,"#! %s \n",LI.outstr);
   if(LI.T==0){fprintf(LI.fout,"# program spins: temperature not found in %s - setting T=1 K\n",argv[2]);LI.T=1;}
  fprintf(LI.fout,"#!%s: nr1=%i nr2=%i nr3=%i \n",LI.outstr,LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc());
  fprintf(LI.fout,"#1           2     3     4     5       6       7       8    9    10\n");
  fprintf(LI.fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ia> <Ib> <Ic> [created by program spins]\n");
  Vector dd3(1,3),dd0(1,3);
  for (int i=1;i<=LI.savmf->na();++i){for(int j=1;j<=LI.savmf->nb();++j){for(int k=1;k<=LI.savmf->nc();++k)
  {for(int ii=1;ii<=LI.inputpars->cs.nofatoms;++ii)
   {// output the positions
    dd3=LI.savmf->pos_dabc(i,j,k,ii, LI.cs);
   //returns position dd3 as components	    with respect to lattice a b c
   
   // dd0 is position as components with respect to primitive crystallographic unit cell
    dd0=LI.savmf->pos_dr123(i,j,k,ii, LI.cs);
    fprintf(LI.fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            (*LI.inputpars->jjj[ii]).sipffilename,dd3(1),dd3(2),dd3(3),dd0(1),dd0(2),dd0(3));
    //output the "magnetic" moment if possible ... actually it outputs Ia Ib Ic 
    for(int nt=1;nt<=3;++nt){fprintf(LI.fout," %4.4f",myround(1e-5,LI.savmf->m(i,j,k)(LI.inputpars->cs.nofcomponents*(ii-1)+nt)));}      
    fprintf(LI.fout,"\n");
   }
  }}}

  fclose(LI.fout);
    LI.gp.spins_scale_moment=1;LI.gp.phonon_scale_static_displacements=1;
    LI.gp.read();
    LI.fin = fopen_errchk ("./results/spins_prim.jvx", "w"); // here draw some graphics for the spinconfiguration
                       // with option -f. Yet this is very limited, because we do not know what is the magnetic
                       // moment, i.e. a guess is made using the components 1 2 3 of the LI.savmf 
     LI.gp.showprim=1;
     LI.gp.scale_density_vectors=0;LI.gp.show_density=0; // do not show charge densities
             Vector hkl1(1,3);hkl1=0;
             Vector gjmbHxc1(1,3);gjmbHxc1=0;
     spincf magmom(LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc(),LI.savmf->nofatoms,3); // the magnetic moment guess
     int i,j,k,l;for (i=1;i<=LI.savmf->na();++i){for (j=1;j<=LI.savmf->nb();++j){for (k=1;k<=LI.savmf->nc();++k){for(l=1;l<=magmom.nofatoms;++l)  
         for(int momdim=1;momdim<=3&&momdim<=LI.savmf->nofcomponents;++momdim)
               {magmom.moment(i,j,k,l)(momdim)=LI.savmf->moment(i,j,k,l)(momdim);}
        }}}

     LI.savmf->jvx_cd(LI.fin,LI.outstr,LI.cs,
                  LI.gp,0.0,(*LI.savmf)*0.0,(*LI.savmf)*0.0,
                  hkl1,LI.T,gjmbHxc1,LI.Hext,LI.cs,magmom,magmom * 0.0,magmom* 0.0,magmom * 0.0);
    fclose (LI.fin);
     break; }
   default: break;
                     }
   return 0;  // we are done with File mode -f -fc exit !
   break;
// .............................................................................                                
case OperationMode::Graphical:
   LI.savmf->print(stdout,LI.nofcomp);
   // no break because also the next lines should be printed
// .............................................................................                                
case OperationMode::Table:
   {fprintf (LI.fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(LI.fout,"#! strain tensor: eps1=%4.4g (epsii) eps2=%4.4g (epsjj) eps3=%4.4g (epskk) eps4=%4.4g (2epsjk) eps5=%4.4g (2epsik) eps6=%4.4g (2epsij)\n",
   myround(LI.savmf->epsilon(1)),myround(LI.savmf->epsilon(2)),myround(LI.savmf->epsilon(3)),myround(LI.savmf->epsilon(4)),myround(LI.savmf->epsilon(5)),myround(LI.savmf->epsilon(6)));
   break;}
default: break;
}

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
//             Column Headers 
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

switch (LI.mode) {
 case OperationMode::Graphical:
 // the following is for the printout of spins.out ...........................
 section4header(LI.fout);
 fprintf(LI.fout,"#! %s : ",LI.outstr);
 fprintf(LI.fout," nr1=%i nr2=%i nr3=%i nat=%i  atoms in primitive magnetic unit cell\n",LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc(),LI.cs4.nofatoms*LI.savmf->na()*LI.savmf->nb()*LI.savmf->nc());
 //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
 //                                        to Sa Sb Sc La Lb Lc
 fprintf(LI.fout,"#1           2     3     4     5       6       7       8    9    10                  11   12   13   14   15   16\n");
 fprintf(LI.fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)\n");
 fprintf(LI.fout,"#          corresponding exchange fields [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)\n");
	break;
// .............................................................................                                
case OperationMode::Table: 
  {fprintf(stdout,"#! nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc(),LI.cs4.nofatoms*LI.savmf->na()*LI.savmf->nb()*LI.savmf->nc());
   tstr[0]='\0';for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col)snprintf(tstr,MAXNOFCHARINLINE,"%s %s",tstr,LI.out.data()[col]);
 snprintf(tstr,MAXNOFCHARINLINE,"%s n nr1 nr2 nr3 {sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] ",tstr);
  switch (LI.tableMode) 
        {
   case TableMode::TMSL:
     snprintf(tstr,MAXNOFCHARINLINE,"%s <Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)]\n",tstr);
     break;
   case TableMode::TI:
     snprintf(tstr,MAXNOFCHARINLINE,"%s <I1> <I2> <I3> ... <Inofcomponents>\n",tstr);
     break;
   case TableMode::THex:
     snprintf(tstr,MAXNOFCHARINLINE,"%s <Hex1> <Hex2> <Hex3> ... <Inofcomponents> (meV)\n",tstr);
     break;
   case TableMode::TL2:
     snprintf(tstr,MAXNOFCHARINLINE,"%s bondlength_to_neighbour [A]  elongation^2 [A^2]\n",tstr);
     nlimits_calc(LI.nmin,LI.nmax, LI.tL2,LI.p);
   default: break ;
        }
    printf("# ");print_col_numbers(stdout,tstr);
   printf("\n#%s",tstr);
   break;}
  default: break; 
 }       


//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
//    M E A N F I E L D - P R O C E S S I N G - loop all atoms
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

spincf spinconf(LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc(),LI.cs4.nofatoms,3);
spincf sc_phonon(LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc(),LI.cs4.nofatoms,3);

//  1. from the meanfieldconfiguration (*LI.savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: Icalc has to be used to calculate all the <Olm>.
 
Vector h(1,LI.inputpars->cs.nofcomponents);
Vector hh(1, LI.savmf->nofcomponents * LI.savmf->nofatoms);
Vector I(1,LI.inputpars->cs.nofcomponents);
spincf densitycf(
    LI.savmf->na(),
    LI.savmf->nb(),
    LI.savmf->nc(),
    LI.savmf->nofatoms,
    LI.dim);
h=0;for(int ii=1;ii<=LI.inputpars->cs.nofatoms;++ii)
{(*LI.inputpars->jjj[ii]).Icalc_parameter_storage_init(h,LI.Hext,LI.T);} // initialize Icalc module parameter storage


bool moretoprint=true;
while(moretoprint){
 for (int i=1;i<=LI.savmf->na();++i)for(int j=1;j<=LI.savmf->nb();++j)for(int k=1;k<=LI.savmf->nc();++k)
 {
    hh=LI.savmf->m(i,j,k);
  densitycf.m(i,j,k)=0;int i4=1,ii4=1;
  for(int ii=1;ii<=LI.inputpars->cs.nofatoms;++ii)
 {  bool doTableLine=(LI.nofat==MAXNOFATOMS||LI.nofat==ii); 
     // if nofat was initialized - print only table line for this atom
    Vector magmom(1,3),mom(1,3);
    Vector Lmom(1,3);
    Vector Smom(1,3);
    Vector moments(1,LI.dim);
      Vector momS(1,SPINDENS_EV_DIM);
  Vector momL(1,ORBMOMDENS_EV_DIM);
  Vector momentsx(1,SPINDENS_EV_DIM);
  Vector momentsy(1,SPINDENS_EV_DIM);
  Vector momentsz(1,SPINDENS_EV_DIM);
  Vector momentlx(1,ORBMOMDENS_EV_DIM);
  Vector momently(1,ORBMOMDENS_EV_DIM);
  Vector momentlz(1,ORBMOMDENS_EV_DIM);
    h=0;
   for(int nt=1;nt<=LI.inputpars->cs.nofcomponents;++nt){h(nt)=hh(nt+LI.inputpars->cs.nofcomponents*(ii-1));}


switch(LI.arrow)
{case 1: (*LI.inputpars->jjj[ii]).Scalc(mom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage);break;
 case 2: (*LI.inputpars->jjj[ii]).Lcalc(mom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage);break;
 case 3: (*LI.inputpars->jjj[ii]).mcalc(mom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage);break;
 case 5: (*LI.inputpars->jjj[ii]).pelcalc(mom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage);break;
 default: break;
}

switch(LI.arrow)
{case 1:
 case 2:
 case 3:
 case 5:
         for(int nt=1;nt<=3;++nt)
		        {spinconf.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };break;
 case 4: int dim4;
         dim4=3;if((*LI.inputpars->jjj[ii]).module_clust==true)
                           dim4=(*(*LI.inputpars->jjj[ii]).clusterpars).cs.nofatoms*3;
         Vector momi(1,dim4);
         (*LI.inputpars->jjj[ii]).micalc(momi,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage);
         // now put the components of momi to the spinconf                                               
         for(int nt=1;nt<=dim4;++nt){spinconf.m(i,j,k)(nt+i4-1)=momi(nt);};i4+=dim4;
         break;
}
   Vector dd0(1,3),dd3(1,3);
    
// output atoms and moments in primitive unit cell to LI.fout  ------------------------------------
if(LI.arrow==4&&(*LI.inputpars->jjj[ii]).module_clust==true){
    for(int nt=1;nt<=(*(*LI.inputpars->jjj[ii]).clusterpars).cs.nofatoms;++nt)
     {dd3=spinconf.pos_dabc(i,j,k,ii4, LI.cs4);
      dd0=spinconf.pos_dr123(i,j,k,ii4, LI.cs4);
      if (LI.mode == OperationMode::Table&&doTableLine){
              //fprintf(LI.fout," %s ", LI.outstr);
       double value;
      for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col)
       {extract(LI.outstr,LI.out.data()[col],value); fprintf(LI.fout,"%g ",myround(value));}
      fprintf(LI.fout," %i %i %i %i ",LI.savmf->n(),LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc());
                                  
                                     }
      fprintf(LI.fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            LI.cs4.sipffilenames[ii4],dd3(1),dd3(2),dd3(3),dd0(1),dd0(2),dd0(3));
      if (LI.mode != OperationMode::Table||(LI.tableMode==TableMode::TMSL&&doTableLine))
         {
            fprintf(LI.fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(1+3*(ii4-1))));
            fprintf(LI.fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(2+3*(ii4-1))));
            fprintf(LI.fout," %4.4f\n",myround(1e-5,spinconf.m(i,j,k)(3+3*(ii4-1))));
          } 
       // do not print out L S or exchange fields for cluster module
      ++ii4;
     }
    }
else{   ++ii4;
    // output the positions
    dd3=LI.savmf->pos(i,j,k,ii, LI.cs);
      // if module allows to calculate position shift of an atom - use this for output 
    if(true==(*LI.inputpars->jjj[ii]).pcalc(mom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage))
     {dd3+=mom; 
     //  fprintf(stderr,"# Attention: atom %i shifted from equilibrium position by (%g %g %g) A \n",ii,mom(1),mom(2),mom(3));
     }
     // if module cannot calculate position shift of an atom - check if the atom 
     // is actually a charge cloud on top of another atom and if yes take position shift
     // from this other atom
      else
     {for(int iii=1;iii<=LI.inputpars->cs.nofatoms;++iii)
      {if(ii!=iii)
       {Vector ddc(1,3); ddc=LI.savmf->pos(i,j,k,iii, LI.cs);
        if(abs(ddc-dd3)<0.01){
Vector hhh(1,LI.inputpars->cs.nofcomponents); hhh=0;
for(int nt=1;nt<=LI.inputpars->cs.nofcomponents;++nt){hhh(nt)=hh(nt+LI.inputpars->cs.nofcomponents*(iii-1));}

 if(true==(*LI.inputpars->jjj[iii]).pcalc(mom,LI.T,hhh,LI.Hext,(*LI.inputpars->jjj[iii]).Icalc_parstorage))
    {dd3+=mom;
 //printf(stderr,"# Attention: atom %i shifted from equilibrium position by (%g %g %g) A taken from atom %i\n",ii,mom(1),mom(2),mom(3),iii);
 }
                               }
       }
      } 
     }

Vector dd(1,3);
    dd0=LI.p.Inverse()*dd3;
dd0(1)*=LI.savmf->na();dd0(2)*=LI.savmf->nb();dd0(3)*=LI.savmf->nc();
    Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,LI.cs.abc);
    dd=abc_in_ijk.Inverse()*dd3;
 if(doTableLine)
 {  if(LI.mode==OperationMode::Table)
    {snprintf(tstr,MAXNOFCHARINLINE,"{%s} %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f ",
            LI.cs.sipffilenames[ii],dd(1),dd(2),dd(3),dd0(1),dd0(2),dd0(3));
        
        double value;
      for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col)
       {extract(LI.outstr,LI.out.data()[col],value); fprintf(LI.fout,"%g ",myround(value));}
      fprintf(LI.fout," %i %i %i %i ",LI.savmf->n(),LI.savmf->na(),LI.savmf->nb(),LI.savmf->nc());
      fprintf(LI.fout,"%s",tstr);
       } 
    
   switch(LI.tableMode) { 
     case TableMode::TMSL :
     case TableMode::None :
       { 
     //ouput the magnetic moment if possible
     if((*LI.inputpars->jjj[ii]).mcalc(magmom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage))
     {     
      for(int nt=1;nt<=3;++nt){fprintf(LI.fout," %4.4f",myround(1e-5,magmom(nt)));}
      // and output the orbital and spin momentum if possible 
      if((*LI.inputpars->jjj[ii]).Lcalc(Lmom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage)&&
        (*LI.inputpars->jjj[ii]).Scalc(Smom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage))
       {        for(int nt=1;nt<=3;++nt){fprintf(LI.fout," %4.4f",myround(1e-5,Smom(nt)));}
               for(int nt=1;nt<=3;++nt){fprintf(LI.fout," %4.4f",myround(1e-5,Lmom(nt)));}
       }
      }
     break;}
     default: break;
    }
   switch(LI.tableMode) {
     case TableMode::THex :
     case TableMode::None :
      {
    // finally output a line with the exchange fields 
    if (LI.mode != OperationMode::Table)fprintf(LI.fout,"\n                 corresponding exchange fields [meV]-->          ");
                      for(int nt=1;nt<=LI.savmf->nofcomponents;++nt)  // printout exchangefields
                        {fprintf(LI.fout," %4.4f",myround(1e-5,h(nt)));}
                        
        break; }
      
     case TableMode::TI:
     {double lnZ;double U;
     (*LI.inputpars->jjj[ii]).Icalc(I,LI.T,h,LI.Hext,lnZ,U,(*LI.inputpars->jjj[ii]).Icalc_parstorage);
                         for(int nt=1;nt<=LI.savmf->nofcomponents;++nt)  // printout I operator expectation values
                        {fprintf(LI.fout," %4.4f",myround(1e-5,I(nt)));}
                         
      break;}
     case TableMode::TL2:
      {// here calculate up to tl2 Angstroem bondlengths and elongation^2 and output to LI.fout
    L2_calc(LI.tL2,i,j,k,ii,LI.fout,LI.outstr,tstr,LI.nmin,LI.nmax,*LI.savmf,LI.p,LI.T,h,LI.Hext,*LI.inputpars);
      break;}
     default: break;
    }

 if(LI.tableMode!=TableMode::TL2)fprintf(LI.fout,"\n");
}
// -----------------------------------------------------------------------------------------------
   

if(LI.density){
switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c':  (*LI.inputpars->jjj[ii]).chargedensity_coeff (moments,LI.T, h, LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage); break;
 case 's':  if(LI.xx!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).spindensity_coeff (momentsx,1,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.yy!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).spindensity_coeff (momentsy,2,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.zz!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).spindensity_coeff (momentsz,3,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.doijk==3){ moments=LI.xx*momentsx+LI.yy*momentsy+LI.zz*momentsz;}
            else{for(int i=1;i<=SPINDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'o':  if(LI.xx!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentsx,1,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.yy!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentsy,2,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.zz!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentsz,3,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.doijk==3){ moments=LI.xx*momentsx+LI.yy*momentsy+LI.zz*momentsz;}
            else{for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+ORBMOMDENS_EV_DIM)=momentsy(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'm':  if(LI.xx!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).spindensity_coeff (momentsx,1,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.yy!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).spindensity_coeff (momentsy,2,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.zz!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).spindensity_coeff (momentsz,3,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            momS=LI.xx*momentsx+LI.yy*momentsy+LI.zz*momentsz;
            if(LI.xx!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentlx,1,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.yy!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momently,2,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            if(LI.zz!=0||LI.doijk<3)(*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentlz,3,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            momL=LI.xx*momentlx+LI.yy*momently+LI.zz*momentlz;
            for(int i=1;i<=SPINDENS_EV_DIM;++i){
            if(LI.doijk==3){moments(i)=momS(i);moments(i+SPINDENS_EV_DIM)=momL(i);}
            else{moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);
                 moments(i+3*SPINDENS_EV_DIM)=momentlx(i);moments(i+4*SPINDENS_EV_DIM)=momently(i);moments(i+5*SPINDENS_EV_DIM)=momentlz(i);
                }
                                               }
            break;
 case 'j':  (*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentlx,1,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            (*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momently,2,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            (*LI.inputpars->jjj[ii]).orbmomdensity_coeff (momentlz,3,LI.T, h,LI.Hext, (*LI.inputpars->jjj[ii]).Icalc_parstorage);
            for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){
             moments(i)=momentlx(i);moments(i+ORBMOMDENS_EV_DIM)=momently(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentlz(i);
             }
            break;
 default: help();return EXIT_FAILURE;
}
                   for(int nt=1;nt<=LI.dim;++nt)
		        {densitycf.m(i,j,k)(nt+LI.dim*(ii-1))=moments(nt);
                    }
} // LI.gp.show_density


if(LI.phonon==1)  // if module allows to calculate position  - use this for graphics ...
{if(true==(*LI.inputpars->jjj[ii]).pcalc(mom,LI.T,h,LI.Hext,(*LI.inputpars->jjj[ii]).Icalc_parstorage))
{for(int nt=1;nt<=3;++nt)
		        {sc_phonon.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };
}
 
}
} // if cluster or atoms - else
} // for ii
} // for ijk

if(LI.longtable){if(check_for_best(LI.fin,LI.aa,*LI.savmf,x,y,LI.T,LI.Hext,LI.outstr,LI.out.data(),LI.inputpars->cs.abc,1))
              {fclose(LI.fin);return 0;}
                printf("#  next HT\n");
                }
else {
      moretoprint=false;}
} // while moretoprint
if (LI.mode==OperationMode::Table){exit(0);}
  fclose (LI.fout);
   
// create plot of spinconfiguration -----------------------------------------------------------
printf("# ************************************************************************\n");
printf("#%s\n",LI.gp.title);
printf("# ************************************************************************\n");
              if(LI.arrow==0)LI.gp.spins_scale_moment=0;
              if(LI.density==0)LI.gp.show_density=0;
if(LI.eps==1){
    LI.fin = fopen_errchk ("./results/spins.eps", "w");
     spinconf.eps(LI.fin,LI.outstr);
    fclose (LI.fin);

// here the 3d file should be created
    LI.fin = fopen_errchk ("./results/spinsab.eps", "w");


     spinconf.eps3d(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,1,spinconf);
    fclose (LI.fin);
    LI.fin = fopen_errchk ("./results/spinsac.eps", "w");
     spinconf.eps3d(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,2,spinconf);
    fclose (LI.fin);
    LI.fin = fopen_errchk ("./results/spinsbc.eps", "w");
     spinconf.eps3d(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,3,spinconf);
    fclose (LI.fin);
    LI.fin = fopen_errchk ("./results/spins3dab.eps", "w");
     spinconf.eps3d(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,4,spinconf);
    fclose (LI.fin);
    LI.fin = fopen_errchk ("./results/spins3dac.eps", "w");
     spinconf.eps3d(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,5,spinconf);
    fclose (LI.fin);
    LI.fin = fopen_errchk ("./results/spins3dbc.eps", "w");
     spinconf.eps3d(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,6,spinconf);
    fclose (LI.fin);
          }
if(LI.fst==1){
    LI.fin = fopen_errchk ("./results/spins.fst", "w");
     spinconf.fst(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,spinconf);
    fclose (LI.fin);

    
   LI.fin = fopen_errchk ("./results/spins_prim.fst", "w");
     spinconf.fstprim(LI.fin,LI.outstr,LI.cs4.abc,LI.cs4.r,LI.cs4.x,LI.cs4.y,LI.cs4.z,spinconf);
    fclose (LI.fin);
          }
const int status =jvx_create(LI, argc, argv, densitycf, spinconf, sc_phonon, h, hh, x, y);
 for (int i = 1; i <= LI.cs.nofatoms; ++i) {
    delete[] LI.cs.sipffilenames[i];
    LI.cs.sipffilenames[i] = nullptr;
  }
 return status;

}



