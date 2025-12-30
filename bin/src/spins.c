/****************************************************
 * spins - abalyse and display spinconfiguration at given htpoint
 * Author: Martin Rotter
 ****************************************************/
#include "../../version"
#include <par.hpp>
#include "spincf.hpp"
#include "martin.h"
#include "graphic_parameters.hpp"
#include "cryststruct.hpp"
#include "densities_func.c"
#include "elements.c"

void help_and_exit()
    { printf ("\n\
program spins - popout spin/exchange field configuration\n"
"              - and/or display 3d animation of spin/moment/densities and animations\n\n\
use as: spins -f[c 1 13 3 0.1] [-n 2] mcphas.sps T Ha Hb Hc\n\
    or: spins -f[c 1 13 3 0.2] [-n 2] mcphas.sps x y\n\
    or: spins -f[c 1 13 3 0.2] [-n 2] mcphas.sps out1 out2 out3 out4 out5 out6 out7\n\
    or: spins -f[c 1 13 3 0.1] [-n 2] mcphas.tst n\n\
    or: spins -tMSL [-prefix 001] T Ha Hb Hc \n\
    or: spins -tHex [-prefix 001]  T Ha Hb Hc \n\
    or: spins -tI  [-prefix 001] T Ha Hb Hc \n\
    or: spins -tL2 4.5 [-prefix 001] T Ha Hb Hc \n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M|-pel] [-P] [-eps|-fst] [-prefix 001] T Ha Hb Hc [h k l E]\n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M|-pel] [-P] [-eps|-fst] [-prefix 001] x y [h k l E] \n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M|-pel] [-P] [-eps|-fst] [-prefix 001] out1 out2 out3 out4 out5 out6 out7 [[-prefix 020]h k l E]\n\
                    \n\
1) if used with -f file T Ha Hb Hc, this file has to be a mcphas.mf or mcphas.sps file,\n\
   the spin configuration at given temperature T[K] and magnetic effective field H[T]\n\
   is read and extracted from this file and printed on screen (stdout),\n\
   results/spins.out is created (with mag moment chosen to be = <Ia> <Ib> <Ic>),\n\
   a simple graphics to represent the configuration is created in results/spins_prim.jvx \n\
   Note: - T Ha Hb Hc stand for the 3rd 5th 6th 7th column in mcphas.* output files\n\
         (may have different meaning if out3 out5 out6 out7 is set in mcphas.ini)\n\
         - out1,...,out7 refers to column 1-7 in mcphas.* output files\n\
  \n\
2) if used with -f file x y, then this file has to be a mcphas.mf or mcphas.sps file,\n\
   the spin configuration at a given x,y point is read and extracted from this file ,\n\
   and printed on screen (stdout) etc. as 1)\n\
   Note: x y stands for the 1st and 2nd column in mcphas.* output files\n\
         (may have different meaning if out1 or out2 is set in mcphas.ini)\n\
3) if used with -f filen n,  this file has to be a mcphas.tst file,\n\
   the spin configuration number n\n\
   is read and extracted from this file and printed on screen (stdout),\n\
   results/spins.out is created (with mag moment chosen to be = <Ia> <Ib> <Ic>)\n\
   a simple graphics to represent the configuration is created in results/spins_prim.jvx \n\
  \n\
1&2&3) - if used with -fc min max n lim a human readable format is output for spin components with index \n\
   from min to max, only n numbers exceeding  \n\
   absolute value of lim, e.g. -fc 1 13 3 0.1 outputs I1,I2,...,I13, however only at maximum 3 numbers \n\
   for each atom which are all larger (absolute value) than 0.1. \n\
      - if option -n 5 ist present only I1,..,I5 are output (mind 5<=nofcomponents)\n\
 \n\
4) if used without a filename, the information is read from results/mcphas.* results/mcdisp.*\n\
   output files and tables or 3d graphical animations are created.\n\
   for table the options are:\n\
       -tMSL ... output to stdout a table with T Ha Hb Hc atom positions and with \n\
              magnetic moments <Mx> <My> <Mz>, orbital moments and spin \n\
              of each atom in the magnetic unitc cell \n\
       -tHex ...  output to stdout table with T Ha Hb Hc atom positions and exchange fields Hex\n\
       -tI   ... a similar table with expectation values of interaction operators <I>\n\
       -tL2 4.5  ... a table with average squared longitudinal bond elongations up to bond\n\
                     length 4.5 Angstroem, only atoms with phonon degrees of freedom are considered\n\
                     (useful for estimating elastic energy contributions of different bonds)\n\
   for graphical animations the options are: \n\
         -c ... calculate chargedensity\n\
         -s ... calculate spindensity\n\
         -o ... calculate angular orbital momentum density\n\
         -m ... calculate magnetic moment density\n\
         -j ... calculate currentdensity\n\
         -p i j k ... calculate projection of spin/orbital/current/magnetic moment density\n\
                  along direction i j k, e.g. 0 0 1\n\
         -div    ... calculate divergence of spin/orbital/current/magnetic moment density  \n\
         -S  ... show arrow indicating spin\n\
         -L  ... show arrow indicating orbital angular momentum\n\
         -M  ... show arrow indicating magnetic moment (for cluster show total moment)\n\
         -Mi ... show arrow indicating magnetic moment (for cluster show individual moments)\n\
         -P  ... calculate phononic displacement\n\
         -pel  ... show arrow indicating electric dipole moment\n\
         -eps ... create eps (postscript) files in addition to jvx (javaview) output\n\
         -fst ... create fst (fullprof viewer) files in addition to jvx (javaview) output\n\
\n\
         note, that in order to animate changes in the above quantities, the corresponding\n\
         switch has to be enabled in the mcdisp calculation (mcdisp.par) and the single ion\n\
         modules have to be capable of calculating the corresponding observables. \n\
\n\
         -prefix 001 ... use input file(s) results/001mc* instead of results/mc*\n\
 \n\
     example:\n\
        spins -c 2 0 0 1\n\
        ...calculates the charge density at T=2K and H=(0,0,1) Tesla\n\
        spins -tI 2 0 0 1 \n\
        ... outputs a table with atomic positions and expectations values <I> \n\
\n\
 This program outputs a magnetic structure (and magnetic excitation)\n\
 graphic/movie in the output files of different format:\n\
 results/spins*.eps (postscript), results/spins*.fst (fp_studio), \n\
 results/spins.out (ascii) and results/spins*.jvx (javaview)\n\n\
 the graphics output format can be fine tuned in .sps and .qev input files\n"
" or results/graphics_parameters.set by show_abc_unitcell,\n"
" show_primitive_crystal_unitcell, spins_scale_moment, spins_wave_amplitude\n\
 show_magnetic_unitcell, show_atoms, scale_view_1,scale_view_2, scale_view_3 ...\n\n\
 jvx files can be viewed by:\n\
 java javaview results/spins.jvx \n\
 java javaview \"model=results/spins.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n"
" gif images stored by javaview can be connected to animated gif by ImageMagick \n"
" magick -delay 1 -size 100x100 -loop 1 geomAnim.*.gif output.gif\n"
" ppm images stored by javaveiw can be connected to mpeg video by ppmtompeg (part of netpbm package)\n"
" ppmtompeg param \n param is a file with parameters - here an example:\n"
"OUTPUT movie.mpeg\n\
INPUT_DIR .\n\
INPUT\n\
*.ppm [001-021]\n\
END_INPUT\n\
BASE_FILE_FORMAT PNM\n\
INPUT_CONVERT *\n\
FRAME_RATE 25\n\
PATTERN IBBPBBPBBPBBPBB\n\
SLICES_PER_FRAME 16\n\
GOP_SIZE 30\n\
PIXEL HALF\n\
IQSCALE 1\n\
PQSCALE 5\n\
BQSCALE 10\n\
RANGE 5\n\
PSEARCH_ALG TWOLEVEL\n\
BSEARCH_ALG CROSS2\n\
REFERENCE_FRAME DECODED\n"
);
 exit (1);
    }


void section4header (FILE * fout)
{
fprintf(fout,"\
#\n\
# %%SECTION 4%% DESCRIPTION OF MAGNETIC UNIT CELL AND LIST OF MAGNETIC ATOMS\n\
#\n\
#\n\
# here follows the description of the magnetic unit cell with respect\n\
# to the primitive crystallographic unit cell:\n\
# 'nr1', 'nr2', 'nr3' ...the crystallographic unit cell has to be taken\n\
#                        nr1 nr2 and nr3 times along r1 r2 and r3,\n\
#                        respectively to get magnetic unit cell\n\
# 'nat' denotes the number of magnetic atoms in magnetic unit cell\n\
#\n\
# It follows a list of nat lines with to describe the magnetic moment configuration\n\
# Notes:\n\
# 'atom-filename' means the single ion property filename of this magnetic atom:\n\
#                 -it must contain the Formfactor Coefficients (e.g. see international tables)\n\
#                                      Lande factor\n\
#                                      Neutron Scattering Length (10^-12 cm) \n\
#                 -it may contain a    Debey Waller Factor\n\
# 'da' 'db' and 'dc' are not used by the program (unless you enter a line #! use_dadbdc=1)\n\
# 'dr1','dr2' and 'dr3' refer to the primitive lattice given below\n\
# 'Ma','Mb','Mc' denote the magnetic moment components in Bohr magnetons\n\
#                in case of non orthogonal lattices instead of Ma Mb Mc the components Mx My Mz\n\
#                have to be given, which refer to an right handed orthogonal coordinate system \n\
#                defined by y||b, z||(a x b) and x normal to y and z\n\
#  <Sa> <Sb> <Sc>  <La> <Lb > <Lc>  (optional) denote the spin and orbital angular momentum components \n\
# 'mf1' 'mf2' 'mf3' (optional line, used to go beyond dipole approx for formfactor)\n\
#                                     denote the corresponding exchange fields in meV\n\
#\n");
}



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
          
/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ 
fprintf(stderr,"# ***********************************************************\n");
fprintf(stderr,"# * spins - analyse mcphas output and display 3d graphics of*\n");
fprintf(stderr,"# *  spins,moments,densities, etc at given H and T          *\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R              *\n");
fprintf(stderr,"# * %s                                      *\n",MCPHASVERSION);
fprintf(stderr,"# ***********************************************************\n");

 FILE * fin, * fout;
 int i,minl=1,maxl=1,eps=0,fst=0;
 cryststruct cs,cs4;
 
 char outstr[MAXNOFCHARINLINE];
 char tstr[MAXNOFCHARINLINE];
 char infilename[MAXNOFCHARINLINE];
 char mcphasjfilename[MAXNOFCHARINLINE];

 char prefix[MAXNOFCHARINLINE];prefix[0]='\0';
 Vector nmin(1,3),nmax(1,3);
 Matrix N(1,3,1,3); // demagnetisation factor
  int dim=28,nofcomp=1000000000;
 int os=0,maxn=0; int doijk=0,arrow=0,density=0,phonon=0;//,arrowdim=3;
 double xx=0,yy=0,zz=0,limit=0,tL2=0;
graphic_parameters gp;
gp.show_abc_unitcell=1.0;
gp.show_primitive_crystal_unitcell=1.0;
gp.show_magnetic_unitcell=1.0;
gp.show_atoms=1.0;
gp.scale_view_1=1.0;
gp.scale_view_2=1.0;
gp.scale_view_3=1.0;
gp.spins_scale_moment=0;
gp.show_density=0;
snprintf(gp.title,MAXNOFCHARINLINE,"output of program spins");

 // check command line
 if (argc < 2){help_and_exit();}
// first: options without graphics just screendump <I> or exchange field configuration at given HT
 if (strncmp(argv[1],"-f",2)==0)
 { os=2;if (strcmp(argv[1],"-fc")==0){os=6;minl=(int)strtod(argv[2],NULL);maxl=(int)strtod(argv[3],NULL);
                                      maxn=(int)strtod(argv[4],NULL);limit=strtod(argv[5],NULL);
                                     }
   if (strncmp(argv[os],"-n",2)==0){nofcomp=(int)strtod(argv[os+1],NULL);os+=2;}
   fin = fopen_errchk (argv[os], "rb");printf("#* program spins ... reading from file %s\n",argv[os]);   
 }
 else { if (strncmp(argv[1],"-t",2)==0){os=1;fout=stdout;
                          if (strncmp(argv[1],"-tL2",4)==0){os=2;tL2=strtod(argv[2],NULL);}
                                       }
       else  // second ... other options with graphics !!
 {if(strcmp(argv[1],"-c")==0){os=1;}
  if(strcmp(argv[1],"-s")==0){os=1;}
  if(strcmp(argv[1],"-o")==0){os=1;} 
  if(strcmp(argv[1],"-m")==0){os=1;}
  if(strcmp(argv[1],"-j")==0){os=1;}
  if(os==1)
  {gp.show_density=1;density=1;

switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c': dim=CHARGEDENS_EV_DIM;
printf("#chargedensity is expanded in tesseral harmonics Zlm\n\
#   ro(r) sum_lm (a(l,m) R^2(r) Zlm(Omega)\n\
#   M. Rotter et al. J Phys: Conf Ser. 325 (2011) 012005\n#\n ");
 snprintf(gp.title,MAXNOFCHARINLINE,"chargedensity ro(r)");
 gp.threshhold=-0.05;
           break;
 case 's': dim=SPINDENS_EV_DIM;
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
  if(doijk==3) snprintf(gp.title,MAXNOFCHARINLINE,"projection of spindensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,MAXNOFCHARINLINE,"divergence of spindensity div Ms(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,MAXNOFCHARINLINE,"abs value  of spindensity |Ms(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'o': dim=ORBMOMDENS_EV_DIM;
printf("#orbital momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm  aL(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
  if(doijk==3) snprintf(gp.title,MAXNOFCHARINLINE,"projection of orbmomdensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,MAXNOFCHARINLINE,"divergence of orbmomdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,MAXNOFCHARINLINE,"abs value  of orbmomdensity |ML(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'm': dim=SPINDENS_EV_DIM+ORBMOMDENS_EV_DIM;
printf("#magnetic momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm (aS(l,m) R^2(r)+ aL(l,m) F(r)) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
  if(doijk==3) snprintf(gp.title,MAXNOFCHARINLINE,"projection of momdensity M(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,MAXNOFCHARINLINE,"divergence of momdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,MAXNOFCHARINLINE,"abs value  of momdensity |ML(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'j': dim=ORBMOMDENS_EV_DIM;
printf("#currdensity is expanded in tesseral harmonics Zlm\n\
#   j(r).(%g,%g,%g)= sum_lm (b(l,m) R^2(r)+ d(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
  if(doijk==3) snprintf(gp.title,MAXNOFCHARINLINE,"projection of currdensity j(r).(i=%g,j=%g,k=%g)(milliAmp/A^2)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,MAXNOFCHARINLINE,"divergence of currdensity div j(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,MAXNOFCHARINLINE,"abs value  of currdensity |j(r)|(milliAmp/A^2)");
  dim*=6;
gp.threshhold=0.05;
break;
 default: help_and_exit();break;
}
  }
  if(strcmp(argv[os+1],"-div")==0){os+=1;doijk=1;}
  else if(strcmp(argv[os+1],"-p")==0){os+=4;
  xx=strtod(argv[3],NULL);
  yy=strtod(argv[4],NULL);
  zz=strtod(argv[5],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(xx*xx+yy*yy+zz*zz);
  xx/=rr;yy/=rr;zz/=rr;
  doijk=3;
                                     }

if(strcmp(argv[1+os],"-S")==0){os+=1;arrow=1;gp.spins_colour=3; gp.spins_scale_moment=1;//arrowdim=SPIN_EV_DIM;
                              snprintf(gp.title+strlen(gp.title),MAXNOFCHARINLINE-strlen(gp.title)," arrows correspond to the spins");}
else if(strcmp(argv[1+os],"-L")==0){os+=1;arrow=2;gp.spins_colour=2; gp.spins_scale_moment=1;//arrowdim=ORBMOM_EV_DIM;
                                   snprintf(gp.title+strlen(gp.title),MAXNOFCHARINLINE-strlen(gp.title)," arrows correspond to the orbital angular momenta");}
else if(strncmp(argv[1+os],"-M",2)==0){os+=1;arrow=3;gp.spins_colour=1; gp.spins_scale_moment=1;//arrowdim=MAGMOM_EV_DIM;
                                   snprintf(gp.title+strlen(gp.title),MAXNOFCHARINLINE-strlen(gp.title)," arrows correspond to the magnetic moments");
                                   if(strcmp(argv[os],"-Mi")==0){arrow=4;}
                                   }
else if(strncmp(argv[1+os],"-pel",4)==0){os+=1;arrow=5;gp.spins_colour=4; gp.spins_scale_moment=1;//arrowdim=MAGMOM_EV_DIM;
                                   snprintf(gp.title+strlen(gp.title),MAXNOFCHARINLINE-strlen(gp.title)," arrows correspond to the electric dipole moments");
                                   }

if(strcmp(argv[1+os],"-P")==0){os+=1;phonon=1;}
if(strcmp(argv[1+os],"-eps")==0){os+=1;eps=1;}
if(strcmp(argv[1+os],"-fst")==0){os+=1;fst=1;}
}
if(strcmp(argv[1+os],"-prefix")==0){strcpy(prefix,argv[2+os]); // read prefix
                                   fprintf(stdout,"# prefix for input filenames: %s\n",prefix);
 				   os+=2;}
 

 strcpy(infilename,"./results/");strcpy(infilename+10,prefix);
 strcpy(infilename+10+strlen(prefix),"mcphas.mf");fin = fopen(infilename, "rb");
 if(fin==NULL){strcpy(infilename+10,"mcphas.mf");fin = fopen_errchk(infilename, "rb");}
 printf("# reading from file %s\n",infilename);
 
 }

// --------------------- load crystal structure information from mcphas.j for further processing  -----------------
//   (needed by check_for_best () to convert Habc into Hijk , then for spins.out crystallographic info ...)
 FILE * fj;
 strcpy(mcphasjfilename,"./");strcpy(mcphasjfilename+2,prefix); // try prefix filename
  strcpy(mcphasjfilename+2+strlen(prefix),"mcphas.j");fj = fopen(mcphasjfilename, "rb");
 if(fj==NULL){strcpy(mcphasjfilename,"./results/");strcpy(mcphasjfilename+10,prefix); // try prefix filename
              strcpy(mcphasjfilename+10+strlen(prefix),"mcphas.j");fj = fopen(mcphasjfilename, "rb");
   if(fj==NULL){strcpy(mcphasjfilename+2,"mcphas.j");fj = fopen_errchk(mcphasjfilename, "rb");}
              }
 fclose(fj);
fprintf(stderr,"# loading crystal structure from %s ... \n",mcphasjfilename);fflush(stderr);
 par inputpars(mcphasjfilename);

char *out[NOF_USERDEF_MCPHAS_COLS+1];
for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col){out[col]=new char[20];}

 if (strncmp(argv[1],"-t",2)!=0&&strcmp(argv[1],"-fc")!=0){
  fout = fopen_errchk ("./results/spins.out", "w"); // unless it is table option
   cs.print_mcdiff_in_header(fout,"spins",0);
fprintf(fout,"\
#   0.73250   0.00000                       -1.43200 -1.43200  0.71600  0.00000 -0.71600  1.43200  0.00000\n\
#   0.73250   0.00000                       -2.56800 -2.56800  1.28400  0.00000 -1.28400  2.56800  0.00000\n\
#\n\
#\n\
# %%SECTION 3%% DESCRIPTION OF THE LATTICE\n\
#\n\
# -----------------------------------------------------------------------------\n");

// input file header and conf------------------------------------------------------------------
   headerinput(fin,fout,gp,cs,out,N);}
else
  {
// input file header and conf------------------------------------------------------------------
   headerinput(fin,stderr,gp,cs,out,N);
  }
   if(cs.nofatoms<1){fclose (fin);fprintf(stderr,"#!!! Error program spins reading nofatoms=%i - must be >0 !!!\n",cs.nofatoms);exit(1);}
   if(cs.nofcomponents<1){fclose (fin);fprintf(stderr,"#!!! Error program spins reading nofcomponents=%i - must be >0 !!!\n",cs.nofcomponents);exit(1);}

   spincf savmf(1,1,1,cs.nofatoms,cs.nofcomponents);

// ------------------------load spinsconfigurations and check which one is nearest -------------------------------   
double aa[NOF_USERDEF_MCPHAS_COLS+1];for(int i=1;i<=NOF_USERDEF_MCPHAS_COLS;++i)aa[i]=1e100;
int apf=0;
for(i=os+1;i<argc;++i){
if(strcmp(argv[i],"-prefix")==0){apf=2;}
                       }
double lnZ,U,x,y; 
if (strncmp(argv[1],"-f",2)==0&&argc<3+os){
                 // in this case we have exactly one argument left,
                 // i.e.  aa[1] becomes a number of a spinconfig in a file
                 aa[0]=strtod(argv[1+os],NULL);printf("# the configuration number %g\n",aa[1]);
                 }
                           
else{if(argc-1==2+os||argc-1==6+os+apf){
               // in this case  x and y in the phasediagram are given in the command line
               // arguments. We have to find out, which columns ix and iy
               //  correspond to x and y and put into corresponding aa[ix] and aa[iy] the
               // remaining arguments
                aa[0]=0;
                i=find_usrdef_out("x",out);if(i>0)aa[i]=strtod(argv[1+os],NULL);
                i=find_usrdef_out("y",out);if(i>0)aa[i]=strtod(argv[2+os],NULL);
               }
     else
     if(argc-1==4+os||argc-1==8+os+apf)
     {// in this case  T, Ha, Hb, Hc in the phasediagram are given in the command line
               // arguments. We have to find out, which columns iT iHa iHb  and iHc
               //  correspond to T Ha Hb and Hc and put into corresponding aa[iT],... and aa[iHc] the
               // remaining arguments
              aa[0]=0;
             i=find_usrdef_out("T",out);if(i>0)aa[i]=strtod(argv[1+os],NULL);
             i=find_usrdef_out("Ha",out);if(i>0)aa[i]=strtod(argv[2+os],NULL);
             i=find_usrdef_out("Hb",out);if(i>0)aa[i]=strtod(argv[3+os],NULL);
             i=find_usrdef_out("Hc",out);if(i>0)aa[i]=strtod(argv[4+os],NULL);
     }
     else
     if(argc-1==NOF_USERDEF_MCPHAS_COLS+os||argc-1==NOF_USERDEF_MCPHAS_COLS+4+os+apf)
     {aa[0]=0;for(int i=1;i<=NOF_USERDEF_MCPHAS_COLS;++i)aa[i]=strtod(argv[i+os],NULL);
     } else
     {
      fprintf(stderr,"#Error program spins - wrong number of arguments %i (optional %i)!\n",argc-1,os);
      exit(1);
     }
    } 
double T=0; Vector Hext(1,HEXT_DIMENSION);

if(check_for_best(fin,aa,savmf,x,y,T,Hext,outstr,out,inputpars.cs.abc))
  {fclose (fin);fprintf(stderr,"#!!! Error program spins - no stable structure found !!!\n");exit(1);}
fclose (fin);
// ----------------------------output configuration ----------------------------------------------------------------
for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col){delete []out[col];}

  printf("#! %s - configuration\n",outstr);
  fprintf(stdout,"# Components of the Demagnetisation Tensor in SI Units\n");
    fprintf(stdout,"# refering to ijk coordinate system\n");
    fprintf(stdout,"# defined by  j||b, k||(a x b) and i normal to k and j\n");

    fprintf(stdout,"#!Nii=%g Nij=%g Nik=%g\n",N(1,1),N(1,2),N(1,3));
    fprintf(stdout,"#!Njj=%g Njk=%g\n",N(2,2),N(2,3));
    fprintf(stdout,"#!Nkk=%g\n",N(3,3));
if(nofcomp>savmf.nofcomponents){nofcomp=savmf.nofcomponents;}
if(nofcomp<savmf.nofcomponents){printf("#! printing only nofcomponents=%i components\n",nofcomp);}
  if (strncmp(argv[1],"-t",2)!=0){
  if(strcmp(argv[1],"-fc")==0){
if(strcmp(argv[os]+strlen(argv[os])-3,".mf")==0)
{savmf.print_commented(stdout,"Hex",minl,maxl,maxn,limit);}
else{savmf.print_commented(stdout,"I",minl,maxl,maxn,limit);}
                       exit(0);} // exit for -fc option
  else {savmf.print(stdout,nofcomp);}}

  int ii,nt,k,j;
// determine primitive magnetic unit cell
Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3),dd3(1,3),dd(1,3);
if (inputpars.cs.r!=cs.r){cs.r=inputpars.cs.r;}
if (inputpars.cs.abc!=cs.abc){cs.abc=inputpars.cs.abc;}

// check sipffilenames -> if these are not present in cs (from headerinput), put also atomic positions from inputpars
for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
{if(cs.sipffilenames[ii]==NULL)cs.sipffilenames[ii]=new char[MAXNOFCHARINLINE];
 if (strcmp((*inputpars.jjj[ii]).sipffilename,cs.sipffilenames[ii])!=0)
             {strcpy(cs.sipffilenames[ii],(*inputpars.jjj[ii]).sipffilename);
              cs.x[ii]=inputpars.cs.x[ii];
              cs.y[ii]=inputpars.cs.y[ii];
              cs.z[ii]=inputpars.cs.z[ii];
// HERE take care about atoms sitting at nearly the same position (a nucleus and a mangetic shell
// from makenn -cfph ... and put the positions exactly at the same correct value 
// (makenn had shifted the magnetic charge cloud by 0.01 A along c in order to enable 
// the correct evaluation of cf-phonon interactions by mcphas and mcdisp. here we correct for
// this shift in order to get the right output of charge densities and in spins.out the
// right positions for doing mcdiff !!
                   for(int i=1;i<ii;++i){//loop all atoms which have been read
                                        if ( fabs((cs.x[ii]-cs.x[i])*cs.abc[1])<0.2 &&
                                             fabs((cs.y[ii]-cs.y[i])*cs.abc[2])<0.2 &&
                                             fabs((cs.z[ii]-cs.z[i])*cs.abc[3])<0.2 )
                                            { // atom i and ii are the same atom therefore check if they
                                              // are displaced along c and move magnetic atom back to nuclear position
                                              if(fabs((cs.x[ii]-cs.x[i])*cs.abc[1])>0.001 ||
                                                 fabs((cs.y[ii]-cs.y[i])*cs.abc[2])>0.001){fprintf(stderr,"Error spins.c: atoms %i (%g %g %g) and %i (%g %g %g) too close\n",i,cs.x[i],cs.y[i],cs.z[i],ii,cs.x[ii],cs.y[ii],cs.z[ii]);exit(EXIT_FAILURE);}
                                              if(cs.z[ii]>cs.z[i]){cs.z[ii]=cs.z[i];}else{cs.z[i]=cs.z[ii];}
                                            }
                                       }
             }
}
cs4.abc=cs.abc;cs4.r=cs.r;cs4.nofatoms=cs.nofatoms;cs4.nofcomponents=cs.nofcomponents;
savmf.calc_prim_mag_unitcell(p,cs.abc,cs.r);
  
  if (strncmp(argv[1],"-f",2)==0) 
 { inputpars.savelattice(fout);
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(fout,"#! strain tensor: eps1=%4.4g (epsii) eps2=%4.4g (epsjj) eps3=%4.4g (epskk) eps4=%4.4g (2epsjk) eps5=%4.4g (2epsik) eps6=%4.4g (2epsij)\n",
    myround(savmf.epsilon(1)),myround(savmf.epsilon(2)),myround(savmf.epsilon(3)),myround(savmf.epsilon(4)),myround(savmf.epsilon(5)),myround(savmf.epsilon(6)));

  fprintf(fout,"#! %s \n",outstr);
   if(T==0){fprintf(fout,"# program spins: temperature not found in %s - setting T=1 K\n",argv[2]);T=1;}
  fprintf(fout,"#!%s: nr1=%i nr2=%i nr3=%i \n",outstr,savmf.na(),savmf.nb(),savmf.nc());
  fprintf(fout,"#1           2     3     4     5       6       7       8    9    10\n");
  fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ia> <Ib> <Ic> [created by program spins]\n");
  for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
  {for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
   {// output the positions
    dd3=savmf.pos_dabc(i,j,k,ii, cs);
   //returns position dd3 as components	    with respect to lattice a b c
   
   // dd0 is position as components with respect to primitive crystallographic unit cell
    dd0=savmf.pos_dr123(i,j,k,ii, cs);
    fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            (*inputpars.jjj[ii]).sipffilename,dd3(1),dd3(2),dd3(3),dd0(1),dd0(2),dd0(3));
    //output the "magnetic" moment if possible ... actually it outputs Ia Ib Ic 
    for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,savmf.m(i,j,k)(inputpars.cs.nofcomponents*(ii-1)+nt)));}      
    fprintf(fout,"\n");
   }
  }}}
  

  fclose(fout);
    gp.spins_scale_moment=1;gp.phonon_scale_static_displacements=1;
    gp.read();
    fin = fopen_errchk ("./results/spins_prim.jvx", "w"); // here draw some graphics for the spinconfiguration
                       // with option -f. Yet this is very limited, because we do not know what is the magnetic
                       // moment, i.e. a guess is made using the components 1 2 3 of the savmf 
     gp.showprim=1;
     gp.scale_density_vectors=0;gp.show_density=0; // do not show charge densities
             Vector hkl1(1,3);hkl1=0;
             Vector gjmbHxc1(1,3);gjmbHxc1=0;
     spincf magmom(savmf.na(),savmf.nb(),savmf.nb(),savmf.nofatoms,3); // the magnetic moment guess 
     int i,j,k,l;for (i=1;i<=savmf.na();++i){for (j=1;j<=savmf.nb();++j){for (k=1;k<=savmf.nc();++k){for(l=1;l<=magmom.nofatoms;++l)  
         for(int momdim=1;momdim<=3&&momdim<=savmf.nofcomponents;++momdim)
               {magmom.moment(i,j,k,l)(momdim)=savmf.moment(i,j,k,l)(momdim);}
        }}}

     savmf.jvx_cd(fin,outstr,cs,
                  gp,0.0,savmf*0.0,savmf*0.0,
                  hkl1,T,gjmbHxc1,Hext,cs,magmom,magmom * 0.0,magmom* 0.0,magmom * 0.0);
    fclose (fin);

  exit(0);
}
// FROM HERE ON IT IS ONLY EXECUTED IF GRAPHICS ARE DESIRED ... 

gp.read();

  Vector hh(1,savmf.nofcomponents*savmf.nofatoms);
  spincf densitycf(savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,dim);
  ii=0; 
  // if individual ions of the cluster are to be shown make nofatoms in spincf larger !
  for (j=1;j<=inputpars.cs.nofatoms;++j){
         if(arrow==4&&(*inputpars.jjj[j]).module_clust==true){for(k=1;k<=(*(*inputpars.jjj[j]).clusterpars).cs.nofatoms;++k)
                                                 {++ii;par inputpars4((*(*inputpars.jjj[j]).clusterpars));
                                                 cs4.x[ii]=cs.x[j]+(*inputpars4.jjj[k]).xyz(1);
                                                 cs4.y[ii]=cs.y[j]+(*inputpars4.jjj[k]).xyz(2);
                                                 cs4.z[ii]=cs.z[j]+(*inputpars4.jjj[k]).xyz(3);cs4.sipffilenames[ii]=new char[MAXNOFCHARINLINE];
                                                 strcpy(cs4.sipffilenames[ii],(*inputpars4.jjj[k]).sipffilename);
                                                 //check if cluster abc  are the same as inputpars
           if(fabs(cs4.abc(1)-inputpars4.cs.abc(1))>1e-5)
              {fprintf(stderr,"Error program spins - a=%g in cluster %s not the same as a=%g in mcphas.j\n",cs4.abc(1),(*inputpars.jjj[j]).sipffilename,inputpars4.cs.abc(1));exit(1);}
           if(fabs(cs4.abc(2)-inputpars4.cs.abc(2))>1e-5)
              {fprintf(stderr,"Error program spins - b=%g in cluster %s not the same as b=%g in mcphas.j\n",cs4.abc(2),(*inputpars.jjj[j]).sipffilename,inputpars4.cs.abc(2));exit(1);}
           if(fabs(cs4.abc(3)-inputpars4.cs.abc(3))>1e-5)
              {fprintf(stderr,"Error program spins - c=%g in cluster %s not the same as c=%g in mcphas.j\n",cs4.abc(3),(*inputpars.jjj[j]).sipffilename,inputpars4.cs.abc(3));exit(1);}
           if(fabs(cs4.abc(4)-inputpars4.cs.abc(4))>1e-5)
              {fprintf(stderr,"Error program spins - alpha=%g in cluster %s not the same as alpha=%g in mcphas.j\n",cs4.abc(4),(*inputpars.jjj[j]).sipffilename,inputpars4.cs.abc(4));exit(1);}
           if(fabs(cs4.abc(5)-inputpars4.cs.abc(5))>1e-5)
              {fprintf(stderr,"Error program spins - beta=%g in cluster %s not the same as beta=%g in mcphas.j\n",cs4.abc(5),(*inputpars.jjj[j]).sipffilename,inputpars4.cs.abc(5));exit(1);}
           if(fabs(cs4.abc(6)-inputpars4.cs.abc(6))>1e-5)
              {fprintf(stderr,"Error program spins - gamma=%g in cluster %s not the same as gamma=%g in mcphas.j\n",cs4.abc(6),(*inputpars.jjj[j]).sipffilename,inputpars4.cs.abc(6));exit(1);}                                                 
                                                 }
                                               }
                                               else
                                              {++ii;cs4.x[ii]=cs.x[j];cs4.y[ii]=cs.y[j];cs4.z[ii]=cs.z[j];cs4.sipffilenames[ii]=cs.sipffilenames[j];}
              } cs4.nofatoms=ii;

  spincf spinconf(savmf.na(),savmf.nb(),savmf.nc(),ii,3);
  spincf sc_phonon(savmf.na(),savmf.nb(),savmf.nc(),ii,3);

fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(fout,"#! strain tensor: eps1=%4.4g (epsii) eps2=%4.4g (epsjj) eps3=%4.4g (epskk) eps4=%4.4g (2epsjk) eps5=%4.4g (2epsik) eps6=%4.4g (2epsij)\n",
       myround(savmf.epsilon(1)),myround(savmf.epsilon(2)),myround(savmf.epsilon(3)),myround(savmf.epsilon(4)),myround(savmf.epsilon(5)),myround(savmf.epsilon(6)));

if (strncmp(argv[1],"-t",2)!=0){
// the following is for the printout of spins.out ...........................
section4header(fout);
fprintf(fout,"#! %s : ",outstr);
fprintf(fout," nr1=%i nr2=%i nr3=%i nat=%i  atoms in primitive magnetic unit cell\n",savmf.na(),savmf.nb(),savmf.nc(),cs4.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
//MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
//                                        to Sa Sb Sc La Lb Lc
fprintf(fout,"#1           2     3     4     5       6       7       8    9    10                  11   12   13   14   15   16\n");
fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)\n");
fprintf(fout,"#          corresponding exchange fields [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)\n");
// .............................................................................                                
	}
else  //now table options
 { fprintf(fout,"#! nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",savmf.na(),savmf.nb(),savmf.nc(),cs4.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
 //fprintf(fout,"# 1 2  3  4  5           6     7     8     9       10      11      12   13   14                  15   16   17   18   19   20\n");
 fprintf(fout,"# external parameters {sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] ");
 if (strcmp(argv[1],"-tMSL")==0)  fprintf(fout,"<Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)]\n");
 else if (strcmp(argv[1],"-tI")==0) fprintf(fout,"<I1> <I2> <I3> ... <Inofcomponents>\n");
 else if (strcmp(argv[1],"-tHex")==0) fprintf(fout,"<Hex1> <Hex2> <Hex3> ... <Inofcomponents> (meV)\n");
 else if (strcmp(argv[1],"-tL2")==0){ fprintf(fout,"bondlength_to_neighbour [A]  elongation^2 [A^2]\n");
                                nlimits_calc(nmin,nmax, tL2,p);
                               }
 else {fprintf(stderr, "Error spins: option %s not known\n",argv[1]);exit(EXIT_FAILURE);}
 }       
//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: Icalc has to be used to calculate all the <Olm>.
 
Vector h(1,inputpars.cs.nofcomponents);
Vector I(1,inputpars.cs.nofcomponents);
h=0;for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
{(*inputpars.jjj[ii]).Icalc_parameter_storage_init(h,Hext,T);} // initialize Icalc module parameter storage

 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
 {
    hh=savmf.m(i,j,k);
  densitycf.m(i,j,k)=0;int i4=1,ii4=1;
  for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
 {  
    Vector magmom(1,3),mom(1,3);
    Vector Lmom(1,3);
    Vector Smom(1,3);
    Vector moments(1,dim);
      Vector momS(1,SPINDENS_EV_DIM);
  Vector momL(1,ORBMOMDENS_EV_DIM);
  Vector momentsx(1,SPINDENS_EV_DIM);
  Vector momentsy(1,SPINDENS_EV_DIM);
  Vector momentsz(1,SPINDENS_EV_DIM);
  Vector momentlx(1,ORBMOMDENS_EV_DIM);
  Vector momently(1,ORBMOMDENS_EV_DIM);
  Vector momentlz(1,ORBMOMDENS_EV_DIM);
    h=0;
   for(nt=1;nt<=inputpars.cs.nofcomponents;++nt){h(nt)=hh(nt+inputpars.cs.nofcomponents*(ii-1));}


switch(arrow)
{case 1: (*inputpars.jjj[ii]).Scalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 2: (*inputpars.jjj[ii]).Lcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 3: (*inputpars.jjj[ii]).mcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 5: (*inputpars.jjj[ii]).pelcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 default: break;
}

switch(arrow)
{case 1:
 case 2:
 case 3:
 case 5:
         for(nt=1;nt<=3;++nt)
		        {spinconf.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };break;
 case 4: int dim4;
         dim4=3;if((*inputpars.jjj[ii]).module_clust==true)
                           dim4=(*(*inputpars.jjj[ii]).clusterpars).cs.nofatoms*3;
         Vector momi(1,dim4);
         (*inputpars.jjj[ii]).micalc(momi,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);
         // now put the components of momi to the spinconf                                               
         for(nt=1;nt<=dim4;++nt){spinconf.m(i,j,k)(nt+i4-1)=momi(nt);};i4+=dim4;
         break;
}

// output atoms and moments in primitive unit cell to fout  ------------------------------------
if(arrow==4&&(*inputpars.jjj[ii]).module_clust==true){
    for(nt=1;nt<=(*(*inputpars.jjj[ii]).clusterpars).cs.nofatoms;++nt)
     {dd3=spinconf.pos_dabc(i,j,k,ii4, cs4);
      dd0=spinconf.pos_dr123(i,j,k,ii4, cs4);
      if (strncmp(argv[1],"-t",2)==0){
              fprintf(fout," %s ", outstr);
                                     }
      fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            cs4.sipffilenames[ii4],dd3(1),dd3(2),dd3(3),dd0(1),dd0(2),dd0(3));
      if (strncmp(argv[1],"-t",2)!=0||strcmp(argv[1],"-tMSL")==0){
            fprintf(fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(1+3*(ii4-1))));
            fprintf(fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(2+3*(ii4-1))));
            fprintf(fout," %4.4f\n",myround(1e-5,spinconf.m(i,j,k)(3+3*(ii4-1))));
                                                                     } 
       // do not print out L S or exchange fields for cluster module
      ++ii4;
     }
    }
else{   ++ii4;
    // output the positions
    dd3=savmf.pos(i,j,k,ii, cs);
      // if module allows to calculate position shift of an atom - use this for output 
    if(true==(*inputpars.jjj[ii]).pcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
     {dd3+=mom; 
     //  fprintf(stderr,"# Attention: atom %i shifted from equilibrium position by (%g %g %g) A \n",ii,mom(1),mom(2),mom(3));
     }
     // if module cannot calculate position shift of an atom - check if the atom 
     // is actually a charge cloud on top of another atom and if yes take position shift
     // from this other atom
      else
     {for(int iii=1;iii<=inputpars.cs.nofatoms;++iii)
      {if(ii!=iii)
       {Vector ddc(1,3); ddc=savmf.pos(i,j,k,iii, cs);
        if(abs(ddc-dd3)<0.01){
Vector hhh(1,inputpars.cs.nofcomponents); hhh=0;
for(nt=1;nt<=inputpars.cs.nofcomponents;++nt){hhh(nt)=hh(nt+inputpars.cs.nofcomponents*(iii-1));}

 if(true==(*inputpars.jjj[iii]).pcalc(mom,T,hhh,Hext,(*inputpars.jjj[iii]).Icalc_parstorage))
    {dd3+=mom;
 //printf(stderr,"# Attention: atom %i shifted from equilibrium position by (%g %g %g) A taken from atom %i\n",ii,mom(1),mom(2),mom(3),iii);
 }
                               }
       }
      } 
     }

    dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
    Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
    dd=abc_in_ijk.Inverse()*dd3; 
    snprintf(tstr,MAXNOFCHARINLINE,"{%s} %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f ",
            cs.sipffilenames[ii],dd(1),dd(2),dd(3),dd0(1),dd0(2),dd0(3));
    if (strcmp(argv[1],"-tL2")!=0) {     
    if (strncmp(argv[1],"-t",2)==0){
       fprintf(fout,"%s ",outstr);
                                   }
      fprintf(fout,"%s",tstr);
    }
    if (strncmp(argv[1],"-t",2)!=0||strcmp(argv[1],"-tMSL")==0)
    {
     //ouput the magnetic moment if possible
     if((*inputpars.jjj[ii]).mcalc(magmom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
     {     for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,magmom(nt)));}
      // and output the orbital and spin momentum if possible 
      if((*inputpars.jjj[ii]).Lcalc(Lmom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage)&&
        (*inputpars.jjj[ii]).Scalc(Smom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
       //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
       //                                        to Sa Sb Sc La Lb Lc
       //{        for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f %4.4f",myround(1e-5,Smom(nt)),myround(1e-5,Lmom(nt)));}
       {        for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,Smom(nt)));}
               for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,Lmom(nt)));}
       }
      }}
    }
   if (strncmp(argv[1],"-t",2)!=0||strcmp(argv[1],"-tHex")==0)
   {
    // finally output a line with the exchange fields 
    if (strncmp(argv[1],"-t",2)!=0)fprintf(fout,"\n                 corresponding exchange fields [meV]-->          ");
                      for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout exchangefields
                        {fprintf(fout," %4.4f",myround(1e-5,h(nt)));}
                        
   }
   if (strcmp(argv[1],"-tI")==0)
   {(*inputpars.jjj[ii]).Icalc(I,T,h,Hext,lnZ,U,(*inputpars.jjj[ii]).Icalc_parstorage);
                         for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout I operator expectation values
                        {fprintf(fout," %4.4f",myround(1e-5,I(nt)));}
                         
   }
   if (strcmp(argv[1],"-tL2")==0)
   {// here calculate up to tl2 Angstroem bondlengths and elongation^2 and output to fout
    L2_calc(tL2,i,j,k,ii,fout,outstr,tstr,nmin,nmax,savmf,p,T,h,Hext,inputpars);
   }
   else   fprintf(fout,"\n");

// -----------------------------------------------------------------------------------------------


if(density){
switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c':  (*inputpars.jjj[ii]).chargedensity_coeff (moments, T, h, Hext, (*inputpars.jjj[ii]).Icalc_parstorage); break;
 case 's':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsy,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=SPINDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'o':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsy,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+ORBMOMDENS_EV_DIM)=momentsy(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'm':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsy,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            momS=xx*momentsx+yy*momentsy+zz*momentsz;
            if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentlx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momently,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentlz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            momL=xx*momentlx+yy*momently+zz*momentlz;
            for(int i=1;i<=SPINDENS_EV_DIM;++i){
            if(doijk==3){moments(i)=momS(i);moments(i+SPINDENS_EV_DIM)=momL(i);}
            else{moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);
                 moments(i+3*SPINDENS_EV_DIM)=momentlx(i);moments(i+4*SPINDENS_EV_DIM)=momently(i);moments(i+5*SPINDENS_EV_DIM)=momentlz(i);
                }
                                               }
            break;
 case 'j':  (*inputpars.jjj[ii]).orbmomdensity_coeff (momentlx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            (*inputpars.jjj[ii]).orbmomdensity_coeff (momently,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            (*inputpars.jjj[ii]).orbmomdensity_coeff (momentlz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){
             moments(i)=momentlx(i);moments(i+ORBMOMDENS_EV_DIM)=momently(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentlz(i);
             }
            break;
 default: help_and_exit();
}
                   for(nt=1;nt<=dim;++nt)
		        {densitycf.m(i,j,k)(nt+dim*(ii-1))=moments(nt);
                    }
} // gp.show_density


if(phonon==1)  // if module allows to calculate position  - use this for graphics ...
{if(true==(*inputpars.jjj[ii]).pcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
{for(nt=1;nt<=3;++nt)
		        {sc_phonon.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };
}
}

  }
}
}}
             
if (strncmp(argv[1],"-t",2)==0){exit(0);}
  fclose (fout);
   
// create plot of spinconfiguration -----------------------------------------------------------
printf("# ************************************************************************\n");
printf("#%s\n",gp.title);
printf("# ************************************************************************\n");
              if(arrow==0)gp.spins_scale_moment=0;
              if(density==0)gp.show_density=0;
if(eps==1){
    fin = fopen_errchk ("./results/spins.eps", "w");
     spinconf.eps(fin,outstr);
    fclose (fin);

// here the 3d file should be created
    fin = fopen_errchk ("./results/spinsab.eps", "w");


     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,1,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spinsac.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,2,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spinsbc.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,3,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spins3dab.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,4,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spins3dac.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,5,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spins3dbc.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,6,spinconf);
    fclose (fin);
          }
if(fst==1){
    fin = fopen_errchk ("./results/spins.fst", "w");
     spinconf.fst(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,spinconf);
    fclose (fin);

    
   fin = fopen_errchk ("./results/spins_prim.fst", "w");
     spinconf.fstprim(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,spinconf);
    fclose (fin);
          }
             Vector hkl(1,3);hkl=0;
             Vector gjmbHxc(1,3);gjmbHxc=0;
             spincf densityev_real(densitycf*0.0);
             spincf densityev_imag(densitycf*0.0);
             spincf spinconfev_real(spinconf*0.0);
             spincf spinconfev_imag(spinconf*0.0);
             spincf spinconfpev_real(spinconf*0.0);
             spincf spinconfpev_imag(spinconf*0.0);
            // to do jvx output of static structure put zeros into these spinconfigurations

// check sipffilenames and put radius= ... in case single ion module is
//  capable of calculating position 
for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
{Vector pos(1,3);   
  for(nt=1;nt<=inputpars.cs.nofcomponents;++nt){h(nt)=hh(nt+inputpars.cs.nofcomponents*(ii-1));}
  if(true==(*inputpars.jjj[ii]).pcalc(pos,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
 {double charge;charge=(*inputpars.jjj[ii]).charge;if(charge==0)charge=0.01;
   snprintf(cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g",
           charge,gp.scale_pointcharges*0.529177*signum(charge)*pow((double)fabs(charge),0.3333));
  // if we can determine from sipffilename the element --> put r g b information
    
  unsigned int r,g,b; char element []  ="E\0";
   element[0]=toupper((*inputpars.jjj[ii]).sipffilename[0]);
   if(isalpha((*inputpars.jjj[ii]).sipffilename[1]))
   element[1]=tolower((*inputpars.jjj[ii]).sipffilename[1]);
  for(int i=0;i<NOFELEMENTS_RGB;++i)
   if(0==strncmp(element,elstr[i],2)){
r=(unsigned int)(unsigned char)elstr[i][2];
g=(unsigned int)(unsigned char)elstr[i][3];
b=(unsigned int)(unsigned char)elstr[i][4];
   snprintf(cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g %s r=%i g=%i b=%i",
           charge,gp.scale_pointcharges*0.529177*signum(charge)*pow((double)fabs(charge),0.3333),
           element,r,g,b);}
   // printf("%s\n",cs.sipffilenames[ii]);

// printf("#! atom %i %s displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,cs.sipffilenames[ii],ii,pos(1),ii,pos(2),ii,pos(3));
  printf("#! atom %i  displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,ii,pos(1),ii,pos(2),ii,pos(3));
 }
}

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin = fopen_errchk ("./results/spins.jvx", "w");
    gp.showprim=0;gp.spins_wave_amplitude=0;
     densitycf.jvx_cd(fin,outstr,cs,gp,0.0,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon);
    fclose (fin);

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin = fopen_errchk ("./results/spins_prim.jvx", "w");
     gp.showprim=1;
     densitycf.jvx_cd(fin,outstr,cs,gp,0.0,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon);
    fclose (fin);

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
// try movie - a spinwave picture  ... including phonons and spindensity changes ...
//***************************************************************************************************************
//***************************************************************************************************************
// check if there is a second 
for(i=os+1;i<argc;++i){
if(strcmp(argv[i],"-prefix")==0){strcpy(prefix,argv[i+1]); // read alternative prefix
                                   fprintf(stdout,"# prefix for qee qep etc input filenames: %s\n",prefix);
 				   os+=2;}
                       }
if(argc-1==NOF_USERDEF_MCPHAS_COLS+4+os)os+=NOF_USERDEF_MCPHAS_COLS-4;
if(argc-1==6+os)os-=2;  // 
//argc-1==8+os  ... ok 
if (argc-1==8+os){
              // double E;
              char outhklstr[MAXNOFCHARINLINE];
             gp.spins_wave_amplitude=1.0;gp.spins_show_ellipses=1.0;gp.spins_show_oscillation=1.0;
             gp.phonon_wave_amplitude=1.0;gp.phonon_scale_static_displacements=1.0;
//----------------------------------------------------------------------------------------------------------
           Vector thkl(1,3);
                
           if(arrow>0){
             strcpy(infilename,"./results/");strcpy(infilename+10,prefix);
             switch(arrow)
             {case 1: strcpy(infilename+10+strlen(prefix),"mcdisp.qes");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qes", "rb");
                      break;
              case 2: strcpy(infilename+10+strlen(prefix),"mcdisp.qeo");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qeo", "rb");
                      break;
              case 3: strcpy(infilename+10+strlen(prefix),"mcdisp.qem");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qem", "rb");
                      break;
              case 4: fprintf(stderr,"mcdisp: output of individual moment oscillation in eigenvector file mcdisp.qemi not yet implemented - thus exiting program spins\n");
                      strcpy(infilename+10+strlen(prefix),"mcdisp.qemi");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qemi", "rb");
                      break;   
              case 5: strcpy(infilename+10+strlen(prefix),"mcdisp.qpe");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qpe", "rb");
                      break;                   
             }
             check_for_best_excitation_and_close(fin,gp,"qev",spinconf,spinconfev_real,spinconfev_imag,
                          x,y,T,Hext,outstr,outhklstr,"moment oscillation",thkl,hkl,argv,os,
                          inputpars.cs.abc,inputpars.cs.nofatoms,doijk, xx,yy,zz,dim);

             }//arrow
     
//----------------------------------------------------------------------------------------------------------
            if(phonon>0){
             strcpy(infilename,"./results/");strcpy(infilename+10,prefix);
             strcpy(infilename+10+strlen(prefix),"mcdisp.qep");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qep", "rb");
             check_for_best_excitation_and_close(fin,gp,"qep",spinconf,spinconfpev_real,spinconfpev_imag,
                          x,y,T,Hext,outstr,outhklstr,"phonon oscillation",thkl,hkl,argv,os,
                          inputpars.cs.abc,inputpars.cs.nofatoms,doijk, xx,yy,zz,dim);

             }//phonon
//----------------------------------------------------------------------------------------------------------
            if(density){
             strcpy(infilename,"./results/");strcpy(infilename+10,prefix);
             switch(argv[1][1]) // dimension definition from jjjpar.hpp
                {case 'c': strcpy(infilename+10+strlen(prefix),"mcdisp.qee");fin = fopen(infilename, "rb");
                           if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qee", "rb");
                           break;
                 case 's': strcpy(infilename+10+strlen(prefix),"mcdisp.qsd");fin = fopen(infilename, "rb");
                           if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qsd", "rb");
                           break;
                 case 'o': strcpy(infilename+10+strlen(prefix),"mcdisp.qod");fin = fopen(infilename, "rb");
                           if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qod", "rb");
                           break;
                 case 'm': fprintf(stderr,"Error spins: magnetic moment density oscillation not yet implemented\n");exit(1);break;
                             // would have to look into qsd and qod files !!
                 case 'j': strcpy(infilename+10+strlen(prefix),"mcdisp.qod");fin = fopen(infilename, "rb");
                           if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qod", "rb");
                           break;
                }
             check_for_best_excitation_and_close(fin,gp,"qee/qsd/qod",densitycf,densityev_real,densityev_imag,
          x,y,T,Hext,outstr,outhklstr,"density oscillation",thkl,hkl,argv,os,
          inputpars.cs.abc,inputpars.cs.nofatoms,doijk, xx,yy,zz,dim);

             }//gp.show_density
//----------------------------------------------------------------------------------------------------------           
              gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
                        // in case he wants to overwrite some default settings
              if(arrow==0)gp.spins_scale_moment=0;
              if(density==0)gp.show_density=0;
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
             // complex <double> im(0,1);
              for(i=0;i<16;++i)
              {phase=2*3.1415*i/15;
               printf("\n********************************************\n");
               printf(" calculating movie sequence %i(16)\n",i+1);
               printf("********************************************\n");
               char filename[MAXNOFCHARINLINE];
               snprintf(filename,MAXNOFCHARINLINE,"./results/spins.%i.jvx",i+1);
               fin = fopen_errchk (filename, "w");gp.showprim=0;
                     densitycf.jvx_cd(fin,outhklstr,cs,gp,
                                  phase,densityev_real,densityev_imag,hkl,T,hh,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon,spinconfpev_real,spinconfpev_imag);
               fclose (fin);
               snprintf(filename,MAXNOFCHARINLINE,"./results/spins_prim.%i.jvx",i+1);
               fin = fopen_errchk (filename, "w");gp.showprim=1;
                     densitycf.jvx_cd(fin,outhklstr,cs,gp,
                                  phase,densityev_real,densityev_imag,hkl,T,hh,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon,spinconfpev_real,spinconfpev_imag);
               fclose (fin);
              }
          printf("# %s\n",outhklstr);
          }
fprintf(stderr,"# %s\n",outstr);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spins running at \n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/spins.jvx\n");
fprintf(stderr,"# * java javaview \"model=results/spins.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n");
fprintf(stderr,"# * saved density mesh in results/spins.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  for(i=1;i<=cs.nofatoms;++i){  delete cs.sipffilenames[i];}
  return 0;
}


