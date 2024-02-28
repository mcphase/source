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

void help_and_exit()
    { printf ("\n\
program spins - popout spin/exchange field configuration\n"
"              - and/or display 3d animation of spin/moment/densities and animations\n\n\
use as: spins -f[c 13 0.1] mcphas.sps T Ha Hb Hc\n\
    or: spins -f[c 13 0.2] mcphas.sps x y\n\
    or: spins -f[c 14 0.1] mcphas.tst n\n\
    or: spins -tMSL [-prefix 001] T Ha Hb Hc \n\
    or: spins -tHex [-prefix 001]  T Ha Hb Hc \n\
    or: spins -tI  [-prefix 001] T Ha Hb Hc \n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M] [-P] [-prefix 001] T Ha Hb Hc [h k l E]\n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M] [-P] [-prefix 001] x y\n\
                    \n\
1) if used with -f file T Ha Hb Hc, this file has to be a mcphas.mf or mcphas.sps file,\n \
   the spin configuration at given temperature T[K] and magnetic effective field H[T]\n \
   is read and extracted from this file and printed on screen (stdout),\n \
   results/spins.out is created (with mag moment chosen to be = <Ia> <Ib> <Ic>),\n\
   a simple graphics to represent the configuration is created in results/spins_prim.jvx \n\
  \n\
2) if used with -f file x y, then this file has to be a mcphas.mf or mcphas.sps file,\n\
   the spin configuration at a given x,y point is read and extracted from this file ,\n\
   and printed on screen (stdout) etc. as 1)\n\
3) if used with -f filen n,  this file has to be a mcphas.tst file,\n \
   the spin configuration number n\n \
   is read and extracted from this file and printed on screen (stdout),\n \
   results/spins.out is created (with mag moment chosen to be = <Ia> <Ib> <Ic>)\n\
   a simple graphics to represent the configuration is created in results/spins_prim.jvx \n\
  \n\
1&2&3) if used with -fc min max n lim a human readable format is output for spin components with index \n\
   from min to max, only n numbers exceeding  \n\
   absolute value of lim, e.g. -fc 1 3 13 0.1 outputs at maximum 13 components which are all larger \n\
   (absolute value) than 0.1 \n\
 \n\
4) if used without a filename, the information is read from results/mcphas.* results/mcdisp.*\n\
   output files and tables or 3d graphical animations are created.\n\
   for table the options are:\n\
       -tMSL ... output to stdout a table with T Ha Hb Hc atom positions and with \n\
              magnetic moments <Mx> <My> <Mz>, orbital moments and spin \n\
              of each atom in the magnetic unitc cell \n\
       -tHex ...  output to stdout table with T Ha Hb Hc atom positions and exchange fields Hex\n\
       -tI   ... a similar table with expectation values of interaction operators <I>\n\
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
 This program outputs a magnetic structure (and magnetic excitation)\n \
 graphic/movie in the output files of different format:\n \
 results/spins*.eps (postscript), results/spins*.fst (fp_studio), \n \
 results/spins.out (ascii) and results/spins*.jvx (javaview)\n\n \
 the graphics output format can be fine tuned in .sps and .qev input files\n"
" or results/graphics_parameters.set by show_abc_unitcell,\n"
" show_primitive_crystal_unitcell, spins_scale_moment, spins_wave_amplitude\n \
 show_magnetic_unitcell, show_atoms, scale_view_1,scale_view_2, scale_view_3 ...\n\n \
 jvx files can be viewed by:\n\
 java javaview results/spins.jvx \n \
 java javaview \"model=results/spins.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n"
" gif images stored by javaview can be connected to movies by ImageMagick \n"
" convert  -delay 1 -size 100x100 -loop 1 geomAnim.*.gif output.gif\n");
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
# 'Hxc1' 'Hxc2' 'Hxc3' (optional line, used to go beyond dipole approx for formfactor)\n\
#                                     denote the corresponding exchange fields in meV\n\
#\n");
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
double T=0; Vector Hext(1,3),Hextijk(1,3);
 int i,n=0,minl=1,maxl=1;//,dophon=0;
 cryststruct cs,cs4;
 //float numbers[13];numbers[9]=1;numbers[10]=3;
 //numbers[0]=13;
 char outstr[MAXNOFCHARINLINE];
 char infilename[MAXNOFCHARINLINE];
 char prefix[MAXNOFCHARINLINE];prefix[0]='\0';
 
  int dim=28;
 char text[1000];
 int os=0,maxn=0; int doijk=0,arrow=0,density=0,phonon=0;//,arrowdim=3;
 double xx=0,yy=0,zz=0,limit=0;
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
snprintf(gp.title,sizeof(gp.title),"output of program spins");

 // check command line
 if (argc < 2){help_and_exit();}
// first: options without graphics just screendump <I> or exchange field configuration at given HT
 if (strncmp(argv[1],"-f",2)==0)
 { os=2;if (strcmp(argv[1],"-fc")==0){os=6;minl=(int)strtod(argv[2],NULL);maxl=(int)strtod(argv[3],NULL);
 maxn=(int)strtod(argv[4],NULL);limit=strtod(argv[5],NULL);
   }
   fin = fopen_errchk (argv[os], "rb");printf("#* program spins ... reading from file %s\n",argv[os]);   
 }
 else { if (strncmp(argv[1],"-t",2)==0){os=1;fout=stdout;}
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
 snprintf(text,sizeof(text),"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , chargedensity ro(r)</title>\n", T,Hext(1),Hext(2),Hext(3));
 snprintf(gp.title,sizeof(gp.title),"chargedensity ro(r)");
 gp.threshhold=-0.05;
           break;
 case 's': dim=SPINDENS_EV_DIM;
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 snprintf(text,sizeof(text),"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , spindensity S(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
  if(doijk==3) snprintf(gp.title,sizeof(gp.title),"projection of spindensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,sizeof(gp.title),"divergence of spindensity div Ms(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,sizeof(gp.title),"abs value  of spindensity |Ms(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'o': dim=ORBMOMDENS_EV_DIM;
printf("#orbital momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm  aL(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 snprintf(text,sizeof(text),"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , orbital momdensity L(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
  if(doijk==3) snprintf(gp.title,sizeof(gp.title),"projection of orbmomdensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,sizeof(gp.title),"divergence of orbmomdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,sizeof(gp.title),"abs value  of orbmomdensity |ML(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'm': dim=SPINDENS_EV_DIM+ORBMOMDENS_EV_DIM;
printf("#magnetic momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm (aS(l,m) R^2(r)+ aL(l,m) F(r)) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 snprintf(text,sizeof(text),"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , magnetic momdensity M(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
  if(doijk==3) snprintf(gp.title,sizeof(gp.title),"projection of momdensity M(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,sizeof(gp.title),"divergence of momdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,sizeof(gp.title),"abs value  of momdensity |ML(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'j': dim=ORBMOMDENS_EV_DIM;
printf("#currdensity is expanded in tesseral harmonics Zlm\n\
#   j(r).(%g,%g,%g)= sum_lm (b(l,m) R^2(r)+ d(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 snprintf(text,sizeof(text),"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , currentdensity j(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
  if(doijk==3) snprintf(gp.title,sizeof(gp.title),"projection of currdensity j(r).(i=%g,j=%g,k=%g)(milliAmp/A^2)",xx,yy,zz);
  if(doijk==1){snprintf(gp.title,sizeof(gp.title),"divergence of currdensity div j(r)");gp.scale_density_vectors=0;}
  if(doijk==0) snprintf(gp.title,sizeof(gp.title),"abs value  of currdensity |j(r)|(milliAmp/A^2)");
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
                              snprintf(gp.title,sizeof(gp.title),"%s arrows correspond to the spins",gp.title);}
else if(strcmp(argv[1+os],"-L")==0){os+=1;arrow=2;gp.spins_colour=2; gp.spins_scale_moment=1;//arrowdim=ORBMOM_EV_DIM;
                                   snprintf(gp.title,sizeof(gp.title),"%s arrows correspond to the orbital angular momenta",gp.title);}
else if(strncmp(argv[1+os],"-M",2)==0){os+=1;arrow=3;gp.spins_colour=1; gp.spins_scale_moment=1;//arrowdim=MAGMOM_EV_DIM;
                                   snprintf(gp.title,sizeof(gp.title),"%s arrows correspond to the magnetic moments",gp.title);
                                   if(strcmp(argv[os],"-Mi")==0){arrow=4;}
                                   }

if(strcmp(argv[1+os],"-P")==0){os+=1;phonon=1;}
if(strcmp(argv[1+os],"-prefix")==0){strcpy(prefix,argv[2+os]); // read prefix
                                   fprintf(stdout,"# prefix for input filenames: %s\n",prefix);
 				   os+=2;}
}
 strcpy(infilename,"./results/");strcpy(infilename+10,prefix);
 strcpy(infilename+10+strlen(prefix),"mcphas.mf");fin = fopen(infilename, "rb");
 if(fin==NULL){strcpy(infilename+10,"mcphas.mf");fin = fopen_errchk(infilename, "rb");}
 printf("# reading from file %s\n",infilename);
  
 }
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
   n=headerinput(fin,fout,gp,cs);}
else
  {
// input file header and conf------------------------------------------------------------------
   n=headerinput(fin,stderr,gp,cs);
  }
   if(cs.nofatoms<1){fclose (fin);fprintf(stderr,"#!!! Error program spins reading nofatoms=%i - must be >0 !!!\n",cs.nofatoms);exit(1);}
   if(cs.nofcomponents<1){fclose (fin);fprintf(stderr,"#!!! Error program spins reading nofcomponents=%i - must be >0 !!!\n",cs.nofcomponents);exit(1);}

   spincf savmf(1,1,1,cs.nofatoms,cs.nofcomponents);

// load spinsconfigurations and check which one is nearest -------------------------------   
double TT=0; TT=strtod(argv[1+os],NULL);
double HHx=0,HHy=0,HHz=0,lnZ,U;
if (strncmp(argv[1],"-f",2)==0&&argc-os<3){TT=-TT;printf("# the configuration number %g\n",-TT);} // here TT becomes a number of a spinconfig in a file
else{if(argc<4+os){TT=0;HHx=strtod(argv[1+os],NULL);HHy=strtod(argv[2+os],NULL);
               }// here Hx and Hy become x and y in the phasediagram and TT=0 indicates this fact
     else
     {HHx=strtod(argv[2+os],NULL);HHy=strtod(argv[3+os],NULL);HHz=strtod(argv[4+os],NULL);}
     }
if(check_for_best(fin,TT,HHx,HHy,HHz,savmf,T,Hext,outstr))
  {fclose (fin);fprintf(stderr,"#!!! Error program spins - no stable structure found !!!\n");exit(1);}
fclose (fin);

  printf("#! %s - configuration\n",outstr);
  if (strncmp(argv[1],"-t",2)!=0){
  if(strcmp(argv[1],"-fc")==0){
if(strcmp(argv[os]+strlen(argv[os])-3,".mf")==0)
{savmf.print_commented(stdout,"Hex",minl,maxl,maxn,limit);}
else{savmf.print_commented(stdout,"I",minl,maxl,maxn,limit);}

exit(0);}
  else {savmf.print(stdout);}}
  par inputpars("./mcphas.j");
  int ii,nt,k,j;
// determine primitive magnetic unit cell
Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3),dd3(1,3),dd(1,3);
if (inputpars.cs.r!=cs.r){cs.r=inputpars.cs.r;}
if (inputpars.cs.abc!=cs.abc){cs.abc=inputpars.cs.abc;}

// check sipffilenames
for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
{if(cs.sipffilenames[ii]==NULL)cs.sipffilenames[ii]=new char[MAXNOFCHARINLINE];
 if (strcmp((*inputpars.jjj[ii]).sipffilename,cs.sipffilenames[ii])!=0){strcpy(cs.sipffilenames[ii],(*inputpars.jjj[ii]).sipffilename);}
}
cs4.abc=cs.abc;cs4.r=cs.r;cs4.nofatoms=cs.nofatoms;cs4.nofcomponents=cs.nofcomponents;
savmf.calc_prim_mag_unitcell(p,cs.abc,cs.r);
  
  if (strncmp(argv[1],"-f",2)==0) 
 { inputpars.savelattice(fout);
   fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(fout,"#! strain tensor: eps1=%4.4g=epsii eps2=%4.4g=epsjj eps3=%4.4g=epskk eps4=%4.4g=2epsjk eps5=%4.4g=2epsik eps6=%4.4g=2epsij\n",
    myround(savmf.epsilon(1)),myround(savmf.epsilon(2)),myround(savmf.epsilon(3)),myround(savmf.epsilon(4)),myround(savmf.epsilon(5)),myround(savmf.epsilon(6)));

  fprintf(fout,"#! %s \n",outstr);
   if(T==0){fprintf(fout,"# program spins: temperature not found in %s - setting T=1 K\n",argv[2]);T=1;}
  fprintf(fout,"#!T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",T,Hext(1),Hext(2),Hext(3),savmf.na(),savmf.nb(),savmf.nc(),cs4.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
  fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ia> <Ib> <Ic> [created by program spins]\n");
  for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
  {for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
   {// output the positions
    dd3=savmf.pos_dabc(i,j,k,ii, cs);
   //returns position dd3 as components with respect to lattice a b c
   
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
    gp.spins_scale_moment=1;
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
 //fprintf(stderr,"savmf.nofcomponents=%i\n",savmf.nofcomponents);

     savmf.jvx_cd(fin,outstr,cs,
                  gp,0.0,savmf*0.0,savmf*0.0,
                  hkl1,T,gjmbHxc1,Hextijk,cs,magmom,magmom * 0.0,magmom* 0.0);
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
         if(arrow==4&&(*inputpars.jjj[j]).module_type==5){for(k=1;k<=(*(*inputpars.jjj[j]).clusterpars).cs.nofatoms;++k)
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

fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(fout,"#! strain tensor: eps1=%4.4g=epsii eps2=%4.4g=epsjj eps3=%4.4g=epskk eps4=%4.4g=2epsjk eps5=%4.4g=2epsik eps6=%4.4g=2epsij\n",
    myround(savmf.epsilon(1)),myround(savmf.epsilon(2)),myround(savmf.epsilon(3)),myround(savmf.epsilon(4)),myround(savmf.epsilon(5)),myround(savmf.epsilon(6)));

if (strncmp(argv[1],"-t",2)!=0){
// the following is for the printout of spins.out ...........................
section4header(fout);
fprintf(fout,"#! %s \n",outstr);
fprintf(fout,"#!T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",T,Hext(1),Hext(2),Hext(3),savmf.na(),savmf.nb(),savmf.nc(),cs4.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
//MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
//                                        to Sa Sb Sc La Lb Lc
//fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] [optional <Sa> <La> <Sb> <Lb> <Sc> <Lc> (hbar)\n");
fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)\n");
fprintf(fout,"#          corresponding exchange fields hxc [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)\n");
// .............................................................................                                
	}
else  //now table options
 { fprintf(fout,"#! nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",savmf.na(),savmf.nb(),savmf.nc(),cs4.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
  fprintf(fout,"# T Ha Hb Hc {sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3]");
 if (strcmp(argv[1],"-tMSL")==0) fprintf(fout,"<Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)\n");
 if (strcmp(argv[1],"-tI")==0) fprintf(fout,"<I1> <I2> <I3> ... <Inofcomponents>\n");
 if (strcmp(argv[1],"-tHex")==0) fprintf(fout,"<Hex1> <Hex2> <Hex3> ... <Inofcomponents> (meV)\n");
 }       
//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: Icalc has to be used to calculate all the <Olm>.
 Vector abc(1,6); abc(1)=1; abc(2)=1; abc(3)=1; // trick to get Habc as components along a,b,c
                  abc(4)=inputpars.cs.alpha(); abc(5)=inputpars.cs.beta(); abc(6)=inputpars.cs.gamma();
 dadbdc2ijk(Hextijk,Hext,abc); // transform Habc=Hext to ijk coordinates ... this is Hextijk
 
Vector h(1,inputpars.cs.nofcomponents);
Vector I(1,inputpars.cs.nofcomponents);
h=0;for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
{(*inputpars.jjj[ii]).Icalc_parameter_storage_init(h,Hextijk,T);} // initialize Icalc module parameter storage

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
{case 1: (*inputpars.jjj[ii]).Scalc(mom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 2: (*inputpars.jjj[ii]).Lcalc(mom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 3: (*inputpars.jjj[ii]).mcalc(mom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage);break;
}

switch(arrow)
{case 1:
 case 2:
 case 3:
         for(nt=1;nt<=3;++nt)
		        {spinconf.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };break;
 case 4: int dim4;
         dim4=3;if((*inputpars.jjj[ii]).module_type==5)
                           dim4=(*(*inputpars.jjj[ii]).clusterpars).cs.nofatoms*3;
         Vector momi(1,dim4);
         (*inputpars.jjj[ii]).micalc(momi,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage);
         // now put the components of momi to the spinconf                                               
         for(nt=1;nt<=dim4;++nt){spinconf.m(i,j,k)(nt+i4-1)=momi(nt);};i4+=dim4;
         break;
}

// output atoms and moments in primitive unit cell to fout  ------------------------------------
if(arrow==4&&(*inputpars.jjj[ii]).module_type==5){
    for(nt=1;nt<=(*(*inputpars.jjj[ii]).clusterpars).cs.nofatoms;++nt)
     {dd3=spinconf.pos_dabc(i,j,k,ii4, cs4);
      dd0=spinconf.pos_dr123(i,j,k,ii4, cs4);
      if (strncmp(argv[1],"-t",2)==0){fprintf(fout,"%4.4f %4.4f %4.4f %4.4f ",T,Hext(1),Hext(2),Hext(3));}
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
    if(true==(*inputpars.jjj[ii]).pcalc(mom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage))
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

 if(true==(*inputpars.jjj[iii]).pcalc(mom,T,hhh,Hextijk,(*inputpars.jjj[iii]).Icalc_parstorage))
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
    if (strncmp(argv[1],"-t",2)==0){fprintf(fout,"%4.4f %4.4f %4.4f %4.4f ",T,Hext(1),Hext(2),Hext(3));}
    fprintf(fout,"{%s} %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f ",
            cs.sipffilenames[ii],dd(1),dd(2),dd(3),dd0(1),dd0(2),dd0(3));
    if (strncmp(argv[1],"-t",2)!=0||strcmp(argv[1],"-tMSL")==0){
    //ouput the magnetic moment if possible
    if((*inputpars.jjj[ii]).mcalc(magmom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage))
    {           for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,magmom(nt)));}
     // and output the orbital and spin momentum if possible 
     if((*inputpars.jjj[ii]).Lcalc(Lmom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage)&&
        (*inputpars.jjj[ii]).Scalc(Smom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage))
      //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
      //                                        to Sa Sb Sc La Lb Lc
      //{        for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f %4.4f",myround(1e-5,Smom(nt)),myround(1e-5,Lmom(nt)));}
      {        for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,Smom(nt)));}
               for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,Lmom(nt)));}
      }
   }}}
if (strncmp(argv[1],"-t",2)!=0||strcmp(argv[1],"-tHex")==0)
  {
    // finally output a line with the exchange fields 
    if (strncmp(argv[1],"-t",2)!=0)fprintf(fout,"\n                 corresponding exchange fields hxc [meV]-->          ");
                      for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout exchangefields
                        {fprintf(fout," %4.4f",myround(1e-5,h(nt)));}
                        
  }
if (strcmp(argv[1],"-tI")==0)
  {(*inputpars.jjj[ii]).Icalc(I,T,h,Hextijk,lnZ,U,(*inputpars.jjj[ii]).Icalc_parstorage);
                         for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout I operator expectation values
                        {fprintf(fout," %4.4f",myround(1e-5,I(nt)));}
                         
  }
fprintf(fout,"\n");

// -----------------------------------------------------------------------------------------------


if(density){
switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c':  (*inputpars.jjj[ii]).chargedensity_coeff (moments, T, h, Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage); break;
 case 's':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsx,1, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsy,2, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsz,3, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=SPINDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'o':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsx,1, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsy,2, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsz,3, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+ORBMOMDENS_EV_DIM)=momentsy(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'm':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsx,1, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsy,2, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsz,3, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            momS=xx*momentsx+yy*momentsy+zz*momentsz;
            if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentlx,1, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momently,2, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentlz,3, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            momL=xx*momentlx+yy*momently+zz*momentlz;
            for(int i=1;i<=SPINDENS_EV_DIM;++i){
            if(doijk==3){moments(i)=momS(i);moments(i+SPINDENS_EV_DIM)=momL(i);}
            else{moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);
                 moments(i+3*SPINDENS_EV_DIM)=momentlx(i);moments(i+4*SPINDENS_EV_DIM)=momently(i);moments(i+5*SPINDENS_EV_DIM)=momentlz(i);
                }
                                               }
            break;
 case 'j':  (*inputpars.jjj[ii]).orbmomdensity_coeff (momentlx,1, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            (*inputpars.jjj[ii]).orbmomdensity_coeff (momently,2, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
            (*inputpars.jjj[ii]).orbmomdensity_coeff (momentlz,3, T, h,Hextijk, (*inputpars.jjj[ii]).Icalc_parstorage);
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
{if(true==(*inputpars.jjj[ii]).pcalc(mom,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage))
{for(nt=1;nt<=3;++nt)
		        {spinconf.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };
}
}

  }}
}}

if (strncmp(argv[1],"-t",2)==0){exit(0);}
  fclose (fout);
   
// create plot of spinconfiguration -----------------------------------------------------------
printf("# ************************************************************************\n");
printf("#%s\n",gp.title);
printf("# ************************************************************************\n");
              if(arrow==0)gp.spins_scale_moment=0;
              if(density==0)gp.show_density=0;

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

    fin = fopen_errchk ("./results/spins.fst", "w");
     spinconf.fst(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,spinconf);
    fclose (fin);

    
   fin = fopen_errchk ("./results/spins_prim.fst", "w");
     spinconf.fstprim(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,spinconf);
    fclose (fin);
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
  if(true==(*inputpars.jjj[ii]).pcalc(pos,T,h,Hextijk,(*inputpars.jjj[ii]).Icalc_parstorage))
 {double charge;charge=(*inputpars.jjj[ii]).charge;if(charge==0)charge=0.01;
  snprintf(cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g",charge,gp.scale_pointcharges*0.529177*signum(charge)*pow((double)fabs(charge),0.3333));
// printf("#! atom %i %s displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,cs.sipffilenames[ii],ii,pos(1),ii,pos(2),ii,pos(3));
  printf("#! atom %i  displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,ii,pos(1),ii,pos(2),ii,pos(3));
 }
}

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin = fopen_errchk ("./results/spins.jvx", "w");
    gp.showprim=0;gp.spins_wave_amplitude=0;
     densitycf.jvx_cd(fin,outstr,cs,gp,0.0,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hextijk,cs4,spinconf,spinconfev_real,spinconfev_imag);
    fclose (fin);

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin = fopen_errchk ("./results/spins_prim.jvx", "w");
     gp.showprim=1;
     densitycf.jvx_cd(fin,outstr,cs,gp,0.0,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hextijk,cs4,spinconf,spinconfev_real,spinconfev_imag);
    fclose (fin);

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
// try a spinwave picture !!!  ... include phonons and spindensity changes ...
//***************************************************************************************************************
//***************************************************************************************************************
if (argc-os>=6){
              // double E;
             long int pos=0;
             int extended_eigenvector_dimension;
              char instr[MAXNOFCHARINLINE];
              float numbers[20];numbers[9]=1;numbers[10]=3;
              numbers[0]=20;
             gp.spins_wave_amplitude=1.0;gp.spins_show_ellipses=1.0;gp.spins_show_oscillation=1.0;
             gp.phonon_wave_amplitude=1.0;gp.phonon_scale_static_displacements=1.0;
//----------------------------------------------------------------------------------------------------------
            if(arrow>0){double checkdd=1e7;
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
             }
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_oscillation",gp.spins_show_oscillation);
              }

               j=fseek(fin,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               double delta,dd,ddT,ddHa,ddHb,ddHc,ddh,ddk,ddl,ddE;
               for (delta=1000.0;feof(fin)==0                      //end of file
                     &&(n=inputline(fin,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin); 
                 spincf ev_real(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,3);
                 spincf ev_imag(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,3);
                 ev_real.load(fin);ev_imag.load(fin);
                 ddT=T-numbers[4];ddT*=ddT;
                 ddHa=Hext(1)-numbers[1];ddHa*=ddHa;
                 ddHb=Hext(2)-numbers[2];ddHb*=ddHb;
                 ddHc=Hext(3)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5+os],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6+os],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7+os],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8+os],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;checkdd=fabs(T-numbers[4])+fabs(Hext(1)-numbers[1])+fabs(Hext(2)-numbers[2])+fabs(Hext(3)-numbers[3]);
                  snprintf(outstr,sizeof(outstr),"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];//E=numbers[9]; 
                  spinconfev_real=ev_real;
                  spinconfev_imag=ev_imag;                  
                 }
                 pos=ftell(fin); 
                 fgets(instr,MAXNOFCHARINLINE,fin); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin);fgets(instr,MAXNOFCHARINLINE,fin);}
                 j=fseek(fin,pos,SEEK_SET);
               }
              fclose (fin);
             if(checkdd>1e-7)
               {fprintf(stderr,"Error program spins - inconsistent output files mcphas.mf and mcdisp.*:  temperature/magnetic field in static configuration\n" 
                               "(from results/mcphas.mf) is different from mcdisp output files results/mcdisp.qes,qem,qeo:\n %s\n"
                               "... probably you need to rerun setup_mcdisp_mf and mcdisp with T=%g Ha=%g Hb=%g Hc=%g",outstr,T,Hext(1),Hext(2),Hext(3));exit(EXIT_FAILURE);}
               fprintf(stdout,"#%s - moment oscillation - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              spinconfev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              spinconfev_imag.print(stdout);
             }//arrow
     
//----------------------------------------------------------------------------------------------------------
            if(phonon>0){double checkdd=1e7;
             strcpy(infilename,"./results/");strcpy(infilename+10,prefix);
             strcpy(infilename+10+strlen(prefix),"mcdisp.qep");fin = fopen(infilename, "rb");
                      if(fin==NULL)fin = fopen_errchk ("./results/mcdisp.qep", "rb");
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"phonon_wave_amplitude",gp.phonon_wave_amplitude);
               extract(instr,"phonon_scale_static_displacements",gp.phonon_scale_static_displacements);
              }
               j=fseek(fin,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qep file format\n");exit (EXIT_FAILURE);}
   
               double delta,dd,ddT,ddHa,ddHb,ddHc,ddh,ddk,ddl,ddE;
               for (delta=1000.0;feof(fin)==0                      //end of file
                     &&(n=inputline(fin,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin); 
                 spincf ev_real(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,3);
                 spincf ev_imag(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,3);
                 ev_real.load(fin);ev_imag.load(fin);
                 ddT=T-numbers[4];ddT*=ddT;
                 ddHa=Hext(1)-numbers[1];ddHa*=ddHa;
                 ddHb=Hext(2)-numbers[2];ddHb*=ddHb;
                 ddHc=Hext(3)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5+os],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6+os],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7+os],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8+os],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;checkdd=fabs(T-numbers[4])+fabs(Hext(1)-numbers[1])+fabs(Hext(2)-numbers[2])+fabs(Hext(3)-numbers[3]);
                  snprintf(outstr,sizeof(outstr),"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];//E=numbers[9]; 
                  spinconfpev_real=ev_real;
                  spinconfpev_imag=ev_imag;                  
                 }
                 pos=ftell(fin); 
                 fgets(instr,MAXNOFCHARINLINE,fin); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin);fgets(instr,MAXNOFCHARINLINE,fin);}
                 j=fseek(fin,pos,SEEK_SET);
               }
              fclose (fin);
             if(checkdd>1e-7)
               {fprintf(stderr,"Error program spins - inconsistent output files mcphas.mf and mcdisp.*:  temperature/magnetic field in static configuration\n" 
                               "(from results/mcphas.mf) is different from mcdisp output files results/mcdisp.qep:\n %s\n"
                               "... probably you need to rerun setup_mcdisp_mf and mcdisp with T=%g Ha=%g Hb=%g Hc=%g",outstr,T,Hext(1),Hext(2),Hext(3));exit(EXIT_FAILURE);}
               fprintf(stdout,"#%s - phonon oscillation - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              spinconfpev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              spinconfpev_imag.print(stdout);
             }//phonon
//----------------------------------------------------------------------------------------------------------
            if(density){double checkdd=1e7;
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
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_oscillation",gp.spins_show_oscillation);
               extract(instr,"extended_eigenvector_dimension",extended_eigenvector_dimension);
              }

               j=fseek(fin,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               double delta,dd,ddT,ddHa,ddHb,ddHc,ddh,ddk,ddl,ddE;
               for (delta=1000.0;feof(fin)==0                      //end of file
                    &&(n=inputline(fin,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin); 
                 spincf ev_real(densitycf.na(),densitycf.nb(),densitycf.nc(),densitycf.nofatoms,extended_eigenvector_dimension);
                 spincf ev_imag(densitycf.na(),densitycf.nb(),densitycf.nc(),densitycf.nofatoms,extended_eigenvector_dimension);
                 ev_real.load(fin);ev_imag.load(fin);
                 ddT=T-numbers[4];ddT*=ddT;
                 ddHa=Hext(1)-numbers[1];ddHa*=ddHa;
                 ddHb=Hext(2)-numbers[2];ddHb*=ddHb;
                 ddHc=Hext(3)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5+os],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6+os],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7+os],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8+os],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;checkdd=fabs(T-numbers[4])+fabs(Hext(1)-numbers[1])+fabs(Hext(2)-numbers[2])+fabs(Hext(3)-numbers[3]);
                  snprintf(outstr,sizeof(outstr),"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];//E=numbers[9]; 

            switch(argv[1][1]) // dimension definition from jjjpar.hpp
            {case 's': 
             case 'o': 
                      if(doijk==3){// moments=xx*momentsx+yy*momentsy+zz*momentsz;
                   for(ii=1;ii<=inputpars.cs.nofatoms;++ii)
                   for (i=1;i<=savmf.na();++i)for(j=1;j<=savmf.nb();++j)for(k=1;k<=savmf.nc();++k)
                  for(nt=1;nt<=dim;++nt){densityev_real.m(i,j,k)(nt+dim*(ii-1))=xx*ev_real.m(i,j,k)(nt+3*dim*(ii-1))+yy*ev_real.m(i,j,k)(nt+dim+3*dim*(ii-1))+zz*ev_real.m(i,j,k)(nt+2*dim+3*dim*(ii-1));
                                         densityev_imag.m(i,j,k)(nt+dim*(ii-1))=xx*ev_imag.m(i,j,k)(nt+3*dim*(ii-1))+yy*ev_imag.m(i,j,k)(nt+dim+3*dim*(ii-1))+zz*ev_imag.m(i,j,k)(nt+2*dim+3*dim*(ii-1));
                                        }
                      }
                      else{densityev_real=ev_real;
                           densityev_imag=ev_imag;
                     }
                       break;
            case 'c': 
            case 'j': densityev_real=ev_real;
                      densityev_imag=ev_imag;
                     break;
            default: help_and_exit();
                     }
                    }
               }
              fclose (fin);
             if(checkdd>1e-7)
               {fprintf(stderr,"Error program spins - inconsistent output files mcphas.mf and mcdisp.*:  temperature/magnetic field in static configuration\n" 
                               "(from results/mcphas.mf) is different from mcdisp output files results/mcdisp.qee,qsd,qod:\n %s\n"
                               "... probably you need to rerun setup_mcdisp_mf and mcdisp with T=%g Ha=%g Hb=%g Hc=%g",outstr,T,Hext(1),Hext(2),Hext(3));exit(EXIT_FAILURE);}
              fprintf(stdout,"#%s - density oscillation - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              densityev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              densityev_imag.print(stdout);
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
               snprintf(filename,sizeof(filename),"./results/spins.%i.jvx",i+1);
               fin = fopen_errchk (filename, "w");gp.showprim=0;
                     densitycf.jvx_cd(fin,outstr,cs,gp,
                                  phase,densityev_real,densityev_imag,hkl,T,hh,Hextijk,cs4,spinconf,spinconfev_real,spinconfev_imag,spinconfpev_real,spinconfpev_imag);
               fclose (fin);
               snprintf(filename,sizeof(filename),"./results/spins_prim.%i.jvx",i+1);
               fin = fopen_errchk (filename, "w");gp.showprim=1;
                     densitycf.jvx_cd(fin,outstr,cs,gp,
                                  phase,densityev_real,densityev_imag,hkl,T,hh,Hextijk,cs4,spinconf,spinconfev_real,spinconfev_imag,spinconfpev_real,spinconfpev_imag);
               fclose (fin);
              }
          printf("# %s\n",outstr);
          }
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spins\n");
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


