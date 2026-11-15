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
#include <memory>
#include <vector>

#define MAXDEG 100    // maximal degeneracy
#define MAXNOFSPINCF 10000000 // maximal number of spinconfigurations to be loaded
#define MAXNOFATOMS  10000

enum class OperationMode {
    Graphical,
    File,
    Table,
    Eigenvector
};

struct LoadedInput {
    OperationMode mode = OperationMode::Graphical;

    FILE* fin = nullptr;
    FILE* fout = nullptr;

    int minl = 1, maxl = 1, eps = 0, fst = 0,dim=28;
    int nomorethanthis = MAXNOFSPINCF;
    int nofcomp = 1000000000; 
    int nofat = MAXNOFATOMS; // nofat: for  -fc   it contains maximum number of atoms
                             //        for  -t* -n it contains atom number to be printed
    int os = 0, maxn = 0, doijk = 0, arrow = 0, density = 0, phonon = 0;
    
    double T=0; 
    double xx = 0, yy = 0, zz = 0, limit = 0, tL2 = 0;

    bool longtable = false;

    char prefix[MAXNOFCHARINLINE]{};
    char outstr[MAXNOFCHARINLINE]{};
    char infilename[MAXNOFCHARINLINE]{};
    double aa[NOF_USERDEF_MCPHAS_COLS+1]{};

    cryststruct cs;
    cryststruct cs4;
    graphic_parameters gp;
    Matrix N{1, 3, 1, 3};// demagnetisation factor
    Matrix p{1, 3, 1, 3};// primitive lattice
    Vector nmin{1, 3}, nmax{1, 3};
    Vector Hext{1,HEXT_DIMENSION};
    // These are constructed after command-line parsing / structure loading.
    std::unique_ptr<par> inputpars;
    std::unique_ptr<spincf> savmf;

    // Own the buffers, and expose their addresses as char* values for
    // functions that expect the existing char** API.
    std::vector<std::unique_ptr<char[]>> outStorage;
    std::vector<char*> out;
};

bool detectOperationMode(int argc, char** argv, LoadedInput& loaded)
{
    loaded.mode = OperationMode::Graphical; // default

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-q") == 0) {
            loaded.mode = OperationMode::Eigenvector;
            return true;
        }
    }

    // Handle other modes using the program's actual command syntax.
    if (argc > 1 && strncmp(argv[1], "-t", 2) == 0) {
        loaded.mode = OperationMode::Table;
    } else if (argc > 1 && strcmp(argv[1], "-f") == 0) {
        loaded.mode = OperationMode::File;
    }

    return true;
}


void help()
    { printf ("\n\
program spins - popout spin/exchange field configuration\n"
"              - and/or display 3d animation of spin/moment/densities and animations\n\n\
use as\n\
  1) file-mode\n\
        spins -f[c 1 13 3 0.1] [-n 2] mcphas.sps T Ha Hb Hc\n\
    or: spins -f[c 1 13 3 0.2] [-n 2] mcphas.sps x y\n\
    or: spins -f[c 1 13 3 0.2] [-n 2] mcphas.sps out1 out2 out3 out4 out5 out6 out7\n\
    or: spins -f[c 1 13 3 0.1] [-n 2] mcphas.tst n\n\
  2) table-mode \n\
    or: spins -tMSL [-n 2] [-prefix 001] [T Ha Hb Hc] \n\
    or: spins -tHex [-n 2] [-prefix 001] [x y] \n\
    or: spins -tI  [-n 2] [-prefix 001] [out1 out2 out3 out4 out5 out6 out7] \n\
    or: spins -tL2 4.5 [-n 2] [-prefix 001] [T Ha Hb Hc] \n\
  3) graphic-mode\n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M|-pel] [-P] [-eps|-fst] [-prefix 001] T Ha Hb Hc [ h k l E]\n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M|-pel] [-P] [-eps|-fst] [-prefix 001] x y [h k l E] \n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M|-pel] [-P] [-eps|-fst] [-prefix 001] out1 out2 out3 out4 out5 out6 out7 [[-prefix 020] h k l E]\n\
                    \n\
  4) eigenvector-mode\n\
    or: spins [-symop|-symoa [-]3 cryststruc.cif] -q mcdisp.qep h k l [E]\n\
                    \n\
1) file-mode:\n\
   a) if used with -f file T Ha Hb Hc, this file has to be a mcphas.mf or mcphas.sps file,\n\
   the spin configuration at given temperature T[K] and magnetic effective field H[T]\n\
   is read and extracted from this file and printed on screen (stdout),\n\
   results/spins.out is created (with mag moment chosen to be = <Ia> <Ib> <Ic>),\n\
   a simple graphics to represent the configuration is created in results/spins_prim.jvx \n\
   Note: - T Ha Hb Hc stand for the 3rd 5th 6th 7th column in mcphas.* output files\n\
         (may have different meaning if out3 out5 out6 out7 is set in mcphas.ini)\n\
         - out1,...,out7 refers to column 1-7 in mcphas.* output files\n\
  \n\
   b) if used with -f file x y, then this file has to be a mcphas.mf or mcphas.sps file,\n\
   the spin configuration at a given x,y point is read and extracted from this file ,\n\
   and printed on screen (stdout) etc. as 1)\n\
   Note: x y stands for the 1st and 2nd column in mcphas.* output files\n\
         (may have different meaning if out1 or out2 is set in mcphas.ini)\n\
   c) if used with -f filen n,  this file has to be a mcphas.tst file,\n\
   the spin configuration number n\n\
   is read and extracted from this file and printed on screen (stdout),\n\
   results/spins.out is created (with mag moment chosen to be = <Ia> <Ib> <Ic>)\n\
   a simple graphics to represent the configuration is created in results/spins_prim.jvx \n\
  \n\
a&b&c) - if used with -fc min max n lim a human readable format is output for spin components with index \n\
   from min to max, only n numbers exceeding  \n\
   absolute value of lim, e.g. -fc 1 13 3 0.1 outputs I1,I2,...,I13, however only at maximum 3 numbers \n\
   for each atom which are all larger (absolute value) than 0.1. \n\
      - if option -n 5 ist present and not -fc only I1,..,I5 are output for each atom (mind: number 5<=nofcomponents)\n\
      - if option -n 5 and option -fc are present - information is output only for atoms 1,2,3,4,5 (mind in this case 5<=nofatoms \n\
 \n\
2) table mode \n\
   in this mode spins is used without a filename, the information is read from \n\
   results/mcphas.*  and output is created for the configuration which matches \n\
   best the external parameters x y or T Ha Hb Hc or out1-ou7 specified in the command line.\n\
   If no external parameters are given, a table is created containing all configurations which \n\
   are available in the output files of mcphas/mcdisp. \n\
   for table the options are:\n\
       -tMSL ... output to stdout a table with T Ha Hb Hc atom positions and with \n\
              magnetic moments <Mx> <My> <Mz>, orbital moments and spin \n\
              of each atom in the magnetic unitc cell \n\
       -tHex ...  output to stdout table with T Ha Hb Hc atom positions and exchange fields Hex\n\
       -tI   ... a similar table with expectation values of interaction operators <I>\n\
       -tL2 4.5  ... a table with average squared longitudinal bond elongations up to bond\n\
                     length 4.5 Angstroem, only atoms with phonon degrees of freedom are considered\n\
                     (useful for estimating elastic energy contributions of different bonds)\n\
       -n 3   ....  only print out for atom number 3 the relevant information, omit other atoms \n\
3) graphic / excitation mode \n\
   also in this mode spins is used without a filename, the information is read from \n\
   results/mcphas.* results/mcdisp.* , for graphical animations the options are: \n\
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
         Note, that in order to animate changes in the above quantities, the corresponding\n\
         switch has to be enabled in the mcdisp calculation (mcdisp.par) and the single ion\n\
         modules have to be capable of calculating the corresponding observables. In this\n\
         case mcdisp produces files results/mcdisp.qev .pev .sev etc. which contain the \n\
         correlation function for the observable (for displacements this corresponds to \n\
         the phonon eigenvector for phonon displacements in results/mcdisp.qep; for magnetic\n\
         moments this corresponds to the spinwave eigenvector  in results/mcdisp.qem, ... )\n\
\n\
         -prefix 001 ... use input file(s) results/001mc* instead of results/mc*\n\
 \n\
4) eigenvector mode:\n\
         also in this mode spins is used without a filename, the information is read from \n\
         results/mcphas.* results/mcdisp.*  for files containing the eigenvectors corresponding to correlation functions of\n\
         observables (see manual Evaluating results - Graphical Visualisation of Spin and Charge \n\
         Oscillations of Excitations) can be evaluated using the eigenvector mode.\n\
         The filename corresponding to the observable (e.g. for magnetic moment, phonon displacement \n\
         it is mcdisp.qem, mcdisp.qep, respectively) and the propagation\n\
         vector q=(hkl) and the energy have to be given and the program looks for the \n\
         closest match and prints out the eigenvector as stored in the file and also \n\
         in human readable form.\n\
         Using option -symop the symmetry operations in a .cif file under the tag\n\
          loop_\n\
          _space_group_symop_id\n\
          _space_group_symop_operation_xyz\n\
          1 'y+1/4, x+3/4, -z+1/4'\n\
          2 'y+1/4, -x+3/4, -z+1/4'\n\
         are read and spins outputs the eigenvector after applying the specified symmetry operation assuming\n\
         that the observable is a polar vector, e.g. phonon displacememnts. Use option -symoa for axial \n\
         vectors such as the magnetic moment.\n\
         If a negative sign is specified for the symmetry operation (e.g. -3), all symmetry operations\n\
         are employed and characters are output (provided that the eigenvector is invariant).\n\
         If in the .cif file the character table of the point group is provided, then \n\
         the irreducible representation is identified. The .cif file has to contain a character table\n\
         as follows:\n\
         _point_group_name D4h\n\
         _point_group_classes                   'E  EI C42 C42I 2C4 2C4I 2C2 2C2I 2C2d 2C2dI '\n\
         _point_group_space_symmetry_operations '16 8  11  3    13  5    1   9    6    14'\n\
         loop_\n\
         _irrep_point_group\n\
         _character_point_group\n\
         A1g ' 1  1  1  1  1  1  1  1  1  1'\n\
         A1u ' 1 -1  1 -1  1 -1  1 -1  1 -1'\n\
         If energy is omitted, then all energies in the qep file for the wave vectors are\n\
         scanned and the characters are output and  the  irrep are determined.\n\
\n\
\n\
     examples:\n\
        spins -c 2 0 0 1\n\
        ...calculates the charge density at T=2K and H=(0,0,1) Tesla and produces jvx file to view\n\
        spins -tI 2 0 0 1 \n\
        ... outputs a table with atomic positions and expectations values <I> \n\
        spins -M 2 0 0 1   0.5 0  0  23\n\
        ...calculates the magnetic moment precession  at T=2K and H=(0,0,1) Tesla for an excitation at\n\
          (hkl)=(0.5 0 0) at 23 meV - produces a jvx file with animation\n\
        \n\
 This program outputs a magnetic structure (and magnetic excitation)\n\
 graphic/movie in the output files of different format:\n\
 results/spins*.eps (postscript), results/spins*.fst (fp_studio), \n\
 results/spins.out (ascii) and results/spins*.jvx (javaview)\n\n\
 the graphics output format can be fine tuned in .sps and .qev input files\n"
" or results/graphic_parameters.set by show_abc_unitcell,\n"
" show_primitive_crystal_unitcell, spins_scale_moment, spins_wave_amplitude\n\
 show_magnetic_unitcell, show_atoms, scale_view_1,scale_view_2, scale_view_3 ...\n\n\
 jvx files can be viewed by:\n\
 java JvxAnimator results/spins.jvx \n\
 java JvxAnimator \"results/spins.*.jvx\" \n"
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
                  
void print_eigenvector_human_legible(par & p,spincf & ev_r,spincf & ev_i,char * outhklstr)
{printf("%s\n#sipf-file at (da[a] db[b] dc[c]) eigenvector\n",outhklstr);
 
  for(int ii=1;ii<=p.cs.nofatoms;++ii)
   for(int i=1;i<=ev_r.na();++i)
   for(int j=1;j<=ev_r.nb();++j)
   for(int k=1;k<=ev_r.nc();++k) //loop supercells
 { Vector dd3(1,3),dd4(1,3);
    dd3=ev_r.pos_dabc(i,j,k,ii, p.cs);
    dd4=ev_r.pos_dr123(i,j,k,ii, p.cs);
    printf("# %s at dabc=(%5.3f %5.3f %5.3f) dr1r2r3=(%5.3f %5.3f %5.3f) [ ",(*p.jjj[ii]).sipffilename,dd3(1),dd3(2),dd3(3),dd4(1),dd4(2),dd4(3));
    for(int a=1;a<=ev_r.nofcomponents;++a) {if( ev_r.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)!=0.0||
                                                ev_i.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)==0.0)
                                            printf("%5.4f",ev_r.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)); 
                                           if( ev_i.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)!=0.0)
                                            printf("%+5.4fi",ev_i.m(i,j,k)((ii-1)*ev_i.nofcomponents+a));
                                           printf(" ");
                                         }     
    printf("]\n");  
 }  
}

// apply symmetry operation sn from ciffile  and return character
// if verbose=true: print eigenvectors
complex<double> symop(par & p,spincf * scr,spincf * sci,int deg,char * outhklstr,char * sn,char * ciffile, bool polar,bool verbose)
{FILE * fin; fin = fopen_errchk (ciffile,"rb");
 char instr[MAXNOFCHARINLINE];bool found=false,begin=false;
 while(feof(fin)==false&&found==false)
 {fgets(instr, MAXNOFCHARINLINE, fin);
  if(instr[strspn(instr," \t")]!='#')
  {if(begin==true)
    {if(strncmp(sn,instr,strlen(sn))==0)found=true;
    }
   if(strncmp(instr,"_space_group_symop_operation_xyz",32)==0)begin=true;
  }
 }
 fclose(fin);
if(found==false){fprintf(stderr,"Error program spins: symmetry operation %s not found\n",sn);exit(EXIT_FAILURE); }

if(verbose)printf("#Applying symmetry operation %s\n",instr);
  
parser ob;
char * xd=mystrnstr(instr,"'",MAXNOFCHARINLINE);++xd;
char * yd=mystrnstr(xd,",",MAXNOFCHARINLINE);*yd='\0';++yd;
char * zd=mystrnstr(yd,",",MAXNOFCHARINLINE);*zd='\0';++zd;
char * d=mystrnstr(zd,"'",MAXNOFCHARINLINE);*d='\0';
//printf("%s\n%s\n%s\n",xd,yd,zd);
char f[MAXNOFCHARINLINE];
Vector t(1,3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",0.0,0.0,0.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",t(1),ob);extract(f,"yd",t(2),ob);extract(f,"zd",t(3),ob);
if(verbose)myPrintVector(t,"corresponding to x'=A*x+t and ev'=B*ev\n# translation t:");
Matrix A(1,3,1,3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",1.0,0.0,0.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",A(1,1),ob);extract(f,"yd",A(2,1),ob);extract(f,"zd",A(3,1),ob);
                A(1,1)-=t(1);A(2,1)-=t(2);A(3,1)-=t(3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",0.0,1.0,0.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",A(1,2),ob);extract(f,"yd",A(2,2),ob);extract(f,"zd",A(3,2),ob);
                A(1,2)-=t(1);A(2,2)-=t(2);A(3,2)-=t(3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",0.0,0.0,1.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",A(1,3),ob);extract(f,"yd",A(2,3),ob);extract(f,"zd",A(3,3),ob);
                A(1,3)-=t(1);A(2,3)-=t(2);A(3,3)-=t(3);
if(verbose)myPrintMatrix(A,"Roto(inversion) Matrix A:");
// determine transformation matrix for vectors (and pseudovectors) from
// spatial transformation matrix A: B should transform eigenvectors like a spatial vector
// in ijk coordinates.
// dabcd=A*dabc  
// dijk=abc_in_ijk*dabc 
// dijkd=abc_in_ijk*dabcd 
//  --> dijkd=abc_in_ijk*A*abc_in_ijk.Inverse()*dijk
// -->
Matrix B(1,3,1,3);Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,p.cs.abc);
B=abc_in_ijk*A*abc_in_ijk.Inverse(); 
if(polar==false&&Det(B)<0){B=-B;} // for axial vectors ignore the inversion symmetry operation
if(verbose)myPrintMatrix(B,"Observable transformation Matrix B:");

// go through all degenerate configurations, transform them and store original and transformed
// vectors in matrices evm and evm'
int dim=scr[1].nofcomponents*scr[1].nofatoms*scr[1].n();
ComplexMatrix evm(1,dim,1,deg),evmd(1,dim,1,deg);
for(int dg=1;dg<=deg;++dg)
{if(scr[dg].nofcomponents!=3){fprintf(stderr,"Error spins: nofcomponents=%i symmetry operations can currently only be performed for vector observables with nofcomponents=3\n",scr[dg].nofcomponents);exit(EXIT_FAILURE);}
spincf evd_r=scr[dg].symop(A,t,B,p.cs);
spincf evd_i=sci[dg].symop(A,t,B,p.cs);
if(verbose){printf("# Original eigenvector:\n");print_eigenvector_human_legible(p,scr[dg],sci[dg],outhklstr);}
if(verbose){printf("# Transformed eigenvector:\n");print_eigenvector_human_legible(p,evd_r,evd_i,outhklstr);}
// store as vector in order to determine character of representation ...
Vector vr=scr[dg].asVector();
Vector vi=sci[dg].asVector();
ComplexVector v(vr,vi);for(int i=v.Lo();i<=v.Hi();++i)evm(i,dg)=v(i);

Vector vdr=evd_r.asVector();
Vector vdi=evd_i.asVector();
ComplexVector vd(vdr,vdi);for(int i=vd.Lo();i<=vd.Hi();++i)evmd(i,dg)=vd(i);
}
// linear regression to determine deg x deg transformation matrix T
// evmd=evm*T
// evm.Transpose * evmd= evm.Transpose * evm * T = C * T
// T=C.inverse * evm.Transpose * evmd
ComplexMatrix C(1,deg,1,deg); C=evm.Hermitean()*evm;
ComplexMatrix T(1,deg,1,deg); T=C.Inverse()*(evm.Hermitean()*evmd);

complex<double> ch=Trace(T);
// detect how good is linear regression in order to see if something is odd
// evmd=evm*T should be tue. 
// evm*T-evmd=0 should be true - check it
double e=NormFro(evm*T-evmd)/dim/deg;
if(e>0.0001) {myPrintMatrix(evmd,"Matrix of transformed eigenvectors:");evmd=evm*T-evmd;
             myPrintMatrix(evmd,"Difference Matrix to Matrix of best fitting linear combination of original eigenvectors:");
             fprintf(stderr,"Error: spin configuration not invariant to symmetry operation %s |s-s'|=%g\n",instr,e);
             exit(EXIT_FAILURE); }

if(verbose)printf("ch%s=%5.3f%+5.3fi ",sn,ch.real(),ch.imag());
return ch;
}


bool identify_irrep(complex<double> * ch,int nofsymop,char * ciffile)
// look if cif file contains point group character table and if yes find out irrep
//_point_group_name D4h
//_point_group_classes                   'E  EI C42 C42I 2C4 2C4I 2C2 2C2I 2C2d 2C2dI '
//_point_group_space_symmetry_operations '16 8  11  3    13  5    1   9    6    14'
//loop_
//_irrep_point_group
//_character_point_group
// A1g ' 1  1  1  1  1  1  1  1  1  1'
// A1u ' 1 -1  1 -1  1 -1  1 -1  1 -1'
{FILE * fin; fin = fopen_errchk (ciffile,"rb");
 char instr[MAXNOFCHARINLINE];
 char pointgroup[MAXNOFCHARINLINE];
 char classes[MAXNOFCHARINLINE];
 float pso[nofsymop+1];pso[0]=nofsymop;
 bool found=false,begin=false,eot=false;int nofc=0,nd=0;char * irrep;
 while(feof(fin)==false&&found==false&&eot==false)
 {fgets(instr, MAXNOFCHARINLINE, fin);
  if(instr[strspn(instr," \t")]!='#')
  {char * t=instr+strspn(instr," \t");
   if(strncmp(t,"_point_group_name",17)==0)snprintf(pointgroup,MAXNOFCHARINLINE,"%s",t+18);
   if(strncmp(t,"_point_group_classes",20)==0)snprintf(classes,MAXNOFCHARINLINE,"%s",t+20);
   if(strncmp(t,"_point_group_space_symmetry_operations",38)==0)
       {t+=38;t+=strspn(t," \t'");
        if((nofc=splitstring(t,pso))>nofsymop){fprintf(stderr,"Error spins reading %s: _point_group_space_symmetry_operations > number of symmetry operations read %i\n",ciffile,nofsymop);exit(EXIT_FAILURE); }
       }
   if(begin==true) // we are in the table - compare irrep characters with irreps
    {if(instr[strspn(instr," \t")]=='_'||strncmp(instr+strspn(instr," \t"),"loop_",5)==0)eot=true;
     else{
     irrep=1+mystrnstr(instr,"'",MAXNOFCHARINLINE);
     float irrepch[nofsymop+1];irrepch[0]=nofsymop;
     int nofd=splitstring(irrep,irrepch);
     if(nofc!=nofd){fprintf(stderr,"Error - irrep contains %i characters not equal to %i in tag _point_group_space_symmetry_operations \n",nofc,nofd);exit(EXIT_FAILURE);}
     int i=0;
     while(i<nofc&&fabs(irrepch[i+1]-ch[(int)pso[i+1]].real())<1e-3){++i;//printf("%8.6f=%8.6f ",irrepch[i],ch[(int)pso[i]].real());
      }
     //printf("\n");
     if(i==nofc)found=true;
     //else printf("i=%i table: %8.6f  calc: %8.6f\n",(int)pso[i+1],irrepch[i+1],ch[(int)pso[i+1]].real());
     }}
   if(strncmp(instr,"_character_point_group",22)==0)begin=true;
   }
 ++nd;
 }
 fclose(fin);
 if(found==true){printf("irrep=%s",instr);return true;}
 else {printf("irrep not found in %s\n",ciffile);return false;}
}
                   

/**********************************************************************/
// programm fuer option -q (koennte auch eigenstaendig compiliert werden)
/**********************************************************************/
int spins_q(int argc,char **argv)
{FILE * fin;int os=0;
 
 graphic_parameters gp;bool deg=true;
 char prefix[MAXNOFCHARINLINE];prefix[0]='\0';
 par inputpars("mcphas.j",prefix); 
 spincf spinconf(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 spincf spinconfev_real(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 spincf spinconfev_imag(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 Vector Hext(1,HEXT_DIMENSION);Hext=0;
 char outstr[MAXNOFCHARINLINE];
 char outhklstr[MAXNOFCHARINLINE];
 char * ciffile;ciffile=NULL;
 Vector hkl(1,3);double E;double T=0;char * sn=NULL;
 bool polar=true;
if(strncmp(argv[os+1],"-symo",5)==0){if(argv[os+1][5]=='a')polar=false;
                                     sn=argv[os+2];// label of symoperation in table
                                     ciffile=argv[os+3];
                                     os+=3;// here we could enter some function to select desired symop from file                                  
                                    } 
os-=2;
spincf  scr[MAXDEG];
spincf  sci[MAXDEG];
double Etable[MAXDEG];
Etable[0]=1;
if(argc==os+8) // no energy is given - probe all energies
{fin = fopen_errchk (argv[os+4], "rb");
 char instr[MAXNOFCHARINLINE];instr[0]='#';
 while (feof(fin)==0)
 {fgets(instr,MAXNOFCHARINLINE,fin);
   if (instr[strspn(instr," \t")]!='#') // pointer to 'ltrimstring' 
    {if(0==extract(instr,"E",Etable[(int)Etable[0]]))Etable[0]++;
    }
 }Etable[0]--;
 fclose(fin);
}
printf("#* program spins ... reading from file %s\n",argv[os+4]); 
for(int ie=1;ie<=Etable[0];++ie)
{E=Etable[ie];T=0;deg=true;
while(deg){fin = fopen_errchk (argv[os+4], "rb");  
 deg=check_for_best_excitation_and_close(fin,gp,"ddd",spinconf,spinconfev_real,spinconfev_imag,
                          0,0,T,Hext,outstr,outhklstr,"best",hkl,E,argc,argv,os,
                          inputpars.cs.abc,inputpars.cs.nofatoms,0,0.0,0.0,0.0,0); 
 T--;
 if(argc>os+8)print_eigenvector_human_legible(inputpars,spinconfev_real,spinconfev_imag,outhklstr);
 scr[-(int)T]=spincf(spinconfev_real);
 sci[-(int)T]=spincf(spinconfev_imag);
 }

 if(ciffile!=NULL){if(sn[0]=='-')
                   {char snn[10];int nofsymop=(int)strtod(sn+1,NULL);
                    complex<double> ch[nofsymop+1];
                    for(int i=1;i<=nofsymop;++i)
                    {snprintf(snn,10,"%i",i);
                     ch[i]=symop(inputpars,scr,sci,(int)-T,outhklstr,snn,ciffile,polar,false);
                    }
                    printf("E=%g ",E);
                   // look if cif file contains point group character table and if yes find out irrep
                   if(identify_irrep(ch,nofsymop,ciffile)==false)
                    {for(int i=1;i<=nofsymop;++i)printf("ch%i=%5.3f%+5.3fi ",i,ch[i].real(),ch[i].imag());
                     printf("\n");
                    }
                   }
                   else symop(inputpars,scr,sci,(int)-T,outhklstr,sn,ciffile,polar,true);
                  }
} 
 return 0; 
}
          
// to be done for table and graphical mode 
void getinfilename(LoadedInput& s, int argc, char** argv)
{if(argc>s.os+1)if(strcmp(argv[1+s.os],"-prefix")==0){strcpy(s.prefix,argv[2+s.os]); // read s.prefix
                                   fprintf(stdout,"# prefix for input filenames: %s\n",s.prefix);
 				   s.os+=2;}
 

 strcpy(s.infilename,"./results/");strcpy(s.infilename+10,s.prefix);
 strcpy(s.infilename+10+strlen(s.prefix),"mcphas.mf");s.fin = fopen(s.infilename, "rb");
 if(s.fin==NULL){strcpy(s.infilename+10,"mcphas.mf");s.fin = fopen_errchk(s.infilename, "rb");}
 printf("# reading from file %s\n",s.infilename);
 
 }


bool parse_options_and_load_structure(
    int argc,
    char** argv,
    LoadedInput& s)
{
    if (argc < 2) {
        help();
        return false;
    }
fprintf(stderr,"# ***********************************************************\n");
fprintf(stderr,"# * spins - analyse mcphas output and display 3d graphics of*\n");
fprintf(stderr,"# *  spins,moments,densities, etc at given H and T          *\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R              *\n");
fprintf(stderr,"# * %s                                      *\n",MCPHASVERSION);
fprintf(stderr,"# ***********************************************************\n");
 


    // Move the command-line parsing block here.
    // Change local variables to state members, for example:
    //   fin       -> s.fin
    //   gp        -> s.gp
    //   prefix    -> s.prefix
    //   nofat     -> s.nofat
    //   os        -> s.os
    //
    // For each error that currently calls exit(...), report it and return
    // false instead. Example:
    //
    //   if (too_few_arguments) {
    //       fprintf(stderr, "Error: too few arguments\n");
    //       return false;
    //   }
// check command line
detectOperationMode(argc,argv,s);

s.prefix[0]='\0';
// options without graphics just screendump <I> or exchange field configuration at given HT
switch (s.mode) {
case OperationMode::Eigenvector:
             if (strncmp(argv[s.os+1],"-symo",5)==0) {s.os+=3;// here we could enter some function to select desired symop from file
                                    }
             if(argc>5+s.os){return spins_q(argc,argv);}else{fprintf(stderr,"Error spins: too few arguments for opiont -q\n");return false;}
                          break;
case OperationMode::File:
  s.os+=2;if (strcmp(argv[s.os-1],"-fc")==0){s.minl=(int)strtod(argv[s.os],NULL);s.maxl=(int)strtod(argv[s.os+1],NULL);
                                      s.maxn=(int)strtod(argv[s.os+2],NULL);s.limit=strtod(argv[s.os+3],NULL);
                                     s.os+=4;}
   if (strncmp(argv[s.os],"-n",2)==0){if(strcmp(argv[1],"-fc")==0)s.nofat=(int)strtod(argv[s.os+1],NULL);
                                    else s.nofcomp=(int)strtod(argv[s.os+1],NULL);
                                    s.os+=2;}
   s.fin = fopen_errchk (argv[s.os], "rb");strcpy(s.infilename,argv[s.os]);printf("#* program spins ... reading from file %s\n",s.infilename);   
 
 case OperationMode::Table:
    // Table-mode work
 s.os+=1;s.fout=stdout;
                          if (strncmp(argv[1],"-tL2",4)==0){s.os+=1;s.tL2=strtod(argv[2],NULL);}
                          if(argc>s.os+2)if (strncmp(argv[s.os+1],"-n",2)==0){s.nofat=(int)strtod(argv[s.os+2],NULL);
                                                                 s.os+=2;
                                           // printf("output only for atom %i\n",s.nofat);
                                              }
  getinfilename(s,argc,argv);
  break;
 // second ... other options with graphics !!
 case OperationMode::Graphical:
    // Graphical-mode work
  if(strcmp(argv[1],"-c")==0){s.os=1;}
  if(strcmp(argv[1],"-s")==0){s.os=1;}
  if(strcmp(argv[1],"-o")==0){s.os=1;} 
  if(strcmp(argv[1],"-m")==0){s.os=1;}
  if(strcmp(argv[1],"-j")==0){s.os=1;}
  if(s.os==1)
  {s.gp.show_density=1;s.density=1;

switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c': s.dim=CHARGEDENS_EV_DIM;
printf("#chargedensity is expanded in tesseral harmonics Zlm\n\
#   ro(r) sum_lm (a(l,m) R^2(r) Zlm(Omega)\n\
#   M. Rotter et al. J Phys: Conf Ser. 325 (2011) 012005\n#\n ");
 snprintf(s.gp.title,MAXNOFCHARINLINE,"chargedensity ro(r)");
 s.gp.threshhold=-0.05;
           break;
 case 's': s.dim=SPINDENS_EV_DIM;
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",s.xx,s.yy,s.zz);
  if(s.doijk==3) snprintf(s.gp.title,MAXNOFCHARINLINE,"projection of spindensity Ms(r).(%g,%g,%g)",s.xx,s.yy,s.zz);
  if(s.doijk==1){snprintf(s.gp.title,MAXNOFCHARINLINE,"divergence of spindensity div Ms(r)");s.gp.scale_density_vectors=0;}
  if(s.doijk==0) snprintf(s.gp.title,MAXNOFCHARINLINE,"abs value  of spindensity |Ms(r)|");
if(s.doijk<3){s.dim*=3;}
s.gp.threshhold=0.05;
break;
 case 'o': s.dim=ORBMOMDENS_EV_DIM;
printf("#orbital momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm  aL(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",s.xx,s.yy,s.zz);
  if(s.doijk==3) snprintf(s.gp.title,MAXNOFCHARINLINE,"projection of orbmomdensity Ms(r).(%g,%g,%g)",s.xx,s.yy,s.zz);
  if(s.doijk==1){snprintf(s.gp.title,MAXNOFCHARINLINE,"divergence of orbmomdensity div ML(r)");s.gp.scale_density_vectors=0;}
  if(s.doijk==0) snprintf(s.gp.title,MAXNOFCHARINLINE,"abs value  of orbmomdensity |ML(r)|");
if(s.doijk<3){s.dim*=3;}
s.gp.threshhold=0.05;
break;
 case 'm': s.dim=SPINDENS_EV_DIM+ORBMOMDENS_EV_DIM;
printf("#magnetic momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm (aS(l,m) R^2(r)+ aL(l,m) F(r)) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",s.xx,s.yy,s.zz);
  if(s.doijk==3) snprintf(s.gp.title,MAXNOFCHARINLINE,"projection of momdensity M(r).(%g,%g,%g)",s.xx,s.yy,s.zz);
  if(s.doijk==1){snprintf(s.gp.title,MAXNOFCHARINLINE,"divergence of momdensity div ML(r)");s.gp.scale_density_vectors=0;}
  if(s.doijk==0) snprintf(s.gp.title,MAXNOFCHARINLINE,"abs value  of momdensity |ML(r)|");
if(s.doijk<3){s.dim*=3;}
s.gp.threshhold=0.05;
break;
 case 'j': s.dim=ORBMOMDENS_EV_DIM;
printf("#currdensity is expanded in tesseral harmonics Zlm\n\
#   j(r).(%g,%g,%g)= sum_lm (b(l,m) R^2(r)+ d(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",s.xx,s.yy,s.zz);
  if(s.doijk==3) snprintf(s.gp.title,MAXNOFCHARINLINE,"projection of currdensity j(r).(i=%g,j=%g,k=%g)(milliAmp/A^2)",s.xx,s.yy,s.zz);
  if(s.doijk==1){snprintf(s.gp.title,MAXNOFCHARINLINE,"divergence of currdensity div j(r)");s.gp.scale_density_vectors=0;}
  if(s.doijk==0) snprintf(s.gp.title,MAXNOFCHARINLINE,"abs value  of currdensity |j(r)|(milliAmp/A^2)");
  s.dim*=6;
s.gp.threshhold=0.05;
break;
 default: help();return false;break;
}
  }
  if(argc>s.os+1)if(strcmp(argv[s.os+1],"-div")==0){s.os+=1;s.doijk=1;}
  else if(strcmp(argv[s.os+1],"-p")==0){s.os+=4;
  s.xx=strtod(argv[3],NULL);
  s.yy=strtod(argv[4],NULL);
  s.zz=strtod(argv[5],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(s.xx*s.xx+s.yy*s.yy+s.zz*s.zz);
  s.xx/=rr;s.yy/=rr;s.zz/=rr;
  s.doijk=3;
                                     }

if(argc>s.os+1)if(strcmp(argv[1+s.os],"-S")==0){s.os+=1;s.arrow=1;s.gp.spins_colour=3; s.gp.spins_scale_moment=1;//arrowdim=SPIN_EV_DIM;
                              snprintf(s.gp.title+strlen(s.gp.title),MAXNOFCHARINLINE-strlen(s.gp.title)," arrows correspond to the spins");}
else if(strcmp(argv[1+s.os],"-L")==0){s.os+=1;s.arrow=2;s.gp.spins_colour=2; s.gp.spins_scale_moment=1;//arrowdim=ORBMOM_EV_DIM;
                                   snprintf(s.gp.title+strlen(s.gp.title),MAXNOFCHARINLINE-strlen(s.gp.title)," arrows correspond to the orbital angular momenta");}
else if(strncmp(argv[1+s.os],"-M",2)==0){s.os+=1;s.arrow=3;s.gp.spins_colour=1; s.gp.spins_scale_moment=1;//arrowdim=MAGMOM_EV_DIM;
                                   snprintf(s.gp.title+strlen(s.gp.title),MAXNOFCHARINLINE-strlen(s.gp.title)," arrows correspond to the magnetic moments");
                                   if(strcmp(argv[s.os],"-Mi")==0){s.arrow=4;}
                                   }
else if(strncmp(argv[1+s.os],"-pel",4)==0){s.os+=1;s.arrow=5;s.gp.spins_colour=4; s.gp.spins_scale_moment=1;//arrowdim=MAGMOM_EV_DIM;
                                   snprintf(s.gp.title+strlen(s.gp.title),MAXNOFCHARINLINE-strlen(s.gp.title)," arrows correspond to the electric dipole moments");
                                   }

if(argc>s.os+1)if(strcmp(argv[1+s.os],"-P")==0){s.os+=1;s.phonon=1;}
if(argc>s.os+1)if(strcmp(argv[1+s.os],"-eps")==0){s.os+=1;s.eps=1;}
if(argc>s.os+1)if(strcmp(argv[1+s.os],"-fst")==0){s.os+=1;s.fst=1;}
getinfilename(s,argc,argv);
break;
}
    // Initialize the graphic parameters here, using s.gp.
s.gp.show_abc_unitcell=1.0;
s.gp.show_primitive_crystal_unitcell=1.0;
s.gp.show_magnetic_unitcell=1.0;
s.gp.show_atoms=1.0;
s.gp.scale_view_1=1.0;
s.gp.scale_view_2=1.0;
s.gp.scale_view_3=1.0;
s.gp.spins_scale_moment=0;
s.gp.show_density=0;
snprintf(s.gp.title,MAXNOFCHARINLINE,"output of program spins");

    // After options have set s.prefix:
// --------------------- load crystal structure information from mcphas.j for further processing  -----------------
//   (needed by check_for_best () to convert Habc into Hijk , then for spins.out crystallographic info ...)
// par inputpars("mcphas.j",prefix);



    s.inputpars.reset(new par("mcphas.j", s.prefix));

//char *out[NOF_USERDEF_MCPHAS_COLS+1];
//for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col){out[col]=new char[20];}
    const int columnCount = NOF_USERDEF_MCPHAS_COLS + 1;
    s.outStorage.reserve(columnCount);
    s.out.resize(columnCount, nullptr);

    for (int col = 1; col < columnCount; ++col) {
        s.outStorage.push_back(std::unique_ptr<char[]>(new char[20]));
        s.out.back(); 
        s.out[col] = s.outStorage.back().get();
    }

    // Move the headerinput(...) call here, updating its arguments to use
    // s.fin, s.gp, s.cs, s.out.data(), and s.N as appropriate.
    //
    // For example, adapt:
    //   headerinput(fin, ..., gp, cs, out, N);
    // to:
    //   headerinput(s.fin, ..., s.gp, s.cs, s.out.data(), s.N);
if (s.mode != OperationMode::Table&&strcmp(argv[1],"-fc")!=0){
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
if(!s.longtable){fclose (s.fin);}


// ----------------------------output configuration ----------------------------------------------------------------

  printf("#! %s - configuration\n",s.outstr);
  fprintf(stdout,"# Components of the Demagnetisation Tensor in SI Units\n");
    fprintf(stdout,"# refering to ijk coordinate system\n");
    fprintf(stdout,"# defined by  j||b, k||(a x b) and i normal to k and j\n");

    fprintf(stdout,"#!Nii=%g Nij=%g Nik=%g\n",s.N(1,1),s.N(1,2),s.N(1,3));
    fprintf(stdout,"#!Njj=%g Njk=%g\n",s.N(2,2),s.N(2,3));
    fprintf(stdout,"#!Nkk=%g\n",s.N(3,3));
if(s.nofcomp>s.savmf->nofcomponents){s.nofcomp=s.savmf->nofcomponents;}
if(s.nofcomp<s.savmf->nofcomponents){printf("#! printing only nofcomponents=%i components\n",s.nofcomp);}
  if (s.mode != OperationMode::Table){ // do not execute this for -t* options
  if(strcmp(argv[1],"-fc")==0){ // for option -fc limit nofat (atoms to be printed)
if(s.nofat>s.savmf->nofatoms){s.nofat=s.savmf->nofatoms;}
if(s.nofat<s.savmf->nofatoms){printf("#! printing only nofatoms=%i atoms\n",s.nofat);}
if(strcmp(argv[s.os]+strlen(argv[s.os])-3,".mf")==0)
{s.savmf->print_commented(stdout,"Hex",s.minl,s.maxl,s.maxn,s.limit,1,s.nofat);}
else{s.savmf->print_commented(stdout,"I",s.minl,s.maxl,s.maxn,s.limit,1,s.nofat);}
                       exit(0);} // exit for -fc option
  else {s.savmf->print(stdout,s.nofcomp);}}

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
                                                 fabs((s.cs.y[ii]-s.cs.y[i])*s.cs.abc[2])>0.001){fprintf(stderr,"Error spins.c: atoms %i (%g %g %g) and %i (%g %g %g) too close\n",i,s.cs.x[i],s.cs.y[i],s.cs.z[i],ii,s.cs.x[ii],s.cs.y[ii],s.cs.z[ii]);exit(EXIT_FAILURE);}
                                              if(s.cs.z[ii]>s.cs.z[i]){s.cs.z[ii]=s.cs.z[i];}else{s.cs.z[i]=s.cs.z[ii];}
                                            }
                                       }
             }
}
s.cs4.abc=s.cs.abc;s.cs4.r=s.cs.r;s.cs4.nofatoms=s.cs.nofatoms;s.cs4.nofcomponents=s.cs.nofcomponents;
s.savmf->calc_prim_mag_unitcell(s.p,s.cs.abc,s.cs.r);
  
  if (strcmp(argv[1],"-f")==0)
 { s.inputpars->savelattice(s.fout);
   fprintf (s.fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(s.fout,"#! strain tensor: eps1=%4.4g (epsii) eps2=%4.4g (epsjj) eps3=%4.4g (epskk) eps4=%4.4g (2epsjk) eps5=%4.4g (2epsik) eps6=%4.4g (2epsij)\n",
    myround(s.savmf->epsilon(1)),myround(s.savmf->epsilon(2)),myround(s.savmf->epsilon(3)),myround(s.savmf->epsilon(4)),myround(s.savmf->epsilon(5)),myround(s.savmf->epsilon(6)));

  fprintf(s.fout,"#! %s \n",s.outstr);
   if(s.T==0){fprintf(s.fout,"# program spins: temperature not found in %s - setting T=1 K\n",argv[2]);s.T=1;}
  fprintf(s.fout,"#!%s: nr1=%i nr2=%i nr3=%i \n",s.outstr,s.savmf->na(),s.savmf->nb(),s.savmf->nc());
  fprintf(s.fout,"#1           2     3     4     5       6       7       8    9    10\n");
  fprintf(s.fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ia> <Ib> <Ic> [created by program spins]\n");
  for (i=1;i<=s.savmf->na();++i){for(j=1;j<=s.savmf->nb();++j){for(k=1;k<=s.savmf->nc();++k)
  {for(ii=1;ii<=s.inputpars->cs.nofatoms;++ii)
   {// output the positions
    dd3=s.savmf->pos_dabc(i,j,k,ii, s.cs);
   //returns position dd3 as components	    with respect to lattice a b c
   
   // dd0 is position as components with respect to primitive crystallographic unit cell
    dd0=s.savmf->pos_dr123(i,j,k,ii, s.cs);
    fprintf(s.fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            (*s.inputpars->jjj[ii]).sipffilename,dd3(1),dd3(2),dd3(3),dd0(1),dd0(2),dd0(3));
    //output the "magnetic" moment if possible ... actually it outputs Ia Ib Ic 
    for(nt=1;nt<=3;++nt){fprintf(s.fout," %4.4f",myround(1e-5,s.savmf->m(i,j,k)(s.inputpars->cs.nofcomponents*(ii-1)+nt)));}      
    fprintf(s.fout,"\n");
   }
  }}}
  

  fclose(s.fout);
    s.gp.spins_scale_moment=1;s.gp.phonon_scale_static_displacements=1;
    s.gp.read();
    s.fin = fopen_errchk ("./results/spins_prim.jvx", "w"); // here draw some graphics for the spinconfiguration
                       // with option -f. Yet this is very limited, because we do not know what is the magnetic
                       // moment, i.e. a guess is made using the components 1 2 3 of the s.savmf 
     s.gp.showprim=1;
     s.gp.scale_density_vectors=0;s.gp.show_density=0; // do not show charge densities
             Vector hkl1(1,3);hkl1=0;
             Vector gjmbHxc1(1,3);gjmbHxc1=0;
     spincf magmom(s.savmf->na(),s.savmf->nb(),s.savmf->nc(),s.savmf->nofatoms,3); // the magnetic moment guess
     int i,j,k,l;for (i=1;i<=s.savmf->na();++i){for (j=1;j<=s.savmf->nb();++j){for (k=1;k<=s.savmf->nc();++k){for(l=1;l<=magmom.nofatoms;++l)  
         for(int momdim=1;momdim<=3&&momdim<=s.savmf->nofcomponents;++momdim)
               {magmom.moment(i,j,k,l)(momdim)=s.savmf->moment(i,j,k,l)(momdim);}
        }}}

     s.savmf->jvx_cd(s.fin,s.outstr,s.cs,
                  s.gp,0.0,(*s.savmf)*0.0,(*s.savmf)*0.0,
                  hkl1,s.T,gjmbHxc1,s.Hext,s.cs,magmom,magmom * 0.0,magmom* 0.0,magmom * 0.0);
    fclose (s.fin);

  exit(0);
}
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

/**********************************************************************/
// hauptprogramm
/**********************************************************************/
int main (int argc, char **argv)
{ double x,y;

  

LoadedInput loaded;

    if (!parse_options_and_load_structure(argc, argv, loaded)) {
        if (loaded.fin) fclose(loaded.fin);
        if (loaded.fout && loaded.fout != stdout) fclose(loaded.fout);
        return EXIT_FAILURE;
    }
// Example conversions:
    // fin       -> loaded.fin
    // gp        -> loaded.gp
    // cs        -> loaded.cs
    // inputpars -> *loaded.inputpars
    // savmf     -> *loaded.savmf
    // out       -> loaded.out.data()


  spincf spinconf(loaded.savmf->na(),loaded.savmf->nb(),loaded.savmf->nc(),loaded.cs4.nofatoms,3);
  spincf sc_phonon(loaded.savmf->na(),loaded.savmf->nb(),loaded.savmf->nc(),loaded.cs4.nofatoms,3);

fprintf (loaded.fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
   fprintf(loaded.fout,"#! strain tensor: eps1=%4.4g (epsii) eps2=%4.4g (epsjj) eps3=%4.4g (epskk) eps4=%4.4g (2epsjk) eps5=%4.4g (2epsik) eps6=%4.4g (2epsij)\n",
       myround(loaded.savmf->epsilon(1)),myround(loaded.savmf->epsilon(2)),myround(loaded.savmf->epsilon(3)),myround(loaded.savmf->epsilon(4)),myround(loaded.savmf->epsilon(5)),myround(loaded.savmf->epsilon(6)));

if (loaded.mode != OperationMode::Table){
// the following is for the printout of spins.out ...........................
section4header(loaded.fout);
fprintf(loaded.fout,"#! %s : ",loaded.outstr);
fprintf(loaded.fout," nr1=%i nr2=%i nr3=%i nat=%i  atoms in primitive magnetic unit cell\n",loaded.savmf->na(),loaded.savmf->nb(),loaded.savmf->nc(),loaded.cs4.nofatoms*loaded.savmf->na()*loaded.savmf->nb()*loaded.savmf->nc());
//MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
//                                        to Sa Sb Sc La Lb Lc
fprintf(loaded.fout,"#1           2     3     4     5       6       7       8    9    10                  11   12   13   14   15   16\n");
fprintf(loaded.fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)\n");
fprintf(loaded.fout,"#          corresponding exchange fields [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)\n");
// .............................................................................                                
	}
else  //now table options
 { fprintf(loaded.fout,"#! nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",loaded.savmf->na(),loaded.savmf->nb(),loaded.savmf->nc(),loaded.cs4.nofatoms*loaded.savmf->na()*loaded.savmf->nb()*loaded.savmf->nc());
 //fprintf(loaded.fout,"# 1 2  3  4  5           6     7     8     9       10      11      12   13   14                  15   16   17   18   19   20\n");
 //fprintf(loaded.fout,"# external parameters {sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] ");
 fprintf(loaded.fout,"#");
 for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col)fprintf(loaded.fout,"%s ",loaded.out.data()[col]);
 fprintf(loaded.fout," n nr1 nr2 nr3 {sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] ");
 if (strcmp(argv[1],"-tMSL")==0)  fprintf(loaded.fout,"<Ma> <Mb> <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> (hbar)]\n");
 else if (strcmp(argv[1],"-tI")==0) fprintf(loaded.fout,"<I1> <I2> <I3> ... <Inofcomponents>\n");
 else if (strcmp(argv[1],"-tHex")==0) fprintf(loaded.fout,"<Hex1> <Hex2> <Hex3> ... <Inofcomponents> (meV)\n");
 else if (strcmp(argv[1],"-tL2")==0){ fprintf(loaded.fout,"bondlength_to_neighbour [A]  elongation^2 [A^2]\n");
                                nlimits_calc(loaded.nmin,loaded.nmax, loaded.tL2,loaded.p);
                               }
 else {fprintf(stderr, "Error spins: option %s not known\n",argv[1]);exit(EXIT_FAILURE);}
 }       
//  1. from the meanfieldconfiguration (*loaded.savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: Icalc has to be used to calculate all the <Olm>.
 
Vector h(1,loaded.inputpars->cs.nofcomponents);
Vector hh(1, loaded.savmf->nofcomponents * loaded.savmf->nofatoms);
Vector I(1,loaded.inputpars->cs.nofcomponents);
spincf densitycf(
    loaded.savmf->na(),
    loaded.savmf->nb(),
    loaded.savmf->nc(),
    loaded.savmf->nofatoms,
    loaded.dim);
h=0;for(int ii=1;ii<=loaded.inputpars->cs.nofatoms;++ii)
{(*loaded.inputpars->jjj[ii]).Icalc_parameter_storage_init(h,loaded.Hext,loaded.T);} // initialize Icalc module parameter storage


bool moretoprint=true;
while(moretoprint){
 for (int i=1;i<=loaded.savmf->na();++i)for(int j=1;j<=loaded.savmf->nb();++j)for(int k=1;k<=loaded.savmf->nc();++k)
 {
    hh=loaded.savmf->m(i,j,k);
  densitycf.m(i,j,k)=0;int i4=1,ii4=1;
  for(int ii=1;ii<=loaded.inputpars->cs.nofatoms;++ii)
 {  
    Vector magmom(1,3),mom(1,3);
    Vector Lmom(1,3);
    Vector Smom(1,3);
    Vector moments(1,loaded.dim);
      Vector momS(1,SPINDENS_EV_DIM);
  Vector momL(1,ORBMOMDENS_EV_DIM);
  Vector momentsx(1,SPINDENS_EV_DIM);
  Vector momentsy(1,SPINDENS_EV_DIM);
  Vector momentsz(1,SPINDENS_EV_DIM);
  Vector momentlx(1,ORBMOMDENS_EV_DIM);
  Vector momently(1,ORBMOMDENS_EV_DIM);
  Vector momentlz(1,ORBMOMDENS_EV_DIM);
    h=0;
   for(int nt=1;nt<=loaded.inputpars->cs.nofcomponents;++nt){h(nt)=hh(nt+loaded.inputpars->cs.nofcomponents*(ii-1));}


switch(loaded.arrow)
{case 1: (*loaded.inputpars->jjj[ii]).Scalc(mom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage);break;
 case 2: (*loaded.inputpars->jjj[ii]).Lcalc(mom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage);break;
 case 3: (*loaded.inputpars->jjj[ii]).mcalc(mom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage);break;
 case 5: (*loaded.inputpars->jjj[ii]).pelcalc(mom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage);break;
 default: break;
}

switch(loaded.arrow)
{case 1:
 case 2:
 case 3:
 case 5:
         for(int nt=1;nt<=3;++nt)
		        {spinconf.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };break;
 case 4: int dim4;
         dim4=3;if((*loaded.inputpars->jjj[ii]).module_clust==true)
                           dim4=(*(*loaded.inputpars->jjj[ii]).clusterpars).cs.nofatoms*3;
         Vector momi(1,dim4);
         (*loaded.inputpars->jjj[ii]).micalc(momi,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage);
         // now put the components of momi to the spinconf                                               
         for(int nt=1;nt<=dim4;++nt){spinconf.m(i,j,k)(nt+i4-1)=momi(nt);};i4+=dim4;
         break;
}
   Vector dd0(1,3),dd3(1,3);
    
// output atoms and moments in primitive unit cell to loaded.fout  ------------------------------------
if(loaded.arrow==4&&(*loaded.inputpars->jjj[ii]).module_clust==true){
    for(int nt=1;nt<=(*(*loaded.inputpars->jjj[ii]).clusterpars).cs.nofatoms;++nt)
     {dd3=spinconf.pos_dabc(i,j,k,ii4, loaded.cs4);
      dd0=spinconf.pos_dr123(i,j,k,ii4, loaded.cs4);
      if (strncmp(argv[1],"-t",2)==0&&(loaded.nofat==MAXNOFATOMS||loaded.nofat==ii)){
              //fprintf(loaded.fout," %s ", loaded.outstr);
       double value;
      for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col)
       {extract(loaded.outstr,loaded.out.data()[col],value); fprintf(loaded.fout,"%g ",myround(value));}
      fprintf(loaded.fout," %i %i %i %i ",loaded.savmf->n(),loaded.savmf->na(),loaded.savmf->nb(),loaded.savmf->nc());
                                  
                                     }
      fprintf(loaded.fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            loaded.cs4.sipffilenames[ii4],dd3(1),dd3(2),dd3(3),dd0(1),dd0(2),dd0(3));
      if (loaded.mode != OperationMode::Table||strcmp(argv[1],"-tMSL")==0)
      if(loaded.nofat==MAXNOFATOMS||loaded.nofat==ii){
            fprintf(loaded.fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(1+3*(ii4-1))));
            fprintf(loaded.fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(2+3*(ii4-1))));
            fprintf(loaded.fout," %4.4f\n",myround(1e-5,spinconf.m(i,j,k)(3+3*(ii4-1))));
                                                                     } 
       // do not print out L S or exchange fields for cluster module
      ++ii4;
     }
    }
else{   ++ii4;
    // output the positions
    dd3=loaded.savmf->pos(i,j,k,ii, loaded.cs);
      // if module allows to calculate position shift of an atom - use this for output 
    if(true==(*loaded.inputpars->jjj[ii]).pcalc(mom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage))
     {dd3+=mom; 
     //  fprintf(stderr,"# Attention: atom %i shifted from equilibrium position by (%g %g %g) A \n",ii,mom(1),mom(2),mom(3));
     }
     // if module cannot calculate position shift of an atom - check if the atom 
     // is actually a charge cloud on top of another atom and if yes take position shift
     // from this other atom
      else
     {for(int iii=1;iii<=loaded.inputpars->cs.nofatoms;++iii)
      {if(ii!=iii)
       {Vector ddc(1,3); ddc=loaded.savmf->pos(i,j,k,iii, loaded.cs);
        if(abs(ddc-dd3)<0.01){
Vector hhh(1,loaded.inputpars->cs.nofcomponents); hhh=0;
for(int nt=1;nt<=loaded.inputpars->cs.nofcomponents;++nt){hhh(nt)=hh(nt+loaded.inputpars->cs.nofcomponents*(iii-1));}

 if(true==(*loaded.inputpars->jjj[iii]).pcalc(mom,loaded.T,hhh,loaded.Hext,(*loaded.inputpars->jjj[iii]).Icalc_parstorage))
    {dd3+=mom;
 //printf(stderr,"# Attention: atom %i shifted from equilibrium position by (%g %g %g) A taken from atom %i\n",ii,mom(1),mom(2),mom(3),iii);
 }
                               }
       }
      } 
     }

Vector dd(1,3);
char  tstr [MAXNOFCHARINLINE];
    dd0=loaded.p.Inverse()*dd3;
dd0(1)*=loaded.savmf->na();dd0(2)*=loaded.savmf->nb();dd0(3)*=loaded.savmf->nc();
    Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,loaded.cs.abc);
    dd=abc_in_ijk.Inverse()*dd3;
 if(loaded.nofat==MAXNOFATOMS||loaded.nofat==ii)
 { snprintf(tstr,MAXNOFCHARINLINE,"{%s} %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f ",
            loaded.cs.sipffilenames[ii],dd(1),dd(2),dd(3),dd0(1),dd0(2),dd0(3));
if (strcmp(argv[1],"-tL2")!=0) {     
    if (strncmp(argv[1],"-t",2)==0){
       //fprintf(loaded.fout,"%s ",loaded.outstr);
       double value;
      for (int col=1;col<=NOF_USERDEF_MCPHAS_COLS;++col)
       {extract(loaded.outstr,loaded.out.data()[col],value); fprintf(loaded.fout,"%g ",myround(value));}
      fprintf(loaded.fout," %i %i %i %i ",loaded.savmf->n(),loaded.savmf->na(),loaded.savmf->nb(),loaded.savmf->nc());
                                   }
      fprintf(loaded.fout,"%s",tstr);

    }
    if (loaded.mode != OperationMode::Table||strcmp(argv[1],"-tMSL")==0)
    {
     //ouput the magnetic moment if possible
     if((*loaded.inputpars->jjj[ii]).mcalc(magmom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage))
     {     for(int nt=1;nt<=3;++nt){fprintf(loaded.fout," %4.4f",myround(1e-5,magmom(nt)));}
      // and output the orbital and spin momentum if possible 
      if((*loaded.inputpars->jjj[ii]).Lcalc(Lmom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage)&&
        (*loaded.inputpars->jjj[ii]).Scalc(Smom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage))
       //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
       //                                        to Sa Sb Sc La Lb Lc
       //{        for(int nt=1;nt<=3;++nt){fprintf(loaded.fout," %4.4f %4.4f",myround(1e-5,Smom(nt)),myround(1e-5,Lmom(nt)));}
       {        for(int nt=1;nt<=3;++nt){fprintf(loaded.fout," %4.4f",myround(1e-5,Smom(nt)));}
               for(int nt=1;nt<=3;++nt){fprintf(loaded.fout," %4.4f",myround(1e-5,Lmom(nt)));}
       }
      }}
    }
   if (loaded.mode != OperationMode::Table||strcmp(argv[1],"-tHex")==0)
   {
    // finally output a line with the exchange fields 
    if (loaded.mode != OperationMode::Table)fprintf(loaded.fout,"\n                 corresponding exchange fields [meV]-->          ");
                      for(int nt=1;nt<=loaded.savmf->nofcomponents;++nt)  // printout exchangefields
                        {fprintf(loaded.fout," %4.4f",myround(1e-5,h(nt)));}
                        
   }
   if (strcmp(argv[1],"-tI")==0)
   {double lnZ;double U;
     (*loaded.inputpars->jjj[ii]).Icalc(I,loaded.T,h,loaded.Hext,lnZ,U,(*loaded.inputpars->jjj[ii]).Icalc_parstorage);
                         for(int nt=1;nt<=loaded.savmf->nofcomponents;++nt)  // printout I operator expectation values
                        {fprintf(loaded.fout," %4.4f",myround(1e-5,I(nt)));}
                         
   }
   if (strcmp(argv[1],"-tL2")==0)
   {// here calculate up to tl2 Angstroem bondlengths and elongation^2 and output to loaded.fout
    L2_calc(loaded.tL2,i,j,k,ii,loaded.fout,loaded.outstr,tstr,loaded.nmin,loaded.nmax,*loaded.savmf,loaded.p,loaded.T,h,loaded.Hext,*loaded.inputpars);
   }
   else   fprintf(loaded.fout,"\n");
 }
// -----------------------------------------------------------------------------------------------
   

if(loaded.density){
switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c':  (*loaded.inputpars->jjj[ii]).chargedensity_coeff (moments,loaded.T, h, loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage); break;
 case 's':  if(loaded.xx!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).spindensity_coeff (momentsx,1,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.yy!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).spindensity_coeff (momentsy,2,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.zz!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).spindensity_coeff (momentsz,3,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.doijk==3){ moments=loaded.xx*momentsx+loaded.yy*momentsy+loaded.zz*momentsz;}
            else{for(int i=1;i<=SPINDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'o':  if(loaded.xx!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentsx,1,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.yy!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentsy,2,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.zz!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentsz,3,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.doijk==3){ moments=loaded.xx*momentsx+loaded.yy*momentsy+loaded.zz*momentsz;}
            else{for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+ORBMOMDENS_EV_DIM)=momentsy(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'm':  if(loaded.xx!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).spindensity_coeff (momentsx,1,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.yy!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).spindensity_coeff (momentsy,2,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.zz!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).spindensity_coeff (momentsz,3,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            momS=loaded.xx*momentsx+loaded.yy*momentsy+loaded.zz*momentsz;
            if(loaded.xx!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentlx,1,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.yy!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momently,2,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            if(loaded.zz!=0||loaded.doijk<3)(*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentlz,3,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            momL=loaded.xx*momentlx+loaded.yy*momently+loaded.zz*momentlz;
            for(int i=1;i<=SPINDENS_EV_DIM;++i){
            if(loaded.doijk==3){moments(i)=momS(i);moments(i+SPINDENS_EV_DIM)=momL(i);}
            else{moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);
                 moments(i+3*SPINDENS_EV_DIM)=momentlx(i);moments(i+4*SPINDENS_EV_DIM)=momently(i);moments(i+5*SPINDENS_EV_DIM)=momentlz(i);
                }
                                               }
            break;
 case 'j':  (*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentlx,1,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            (*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momently,2,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            (*loaded.inputpars->jjj[ii]).orbmomdensity_coeff (momentlz,3,loaded.T, h,loaded.Hext, (*loaded.inputpars->jjj[ii]).Icalc_parstorage);
            for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){
             moments(i)=momentlx(i);moments(i+ORBMOMDENS_EV_DIM)=momently(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentlz(i);
             }
            break;
 default: help();return false;
}
                   for(int nt=1;nt<=loaded.dim;++nt)
		        {densitycf.m(i,j,k)(nt+loaded.dim*(ii-1))=moments(nt);
                    }
} // loaded.gp.show_density


if(loaded.phonon==1)  // if module allows to calculate position  - use this for graphics ...
{if(true==(*loaded.inputpars->jjj[ii]).pcalc(mom,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage))
{for(int nt=1;nt<=3;++nt)
		        {sc_phonon.m(i,j,k)(nt+3*(ii-1))=mom(nt); // here we set moment to be output as arrow
                    };
}
 
}// if cluster or atoms
} // for ii
} // for ijk
if(loaded.longtable){if(check_for_best(loaded.fin,loaded.aa,*loaded.savmf,x,y,loaded.T,loaded.Hext,loaded.outstr,loaded.out.data(),loaded.inputpars->cs.abc,1))
              {fclose(loaded.fin);exit(0);}
                           }
else {
      moretoprint=false;}
} // while moretoprint
if (loaded.mode==OperationMode::Table){exit(0);}
  fclose (loaded.fout);
   
// create plot of spinconfiguration -----------------------------------------------------------
printf("# ************************************************************************\n");
printf("#%s\n",loaded.gp.title);
printf("# ************************************************************************\n");
              if(loaded.arrow==0)loaded.gp.spins_scale_moment=0;
              if(loaded.density==0)loaded.gp.show_density=0;
if(loaded.eps==1){
    loaded.fin = fopen_errchk ("./results/spins.eps", "w");
     spinconf.eps(loaded.fin,loaded.outstr);
    fclose (loaded.fin);

// here the 3d file should be created
    loaded.fin = fopen_errchk ("./results/spinsab.eps", "w");


     spinconf.eps3d(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,1,spinconf);
    fclose (loaded.fin);
    loaded.fin = fopen_errchk ("./results/spinsac.eps", "w");
     spinconf.eps3d(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,2,spinconf);
    fclose (loaded.fin);
    loaded.fin = fopen_errchk ("./results/spinsbc.eps", "w");
     spinconf.eps3d(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,3,spinconf);
    fclose (loaded.fin);
    loaded.fin = fopen_errchk ("./results/spins3dab.eps", "w");
     spinconf.eps3d(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,4,spinconf);
    fclose (loaded.fin);
    loaded.fin = fopen_errchk ("./results/spins3dac.eps", "w");
     spinconf.eps3d(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,5,spinconf);
    fclose (loaded.fin);
    loaded.fin = fopen_errchk ("./results/spins3dbc.eps", "w");
     spinconf.eps3d(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,6,spinconf);
    fclose (loaded.fin);
          }
if(loaded.fst==1){
    loaded.fin = fopen_errchk ("./results/spins.fst", "w");
     spinconf.fst(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,spinconf);
    fclose (loaded.fin);

    
   loaded.fin = fopen_errchk ("./results/spins_prim.fst", "w");
     spinconf.fstprim(loaded.fin,loaded.outstr,loaded.cs4.abc,loaded.cs4.r,loaded.cs4.x,loaded.cs4.y,loaded.cs4.z,spinconf);
    fclose (loaded.fin);
          }
             Vector hkl(1,3);hkl=0;double E;E=0;
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
for(int ii=1;ii<=loaded.inputpars->cs.nofatoms;++ii)
{Vector pos(1,3);   
  for(int nt=1;nt<=loaded.inputpars->cs.nofcomponents;++nt){h(nt)=hh(nt+loaded.inputpars->cs.nofcomponents*(ii-1));}
  if(true==(*loaded.inputpars->jjj[ii]).pcalc(pos,loaded.T,h,loaded.Hext,(*loaded.inputpars->jjj[ii]).Icalc_parstorage))
 {double charge;charge=(*loaded.inputpars->jjj[ii]).charge;if(charge==0)charge=0.01;
   snprintf(loaded.cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g",
           charge,loaded.gp.scale_pointcharges*0.529177*signum(charge)*pow((double)fabs(charge),0.3333));
  // if we can determine from sipffilename the element --> put r g b information
    
  unsigned int r,g,b; char element []  ="E\0";double radius=0.3333;
   element[0]=toupper((*loaded.inputpars->jjj[ii]).sipffilename[0]);
   if(isalpha((*loaded.inputpars->jjj[ii]).sipffilename[1]))
   element[1]=tolower((*loaded.inputpars->jjj[ii]).sipffilename[1]);
  for(int i=0;i<NOFELEMENTS_RGB;++i)
   if(0==strncmp(element,elstr[i],2)){
r=(unsigned int)(unsigned char)elstr[i][2];
g=(unsigned int)(unsigned char)elstr[i][3];
b=(unsigned int)(unsigned char)elstr[i][4];
radius=strtod(elstr[i]+5,NULL);
   snprintf(loaded.cs.sipffilenames[ii],MAXNOFCHARINLINE,"pointcharge %g |e| radius=%g %s r=%i g=%i b=%i",
           charge,loaded.gp.scale_pointcharges*radius*signum(charge)*pow((double)fabs(charge),0.3333),
           element,r,g,b);
//    printf("%s scale=%g radius=%g charge=%g sig(charge)=%g pow=%g  estr=%s\n",loaded.cs.sipffilenames[ii],
//loaded.gp.scale_pointcharges,radius,charge,signum(charge),pow((double)fabs(charge),0.3333),elstr[i]+5);

}

// printf("#! atom %i %s displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,loaded.cs.sipffilenames[ii],ii,pos(1),ii,pos(2),ii,pos(3));
  printf("#! atom %i  displacement u%ix=%g u%iy=%g u%iz=%g A\n",ii,ii,pos(1),ii,pos(2),ii,pos(3));
 }
}

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   loaded.fin = fopen_errchk ("./results/spins.jvx", "w");
    loaded.gp.showprim=0;loaded.gp.spins_wave_amplitude=0;
     densitycf.jvx_cd(loaded.fin,loaded.outstr,loaded.cs,loaded.gp,0.0,densityev_real,densityev_imag,hkl,loaded.T,gjmbHxc,loaded.Hext,loaded.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon);
    fclose (loaded.fin);

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   loaded.fin = fopen_errchk ("./results/spins_prim.jvx", "w");
     loaded.gp.showprim=1;
     densitycf.jvx_cd(loaded.fin,loaded.outstr,loaded.cs,loaded.gp,0.0,densityev_real,densityev_imag,hkl,loaded.T,gjmbHxc,loaded.Hext,loaded.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon);
    fclose (loaded.fin);

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
// try movie - a spinwave picture  ... including phonons and spindensity changes ...
//***************************************************************************************************************
//***************************************************************************************************************
// check if there is a second 
for(int i=loaded.os+1;i<argc;++i){
if(strcmp(argv[i],"-prefix")==0){strcpy(loaded.prefix,argv[i+1]); // read alternative loaded.prefix
                                   fprintf(stdout,"# prefix for qee qep etc input filenames: %s\n",loaded.prefix);
 				   loaded.os+=2;}
                       }
if(argc-1==NOF_USERDEF_MCPHAS_COLS+4+loaded.os)loaded.os+=NOF_USERDEF_MCPHAS_COLS-4;
if(argc-1==6+loaded.os)loaded.os-=2;  // 
//argc-1==8+loaded.os  ... ok 
if (argc-1==8+loaded.os){
              // double E;
              char outhklstr[MAXNOFCHARINLINE];
             loaded.gp.spins_wave_amplitude=1.0;loaded.gp.spins_show_ellipses=1.0;loaded.gp.spins_show_oscillation=1.0;
             loaded.gp.phonon_wave_amplitude=1.0;loaded.gp.phonon_scale_static_displacements=1.0;
//----------------------------------------------------------------------------------------------------------
 
           if(loaded.arrow>0){
             strcpy(loaded.infilename,"./results/");strcpy(loaded.infilename+10,loaded.prefix);
             switch(loaded.arrow)
             {case 1: strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qes");loaded.fin = fopen(loaded.infilename, "rb");
                      if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qes", "rb");
                      break;
              case 2: strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qeo");loaded.fin = fopen(loaded.infilename, "rb");
                      if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qeo", "rb");
                      break;
              case 3: strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qem");loaded.fin = fopen(loaded.infilename, "rb");
                      if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qem", "rb");
                      break;
              case 4: fprintf(stderr,"mcdisp: output of individual moment oscillation in eigenvector file mcdisp.qemi not yet implemented - thus exiting program spins\n");
                      strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qemi");loaded.fin = fopen(loaded.infilename, "rb");
                      if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qemi", "rb");
                      break;   
              case 5: strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qpe");loaded.fin = fopen(loaded.infilename, "rb");
                      if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qpe", "rb");
                      break;                   
             }
             check_for_best_excitation_and_close(loaded.fin,loaded.gp,"qev",spinconf,spinconfev_real,spinconfev_imag,
                          x,y,loaded.T,loaded.Hext,loaded.outstr,outhklstr,"moment oscillation",hkl,E,argc,argv,loaded.os,
                          loaded.inputpars->cs.abc,loaded.inputpars->cs.nofatoms,loaded.doijk, loaded.xx,loaded.yy,loaded.zz,loaded.dim);

             }//arrow
     
//----------------------------------------------------------------------------------------------------------
            if(loaded.phonon>0){
             strcpy(loaded.infilename,"./results/");strcpy(loaded.infilename+10,loaded.prefix);
             strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qep");loaded.fin = fopen(loaded.infilename, "rb");
                      if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qep", "rb");
             check_for_best_excitation_and_close(loaded.fin,loaded.gp,"qep",spinconf,spinconfpev_real,spinconfpev_imag,
                          x,y,loaded.T,loaded.Hext,loaded.outstr,outhklstr,"phonon loaded.oscillation",hkl,E,argc,argv,loaded.os,
                          loaded.inputpars->cs.abc,loaded.inputpars->cs.nofatoms,loaded.doijk, loaded.xx,loaded.yy,loaded.zz,loaded.dim);

             }//phonon
//----------------------------------------------------------------------------------------------------------
            if(loaded.density){
             strcpy(loaded.infilename,"./results/");strcpy(loaded.infilename+10,loaded.prefix);
             switch(argv[1][1]) // dimension deloaded.finition from jjjpar.hpp
                {case 'c': strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qee");loaded.fin = fopen(loaded.infilename, "rb");
                           if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qee", "rb");
                           break;
                 case 's': strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qsd");loaded.fin = fopen(loaded.infilename, "rb");
                           if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qsd", "rb");
                           break;
                 case 'o': strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qod");loaded.fin = fopen(loaded.infilename, "rb");
                           if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qod", "rb");
                           break;
                 case 'm': fprintf(stderr,"Error spins: magnetic moment density oscillation not yet implemented\n");exit(1);break;
                             // would have to look into qsd and qod files !!
                 case 'j': strcpy(loaded.infilename+10+strlen(loaded.prefix),"mcdisp.qod");loaded.fin = fopen(loaded.infilename, "rb");
                           if(loaded.fin==NULL)loaded.fin = fopen_errchk ("./results/mcdisp.qod", "rb");
                           break;
                }
             check_for_best_excitation_and_close(loaded.fin,loaded.gp,"qee/qsd/qod",densitycf,densityev_real,densityev_imag,
          x,y,loaded.T,loaded.Hext,loaded.outstr,outhklstr,"density oscillation",hkl,E,argc,argv,loaded.os,
          loaded.inputpars->cs.abc,loaded.inputpars->cs.nofatoms,loaded.doijk, loaded.xx,loaded.yy,loaded.zz,loaded.dim);

             }//loaded.gp.show_density
//----------------------------------------------------------------------------------------------------------           
              loaded.gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
                        // in case he wants to overwrite some default settings
              if(loaded.arrow==0)loaded.gp.spins_scale_moment=0;
              if(loaded.density==0)loaded.gp.show_density=0;
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
             // complex <double> im(0,1);
              for(int i=0;i<loaded.gp.nof_frames;++i)
              {phase=2*3.1415*i/(loaded.gp.nof_frames-1);
               printf("\n********************************************\n");
               printf(" calculating movie sequence %i(%i)\n",i+1,loaded.gp.nof_frames);
               printf("********************************************\n");
               char filename[MAXNOFCHARINLINE];
               snprintf(filename,MAXNOFCHARINLINE,"./results/spins.%02d.jvx",i+1);
               loaded.fin = fopen_errchk (filename, "w");loaded.gp.showprim=0;
                     densitycf.jvx_cd(loaded.fin,outhklstr,loaded.cs,loaded.gp,
                                  phase,densityev_real,densityev_imag,hkl,loaded.T,hh,loaded.Hext,loaded.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon,spinconfpev_real,spinconfpev_imag);
               fclose (loaded.fin);
               snprintf(filename,MAXNOFCHARINLINE,"./results/spins_prim.%02d.jvx",i+1);
               loaded.fin = fopen_errchk (filename, "w");loaded.gp.showprim=1;
                     densitycf.jvx_cd(loaded.fin,outhklstr,loaded.cs,loaded.gp,
                                  phase,densityev_real,densityev_imag,hkl,loaded.T,hh,loaded.Hext,loaded.cs4,spinconf,spinconfev_real,spinconfev_imag,sc_phonon,spinconfpev_real,spinconfpev_imag);
               fclose (loaded.fin);
              }
          printf("# %s\n",outhklstr);
          }
fprintf(stderr,"# %s\n",loaded.outstr);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spins running at \n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/spins.jvx\n");
fprintf(stderr,"# * java JvxAnimator \"results/spins.*.jvx\" \n");
fprintf(stderr,"# ************************************************************************\n");
for (int i = 1; i <= loaded.cs.nofatoms; ++i) {
    delete[] loaded.cs.sipffilenames[i];
    loaded.cs.sipffilenames[i] = nullptr;
}
  return 0;
}


