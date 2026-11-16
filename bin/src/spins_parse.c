// parsing command line functions for program spins
#define MAXDEG 100    // maximal degeneracy
#define MAXNOFSPINCF 10000000 // maximal number of spinconfigurations to be loaded

enum class OperationMode {
    Graphical,
    File,
    Table,
    Eigenvector
};

enum class TableMode {
    None,
    TI,
    THex,
    TL2,
    TMSL
};

enum class FileMode {
    None,
    F,
    FC
};


enum class SymmetryMode {
    None,
    Symop,
    Symoa
};

struct EigenvectorOptions {
    SymmetryMode symmetryMode;
    int symmetryIndex;             // signed value, e.g. 3 or -3
    std::string cif;  // CIF filename
    std::string qFile;
    std::array<double, 3> hkl;
    bool hasEnergy;
    double energy;

    EigenvectorOptions()
        : symmetryMode(SymmetryMode::None),
          symmetryIndex(0),
          hkl{{0.0, 0.0, 0.0}},
          hasEnergy(false),
          energy(0.0)
    {}
};



struct LoadedInput {
    OperationMode mode = OperationMode::Graphical;
    TableMode tableMode;
    FileMode fileMode;
    EigenvectorOptions eigenvector;

    LoadedInput()
        : mode(OperationMode::Graphical),
          tableMode(TableMode::None),
          fileMode(FileMode::None) {}

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
" or results/graphic_parameters.gp by show_abc_unitcell,\n"
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



static bool parseDoubleArg(const char* text, double& value)
{
    if (text == NULL || *text == '\0')
        return false;

    errno = 0;
    char* end = NULL;
    const double parsed = std::strtod(text, &end);

    if (errno != 0 || end == text || *end != '\0')
        return false;

    value = parsed;
    return true;
}

static bool parseIntArg(const char* text, int& value)
{
    if (text == NULL || *text == '\0')
        return false;

    errno = 0;
    char* end = NULL;
    const long parsed = std::strtol(text, &end, 10);

    if (errno != 0 || end == text || *end != '\0' ||
        parsed < static_cast<long>(-2147483647 - 1) ||
        parsed > static_cast<long>(2147483647))
        return false;

    value = static_cast<int>(parsed);
    return true;
}

bool parseEigenvectorOptions(
    int argc,
    char** argv,
    int firstArg,
    LoadedInput& s,
    std::string& error)
{
    EigenvectorOptions options;
    int i = firstArg;

    // Optional: -symop <signed-index> <crystal-structure.cif>
    //        or -symoa <signed-index> <crystal-structure.cif>
    if (i < argc &&
        (std::string(argv[i]) == "-symop" ||
         std::string(argv[i]) == "-symoa"))
    {
        options.symmetryMode =
            (std::string(argv[i]) == "-symop")
                ? SymmetryMode::Symop
                : SymmetryMode::Symoa;
        ++i;

        if (i >= argc || !parseIntArg(argv[i], options.symmetryIndex)) {
            error = "Expected a signed symmetry index after -symop/-symoa";
            return false;
        }
        ++i;

        if (i >= argc) {
            error = "Expected a crystal-structure CIF filename after the symmetry index";
            return false;
        }
        options.cif = argv[i++];
    }

    if (i >= argc || std::string(argv[i]) != "-q") {
        error = "Expected -q followed by a q-file and h k l values";
        return false;
    }
    ++i;

    if (i >= argc) {
        error = "Expected a q-file after -q";
        return false;
    }
    options.qFile = argv[i++];
    s.os=i-5; // needed in spins_q to call check_for_best_excitation_and_close
    for (int axis = 0; axis < 3; ++axis) {
        if (i >= argc || !parseDoubleArg(argv[i], options.hkl[axis])) {
            error = "Expected three numeric h k l values after the q-file";
            return false;
        }
        ++i;
    }

    if (i < argc) {
        if (!parseDoubleArg(argv[i], options.energy)) {
            error = "Optional energy must be numeric";
            return false;
        }
        options.hasEnergy = true;
        ++i;
    }

    if (i != argc) {
        error = "Unexpected extra argument in eigenvector mode";
        return false;
    }

    s.mode = OperationMode::Eigenvector;
    s.eigenvector = options;
    return true;
}

bool detectOperationMode(int argc, char** argv, LoadedInput& s)
{ std::string error;
const int firstEigenvectorArg = s.os + 1;
    s.mode = OperationMode::Graphical; // default

   

    // Handle other modes using the program's actual command syntax.
    

if (argc > 1 && strcmp(argv[1],"-fc") == 0) {
    s.mode = OperationMode::File;
    s.fileMode = FileMode::FC;
    s.os = 1;
}
else if (argc > 1 && strcmp(argv[1],"-f") == 0) {
    s.mode = OperationMode::File;
    s.fileMode = FileMode::F;
    s.os = 1;
}
else if (argc > 1 && strcmp(argv[1],"-tI") == 0) {
    s.mode = OperationMode::Table;
    s.tableMode = TableMode::TI;
    s.os = 1;
}
else if (argc > 1 && strcmp(argv[1],"-tHex") == 0) {
    s.mode = OperationMode::Table;
    s.tableMode = TableMode::THex;
    s.os = 1;
}
else if (argc > 1 && strcmp(argv[1],"-tL2") == 0) {
    s.mode = OperationMode::Table;
    s.tableMode = TableMode::TL2;
    s.os = 1;
}
else if (argc > 1 && strcmp(argv[1],"-tMSL") == 0) {
    s.mode = OperationMode::Table;
    s.tableMode = TableMode::TMSL;
    s.os = 1;
}
else 
{ for( int i = 1 ; i < argc ; ++i)
  { if(strcmp(argv[i],"-q") == 0)
    {s.mode = OperationMode::Eigenvector;
 if (!parseEigenvectorOptions(argc, argv, firstEigenvectorArg, s, error)) {
    std::cerr << "Error: " << error << '\n';
    return false;
      } else  {return true;}
    }
   }
}
s.mode = OperationMode::Graphical;
    return true;
}

// to be done for table and graphical mode - get the correct .mf filename
void getinfilename(LoadedInput& s, int argc, char** argv)
{if(argc>s.os+1)if(strcmp(argv[1+s.os],"-prefix")==0){strcpy(s.prefix,argv[2+s.os]); // read s.prefix
                                   fprintf(stdout,"# prefix for input filenames: %s\n",s.prefix);
 				   s.os+=2;}
 

 strcpy(s.infilename,"./results/");strcpy(s.infilename+10,s.prefix);
 strcpy(s.infilename+10+strlen(s.prefix),"mcphas.mf");s.fin = fopen(s.infilename, "rb");
 if(s.fin==NULL){strcpy(s.infilename+10,"mcphas.mf");s.fin = fopen_errchk(s.infilename, "rb");}
 printf("# reading from file %s\n",s.infilename);
 
 }
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

bool parse_options(
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
if(!detectOperationMode(argc,argv,s)){return false;}

s.prefix[0]='\0';
// options without graphics just screendump <I> or exchange field configuration at given HT
switch (s.mode) {
case OperationMode::File:
s.os+=1;
if (s.fileMode == FileMode::FC)
  {s.minl=(int)strtod(argv[s.os],NULL);
   s.maxl=(int)strtod(argv[s.os+1],NULL);
   s.maxn=(int)strtod(argv[s.os+2],NULL);
   s.limit=strtod(argv[s.os+3],NULL);
   s.os+=4;}
   if (strncmp(argv[s.os],"-n",2)==0){
    if(s.fileMode == FileMode::FC)s.nofat=(int)strtod(argv[s.os+1],NULL);
    else s.nofcomp=(int)strtod(argv[s.os+1],NULL);
                              s.os+=2;}
   s.fin = fopen_errchk (argv[s.os], "rb");
   strcpy(s.infilename,argv[s.os]);
   printf("#* program spins ... reading from file %s\n",s.infilename);   
break;
 case OperationMode::Table:
    // Table-mode work
    s.fout=stdout;
                          if (strncmp(argv[1],"-tL2",4)==0){s.os+=1;s.tL2=strtod(argv[2],NULL);}
                          if(argc>s.os+2)if (strncmp(argv[s.os+1],"-n",2)==0){s.nofat=(int)strtod(argv[s.os+2],NULL);
                                                                 s.os+=2;
                                           // printf("output only for atom %i\n",s.nofat);
                                              }
  getinfilename(s,argc,argv); //  get the correct .mf filename
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
 default: help();return 1;break;
}
  }
  

if (argc > s.os + 1) {
    if (std::strcmp(argv[s.os + 1], "-div") == 0) {
        s.os += 1;
        s.doijk = 1;
    } else if (std::strcmp(argv[s.os + 1], "-p") == 0) {
        // -p plus three direction values must be present.
        if (argc <= s.os + 4) {
            // Handle missing direction values.
        } else {
            s.xx = std::strtod(argv[s.os + 2], NULL);
            s.yy = std::strtod(argv[s.os + 3], NULL);
            s.zz = std::strtod(argv[s.os + 4], NULL);

            const double rr = std::sqrt(s.xx * s.xx +
                                        s.yy * s.yy +
                                        s.zz * s.zz);
            if (rr == 0.0) {
                // Handle a zero-length direction vector.
            } else {
                s.xx /= rr;
                s.yy /= rr;
                s.zz /= rr;
                s.os += 4; // Point to the last consumed argument.
                s.doijk = 3;
            }
        }
    }
}

/*
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
*/
if (argc > s.os + 1) {
    const char* option = argv[s.os + 1];

    if (std::strcmp(option, "-S") == 0) {
        ++s.os;
        s.arrow = 1;
        s.gp.spins_colour = 3;
        s.gp.spins_scale_moment = 1;
        std::snprintf(s.gp.title + std::strlen(s.gp.title),
                      MAXNOFCHARINLINE - std::strlen(s.gp.title),
                      " arrows correspond to the spins");

    } else if (std::strcmp(option, "-L") == 0) {
        ++s.os;
        s.arrow = 2;
        s.gp.spins_colour = 2;
        s.gp.spins_scale_moment = 1;
        std::snprintf(s.gp.title + std::strlen(s.gp.title),
                      MAXNOFCHARINLINE - std::strlen(s.gp.title),
                      " arrows correspond to the orbital angular momenta");

    } else if (std::strncmp(option, "-M", 2) == 0) {
        ++s.os;
        s.arrow = (std::strcmp(option, "-Mi") == 0) ? 4 : 3;
        s.gp.spins_colour = 1;
        s.gp.spins_scale_moment = 1;
        std::snprintf(s.gp.title + std::strlen(s.gp.title),
                      MAXNOFCHARINLINE - std::strlen(s.gp.title),
                      " arrows correspond to the magnetic moments");

    } else if (std::strncmp(option, "-pel", 4) == 0) {
        ++s.os;
        s.arrow = 5;
        s.gp.spins_colour = 4;
        s.gp.spins_scale_moment = 1;
        std::snprintf(s.gp.title + std::strlen(s.gp.title),
                      MAXNOFCHARINLINE - std::strlen(s.gp.title),
                      " arrows correspond to the electric dipole moments");
    }
}


if(argc>s.os+1)if(strcmp(argv[1+s.os],"-P")==0){s.os+=1;s.phonon=1;}
if(argc>s.os+1)if(strcmp(argv[1+s.os],"-eps")==0){s.os+=1;s.eps=1;}
if(argc>s.os+1)if(strcmp(argv[1+s.os],"-fst")==0){s.os+=1;s.fst=1;}
getinfilename(s,argc,argv);
break;
default: break;
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
return true;
}