// methods for class inimcdis 
#include "inimcdiff.hpp"
#include <martin.h>
#include <mfcf.hpp>
#include <mcdiff.h>
#include "../../version"

#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__FreeBSD__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <windows.h>
#endif

#define NOF_OUT_VARIABLES  nofoutputcolumns+1


int cc[]= {2,31,32,33,34,35,36,27,28,29,1,0,30}; // field to store code for assigning type of data to columns of output,
                                           // set default values here (see list below for different types)
                                           // using the out0 ... out12 variable in mcdiff.in these codes can be modified

#define COLHEADDIM 47

// different output data for columns 1 to and 12
const char * ch []= {  "LF ",          //  0
                            "|NSF|[b] ",    //  1
                            "Re(NSF)[b] ",    //  2
                            "Im(NSF)[b] ",    //  3
                            "|MSF| ",    //  4
                            "|MSF.P| ",    //  5
                            "Re(MSF.P) ",    //  6
                            "Im(MSF.P) ",    //  7
                            "|MSFdip| ",    //  8
                            "|MSFdip.P| ",    //  9
                            "Re(MSFdip.P) ",    //  10
                            "Im(MSFdip.P) ",    //  11
                            "angl(Q,P)[°] ",    //  12
                            "i(MSFxMSF*).P ",    //  13
                            "I+ ",    //  14
                            "I- ",    //  15
                            "I+/I- ",    // 16
                            "i(MSFxMSF*)dip.P ",    //  17
                            "Idip+ ",    //  18
                            "Idip- ",    //  19
                            "Idip+/Idip- ",    //  20
                            "2*|MSF.P|/sin^2(angl(Q,P) ",    //  21
                            "2*|MSFdip.P|/sin^2(angl(Q,P) ",    //  22
                            "2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-) ",    //  23
                            "2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-) ",    //  24
                            "2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-) ",    //  25
                            "2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-) ",    //  26
                            "Inuc(2t)   ",    //  27
                            "Imag(2t)   ",    //  28
                            "Itot(2t)   ",    //  29
                            "Imag_dip(2t) ",    //  30
                            "h      ",    //  31
                            "k      ",    //  32
                            "l      ",    //  33
                            "d[A]       ",    //  34             
                            "|Q|[1/A]   ",    //  35
                            "2theta     ",    //  36
                            "Qi[1/A]    ",    //  37    euclidean components of scattering vector 
                            "Qj[1/A]    ",    //  38    with j||b, k||(a x b) and i normal to k and j
                            "Qk[1/A]    ",    //  39
                            "T[K]       ",    //  40
                            "Ha[T]      ",    //  41
                            "Hb[T]      ",    //  42
                            "Hc[T]      ",    //  43
                            "hprim      ",    //  44
                            "kprim      ",    //  45
                            "lprim      ",    //  46
			    "Itotdip(2t) ",    //  47
                                           };


 // *************************************************************************
 // ************************ inimcdiff *************************************
 // *************************************************************************
 // class of initial parameters for program mcdiff
 // *************************************************************************

// print user defined column headers
void inimcdiff::print_usrdefcolhead(FILE *fout)
{fprintf(fout,"#");
 for(int i=1;i<= nofoutputcolumns;++i){fprintf(fout,"%i%*s",i,(int)strlen(colhead[colcod[i]])-1,"");}
fprintf(fout,"\n#");
 for(int i=1;i<= nofoutputcolumns;++i)fprintf(fout,"%s",colhead[colcod[i]]);
}

// print user defined columns
void inimcdiff::print_usrdefcols(FILE *fout,float ** out,int i)
{
for(int j=1;j<= nofoutputcolumns;++j){switch(colcod[j]){case 31: case 32: case 33: case 44: case 45:case 46:{fprintf(fout,"%6.3f ",myround(out[j][i]));break;}default: fprintf(fout,"%5.4E ",out[j][i]);}}
      
}






// save parameters (which were read from mcdiff.in)
void inimcdiff::save()
{save(savfilename);
}

// save parameters (which were read from mcdiff.in)
void inimcdiff::save(const char * filename)
{  FILE * fout;int i;
  fout=fopen(filename,"w");if (fout==NULL) {fprintf(stderr,"ERROR - file %s cannot be opened \n",filename);exit(EXIT_FAILURE);} 
  fprintf(fout,"# this file is the input file read by program %s ",MCDIFFVERSION);
  time_t curtime;
  struct tm * loctime;
  curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
  fprintf(fout,"#<!--mcdiff.mcdiff.in>\n");
fprintf(fout,"#***************************************************************\n");
fprintf(fout,"#      mcdiff is a program for the calculation of elastic\n");
fprintf(fout,"#   neutron diffraction and resonant magnetic Xray scattering \n");
fprintf(fout,"#  reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
fprintf(fout,"#*************************************************************** \n");
fprintf(fout,"# this input file contains 4 sections corresponding to different\n");
fprintf(fout,"# groups of parameters\n");
fprintf(fout,"#\n");
fprintf(fout,"# - all lines have to start with a # sign with the  exception of \n");
fprintf(fout,"#   the lines containing atomic positional parameters\n");
fprintf(fout,"# - the other parameters have to be defined in the corresponding \n");
fprintf(fout,"#   section by statements such as parameter=value\n");
fprintf(fout,"# - the sequence of the parameters within a section is arbitrary\n");
fprintf(fout,"# \n");
fprintf(fout,"#\n");

fprintf(fout,"# %%SECTION 1%%  OVERALL PARAMETERS\n");
fprintf(fout,"#\n");
fprintf(fout,"#! lambda   = %g  wavelength (A)\n",lambda);
fprintf(fout,"#\n");
fprintf(fout,"#! thetamax = %g   maximum bragg angle (deg)\n",thetamax);
fprintf(fout,"#\n");
fprintf(fout,"#! ovalltemp= %g  overall temperature factor (A^2) \n",ovalltemp);
fprintf(fout,"#           ...I ~ EXP(-2 * ovalltemp * sintheta^2 / lambda^2) \n");
fprintf(fout,"#                  relation to other notations:\n");
fprintf(fout,"#                  ovalltemp = Biso = 8 pi^2 Uiso^2\n");
fprintf(fout,"#\n");
fprintf(fout,"#! lorentz=%i  type of lorentzfactor to be used\n",lorenz);
fprintf(fout,"#            0.....no lorentzfactor \n");
fprintf(fout,"#            1.....neutron powder flat sample\n");
fprintf(fout,"#            2.....neutron powder cylindrical sample\n");
fprintf(fout,"#            3.....neutron single crystal\n");
fprintf(fout,"#            4.....neutron TOF powder cyl. sample - d-pattern log scaled\n");
fprintf(fout,"#            5.....neutron TOF powder cyl. sample - d-pattern normal scaled\n#\n");
fprintf(fout,"# out*  controls the type of output in  mcdiff.out \n");
 fprintf(fout,"#!");for(int i=0;i<= nofoutputcolumns;++i){fprintf(fout,"out%i=%i ",i,colcod[i]);
if(i==0)fprintf(fout," 0..short header, 1...standard header, 2... long header in mcdiff.out\n" \
                     "# (a negative value triggers calculation of magnetic xray scattering)\n" \
                     "#\n#!  nofoutputcolumns=%i  number of columns in output file mcdiff.out\n" \
                     "# choose out* to set type of desired output in column 1 to %i\n" \
                     "# (default is h k l d Q 2theta Inuc Imag Itot)\n#\n#!", nofoutputcolumns, nofoutputcolumns);
  }fprintf(fout,"\n#\n");
  fprintf(fout,"#     ... in out*=n the numbers n have the following meaning:\n");
  for(i=0;i<=COLHEADDIM;++i){
  fprintf(fout,"#            %i....%s",i,colhead[i]);
if(i==38)fprintf (fout," Qi Qj Qk are euclidean components of scattering vector");
if(i==39)fprintf (fout," with j||b, k||(a x b) and i normal to k and j");
  fprintf(fout,"#\n");
                   }
  fprintf(fout,"#\n");
fprintf(fout,"#           In the above the intensities I+ and I- are the intensities in a polarised neutron\n");
fprintf(fout,"#           experiment with incident polarisation up (+) or down (-):\n");
fprintf(fout,"#            I+-=LF exp(-OTF Q^2/8pi^2) \n");
fprintf(fout,"#                    [ |NSF/NB|^2 + 3.65/4pi (|MSF|^2-+i(MSF x MSF*).P)/NB^2 \n");
fprintf(fout,"#                        +-  sqrt(3.65/4pi)/NB^2 (NSF (MSF*.P) + NSF* (MSF.P)]\n"
             "#           LF  ..... Lorentzfactor\n"
             "#           MSF ..... magnetic structure factor\n"
             "#           NSF ..... nuclear structure factor\n");
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"#             For some of the above options we need the\n");
fprintf(fout,"#! Pa=%8.4f   Components of Polarisation Vector in terms of lattice vectors P=(Pa * a + Pb * b + Pc *c)\n",P(1));
fprintf(fout,"#! Pb=%8.4f   Note: the length of P, i.e. |P| indicates the degree of beam polarisation (|P|<=1)\n",P(2));
fprintf(fout,"#! Pc=%8.4f\n",P(3));
fprintf(fout,"#\n");
fprintf(fout,"#\n");

 fprintf(fout,"# %%SECTION 2%% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL\n");
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"#! natcryst=%i      number of nonmagnetic atoms in primitive crystalographic unit cell\n",nat);
fprintf(fout,"#\n");
fprintf(fout,"# it follows a list of natcryst lines with nonmagnetic atoms\n");
fprintf(fout,"# ... notes: - if an occupancy other than 1.0 is needed, just reduce \n");
fprintf(fout,"#              the scattering length linear accordingly\n");
fprintf(fout,"#            - Debye Waller Factor notation: sqr(Intensity) ~ structure factor ~ \n");
fprintf(fout,"#              ~sum_n ()n exp(-2 DWFn sin^2(theta) / lambda^2)=EXP (-Wn),  \n");
fprintf(fout,"#              relation to other notations: 2*DWF = B = 8 pi^2 <u^2>, units DWF (A^2)\n");
fprintf(fout,"#\n");
fprintf(fout,"#! use_dadbdc=%i\n",use_dadbdc);
fprintf(fout,"#            - 0 means: da db and dc are not used by the program (unless you enter a line #! use_dadbdc=1),\n");
fprintf(fout,"#               dr1,dr2 and dr3 refer to the primitive lattice given below\n");
fprintf(fout,"# Real Imag[scattering length(10^-12cm)]   da(a)    db(b)    dc(c)    dr1(r1)  dr2(r2)  dr3(r3)  DWF(A^2)\n");
 if (nat!=0){ 
               for(i=1;i<=nat;++i) {
   fprintf(fout,"  %8.5f  %8.5f                       %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",sl1r[i],sl1i[i],da[i],db[i],dc[i],x1[i],y1[i],z1[i],dwf1[i]);
                                   }
             }
 fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"# %%SECTION 3%% DESCRIPTION OF THE LATTICE\n");
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"# Note: what follows here may directly be taken from the output of program spins \n");
fprintf(fout,"#       (file spins.out) or charges (file charges.out)\n");
fprintf(fout,"# -----------------------------------------------------------------------------\n");
fprintf(fout,"#\n");
fprintf(fout,"# lattice constants (A) and angles \n");
fprintf(fout,"#! a=%g b=%g c=%g alpha=  %g beta=  %g gamma=  %g\n",a,b,c,alpha,beta,gamma);
fprintf(fout,"#\n");
fprintf(fout,"# primitive lattice vectors \n");
fprintf(fout,"#! r1a= %7f r2a= %7f r3a= %7f\n",r1(1),r2(1),r3(1));
fprintf(fout,"#! r1b= %7f r2b= %7f r3b= %7f   primitive lattice vectors (a)(b)(c)\n",r1(2),r2(2),r3(2));
fprintf(fout,"#! r1c= %7f r2c= %7f r3c= %7f\n",r1(3),r2(3),r3(3));
fprintf(fout,"#\n");
fprintf (fout, "#      - coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
fprintf(fout,"#! strain tensor (optional): eps1=%4.4g=epsii eps2=%4.4g=epsjj eps3=%4.4g=epskk eps4=%4.4g=2epsjk eps5=%4.4g=2epsik eps6=%4.4g=2epsij\n",
        eps(1,1),eps(2,2),eps(3,3),eps(2,3),eps(1,3),eps(1,2));
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"# %%SECTION 4%% DESCRIPTION OF MAGNETIC UNIT CELL AND LIST OF MAGNETIC ATOMS\n");
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"# here follows the description of the magnetic unit cell with respect\n");
fprintf(fout,"# to the primitive crystallographic unit cell:\n");
fprintf(fout,"# 'nr1', 'nr2', 'nr3' ...the crystallographic unit cell has to be taken \n");
fprintf(fout,"#                        nr1 nr2 and nr3 times along r1 r2 and r3,\n");
fprintf(fout,"#                        respectively to get magnetic unit cell\n");
fprintf(fout,"# 'nat' denotes the number of magnetic atoms in magnetic unit cell\n");
fprintf(fout,"#\n");
fprintf(fout,"# Temperature,  External Magnetic Field: Magnetic Unit Cell\n");
fprintf(fout,"#! T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i \n",T,H(1),H(2),H(3),nr1,nr2,nr3,natmagnetic);
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"# It follows a list of nat lines with to describe the magnetic moment configuration\n");
fprintf(fout,"# Notes:\n");
fprintf(fout,"# 'atom-filename' means the single ion property filename of this magnetic atom:\n");
fprintf(fout,"#                 -it must contain the Formfactor Coefficients (e.g. see international tables)\n");
fprintf(fout,"#                                      Lande factor\n");
fprintf(fout,"#                                      Neutron Scattering Length (10^-12 cm) \n");
fprintf(fout,"#                 -it may contain a    Debey Waller Factor\n");
fprintf(fout,"# 'da' 'db' and 'dc' are not used by the program (unless you enter a line #! use_dadbdc=1)\n");
fprintf(fout,"# 'dr1','dr2' and 'dr3' refer to the primitive lattice given below\n");
fprintf(fout,"# 'Ma','Mb','Mc' denote the magnetic moment components in Bohr magnetons\n");
fprintf(fout,"#                in case of non orthogonal lattices instead of Ma Mb Mc the components Mx My Mz\n");
fprintf(fout,"#                have to be given, which refer to an right handed orthogonal coordinate system \n");
fprintf(fout,"#                defined by y||b, z||(a x b) and x normal to y and z\n");
//MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
//                                        to Sa Sb Sc La Lb Lc
//fprintf(fout,"#  <Sa>  <La> <Sb> <Lb >  <Sc> <Lc>  (optional) denote the spin and orbital angular momentum components \n");
fprintf(fout,"#  <Sa> <Sb> <Sc>  <La> <Lb > <Lc>  (optional) denote the spin and orbital angular momentum components \n");
fprintf(fout,"# 'mf1' 'mf2' 'mf3' (optional line, used to go beyond dipole approx for formfactor)\n");
fprintf(fout,"#                                     denote the corresponding exchange fields in meV\n");
fprintf(fout,"#\n");
//MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
//fprintf(fout,"#{atom-file} da[a]  db[b]    dc[c]     dr1[r1]  dr2[r2]  dr3[r3]   <Ma>     <Mb>     <Mc> [mb] [optional <Sa> <La> <Sb> <Lb> <Sc> <Lc> ]\n");
fprintf(fout,"#{atom-file} da[a]  db[b]    dc[c]     dr1[r1]  dr2[r2]  dr3[r3]   <Ma>     <Mb>     <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> ]\n");
fprintf(fout,"#{corresponding exchange fields [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)}\n");

for(i=1;i<=natmagnetic;++i){
fprintf(fout,"{%s} %8.5f %8.5f %8.5f  ",(*jjjpars[i]).sipffilename,(*jjjpars[i]).xyz(1)*nr1*r1s(1)+(*jjjpars[i]).xyz(2)*nr2*r2s(1)+(*jjjpars[i]).xyz(3)*nr3*r3s(1),
                                                                   (*jjjpars[i]).xyz(1)*nr1*r1s(2)+(*jjjpars[i]).xyz(2)*nr2*r2s(2)+(*jjjpars[i]).xyz(3)*nr3*r3s(2),
                                                                   (*jjjpars[i]).xyz(1)*nr1*r1s(3)+(*jjjpars[i]).xyz(2)*nr2*r2s(3)+(*jjjpars[i]).xyz(3)*nr3*r3s(3));
fprintf(fout,"%8.5f %8.5f %8.5f  ",(*jjjpars[i]).xyz(1)*nr1,(*jjjpars[i]).xyz(2)*nr2,(*jjjpars[i]).xyz(3)*nr3); // positions
fprintf(fout," %+8.5f %+8.5f %+8.5f ",(*jjjpars[i]).mom(1),(*jjjpars[i]).mom(2),(*jjjpars[i]).mom(3)); // magmoments
fprintf(fout," %+8.5f",(*jjjpars[i]).mom(4));
fprintf(fout," %+8.5f",(*jjjpars[i]).mom(6));
fprintf(fout," %+8.5f",(*jjjpars[i]).mom(8));
fprintf(fout," %+8.5f",(*jjjpars[i]).mom(5));
fprintf(fout," %+8.5f",(*jjjpars[i]).mom(7));
fprintf(fout," %+8.5f",(*jjjpars[i]).mom(9));
fprintf(fout,"\n");
if((*jjjpars[i]).FF_type<0) // i.e. go beyond 
          {                  
fprintf(fout,"                    corresponding exchange fields [meV]-->");
for(int k=1;k<=(*jjjpars[i]).MF.Hi();++k){fprintf(fout," %+8.5f",(*jjjpars[i]).MF(k));}
fprintf(fout,"\n");
          }
 	
                             }
  fclose (fout);
}
// *************************************************************************
//constructor ... load initial parameters from file
inimcdiff::inimcdiff (const char * file,char * pref,int verb)
{ errno=1;long int pos=0;int j,k;nofoutputcolumns=12;
  verbose=verb;
  char instr[MAXNOFCHARINLINE],somestring[MAXNOFCHARINLINE],sipffilename[MAXNOFCHARINLINE],infile[MAXNOFCHARINLINE]; 
  //,hklline[MAXNOFCHARINLINE];
 // Hext=Vector(1,3);
  P=Vector(1,3);
  FILE *fin; //,*finhkl;float N,M,h0,k0,l0,h1,k1,l1,hN,kN,lN,hM,kM,lM;
  prefix= new char [strlen(pref)+1]; strcpy(prefix,pref); // set prefix
  //********************************  
  savfilename= new char [strlen(file)+strlen(prefix)+11];
  outfilename= new char [strlen(prefix)+strlen("./results/mcdiff.out")+1];
  snprintf(outfilename,MAXNOFCHARINLINE,"./results/%smcdiff.out",prefix);
  // check if directory results exists and can be written to ...
  fin = fopen_errchk (outfilename, "a");fclose(fin);

  errno = 0;
  // **************** initialize parameters to default values ********************************************
  T=0;
 // ******************************** reading parameters  from mcdiff.in ****************************************************
  //int i=0; //,hklblock=0,QxQyQzblock=0,j;
  snprintf(infile,MAXNOFCHARINLINE,"%s%s",prefix,file);// use for the moment savfilename to store input filename
  fin = fopen(infile, "rb"); 
  if (fin==NULL) { // try standard file without prefix
                    printf("input file %s not found - try reading file %s\n",infile,file);
                    snprintf(infile,MAXNOFCHARINLINE,"%s",file);
                    fin = fopen(infile, "rb"); 
  if (fin==NULL)  {fprintf(stderr,"# Error -  mcdiff input file %s not found\n",infile);exit(1);}
                 }
  else { printf("reading file %s\n",infile);}

  // ...  set savfilename correctly to save mcdiff.in
  snprintf(savfilename,MAXNOFCHARINLINE,"results/_%s%s",prefix,file);

// input section 1 *******************************************************

  instr[0]='#';pos=ftell(fin);
// look if nofoutputcolumns is there ....
 while (instr[strspn(instr," \t")]=='#'&&strstr (instr, "%SECTION 2%")==NULL) // pointer to 'ltrimstring' 
  { if (pos==-1) {fprintf(stderr,"Error mcdiff: wrong mcdiff.in file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin); 
   extract(instr,"nofoutputcolumns",nofoutputcolumns);  
  }
colcod=new int[NOF_OUT_VARIABLES+1];for(int i=0;i<NOF_OUT_VARIABLES;++i){if(i<=12){colcod[i]=cc[i];}else{colcod[i]=0;}}
colhead=new char *[COLHEADDIM+1];for(int i=0;i<=COLHEADDIM;++i){colhead[i]=new char [strlen(ch[i])+1];strcpy(colhead[i],ch[i]);}

fseek(fin,pos,SEEK_SET); // go back to beginning of file and load other pars of section 1
instr[0]='#';instr[1]='\0';
while (instr[strspn(instr," \t")]=='#'&&strstr (instr, "%SECTION 2%")==NULL) // pointer to 'ltrimstring' 
  { pos=ftell(fin); 
    if (pos==-1) {fprintf(stderr,"Error mcdiff: wrong mcdiff.in file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin); 
   extract(instr,"nofatoms",nofatoms);  
   extract(instr,"lambda", lambda);
   extract(instr, "thetamax", thetamax);
   extract(instr, "nat", nat);
   extract(instr, "natcryst", nat);
   extract(instr, "ovalltemp", ovalltemp);
   extract(instr, "lorentz", lorenz);
   for(int i=0;i<= nofoutputcolumns;++i) // extract user defined output columns
   {snprintf(somestring,MAXNOFCHARINLINE,"out%i",i);
    extract(instr, somestring,colcod[i]);
   }
   extract(instr, "Pa",P(1));
   extract(instr, "Pb",P(2));
   extract(instr, "Pc",P(3));
   extract(instr, "nofthreads",P(3));
  }

  fseek(fin,pos,SEEK_SET); 
if (lorenz == 0){fprintf(stderr,"Warning mcdiff: read lorentz=0, will calculate no Lorentzfactor.\n");}
if (lambda == 0){fprintf(stderr,"ERROR mcdiff: no wavelength lambda given or line does not start with # in section 1\n");exit(EXIT_FAILURE);}
if (thetamax == 0){fprintf(stderr,"ERROR mcdiff: no thetamax given or line does not start with # in section 1\n");exit(EXIT_FAILURE);}
 // Checks if nofthreads set in mcdiff.in, if not check environment or use system calls
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
 // printf("# nofthreads=%i\n",nofthreads);
 
printf("     section 1 - lambda=%g A thetamax= %g deg\n",lambda, thetamax);
printf("                 ovalltemp=%g A^2 lorentz-type=%i\n",ovalltemp,lorenz);
printf("                 output: column ");
for(int i=2;i<= nofoutputcolumns;++i){printf(" %i=%s",i,colhead[colcod[i]]);
                                  if(!(i%5))printf("\n                                ");}
printf("\n");

// input section 2  *******************************************************

 instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&a==0&&b==0&&c==0) // pointer to 'ltrimstring' 
  { pos=ftell(fin); 
    if (pos==-1)  {fprintf(stderr,"Error mcdiff: wrong mcdiff.in file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin); 
   extract(instr, "nat", nat);
   extract(instr, "natcryst", nat);
   extract(instr, "a", a);
   extract(instr, "b", b);
   extract(instr, "c", c);
   extract(instr, "alpha", alpha);
   extract(instr, "beta", beta);
   extract(instr, "gamma", gamma);
   extract(instr, "use_dadbdc",use_dadbdc);
  }
  fseek(fin,pos,SEEK_SET); 


  printf ("     section 2 - natcryst=%i\n",nat);
  x1=new float[nat+1];y1=new float[nat+1];z1=new float[nat+1];
  da=new float[nat+1];db=new float[nat+1];dc=new float[nat+1];
  sl1r=new float[nat+1];sl1i=new float[nat+1];;dwf1=new float[nat+1];

  float numbers[70];numbers[0]=70;
  float numbers1[70];numbers1[0]=70;
  if (nat!=0){ for(int i=1;i<=nat;++i) { pos=ftell(fin); 
                                     int n=inputline(fin,numbers);
                                     if (n==0) {if(feof(fin)==true){fprintf(stderr,"Error mcdiff: end of input file in section 2\n");exit (EXIT_FAILURE);}
                                                fseek(fin,pos,SEEK_SET); 
                                                fgets(instr,MAXNOFCHARINLINE,fin); 
                                                if(strstr (instr, "%%SECTION 3%%")!=NULL){fprintf (stderr,"ERROR mcdiff: Section 3 started before all nat=%i atoms of crystallographic unit cell were listed !\n",nat);exit (EXIT_FAILURE);}
                                               --i;}
                                     else      {if (n<9) {fprintf (stderr,"ERROR mcdiff: Section 2 - Nonmagnetic Atoms: too few positional parameters for atom %i!\n",i);exit (EXIT_FAILURE);}
                                                sl1r[i]=numbers[1];sl1i[i]=numbers[2]; x1[i] = numbers[6]; y1[i] = numbers[7]; z1[i] = numbers[8];dwf1[i]=numbers[9];
                                                                                       da[i] = numbers[3]; db[i] = numbers[4]; dc[i] = numbers[5];
                                                printf("                 sl=%g%+gi 10^-12cm at %g*r1%+g*r2%+g*r3 DWF=%g A^2\n",sl1r[i],sl1i[i],x1[i],y1[i],z1[i],dwf1[i]);
                                               }
                                    }
              }


// input section 3 *********************************************************
  instr[0]='#';
 eps=Matrix(1,3,1,3);r1=Vector(1,3);r2=Vector(1,3);r3=Vector(1,3);H=Vector(1,3);

 nr1=0;nr2=0; nr3=0;
 while (instr[strspn(instr," \t")]=='#'&&nr1*nr2*nr3==0) 
  { pos=ftell(fin); 
   fgets(instr,MAXNOFCHARINLINE,fin); 
   if(a==0)extract(instr, "a", a);
   if(b==0)extract(instr, "b", b);
   if(c==0)extract(instr, "c", c);
   if(alpha==0)extract(instr, "alpha", alpha);
   if(beta==0)extract(instr, "beta", beta);
   if(gamma==0)extract(instr, "gamma", gamma);
    extract(instr,"eps1",eps(1,1));
    extract(instr,"eps2",eps(2,2));
    extract(instr,"eps3",eps(3,3));
    extract(instr,"eps4",eps(2,3));
    extract(instr,"eps5",eps(1,3));
    extract(instr,"eps6",eps(1,2));

    extract(instr, "r1x", r1(1));
    extract(instr, "r1y", r1(2));
    extract(instr, "r1z", r1(3));
    extract(instr, "r2x", r2(1));
    extract(instr, "r2y", r2(2));
    extract(instr, "r2z", r2(3));
    extract(instr, "r3x", r3(1));
    extract(instr, "r3y", r3(2));
    extract(instr, "r3z", r3(3));
    extract(instr, "r1a", r1(1));
    extract(instr, "r1b", r1(2));
    extract(instr, "r1c", r1(3));
    extract(instr, "r2a", r2(1));
    extract(instr, "r2b", r2(2));
    extract(instr, "r2c", r2(3));
    extract(instr, "r3a", r3(1));
    extract(instr, "r3b", r3(2));
    extract(instr, "r3c", r3(3));
    extract(instr, "nr1", nr1);
    extract(instr, "nr2", nr2);
    extract(instr, "nr3", nr3);
    extract(instr, "nat", natmagnetic);
    extract(instr, "T", T);
    extract(instr, "Ha", H(1));
    extract(instr, "Hb", H(2));
    extract(instr, "Hc", H(3));
  }

  fseek(fin,pos,SEEK_SET); 
if (a == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant a given in section 3 or line does not start with # or nat too small: \n%s\n",instr);exit (EXIT_FAILURE);}
if (b == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant b given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
if (c == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant c given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
printf("     section 3 - a=%g A  b=%g A c=%g A alpha=%g  beta=%g gamma=%g\n",a,b,c,alpha,beta,gamma);
 unitcellstr=new char[MAXNOFCHARINLINE+1];
snprintf(unitcellstr,sizeof(unitcellstr)," a= %g A  b= %g A c= %g A  alpha=%g  beta=%g gamma=%g\n",a,b,c,alpha,beta,gamma);
printf("                 r1= %5.3ga + %5.3gb + %5.3gc\n", r1(1), r1(2), r1(3));
printf("                 r2= %5.3ga + %5.3gb + %5.3gc\n", r2(1), r2(2), r2(3));
printf("                 r3= %5.3ga + %5.3gb + %5.3gc\n", r3(1), r3(2), r3(3));

//printf("                    / %5.3ga \\     / %5.3ga \\     / %5.3ga \\    x||c \n", r1(1), r2(1), r3(1));
//printf("                 r1=| %5.3gb |  r2=| %5.3gb |  r3=| %5.3gb |    y||a\n", r1(2), r2(2), r3(2));
//printf("                    \\ %5.3gc /     \\ %5.3gc /     \\ %5.3gc /    z||b\n", r1(3), r2(3), r3(3));

  double da1,db1,dc1;
  rez1=Vector(1,3);rez2=Vector(1,3);rez3=Vector(1,3);
  rezcalc (r1, r2, r3, rez1, rez2, rez3);
  if (nat!=0){ for(int i=1;i<=nat;++i) {
   // calculate da db dc from dr1 dr2 dr3 and print to results/_mcdiff.in
       da1= x1[i]*r1(1)+y1[i]*r2(1)+z1[i]*r3(1)-da[i];
       db1= x1[i]*r1(2)+y1[i]*r2(2)+z1[i]*r3(2)-db[i];
       dc1= x1[i]*r1(3)+y1[i]*r2(3)+z1[i]*r3(3)-dc[i];
       double dd=sqrt(da1*da1+db1*db1+dc1*dc1);
                                 da1=x1[i]- (da[i]*rez1(1)+db[i]*rez1(2)+dc[i]*rez1(3))/2/PI;
                                 db1=y1[i]- (da[i]*rez2(1)+db[i]*rez2(2)+dc[i]*rez2(3))/2/PI;
                                 dc1=z1[i]- (da[i]*rez3(1)+db[i]*rez3(2)+dc[i]*rez3(3))/2/PI;
       dd+=sqrt(da1*da1+db1*db1+dc1*dc1);
       if(dd>SMALLPOSITIONDEVIATION){fprintf (stderr,"Warning: atomic positions da db dc and dr1 dr2 dr3 inconsistent !\n");
                    fprintf (stderr,"         use_dadbdc=%i\n",use_dadbdc);
                    if(use_dadbdc==0){ fprintf (stderr,"using dr1 dr2 dr3 and recalculating da db dc...\n");}
                    else {fprintf (stderr,"using da db dc and recalculating dr1 dr2 dr3...\n");}
                   i=nat;}
                                   }

               for(int i=1;i<=nat;++i) {
if(
    fabs(x1[i]*r1(1)+y1[i]*r2(1)+z1[i]*r3(1)-da[i])>SMALLPOSITIONDEVIATION
||  fabs(x1[i]*r1(2)+y1[i]*r2(2)+z1[i]*r3(2)-db[i])>SMALLPOSITIONDEVIATION
||  fabs(x1[i]*r1(3)+y1[i]*r2(3)+z1[i]*r3(3)-dc[i])>SMALLPOSITIONDEVIATION
  )
{fprintf(stderr,"Warning mcdiff: da db dc and dr1 dr2 dr3 inconsistent for nonmagnetic ion number %i \n",i);
 if(use_dadbdc==0){ fprintf (stderr,"using dr1 dr2 dr3 and recalculating da db dc...\n");}
             else {fprintf (stderr,"using da db dc and recalculating dr1 dr2 dr3...\n");}                
 }
        if(use_dadbdc==0){       da[i]= x1[i]*r1(1)+y1[i]*r2(1)+z1[i]*r3(1);
                                 db[i]= x1[i]*r1(2)+y1[i]*r2(2)+z1[i]*r3(2);
                                 dc[i]= x1[i]*r1(3)+y1[i]*r2(3)+z1[i]*r3(3);
                         }
        else
                         {       x1[i]= (da[i]*rez1(1)+db[i]*rez1(2)+dc[i]*rez1(3))/2/PI;
                                 y1[i]= (da[i]*rez2(1)+db[i]*rez2(2)+dc[i]*rez2(3))/2/PI;
                                 z1[i]= (da[i]*rez3(1)+db[i]*rez3(2)+dc[i]*rez3(3))/2/PI;
                         }
                                   }
             }


eps=0.5*(eps+eps.Transpose());

printf("  Strain Tensor epsilon (coordinate system ijk defined by  j||b, k||(a x b) and i normal to k and j\n");
myPrintMatrix(stdout,eps);   

r1s=Vector(1,3);r2s=Vector(1,3);r3s=Vector(1,3);
r1s=r1;r2s=r2;r3s=r3;
rtoijk=Matrix(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the ijk coordinate system

// formulas
//  j||b, k||(a x b) and i normal to j and k
// (ijk form an euclidian righthanded coordinate system)
//  ri=ris(1)*a+ris(2)*b+ris(3)*c = ri(1)*i+ri(2)*j+ri(3)*k    (1)
//  (and then finally consider magnetic supercell and transform ri=ri*nri)
// so how to get ri(1,..,3): calculate the matrix rtoijk, which should obey
// ri()=rtoijk*ris() for i=1,2,3
// to get the components of this matrix we multiply (1) by i,j,k and get:
// ri(1)=ris(1)*(a.i)+ris(2) (b.i)+ ris(3)*(c.i)
// ri(2)=ris(1)*(a.j)+ris(2) (b.j)+ ris(3)*(c.j)
// ri(3)=ris(1)*(a.k)+ris(2) (b.k)+ ris(3)*(c.k)
// using j||b, k||(a x b) and i normal to j and k  we get
//i=(bx(axb))/(|a||b|^2sin(gamma)sin(angl(b,axb))
// ri(1)=ris(1)*(a.i)                 + ris(3)*(c.i)
// ri(2)=ris(1)*(a.b)/|b|+ris(2) |b|  + ris(3)*(c.b)/|b|
// ri(3)=                             + ris(3)*(c.(axb))/(|a||b|sin(gamma))
// note (a.i)=(a.(bx(axb))/|bx(axb)|=|a|sin(gamma)
// i.e.
//         | |a|sin(gamma) 0         (c.i)                         |
// rtoijk= | |a|cos(gamma) |b|       |c|cos(alpha)                 |
//         | 0             0         (c.(axb))/(|a||b|sin(gamma))  |
//
// to get (c.i) we write in components
// a=|a|(sin(gamma),cos(gamma),0)
// c=|c|(eps,cos(alpha),delt)
// (a.c)=|a||c|cos(beta)=|a||c|(eps*sin(gamma)+cos(gamma)*cos(alpha)
// --> eps=(cos(beta)-cos(gamma)cos(alpha))/sin(gamma)
//  (c.i)=|c|*eps
//  delta and (c.k) we get from the condition that length of c is |c|.

if (gamma>180||gamma<=0){fprintf(stderr,"ERROR mcdiff: gamma must be between 0 and 180 degrees\n");exit(EXIT_FAILURE);}
rtoijk(1,1)=a*sin(gamma*PI/180);
rtoijk(2,1)=a*cos(gamma*PI/180);
rtoijk(3,1)=0;

rtoijk(1,2)=0;
rtoijk(2,2)=b;
rtoijk(3,2)=0;

rtoijk(1,3)=c*(cos(beta*PI/180)-cos(gamma*PI/180)*cos(alpha*PI/180))/sin(gamma*PI/180);
if (fabs(rtoijk(1,3))>c){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoijk(2,3)=c*cos(alpha*PI/180);
rtoijk(3,3)=c*c-rtoijk(1,3)*rtoijk(1,3)-rtoijk(2,3)*rtoijk(2,3);
if (rtoijk(3,3)<=0){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoijk(3,3)=sqrt(rtoijk(3,3));
// rtoijk columns are a b c lattice vectors in terms of ijk euclidean system

// calculate reciprocal lattice of rtoijk
rtoijk_rez=Matrix(1,3,1,3);
rtoijk_rez=rezcalc(rtoijk);

// finally compute the components of r1 r2 r3 (magnetic unit cell) in terms of ijk coordinate system
r1=(rtoijk*r1)*(double)nr1;
r2=(rtoijk*r2)*(double)nr2;
r3=(rtoijk*r3)*(double)nr3;

// consider the strain tensor
r1+=eps*r1;
r2+=eps*r2;
r3+=eps*r3;


// transform also Projection vector
Pxyz=Vector (1,3);
Pxyz=(rtoijk*P);
// P/=Norm(Pxyz);Pxyz/=Norm(Pxyz); // normalise to length 1 : removed 14.10.2011 to be able to calculate different degrees of beam polarisation
if(Norm(Pxyz)>1){fprintf(stderr,"Warning mcdiff: length of polarization vector |P|>1 ... taking full polarised beam, i.e. normalising length of P to |P|=1\n");
                 P/=Norm(Pxyz);Pxyz/=Norm(Pxyz); }  // normalize only if |P|>1
printf("#Length of Polarization Vector P (beam polarisation for calculation of I+,I-,MSF.P): |P|=%g \n",Norm(Pxyz));

//r1(1) = a * r1(1) * nr1;
//r2(1) = a * r2(1) * nr2;
//r3(1) = a * r3(1) * nr3; 
//r1(2) = b * r1(2) * nr1;
//r2(2) = b * r2(2) * nr2;
//r3(2) = b * r3(2) * nr3;
//r1(3) = c * r1(3) * nr1;
//r2(3) = c * r2(3) * nr2;
//r3(3) = c * r3(3) * nr3;
// input section 4 *********************************************************

if (nr1 == 0){fprintf(stderr,"ERROR mcdiff: nr1 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (nr2 == 0){fprintf(stderr,"ERROR mcdiff: nr2 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (nr3 == 0){fprintf(stderr,"ERROR mcdiff: nr3 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (T <= 0){fprintf(stderr,"ERROR mcdiff: Temperature read from input file T=%g < 0\n",T);exit(EXIT_FAILURE);}
printf ("     section 4 - nr1=%i nr2=%i nr3=%i\n",nr1,nr2,nr3);
printf ("                 nat=%i magnetic atoms\n",natmagnetic);

n = nr1 * nr2 * nr3 * nat + natmagnetic; //atoms in der magnetic unit cell
jjjpars = new jjjpar * [n+1];

printf("                 reading magnetic atoms and moments ...\n");

mfcf mfields(1,1,1,natmagnetic,MAX_NOF_MF_COMPONENTS); // MAX_NOF_MF_COMPONENTS is maximum of nofmfcomponents - we take it here !
mfields.clear();
int maxmfcomponents=3;// for printout of mcdiff.mf we check what is the largest nofmfcomponents in mcdiff.in
for(int i=1;i<=natmagnetic;++i){
                            instr[0]='#';
                            while(instr[strspn(instr," \t")]=='#'){pos=ftell(fin);
                                                                   if(feof(fin)==1){fprintf(stderr,"mcdiff Error: end of file before all magnetic atoms could be read\n");exit(EXIT_FAILURE);}
                                                                  fgets(instr,MAXNOFCHARINLINE,fin);
                                                                  }
			     // get sipffilename out of "{filename}   ..."

                            if(instr[strspn(instr," \t")]!='{'){fprintf(stderr,"ERROR mcdiff: magnetic atom line has to start with '{'\n");exit (EXIT_FAILURE);}
			    if (strchr(instr,'}')==NULL){fprintf(stderr,"ERROR mcdiff: no '}' found after filename for magnetic atom %s\n",instr);exit (EXIT_FAILURE);}

                            instr[strspn(instr," \t")]='=';
			    extract(instr,"",sipffilename,(size_t)MAXNOFCHARINLINE,1);
                            if(strchr(sipffilename,'}')!=NULL){*strchr(sipffilename,'}')='\0';}
                            if(strchr(sipffilename,' ')!=NULL){*strchr(sipffilename,' ')='\0';}
                            if(strchr(sipffilename,'\t')!=NULL){*strchr(sipffilename,'\t')='\0';}
                            //printf("%s\n",sipffilename);

                             // read the rest of the line and split into numbers
                            fseek(fin,pos+strchr(instr,'}')-instr+1,SEEK_SET); 
                            j=inputline(fin,numbers);
   //MR23.10.2022 change operator sequence from Sa La Sb Lb Sc Lc --------
   //                                        to Sa Sb Sc La Lb Lc
   //                                           
   double dum; dum=numbers[9+2];numbers[9+2]=numbers[9+4];numbers[9+4]=numbers[9+5];numbers[9+5]=numbers[9+3];numbers[9+3]=dum;

if(
 fabs( numbers[4]- (numbers[1]*rez1(1)+numbers[2]*rez1(2)+numbers[3]*rez1(3))/2/PI)>SMALLPOSITIONDEVIATION
||  fabs(numbers[5]- (numbers[1]*rez2(1)+numbers[2]*rez2(2)+numbers[3]*rez2(3))/2/PI)>SMALLPOSITIONDEVIATION
||  fabs(numbers[6]- (numbers[1]*rez3(1)+numbers[2]*rez3(2)+numbers[3]*rez3(3))/2/PI)>SMALLPOSITIONDEVIATION
  )
{fprintf(stderr,"Warning mcdiff: da db dc and dr1 dr2 dr3 inconsistent for magnetic ion number %i \n",i);
 if(use_dadbdc==0){ fprintf (stderr,"using dr1 dr2 dr3 and recalculating da db dc...\n");}
                    else {fprintf (stderr,"using da db dc and recalculating dr1 dr2 dr3...\n");}                
 }

                      
if(use_dadbdc!=0)        {       numbers[4]= (numbers[1]*rez1(1)+numbers[2]*rez1(2)+numbers[3]*rez1(3))/2/PI;
                                 numbers[5]= (numbers[1]*rez2(1)+numbers[2]*rez2(2)+numbers[3]*rez2(3))/2/PI;
                                 numbers[6]= (numbers[1]*rez3(1)+numbers[2]*rez3(2)+numbers[3]*rez3(3))/2/PI;
                         }
                            if (j<9) {fprintf(stderr,"ERROR mcdiff: too few parameters for magnetic atom %i: %s\n",i,instr);exit(EXIT_FAILURE);}
                             // determine jxc .... dimension of exchange field if present >>>>>>>>>>>>>>>>
                            long int currentpos=ftell(fin);instr[0]='#';int jxc;
                            while(instr[strspn(instr," \t")]=='#'&&feof(fin)==0){pos=ftell(fin);fgets(instr,MAXNOFCHARINLINE,fin);}
                            if (strchr(instr,'>')==NULL||instr[strspn(instr," \t")]=='#')
                             {jxc=1;} // no ">" found --> do dipole approx
                             else          
                             {fseek(fin,pos+strchr(instr,'>')-instr+1,SEEK_SET); 
                              jxc=inputline(fin,numbers1);if(verbose)printf("dimension of mf = %i\n",jxc);
                             }
                             fseek(fin,currentpos,SEEK_SET); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             jjjpars[i]=new jjjpar((double)numbers[4] / nr1,(double)numbers[5] / nr2,(double)numbers[6] / nr3, sipffilename,jxc);
                             //J[i]=-1;
                             //J[i]         1   0    -1   -2   -3
                             //FF_type      1  -2    +2   +3   -3
                             (*jjjpars[i]).FF_type=+2;
                             (*jjjpars[i]).save_sipf("./results/_");// save read single ion parameter file
                              // store moment and components of S and L (if given)
                              for(k=7;k<=j&&k<=15;++k){(*jjjpars[i]).mom(k-6) = numbers[k];}
                              if((*jjjpars[i]).gJ==0){if(j>=15){(*jjjpars[i]).FF_type=+3;//J[i]=-2; // do not use input moment but spin and angular momentum for calculation
                                                                // do some consistency checks
                                                                if (fabs((*jjjpars[i]).mom(1)-2*(*jjjpars[i]).mom(4)-(*jjjpars[i]).mom(5))/(fabs((*jjjpars[i]).mom(1))+1.0)>0.001){fprintf(stderr,"Warning mcdiff: a-component magnetic moment=%g and <La>+2<Sa>=%g not consistent for atom %i - setting moment=<La>+2<Sa>\n",(*jjjpars[i]).mom(1),2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5),i);}
                                                                if (fabs((*jjjpars[i]).mom(2)-2*(*jjjpars[i]).mom(6)-(*jjjpars[i]).mom(7))/(fabs((*jjjpars[i]).mom(2))+1.0)>0.001){fprintf(stderr,"Warning mcdiff: b-component magnetic moment=%g and <Lb>+2<Sb>=%g not consistent for atom %i - setting moment=<Lb>+2<Sb>\n",(*jjjpars[i]).mom(2),2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7),i);}
                                                                if (fabs((*jjjpars[i]).mom(3)-2*(*jjjpars[i]).mom(8)-(*jjjpars[i]).mom(9))/(fabs((*jjjpars[i]).mom(3))+1.0)>0.001){fprintf(stderr,"Warning mcdiff: c-component magnetic moment=%g and <Lc>+2<Sc>=%g not consistent for atom %i - setting moment=<Lc>+2<Sc>\n",(*jjjpars[i]).mom(3),2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9),i);}
                                                                (*jjjpars[i]).mom(1)=2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5);
                                                                (*jjjpars[i]).mom(2)=2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7);
                                                                (*jjjpars[i]).mom(3)=2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9);
                                                                }
                                                           else {(*jjjpars[i]).FF_type=+2;}//J[i]=-1;} // just use spin formfactor
                                                      }
                            instr[0]='#';
                            while(instr[strspn(instr," \t")]=='#'&&feof(fin)==0){pos=ftell(fin);fgets(instr,MAXNOFCHARINLINE,fin);}
                            if (strchr(instr,'>')==NULL||instr[strspn(instr," \t")]=='#')
                             {fseek(fin,pos,SEEK_SET);} // no ">" found --> do dipole approx
                             else          
                             {Vector Qvec(1,3);Qvec=0;ComplexVector Mq(1,3);
                              fseek(fin,pos+strchr(instr,'>')-instr+1,SEEK_SET); 
                              j=inputline(fin,numbers);if(verbose)printf("dimension of mf = %i\n",j);
                              if(j>maxmfcomponents){maxmfcomponents=j;}
                              if(j>mfields.nofcomponents){fprintf(stderr,"ERROR mcdiff: number of exchange field components too large (%i>%i) recompile with larger MAX_NOF_MF_COMPONENTS\n",j,mfields.nofcomponents);exit(EXIT_FAILURE);}
                              Vector gjmbHxc(1,j);for(k=1;k<=j;++k){gjmbHxc(k)=numbers[k];mfields.mf(1,1,1)(mfields.nofcomponents*(i-1)+k)=gjmbHxc(k);}
                              (*jjjpars[i]).eigenstates(gjmbHxc,H,T); // calculate eigenstates
                              (*jjjpars[i]).Icalc_parameter_storage_init(gjmbHxc,H,T);// initialise parameter storage for Icalc
                              // do some consistency checks
                               ComplexMatrix Icalcpars((*jjjpars[i]).Icalc_parstorage.Rlo(),(*jjjpars[i]).Icalc_parstorage.Rhi(),(*jjjpars[i]).Icalc_parstorage.Clo(),(*jjjpars[i]).Icalc_parstorage.Chi());
                                             Icalcpars=(*jjjpars[i]).Icalc_parstorage;
                              // check if M(Q) works for this sipf module - if yes do beyond calculation for this
                              // ion 
                              if((*jjjpars[i]).MQ(Mq,Qvec))
                              {Vector moment(1,3),L(1,3),S(1,3);
                               // check if mcalc is present !!!, if yes:
                               if((*jjjpars[i]).mcalc(moment,T,gjmbHxc,H,Icalcpars))
                               {if (fabs((*jjjpars[i]).mom(1)-moment(1))>0.001){fprintf(stderr,"Warning mcdiff: a-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                                if (fabs((*jjjpars[i]).mom(2)-moment(2))>0.001){fprintf(stderr,"Warning mcdiff: b-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                                if (fabs((*jjjpars[i]).mom(3)-moment(3))>0.001){fprintf(stderr,"Warning mcdiff: c-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                                (*jjjpars[i]).mom(1)=moment(1);
                                (*jjjpars[i]).mom(2)=moment(2);
                                (*jjjpars[i]).mom(3)=moment(3);
                               }
                               // check if <L> and <S> are present, if yes   
                              if((*jjjpars[i]).FF_type==+3)//if(J[i]==-2)
   		              {(*jjjpars[i]).FF_type=-3;//J[i]=-3;
                              // check if Lcalc and Scalc are present
                               if((*jjjpars[i]).Lcalc(L,T,gjmbHxc,H,Icalcpars)&&
                              (*jjjpars[i]).Scalc(S,T,gjmbHxc,H,Icalcpars))
                               {if (fabs((*jjjpars[i]).mom(4)-S(1))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <Sa> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(4)=S(1);
   			       if (fabs((*jjjpars[i]).mom(5)-L(1))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <La> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(5)=L(1);
   			       if (fabs((*jjjpars[i]).mom(6)-S(2))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <Sb> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(6)=S(2);
   			       if (fabs((*jjjpars[i]).mom(7)-L(2))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <Lb> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(7)=L(2);
   			       if (fabs((*jjjpars[i]).mom(8)-S(3))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <Sc> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(8)=S(3);
   			       if (fabs((*jjjpars[i]).mom(9)-L(3))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <Lc> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(9)=L(3);
                                }
			      }
			      else
			      { (*jjjpars[i]).FF_type=-2;   
                              //J[i]=0; // J=0 tells that full calculation should be done for this ion using 
                                      // for dip intensities Ma Mb and Mc
                               }

                             (*jjjpars[i]).checkFFcoeffnonzero(4);
                             (*jjjpars[i]).checkFFcoeffnonzero(6);

 			      }else{fprintf(stderr,"WARNING mcdiff: exchange fields given in mcdiff.in (probably to go beyond dipole approximation) but MQ function not implemented for ion in file %s - switching to dipole approximation\n",sipffilename);}
                              }

                             if((*jjjpars[i]).SLR==0){fprintf(stderr,"WARNING mcdiff: SCATTERINGLENGTHREAL not found or zero in file %s\n",sipffilename);}
//                             if((*jjjpars[i]).gJ==0){fprintf(stderr,"WARNING mcdiff: GJ not found or zero in file %s - gJ=0 means Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc !\n",sipffilename);}
                             (*jjjpars[i]).checkFFcoeffnonzero(0);
                             (*jjjpars[i]).checkFFcoeffnonzero(2);

                           }
  fclose(fin);

 mfields.resetnofc(maxmfcomponents);

// print spinconfiguration to mcdiff.sps  (useful for viewing)
print_sps(*this);
print_mf(*this,mfields);

//now insert also nonmagnetic elements into the unit cell
int ncryst,na,nb,nc;
ncryst = natmagnetic;
for(na = 1;na<=nr1;++na){
 for(nb = 1;nb<=nr2;++nb){
  for(nc = 1;nc<=nr3;++nc){
   if(nat!=0){
    for(int i=1;i<=nat;++i){
      ++ncryst;
      jjjpars[ncryst]=new jjjpar((na + x1[i] - 1) / nr1,(nb + y1[i] - 1) / nr2,(nc + z1[i] - 1) / nr3,sl1r[i],sl1i[i],dwf1[i]);
      (*jjjpars[ncryst]).FF_type=+1;//J[ncryst]=1;
      (*jjjpars[ncryst]).mom=0;
      (*jjjpars[ncryst]).gJ=0;
      }
    }
}}}
// finally do reziprocal lattice from r1 r2 r3  (which were recomputed to yield magnetic unit cell in terms of ijk)
   rezcalc (r1, r2, r3, rez1, rez2, rez3);

  save();
}

//kopier-konstruktor 
inimcdiff::inimcdiff (const inimcdiff & p)
{ colcod=new int[NOF_OUT_VARIABLES+1];for(int i=0;i<NOF_OUT_VARIABLES;++i){colcod[i]=cc[i];}
  colhead=new char *[COLHEADDIM+1];for(int i=0;i<=COLHEADDIM;++i){colhead[i]=new char [strlen(ch[i])+1];strcpy(colhead[i],ch[i]);}
  verbose=p.verbose;
  savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);
  outfilename= new char [strlen(p.outfilename)+1];
  strcpy(outfilename,p.outfilename);
  infile= new char [strlen(p.infile)+1];
  strcpy(infile,p.infile);
    prefix= new char [strlen(p.prefix)+1]; strcpy(prefix,p.prefix);
  unitcellstr= new char [strlen(p.unitcellstr)+1]; strcpy(unitcellstr,p.unitcellstr);
  T=p.T;  
  lambda=p.lambda;
  thetamax=p.thetamax;
  ovalltemp=p.ovalltemp;
  lorenz=p.lorenz;
  nat=p.nat;
  x1=new float[nat+1];y1=new float[nat+1];z1=new float[nat+1];
  da=new float[nat+1];db=new float[nat+1];dc=new float[nat+1];
  sl1r=new float[nat+1];sl1i=new float[nat+1];;dwf1=new float[nat+1];
  for (int i=1;i<=nat;++i){
  x1[i]=p.x1[i];y1[i]=p.y1[i];z1[i]=p.z1[i];
  da[i]=p.da[i];db[i]=p.db[i];dc[i]=p.dc[i];
  sl1r[i]=p.sl1r[i];sl1i[i]=p.sl1i[i];dwf1[i]=p.dwf1[i];
                          } 
  n=p.n;
  nofatoms=p.nofatoms;
  natmagnetic=p.natmagnetic;
  P=p.P;Pxyz=p.Pxyz;eps=p.eps;r1=p.r1;r2=p.r2;r3=p.r3;rez1=p.rez1;rez2=p.rez2;rez3=p.rez3;
  a=p.a;b=p.b;c=p.c;alpha=p.alpha;beta=p.beta;gamma=p.gamma;
  use_dadbdc=p.use_dadbdc;
  nr1=p.nr1;nr2=p.nr2;nr3=p.nr3;H=p.H;
  jjjpars = new jjjpar * [n+1];
  r1s=Vector(1,3);r2s=Vector(1,3);r3s=Vector(1,3);
  r1s=p.r1s;r2s=p.r2s;r3s=p.r3s;
  rtoijk=Matrix(1,3,1,3);
  rtoijk_rez=Matrix(1,3,1,3);
  rtoijk=p.rtoijk;
  rtoijk_rez=p.rtoijk_rez;
 //mf=mfcf(1,1,1,nofatoms,nofcomponents);mf=p.mf;T=p.T;
  //nofhkls=p.nofhkls;
  
}

//destruktor
inimcdiff::~inimcdiff ()
{delete []savfilename;delete []infile;delete []unitcellstr;
 delete []prefix;delete []x1;delete []y1;delete []z1;
 delete []da;delete []db;delete []dc;
 delete []sl1r;delete []sl1i;delete []dwf1;
//  for (i=1;i<=n;++i){delete jjjpars[i];}
//  delete []jjjpars*;
delete []colcod;
for(int i=0;i<=COLHEADDIM;++i){delete []colhead[i];}
delete []colhead;

/* int i;
  for (i=1;i<=nofhkls;++i) 
   { delete []hkls[i];}
   delete []hkls;
   delete  []hklfile_start_index;
 */
}
