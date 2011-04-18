/***********************************************************************
 *
 * mcdiff - program to calculate neutron and magnetic xray diffraction
 *
 * reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R
 ***********************************************************************/
#include <mcdiff.h>
#include "mcdiff_intcalc.c"
#include "mcdiff_output.c"


// hauptprogramm
int main (int argc, char **argv)
{ FILE * fin, * fin_coq, * fout;
  float ovalltemp,thetamax,lambda,a=0,b=0,c=0,alpha=0,beta=0,gamma=0;
  double T=0;
  int i,j,k,n,lorenz,nat, nofatoms,nr1=0,nr2=0,nr3=0,natmagnetic;
  int colcode[40]; // field to store code for assigning type of data to columns of output
  long int pos=0;
  int use_dadbdc=0;
  char instr[MAXNOFCHARINLINE+1];
  char cffilename[MAXNOFCHARINLINE+1];
  char unitcellstr[MAXNOFCHARINLINE+1];
  float numbers[70];numbers[0]=70;
  Vector r1(1,3),r2(1,3),r3(1,3),H(1,3),P(1,3);
  Vector rez1(1,3),rez2(1,3),rez3(1,3);
fprintf(stderr,"***********************************************************************\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* mcdiff - program to calculate neutron and magnetic xray diffraction\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
fprintf(stderr,"***********************************************************************\n");

  // check command line
   if (argc > 1)
    { if (strcmp(argv[1],"-h")==0)
      {printf (" \n \
                use as: mcdiff [hkllistfilename]\n \
		- for format of input file mcdiff.in see mcphase manual\n \
                - optional an hkl list can be given in file hkllistfilename-in\n \
                  this case the program computes reflections in this list.\n \
                  if it is a 3column list,the azimuth dependence of\n \
                  magnetic scattering is calculated.If neutron intensities\n \
                  are listed in column 4 the program computes\n \
                  rpvalue=100 * sum_i |Icalc_i-Iexp_i|/sum_i |Iobs_i|\n \
                  and does not output the azimuth dependence \n \
                  if in addition experimental errors are given in column 5,\n \
                  chisquared=1/N *sum_i (Icalc_i - Iexp_i)^2/err_i^2 is\n \
                  calculated also.\n \
                - results are saved in results/mcdiff.out\n");
        exit (1);}
    }

// check if directory results exists and can be written to ...
   fout = fopen_errchk ("./results/mcdiff.out", "a");fclose(fout);

// test to test threej function
/* if (argc>6) {printf ("cint(%g)=%i\n", strtod(argv[6],NULL),cint(strtod(argv[6],NULL)));
		// test matpack routine  
                    int n,ndim=20; double thrcof[20];int errflag;
		    double min,max;
		     ThreeJSymbolM	(strtod(argv[1],NULL),
	             strtod(argv[2],NULL),
	             strtod(argv[3],NULL),
	             strtod(argv[4],NULL),
	             min,max, thrcof, ndim, 
			 errflag);
		     n=(int)(strtod(argv[5],NULL)-min);	 

              printf ("threej symbol=%g=%g\n",
              threej(strtod(argv[1],NULL),
	             strtod(argv[2],NULL),
	             strtod(argv[3],NULL),
	             strtod(argv[4],NULL),
	             strtod(argv[5],NULL),
	             strtod(argv[6],NULL)
	             ),thrcof[n]);
		     return 0;}
*/
// test spherical harmonics
/*if (argc>3) {
  int l,m; double theta,phi;
  l=(int)strtod(argv[1],NULL);
  m=(int)strtod(argv[2],NULL);
  theta=strtod(argv[3],NULL);
  phi=strtod(argv[4],NULL);
  
  printf("Y%i%i(%g,%g)=%g %+g i\n",l,m,theta,phi,real(SphericalHarmonicY (l,m,theta,phi)),imag(SphericalHarmonicY (l,m,theta,phi)));
 return 0;}
*/


fin_coq = fopen_errchk ("./mcdiff.in", "rb");
 fprintf(stdout,"\n reading file mcdiff.in\n\n");
fout = fopen_errchk ("./results/_mcdiff.in", "w"); //copy input file parameters to results
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

// input section 1 *******************************************************

  instr[0]='#';colcode[10]=1;colcode[11]=0;
 while (instr[strspn(instr," \t")]=='#'&&strstr (instr, "%SECTION 2%")==NULL) // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
    if (pos==-1) 
       {fprintf(stderr,"Error mcdiff: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   extract(instr,"nofatoms",nofatoms);  
   extract(instr,"lambda", lambda);
   extract(instr, "thetamax", thetamax);
   extract(instr, "nat", nat);
   extract(instr, "ovalltemp", ovalltemp);
   extract(instr, "lorentz", lorenz);
   extract(instr, "out10",colcode[10]);
   extract(instr, "out11",colcode[11]);
   extract(instr, "Pa",P(1));
   extract(instr, "Pb",P(2));
   extract(instr, "Pc",P(3));
  }
  fseek(fin_coq,pos,SEEK_SET); 

if (lorenz == 0){fprintf(stderr,"Warning mcdiff: read lorentz=0, will calculate no Lorentzfactor.\n");}
if (lambda == 0){fprintf(stderr,"ERROR mcdiff: no wavelength lambda given or line does not start with # in section 1\n");exit(EXIT_FAILURE);}
if (thetamax == 0){fprintf(stderr,"ERROR mcdiff: no thetamax given or line does not start with # in section 1\n");exit(EXIT_FAILURE);}
printf("     section 1 - lambda=%g A thetamax= %g deg\n",lambda, thetamax);
printf("                 ovalltemp=%g A^2 lorentz-type=%i\n",ovalltemp,lorenz);
printf("                 output: column 10=%s column 11=%s\n",colheader[colcode[10]],colheader[colcode[11]]);

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
fprintf(fout,"#            5.....neutron TOF powder cyl. sample - d-pattern normal scaled\n");
fprintf(fout,"#! out10=%i    type of desired output in column 10 and 11 of mcdiff.out\n",colcode[10]);
fprintf(fout,"#! out11=%i    (optional) default is NSF in column 10 and LF in column 11\n",colcode[11]);
for(i=0;i<=20;++i){
fprintf(fout,"#            %i....%s\n",i,colheader[i]);
                   }
fprintf(fout,"#\n");
fprintf(fout,"#           In the above the intensities I+ and I- are the spinflip and nonspinflip intensities\n");
fprintf(fout,"#           in a polarised neutron experiment:\n");
fprintf(fout,"#            I+-=LF exp(-OTF Q^2/8pi^2) \n");
fprintf(fout,"#                    [ |NSF/NB|^2 + 3.65/4pi (|MSF/NB|^2-i(MSF x MSF*).P) \n");
fprintf(fout,"#                        +-  sqrt(3.65/4pi)/NB^2 (NSF (MSF*.P) + NSF* (MSF.P)]\n");
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"#             For some of the above options we need the\n");
fprintf(fout,"#! Pa=%8.4f   Components of Projection Vector P=(Pa * a + Pb * b + Pc *c)/Norm(Pa * a + Pb * b + Pc *c)\n",P(1));
fprintf(fout,"#! Pb=%8.4f\n",P(2));
fprintf(fout,"#! Pc=%8.4f\n",P(3));
fprintf(fout,"#\n");
fprintf(fout,"#\n");

// input section 2 *********************************************************

  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&a==0&&b==0&&c==0) // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
    if (pos==-1) 
       {fprintf(stderr,"Error mcdiff: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   extract(instr, "nat", nat);
   extract(instr, "a", a);
   extract(instr, "b", b);
   extract(instr, "c", c);
   extract(instr, "alpha", alpha);
   extract(instr, "beta", beta);
   extract(instr, "gamma", gamma);
   extract(instr, "use_dadbdc",use_dadbdc);
  }
  fseek(fin_coq,pos,SEEK_SET); 
  printf ("     section 2 - nat=%i\n",nat);
  float *x1;x1=new float[nat+1];float*y1;y1=new float[nat+1];float*z1;z1=new float[nat+1];
  float *da;da=new float[nat+1];float*db;db=new float[nat+1];float*dc;dc=new float[nat+1];
  float *sl1r;sl1r=new float[nat+1];float*sl1i;sl1i=new float[nat+1];float *dwf1;dwf1=new float[nat+1];

fprintf(fout,"# %%SECTION 2%% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL\n");
fprintf(fout,"#\n");
fprintf(fout,"#\n");
fprintf(fout,"#! nat=%i      number of nonmagnetic atoms in primitive crystalographic unit cell\n",nat);
fprintf(fout,"#\n");
fprintf(fout,"# it follows a list of nat lines with nonmagnetic atoms\n");
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

  if (nat!=0){ for(i=1;i<=nat;++i) { pos=ftell(fin_coq); 
                                     n=inputline(fin_coq,numbers);
                                     if (n==0) {if(feof(fin_coq)==true){fprintf(stderr,"Error mcdiff: end of input file in section 2\n");exit (EXIT_FAILURE);}
                                                fseek(fin_coq,pos,SEEK_SET); 
                                                fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                                                if(strstr (instr, "%%SECTION 3%%")!=NULL){fprintf (stderr,"ERROR mcdiff: Section 3 started before all nat=%i atoms of crystallographic unit cell were listed !\n",nat);exit (EXIT_FAILURE);}
                                               --i;}
                                     else      {if (n<9) {fprintf (stderr,"ERROR mcdiff: Section 2 - Nonmagnetic Atoms: too few positional parameters for atom %i!\n",i);exit (EXIT_FAILURE);}
                                                sl1r[i]=numbers[1];sl1i[i]=numbers[2]; x1[i] = numbers[6]; y1[i] = numbers[7]; z1[i] = numbers[8];dwf1[i]=numbers[9];
                                                                                       da[i] = numbers[3]; db[i] = numbers[4]; dc[i] = numbers[8];
                                                printf("                 sl=%g%+gi 10^-12cm at %g*r1%+g*r2%+g*r3 DWF=%g A^2\n",sl1r[i],sl1i[i],x1[i],y1[i],z1[i],dwf1[i]);
                                               }
                                    }
              }

// input section 3 *********************************************************
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&nr1*nr2*nr3==0) 
  { pos=ftell(fin_coq); 
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   if(a==0)extract(instr, "a", a);
   if(b==0)extract(instr, "b", b);
   if(c==0)extract(instr, "c", c);
   if(alpha==0)extract(instr, "alpha", alpha);
   if(beta==0)extract(instr, "beta", beta);
   if(gamma==0)extract(instr, "gamma", gamma);
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
  fseek(fin_coq,pos,SEEK_SET); 
if (a == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant a given in section 3 or line does not start with # or nat too small: \n%s\n",instr);exit (EXIT_FAILURE);}
if (b == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant b given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
if (c == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant c given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
printf("     section 3 - a=%g A  b=%g A c=%g A alpha=%g  beta=%g gamma=%g\n",a,b,c,alpha,beta,gamma);
sprintf(unitcellstr," a= %g A  b= %g A c= %g A  alpha=%g  beta=%g gamma=%g\n",a,b,c,alpha,beta,gamma);
printf("                 r1= %5.3ga + %5.3gb + %5.3gc\n", r1(1), r1(2), r1(3));
printf("                 r2= %5.3ga + %5.3gb + %5.3gc\n", r2(1), r2(2), r2(3));
printf("                 r3= %5.3ga + %5.3gb + %5.3gc\n", r3(1), r3(2), r3(3));

//printf("                    / %5.3ga \\     / %5.3ga \\     / %5.3ga \\    x||c \n", r1(1), r2(1), r3(1));
//printf("                 r1=| %5.3gb |  r2=| %5.3gb |  r3=| %5.3gb |    y||a\n", r1(2), r2(2), r3(2));
//printf("                    \\ %5.3gc /     \\ %5.3gc /     \\ %5.3gc /    z||b\n", r1(3), r2(3), r3(3));

  double da1,db1,dc1,dd;
  rezcalc (r1, r2, r3, rez1, rez2, rez3);
  if (nat!=0){ for(i=1;i<=nat;++i) {
   // calculate da db dc from dr1 dr2 dr3 and print to results/_mcdiff.in
       da1= x1[i]*r1(1)+y1[i]*r2(1)+z1[i]*r3(1)-da[i];
       db1= x1[i]*r1(2)+y1[i]*r2(2)+z1[i]*r3(2)-db[i];
       dc1= x1[i]*r1(3)+y1[i]*r2(3)+z1[i]*r3(3)-dc[i];
       dd=sqrt(da1*da1+db1*db1+dc1*dc1);
                                 da1=x1[i]- (da[i]*rez1(1)+db[i]*rez1(2)+dc[i]*rez1(3))/2/PI;
                                 db1=y1[i]- (da[i]*rez2(1)+db[i]*rez2(2)+dc[i]*rez2(3))/2/PI;
                                 dc1=z1[i]- (da[i]*rez3(1)+db[i]*rez3(2)+dc[i]*rez3(3))/2/PI;
       dd+=sqrt(da1*da1+db1*db1+dc1*dc1);
       if(dd>SMALL){fprintf (stderr,"Warning: atomic positions da db dc and dr1 dr2 dr3 inconsistent !\n");
                    fprintf (stderr,"         use_dadbdc=%i\n",use_dadbdc);
                    if(use_dadbdc==0){ fprintf (stderr,"using dr1 dr2 dr3 and recalculating da db dc...\n");}
                    else {fprintf (stderr,"using da db dc and recalculating dr1 dr2 dr3...\n");}
                   i=nat;}
                                   }

               for(i=1;i<=nat;++i) {
        if(use_dadbdc==0){       da[i]= x1[i]*r1(1)+y1[i]*r2(1)+z1[i]*r3(1);
                                 db[i]= x1[i]*r1(2)+y1[i]*r2(2)+z1[i]*r3(2);
                                 dc[i]= x1[i]*r1(3)+y1[i]*r2(3)+z1[i]*r3(3);
                         }
        else
                         {       x1[i]= (da[i]*rez1(1)+db[i]*rez1(2)+dc[i]*rez1(3))/2/PI;
                                 y1[i]= (da[i]*rez2(1)+db[i]*rez2(2)+dc[i]*rez2(3))/2/PI;
                                 z1[i]= (da[i]*rez3(1)+db[i]*rez3(2)+dc[i]*rez3(3))/2/PI;
                         }
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
fprintf(fout,"#\n");
fprintf(fout,"#\n");


Vector r1s(1,3),r2s(1,3),r3s(1,3);
r1s=r1;r2s=r2;r3s=r3;
Matrix rtoxyz(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the xyz coordinate system
rtoxyz(1,1)=0;
rtoxyz(2,1)=a*sin(gamma*PI/180);
rtoxyz(3,1)=a*cos(gamma*PI/180);

rtoxyz(1,2)=0;
rtoxyz(2,2)=0;
rtoxyz(3,2)=b;

rtoxyz(3,3)=c*cos(alpha*PI/180);
rtoxyz(2,3)=(a*c*cos(beta*PI/180)-rtoxyz(3,3)*rtoxyz(3,1))/rtoxyz(2,1);
if (fabs(rtoxyz(2,3))>c){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoxyz(1,3)=c*c-rtoxyz(2,3)*rtoxyz(2,3)-rtoxyz(3,3)*rtoxyz(3,3);
if (rtoxyz(1,3)<=0){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoxyz(1,3)=sqrt(rtoxyz(1,3));

r1=(rtoxyz*r1)*(double)nr1;
r2=(rtoxyz*r2)*(double)nr2;
r3=(rtoxyz*r3)*(double)nr3;

// transform also Projection vector
Vector Pxyz (1,3);
Pxyz=(rtoxyz*P);
P/=Norm(Pxyz);Pxyz/=Norm(Pxyz); // normalise to length 1

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


int *J;J=new int[n+1]; // code for indicating if ion is nonmagnetic (J=1),
            // go beyond dipole approx for rare earth (J=0)
            // use magnetic moment with dipole approx (J=-1)
            // use L and S values with dipole approx (J=-2) 
            // go beyond dipole approx for gJ=0 (L and S separately)
jjjpar ** jjjpars = new jjjpar * [n+1];
double lnZ,U;
printf("                 reading magnetic atoms and moments ...\n");

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
fprintf(fout,"# Temperature,    Magnetic Field: Magnetic Unit Cell\n");
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
fprintf(fout,"#                in case of non orthogonal lattices instead of Ma Mb Mc the components Mi Mj Mk\n");
fprintf(fout,"#                have to be given, which refer to an right handed orthogonal coordinate system \n");
fprintf(fout,"#                defined by j||b, k||(a x b) and i normal to k and j\n");
fprintf(fout,"# '<Ja>' '<Jb>' '<Jc>' (optional) denote the momentum components \n");
fprintf(fout,"# 'gjmbHeffa' 'gjmbHeffb' 'gjmbHeffc' (optional line, used to go beyond dipole approx for formfactor)\n");
fprintf(fout,"#                                     denote the corresponding meanfields multiplied by \n");
fprintf(fout,"#                                     Lande factor and Bohr magneton \n");
fprintf(fout,"#\n");
fprintf(fout,"#{atom-file} da[a]  db[b]    dc[c]     dr1[r1]  dr2[r2]  dr3[r3]   <Ma>     <Mb>     <Mc> [mb] <Ja>     <Jb>     <Jc> ...\n");
fprintf(fout,"#{corresponding effective fields gjmbHeff [meV]- if passed to mcdiff only these are used for caculation (not the magnetic moments)}\n");

mfcf mfields(1,1,1,natmagnetic,51); // 51 is maximum of nofmfcomponents - we take it here !
mfields.clear();

for(i=1;i<=natmagnetic;++i){
                            instr[0]='#';J[i]=-1;
                            while(instr[strspn(instr," \t")]=='#'){pos=ftell(fin_coq);
                                                                   if(feof(fin_coq)==1){fprintf(stderr,"mcdiff Error: end of file before all magnetic atoms could be read\n");exit(EXIT_FAILURE);}
                                                                  fgets(instr,MAXNOFCHARINLINE,fin_coq);
                                                                  }
			     // get cffilename out of "{filename}   ..."

                            if(instr[strspn(instr," \t")]!='{'){fprintf(stderr,"ERROR mcdiff: magnetic atom line has to start with '{'\n");exit (EXIT_FAILURE);}
			    if (strchr(instr,'}')==NULL){fprintf(stderr,"ERROR mcdiff: no '}' found after filename for magnetic atom %s\n",instr);exit (EXIT_FAILURE);}

                            instr[strspn(instr," \t")]='=';
			    extract(instr,"",cffilename,(size_t)MAXNOFCHARINLINE);
                            if(strchr(cffilename,'}')!=NULL){*strchr(cffilename,'}')='\0';}
                            if(strchr(cffilename,' ')!=NULL){*strchr(cffilename,' ')='\0';}
                            if(strchr(cffilename,'\t')!=NULL){*strchr(cffilename,'\t')='\0';}
                            //printf("%s\n",cffilename);

                             // read the rest of the line and split into numbers
                            fseek(fin_coq,pos+strchr(instr,'}')-instr+1,SEEK_SET); 
                            j=inputline(fin_coq,numbers);
if(use_dadbdc!=0)        {       numbers[4]= (numbers[1]*rez1(1)+numbers[2]*rez1(2)+numbers[3]*rez1(3))/2/PI;
                                 numbers[5]= (numbers[1]*rez2(1)+numbers[2]*rez2(2)+numbers[3]*rez2(3))/2/PI;
                                 numbers[6]= (numbers[1]*rez3(1)+numbers[2]*rez3(2)+numbers[3]*rez3(3))/2/PI;
                         }
                            if (j<9) {fprintf(stderr,"ERROR mcdiff: too few parameters for magnetic atom %i: %s\n",i,instr);exit(EXIT_FAILURE);}
                             jjjpars[i]=new jjjpar((double)numbers[4] / nr1,(double)numbers[5] / nr2,(double)numbers[6] / nr3, cffilename);
                             (*jjjpars[i]).save_sipf("./results/_");// save read single ion parameter file
                              // store moment and components of S and L (if given)
                              for(k=7;k<=j&&k<=15;++k){(*jjjpars[i]).mom(k-6) = numbers[k];}
                              if((*jjjpars[i]).gJ==0){if(j>=15){J[i]=-2; // do not use input moment but spin and angular momentum for calculation
                                                                // do some consistency checks
                                                                if (fabs((*jjjpars[i]).mom(1)-2*(*jjjpars[i]).mom(4)-(*jjjpars[i]).mom(5))/(fabs((*jjjpars[i]).mom(1))+1.0)>0.001){fprintf(stderr,"Warning mcdiff: a-component magnetic moment=%g and <La>+2<Sa>=%g not consistent for atom %i - setting moment=<La>+2<Sa>\n",(*jjjpars[i]).mom(1),2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5),i);}
                                                                if (fabs((*jjjpars[i]).mom(2)-2*(*jjjpars[i]).mom(6)-(*jjjpars[i]).mom(7))/(fabs((*jjjpars[i]).mom(2))+1.0)>0.001){fprintf(stderr,"Warning mcdiff: b-component magnetic moment=%g and <Lb>+2<Sb>=%g not consistent for atom %i - setting moment=<Lb>+2<Sb>\n",(*jjjpars[i]).mom(2),2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7),i);}
                                                                if (fabs((*jjjpars[i]).mom(3)-2*(*jjjpars[i]).mom(8)-(*jjjpars[i]).mom(9))/(fabs((*jjjpars[i]).mom(3))+1.0)>0.001){fprintf(stderr,"Warning mcdiff: c-component magnetic moment=%g and <Lc>+2<Sc>=%g not consistent for atom %i - setting moment=<Lc>+2<Sc>\n",(*jjjpars[i]).mom(3),2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9),i);}
                                                                (*jjjpars[i]).mom(1)=2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5);
                                                                (*jjjpars[i]).mom(2)=2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7);
                                                                (*jjjpars[i]).mom(3)=2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9);
                                                                }
                                                           else {J[i]=-1;} // just use spin formfactor
                                                      }
fprintf(fout,"{%s} %8.5f %8.5f %8.5f  ",cffilename,numbers[4]*r1s(1)+numbers[5]*r2s(1)+numbers[6]*r3s(1),numbers[4]*r1s(2)+numbers[5]*r2s(2)+numbers[6]*r3s(2),numbers[4]*r1s(3)+numbers[5]*r2s(3)+numbers[6]*r3s(3));
fprintf(fout,"%8.5f %8.5f %8.5f  ",numbers[4],numbers[5],numbers[6]); // positions
fprintf(fout," %+8.5f %+8.5f %+8.5f ",(*jjjpars[i]).mom(1),(*jjjpars[i]).mom(2),(*jjjpars[i]).mom(3)); // magmoments
for(k=10;k<=j;++k){fprintf(fout," %+8.5f",numbers[k]);}
fprintf(fout,"\n");
                            instr[0]='#';
                            while(instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0){pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
                            if (strchr(instr,'>')==NULL)
                             {fseek(fin_coq,pos,SEEK_SET);} // no ">" found --> do dipole approx
                             else          
                             {J[i]=0; // J=0 tells that full calculation should be done for this ion
                              fseek(fin_coq,pos+strchr(instr,'>')-instr+1,SEEK_SET); 
                              j=inputline(fin_coq,numbers);printf("dimension of mf = %i\n",j);
                              Vector heff(1,j);for(k=1;k<=j;++k){heff(k)=numbers[k];mfields.mf(1,1,1)(51*(i-1)+k)=heff(k);}
                              if ((*jjjpars[i]).gJ==0)
 			      {J[i]=-3;fprintf(stderr,"mcdiff: gJ=0 - going beyond dipolar approximation for intermediate coupling");
   			             (*jjjpars[i]).eigenstates(heff,T); // calculate eigenstates
                                     (*jjjpars[i]).mcalc_parameter_storage_init(heff,T);// initialise parameter storage for mcalc
                               // do some consistency checks
                               ComplexMatrix mcalcpars((*jjjpars[i]).mcalc_parstorage.Rlo(),(*jjjpars[i]).mcalc_parstorage.Rhi(),(*jjjpars[i]).mcalc_parstorage.Clo(),(*jjjpars[i]).mcalc_parstorage.Chi());
                                             mcalcpars=(*jjjpars[i]).mcalc_parstorage;
                               Vector moment(1,j);(*jjjpars[i]).mcalc(moment,T,heff,lnZ,U,mcalcpars);
                               for(k=1;k<=j;++k){if (fabs((*jjjpars[i]).mom(k+3)-moment(k))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <J> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                                                  (*jjjpars[i]).mom(3+k)=moment(k);//printf("m(%i)=%g ",k,moment(k));
                                                 }
                               (*jjjpars[i]).mom(1)=2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5);
                               (*jjjpars[i]).mom(2)=2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7);
                               (*jjjpars[i]).mom(3)=2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9);
			      }
			      else
			      {// beyond formalism for rare earth		    
                               (*jjjpars[i]).eigenstates(heff,T); //calculate some eigenstates
                               (*jjjpars[i]).mcalc_parameter_storage_init(heff,T);// initialise parameter storage for mcalc
                                     // do some consistency checks
                               ComplexMatrix mcalcpars((*jjjpars[i]).mcalc_parstorage.Rlo(),(*jjjpars[i]).mcalc_parstorage.Rhi(),(*jjjpars[i]).mcalc_parstorage.Clo(),(*jjjpars[i]).mcalc_parstorage.Chi());
                                             mcalcpars=(*jjjpars[i]).mcalc_parstorage;
                               Vector moment(1,j);(*jjjpars[i]).mcalc(moment,T,heff,lnZ,U,mcalcpars);
                               if (fabs((*jjjpars[i]).mom(1)-(*jjjpars[i]).gJ*moment(1))>0.001){fprintf(stderr,"Warning mcdiff: a-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                               if (fabs((*jjjpars[i]).mom(2)-(*jjjpars[i]).gJ*moment(2))>0.001){fprintf(stderr,"Warning mcdiff: b-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                               if (fabs((*jjjpars[i]).mom(3)-(*jjjpars[i]).gJ*moment(3))>0.001){fprintf(stderr,"Warning mcdiff: c-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(1)=(*jjjpars[i]).gJ*moment(1);
                               (*jjjpars[i]).mom(2)=(*jjjpars[i]).gJ*moment(2);
                               (*jjjpars[i]).mom(3)=(*jjjpars[i]).gJ*moment(3);
                            
                     
                               if(Norm((*jjjpars[i]).Zc)==0){fprintf(stderr,"WARNING mcdiff: Z(K) coefficients not found or zero in file %s\n",cffilename);}
                               }
                               if(Norm((*jjjpars[i]).magFFj4)==0){fprintf(stderr,"WARNING mcdiff: <j4(Q)> coefficients not found or zero in file %s\n",cffilename);}
                               if(Norm((*jjjpars[i]).magFFj6)==0){fprintf(stderr,"WARNING mcdiff: <j6(Q)> coefficients not found or zero in file %s\n",cffilename);}

fprintf(fout,"                  corresponding effective fields gjmbHeff [meV]-->");
for(k=1;k<=j;++k){fprintf(fout," %+8.5f",heff(k));}
fprintf(fout,"\n");
 			      }
                             if((*jjjpars[i]).SLR==0){fprintf(stderr,"WARNING mcdiff: SCATTERINGLENGTHREAL not found or zero in file %s\n",cffilename);}
//                             if((*jjjpars[i]).gJ==0){fprintf(stderr,"WARNING mcdiff: GJ not found or zero in file %s - gJ=0 means Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc !\n",cffilename);}
                             if(Norm((*jjjpars[i]).magFFj0)==0){fprintf(stderr,"WARNING mcdiff: <j0(Q)> coefficients not found or zero in file %s\n",cffilename);}
                             if(Norm((*jjjpars[i]).magFFj2)==0){fprintf(stderr,"WARNING mcdiff: <j2(Q)> coefficients not found or zero in file %s\n",cffilename);}
                           }
  fclose(fin_coq);
  fclose(fout);



// print spinconfiguration to mcdiff.sps  (useful for viewing)
print_sps("./results/mcdiff.sps",natmagnetic,a,b,c,alpha,beta,gamma,nr1,nr2,nr3,r1s,r2s,r3s,jjjpars,T,H);
print_mf("./results/mcdiff.sps",mfields,natmagnetic,a,b,c,alpha,beta,gamma,nr1,nr2,nr3,r1s,r2s,r3s,jjjpars,T,H);

printf ("calculating ...\n");  

//now insert also nonmagnetic elements into the unit cell
int ncryst,na,nb,nc;
ncryst = natmagnetic;
for(na = 1;na<=nr1;++na){
 for(nb = 1;nb<=nr2;++nb){
  for(nc = 1;nc<=nr3;++nc){
   if(nat!=0){
    for(i=1;i<=nat;++i){
      ++ncryst;
      J[ncryst]=1;
      jjjpars[ncryst]=new jjjpar((na + x1[i] - 1) / nr1,(nb + y1[i] - 1) / nr2,(nc + z1[i] - 1) / nr3,sl1r[i],sl1i[i],dwf1[i]);
      (*jjjpars[ncryst]).mom=0;
      (*jjjpars[ncryst]).gJ=0;
      }
    }
}}}
delete []x1;delete []y1;delete []z1;delete []da;delete []db;delete[]dc;
delete []sl1r;delete[]sl1i;delete[]dwf1;

int m=0;
Vector * hkl = new Vector[MAXNOFREFLECTIONS+1];for(i=0;i<=MAXNOFREFLECTIONS;++i){hkl[i]=Vector(1,3);}
float D[MAXNOFREFLECTIONS+1];
float theta[MAXNOFREFLECTIONS+1];
float intmag[MAXNOFREFLECTIONS+1];
float intmagdip[MAXNOFREFLECTIONS+1];
float ikern[MAXNOFREFLECTIONS+1];
float out10[MAXNOFREFLECTIONS+1];
float out11[MAXNOFREFLECTIONS+1];
complex <double> mx[MAXNOFREFLECTIONS+1];
complex <double> my[MAXNOFREFLECTIONS+1];
complex <double> mz[MAXNOFREFLECTIONS+1];
complex <double> mxmy[MAXNOFREFLECTIONS+1];
complex <double> mxmz[MAXNOFREFLECTIONS+1];
complex <double> mymz[MAXNOFREFLECTIONS+1];
complex <double> mx2[MAXNOFREFLECTIONS+1];
complex <double> my2[MAXNOFREFLECTIONS+1];
complex <double> mz2[MAXNOFREFLECTIONS+1];

rezcalc (r1, r2, r3, rez1, rez2, rez3);
Vector hhkkll(1,3);
int code=0;
// if hkllist is given, read the file and put hkls to hkl[i], m is number of reflections to be considered
if (argc>1){int nr;
      float nn[20];nn[0]=19;
     // open hkllist file
     fprintf(stdout,"reading hkl list from file %s\n",argv[1]);
       fin = fopen_errchk (argv[1], "rb");
       while(feof(fin)==false){nr=inputline(fin,nn);
                               if(nr>2)
                               {hhkkll(1)=nn[1];hhkkll(2)=nn[2];hhkkll(3)=nn[3];++m;
                                code=1;                                
                               // transformieren der millerindizes auf magnetische einheitszelle
                                  hkl[m](1)=hhkkll*(rtoxyz.Inverse()*r1);
                                  hkl[m](2)=hhkkll*(rtoxyz.Inverse()*r2);
                                  hkl[m](3)=hhkkll*(rtoxyz.Inverse()*r3);
                              // check if magnetic reflection is indeed on magnetic reciprocal lattice
                              if(fabs(rint(hkl[m](1))-hkl[m](1))>SMALL||fabs(rint(hkl[m](2))-hkl[m](2))>SMALL||fabs(rint(hkl[m](3))-hkl[m](3))>SMALL)
                                {fprintf(stderr,"Warning - reading (%g %g %g): calculation impossible, because this corresponds to ", hhkkll(1),hhkkll(2),hhkkll(3));
                                 fprintf(stderr,"non integer magnetic reciprocal lattice (%g %g %g)\n\n", hkl[m](1), hkl[m](2), hkl[m](3));
                                }
                                if(nr>3){mx[m]=complex <double> (nn[4],0);code=2;}// intensities given                                
                                if(nr>4){my[m]=complex <double> (nn[5],0);code=3;}// errors given                                
                              }}
       fclose(fin);      
           }

// transformieren der millerindizes auf kristallographische einheitszelle


neutint(jjjpars,code,T,lambda, thetamax, ovalltemp, lorenz, r1, r2, r3, n,  J, m, hkl, D, theta, intmag,intmagdip, ikern, out10, out11,mx,my,mz,mxmy,mxmz,mymz,mx2,my2,mz2,colcode,Pxyz);



// transformieren der millerindizes auf kristallographische einheitszelle

for(i=1;i<=m;++i){hhkkll=hkl[i];
                  hkl[i]=hhkkll(1)*rez1+hhkkll(2)*rez2+hhkkll(3)*rez3;
                  hkl[i]/=2.0*PI;
                  hhkkll=hkl[i];
                  hkl[i](1)=hhkkll*rtoxyz.Column(1);
                  hkl[i](2)=hhkkll*rtoxyz.Column(2);
                  hkl[i](3)=hhkkll*rtoxyz.Column(3);
                 }


printeln(jjjpars,code,"./results/mcdiff.out","mcdiff.in", unitcellstr,T,H, lambda, ovalltemp, lorenz, r1, 
          r2, r3, n,  J, m, hkl, ikern, intmag,intmagdip, D, theta, out10, out11,mx,my,mz,mxmy,mxmz,mymz,
          mx2,my2,mz2,a,b,c,colcode,P);

fprintf (stderr,"...results written to ./results/mcdiff.out\n");
fprintf (stderr,"***********************************************************\n");
fprintf (stderr,"                   End of Program mcdiff\n");
fprintf (stderr,"reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
fprintf (stderr,"***********************************************************\n");

//  for (i=1;i<=n;++i){delete jjjpars[i];}
//  delete []jjjpars;
  delete []hkl;
  delete []J;
 return 0;
}


