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
{ std::clock_t startcputime = std::clock();
  FILE * fin, * fin_coq, * fout;
  float ovalltemp,thetamax,lambda,a=0,b=0,c=0,alpha=0,beta=0,gamma=0;
  double T=0;
  int i,j,k,n,lorenz,nat, nofatoms,nr1=0,nr2=0,nr3=0,natmagnetic;
  long int pos=0;
  int use_dadbdc=0;
  char instr[MAXNOFCHARINLINE+1];
  char somestring[MAXNOFCHARINLINE+1];
  char sipffilename[MAXNOFCHARINLINE+1];
  char unitcellstr[MAXNOFCHARINLINE+1];
  float numbers[70];numbers[0]=70;
  float numbers1[70];numbers1[0]=70;
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

  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&strstr (instr, "%SECTION 2%")==NULL) // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
    if (pos==-1) 
       {fprintf(stderr,"Error mcdiff: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   extract(instr,"nofatoms",nofatoms);  
   extract(instr,"lambda", lambda);
   extract(instr, "thetamax", thetamax);
   extract(instr, "nat", nat);
   extract(instr, "natcryst", nat);
   extract(instr, "ovalltemp", ovalltemp);
   extract(instr, "lorentz", lorenz);
   for(int i=1;i<=usrdefoutcols[0];++i) // extract user defined output columns
   {sprintf(somestring,"out%i",usrdefoutcols[i]);
    extract(instr, somestring,colcode[usrdefoutcols[i]]);
   }
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
printf("                 output:");
for(int i=1;i<=usrdefoutcols[0];++i)printf("column %i=%s",usrdefoutcols[i],colheader[colcode[usrdefoutcols[i]]]);
printf("\n");
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
fprintf(fout,"#     out* controls the type of output in user defined column * of mcdiff.out (optional)\n");
for(int i=1;i<=usrdefoutcols[0];++i)
fprintf(fout,"#! out%i=%i \n",usrdefoutcols[i],colcode[usrdefoutcols[i]]);
fprintf(fout,"#     ... in out*=n the numbers n have the following meaning:\n");
for(i=0;i<=COLHEADERDIM;++i){
fprintf(fout,"#            %i....%s\n",i,colheader[i]);
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

// input section 2 *********************************************************

  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&a==0&&b==0&&c==0) // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
    if (pos==-1) 
       {fprintf(stderr,"Error mcdiff: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
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
  fseek(fin_coq,pos,SEEK_SET); 
  printf ("     section 2 - natcryst=%i\n",nat);
  float *x1;x1=new float[nat+1];float*y1;y1=new float[nat+1];float*z1;z1=new float[nat+1];
  float *da;da=new float[nat+1];float*db;db=new float[nat+1];float*dc;dc=new float[nat+1];
  float *sl1r;sl1r=new float[nat+1];float*sl1i;sl1i=new float[nat+1];float *dwf1;dwf1=new float[nat+1];

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
Matrix rtoijk(1,3,1,3); // define transformation matrix to calculate components of
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

// --------------------- old internal coordinates - removed and changed to the above 30.10.2011 MR
/* rtoijk(1,1)=0;
rtoijk(2,1)=a*sin(gamma*PI/180);
rtoijk(3,1)=a*cos(gamma*PI/180);

rtoijk(1,2)=0;
rtoijk(2,2)=0;
rtoijk(3,2)=b;

rtoijk(3,3)=c*cos(alpha*PI/180);
rtoijk(2,3)=c*(cos(beta*PI/180)-cos(alpha*PI/180)*cos(gamma*PI/180))/sin(gamma*PI/180);
if (fabs(rtoijk(2,3))>c){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoijk(1,3)=c*c-rtoijk(2,3)*rtoijk(2,3)-rtoijk(3,3)*rtoijk(3,3);
if (rtoijk(1,3)<=0){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoijk(1,3)=sqrt(rtoijk(1,3));
*/
r1=(rtoijk*r1)*(double)nr1;
r2=(rtoijk*r2)*(double)nr2;
r3=(rtoijk*r3)*(double)nr3;

// transform also Projection vector
Vector Pxyz (1,3);
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

jjjpar ** jjjpars = new jjjpar * [n+1];
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
fprintf(fout,"#  <Sa>  <La> <Sb> <Lb >  <Sc> <Lc>  (optional) denote the spin and orbital angular momentum components \n");
fprintf(fout,"# 'Hxc1' 'Hxc2' 'Hxc3' (optional line, used to go beyond dipole approx for formfactor)\n");
fprintf(fout,"#                                     denote the corresponding exchange fields in meV\n");
fprintf(fout,"#\n");
fprintf(fout,"#{atom-file} da[a]  db[b]    dc[c]     dr1[r1]  dr2[r2]  dr3[r3]   <Ma>     <Mb>     <Mc> [mb] [optional <Sa> <La> <Sb> <Lb> <Sc> <Lc> ]\n");
fprintf(fout,"#{corresponding exchange fields Hxc [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)}\n");

mfcf mfields(1,1,1,natmagnetic,MAX_NOF_MF_COMPONENTS); // MAX_NOF_MF_COMPONENTS is maximum of nofmfcomponents - we take it here !
mfields.clear();
int maxmfcomponents=3;// for printout of mcdiff.mf we check what is the largest nofmfcomponents in mcdiff.in
for(i=1;i<=natmagnetic;++i){
                            instr[0]='#';
                            while(instr[strspn(instr," \t")]=='#'){pos=ftell(fin_coq);
                                                                   if(feof(fin_coq)==1){fprintf(stderr,"mcdiff Error: end of file before all magnetic atoms could be read\n");exit(EXIT_FAILURE);}
                                                                  fgets(instr,MAXNOFCHARINLINE,fin_coq);
                                                                  }
			     // get sipffilename out of "{filename}   ..."

                            if(instr[strspn(instr," \t")]!='{'){fprintf(stderr,"ERROR mcdiff: magnetic atom line has to start with '{'\n");exit (EXIT_FAILURE);}
			    if (strchr(instr,'}')==NULL){fprintf(stderr,"ERROR mcdiff: no '}' found after filename for magnetic atom %s\n",instr);exit (EXIT_FAILURE);}

                            instr[strspn(instr," \t")]='=';
			    extract(instr,"",sipffilename,(size_t)MAXNOFCHARINLINE);
                            if(strchr(sipffilename,'}')!=NULL){*strchr(sipffilename,'}')='\0';}
                            if(strchr(sipffilename,' ')!=NULL){*strchr(sipffilename,' ')='\0';}
                            if(strchr(sipffilename,'\t')!=NULL){*strchr(sipffilename,'\t')='\0';}
                            //printf("%s\n",sipffilename);

                             // read the rest of the line and split into numbers
                            fseek(fin_coq,pos+strchr(instr,'}')-instr+1,SEEK_SET); 
                            j=inputline(fin_coq,numbers);
if(use_dadbdc!=0)        {       numbers[4]= (numbers[1]*rez1(1)+numbers[2]*rez1(2)+numbers[3]*rez1(3))/2/PI;
                                 numbers[5]= (numbers[1]*rez2(1)+numbers[2]*rez2(2)+numbers[3]*rez2(3))/2/PI;
                                 numbers[6]= (numbers[1]*rez3(1)+numbers[2]*rez3(2)+numbers[3]*rez3(3))/2/PI;
                         }
                            if (j<9) {fprintf(stderr,"ERROR mcdiff: too few parameters for magnetic atom %i: %s\n",i,instr);exit(EXIT_FAILURE);}
                             // determine jxc .... dimension of exchange field if present >>>>>>>>>>>>>>>>
                            long int currentpos=ftell(fin_coq);instr[0]='#';int jxc;
                            while(instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0){pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
                            if (strchr(instr,'>')==NULL||instr[strspn(instr," \t")]=='#')
                             {jxc=1;} // no ">" found --> do dipole approx
                             else          
                             {fseek(fin_coq,pos+strchr(instr,'>')-instr+1,SEEK_SET); 
                              jxc=inputline(fin_coq,numbers1);printf("dimension of mf = %i\n",jxc);
                             }
                             fseek(fin_coq,currentpos,SEEK_SET); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
fprintf(fout,"{%s} %8.5f %8.5f %8.5f  ",sipffilename,numbers[4]*r1s(1)+numbers[5]*r2s(1)+numbers[6]*r3s(1),numbers[4]*r1s(2)+numbers[5]*r2s(2)+numbers[6]*r3s(2),numbers[4]*r1s(3)+numbers[5]*r2s(3)+numbers[6]*r3s(3));
fprintf(fout,"%8.5f %8.5f %8.5f  ",numbers[4],numbers[5],numbers[6]); // positions
fprintf(fout," %+8.5f %+8.5f %+8.5f ",(*jjjpars[i]).mom(1),(*jjjpars[i]).mom(2),(*jjjpars[i]).mom(3)); // magmoments
for(k=10;k<=j;++k){fprintf(fout," %+8.5f",numbers[k]);}
fprintf(fout,"\n");
                            instr[0]='#';
                            while(instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0){pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
                            if (strchr(instr,'>')==NULL||instr[strspn(instr," \t")]=='#')
                             {fseek(fin_coq,pos,SEEK_SET);} // no ">" found --> do dipole approx
                             else          
                             {Vector Qvec(1,3);Qvec=0;ComplexVector Mq(1,3);
                              fseek(fin_coq,pos+strchr(instr,'>')-instr+1,SEEK_SET); 
                              j=inputline(fin_coq,numbers);printf("dimension of mf = %i\n",j);
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

fprintf(fout,"                    corresponding exchange fields gjmbHxc [meV]-->");
for(k=1;k<=j;++k){fprintf(fout," %+8.5f",gjmbHxc(k));}
fprintf(fout,"\n");
 			      }}
                             if((*jjjpars[i]).SLR==0){fprintf(stderr,"WARNING mcdiff: SCATTERINGLENGTHREAL not found or zero in file %s\n",sipffilename);}
//                             if((*jjjpars[i]).gJ==0){fprintf(stderr,"WARNING mcdiff: GJ not found or zero in file %s - gJ=0 means Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc !\n",sipffilename);}
                             (*jjjpars[i]).checkFFcoeffnonzero(0);
                             (*jjjpars[i]).checkFFcoeffnonzero(2);
                           }
  fclose(fin_coq);
  fclose(fout);

mfields.resetnofc(maxmfcomponents);

// print spinconfiguration to mcdiff.sps  (useful for viewing)
print_sps(natmagnetic,a,b,c,alpha,beta,gamma,nr1,nr2,nr3,r1s,r2s,r3s,jjjpars,T,H);
print_mf(mfields,natmagnetic,a,b,c,alpha,beta,gamma,nr1,nr2,nr3,r1s,r2s,r3s,jjjpars,T,H);

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
      jjjpars[ncryst]=new jjjpar((na + x1[i] - 1) / nr1,(nb + y1[i] - 1) / nr2,(nc + z1[i] - 1) / nr3,sl1r[i],sl1i[i],dwf1[i]);
      (*jjjpars[ncryst]).FF_type=+1;//J[ncryst]=1;
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
float intmag[MAXNOFREFLECTIONS+1];
float intmagdip[MAXNOFREFLECTIONS+1];
float ikern[MAXNOFREFLECTIONS+1];
float * out[NOFOUTPUTCOLUMNS+1];
for(int i=1;i<=usrdefoutcols[0];++i)
 {out[usrdefoutcols[i]]=new float [MAXNOFREFLECTIONS+1];
  if(out[usrdefoutcols[i]]==NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
 }

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
                                  hkl[m](1)=hhkkll*(rtoijk.Inverse()*r1);
                                  hkl[m](2)=hhkkll*(rtoijk.Inverse()*r2);
                                  hkl[m](3)=hhkkll*(rtoijk.Inverse()*r3);
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
printheader(jjjpars,code,"./results/mcdiff.out","mcdiff.in", unitcellstr,T,H, lambda, ovalltemp, lorenz, r1, 
          r2, r3, n,  m,a,b,c,P,Pxyz);


neutint(jjjpars,code,T,lambda, thetamax, ovalltemp, lorenz, r1, r2, r3, n,  m, hkl, D, intmag,intmagdip, ikern, out,mx,my,mz,mxmy,mxmz,mymz,mx2,my2,mz2,Pxyz);



// transformieren der millerindizes auf kristallographische einheitszelle

for(i=1;i<=m;++i){hhkkll=hkl[i];
                  hkl[i]=hhkkll(1)*rez1+hhkkll(2)*rez2+hhkkll(3)*rez3;
                  hkl[i]/=2.0*PI;
                  hhkkll=hkl[i];
                  hkl[i](1)=hhkkll*rtoijk.Column(1);
                  hkl[i](2)=hhkkll*rtoijk.Column(2);
                  hkl[i](3)=hhkkll*rtoijk.Column(3);
                 }


printreflist(jjjpars,code,"./results/mcdiff.out","mcdiff.in", unitcellstr,T,H, lambda, ovalltemp, lorenz, r1, 
          r2, r3, n,  m, hkl, ikern, intmag,intmagdip, D, out,mx,my,mz,mxmy,mxmz,mymz,
          mx2,my2,mz2,a,b,c,P,Pxyz);

double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
std::cout << "#! nofhkls=" << m << " different q vectors generated " << std::endl;

fprintf (stderr,"...results written to ./results/mcdiff.out\n");
fprintf (stderr,"***********************************************************\n");
fprintf (stderr,"                   End of Program mcdiff\n");
fprintf (stderr,"reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
fprintf (stderr,"***********************************************************\n");

//  for (i=1;i<=n;++i){delete jjjpars[i];}
//  delete []jjjpars;
  delete []hkl;
  for(i=1;i<=usrdefoutcols[0];++i)delete out[usrdefoutcols[i]];  
 return 0;
}


