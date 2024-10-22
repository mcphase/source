// class to store a crystal structure


#ifndef CRYSTSTRUC
#define CRYSTSTRUC

#define MAXNOFATOMS 5000

class cryststruct
{ public:
int nofatoms,nofcomponents,maxnofatoms;
   float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS];
  char * sipffilenames[MAXNOFATOMS];
   Matrix r;
   Vector abc;

cryststruct()
{nofatoms=0;nofcomponents=3,maxnofatoms=MAXNOFATOMS;
 r=Matrix(1,3,1,3);
 abc=Vector(1,6);
 abc(4)=90;
 abc(5)=90;
 abc(6)=90;

}

double a(){return abc(1);}
double b(){return abc(2);}
double c(){return abc(3);}
double alpha(){return abc(4);}
double beta(){return abc(5);}
double gamma(){return abc(6);}


void print_mcdiff_in_header(FILE * fout,const char * program,int natcryst)
{
fprintf(fout,"\
#\n\
#<!--mcdiff.mcdiff.in>\n\
#***************************************************************\n\
#    Input file for mcdiff created by program %s \n\
#      mcdiff is a program for the calculation of elastic\n\
#   neutron diffraction and resonant magnetic Xray scattering\n\
#  reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n\
#***************************************************************\n\
# this input file contains 4 sections corresponding to different\n\
# groups of parameters\n\
#\n\
# - all lines have to start with a # sign with the  exception of\n\
#   the lines containing atomic positional parameters\n\
# - the other parameters have to be defined in the corresponding\n\
#   section by statements such as parameter=value\n\
# - the sequence of the parameters within a section is arbitrary\n\
#\n\
#\n\
# %%SECTION 1%%  OVERALL PARAMETERS\n\
#\n\
#! lambda   = 2.3587  wavelength (A)\n\
#\n\
#! thetamax = 10   maximum bragg angle (deg)\n\
#\n\
#! ovalltemp= 0  overall temperature factor (A^2)\n\
#           ...I ~ EXP(-2 * ovalltemp * sintheta^2 / lambda^2)\n\
#                  relation to other notations:\n\
#                  ovalltemp = Biso = 8 pi^2 Uiso^2\n\
#\n\
#! lorentz=0  type of lorentzfactor to be used\n\
#            0.....no lorentzfactor \n\
#            1.....neutron powder flat sample\n\
#            2.....neutron powder cylindrical sample\n\
#            3.....neutron single crystal\n\
#            4.....neutron TOF powder cyl. sample - d-pattern log scaled\n\
#            5.....neutron TOF powder cyl. sample - d-pattern normal scaled\n\
#\n\
# out*  controls the type of output in  mcdiff.out \n\
#!out0=1  0..short header, 1...standard header, 2... long header in mcdiff.out\n\
# (a negative value triggers calculation of magnetic xray scattering)\n\
#\n\
#!  nofoutputcolumns=12  number of columns in output file mcdiff.out\n\
# choose out* to set type of desired output in column 1 to 12\n\
# (default is h k l d Q 2theta Inuc Imag Itot)\n\
#\n\
#!out1=31 out2=32 out3=33 out4=37 out5=38 out6=39 out7=27 out8=28 out9=29 out10=1 out11=0 out12=30 \n\
#\n\
#     ... in out*=n the numbers n have the following meaning:\n\
#            0....LF #\n\
#            1....|NSF|[b] #\n\
#            2....Re(NSF)[b] #\n\
#            3....Im(NSF)[b] #\n\
#            4....|MSF| #\n\
#            5....|MSF.P| #\n\
#            6....Re(MSF.P) #\n\
#            7....Im(MSF.P) #\n\
#            8....|MSFdip| #\n\
#            9....|MSFdip.P| #\n\
#            10....Re(MSFdip.P) #\n\
#            11....Im(MSFdip.P) #\n\
#            12....angl(Q,P)[°] #\n\
#            13....i(MSFxMSF*).P #\n\
#            14....I+ #\n\
#            15....I- #\n\
#            16....I+/I- #\n\
#            17....i(MSFxMSF*)dip.P #\n\
#            18....Idip+ #\n\
#            19....Idip- #\n\
#            20....Idip+/Idip- #\n\
#            21....2*|MSF.P|/sin^2(angl(Q,P) #\n\
#            22....2*|MSFdip.P|/sin^2(angl(Q,P) #\n\
#            23....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-) #\n\
#            24....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-) #\n\
#            25....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-) #\n\
#            26....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-) #\n\
#            27....Inuc(2t)	#\n\
#            28....Imag(2t)   #\n\
#            29....Itot(2t)   #\n\
#            30....Imag_dip(2t) #\n\
#            31....h     #\n\
#            32.... k     #\n\
#            33.... l     #\n\
#            34....d[A]       #\n\
#            35....|Q|[1/A]   #\n\
#            36....2theta     #\n\
#            37....Qi[1/A]    #\n\
#            38....Qj[1/A]     Qi Qj Qk are euclidean components of scattering vector#\n\
#            39....Qk[1/A]     with j||b, k||(a x b) and i normal to k and j#\n\
#            40....T[K]       #\n\
#            41....Ha[T]      #\n\
#            42....Hb[T]      #\n\
#            43....Hc[T]      #\n\
#            44....hprim      #\n\
#            45....kprim      #\n\
#            46....lprim      #\n\
#            47....Itotdip(2t) #\n\
#\n\
#           In the above the intensities I+ and I- are the intensities in a polarised neutron\n\
#           experiment with incident polarisation up (+) or down (-):\n\
#            I+-=LF exp(-OTF Q^2/8pi^2) \n\
#                    [ |NSF/NB|^2 + 3.65/4pi (|MSF|^2-+i(MSF x MSF*).P)/NB^2 \n\
#                        +-  sqrt(3.65/4pi)/NB^2 (NSF (MSF*.P) + NSF* (MSF.P)]\n\
#           LF  ..... Lorentzfactor\n\
#           MSF ..... magnetic structure factor\n\
#           NSF ..... nuclear structure factor\n\
#\n\
#\n\
#             For some of the above options we need the\n\
#! Pa=  0.0000   Components of Polarisation Vector in terms of lattice vectors P=(Pa * a + Pb * b + Pc *c)\n\
#! Pb=  0.0000   Note: the length of P, i.e. |P| indicates the degree of beam polarisation (|P|<=1)\n\
#! Pc=  1.0000\n\
#\n\
#\n\
#\n\
# %%SECTION 2%% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL\n\
#\n",program);
if(natcryst==0)fprintf(fout,"# !!!  NONMAGNETIC ATOMS HAVE NOT BEEN GENERATED by prgram %s!!!\n",program);

fprintf(fout,"\
#\n\
#! natcryst=%i      number of nonmagnetic atoms in primitive crystalographic unit cell\n\
#\n\
# it follows a list of natcryst lines with nonmagnetic atoms\n\
# ... notes: - if an occupancy other than 1.0 is needed, just reduce\n\
#              the scattering length linear accordingly\n\
#            - Debye Waller Factor notation: sqr(Intensity) ~ structure factor ~\n\
#              ~sum_n ()n exp(-2 DWFn sin^2(theta) / lambda^2)=EXP (-Wn),\n\
#              relation to other notations: 2*DWF = B = 8 pi^2 <u^2>, units DWF (A^2)\n\
#\n\
#! use_dadbdc=1\n\
#            - 0 means: da db and dc are not used by the program (unless you enter a line #! use_dadbdc=1),\n\
#               dr1,dr2 and dr3 refer to the primitive lattice given below\n\
#\n\
# Real Imag[scattering length(10^-12cm)]   da(a)    db(b)    dc(c)    dr1(r1)  dr2(r2)  dr3(r3)  DWF(A^2)\n",natcryst);


}



};

#endif
