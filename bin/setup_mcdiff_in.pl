#!/usr/bin/perl
use Getopt::Long;

GetOptions("help"=>\$helpflag);
usage() if $helpflag||$#ARGV<2;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$ARGV[0]=eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;$ARGV[1]=eval $ARGV[1];
$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;$ARGV[2]=eval $ARGV[2];
$ARGV[3]=~s/exp/essp/g;$ARGV[3]=~s/x/*/g;$ARGV[3]=~s/essp/exp/g;$ARGV[3]=eval $ARGV[3];

print STDOUT << "EOF";
*******************************************************
setting up mcdisp.mf to be used by mcdisp
T=$ARGV[0] K Ha=$ARGV[1] T Hb=$ARGV[2] T Hc=$ARGV[3] T
*******************************************************
reading results/mcphas.mf
.... trying to calculate results/spins.*
EOF

print STDOUT << "EOF";
reading results/mcphas.mf
....writing results/spins.*
EOF

system ("spins $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]");

print STDOUT << "EOF";
generating mcdiff.in ...
EOF

open (Fout,">mcdiff.in");
print Fout << "EOF";
#
#<!--mcdiff.mcdiff.in>
#***************************************************************
#      mcdiff is a program for the calculation of elastic
#   neutron diffraction and resonant magnetic Xray scattering
#  reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R
#***************************************************************
# this input file contains 4 sections corresponding to different
# groups of parameters
#
# - all lines have to start with a # sign with the  exception of
#   the lines containing atomic positional parameters
# - the other parameters have to be defined in the corresponding
#   section by statements such as parameter=value
# - the sequence of the parameters within a section is arbitrary
#
#
# %SECTION 1%  OVERALL PARAMETERS
#
#! lambda   = 2.3587  wavelength (A)
#
#! thetamax = 10   maximum bragg angle (deg)
#
#! ovalltemp= 0  overall temperature factor (A^2)
#           ...I ~ EXP(-2 * ovalltemp * sintheta^2 / lambda^2)
#                  relation to other notations:
#                  ovalltemp = Biso = 8 pi^2 Uiso^2
#
#! lorentz=0  type of lorentzfactor to be used
#            0.....no lorentzfactor 
#            1.....neutron powder flat sample
#            2.....neutron powder cylindrical sample
#            3.....neutron single crystal
#            4.....neutron TOF powder cyl. sample - d-pattern log scaled
#            5.....neutron TOF powder cyl. sample - d-pattern normal scaled
#
#     out* controls the type of output in user defined column * of mcdiff.out
#! out4=30    (optional)
#! out5=31    (optional)
#! out6=32    (optional)
#! out10=1    (optional)
#! out11=0    (optional)
#     ... in out*=n the numbers n have the following meaning:
#            0....LF          
#            1....|NSF|[b]    
#            2....Re(NSF)[b]  
#            3....Im(NSF)[b]  
#            4....|MSF|       
#            5....|MSF.P|     
#            6....Re(MSF.P)   
#            7....Im(MSF.P)   
#            8....|MSFdip|    
#            9....|MSFdip.P|  
#            10....Re(MSFdip.P)
#            11....Im(MSFdip.P)
#            12....angl(Q,P)[°]
#            13....i(MSFxMSF*).P
#            14....I+          
#            15....I-          
#            16....I+/I-       
#            17....i(MSFxMSF*)dip.P
#            18....Idip+       
#            19....Idip-       
#            20....Idip+/Idip- 
#            21....2*|MSF.P|/sin^2(angl(Q,P)
#            22....2*|MSFdip.P|/sin^2(angl(Q,P)
#            23....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
#            24....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
#            25....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
#            26....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
#            27....Qx[1/A]     
#            28....Qy[1/A]     
#            29....Qz[1/A]     
#            30....d[A]        
#            31....|Q|[1/A]    
#            32....2theta      
#
#           In the above the intensities I+ and I- are the intensities in apolarised neutron experiment
#           with incident polarisation up (+) or down (-):
#            I+-=LF exp(-OTF Q^2/8pi^2) 
#                    [ |NSF/NB|^2 + 3.65/4pi (|MSF|^2-+i(MSF x MSF*).P)/NB^2 
#                        +-  sqrt(3.65/4pi)/NB^2 (NSF (MSF*.P) + NSF* (MSF.P)]
#           LF  ..... Lorentzfactor
#           MSF ..... magnetic structure factor
#           NSF ..... nuclear structure factor
#
#
#             For some of the above options we need the
#! Pa=  0.0000   Components of Polarisation Vector in terms of lattice vectors P=(Pa * a + Pb * b + Pc *c)
#! Pb=  0.0000   Note: the length of P, i.e. |P| indicates the degree of beam polarisation (|P|<=1)
#! Pc=  1.0000
#
#
#
# %SECTION 2% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL
#
# !!!  NONMAGNETIC ATOMS HAVE NOT BEEN GENERATED BY setup_mcdiff_in !!!
#
#! natcryst=0      number of nonmagnetic atoms in primitive crystalographic unit cell
#
# it follows a list of natcryst lines with nonmagnetic atoms
# ... notes: - if an occupancy other than 1.0 is needed, just reduce
#              the scattering length linear accordingly
#            - Debye Waller Factor notation: sqr(Intensity) ~ structure factor ~
#              ~sum_n ()n exp(-2 DWFn sin^2(theta) / lambda^2)=EXP (-Wn),
#              relation to other notations: 2*DWF = B = 8 pi^2 <u^2>, units DWF (A^2)
#
#! use_dadbdc=0
#            - 0 means: da db and dc are not used by the program (unless you enter a line #! use_dadbdc=1),
#               dr1,dr2 and dr3 refer to the primitive lattice given below
#
# Real Imag[scattering length(10^-12cm)]   da(a)    db(b)    dc(c)    dr1(r1)  dr2(r2)  dr3(r3)  DWF(A^2)
#   0.73250   0.00000                       -1.43200 -1.43200  0.71600  0.00000 -0.71600  1.43200  0.00000
#   0.73250   0.00000                       -2.56800 -2.56800  1.28400  0.00000 -1.28400  2.56800  0.00000
#
#
# %SECTION 3% DESCRIPTION OF THE LATTICE
#
#
# Note: what follows here may directly be taken from the output of program spins
#       (file spins.out) or charges (file charges.out)
# -----------------------------------------------------------------------------
EOF

{print "\n\nSetting up mcdiff.in ...\n";
open (Fin,"results/spins.out");
while(<Fin>) {print Fout $_;}
close Fout,Fin;
}

print STDOUT << "EOF";

    mcdiff.in generated: you can start now mcdiff
    However, please remember to set wavelength, lorentzfactor etc.
    in input file mcdiff.in

    You can view the magnetic structure in postscriptfiles
    results/spins*.eps, by fp_studio results/spins.fst and
    by javaview results/spins.jvx

EOF
exit;

sub usage() {

  print STDERR << "EOF";

    setup_mcdiff_in: program to setup mcdiff_in with information on spinconfiguration
                    to be used by program mcdiff. Note, you must
                    have done a mcphas calculation to stabilise
                    a magnetic structure at the desired Temperature/Field.
                    setup_mcdiff_in reads the results of this calculation
                    from results/mcphas.sps and results/mcphas.mf generates an
                    input file mcdiff.in

    usage: setup_mcdiff_in T Ha Hb Hc

     -h          : this (help) message
      T          : Temperature (K)
      Ha,Hb,Hc   : Magnetic Field (T)

    required input files:

    results/mcphas.sps
                 :  result of a mcphas calculation

    output files:

    mcdiff.in    :  required input file for mcdiff

    - after running this program you can start mcdiff to do the calculation
      magnetic diffraction pattern
EOF

  exit;

}
# **********************************************************************************************
# extracts variable from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# ($standarddeviation)=extract("sta","data.dat");
#
# ... it stores 0.24 in the variable $standarddeviation
#
sub extract {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";$value="";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);}}
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
