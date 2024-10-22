#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use File::Copy;
#function make_input_mcphas
#% Lucian G Pascut
#% pascutlucian@yahoo.com

unless ($#ARGV >=0)
{print STDERR << "EOF";

    $0:

    program $0  used to create mcphas.j type file from output of powdercell

    usage: $0 *.*  [*.cel]

    *.* .. filenname of powdercell file with atomic positions
    *.cel ... optional, filename of powdercell .CEL file

    output is written to mcphas.j, mcdiff.in

 In case of any error please check the format of the input file *.*

 Example:
  No   name       crystal coordinates          cartesian coordinates
                 x        y        z           x        y        z
  ------------------------------------------------------------------
  1     Sr1    0.3644   0.0000   0.2500     1.0962  -4.1497  -2.7991
  ...

EOF
 exit 0;
}else{print STDERR "#* $0 *\n";}

  if($ARGV[1])
  {# read .cel file if given
     $file=$ARGV[1];
    unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
    while($line=<Fin>)
     {
      if ($line=~/^\s*CELL/) {@c=split(" ",$line);
                          $alattice=$c[1];
                          $blattice=$c[2];
                          $clattice=$c[3];
                          $alpha=$c[4];
                          $beta=$c[5];
                          $gamma=$c[6];
                   print "\nLATTICE from $file\n";
                   print " a=$alattice A b=$blattice A c=$clattice A\n alpha=$alpha deg beta=$beta deg gamma=$gamma deg\n\n";
                              }
     }
  }
  else
  {
  print "Give the lattice constants a(A) and then press enter\n";
  $alattice= <STDIN>;$alattice=~s/\n//;
  print "Give the lattice constants b(A) and then press enter\n";
  $blattice=<STDIN>;$blattice=~s/\n//;
  print "Give the lattice constants c(A) and then press enter\n";
  $clattice=<STDIN>;$clattice=~s/\n//;
  print "Give the angle alpha(deg) and then press enter\n";
  $alpha=<STDIN>;$alpha=~s/\n//;
  print "Give the angle beta(deg) and then press enter\n";
  $beta=<STDIN>;$beta=~s/\n//;
  print "Give the angle gamma(deg)  and then press enter\n";
  $gamma=<STDIN>;$gamma=~s/\n//;
  }


 scattlengthoutput();
   $p=0;
 
     $file=$ARGV[0];
 #   %
 #   %begin:"find how many different types of atoms do we have in the structure"
 #   %
 #   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
    while($line=<Fin>)
     {
      if ($line=~/^\s*\d/) {@c=split(" ",$line);
                            $atomtype[$p]=$c[1];
                            $x[$p]=$c[2];$y[$p]=$c[3];$z[$p]=$c[4];
                            while($x[$p]<0.0){$x[$p]+=1.0;};while($x[$p]>1.0){$x[$p]-=1.0;};
                            while($y[$p]<0.0){$y[$p]+=1.0;};while($y[$p]>1.0){$y[$p]-=1.0;};
                            while($z[$p]<0.0){$z[$p]+=1.0;};while($z[$p]>1.0){$z[$p]-=1.0;};
                            $p++;
                           }
     }
    $pp=0;
    for($i=0;$i<$p;++$i)
    { $ppp=1;
      foreach(@atomtypetempp)
       {if($_=~$atomtype[$i]){$ppp=0;}
       }
      if($ppp==1){$atomtypetempp[$pp]=$atomtype[$i];$pp+=$ppp;}
    }
    print "\n   There are ".($pp)." different atom types:\n";
#    %find the number of atoms in the same group
    $noofatom=0;$i=0;
      foreach(@atomtypetempp)
       {foreach(@atomtype){
        if($_=~$atomtypetempp[$i])        {++$nofatoms[$i];}
        }
        ++$i;
       }
    print "@atomtypetempp\n";
    print "    Number of different atoms for each type:\n";
    print "@nofatoms\n";
    @nofmagatoms=@nofatoms;

$i=0;
  foreach(@atomtypetempp)
  {  open (Fout1, ">$_.sipf");
     print "For $_ atom the oxidation state is? \nPlease give the charge in units of elementary charge |e|; e.g. -2 or +1 or +3:\n";
     $charge=<STDIN>;$charge=~s/\n//;
     print " ... is it magnetic (1) or nonmagnetic (0) ? "; $magnetic=<STDIN>;$magnetic=~s/\n//;
     if ($magnetic==0){$nofmagatoms[$i]=0;
     print Fout1 << "EOF";
#!MODULE=PHONON
#<!--mcphase.sipf-->
#IONTYPE=$_
CHARGE=$charge
MAGNETIC=$magnetic

MODPAR1=1  #mass in(m0) !powdercell defaults to hydrogen mass=1 please insert correct value !!
MODPAR2=0   # Kxx
MODPAR3=0   # Kyy
MODPAR4=0   # Kzz
MODPAR5=0  # Kxy  in (meV)
MODPAR6=0  # Kxz
MODPAR7=0  # Kyz
MODPAR8=0.4 # umax          maximum (cutoff) for displacement [a0=0.5219 A]
MODPAR9=0   #                0       umax restriction in all directions
            #                1,2,3   umax restriction in x y z direction only
            #                4       umax restriction in x and y direction
            #                5       umax restriction in x and z direction
            #                6       umax restriction in y and z direction#-------------------------------------------------------
# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)
#                      with s=sin(theta)/lambda=Q/4pi
# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>
# unit of DWF is [A^2]
#-------------------------------------------------------
DWF=0

EOF
}
else
{print Fout1 << "EOF";
#!MODULE=kramer
#<!--mcphase.sipf-->
#***************************************************************
# Single Ion Parameter File for Module Kramer for
# mcphas version 3.0
# - program to calculate static magnetic properties
# reference: M. Rotter JMMM 272-276 (2004) 481
# mcdisp version 3.0
# - program to calculate the dispersion of magnetic excitations
# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
# mcdiff version 3.0
# - program to calculate neutron and magnetic xray diffraction
# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R
#***************************************************************
#IONTYPE=$_
CHARGE=$charge
MAGNETIC=$magnetic

# this is a crystal field ground state doublet
# module, parameters are the following 3 matrix
# elements
#
# A=|<+-|Ja|-+>| B=|<+-|Jb|-+>| C=|<+-|Jc|+->|
A = 2.000000
B = 2.543750
C = 1.600000
#----------------
# Lande factor gJ
#----------------
GJ=2

#-------------------------------------------------------
# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)
#                      with s=sin(theta)/lambda=Q/4pi
# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>
# unit of DWF is [A^2]
#-------------------------------------------------------
DWF=0

EOF
}
     ++$i;

  unless (open (Fin, "scatteringlengths.txt")){die "\n error:unable to open $file\n";}

$element=$_;
$element=~s/p\./\./;
$element=~s/\.sipf//;
$element=~s/\d//;
$slrn=0;
$slin=0;
while($line=<Fin>)
      {if($line=~/^$element\s/)
        {@numbers=split(" ",$line);#print @numbers[2]."\n";
         $sl=@numbers[2];
         if($sl=~/i/) {($slrn,$slin)=($sl=~m/([\+\-]?[\d\.]+)([\+\-]?[\d\.]+)i/);}
         else
         {$slrn=$sl;$slin=0;}

        }
      }


      close Fin;


if(abs($slrn)+abs($slin)>0.00000001)
 {
 print  "# Neutron Scattering Length (10^-12 cm) for element $element from internal table:\n";
 print  "SCATTERINGLENGTHREAL=".(0.1*$slrn)."\n";
 print  "SCATTERINGLENGTHIMAG=".(0.1*$slin)."\n";

 print Fout1 "\n\n\n\n#-------------------------------------------------------\n";
 print Fout1 "# Neutron Scattering Length (10^-12 cm) (can be complex)\n";
 print Fout1 "#-------------------------------------------------------\n";
 print Fout1 "SCATTERINGLENGTHREAL=".(0.1*$slrn)."\n";
 print Fout1 "SCATTERINGLENGTHIMAG=".(0.1*$slin)."\n";
 print Fout1 "#  ... note: - if an occupancy other than 1.0 is needed, just reduce \n";
 print Fout1 "#              the scattering length linear accordingly\n\n";
 }
 else
 {print "scatteringlength for element $element not found in internal table";}
     close Fout1;
}

 my $total = 0;
($total+=$_) for @nofatoms;
 my $totalmag = 0;
($totalmag+=$_) for @nofmagatoms;
 $totalnonmag=$total-$totalmag;

open (Fout2, ">mcdiff.in");
print Fout2 << "EOF";
# this file is the input file created by program powdercell2j
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
#! lambda   = 2.4  wavelength (A)
#
#! thetamax = 60   maximum bragg angle (deg)
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
# out*  controls the type of output in  mcdiff.out 
#!out0=1  0..short header, 1...standard header, 2... long header in mcdiff.out
#
# type of desired output in column 1 to 12 (default is h k l d Q 2theta Inuc Imag Itot
#
#!out1=31 out2=32 out3=33 out4=34 out5=35 out6=36 out7=27 out8=28 out9=29 out10=1 out11=0 out12=30 
#
#     ... in out*=n the numbers n have the following meaning:
#            0....LF #
#            1....|NSF|[b] #
#            2....Re(NSF)[b] #
#            3....Im(NSF)[b] #
#            4....|MSF| #
#            5....|MSF.P| #
#            6....Re(MSF.P) #
#            7....Im(MSF.P) #
#            8....|MSFdip| #
#            9....|MSFdip.P| #
#            10....Re(MSFdip.P) #
#            11....Im(MSFdip.P) #
#            12....angl(Q,P)[°] #
#            13....i(MSFxMSF*).P #
#            14....I+ #
#            15....I- #
#            16....I+/I- #
#            17....i(MSFxMSF*)dip.P #
#            18....Idip+ #
#            19....Idip- #
#            20....Idip+/Idip- #
#            21....2*|MSF.P|/sin^2(angl(Q,P) #
#            22....2*|MSFdip.P|/sin^2(angl(Q,P) #
#            23....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-) #
#            24....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-) #
#            25....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-) #
#            26....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-) #
#            27....Inuc(2t)	#
#            28....Imag(2t)   #
#            29....Itot(2t)   #
#            30....Imag_dip(2t) #
#            31....h     #
#            32.... k     #
#            33.... l     #
#            34....d[A]       #
#            35....|Q|[1/A]   #
#            36....2theta     #
#            37....Qi[1/A]    #
#            38....Qj[1/A]     Qi Qj Qk are euclidean components of scattering vector#
#            39....Qk[1/A]     with j||b, k||(a x b) and i normal to k and j#
#            40....T[K]       #
#            41....Ha[T]      #
#            42....Hb[T]      #
#            43....Hc[T]      #
#            44....hprim      #
#            45....kprim      #
#            46....lprim      #
#            47....Itotdip(2t) #
#
#           In the above the intensities I+ and I- are the intensities in a polarised neutron
#           experiment with incident polarisation up (+) or down (-):
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
#
# %SECTION 2% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL
#
#
#! natcryst=$totalnonmag      number of nonmagnetic atoms in primitive crystalographic unit cell
#
# it follows a list of nat lines with nonmagnetic atoms
# ... notes: - if an occupancy other than 1.0 is needed, just reduce
#              the scattering length linear accordingly
#            - da db and dc are not used by the program, dr1,dr2 and dr3
#              refer to the primitive lattice given below
#            - Debye Waller Factor notation: sqr(Intensity) ~ structure factor ~
#              ~sum_n ()n exp(-2 DWFn sin^2(theta) / lambda^2)=EXP (-Wn),
#              relation to other notations: 2*DWF = B = 8 pi^2 <u^2>, units DWF (A^2)
#
#! use_dadbdc=1
#
# Real Imag[scattering length(10^-12cm)]   da(a)    db(b)    dc(c)    dr1(r1)  dr2(r2)  dr3(r3)  DWF(A^2)
EOF

    open (Fout, ">mcphas_all_atoms.j");
    open (Foutmag, ">mcphas_magnetic_atoms.j");
  print "\n\n\n\n\n\n";
  print Fout "#<!--mcphase.mcphas.j-->\n";
  print Fout "#***************************************************************\n";
  print Fout "# Lattice and Exchange Parameter file for\n";
  print Fout "# mcphas version 3.0\n";
  print Fout "# - program to calculate static magnetic properties\n";
  print Fout "# reference: M. Rotter JMMM 272-276 (2004) 481\n";
  print Fout "# mcdisp version 3.0\n";
  print Fout "# - program to calculate the dispersion of magnetic excitations\n";
  print Fout "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n";
  print Fout "#***************************************************************\n";
  print Fout "#\n";
  print Foutmag "# ";
  print Foutmag "#<!--mcphase.mcphas.j-->\n";
  print Foutmag "#***************************************************************\n";
  print Foutmag "# Lattice and Exchange Parameter file for\n";
  print Foutmag "# mcphas version 3.0\n";
  print Foutmag "# - program to calculate static magnetic properties\n";
  print Foutmag "# reference: M. Rotter JMMM 272-276 (2004) 481\n";
  print Foutmag "# mcdisp version 3.0\n";
  print Foutmag "# - program to calculate the dispersion of magnetic excitations\n";
  print Foutmag "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n";
  print Foutmag "#***************************************************************\n";
  print Foutmag "#\n";
  print Foutmag "# ";
  $i=0;foreach(@atomtypetempp)
  {print Fout " $_ (".($nofatoms[$i]/min(@nofatoms)).")";
   print " $_ (".($nofatoms[$i]/min(@nofatoms)).")";++$i;
  }

print STDOUT << "EOF";


 If this is not the formula of your system
 please check the format of the input file

Example:
 No   name       crystal coordinates          cartesian coordinates
                x        y        z           x        y        z
 ------------------------------------------------------------------
 1     Sr1    0.3644   0.0000   0.2500     1.0962  -4.1497  -2.7991
 ...
EOF

  print Fout "\n";
  print Fout "#\n";
  print Fout "# Lattice Constants (A)\n";
  print Fout "#\n";
  print Fout "#! a= $alattice b= $blattice c= $clattice  alpha= $alpha beta= $beta gamma= $gamma\n";
  print Fout "#\n";
  print Fout "#! r1a=   1 r2a= 0 r3a=  0\n";
  print Fout "#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]\n";
  print Fout "#! r1c=   0 r2c= 0 r3c=  1\n";
  print Fout "#\n";
  print Fout "#! nofatoms= $total  nofcomponents=3  number of atoms in primitive unit cell/number of components of each spin\n";
  print Foutmag "\n";
  print Foutmag "#\n";
  print Foutmag "# Lattice Constants (A)\n";
  print Foutmag "#\n";
  print Foutmag "#! a= $alattice b= $blattice c= $clattice  alpha= $alpha beta= $beta gamma= $gamma\n";
  print Foutmag "#\n";
  print Foutmag "#! r1a=   1 r2a= 0 r3a=  0\n";
  print Foutmag "#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]\n";
  print Foutmag "#! r1c=   0 r2c= 0 r3c=  1\n";
  print Foutmag "#\n";
  print Foutmag "#! nofatoms= $totalmag  nofcomponents=3  number of atoms in primitive unit cell/number of components of each spin\n";

$i=0;$nr=1;$nrmag=1;
  foreach(@atomtypetempp)
   {$ii=0;
    foreach(@atomtype)
      {  if($_=~$atomtypetempp[$i])
           {
           print Fout "#********************************************************************* \n";
           print Fout "#ATOM TYPE $_ ; number of the atom in the UNIT CELL = $nr ; number of the atom within this type = $nofatoms[$i]\n";
           print Fout "#! da= $x[$ii] [a] db= $y[$ii] [b] dc= $z[$ii] [c] nofneighbours=0 diagonalexchange=1 sipffilename= $_.sipf\n";

           ($SR)=extract("SCATTERINGLENGTHREAL","$_.sipf");
           ($SI)=extract("SCATTERINGLENGTHIMAG","$_.sipf");
           ($magnetic)=extract("MAGNETIC","$_.sipf");
           if($magnetic)
           {print Fout2 "# $SR $SI $x[$ii] $y[$ii] $z[$ii] $x[$ii] $y[$ii] $z[$ii] 0  # $_.sipf\n";
           print Foutmag "#********************************************************************* \n";
           print Foutmag "#ATOM TYPE $_ ; number of the atom in the UNIT CELL = $nrmag ; number of the atom within this type = $nofatoms[$i]\n";
           print Foutmag "#! da= $x[$ii] [a] db= $y[$ii] [b] dc= $z[$ii] [c] nofneighbours=0 diagonalexchange=1 sipffilename= $_.sipf\n";
           ++$nrmag;
           }
           else
           {print Fout2 "$SR $SI $x[$ii] $y[$ii] $z[$ii] $x[$ii] $y[$ii] $z[$ii] 0  # $_.sipf\n";}

           ++$nr; }
       ++$ii;
      }
    ++$i;
   }
print Fout "#********************************************************************* \n";
close Fout;
print Foutmag "#********************************************************************* \n";
close Foutmag;

 
print "\n\n";
print "The folowing file have been created by this program:\n\n";
print "mcphas.j - the input file for McPhase program\n";
print "\n";
foreach(@atomtypetempp) {print "$_.sipf - contain the oxidation state of $_\n";}


print STDOUT << "EOF";

 mcphas_all.j and mcphas_magnetic_atoms.j=mcphas.j created

EOF
print STDOUT "press enter to generate mcdiff.in file, too\n";
<STDIN>;
print " Please enter the magnetic supercell dimension na nb nc and the propagation vector (h,k,l):\n";
print " na (>=1)? (1)\n";$ans = <STDIN>; $ans =~ s/[\R\n\r]*//g;  
        if($ans=~/[yYtT]?/) {  $na = 1;        } else { $na = $ans;}
print " nb (>=1) ? (1)\n";$ans = <STDIN>; $ans =~ s/[\R\n\r]*//g;  
         if($ans=~/[yYtT]?/) {  $nb = 1;        } else { $nb = $ans;}
print " nc (>=1) ? (1) \n";$ans = <STDIN>; $ans =~ s/[\R\n\r]*//g;  
         if($ans=~/[yYtT]?/) {  $nc = 1;        } else { $nc = $ans;}
print " h ? (0)\n";$ans = <STDIN>; $ans =~ s/[\R\n\r]*//g;  
         if($ans=~/[yYtT]?/) {  $h = 0;        } else { $h = $ans;}
print " k ? (0)\n";$ans = <STDIN>; $ans =~ s/[\R\n\r]*//g;  
         if($ans=~/[yYtT]?/) {  $k = 0;        } else { $k = $ans;}
print " l ? (0)\n";$ans = <STDIN>; $ans =~ s/[\R\n\r]*//g;  
         if($ans=~/[yYtT]?/) {  $l = 0;        } else { $l = $ans;}

copy("mcphas_magnetic_atoms.j","mcphas.j");
system ("spinsfromq $na $nb $nc $h $k $l > results/mcphas.mf");
#print "spinsfromq $na $nb $nc $h $k $l > results/mcphas.mf\n";
#system ("fact 1 0 powdercell2j.sps");
#system ("spins 0 0 0 0 powdercell2j.sps > powdercell2j.spo");
#system ("javaview results/spins.jvx");
system("display_densities -M 1 0 0 0 > results/powdercell2j.spo");
mydel  ("results/powdercell2j.spo");
mydel  ("results/mcphas.mf");
mydel ("scatteringlengths.txt");

print Fout2 << "EOF";
#
#
# %SECTION 3% DESCRIPTION OF THE LATTICE
#
#
# Note: what follows here may directly be taken from the output of program spins
#       (file spins.out) or charges (file charges.out)
# -----------------------------------------------------------------------------
EOF
open (Fin,"results/spins.out");
while($line=<Fin>) {print Fout2 $line;}
close Fin, Fout2;

print STDOUT << "EOF";

END OF PROGRAM powdercell2j

Created files:

 mcphas_all.j                      ... contains all atoms (magnetic and nonmag)
 mcphas_magnetic_atoms.j=mcphas.j  ... contains only magnetic atoms

   running \"makenn R\" command  will create from mcphas.j a number of files
   \"makenn.aN.pc\" equal with the number of the atoms in the unit cell.
   in this file you find the neighbours of the N atom in the unit cell.
   It is important to know what number corresponds to each atom type!


 mcdiff.in                              ... input file for mcdiff
   running javaview results/spins.jvx displays the corresponding
   spinstructure (collinear), !! mind to edit formfactor information
   in the magnetic ions .sipf files before running mcdiff !!



EOF



sub min
{
    my (@numbers);

    @numbers = @_;

    my ($min);

    $min = $numbers[0];

    foreach my $i (@numbers)
    {
        if ($i < $min)
        {
            $min = $i;
        }
    }

    return ($min);
}

sub sum
{
    my (@numbers);

    @numbers = @_;

    my ($sum);



    foreach my $i (@numbers)
    { $sum+=$i;
    }

    return ($min);
}

sub mydel  { my ($file1)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;
                               return system("del ".$file1);
                              }
                 else
                              {return system("rm ".$file1);
                              }

           }

sub scattlengthoutput {
open (Fout, ">scatteringlengths.txt");
print Fout << "EOF";
#Column 	Unit 	Quantity
#1 	--- 	Isotope
#2 	--- 	Natural abundance (For radioisotopes the half-life is given instead)
#3 	fm 	bound coherent scattering length
#4 	fm 	bound incoherent scattering length
#5 	barn 	bound coherent scattering cross section
#6 	barn 	bound incoherent scattering cross section
#7 	barn 	total bound scattering cross section
#8 	barn 	absorption cross section for 2200 m/s neutrons
###
##
#
#Note: 1fm=1E-15 m, 1barn=1E-24 cm^2, scattering lengths and cross sections in parenthesis are uncertainties.
#Neutron scattering lengths and cross sections
#Isotope conc 	Coh b 	Inc b 	Coh xs 	Inc xs 	Scatt xs 	Abs xs
H 	--- 	-3.7390 	--- 	1.7568 	80.26 	82.02 	0.3326
1H 	99.985 	-3.7406 	25.274 	1.7583 	80.27 	82.03 	0.3326
2H 	0.015 	6.671 	4.04 	5.592 	2.05 	7.64 	0.000519
3H 	(12.32  	4.792 	-1.04 	2.89 	0.14 	3.03 	0
He 	--- 	3.26(3) 	--- 	1.34 	0 	1.34 	0.00747
3He 	0.00014 	5.74-1.483i 	-2.5+2.568i 	4.42 	1.6 	6 	5333.(7.)
4He 	99.99986 	3.26 	0 	1.34 	0 	1.34 	0
Li 	--- 	-1.90 	--- 	0.454 	0.92 	1.37 	70.5
6Li 	7.5 	2.00-0.261i 	-1.89+0.26i 	0.51 	0.46 	0.97 	940.(4.)
7Li 	92.5 	-2.22 	-2.49 	0.619 	0.78 	1.4 	0.0454
Be 	100 	7.79 	0.12 	7.63 	0.0018 	7.63 	0.0076
B 	--- 	5.30-0.213i 	--- 	3.54 	1.7 	5.24 	767.(8.)
10B 	20 	-0.1-1.066i 	-4.7+1.231i 	0.144 	3 	3.1 	3835.(9.)
11B 	80 	6.65 	-1.3 	5.56 	0.21 	5.77 	0.0055
C 	--- 	6.6460 	--- 	5.551 	0.001 	5.551 	0.0035
12C 	98.9 	6.6511 	0 	5.559 	0 	5.559 	0.00353
13C 	1.1 	6.19 	-0.52 	4.81 	0.034 	4.84 	0.00137
N 	--- 	9.36 	--- 	11.01 	0.5 	11.51 	1.9
14N 	99.63 	9.37 	2.0 	11.03 	0.5 	11.53 	1.91
15N 	0.37 	6.44 	-0.02 	5.21 	0.00005 	5.21 	0.000024
O 	--- 	5.803 	--- 	4.232 	0.0008 	4.232 	0.00019
16O 	99.762 	5.803 	0 	4.232 	0 	4.232 	0.0001
17O 	0.038 	5.78 	0.18 	4.2 	0.004 	4.2 	0.236
18O 	0.2 	5.84 	0 	4.29 	0 	4.29 	0.00016
F 	100 	5.654 	-0.082 	4.017 	0.0008 	4.018 	0.0096
Ne 	--- 	4.566 	--- 	2.62 	0.008 	2.628 	0.039
20Ne 	90.51 	4.631 	0 	2.695 	0 	2.695 	0.036
21Ne 	0.27 	6.66 	(+/-)0.6 	5.6 	0.05 	5.7 	0.67
22Ne 	9.22 	3.87 	0 	1.88 	0 	1.88 	0.046
Na 	100 	3.63 	3.59 	1.66 	1.62 	3.28 	0.53
Mg 	--- 	5.375 	--- 	3.631 	0.08 	3.71 	0.063
24Mg 	78.99 	5.66 	0 	4.03 	0 	4.03 	0.05
25Mg 	10 	3.62 	1.48 	1.65 	0.28 	1.93 	0.19
26Mg 	11.01 	4.89 	0 	3 	0 	3 	0.0382
Al 	100 	3.449 	0.256 	1.495 	0.0082 	1.503 	0.231
Si 	--- 	4.1491 	--- 	2.163 	0.004 	2.167 	0.171
28Si 	92.23 	4.107 	0 	2.12 	0 	2.12 	0.177
29Si 	4.67 	4.70 	0.09 	2.78 	0.001 	2.78 	0.101
30Si 	3.1 	4.58 	0 	2.64 	0 	2.64 	0.107
P 	100 	5.13 	0.2 	3.307 	0.005 	3.312 	0.172
S 	--- 	2.847 	--- 	1.0186 	0.007 	1.026 	0.53
32S 	95.02 	2.804 	0 	0.988 	0 	0.988 	0.54
33S 	0.75 	4.74 	1.5 	2.8 	0.3 	3.1 	0.54
34S 	4.21 	3.48 	0 	1.52 	0 	1.52 	0.227
36S 	0.02 	3.(1.) 	0 	1.1 	0 	1.1 	0.15
Cl 	--- 	9.5770 	--- 	11.5257 	5.3 	16.8 	33.5
35Cl 	75.77 	11.65 	6.1 	17.06 	4.7 	21.8 	44.1
37Cl 	24.23 	3.08 	0.1 	1.19 	0.001 	1.19 	0.433
Ar 	--- 	1.909 	--- 	0.458 	0.225 	0.683 	0.675
36Ar 	0.337 	24.90 	0 	77.9 	0 	77.9 	5.2
38Ar 	0.063 	3.5 	0 	1.5(3.1) 	0 	1.5(3.1) 	0.8
40Ar 	99.6 	1.830 	0 	0.421 	0 	0.421 	0.66
K 	--- 	3.67 	--- 	1.69 	0.27 	1.96 	2.1
39K 	93.258 	3.74 	1.4 	1.76 	0.25 	2.01 	2.1
40K 	0.012 	3.(1.) 	--- 	1.1 	0.5 	1.6 	35.(8.)
41K 	6.73 	2.69 	1.5 	0.91 	0.3 	1.2 	1.46
Ca 	--- 	4.70 	--- 	2.78 	0.05 	2.83 	0.43
40Ca 	96.941 	4.80 	0 	2.9 	0 	2.9 	0.41
42Ca 	0.647 	3.36 	0 	1.42 	0 	1.42 	0.68
43Ca 	0.135 	-1.56 	--- 	0.31 	0.5 	0.8 	6.2
44Ca 	2.086 	1.42 	0 	0.25 	0 	0.25 	0.88
46Ca 	0.004 	3.6 	0 	1.6 	0 	1.6 	0.74
48Ca 	0.187 	0.39 	0 	0.019 	0 	0.019 	1.09
Sc 	100 	12.29 	-6.0 	19 	4.5 	23.5 	27.5
Ti 	--- 	-3.438 	--- 	1.485 	2.87 	4.35 	6.09
46Ti 	8.2 	4.93 	0 	3.05 	0 	3.05 	0.59
47Ti 	7.4 	3.63 	-3.5 	1.66 	1.5 	3.2 	1.7
48Ti 	73.8 	-6.08 	0 	4.65 	0 	4.65 	7.84
49Ti 	5.4 	1.04 	5.1 	0.14 	3.3 	3.4 	2.2
50Ti 	5.2 	6.18 	0 	4.8 	0 	4.8 	0.179
V 	--- 	-0.3824 	--- 	0.0184 	5.08 	5.1 	5.08
50V 	0.25 	7.6 	--- 	7.3(1.1) 	0.5 	7.8(1.0) 	60.(40.)
51V 	99.75 	-0.402 	6.35 	0.0203 	5.07 	5.09 	4.9
Cr 	--- 	3.635 	--- 	1.66 	1.83 	3.49 	3.05
50Cr 	4.35 	-4.50 	0 	2.54 	0 	2.54 	15.8
52Cr 	83.79 	4.920 	0 	3.042 	0 	3.042 	0.76
53Cr 	9.5 	-4.20 	6.87 	2.22 	5.93 	8.15 	18.1(1.5)
54Cr 	2.36 	4.55 	0 	2.6 	0 	2.6 	0.36
Mn 	100 	-3.73 	1.79 	1.75 	0.4 	2.15 	13.3
Fe 	--- 	9.45 	--- 	11.22 	0.4 	11.62 	2.56
54Fe 	5.8 	4.2 	0 	2.2 	0 	2.2 	2.25
56Fe 	91.7 	9.94 	0 	12.42 	0 	12.42 	2.59
57Fe 	2.2 	2.3 	--- 	0.66 	0.3 	1 	2.48
58Fe 	0.3 	15.(7.) 	0 	28 	0 	28.(26.) 	1.28
Co 	100 	2.49 	-6.2 	0.779 	4.8 	5.6 	37.18
Ni 	--- 	10.3 	--- 	13.3 	5.2 	18.5 	4.49
58Ni 	68.27 	14.4 	0 	26.1 	0 	26.1 	4.6
60Ni 	26.1 	2.8 	0 	0.99 	0 	0.99 	2.9
61Ni 	1.13 	7.60 	(+/-)3.9 	7.26 	1.9 	9.2 	2.5
62Ni 	3.59 	-8.7 	0 	9.5 	0 	9.5 	14.5
64Ni 	0.91 	-0.37 	0 	0.017 	0 	0.017 	1.52
Cu 	--- 	7.718 	--- 	7.485 	0.55 	8.03 	3.78
63Cu 	69.17 	6.43 	0.22 	5.2 	0.006 	5.2 	4.5
65Cu 	30.83 	10.61 	1.79 	14.1 	0.4 	14.5 	2.17
Zn 	--- 	5.680 	--- 	4.054 	0.077 	4.131 	1.11
64Zn 	48.6 	5.22 	0 	3.42 	0 	3.42 	0.93
66Zn 	27.9 	5.97 	0 	4.48 	0 	4.48 	0.62
67Zn 	4.1 	7.56 	-1.50 	7.18 	0.28 	7.46 	6.8
68Zn 	18.8 	6.03 	0 	4.57 	0 	4.57 	1.1
70Zn 	0.6 	6.(1.) 	0 	4.5 	0 	4.5(1.5) 	0.092
Ga 	--- 	7.288 	--- 	6.675 	0.16 	6.83 	2.75
69Ga 	60.1 	7.88 	-0.85 	7.8 	0.091 	7.89 	2.18
71Ga 	39.9 	6.40 	-0.82 	5.15 	0.084 	5.23 	3.61
Ge 	--- 	8.185 	--- 	8.42 	0.18 	8.6 	2.2
70Ge 	20.5 	10.0 	0 	12.6 	0 	12.6 	3
72Ge 	27.4 	8.51 	0 	9.1 	0 	9.1 	0.8
73Ge 	7.8 	5.02 	3.4 	3.17 	1.5 	4.7 	15.1
74Ge 	36.5 	7.58 	0 	7.2 	0 	7.2 	0.4
76Ge 	7.8 	8.2 	0 	8.(3.) 	0 	8.(3.) 	0.16
As 	100 	6.58 	-0.69 	5.44 	0.06 	5.5 	4.5
Se 	--- 	7.970 	--- 	7.98 	0.32 	8.3 	11.7
74Se 	0.9 	0.8 	0 	0.1 	0 	0.1 	51.8(1.2)
76Se 	9 	12.2 	0 	18.7 	0 	18.7 	85.(7.)
77Se 	7.6 	8.25 	(+/-)0.6(1.6) 	8.6 	0.05 	8.65 	42.(4.)
78Se 	23.5 	8.24 	0 	8.5 	0 	8.5 	0.43
80Se 	49.6 	7.48 	0 	7.03 	0 	7.03 	0.61
82Se 	9.4 	6.34 	0 	5.05 	0 	5.05 	0.044
Br 	--- 	6.795 	--- 	5.8 	0.1 	5.9 	6.9
79Br 	50.69 	6.80 	-1.1 	5.81 	0.15 	5.96 	11
81Br 	49.31 	6.79 	0.6 	5.79 	0.05 	5.84 	2.7
Kr 	--- 	7.81 	--- 	7.67 	0.01 	7.68 	25.(1.)
78Kr 	0.35 	--- 	0 	--- 	0 	--- 	6.4
80Kr 	2.25 	--- 	0 	--- 	0 	--- 	11.8
82Kr 	11.6 	--- 	0 	--- 	0 	--- 	29.(20.)
83Kr 	11.5 	--- 	--- 	--- 	--- 	--- 	185.(30.)
84Kr 	57 	--- 	0 	--- 	0 	6.6 	0.113
86Kr 	17.3 	8.1 	0 	8.2 	0 	8.2 	0.003
Rb 	--- 	7.09 	--- 	6.32 	0.5 	6.8 	0.38
85Rb 	72.17 	7.03 	--- 	6.2 	0.5 	6.7 	0.48
87Rb 	27.83 	7.23 	--- 	6.6 	0.5 	7.1 	0.12
Sr 	--- 	7.02 	--- 	6.19 	0.06 	6.25 	1.28
84Sr 	0.56 	7.(1.) 	0 	6.(2.) 	0 	6.(2.) 	0.87
86Sr 	9.86 	5.67 	0 	4.04 	0 	4.04 	1.04
87Sr 	7 	7.40 	--- 	6.88 	0.5 	7.4 	16.(3.)
88Sr 	82.58 	7.15 	0 	6.42 	0 	6.42 	0.058
Y 	100 	7.75 	1.1 	7.55 	0.15 	7.7 	1.28
Zr 	--- 	7.16 	--- 	6.44 	0.02 	6.46 	0.185
90Zr 	51.45 	6.4 	0 	5.1 	0 	5.1 	0.011
91Zr 	11.32 	8.7 	-1.08 	9.5 	0.15 	9.7 	1.17
92Zr 	17.19 	7.4 	0 	6.9 	0 	6.9 	0.22
94Zr 	17.28 	8.2 	0 	8.4 	0 	8.4 	0.0499
96Zr 	2.76 	5.5 	0 	3.8 	0 	3.8 	0.0229
Nb 	100 	7.054 	-0.139 	6.253 	0.0024 	6.255 	1.15
Mo 	--- 	6.715 	--- 	5.67 	0.04 	5.71 	2.48
92Mo 	14.84 	6.91 	0 	6 	0 	6 	0.019
94Mo 	9.25 	6.80 	0 	5.81 	0 	5.81 	0.015
95Mo 	15.92 	6.91 	--- 	6 	0.5 	6.5 	13.1
96Mo 	16.68 	6.20 	0 	4.83 	0 	4.83 	0.5
97Mo 	9.55 	7.24 	--- 	6.59 	0.5 	7.1 	2.5
98Mo 	24.13 	6.58 	0 	5.44 	0 	5.44 	0.127
100Mo 	9.63 	6.73 	0 	5.69 	0 	5.69 	0.4
Tc 	(2.3E5  	6.8 	--- 	5.8 	0.5 	6.3 	20.(1.)
Ru 	--- 	7.03 	--- 	6.21 	0.4 	6.6 	2.56
96Ru 	5.5 	--- 	0 	--- 	0 	--- 	0.28
98Ru 	1.9 	--- 	0 	--- 	0 	--- 	<8.
99Ru 	12.7 	--- 	--- 	--- 	--- 	--- 	6.9(1.0)
100Ru 	12.6 	--- 	0 	--- 	0 	--- 	4.8
101Ru 	17 	--- 	--- 	--- 	--- 	--- 	3.3
102Ru 	31.6 	--- 	0 	--- 	0 	144.8 	1.17
104Ru 	18.7 	--- 	0 	--- 	0 	4.483 	0.31
Rh 	100 	5.88 	--- 	4.34 	0.3 	4.6 	144.8
Pd 	--- 	5.91 	--- 	4.39 	0.093 	4.48 	6.9
102Pd 	1.02 	7.7(7) 	0 	7.5(1.4) 	0 	7.5(1.4) 	3.4
104Pd 	11.14 	7.7(7) 	0 	7.5(1.4) 	0 	7.5(1.4) 	0.6
105Pd 	22.33 	5.5 	-2.6(1.6) 	3.8 	0.8 	4.6(1.1) 	20.(3.)
106Pd 	27.33 	6.4 	0 	5.1 	0 	5.1 	0.304
108Pd 	26.46 	4.1 	0 	2.1 	0 	2.1 	8.55
110Pd 	11.72 	7.7(7) 	0 	7.5(1.4) 	0 	7.5(1.4) 	0.226
Ag 	--- 	5.922 	--- 	4.407 	0.58 	4.99 	63.3
107Ag 	51.83 	7.555 	1.00 	7.17 	0.13 	7.3 	37.6(1.2)
109Ag 	48.17 	4.165 	-1.60 	2.18 	0.32 	2.5 	91.0(1.0)
Cd 	--- 	4.87-0.70i 	--- 	3.04 	3.46 	6.5 	2520.(50.)
106Cd 	1.25 	5.(2.) 	0 	3.1 	0 	3.1(2.5) 	1
108Cd 	0.89 	5.4 	0 	3.7 	0 	3.7 	1.1
110Cd 	12.51 	5.9 	0 	4.4 	0 	4.4 	11
111Cd 	12.81 	6.5 	--- 	5.3 	0.3 	5.6 	24
112Cd 	24.13 	6.4 	0 	5.1 	0 	5.1 	2.2
113Cd 	12.22 	-8.0-5.73i 	--- 	12.1 	0.3 	12.4 	20600.(400.)
114Cd 	28.72 	7.5 	0 	7.1 	0 	7.1 	0.34
116Cd 	7.47 	6.3 	0 	5 	0 	5 	0.075
In 	--- 	4.065-0.0539i 	--- 	2.08 	0.54 	2.62 	193.8(1.5)
113In 	4.3 	5.39 	(+/-)0.017 	3.65 	0.000037 	3.65 	12.0(1.1)
115In 	95.7 	4.01-0.0562i 	-2.1 	2.02 	0.55 	2.57 	202.(2.)
Sn 	--- 	6.225 	--- 	4.871 	0.022 	4.892 	0.626
112Sn 	1 	6.(1.) 	0 	4.5(1.5) 	0 	4.5(1.5) 	1
114Sn 	0.7 	6.2 	0 	4.8 	0 	4.8 	0.114
115Sn 	0.4 	6.(1.) 	--- 	4.5(1.5) 	0.3 	4.8(1.5) 	30.(7.)
116Sn 	14.7 	5.93 	0 	4.42 	0 	4.42 	0.14
117Sn 	7.7 	6.48 	--- 	5.28 	0.3 	5.6 	2.3
118Sn 	24.3 	6.07 	0 	4.63 	0 	4.63 	0.22
119Sn 	8.6 	6.12 	--- 	4.71 	0.3 	5 	2.2
120Sn 	32.4 	6.49 	0 	5.29 	0 	5.29 	0.14
122Sn 	4.6 	5.74 	0 	4.14 	0 	4.14 	0.18
124Sn 	5.6 	5.97 	0 	4.48 	0 	4.48 	0.133
Sb 	--- 	5.57 	--- 	3.9 	0.007 	3.9 	4.91
121Sb 	57.3 	5.71 	-0.05 	4.1 	0.0003 	4.1 	5.75
123Sb 	42.7 	5.38 	-0.10 	3.64 	0.001 	3.64 	3.8
Te 	--- 	5.80 	--- 	4.23 	0.09 	4.32 	4.7
120Te 	0.096 	5.3 	0 	3.5 	0 	3.5 	2.3
122Te 	2.6 	3.8 	0 	1.8 	0 	1.8 	3.4
123Te 	0.908 	-0.05-0.116i 	-2.04 	0.002 	0.52 	0.52 	418.(30.)
124Te 	4.816 	7.96 	0 	8 	0 	8 	6.8(1.3)
125Te 	7.14 	5.02 	-0.26 	3.17 	0.008 	3.18 	1.55
126Te 	18.95 	5.56 	0 	3.88 	0 	3.88 	1.04
128Te 	31.69 	5.89 	0 	4.36 	0 	4.36 	0.215
130Te 	33.8 	6.02 	0 	4.55 	0 	4.55 	0.29
I 	100 	5.28 	1.58 	3.5 	0.31 	3.81 	6.15
Xe 	--- 	4.92 	3.04 	2.96 	0 	--- 	23.9(1.2)
124Xe 	0.1 	--- 	0 	--- 	0 	--- 	165.(20.)
126Xe 	0.09 	--- 	0 	--- 	0 	--- 	3.5
128Xe 	1.91 	--- 	0 	--- 	0 	--- 	<8
129Xe 	26.4 	--- 	--- 	--- 	--- 	--- 	21.(5.)
130Xe 	4.1 	--- 	0 	--- 	0 	--- 	<26.
131Xe 	21.2 	--- 	--- 	--- 	--- 	--- 	85.(10.)
132Xe 	26.9 	--- 	0 	--- 	0 	--- 	0.45
134Xe 	10.4 	--- 	0 	--- 	0 	--- 	0.265
136Xe 	8.9 	--- 	0 	--- 	0 	--- 	0.26
Cs 	100 	5.42 	1.29 	3.69 	0.21 	3.9 	29.0(1.5)
Ba 	--- 	5.07 	--- 	3.23 	0.15 	3.38 	1.1
130Ba 	0.11 	-3.6 	0 	1.6 	0 	1.6 	30.(5.)
132Ba 	0.1 	7.8 	0 	7.6 	0 	7.6 	7
134Ba 	2.42 	5.7 	0 	4.08 	0 	4.08 	2.0(1.6)
135Ba 	6.59 	4.67 	--- 	2.74 	0.5 	3.2 	5.8
136Ba 	7.85 	4.91 	0 	3.03 	0 	3.03 	0.68
137Ba 	11.23 	6.83 	--- 	5.86 	0.5 	6.4 	3.6
138Ba 	71.7 	4.84 	0 	2.94 	0 	2.94 	0.27
La 	--- 	8.24 	--- 	8.53 	1.13 	9.66 	8.97
138La 	0.09 	8.(2.) 	--- 	8.(4.) 	0.5 	8.5(4.0) 	57.(6.)
139La 	99.91 	8.24 	3.0 	8.53 	1.13 	9.66 	8.93
Ce 	--- 	4.84 	--- 	2.94 	0.001 	2.94 	0.63
136Ce 	0.19 	5.80 	0 	4.23 	0 	4.23 	7.3(1.5)
138Ce 	0.25 	6.70 	0 	5.64 	0 	5.64 	1.1
140Ce 	88.48 	4.84 	0 	2.94 	0 	2.94 	0.57
142Ce 	11.08 	4.75 	0 	2.84 	0 	2.84 	0.95
Pr 	100 	4.58 	-0.35 	2.64 	0.015 	2.66 	11.5
Nd 	--- 	7.69 	--- 	7.43 	9.2 	16.6 	50.5(1.2)
142Nd 	27.16 	7.7 	0 	7.5 	0 	7.5 	18.7
143Nd 	12.18 	14.(2.) 	(+/-)21.(1.) 	25.(7.) 	55.(7.) 	80.(2.) 	337.(10.)
144Nd 	23.8 	2.8 	0 	1 	0 	1 	3.6
145Nd 	8.29 	14.(2.) 	--- 	25.(7.) 	5.(5.) 	30.(9.) 	42.(2.)
146Nd 	17.19 	8.7 	0 	9.5 	0 	9.5 	1.4
148Nd 	5.75 	5.7 	0 	4.1 	0 	4.1 	2.5
150Nd 	5.63 	5.3 	0 	3.5 	0 	3.5 	1.2
Pm 	(2.62  	12.6 	(+/-)3.2(2.5) 	20.0(1.3) 	1.3(2.0) 	21.3(1.5) 	168.4(3.5)
Sm 	--- 	0.80-1.65i 	--- 	0.422 	39.(3.) 	39.(3.) 	5922.(56.)
144Sm 	3.1 	-3.(4.) 	0 	1.(3.) 	0 	1.(3.) 	0.7
147Sm 	15.1 	14.(3.) 	(+/-)11.(7.) 	25.(11.) 	143(19.) 	39.(16.) 	57.(3.)
148Sm 	11.3 	-3.(4.) 	0 	1.(3.) 	0 	1.(3.) 	2.4
149Sm 	13.9 	-19.2-11.7i 	(+/-)31.4-10.3i 	63.5 	137.(5.) 	200.(5.) 	42080.(400.)
150Sm 	7.4 	14.(3.) 	0 	25.(11.) 	0 	25.(11.) 	104.(4.)
152Sm 	26.6 	-5.0 	0 	3.1 	0 	3.1 	206.(6.)
154Sm 	22.6 	9.3 	0 	11.(2.) 	0 	11.(2.) 	8.4
Eu 	--- 	7.22-1.26i 	--- 	6.57 	2.5 	9.2 	4530.(40.)
151Eu 	47.8 	6.13-2.53i 	(+/-)4.5-2.14i 	5.5 	3.1 	8.6 	9100.(100.)
153Eu 	52.2 	8.22 	(+/-)3.2 	8.5 	1.3 	9.8 	312.(7.)
Gd 	--- 	6.5-13.82i 	--- 	29.3 	151.(2.) 	180.(2.) 	49700.(125.)
152Gd 	0.2 	10.(3.) 	0 	13.(8.) 	0 	13.(8.) 	735.(20.)
154Gd 	2.1 	10.(3.) 	0 	13.(8.) 	0 	13.(8.) 	85.(12.)
155Gd 	14.8 	6.0-17.0i 	(+/-)5.(5.)-13.16i 	40.8 	25.(6.) 	66.(6.) 	61100.(400.)
156Gd 	20.6 	6.3 	0 	5 	0 	5 	1.5(1.2)
157Gd 	15.7 	-1.14-71.9i 	(+/-)5.(5.)-55.8i 	650.(4.) 	394.(7.) 	1044.(8.) 	259000.(700.)
158Gd 	24.8 	9.(2.) 	0 	10.(5.) 	0 	10.(5.) 	2.2
160Gd 	21.8 	9.15 	0 	10.52 	0 	10.52 	0.77
Tb 	100 	7.38 	-0.17 	6.84 	0.004 	6.84 	23.4
Dy 	--- 	16.9-0.276i 	--- 	35.9 	54.4(1.2) 	90.3 	994.(13.)
156Dy 	0.06 	6.1 	0 	4.7 	0 	4.7 	33.(3.)
158Dy 	0.1 	6.(4.) 	0 	5.(6.) 	0 	5.(6.) 	43.(6.)
160Dy 	2.34 	6.7 	0 	5.6 	0 	5.6 	56.(5.)
161Dy 	19 	10.3 	(+/-)4.9 	13.3 	3.(1.) 	16.(1.) 	600.(25.)
162Dy 	25.5 	-1.4 	0 	0.25 	0 	0.25 	194.(10.)
163Dy 	24.9 	5.0 	1.3 	3.1 	0.21 	3.3 	124.(7.)
164Dy 	28.1 	49.4-0.79i 	0 	307.(3.) 	0 	307.(3.) 	2840.(40.)
Ho 	100 	8.01 	-1.70 	8.06 	0.36 	8.42 	64.7(1.2)
Er 	--- 	7.79 	--- 	7.63 	1.1 	8.7 	159.(4.)
162Er 	0.14 	8.8 	0 	9.7 	0 	9.7 	19.(2.)
164Er 	1.56 	8.2 	0 	8.4 	0 	8.4 	13.(2.)
166Er 	33.4 	10.6 	0 	14.1 	0 	14.1 	19.6(1.5)
167Er 	22.9 	3.0 	1.0 	1.1 	0.13 	1.2 	659.(16.)
168Er 	27.1 	7.4 	0 	6.9 	0 	6.9 	2.74
170Er 	14.9 	9.6 	0 	11.6 	0 	11.6(1.2) 	5.8
Tm 	100 	7.07 	0.9 	6.28 	0.1 	6.38 	100.(2.)
Yb 	--- 	12.43 	--- 	19.42 	4 	23.4 	34.8
168Yb 	0.14 	-4.07-0.62i 	0 	2.13 	0 	2.13 	2230.(40.)
170Yb 	3.06 	6.77 	0 	5.8 	0 	5.8 	11.4(1.0)
171Yb 	14.3 	9.66 	-5.59 	11.7 	3.9 	15.6 	48.6(2.5)
172Yb 	21.9 	9.43 	0 	11.2 	0 	11.2 	0.8
173Yb 	16.1 	9.56 	-5.3 	11.5 	3.5 	15 	17.1(1.3)
174Yb 	31.8 	19.3 	0 	46.8 	0 	46.8 	69.4(5.0)
176Yb 	12.7 	8.72 	0 	9.6 	0 	9.6 	2.85
Lu 	--- 	7.21 	--- 	6.53 	0.7 	7.2 	74.(2.)
175Lu 	97.39 	7.24 	(+/-)2.2 	6.59 	0.6 	7.2 	21.(3.)
176Lu 	2.61 	6.1-0.57i 	(+/-)3.0+0.61i 	4.7 	1.2 	5.9 	2065.(35.)
Hf 	--- 	7.7 	--- 	7.6 	2.6 	10.2 	104.1
174Hf 	0.2 	10.9(1.1) 	0 	15.(3.) 	0 	15.(3.) 	561.(35.)
176Hf 	5.2 	6.61 	0 	5.5 	0 	5.5 	23.5(3.1)
177Hf 	18.6 	0.8(1.0) 	(+/-)0.9(1.3) 	0.1 	0.1 	0.2 	373.(10.)
178Hf 	27.1 	5.9 	0 	4.4 	0 	4.4 	84.(4.)
179Hf 	13.7 	7.46 	(+/-)1.06 	7 	0.14 	7.1 	41.(3.)
180Hf 	35.2 	13.2 	0 	21.9 	0 	21.9(1.0) 	13.04
Ta 	--- 	6.91 	--- 	6 	0.01 	6.01 	20.6
180Ta 	0.012 	7.(2.) 	--- 	6.2 	0.5 	7.(4.) 	563.(60.)
181Ta 	99.988 	6.91 	-0.29 	6 	0.011 	6.01 	20.5
W 	--- 	4.86 	--- 	2.97 	1.63 	4.6 	18.3
180W 	0.1 	5.(3.) 	0 	3.(4.) 	0 	3.(4.) 	30.(20.)
182W 	26.3 	6.97 	0 	6.1 	0 	6.1 	20.7
183W 	14.3 	6.53 	--- 	5.36 	0.3 	5.7 	10.1
184W 	30.7 	7.48 	0 	7.03 	0 	7.03 	1.7
186W 	28.6 	-0.72 	0 	0.065 	0 	0.065 	37.9
Re 	--- 	9.2 	--- 	10.6 	0.9 	11.5 	89.7(1.)
185Re 	37.4 	9.0 	(+/-)2.0 	10.2 	0.5 	10.7 	112.(2.)
187Re 	62.6 	9.3 	(+/-)2.8 	10.9 	1 	11.9 	76.4(1.)
Os 	--- 	10.7 	--- 	14.4 	0.3 	14.7 	16
184Os 	0.02 	10.(2.) 	0 	13.(5.) 	0 	13.(5.) 	3000.(150.)
186Os 	1.58 	11.6(1.7) 	0 	17.(5.) 	0 	17.(5.) 	80.(13.)
187Os 	1.6 	10.(2.) 	--- 	13.(5.) 	0.3 	13.(5.) 	320.(10.)
188Os 	13.3 	7.6 	0 	7.3 	0 	7.3 	4.7
189Os 	16.1 	10.7 	--- 	14.4 	0.5 	14.9 	25.(4.)
190Os 	26.4 	11.0 	0 	15.2 	0 	15.2 	13.1
192Os 	41 	11.5 	0 	16.6 	0 	16.6(1.2) 	2
Ir 	--- 	10.6 	--- 	14.1 	0.(3.) 	14.(3.) 	425.(2.)
191Ir 	37.3 	--- 	--- 	--- 	--- 	--- 	954.(10.)
193Ir 	62.7 	--- 	--- 	--- 	--- 	--- 	111.(5.)
Pt 	--- 	9.60 	--- 	11.58 	0.13 	11.71 	10.3
190Pt 	0.01 	9.0 	0 	10.(2.) 	0 	10.(2.) 	152.(4.)
192Pt 	0.79 	9.9 	0 	12.3(1.2) 	0 	12.3(1.2) 	10.0(2.5)
194Pt 	32.9 	10.55 	0 	14 	0 	14 	1.44
195Pt 	33.8 	8.83 	-1.00 	9.8 	0.13 	9.9 	27.5(1.2)
196Pt 	25.3 	9.89 	0 	12.3 	0 	12.3 	0.72
198Pt 	7.2 	7.8 	0 	7.6 	0 	7.6 	3.66
Au 	100 	7.63 	-1.84 	7.32 	0.43 	7.75 	98.65
Hg 	--- 	12.692 	--- 	20.24 	6.6 	26.8 	372.3(4.0)
196Hg 	0.2 	30.3(1.0) 	0 	115.(8.) 	0 	115.(8.) 	3080.(180.)
198Hg 	10.1 	--- 	0 	--- 	0 	--- 	2
199Hg 	17 	16.9 	(+/-)15.5 	36.(2.) 	30.(3.) 	66.(2.) 	2150.(48.)
200Hg 	23.1 	--- 	0 	--- 	0 	--- 	<60.
201Hg 	13.2 	--- 	--- 	--- 	--- 	--- 	7.8(2.0)
202Hg 	29.6 	--- 	0 	--- 	0 	9.828 	4.89
204Hg 	6.8 	--- 	0 	--- 	0 	--- 	0.43
Tl 	--- 	8.776 	--- 	9.678 	0.21 	9.89 	3.43
203Tl 	29.524 	6.99 	1.06 	6.14 	0.14 	6.28 	11.4
205Tl 	70.476 	9.52 	-0.242 	11.39 	0.007 	11.4 	0.104
Pb 	--- 	9.405 	--- 	11.115 	0.003 	11.118 	0.171
204Pb 	1.4 	9.90 	0 	12.3 	0 	12.3 	0.65
206Pb 	24.1 	9.22 	0 	10.68 	0 	10.68 	0.03
207Pb 	22.1 	9.28 	0.14 	10.82 	0.002 	10.82 	0.699
208Pb 	52.4 	9.50 	0 	11.34 	0 	11.34 	0.00048
Bi 	100 	8.532 	--- 	9.148 	0.0084 	9.156 	0.0338
Po 	--- 	--- 	0.259 	0 	--- 	--- 	---
At 	--- 	--- 	--- 	0 	--- 	--- 	---
Rn 	--- 	--- 	--- 	0 	--- 	12.6 	---
Fr 	--- 	--- 	--- 	0 	--- 	--- 	---
Ra 	(1.60E3  	10.0(1.0) 	0 	13.(3.) 	0 	13.(3.) 	12.8(1.5)
Ac 	--- 	--- 	--- 	0 	--- 	--- 	---
Th 	100 	10.31 	0 	13.36 	0 	13.36 	7.37
Pa 	(3.28E4  	9.1 	--- 	10.4 	0.1(3.3) 	10.5(3.2) 	200.6(2.3)
U 	--- 	8.417 	--- 	8.903 	0.005 	8.908 	7.57
233U 	(1.59E5  	10.1 	(+/-)1.(3.) 	12.8 	0.1 	12.9 	574.7(1.0)
234U 	0.005 	12.4 	0 	19.3 	0 	19.3 	100.1(1.3)
235U 	0.72 	10.47 	(+/-)1.3 	13.78 	0.2 	14 	680.9(1.1)
238U 	99.275 	8.402 	0 	8.871 	0 	8.871 	2.68
Np 	(2.14E6  	10.55 	--- 	14 	0.5 	14.5 	175.9(2.9)
Pu 	--- 	--- 	--- 	--- 	--- 	--- 	---
238Pu 	(87.74  	14.1 	0 	25.0(1.8) 	0 	25.0(1.8) 	558.(7.)
239Pu 	(2.41E4  	7.7 	(+/-)1.3(1.9) 	7.5 	0.2 	7.7 	1017.3(2.1)
240Pu 	(6.56E3  	3.5 	0 	1.54 	0 	1.54 	289.6(1.4)
242Pu 	(3.76E5  	8.1 	0 	8.2 	0 	8.2 	18.5
Am 	(7.37E3  	8.3 	(+/-)2.(7.) 	8.7 	0.3 	9.0(2.6) 	75.3(1.8)
Cm 	--- 	--- 	--- 	0 	--- 	--- 	---
244Cm 	(18.10  	9.5 	0 	11.3 	0 	11.3 	16.2(1.2)
246Cm 	(4.7E3  	9.3 	0 	10.9 	0 	10.9 	1.36
248Cm 	(3.5E5  	7.7 	0 	7.5 	0 	7.5 	3
EOF
close Fout;
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
