#!/usr/bin/perl

use FileHandle;
use PDL;
use File::Copy;

print "#********************************************************\n";
print "# makenn 220526 - create table with neighbors and interactions\n";
print "# References: M. Rotter et al. PRB 68 (2003) 144418\n";
print "#********************************************************\n";
$PI=3.14159265358979323846;
 $bvkA=25;
 $bvkalpha=0.1;
 # born van karman longitudinal springs:$bvkA*exp(-$bvkalpha*$r*$r);*$r*$r);

unless ($#ARGV>=0) 

{
print STDOUT << "EOF";

 usage: makenn 23.3 [options] [-d]

 meaning take mcphas.j, generate all neighbors within sphere of 23.3A 
 and put them into makenn.j,the output values are sorted by ascending distance

 in interaction columns put by default the classical dipole interaction (meV): this is 
 done assuming the operator sequence 
  I1=Sa I2=La I3=Sb I4=Lb I5=Sc I6=Lc   for sipf files with gJ=0
  I1=Ja I2=Jb I3=Jc    for sipf files with gJ<>0
 (S=spin,L=orbital momentum, J=total angular momentum)

 formula for classical dipole interaction tensor:
 
 Jalphabeta(R)=(mu0/4pi)(gJ muB)^2 (3 Ralpha Rbeta- delta_alphabeta R^2)/R^5

 in the Hamiltonian 

 H= -1/2 sum_ij,alphabeta Jialpha Jalphabeta(Rij) Jjbeta

  Note that in order to use makenn you have to set up a 
 working  mcphas.j file with the crystal structure. 

EOF
print " option -rkky A(meV) kf(1/A) calculates the rkky interaction\n";
print "              according to J(R)=A.cos(2.kf.R)/(2.kf.R)^3\n";
print "              scaling A<0, kf should be the Fermi wavevector (usually\n";
print "              between 0.3-2.5 A^-1 depending on the electrondensity^0.333)\n";
print " option -rkky3d A(meV) ka(1/A) kb(1/A) kc(1/A) calculates the rkky interaction\n";
print "              according to J(R)=A.cos(2.kfR)/(2.kfR)^3\n";
print "              scaling A<0, kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2)\n";
print " option -rkkz A(meV) kf(1/A) calculates the rkky interaction\n";
print "              according to J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4\n";
print "              scaling A>0, kf should be the Fermi wavevector\n";
print " option -rkkz3d A(meV) ka(1/A) kb(1/A) kc(1/A)  calculates the rkky interaction\n";
print "              according to J(R)=A [sin(2.kfR)-2.kfR.cos(2.kfR)]/(2.kfR)^4\n";
print "              scaling A>0, kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2)\n";
print " option -kaneyoshi A(meV) D(A) alpha  calculates the kaneyoshi\n";
print  "             parametrization for the Bethe-Slater\n";
print "              curve: J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2]\n";
print "              with D corresponding to the orbital radius\n";
print "              the exponential alpha is conveniently put to  about 1\n";
print " option -kaneyoshi3d A(meV) Da(A) Db(A) Dc(A) alpha  calculates the 3d-kaneyoshi\n";
print  "             parametrization for the Bethe-Slater\n";
print "              curve: J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2]\n";
print "              with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2)\n";
print "              the exponential alpha is conveniently put to  about 1\n";
print " option -bvk filename\n";
print "              for phonons: take Born van Karman model with longitudinal and\n";
print "              transversal spring constants from file - file format, columns:\n";
print "              #   atom_n_sipf atom_n'_sipf bondlength(A) Clong(N/m) Ctrans(N/m)\n";
print "              mind: into MODPAR2-6 in *.sipf the Einstein-oscillator paramters\n"; 
print "              are written, too. Omit filename to create a sample file with\n";
print "              longitudinal springs:Clong=$bvkA*exp(-$bvkalpha*r/A*r/A) N/m\n";

print STDOUT << "EOF";
 option -cfph 
              calculate crystal field phonon interaction: mcphas.j lists 
              magnetic and non magnetic atoms with charges defined in the 
              sipf files by CHARGE= variable. For magnetic atoms the sipf 
              file the variable MAGNETIC=1 has to be set and information 
              about the ion has to be present (IONTYPE etc.). 
              Foreach magnetic ion a new site is created and shifted
              0.1 A along c in order to not overlap with the original site.
              It is assumed, that the original site will be using an sipf
              file with the MODULE=phonon as well as all the other
              nonmagnetic sites. For the new magnetic site the program
              pointc is used by makenn with option -d to calculate derivatives
              dBlm/du which are inserted as interaction 
              parameters between MODULE=phonon and MODULE=so1ion sites.
              In order to use the resulting file results/makenn.j a phonon
              model has to be set up, the original magnetic atom sites 
              sipffilename has to be changed to the phonon model filename 
              and the phonon model has to be added to makenn.j,  e.g. by
              program addj, moreover magnetic sites sipf files are required, 
              e.g. such as created in results/makenn.a*.sipf. 
 option -e [filename]
 option -f [filename]
 option -dm [filename]
 option -jp [filename]
              read interaction constants from table in file. 
              Use -e for isotropic interactions between momentum Ji and Jj 
                         which only depend on distance (distance vs coupling J)
              Use -f for isotropic interactions between momentum Ji and Jj
                         (neighbour position vs coupling J)
              at positions i and j  
                               ( J   0   0 )
              J Ji.Jj  =    Ji.( 0   J   0 ).Jj  with  H= -1/2 sum_ij J Ji.Jj 
                               ( 0   0   J ) 
              Use -dm for Dzyaloshinski Moriya interactions:
                         (neighbour position vs Dx Dy Dz)
                               ( 0   Dz  -Dy )
              D.(Ji x Jj) = Ji.(-Dz  0    Dx ).Jj  with H= -1/2 sum_ij D.(Ji x Jj) 
                               ( Dy  -Dx  0  )

              Use -jp for jparallel interaction:
                          (neighbour position vs Jp)
                                               ( Rx.Rx   Rx.Ry  Rx.Rz)
              Jp (Ji.R)(Jj.R)/R^2 = Jp/R^2  Ji ( Ry.Rx   Ry.Ry  Ry.Rz) Jj  
                                               ( Rz.Rx   Rz.Ry  Rz.Rz)
 
                           with H= -1/2 sum_ij Jp (Ji.R)(Jj.R)/R^2

              To get a sample file use option -e -f or -dm -jp without a filename.      
EOF
print "        -d puts to the last column the distance of the neighbors (A)\n\n";
print " The neigbours of each atom are also stored in separate files\n";
print " results\/makenn.a*.pc, which can be used with the program pointc to evaluate\n";
print " the pointcharge model and calculate crystal field parameters.\n\n";
 exit 0;}
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($rmax) = eval $ARGV[0];
$rkky=0;$calcdist=0;$bvk=0;$readtable=0;
shift @ARGV; 
$_=$ARGV[0];
if(/-rkky3d/)
  {$rkky=4;shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=$ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$ka=eval $ARGV[0];shift @ARGV;  
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kb=eval $ARGV[0];shift @ARGV;  
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kc=eval $ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A.cos(2.kfR)/(2.kfR)^3 for scale A=$scale meV and\n";
   print "kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n";}
elsif(/-kaneyoshi3d/)
  {$rkky=5;shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=$ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Da=eval $ARGV[0];shift @ARGV;   
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Db=eval $ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Dc=eval $ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$aa=eval $ARGV[0];shift @ARGV;
   print "calculating kaneyoshi parametrization for the Bethe-Slater curve\n";
   print "J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2] for scale A=$scale meV \n";
   print "with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2) with Da=$Da A Db=$Db A Dc=$Dc A  and alpha=$aa\n";}
elsif(/-rkkz3d/)
  {$rkky=6;shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=$ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$ka=eval $ARGV[0];shift @ARGV;  
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kb=eval $ARGV[0];shift @ARGV;  
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kc=eval $ARGV[0];shift @ARGV;   
   print "calculating RKKY interaction J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 for scale A=$scale meV\n";
   print "kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n";}
elsif(/-rkky/)
  {$rkky=1;shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=eval $ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kf=eval $ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A.cos(2.kf.R)/(2.kf.R)^3 for scale A=$scale meV and kf=$kf A^-1\n";}
elsif(/-kaneyoshi/)
  {$rkky=2;shift @ARGV; 
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$scale=eval $ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$D=eval $ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$aa=eval $ARGV[0];shift @ARGV;
   print "calculating kaneyoshi parametrization for the Bethe-Slater curve\n";
   print "J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2] for scale A=$scale meV D=$D A alpha=$aa\n";}
elsif(/-rkkz/)
  {$rkky=3;shift @ARGV; 
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$scale=eval $ARGV[0];shift @ARGV;
   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kf=eval $ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 for scale A=$scale meV and kf=$kf A^-1\n";}
elsif(/-bvk/)
  {$bvk=1;shift @ARGV;
   unless ($#ARGV>=0) # if no filename is given print help
   { $tabout=1; $readtable=$_;}
  else
   {
   $bfk_file=$ARGV[0];shift @ARGV;
   print "creating phononic interactions from Born van Karman model in file $bfk_file\n";
   # read interaction constants from file
   unless(open(Fin,$bfk_file)){ die "could not open $bfk_file\n";}
             $nof_springs=0;
             {while($line=<Fin>){next if $line=~/^\s*#/;
                                 my @numbers=split(" ",$line);
                                 if($#numbers>=4)
                                 {++$nof_springs;print $line;
                                  $atom_n[$nof_springs]=$numbers[0];
                                  $atom_m[$nof_springs]=$numbers[1];
                                  $bondlength[$nof_springs]=$numbers[2];
                                  $long_spring[$nof_springs]=$numbers[3];
                                  $trans_spring[$nof_springs]=$numbers[4];
                                 }
                                }
              close Fin;
       	     }
  }
  }
elsif(/-cfph/)
{$cfph=1; shift @ARGV;
 print "creating crystal field phonon interactions from pointcharge model using program pointc\n";
}
elsif(/-e/||/-f/||/-dm/||/-jp/)
  {$readtable=$_;shift @ARGV;
   unless ($#ARGV>=0) # if no filename is given print help
   {
 $tabout=1; }
  else
 {
   $table_file=$ARGV[0];shift @ARGV;
   $ignore_neihgbours_behind=0;$rtab=$readtable;
   print "reading "; if ($readtable=~/-e/) {print " isotropic ";$DM=0;$Jp=0;} 
                  elsif ($readtable=~/-f/) {print " isotropic ";$DM=0;$Jp=0;} 
                  elsif ($readtable=~/-jp/) {print " Jp ";$DM=0;$Jp=1;} 
                   else {print " DM ";$DM=1;$Jp=0;}
   $readtable=1;
   print" interactions from table in file $table_file\n";
   # read interaction constants from file
   unless(open(Fin,$table_file)){ die "could not open $table_file\n";}
             $n_table=0;
             while(<Fin>){      
               if (/^(#!|[^#])*nr1min\s*=\s*/){$readtable=2;($n1min)=extract("nr1min",$_);}
               if (/^(#!|[^#])*nr2min\s*=\s*/){$readtable=2;($n2min)=extract("nr2min",$_);}
               if (/^(#!|[^#])*nr3min\s*=\s*/){$readtable=2;($n3min)=extract("nr3min",$_);}
               if (/^(#!|[^#])*nr1max\s*=\s*/){$readtable=2;($n1max)=extract("nr1max",$_);}
               if (/^(#!|[^#])*nr2max\s*=\s*/){$readtable=2;($n2max)=extract("nr2max",$_);}
               if (/^(#!|[^#])*nr3max\s*=\s*/){$readtable=2;($n3max)=extract("nr3max",$_);}
               if (/^(#!|[^#])*ignore_neihgbours_behind\s*=\s*/){($ignore_neihgbours_behind)=extract("ignore_neihgbours_behind",$_);}
                                 next if /^\s*#/;$line=$_;
                                 my @numbers=split(" ",$line);
                                 if($rtab=~/-e/||$rtab=~/-jp/)
                                 {++$n_table;
                                  $da[$n_table]=$numbers[0];
                                  $Jex[$n_table]=$numbers[1];
                                 }
                                 if($rtab=~/-f/||$rtab=~/-dm/)
                                 {++$n_table;
                                  $da[$n_table]=$numbers[0];
                                  $db[$n_table]=$numbers[1];
                                  $dc[$n_table]=$numbers[2];
                                  $Jex[$n_table]=$numbers[3];
                                  if($DM==1){
                                  $Jey[$n_table]=$numbers[4];
                                  $Jez[$n_table]=$numbers[5];
                                            }
                                 }
                                }
              close Fin;
       	     }
  }
 
$_=$ARGV[0];
if(/-d/)
  {$calcdist=1;print "putting distance of neighbors (A) to last column of makenn.j\n";}

my ($latt,$p) = getlattice("./mcphas.j"); # gets lattice and atomic positions
my ($a,$b,$c,$alpha,$beta,$gamma,$nofatoms,$nofcomponents) = @{$latt};
unless($tabout) {
 print "rmax=".$rmax." A\n";
 print "a=".$a." b=".$b." c=".$c." alpha=".$alpha." beta=".$beta." gamma=".$gamma."\n";
 print "primitive lattice[abc]:".$p."\n";
}
# define transformation matrix to calculate components of
# r1 r2 and r3 with respect to the ijk coordinate system
# defined by j||b, k||(a x b) and i normal to k and j
$rtoijk = pdl [ [$a*sin($gamma*$PI/180),$a*cos($gamma*$PI/180),0],   # vector a in ijk coordinates
                [0,$b, 0],                                           # vector b in ijk coordinates
                [$c*(cos($beta*$PI/180)-cos($alpha*$PI/180)*cos($gamma*$PI/180))/sin($gamma*$PI/180),$c*cos($alpha*$PI/180),0]
                                                                    # not yet vector c in ijk coordinates
                  ];
if (abs($rtoijk->at(2,0))>$c){die "ERROR makenn: alpha beta and gamma geometrically inconsistent\n";}
$t=$rtoijk->slice("2:2,2:2");
$t.=pdl[[$c*$c-$rtoijk->at(0,2)*$rtoijk->at(0,2)-$rtoijk->at(1,2)*$rtoijk->at(1,2)]];
if ($rtoijk->at(2,2)<=0){die "ERROR makenn: alpha beta and gamma geometrically inconsistent\n";}
$t.=sqrt($t);
#print $t;
#print $rtoijk;
#$invrtoijk=matinv($rtoijk); #invert this matrix for use later
$invrtoijk=inv($rtoijk); #invert this matrix for use later
#print $invrtoijk;
#'x' is hijacked as the matrix multiplication operator. e.g. $c = $a x $b;#
#
#perlDL is row-major not column major so this is actually c(i,j) = sum_k a(k,j) b(i,k) - but when matrices are printed the results will look right. Just remember the indices are reversed. e.g.:#
#
# $a = [                   $b = [
#       [ 1  2  3  0]            [1 1]
#       [ 1 -1  2  7]            [0 2]
#       [ 1  0  0  1]            [0 2]
#      ]                         [1 1]
#                               ]#
#
# gives $c = [
#             [ 1 11]
#             [ 8 10]
#             [ 2  2]
#            ]Note: transpose() does what it says and is a convenient way to turn row vectors into column vectors.

$p= $p x $rtoijk;
unless($tabout){ print "primitive lattice[A]:".$p."\n";}
    $r=0;

# first determine maximum distance of a basis atom to origin of unit cell
    $distmax=0;
  for ($nz=1;$nz<=$nofatoms;++$nz){
   $dabc=pdl [($x[$nz]),($y[$nz]),($z[$nz])];
   $rvec= $dabc x $rtoijk;$rvec=$rvec->slice(":,(0)");
   $rr=inner($rvec, $rvec);
   $r=sqrt($rr->at());
   if($r>$distmax){$distmax=$r;}
  }
  unless($readtable>1)
 {
# determine $nmin,$nmax by looking at a cube with side 3rmax
#     $inv=matinv(transpose($p)); #invert primitive lattice
     $inv=inv(transpose($p)); #invert primitive lattice
# print "inverted primitive lattice[A]:".$inv."\n";
     #loop all corner points
  for ($i1=-1;$i1<=1;$i1+=2){
  for ($i2=-1;$i2<=1;$i2+=2){
  for ($i3=-1;$i3<=1;$i3+=2){
    $n=inner($inv , pdl[$i1*($rmax+$distmax)*1.5,$i2*($rmax+$distmax)*1.5,$i3*($rmax+$distmax)*1.5]);
    if (($n->at(0))<$n1min) {$n1min=int($n->at(0))-1;}
    if (($n->at(1))<$n2min) {$n2min=int($n->at(1))-1;}
    if (($n->at(2))<$n3min) {$n3min=int($n->at(2))-1;}
    if (($n->at(0))>$n1max) {$n1max=int($n->at(0))+1;}
    if (($n->at(1))>$n2max) {$n2max=int($n->at(1))+1;}
    if (($n->at(2))>$n3max) {$n3max=int($n->at(2))+1;}
#   print"corner $i1 $i2 $i3 coordinates in prim bases:$n\n";
#   print "$n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";
  }}}
  }
if($tabout)
 {if($bvk)
  {print STDOUT << "EOF";
#  atom_n_sipf atom_n'_sipf bondlength(A) Clong(N/m) Ctrans(N/m)  generated with Clong=$bvkA*exp(-$bvkalpha*r/A*r/A) N/m                 
EOF
  }
else {    print STDOUT << "EOF";
# Sample table of interaction constants for program makenn
#
# primitive unit cell range to probe for neighbours (optional):
#! nr1min=$n1min  nr2min=$n2min  nr3min=$n3min 
#! nr1max=$n1max  nr2max=$n2max  nr3max=$n3max
#
#! ignore_neihgbours_behind=0   # if set to 1 then only nearest neighbour in each 
#                               # direction is taken
#  
EOF
$_=$readtable;
if(/-e/)
{print STDOUT << "EOF";
# Table of exchange interaction constants J - assumed to be isotropic according to in H= -1/2 sum_ij J Ji.Jj 
# |Rij| [A]  J [meV]  atom_i  atom_j 
EOF

}
elsif(/-jp/)
{print STDOUT << "EOF";
# Table of exchange interaction constants Jp  in H= -1/2 sum_ij Jp (Ji.R)(Jj.R)/R^2 
# |Rij| [A]  Jp [meV]  atom_i  atom_j 
EOF

}
elsif(/-f/)
{print STDOUT << "EOF";
# Table of exchange interaction constants J - assumed to be isotropic according to in H= -1/2 sum_ij J Ji.Jj 
# da [a]    db [b]    dc [c]    J [meV]   atom_i  atom_j  |Rij| [A] 
EOF

}
else
{print STDOUT << "EOF";
# Table of Dzyaloshinska Moriya exchange interaction constants DM  in H= -1/2 sum_ij D.(Ji x Jj) 
#  (default values are neighbour pos in [A]- please modify)
# da [a]    db [b]    dc [c]    Dx [meV] Dy [meV] Dz [meV] atom_i  atom_j  |Rij| [A] 
EOF

}
}}
else
{
print "# $n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";

     # initialize output file results/makenn.j
  ($h,$l)=printlattice("./mcphas.j",">./results/makenn.j");
print "# number of atoms = $nofatoms\n calculating ...\n";
}               
 for ($nnn=1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)    
 {unless($tabout){ if($nnn>$nofatoms){print "# new magnetic ";}
   print "# atom $nnn ...";}
     my $gJ=$gJ[$nnn];
     my $sipffilename=$sipf_file[$nnn];
     my ($rn)=new PDL ();
     my ($an)=new PDL ();

     my ($xn)=new PDL();    
     my ($yn)=new PDL();    
     my ($zn)=new PDL();    

     my ($in)=new PDL();
     my ($jn)=new PDL();
     my ($kn)=new PDL();

     my ($Jaa)=new PDL();   
     my ($Jab)=new PDL();    
     my ($Jac)=new PDL();    

     my ($Jba)=new PDL();    
     my ($Jbb)=new PDL();    
     my ($Jbc)=new PDL();    

     my ($Jca)=new PDL();    
     my ($Jcb)=new PDL();    
     my ($Jcc)=new PDL(); 
     
    
  for ($n1=$n1min;$n1<=$n1max;++$n1){ 
  for ($n2=$n2min;$n2<=$n2max;++$n2){ 
  for ($n3=$n3min;$n3<=$n3max;++$n3){  
  for ($nz=1;$nz<=$nofatoms;++$nz){  
   $dabc=pdl [($x[$nz]-$x[$nnn]),($y[$nz]-$y[$nnn]),($z[$nz]-$z[$nnn])];
   $rvec= $dabc x $rtoijk;$rvec=$rvec->slice(":,(0)");
   $rvec+=$n1*$p->slice(",(0)")+$n2*$p->slice(",(1)")+$n3*$p->slice(",(2)");

   $rr=inner($rvec, $rvec);
   $r=sqrt($rr->at());
   $aabbcc=$rvec x $invrtoijk;$aabbcc=$aabbcc->slice(":,(0)");
   $xx=$aabbcc->at(0);
   $yy=$aabbcc->at(1);
   $zz=$aabbcc->at(2);

   if ($r<=$rmax && $r>0){#save neighbour j format
   if($tabout){
# check if neighbour is already in table
$ff=0;
 for($ntbl=1;$ntbl<=$n_table;++$ntbl){$_=$readtable;
if((/-e/||/-jp/)&&abs($da[$ntbl]-$r)<0.001){$ff=1;}
elsif((/-f/||/-dm/)&&
      abs($da[$ntbl]-$xx)<0.001&&
      abs($db[$ntbl]-$yy)<0.001&&
      abs($dc[$ntbl]-$zz)<0.001){$ff=1;}
elsif((/-bvk/)&&(abs($da[$ntbl]-$r)<0.001)){
       $a1=$sipf_file[$nnn];$a2=$sipf_file[$nz];
       if(($db[$ntbl]=~/$a1/&&$dc[$ntbl]=~/$a2/)
       ||($db[$ntbl]=~/$a2/&&$dc[$ntbl]=~/$a1/)){$ff=1;}}
                                  }

if($ff==0){++$n_table;$ntbl=$n_table;
          $_=$readtable;
if(/-e/||/-jp/)
{$da[$ntbl]=$r;
 print sprintf("%+10.6f     0       a%i a%i \n",$r,$nnn,$nz);
}
elsif(/-f/)
{$da[$ntbl]=$xx;$db[$ntbl]=$yy;$dc[$ntbl]=$zz;
 print sprintf("%+10.6f %+10.6f %+10.6f 0        a%i a%i %+10.6f\n",$xx, $yy ,$zz,$nnn,$nz,$r);
}
elsif(/-dm/)
{$da[$ntbl]=$xx;$db[$ntbl]=$yy;$dc[$ntbl]=$zz;
  print sprintf("%+10.6f %+10.6f %+10.6f    %+10.6f %+10.6f %+10.6f     a%i a%i %+10.6f\n",$xx, $yy ,$zz,$rvec->at(0),$rvec->at(1),$rvec->at(2),$nnn,$nz,$r);
}
elsif(/-bvk/)
{$da[$ntbl]=$r; $db[$ntbl]=$sipf_file[$nnn];$dc[$ntbl]=$sipf_file[$nz];
 $spring= $bvkA*exp(-$bvkalpha*$r*$r);
 print sprintf("%s   %s    %+10.6f   %+10.6f   0 \n",$sipf_file[$nnn],$sipf_file[$nz],$r,$spring);
}   
else{die "Error makenn - creating table for option $_ \n";}
          }
              }

          unless($readtable>0){
    $an=$an->append( pdl ([$nz]));
    $rn=$rn->append( pdl ([$r]));
    $xn=$xn->append( pdl ([$xx]));
    $yn=$yn->append( pdl ([$yy]));
    $zn=$zn->append( pdl ([$zz]));
    $in=$in->append( pdl ([$rvec->at(0)]));
    $jn=$jn->append( pdl ([$rvec->at(1)]));
    $kn=$kn->append( pdl ([$rvec->at(2)]));

    my ($interaction) = getinteraction($gJ,$gJ[$nz],$sipffilename,$sipf_file[$nz],$r,$rvec->at(0),$rvec->at(1),$rvec->at(2));
    my ($jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc) = @{$interaction};
    $Jaa=$Jaa->append( pdl ([$jaa]));
    $Jab=$Jab->append( pdl ([$jab]));
    $Jac=$Jac->append( pdl ([$jac]));
    $Jba=$Jba->append( pdl ([$jba]));
    $Jbb=$Jbb->append( pdl ([$jbb]));
    $Jbc=$Jbc->append( pdl ([$jbc]));
    $Jca=$Jca->append( pdl ([$jca]));
    $Jcb=$Jcb->append( pdl ([$jcb]));
    $Jcc=$Jcc->append( pdl ([$jcc]));
                               } else
                               {# check if neighbour is in readtable - if yes, save it
                               for($ntbl=1;$ntbl<=$n_table;++$ntbl){
                                 if(($rtab=~/-f/&&
                                    abs($da[$ntbl]-$xx)<0.001&&
                                    abs($db[$ntbl]-$yy)<0.001&&
                                    abs($dc[$ntbl]-$zz)<0.001&&
                                    $Jex[$ntbl]!=0.0)||
                                    (($rtab=~/-e/||$rtab=~/-jp/)&&
                                     abs($da[$ntbl]-$r)<0.001&&
                                     $Jex[$ntbl]!=0.0)||
                                   ($rtab=~/-dm/&&
                                    abs($da[$ntbl]-$xx)<0.001&&
                                    abs($db[$ntbl]-$yy)<0.001&&
                                    abs($dc[$ntbl]-$zz)<0.001&&
                                    ($Jex[$ntbl]!=0.0||
                                    $Jey[$ntbl]!=0.0)||
                                    $Jez[$ntbl]!=0.0)
                                   ) 
                                  {# ok save it
    $an=$an->append( pdl ([$nz]));
    $rn=$rn->append( pdl ([$r]));
    $xn=$xn->append( pdl ([$xx]));
    $yn=$yn->append( pdl ([$yy]));
    $zn=$zn->append( pdl ([$zz]));
    $in=$in->append( pdl ([$rvec->at(0)]));
    $jn=$jn->append( pdl ([$rvec->at(1)]));
    $kn=$kn->append( pdl ([$rvec->at(2)]));
unless($DM>0){
 unless($Jp>0){
    $Jaa=$Jaa->append( pdl ([$Jex[$ntbl]]));
    $Jab=$Jab->append( pdl ([0.0]));
    $Jac=$Jac->append( pdl ([0.0]));
    $Jba=$Jba->append( pdl ([0.0]));
    $Jbb=$Jbb->append( pdl ([$Jex[$ntbl]]));
    $Jbc=$Jbc->append( pdl ([0.0]));
    $Jca=$Jca->append( pdl ([0.0]));
    $Jcb=$Jcb->append( pdl ([0.0]));
    $Jcc=$Jcc->append( pdl ([$Jex[$ntbl]]));              
              }
               else  # option $Jp=1
               {$rhx=$rvec->at(0)/$r;$rhy=$rvec->at(1)/$r;$rhz=$rvec->at(2)/$r;
    $Jaa=$Jaa->append( pdl ([$Jex[$ntbl]*$rhx*$rhx]));
    $Jab=$Jab->append( pdl ([$Jex[$ntbl]*$rhx*$rhy]));
    $Jac=$Jac->append( pdl ([$Jex[$ntbl]*$rhx*$rhz]));
    $Jba=$Jba->append( pdl ([$Jex[$ntbl]*$rhy*$rhx]));
    $Jbb=$Jbb->append( pdl ([$Jex[$ntbl]*$rhy*$rhy]));
    $Jbc=$Jbc->append( pdl ([$Jex[$ntbl]*$rhy*$rhz]));
    $Jca=$Jca->append( pdl ([$Jex[$ntbl]*$rhz*$rhx]));
    $Jcb=$Jcb->append( pdl ([$Jex[$ntbl]*$rhz*$rhy]));
    $Jcc=$Jcc->append( pdl ([$Jex[$ntbl]*$rhz*$rhz]));
   

               }
              } else
              {
    $Jaa=$Jaa->append( pdl ([0.0]));
    $Jab=$Jab->append( pdl ([$Jez[$ntbl]]));
    $Jac=$Jac->append( pdl ([-$Jey[$ntbl]]));
    $Jba=$Jba->append( pdl ([-$Jez[$ntbl]]));
    $Jbb=$Jbb->append( pdl ([0.0]));
    $Jbc=$Jbc->append( pdl ([$Jex[$ntbl]]));
    $Jca=$Jca->append( pdl ([$Jey[$ntbl]]));
    $Jcb=$Jcb->append( pdl ([-$Jex[$ntbl]]));
    $Jcc=$Jcc->append( pdl ([0.0]));
              }
                                  }    

                         } 
                               }
                      }
    }}}}  

   $n= qsorti($rn); 
   $nofneighbours[$nnn]=(($rn->dims)[0]-1); 
   if ($readtable>0&&$ignore_neihgbours_behind==1)
   {# check if there is a closer neighbour and if there is delete neighbour
    for ($n1=2;$n1<(($rn->dims)[0]);++$n1)
     {for ($n2=1;$n2<$n1;++$n2){  $innr=$xn->index($n)->at($n1)*$xn->index($n)->at($n2)+
                                   $yn->index($n)->at($n1)*$yn->index($n)->at($n2)+
                                   $zn->index($n)->at($n1)*$zn->index($n)->at($n2);
                                  $rrn1=$xn->index($n)->at($n1)*$xn->index($n)->at($n1)+
                                   $yn->index($n)->at($n1)*$yn->index($n)->at($n1)+
                                   $zn->index($n)->at($n1)*$zn->index($n)->at($n1);
                                  $rrn2=$xn->index($n)->at($n2)*$xn->index($n)->at($n2)+
                                   $yn->index($n)->at($n2)*$yn->index($n)->at($n2)+
                                   $zn->index($n)->at($n2)*$zn->index($n)->at($n2);
                                  
                                if(abs($innr*$innr-$rrn1*$rrn2)<0.0001&&$innr>0)
                                  {# remove neighbour n1
                                   set $rn->index($n),$n1,0;--$nofneighbours[$nnn];
                                   $n2=$n1;
                                  } 
                               }
     }
   }
   if($nofneighbours[$nnn]==-1){$nofneighbours[$nnn]=0;}
   unless($tabout){print $nofneighbours[$nnn]." neighbours found\n";}
   $n= qsorti($rn); 
unless($tabout){printneighbourlist($h,$l,$nofneighbours[$nnn],$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb);}
 }



 endprint($h,$l);   
if($cfph==1){
# for cf phonon interaction recreate makenn.j
my ($h,$l)=printlattice("./mcphas.j",">./results/makenn.j");
for($nnn=$nofatoms+1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)
 { # run pointc to create derivatives of Blm from pointcharge model for magnetic atoms
  system("pointc -d ".$sipf_file[$nnn]." results/makenn.a$nnn.pc > results/makenn.a$nnn.sipf");
  $sipf_file[$nnn]="results/makenn.a$nnn.sipf";
  copy("results/pointc.dBlm","results/makenn.a$nnn.dBlm");
 }
for($nnn=1;$nnn<=$nofatoms;++$nnn){$nofneighbours[$nnn]=0;}
for($nnn=$nofatoms+1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)
 {

# $nph[$nnn] contains the index of the PHONON atom corresponing to the magnetic atom $nnn
# - load the indices of the neighbours from makenn.a$nnn.pc and combine these into a field
# ... using this information push outstring in the next loop into a piddle
# such that it can be output subsequently 
unless (open(Fin,"results/makenn.a$nnn.pc")){die "cannot open file results/makenn.a$nnn.pc\n";}
  $nofn=0;
  while(<Fin>){next if /^\s*#/;
                   my @numbers=split(" ",$_); 
                   if($#numbers>=4)
                    {++$nofn; $ni[$nofn]=$numbers[8];
   }}
close Fin;
  ++$nofn;$ni[$nofn]=$nph[$nnn];$ni[0]=$nofn;
  $nofn=0;
  unless (open(Fin,"results/makenn.a$nnn.dBlm")){die "cannot open file results/makenn.a$nnn.dBlm\n";}
  while(<Fin>){next if /^\s*#/;
                   my @numbers=split(" ",$_); 
                   if($#numbers>=4)
                    {++$nofn;

   $rvec=pdl [$numbers[0],$numbers[1],$numbers[2]];
   $aabbcc=$rvec x $invrtoijk;$aabbcc=$aabbcc->slice(":,(0)");
   $da=$aabbcc->at(0);
   $db=$aabbcc->at(1);
   $dc=$aabbcc->at(2);
  $outstring=sprintf("%+10.6f %+10.6f %+10.6f ",$da, $db ,($dc-0.1/$c));
  $outstringp=sprintf("%+10.6f %+10.6f %+10.6f ",-$da, -$db ,-$dc+0.1/$c);
  # now the off diagonal Elements: Jab Jba Jac Jca Jad Jda .. Jbc Jcb Jbd Jdb ... Jcd Jdc ... etc.
  # a b c d e ... = 11/x 11s/y 10/z  22s 21s 20 21 22  33s 32s ...
  #                 I1    I2    I3    I4  I5
  # for PHONON      ux    uy    uz    dum dum dum ...
  # for so1ion      O11  O11s   O10   O22s O21s ...
  # Hamiltonian Hcfph= - sum_i<j,lmgamma u_gamma(i) (-)dBlm/du_gamma  Olm(j) = 
  #                  =-sumi<j_alphabeta  Ialpha Jalphabeta Ibeta
  #                  =-1/2sumij_alphabeta  Ialpha Jalphabeta Ibeta
  # note: the factor 1/2 comes from the fact, that the pairs are counted twice in the last expression 
  # i.e. J12=0 J13=0 J41=-2 dB22s/dux J14=0 J51=-2 dB21s/dux J15=0 ...
  # 5+9+13=27 ... 27x3=81  5+81=86
  for($i=6;$i<=86;++$i){$outstring.= " ".(-$numbers[$i]);$outstringp.= " ".(-$numbers[$i]);}
  $outstring.= "\n";
  $outstringp.= "\n";
   $ph[$nnn].=$outstring;
  # push the outstring to the appropriate PHONON ions neighbour list
  $ph[$ni[$nofn]].=$outstringp;
  ++$nofneighbours[$ni[$nofn]];
                    }
                   }
  close Fin;
 }
for($nnn=1;$nnn<=$nofatoms;++$nnn)
 {print $l ("#*************************************************************************\n");
  print $l ("#! da=".$x[$nnn]." [a] db=".$y[$nnn]." [b] dc=".$z[$nnn]." nofneighbours=".$nofneighbours[$nnn]." diagonalexchange=2  sipffilename=".$sipf_file[$nnn]."\n");
  print $l ("# crystal field phonon interaction parameters from pointcharge calculation\n");
  print $l ("# da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
  print $l ("#! symmetricexchange=0 indexexchange= 1,4 2,4 3,4 1,5 2,5 3,5 1,6 2,6 3,6 1,7 2,7 3,7 1,8 2,8 3,8 ");
  #O4m
  print $l (" 1,16 2,16 3,16 1,17 2,17 3,17 1,18 2,18 3,18 1,19 2,19 3,19 1,20 2,20 3,20 1,21 2,21 3,21 1,22 2,22 3,22 1,23 2,23 3,23 1,24 2,24 3,24 ");
  #O6m
  print $l (" 1,36 2,36 3,36 1,37 2,37 3,37 1,38 2,38 3,38 1,39 2,39 3,39 1,40 2,40 3,40 1,41 2,41 3,41 1,42 2,42 3,42 1,43 2,43 3,43 1,44 2,44 3,44 1,45 2,45 3,45 1,46 2,46 3,46 1,47 2,47 3,47 1,48 2,48 3,48 \n");
# here we should enter the cf- phonon interactions for the MODULE=phonon oscillators   
#  1) calculate  number of neighbours (only the magnetic atoms) 
#  2) fill values from results/makenn.a$nnn.dBlm  [all done above for magnetic atoms and filled in here]
  print $l $ph[$nnn];
 }
for($nnn=$nofatoms+1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)
 {print $l ("#*************************************************************************\n");
  print $l ("#! da=".$x[$nnn]." [a] db=".$y[$nnn]." [b] dc=".($z[$nnn]+0.1/$c)." nofneighbours=".($nofneighbours[$nnn]+1)." diagonalexchange=2  sipffilename=".$sipf_file[$nnn]."\n");
  print $l ("# crystal field phonon interaction parameters from pointcharge calculation\n");
  print $l ("# da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
  print $l ("#! symmetricexchange=0 indexexchange= 4,1 4,2 4,3 5,1 5,2 5,3 6,1 6,2 6,3 7,1 7,2 7,3 8,1 8,2 8,3 ");
  # O4m
  print $l (" 16,1 16,2 16,3 17,1 17,2 17,3 18,1 18,2 18,3 19,1 19,2 19,3 20,1 20,2 20,3 21,1 21,2 21,3 22,1 22,2 22,3 23,1 23,2 23,3 24,1 24,2 24,3 ");
  # O6m
  print $l (" 36,1 36,2 36,3 37,1 37,2 37,3 38,1 38,2 38,3 39,1 39,2 39,3 40,1 40,2 40,3 41,1 41,2 41,3 42,1 42,2 42,3 43,1 43,2 43,3 44,1 44,2 44,3 45,1 45,2 45,3 46,1 46,2 46,3 47,1 47,2 47,3 48,1 48,2 48,3 \n"); 
  # here come the cf-phonon interactions for the MODULE=so1ion magnetic ions 
  # open pointcharge file 
  print $l $ph[$nnn];
}
endprint($h,$l);
}

if($tabout){exit;}
if($bvk){ # check consistency of files
foreach(@sipf_file) {
unless(/results/){$sf=$_;
 foreach(@sipf_file) {
   if(/results\/$sf\.\d/){ # compare MODPAR1-7
$sfc=$_;
$M[2]=extractfromfile("MODPAR2",$sf);
$M[3]=extractfromfile("MODPAR3",$sf);
$M[4]=extractfromfile("MODPAR4",$sf);
$M[5]=extractfromfile("MODPAR5",$sf);
$M[6]=extractfromfile("MODPAR6",$sf);
$M[7]=extractfromfile("MODPAR7",$sf);

$Mc[2]=extractfromfile("MODPAR2",$sfc);
$Mc[3]=extractfromfile("MODPAR3",$sfc);
$Mc[4]=extractfromfile("MODPAR4",$sfc);
$Mc[5]=extractfromfile("MODPAR5",$sfc);
$Mc[6]=extractfromfile("MODPAR6",$sfc);
$Mc[7]=extractfromfile("MODPAR7",$sfc);
for($i=2;$i<=7;++$i){ if(abs($Mc[$i]-$M[$i])>1e-4){die("ERROR makenn: input file mcphas.j inconsistent ! MODPAR$i in $sfc (".$Mc[$i].")  does not match $sf (".$M[$i].") - sites with same sipf file $sf are not equivalent \n");}
                     }          }
                 }
                        } 
                     }
}

print "created files: results/makenn.j     (interaction parameters)\n";
print "               results/makenn.a*.pc (pointcharge environment files)\n";
print "********************************************************\n";
print "               end of program makenn\n";
print " Reference: M. Rotter et al. PRB 68 (2003) 144418\n";
print "********************************************************\n";

   exit;

#-----------------------------------------------------------------------
sub clearMP{
my ($file)=@_;
my $i,$j;$i=0;
             if(open (Fin,$file))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\bMODPAR2\s*=\s*/) {($Kxxread)=($line=~m/^(?:#!|[^#])*\bMODPAR2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxxread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR2\s*=\s*[^\s\;\n\t\*]+/$1MODPAR2=$Kxxread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR3\s*=\s*/) {($Kyyread)=($line=~m/^(?:#!|[^#])*\bMODPAR3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kyyread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR3\s*=\s*[^\s\;\n\t\*]+/$1MODPAR3=$Kyyread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR4\s*=\s*/) {($Kzzread)=($line=~m/^(?:#!|[^#])*\bMODPAR4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kzzread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR4\s*=\s*[^\s\;\n\t\*]+/$1MODPAR4=$Kzzread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR5\s*=\s*/) {($Kxyread)=($line=~m/^(?:#!|[^#])*\bMODPAR5\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);                                                           $Kxyread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR5\s*=\s*[^\s\;\n\t\*]+/$1MODPAR5=$Kxyread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR6\s*=\s*/) {($Kyzread)=($line=~m/^(?:#!|[^#])*\bMODPAR6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kyzread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR6\s*=\s*[^\s\;\n\t\*]+/$1MODPAR6=$Kyzread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR7\s*=\s*/) {($Kxzread)=($line=~m/^(?:#!|[^#])*\bMODPAR7\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxzread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR7\s*=\s*[^\s\;\n\t\*]+/$1MODPAR7=$Kxzread/g;}
               $line_store[$i]=$line;++$i;}
              close Fin;
              open(Fout,">$file");
              for($j=0;$j<=$i;++$j){print Fout $line_store[$j];}
              close Fout;              
       	     }
             else
             {
             print STDERR "Warning: failed to read/write MODPARS from/to data file \"$filename\"\n";
             }


}


sub addK{
my ($file)=@_;
my $i,$j;$i=0;
             if(open (Fin,$file))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\bMODPAR2\s*=\s*/) {($Kxxread)=($line=~m/^(?:#!|[^#])*\bMODPAR2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxxread+=$Kxx/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR2\s*=\s*[^\s\;\n\t\*]+/$1MODPAR2=$Kxxread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR3\s*=\s*/) {($Kyyread)=($line=~m/^(?:#!|[^#])*\bMODPAR3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kyyread+=$Kyy/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR3\s*=\s*[^\s\;\n\t\*]+/$1MODPAR3=$Kyyread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR4\s*=\s*/) {($Kzzread)=($line=~m/^(?:#!|[^#])*\bMODPAR4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kzzread+=$Kzz/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR4\s*=\s*[^\s\;\n\t\*]+/$1MODPAR4=$Kzzread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR5\s*=\s*/) {($Kxyread)=($line=~m/^(?:#!|[^#])*\bMODPAR5\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxyread+=$Kxy/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR5\s*=\s*[^\s\;\n\t\*]+/$1MODPAR5=$Kxyread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR6\s*=\s*/) {($Kyzread)=($line=~m/^(?:#!|[^#])*\bMODPAR6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kyzread+=$Kyz/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR6\s*=\s*[^\s\;\n\t\*]+/$1MODPAR6=$Kyzread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR7\s*=\s*/) {($Kxzread)=($line=~m/^(?:#!|[^#])*\bMODPAR7\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxzread+=$Kxz/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR7\s*=\s*[^\s\;\n\t\*]+/$1MODPAR7=$Kxzread/g;}
               $line_store[$i]=$line;++$i;}
              close Fin;
              open(Fout,">$file");
              for($j=0;$j<=$i;++$j){print Fout $line_store[$j];}
              close Fout;              
       	     }
             else
             {
             print STDERR "Warning: failed to read/write MODPARS from/to data file \"$filename\"\n";
             }


}

sub getinteraction {
   my ($gJthis,$gJ,$sipffilethis,$sipffile,$r,$rx,$ry,$rz)=@_;
   my $n;

 if($bvk==1)
 { # here do the Born von Karman calculation using spring constants
  $a0 = .5292e-10;#(m)
  $J2meV=1/1.60217646e-22; #1 millielectron volt = 1.60217646 . 10-22 joules
  $jaa=0;$jbb=0;$jcc=0;$jab=0;$jbc=0;$jac=0;
  for($n=1;$n<=$nof_springs;++$n){$a1=$atom_m[$n];$a2=$atom_n[$n];
         if(abs($r-$bondlength[$n])<0.01
            &&(($sipffilethis=~/(results\/)?$a2/&&$sipffile=~/(results\/)?$a1/)
             ||($sipffilethis=~/(results\/)?$a1/&&$sipffile=~/(results\/)?$a2/)
              )
            )
             {# bond found - do something
  $cL=$long_spring[$n]*$a0*$a0*$J2meV;
  $cT=$trans_spring[$n]*$a0*$a0*$J2meV; 
  $Kxx= -(( $cL-$cT)* $rx * $rx + $cT* $r *$r)  /$r /$r;$jaa-=$Kxx;
  $Kyy= -(( $cL-$cT)* $ry * $ry + $cT* $r *$r)  /$r /$r;$jbb-=$Kyy;
  $Kzz= -(( $cL-$cT)* $rz * $rz + $cT* $r *$r)  /$r /$r;$jcc-=$Kzz;

  $Kxy= -($cL-$cT) * ( $rx * $ry) /$r /$r ;$jab-=$Kxy;
  $Kyz= -($cL-$cT) * ( $ry * $rz) /$r /$r ;$jbc-=$Kyz;
  $Kxz= -($cL-$cT) * ( $rx * $rz) /$r /$r ;$jac-=$Kxz;
   # add also something to the Knn and Kmm in $sipffilethis and $sipffile !!! 
   addK($sipffilethis);
   addK($sipffile);
             }
                                 }
   $jba=$jab;   $jcb=$jbc;   $jca=$jac;
 } else{
   if ($gJthis==0){$gJthis=1;} # set gJ=1 for ions with intermediate coupling
   if ($gJ==0){$gJ=1;}

   # calculate classical dipole interaction
  if ($rkky==1)
  {$jaa = $scale*cos(2*$kf*$r)/8/$kf/$kf/$kf/$r/$r/$r;
  $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0; 
  }

  elsif ($rkky==4)
  {$kfr=sqrt($ka*$ka*$rx*$rx+$kb*$kb*$ry*$ry+$kc*$kc*$rz*$rz);
  $jaa = $scale*cos(2*$kfr)/8/$kfr/$kfr/$kfr;
  $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
  }

  elsif ($rkky==3)

  {$jaa = $scale*(sin(2*$kf*$r)-2*$kf*$r*cos(2*$kf*$r))/16/$kf/$kf/$kf/$kf/$r/$r/$r/$r;
  $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
  }

  elsif ($rkky==6)
  {$kfr=sqrt($ka*$ka*$rx*$rx+$kb*$kb*$ry*$ry+$kc*$kc*$rz*$rz);
  $jaa = $scale*(sin(2*$kfr)-2*$kfr*cos(2*$kfr))/16/$kfr/$kfr/$kfr/$kfr;
  $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
  }

  elsif($rkky==2)

  {my ($xx)=$r*$r/$D/$D;
  $jaa= $scale*(-$xx+$xx*$xx)*exp(-$aa*$xx);
  $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;  $jcb=0;$jca=0;
  }

  elsif($rkky==5)

  {my ($xx)=$rx*$rx/$Da/$Da+$ry*$ry/$Db/$Db+$rz*$rz/$Dc/$Dc;
  $jaa= $scale*(-$xx+$xx*$xx)*exp(-$aa*$xx);
  $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;  $jcb=0;$jca=0;
  }

  else

  {
# muB=0.927405e-23 Ampere m^2
# mu0/4 pi=1e-7 kgm s^-2 Amp^-2
# m^3=10^30Angstroem^3
# 1meV= 16.0218e-23 J
#
#c=(mu0/4pi)(gJ muB)^2=0.92740^2  Angstroem^3 meV/16.0218

  my $c = $gJthis * $gJ * .927405 * .927405 / 16.02183;  #[meV A^3]
  $jaa = $c * (3 * $rx * $rx -$r *$r) /$r /$r /$r /$r /$r;
  $jbb = $c * (3 * $ry * $ry -$r *$r) /$r /$r /$r /$r /$r;
  $jcc = $c * (3 * $rz * $rz -$r *$r) /$r /$r /$r /$r /$r; 

  $jab = $c * (3 * $rx * $ry) /$r /$r /$r /$r /$r;
  $jbc = $c * (3 * $ry * $rz) /$r /$r /$r /$r /$r;
  $jac = $c * (3 * $rx * $rz) /$r /$r /$r /$r /$r;
   $jba=$jab;   $jcb=$jbc;   $jca=$jac;
  }     
 }
 return ([$jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc]);  
}

# Get lattic data, reading it from file 

sub getlattice {
    my ($file) = @_;
    my $h = new FileHandle;
    my $n = 0;
    $nofcomponents=0;$nofmagneticatoms=0;
  # input data int piddle
  if(open($h,$file))
  {      while(<$h>)
     {#next if /^\s*#/;
      # detect a= b= c= ...
      if ($a==0){($a)=extract("a",$_);}
      if ($b==0){($b)=extract("b",$_);}
      if ($c==0){($c)=extract("c",$_);}
      if ($alpha==0){($alpha)=extract("alpha",$_);}
      if ($beta==0){($beta)=extract("beta",$_);}
      if ($gamma==0){($gamma)=extract("gamma",$_);}

      if ($r1x==0){($r1x)=extract("r1a",$_);}
      if ($r1y==0){($r1y)=extract("r1b",$_);}
      if ($r1z==0){($r1z)=extract("r1c",$_);}
      if ($r2x==0){($r2x)=extract("r2a",$_);}
      if ($r2y==0){($r2y)=extract("r2b",$_);}
      if ($r2z==0){($r2z)=extract("r2c",$_);}
      if ($r3x==0){($r3x)=extract("r3a",$_);}
      if ($r3y==0){($r3y)=extract("r3b",$_);}
      if ($r3z==0){($r3z)=extract("r3c",$_);}

      if ($r1x==0){($r1x)=extract("r1x",$_);}
      if ($r1y==0){($r1y)=extract("r1y",$_);}
      if ($r1z==0){($r1z)=extract("r1z",$_);}
      if ($r2x==0){($r2x)=extract("r2x",$_);}
      if ($r2y==0){($r2y)=extract("r2y",$_);}
      if ($r2z==0){($r2z)=extract("r2z",$_);}
      if ($r3x==0){($r3x)=extract("r3x",$_);}
      if ($r3y==0){($r3y)=extract("r3y",$_);}
      if ($r3z==0){($r3z)=extract("r3z",$_);}

      if ($nofatoms==0){($nofatoms)=extract("nofatoms",$_);}
      if ($nofcomponents==0){($nofcomponents)=extract("nofcomponents",$_);}


      if (/^(#!|[^#])*nofneighbours\s*=\s*/){++$n;

                                   ($nofneighbours[$n])=extract("nofneighbours",$_);
                                   ($diagonalexchange)=extract("diagonalexchange",$_);
                                   if($diagonalexchange>1){print "Warning program makenn: diagonalexchange=$diagonalexchange not implemented setting diagonalexchange=0\n";$diagonalexchange=0;}
                                   ($x[$n])=extract("da",$_);
                                   ($y[$n])=extract("db",$_);
                                   ($z[$n])=extract("dc",$_);
                                   if (/^.*\Qx=\E/){($x[$n])=extract("x",$_);}
				   if (/^.*\Qy=\E/){($y[$n])=extract("y",$_);}
				   if (/^.*\Qz=\E/){($z[$n])=extract("z",$_);}
				     ($sipffilename)=extractstring("sipffilename",$_);
                                     ($charge[$n])=extractfromfile("CHARGE",$sipffilename);
                                     if($charge[$n]==""){$charge[$n]=$sipffilename;}                                 
                                             #               print "$sipffilename  charge=".$charge[$n]."\n";
                                    ($gJ[$n])=extractfromfile("GJ",$sipffilename); 
                                      unless(open(Fin,$sipffilename)) {die"Error opening $sipffilename\n";}
                                     $line=<Fin>;close Fin;($module[$n])=extractstring("MODULE",$line);
                                    if($cfph!=0){($magnetic)=extractfromfile("MAGNETIC",$sipffilename); 
                                     if($magnetic!=0){unless($module[$n]=~/so1ion/){die "Error makenn:  ion $n in mcphas.j sipf=$sipffilename for option -cfph MODULE must be so1ion\n"; }
                                                      ++$nofmagneticatoms;$nph[$nofmagneticatoms+$nofatoms]=$n;
                                                      $sipf_file[$nofmagneticatoms+$nofatoms]=$sipffilename;
                                                      $gJ[$nofmagneticatoms+$nofatoms]=$gJ[$n];
                                                      $x[$nofmagneticatoms+$nofatoms]=$x[$n];
                                                      $y[$nofmagneticatoms+$nofatoms]=$y[$n];
                                                      $z[$nofmagneticatoms+$nofatoms]=$z[$n];
                                                      $charge[$nofmagneticatoms+$nofatoms]=$charge[$n];
                                                     }  }   
                                   if ($bvk){clearMP($sipffilename);
                                             foreach(@sipf_file) {if(/$sipffilename/&&!$tabout){print "# $sipffilename found twice or more in mcphas.j - storing bvk parameters for ion $n in results/$sipffilename.$n \n";
                                                                                      copy($sipffilename,"results/$sipffilename.$n");
                                                                                      $sipffilename="results/".$sipffilename.".".$n;
                                                                                     }
                                                                 } 
                                            }                 
                                     $sipf_file[$n]=$sipffilename;
				  }

     }

     close $h; 

     if ($n!=$nofatoms) {print STDOUT "# Failed to read data file \"$file\": wrong number of atoms\n";

                         return undef;}
     if ($alpha<=0) { die "ERROR makenn: reading unit cell angle alpha=$alpha <=0\n";}
     if ($beta<=0) { die "ERROR makenn: reading unit cell angle beta=$beta <=0\n";}
     if ($gamma<=0) { die "ERROR makenn: reading unit cell angle gamma=$gamma <=0\n";}

     if($nofcomponents==0) {$nofcomponents=3;}

     return ([$a,$b,$c,$alpha,$beta,$gamma,$nofatoms,$nofcomponents],pdl [[$r1x,$r1y,$r1z],[$r2x,$r2y,$r2z],[$r3x,$r3y,$r3z]]);

    } else {
	print STDOUT "Warning: failed to read lattice from file \"$file\"\n";
	return undef;
    }
}

sub printlattice {
  my ($filein,$fileout)=@_;   
    my $h = new FileHandle;
    my $l = new FileHandle;
    unless (open($l,$fileout)){die "cannot open file $fileout\n";}
    unless (open($h,$filein)) {die "cannot open file $fileout\n";}

     $nofatoms=0;
     while(<$h>)
     {#next if /^\s*#/;
      $text=$_;
      if ($nofatoms==0){($nofatoms)=extract("nofatoms",$_);
          
            if($cfph==1&&$nofatoms!=0){$nofatomsnew=$nofatoms+$nofmagneticatoms;
            $text="#! nofatoms= $nofatomsnew  nofcomponents=48  number of atoms in primitive unit cell/number of components of each spin\n";
                        }
            }
# removed because iterative call of makenn will make fileheader longer
#if ($nofatoms!=0){
#        print $l "#-------------------------------------------------------------------------------------\n";
#        print $l "# output of program makenn $rmax - table with neighbors and interactions\n";
#        print $l "# Reference: M. Rotter et al. PRB 68 (2003) 144418\n";
#                 }
      print $l ($text);
      last if ($nofatoms!=0); # the line nofatoms= must be the last line of the file header !!!!
     }
     if ($nofatoms==0){die "ERROR makenn: unable to find 'nofatoms=' in file $filein\n";}
     
 return ($h,$l);
}



sub printneighbourlist {
  my ($h,$l,$nofn,$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb)=@_;
     
     print $l ("#*************************************************************************\n");
     while(<$h>)
     {$text=$_;
     if (/^(#!|[^#])*nofneighbours\s*=\s*/){($nn0)=extract("nofneighbours",$text);
                                            $text=~s!nofneighbours\s*=\s*\d+!nofneighbours=$nofn!;}
     if (/^(#!|[^#])*diagonalexchange\s*=\s*/){

      if ($rkky>=1||($readtable>0&&$DM<1&&$Jp<1))
       {$text=~s!diagonalexchange\s*=\s*\d+!diagonalexchange=1!;}
      else
       {$text=~s!diagonalexchange\s*=\s*\d+!diagonalexchange=0!;}
      }
      last if (/^(#!|[^#])*diagonalexchange\s*=\s*/);
     }
      print $l ($text);
# the next lines are to advance $h to the end of the numeric table

      while(<$h>)
      {last if ($nn0==0);
       unless (/^\s*#/){--$nn0;}
      }



if ($bvk==1)
{print $l ("# it follows the Born von Karman model according to  springs read from file $bvk_file\n");}
elsif ($rkky==1)
{print $l ("# it follows output of RKKY interaction according to J(R)=A.cos(2.kf.R)/(2.kf.R)^3 with A=$scale meV and kf=$kf A^-1 generated by makenn\n");}
elsif ($rkky==3)
{print $l ("# it follows output of RKKY interaction according to J(R)=A.[sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 with A=$scale meV and kf=$kf A^-1 generated by makenn\n");}
elsif ($rkky==2)
{print $l ("# kaneyoshi parametrization for the Bethe-Slater curve J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2] for scale A=$scale meV D=$D A alpha=$aa\n");}
elsif ($rkky==4)
{print $l ("# it follows output of RKKY interaction J(R)=A.cos(2.kfR)/(2.kfR)^3 for scale A=$scale meV and\n");
 print $l ("# kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n");}
elsif($rkky==5)
{print $l ("# kaneyoshi parametrization for the Bethe-Slater curve J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2] for scale A=$scale meV\n");
 print $l ("# with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2) with Da=$Da A Db=$Db A Dc=$Dc A and alpha=$aa\n");}
elsif($rkky==6)
 {print $l ("# it follows output of RKKY interaction  J(R)=A [sin(2.kfR)-2.kfR.cos(2.kfR)]/(2.kfR)^4 for scale A=$scale meV\n");
  print $l ("# kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n");}
elsif($readtable>0)
  {print $l ("# it follows output generated from interactions read from table $table_file\n");
  }
else
{print $l ("# it follows output of classical DD interaction generated by makenn\n");}

    if($alpha!=90||$beta!=90||$gamma!=90)
     {print $l ("#da[a]    db[b]     dc[c]       Jii[meV]  Jjj[meV]  Jkk[meV]  Jij[meV]  Jji[meV]  Jik[meV]  Jki[meV]  Jjk[meV]  Jkj[meV] with j||b, k||(a x b) and i normal to k and j\n");}
    else
     {print $l ("#da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");}

my $pcout=">./results/makenn.a".$nnn.".pc";
unless (open($l1,$pcout)){die "cannot open file $pcout\n";}
print $l1 "#-------------------------------------------------------------------------------------\n";
print $l1 "#!  table with neighbors and charges for atom n=$nnn at da=".$x[$nnn]." db=".$y[$nnn]." dc=".$z[$nnn]." sipffilename=".$sipf_file[$nnn]."\n";
print $l1 "# output of program makenn:, Reference: M. Rotter et al. PRB 68 (2003) 144418\n";
print $l1 "#-------------------------------------------------------------------------------------\n";
    if($alpha!=90||$beta!=90||$gamma!=90)
     {print $l1 "#orthonormal coordinate system ijk is defined with respect to abc as j||b, k||(a x b) and i normal to k and j\n#charge[|e|]  di[A]   dj[A]   dk[A]        da[a]    db[b]    dc[c]   distance[A] atomnr\n";}
     else
     {print $l1 "#charge[|e|]  da[A]     db[A]     dc[A]          da[a]      db[b]      dc[c]     distance[A]   atomnr\n";}

 for ($n1=1;$n1<(($rn->dims)[0]);++$n1)
 {next if($rn->index($n)->at($n1)==0);
  # the position xyz is relative position (not absolute coordinate of neighbour)
 print $l sprintf("%+10.6f %+10.6f %+10.6f ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
 $ddd=$an->index($n)->at($n1);
  print $l1 sprintf("%8s   %+10.6f %+10.6f %+10.6f     ",$charge[$ddd],$in->index($n)->at($n1),$jn->index($n)->at($n1),$kn->index($n)->at($n1));
  print $l1 sprintf("%+10.6f %+10.6f %+10.6f ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
  print $l1 sprintf("%+10.6f     %s\n",$rn->index($n)->at($n1),$ddd);

  if (($gJ!=0&&$gJ[$ddd]==0)||($gJ==0&&$gJ[$ddd]!=0)){unless($cfph==1){ die "error makenn. mixing of atoms with gJ=0 (intermediate coupling) and gJ>0 not implemented\n";}}
  if (($rkky==0&&$readtable==0)||$bvk==1||($readtable>0&&$DM>0)||($readtable>0&&$Jp>0)) # here anisotropic interaction comes in
   { 
       if($bvk==1||($gJ!=0&&$gJ[$ddd]!=0))
         {print $l sprintf("%+10.9e %+10.9e %+10.9e ",$Jaa->index($n)->at($n1),$Jbb->index($n)->at($n1),$Jcc->index($n)->at($n1));
          for($i=4;$i<=$nofcomponents;++$i){print $l "0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",$Jab->index($n)->at($n1),$Jba->index($n)->at($n1));
          print $l sprintf("%+10.9e %+10.9e ",$Jac->index($n)->at($n1),$Jca->index($n)->at($n1));
          for($i=4;$i<=$nofcomponents;++$i){print $l "0 0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",$Jbc->index($n)->at($n1),$Jcb->index($n)->at($n1));
          for($i=4;$i<=$nofcomponents;++$i){for($ii=$i+1;$ii<=$nofcomponents;++$ii){print $l "0 0 ";}} # add other components
         }
       elsif ($gJ==0&&$gJ[$ddd]==0)
         {print $l sprintf("%+10.9e %+10.9e %+10.9e %+10.9e %+10.9e %+10.9e ",4*$Jaa->index($n)->at($n1),$Jaa->index($n)->at($n1),4*$Jbb->index($n)->at($n1),$Jbb->index($n)->at($n1),4*$Jcc->index($n)->at($n1),$Jcc->index($n)->at($n1));
          for($i=7;$i<=$nofcomponents;++$i){print $l "0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jaa->index($n)->at($n1)),(2*$Jaa->index($n)->at($n1))); #SaLa LaSa
          print $l sprintf("%+10.9e %+10.9e ",(4*$Jab->index($n)->at($n1)),(4*$Jba->index($n)->at($n1))); #SaSb SbSa
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jab->index($n)->at($n1)),(2*$Jba->index($n)->at($n1))); #SaLb LbSa
          print $l sprintf("%+10.9e %+10.9e ",(4*$Jac->index($n)->at($n1)),(4*$Jca->index($n)->at($n1))); #SaSc ScSa
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jac->index($n)->at($n1)),(2*$Jca->index($n)->at($n1))); #SaLc LcSa
          for($i=7;$i<=$nofcomponents;++$i){print $l "0 0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jab->index($n)->at($n1)),(2*$Jba->index($n)->at($n1))); #LaSb SbLa
          print $l sprintf("%+10.9e %+10.9e ",$Jab->index($n)->at($n1),$Jba->index($n)->at($n1));         #LaLb LbLa
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jac->index($n)->at($n1)),(2*$Jca->index($n)->at($n1))); #LaSc ScLa
          print $l sprintf("%+10.9e %+10.9e ",$Jac->index($n)->at($n1),$Jca->index($n)->at($n1));         #LaLc LcLa
          for($i=7;$i<=$nofcomponents;++$i){print $l "0 0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jbb->index($n)->at($n1)),(2*$Jbb->index($n)->at($n1))); #SbLb LbSb
          print $l sprintf("%+10.9e %+10.9e ",(4*$Jbc->index($n)->at($n1)),(4*$Jbc->index($n)->at($n1))); #SbSc ScSb
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jbc->index($n)->at($n1)),(2*$Jcb->index($n)->at($n1))); #SbLc LcSb
          for($i=7;$i<=$nofcomponents;++$i){print $l "0 0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jbc->index($n)->at($n1)),(2*$Jcb->index($n)->at($n1))); #LbSc ScLb
          print $l sprintf("%+10.9e %+10.9e ",$Jbc->index($n)->at($n1),$Jcb->index($n)->at($n1));         #LbLc LcLb
          for($i=7;$i<=$nofcomponents;++$i){print $l "0 0 ";} # add other components
          print $l sprintf("%+10.9e %+10.9e ",(2*$Jcc->index($n)->at($n1)),(2*$Jcc->index($n)->at($n1))); #ScLc LcSc
          for($i=7;$i<=$nofcomponents;++$i){for($ii=$i+1;$ii<=$nofcomponents;++$ii){print $l "0 0 ";}} # add other components
         }
   }
   else  #here the isotropic interaction is written
   {
    if($gJ!=0&&$gJ[$ddd]!=0) 
    {print $l sprintf("%+10.9e %+10.9e %+10.9e ",$Jaa->index($n)->at($n1),$Jbb->index($n)->at($n1),$Jcc->index($n)->at($n1));
     for($i=4;$i<=$nofcomponents;++$i){print $l "0 ";} # add other components
    }
  
    if ($gJ==0&&$gJ[$ddd]==0) # spin - spin interactions only
    {print $l sprintf("%+10.9e 0.0 %+10.9e 0.0 %+10.9e 0.0 ",$Jaa->index($n)->at($n1),$Jbb->index($n)->at($n1),$Jcc->index($n)->at($n1));
     for($i=7;$i<=$nofcomponents;++$i){print $l "0 ";} # add other components
    }
   }

 if ($calcdist==1) {print $l sprintf("%+10.9e a%s",$rn->index($n)->at($n1),$ddd);}

  print $l "\n";

 }

 close $l1;
}



sub endprint {

  my ($h,$l)=@_;   

     close $h; 

     close $l;

}

# **********************************************************************************************
# extracts number from string
# 
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
# 
# ... it stores 0.3 in the variable $standarddeviation
#
sub extract { 
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************

# **********************************************************************************************
# extracts string from string
# 
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
# 
# ... it stores 0.3 in the variable $standarddeviation
#
sub extractstring { 
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*\b([^\n\s]+)[\s\n]/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************




# **********************************************************************************************
# extracts number from file
# 
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:         
#
# ($standarddeviation)=extractfromfile("sta","data.dat");
# 
# ... it stores 0.24 in the variable $standarddeviation
#
sub extractfromfile {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";$value="";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);}}
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\" to extract variable $variable\n";
             }
             return $value;
            }
# **********************************************************************************************
