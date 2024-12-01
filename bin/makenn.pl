#!/usr/bin/perl

use FileHandle;
use PDL;
use File::Copy;
use Getopt::Long;

# Loads the tables of symmetry equivalent positions, and atomic/ionic information
push @INC, $ENV{'MCPHASE_DIR'}.'/bin/';
require 'elements.pl';

print "#********************************************************\n";
print "# makenn 230406 - create table with neighbors and interactions\n";
print "# References: M. Rotter et al. PRB 68 (2003) 144418\n";
print "#********************************************************\n";
$PI=3.14159265358979323846;
 $bvkA=25;
 $bvkalpha=0.1;
 $SMALLdabc=0.00001; # small difference in da , db or dc ...
 $SMALLdr=0.00001; # small difference in distance so that program believes when
                #  comparing table values that this is the same bond
 $delta=0.0001; # for numerical derivative djdx djdy djdz displacement in Angstroem
 $Cel= zeroes (7,7); # for storage of elastic constants if needed
  # born van karman longitudinal springs:$bvkA*exp(-$bvkalpha*$r*$r);*$r*$r);
exit usage() if ($#ARGV<0);
 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($rmax) = eval $ARGV[0];
$rkky=0;$calcdist=0;$readtable=0;$classdip=0;
shift @ARGV; 
@storeargv=@ARGV;
# Parses command line options
GetOptions("rkky3d=s{4}"=>\@rkky3d,
           "kaneyoshi3d=s{5}"=>\@kaneyoshi3d,
           "rkkz3d=s{4}"=>\@rkkz3d,
           "rkky=s{2}"=>\@rkkydummy,
           "kaneyoshi=s{3}"=>\@kaneyoshi,
           "rkkz=s{2}"=>\@rkkz,
           "bvk:s"=>\$bvk,  # : is for optional arguments to follow
           "cfph:s"=>\$cfph,
           "e:s"=>\$e,
           "f:s"=>\$f,
           "jp:s"=>\$jp,
           "dm:s"=>\$dm,
           "d"=>\$d,
	   "npc"=>\$npc,
           "nm"=>\$nm,
           "r"=>\$cfphr,
           "rfunc"=>\$rfunc,
           "djdx"=>\$djdx,
           "djdy"=>\$djdy,
           "djdz"=>\$djdz);

# @ARGV=@storeargv;

die "djdx djdy djdz are exclusive options and cannot be used together\n" if defined ($djdx and $djdy) or  ($djdz and $djdy) or ($djdx and $djdz);
if ($bvk||$cfph){die "djdx djdy djdz cannot be used with bvk and cfph\n" if ($djdx||$djdy||$djdz);}

$ext=".j"; 
if ($djdx) {$ext=".djdx"; }
if ($djdy) {$ext=".djdy"; }
if ($djdz) {$ext=".djdz"; }

#******************************************************** treat options 
#$_=$ARGV[0];
#if(/-nm/){shift @ARGV;$_=$ARGV[0];}   # no new method of finding neighbours
#if(/-rfunc/){shift @ARGV;$_=$ARGV[0];}   # no new method of finding neighbours
#if(/-npc/){shift @ARGV;$_=$ARGV[0];}  # do not keep pointcharge files makenn.a*.pc
if(@rkky3d)
  {$rkky=4;
   $s=$rkky3d[0];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $scale=eval $s;
   $s=$rkky3d[1];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $ka=eval $s;
   $s=$rkky3d[2];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $kb=eval $s;
   $s=$rkky3d[3];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $kc=eval $s;
 #  shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=$ARGV[0];shift @ARGV;
 #  $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$ka=eval $ARGV[0];shift @ARGV;  
 #  $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kb=eval $ARGV[0];shift @ARGV;  
 #  $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kc=eval $ARGV[0];shift @ARGV;
   print "#calculating RKKY interaction J(R)=A.cos(2.kfR)/(2.kfR)^3 for scale A=$scale meV and\n";
   print "#kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n";}
elsif(@kaneyoshi3d)
  {$rkky=5;
   $s=$kaneyoshi3d[0];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $scale=eval $s;
   $s=$kaneyoshi3d[1];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $Da=eval $s;
   $s=$kaneyoshi3d[2];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $Db=eval $s;
   $s=$kaneyoshi3d[3];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $Dc=eval $s;
   $s=$kaneyoshi3d[4];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $aa=eval $s;
#shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=$ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Da=eval $ARGV[0];shift @ARGV;   
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Db=eval $ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Dc=eval $ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$aa=eval $ARGV[0];shift @ARGV;
   print "#calculating kaneyoshi parametrization for the Bethe-Slater curve\n";
   print "# J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2] for scale A=$scale meV \n";
   print "#with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2) with Da=$Da A Db=$Db A Dc=$Dc A  and alpha=$aa\n";}
elsif(@rkkz3d)
  {$rkky=6;
   $s=$rkkz3d[0];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $scale=eval $s;
   $s=$rkkz3d[1];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $ka=eval $s;
   $s=$rkkz3d[2];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $kb=eval $s;
   $s=$rkkz3d[3];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $kc=eval $s;

#   shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=$ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$ka=eval $ARGV[0];shift @ARGV;  
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kb=eval $ARGV[0];shift @ARGV;  
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kc=eval $ARGV[0];shift @ARGV;   
   print "#calculating RKKY interaction J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 for scale A=$scale meV\n";
   print "#kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n";}
elsif(@rkkydummy)
  {$rkky=1;
   $s=$rkkydummy[0];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $scale=eval $s;
   $s=$rkkydummy[1];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $kf=eval $s;
#shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g; $scale=eval $ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kf=eval $ARGV[0];shift @ARGV;
  $rfuncheader="RKKY interaction J(R)=A.cos(2.kf.R)/(2.kf.R)^3 for scale A=$scale meV and kf=$kf A^-1\n";
   print "#calculating".$rfuncheader;
  }
elsif(@kaneyoshi)
  {$rkky=2;
  $s=$kaneyoshi[0];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $scale=eval $s;
  $s=$kaneyoshi[1];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $D=eval $s;
  $s=$kaneyoshi[2];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $aa=eval $s;

#shift @ARGV; 
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$scale=eval $ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$D=eval $ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$aa=eval $ARGV[0];shift @ARGV;
   $rfuncheader="kaneyoshi parametrization for the Bethe-Slater curve\n#J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2] for scale A=$scale meV D=$D A alpha=$aa\n";
   print "#calculating ".$rfuncheader;}
elsif(@rkkz)
  {$rkky=3;
    $s=$rkkz[0];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $scale=eval $s;
  $s=$rkkz[1];  $s=~s/exp/essp/g;$s=~s/x/*/g;$s=~s/essp/exp/g;    $kf=eval $s;

#shift @ARGV; 
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$scale=eval $ARGV[0];shift @ARGV;
#   $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$kf=eval $ARGV[0];shift @ARGV;
$rfuncheader="RKKY interaction J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 for scale A=$scale meV and kf=$kf A^-1\n";
   print "#calculating ".$rfuncheader;
  }
elsif(defined $bvk)
  {#$bvk=1;
  #shift @ARGV;
 #  unless ($#ARGV>=0) # if no filename is given print help
   unless ($bvk) # if no filename is given print help
   { $tabout=1;$bvk=1;}
  else
   { # if filename is given read file ...
     # $bfk_file=$ARGV[0];shift @ARGV;
   $bfk_file=$bvk; 
   $bvk=1;
   print "creating phononic interactions from Born van Karman model in file $bfk_file\n";
   # read interaction constants from file
   unless(open(Fin,$bfk_file)){ die "could not open $bfk_file\n";}
             $nof_springs=0;
             {while($line=<Fin>){next if $line=~/^\s*#/;
                                 my @numbers=split(" ",$line);
                                 if($#numbers>=4)
                                 {++$nof_springs;#print $line;
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
elsif(defined $cfph)
{ 
 # shift @ARGV;
 #$_=$ARGV[0];if(/-r/){shift @ARGV;}
 if($cfph){$screeningfile=$cfph;}
 print "creating crystal field phonon interactions from pointcharge model using program pointc\n";
 #$screeningfile=$ARGV[0];
 $cfph=1;
}
elsif(defined $e||defined $f||defined $dm||defined $jp)
  {
   #shift @ARGV;
   unless ($e||$f||$dm||$jp) # if no filename is given print help
   {
 $tabout=1; }
  else
 { $table_file=$e.$f.$dm.$jp;
   #$table_file=$ARGV[0];shift @ARGV;
   $ignore_neihgbours_behind=0;
   print "reading "; if ($e) {print " isotropic ";$DM=0;$Jp=0;} 
                  elsif ($f) {print " isotropic ";$DM=0;$Jp=0;} 
                  elsif ($jp) {print " Jp ";$DM=0;$Jp=1;} 
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
                                 if($e||$jp)
                                 {++$n_table;
                                  $da[$n_table]=$numbers[0];
                                  $Jex[$n_table]=$numbers[1];
                                 }
                                 if($f||$dm)
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
 else {$classdip=1;} # without options calculate classical dipole interactions
 
if($d) # for option -d
  {$calcdist=1;print "putting distance of neighbors (A) to last column of makenn.j\n";}

#******************************************************** get crystal structure 
my ($latt,$p) = getlattice("./mcphas.j"); # gets lattice and atomic positions
my ($a,$b,$c,$alpha,$beta,$gamma,$nofatoms,$nofcomponents) = @{$latt};
unless($tabout) {
 print "rmax=".$rmax." A\n";
 print "a=".$a." b=".$b." c=".$c." alpha=".$alpha." beta=".$beta." gamma=".$gamma."\n";
 print "Primitive Lattice Vectors (rows) in units of lattice vectors [abc] :".$p."\n";
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


$p= $p x $rtoijk; # primitive lattice in Euclidean ijk coordinates
$pinv=inv($p); #invert this matrix for use later
$pVolume=inner($p->slice(":,(0)"),crossp($p->slice(":,(1)"),$p->slice(":,(2)")));
if($pVolume<0){die "Error makenn - primitive lattice not right handed \n";}
unless($tabout){ print "Primitive Lattice Vectors (rows) in Euclidean Coordinates [A]:".$p."\nPrimitive Unit Cell Volume [A^3]: ".$pVolume."\n";}
    $r=0;
# select a Vector from primitive lattice print $p->slice(":,(0)");
# calculate Volume of primitive unit cell in A^3

#********************************************************
# first determine maximum distance of a basis atom to origin of unit cell
# and check if atoms are in primitive unit cell - if no - move them there
    $distmax=0;
  for ($nz=1;$nz<=$nofatoms;++$nz){
   $dabc=pdl [($x[$nz]),($y[$nz]),($z[$nz])];	
   $rvec= $dabc x $rtoijk;$rvec=$rvec->slice(":,(0)");
   $rr=inner($rvec, $rvec);
   $r=sqrt($rr->at());
   if($r>$distmax){$distmax=$r;}
  # check if atom is in primitive unit cell
  # determine atomic coordinates
    $dr1r2r3=$rvec x $pinv; $dr1r2r3=$dr1r2r3->slice(":,(0)");
  # if not translate it into primitive unit cell
  $moved=0;
    while(($dr1r2r3->at(0))>=1){--$dr1r2r3->slice("(0)");$moved=1;}
    while(($dr1r2r3->at(1))>=1){--$dr1r2r3->slice("(1)");$moved=1;}
    while(($dr1r2r3->at(2))>=1){--$dr1r2r3->slice("(2)");$moved=1;}
    while(($dr1r2r3->at(0))<-$SMALLdabc){++$dr1r2r3->slice("(0)");$moved=1;}
    while(($dr1r2r3->at(1))<-$SMALLdabc){++$dr1r2r3->slice("(1)");$moved=1;}
    while(($dr1r2r3->at(2))<-$SMALLdabc){++$dr1r2r3->slice("(2)");$moved=1;}
   $rvec=$dr1r2r3 x $p; $rvec=$rvec->slice(":,(0)");
   $dabc=$rvec x $invrtoijk; $dabc=$dabc->slice(":,(0)");
  if($moved==1){print STDERR "Warning: atom $nz at da=$x[$nz] db=$y[$nz] dc=$z[$nz] not in primitive unit cell - translating"; 
  $x[$nz]=$dabc->at(0);
  $y[$nz]=$dabc->at(1);
  $z[$nz]=$dabc->at(2);
                print STDERR " it to da=$x[$nz] db=$y[$nz] dc=$z[$nz]\n";
               }
  }

#********************************************************
# Determine the nmin and nmax values in order to cover with the lattice vectors
# the desired sphere including all neighbours up to rmax Angstroem
  unless($readtable>1)
 {if(defined $nm){
# OLD : determine $nmin,$nmax by looking at a cube with side 3rmax
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
}else{
# NEW 8.4.2024 FASTER: determine $nmin,$nmax by looking at condition that
# plane of parallelepiped must be more distant than radius of sphere given in command line
# e.g.
$rvec=crossp($p->slice(":,(1)"),$p->slice(":,(2)")); # $rvec = $r2 x $r3
$rr=inner($rvec, $rvec);
$rvec/=sqrt($rr);

# (n1max-1) * $r1 . $rvec )> $rmax  ... i.e. n1max > 1.0 + $rmax / ($r1 . $rvec)
# (n1min+1) * $r1. $rvec) < -$rmax           n1min < -1.0 - $rmax / ($r1 . $rvec)
$n1max=my_ceil(1.0+$rmax / inner($p->slice(":,(0)"),$rvec));
$n1min=my_floor(-1.0-$rmax / inner($p->slice(":,(0)"),$rvec));
# similar ...
$rvec=crossp($p->slice(":,(2)"),$p->slice(":,(0)")); # $rvec = $r3 x $r1
$rr=inner($rvec, $rvec);
$rvec/=sqrt($rr);
$n2max=my_ceil(1.0+$rmax / inner($p->slice(":,(1)"),$rvec));
$n2min=my_floor(-1.0-$rmax / inner($p->slice(":,(1)"),$rvec));
$rvec=crossp($p->slice(":,(0)"),$p->slice(":,(1)")); # $rvec = $r1 x $r2
$rr=inner($rvec, $rvec);
$rvec/=sqrt($rr);
$n3max=my_ceil(1.0+$rmax / inner($p->slice(":,(2)"),$rvec));
$n3min=my_floor(-1.0-$rmax / inner($p->slice(":,(2)"),$rvec));
}
  }

#******************************************************** put out a sample table if desired
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

if(defined $e)
{print STDOUT << "EOF";
# Table of exchange interaction constants J - assumed to be isotropic according to in H= -1/2 sum_ij J Ji.Jj 
# |Rij| [A]  J [meV]  atom_i  atom_j 
EOF

}
elsif(defined $jp)
{print STDOUT << "EOF";
# Table of exchange interaction constants Jp  in H= -1/2 sum_ij Jp (Ji.R)(Jj.R)/R^2 
# |Rij| [A]  Jp [meV]  atom_i  atom_j 
EOF

}
elsif(defined $f)
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
}} #fi tabout 
else
{
print "# $n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";

     # initialize output file results/makenn$ext
  ($h,$l)=printlattice("./mcphas.j",">./results/makenn$ext");
print "# number of atoms = $nofatoms\n calculating ...\n";
}     
print_time_estimate_until_end(0.0);
#********************************************************        
@atoms=();   # initialize array to push output into ...
 for ($nnn=1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)    
 {unless($tabout){ if($nnn>$nofatoms){print "# splitting atom ".$nph[$nnn]." into nucleus and magnetic electron shell creating a new magnetic \n";}
   print "# atom $nnn ...";} 
  if($nnn>$nofatoms){print "with shifted dc=".($z[$nnn]+0.1/$c)." to distinguish it from the nucleus at dc=".$z[$nnn]."...";}
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
     my $Gmix= zeroes (7,4); # for storage of mixing term phonon strain interaction constants

    # go through neighbours and determine interaction
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
 for($ntbl=1;$ntbl<=$n_table;++$ntbl){
if((defined $e||defined $jp)&&abs($da[$ntbl]-$r)<$SMALLdr){$ff=1;}
elsif((defined $f||defined $dm)&&
      abs($da[$ntbl]-$xx)<$SMALLdabc&&
      abs($db[$ntbl]-$yy)<$SMALLdabc&&
      abs($dc[$ntbl]-$zz)<$SMALLdabc){$ff=1;}
elsif((defined $bvk)&&(abs($da[$ntbl]-$r)<$SMALLdr)){
       $a1=$sipf_file[$nnn];$a2=$sipf_file[$nz];
       if(($db[$ntbl]=~/$a1/&&$dc[$ntbl]=~/$a2/)
       ||($db[$ntbl]=~/$a2/&&$dc[$ntbl]=~/$a1/)){$ff=1;}}
                                  }
if($ff==0){++$n_table;$ntbl=$n_table;
        
if(defined $e||defined $jp)
{$da[$ntbl]=$r;
 print sprintf("%+10.6f     0       a%i a%i \n",$r,$nnn,$nz);
}
elsif(defined $f)
{$da[$ntbl]=$xx;$db[$ntbl]=$yy;$dc[$ntbl]=$zz;
 print sprintf("%+10.6f %+10.6f %+10.6f 0        a%i a%i %+10.6f\n",$xx, $yy ,$zz,$nnn,$nz,$r);
}
elsif(defined $dm)
{$da[$ntbl]=$xx;$db[$ntbl]=$yy;$dc[$ntbl]=$zz;
  print sprintf("%+10.6f %+10.6f %+10.6f    %+10.6f %+10.6f %+10.6f     a%i a%i %+10.6f\n",$xx, $yy ,$zz,$rvec->at(0),$rvec->at(1),$rvec->at(2),$nnn,$nz,$r);
}
elsif(defined $bvk)
{$da[$ntbl]=$r; $db[$ntbl]=$sipf_file[$nnn];$dc[$ntbl]=$sipf_file[$nz];
 $spring= $bvkA*exp(-$bvkalpha*$r*$r);
 print sprintf("%s   %s    %+10.6f   %+10.6f   0 \n",$sipf_file[$nnn],$sipf_file[$nz],$r,$spring);
}   
else{die "Error makenn - creating sample table for options -e -f -dm -jp \n";}
          }
              } #tabout

          unless($readtable>0){
    $an=$an->append( pdl ([$nz]));
    $rn=$rn->append( pdl ([$r]));
    $xn=$xn->append( pdl ([$xx]));
    $yn=$yn->append( pdl ([$yy]));
    $zn=$zn->append( pdl ([$zz]));
    $in=$in->append( pdl ([$rvec->at(0)]));
    $jn=$jn->append( pdl ([$rvec->at(1)]));
    $kn=$kn->append( pdl ([$rvec->at(2)]));

if ($rfunc&&$rfuncdone!=1&&$rfuncheader)  
{
 open(ff,">results/makenn.func.r");
 print ff "#".$rfuncheader;
 print ff "#r[A] J(meV) interaction-present\n";
 for($myr=0.1;$myr<$rmax;$myr+=0.1){
 my ($interaction) = getinteraction($Gmix,$gJ,$gJ[$nz],$sipffilename,$sipf_file[$nz],$myr,0,0,$myr);
 my ($Gmix,$jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc) = @{$interaction};
 print ff $myr." ".$jaa."  0.1\n";
    }
$rfuncdone=1;
}
# calculate the interaction
    my ($interaction) = getinteraction($Gmix,$gJ,$gJ[$nz],$sipffilename,$sipf_file[$nz],$r,$rvec->at(0),$rvec->at(1),$rvec->at(2));
    my ($Gmix,$jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc) = @{$interaction};
if ($rfunc&&$rfuncheader){print ff $r." ".$jaa." 1\n";}

    if($djdx||$djdy||$djdz){my $dx=0; my $dy=0; my $dz=0;
              if($djdx){$dx=$delta;}
              if($djdy){$dy=$delta;}
              if($djdz){$dz=$delta;}
              my $rr=sqrt(($rvec->at(0)+$dx)*($rvec->at(0)+$dx)+($rvec->at(1)+$dy)*($rvec->at(1)+$dy)+($rvec->at(2)+$dz)*($rvec->at(2)+$dz));
              my ($intd) = getinteraction($Gmix,$gJ,$gJ[$nz],$sipffilename,$sipf_file[$nz],$rr,$rvec->at(0)+$dx,$rvec->at(1)+$dy,$rvec->at(2)+$dz);
              my ($Gmix,$djaa,$djab,$djac,$djba,$djbb,$djbc,$djca,$djcb,$djcc) = @{$intd};
               $jaa=($djaa-$jaa)/$delta;
               $jab=($djab-$jab)/$delta;
               $jac=($djac-$jac)/$delta;

               $jba=($djba-$jba)/$delta;
               $jbb=($djbb-$jbb)/$delta;
               $jbc=($djbc-$jbc)/$delta;

               $jca=($djca-$jca)/$delta;
               $jcb=($djcb-$jcb)/$delta;
               $jcc=($djcc-$jcc)/$delta;
} # djdx djdy djdz

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
                                 if((($f)&&
                                    abs($da[$ntbl]-$xx)<$SMALLdabc&&
                                    abs($db[$ntbl]-$yy)<$SMALLdabc&&
                                    abs($dc[$ntbl]-$zz)<$SMALLdabc&&
                                    $Jex[$ntbl]!=0.0)||
                                    (($e||$jp)&&
                                     abs($da[$ntbl]-$r)<$SMALLdr&&
                                     $Jex[$ntbl]!=0.0)||
                                   (($dm)&&
                                    abs($da[$ntbl]-$xx)<$SMALLdabc&&
                                    abs($db[$ntbl]-$yy)<$SMALLdabc&&
                                    abs($dc[$ntbl]-$zz)<$SMALLdabc&&
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
              } #DM
                                  } # save neighbour   

                         } # neighbour in readtable  
                               } #readtable
                      } # 0<r<rmax
    }}}} # n1  n2 n3 nz

if($rfuncdone==1){close ff;}

   $n= qsorti($rn); # sort neighbors according to r
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

unless($tabout){printneighbourlist($Gmix,$h,$nofneighbours[$nnn],$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb);}
print_time_estimate_until_end(($nofatoms+$nofmagneticatoms-$nnn)/$nnn);
 } # next atom nnn

unless($tabout){ endprint($h,$l);  }
 
if($cfph!=0){ # for cf phonon interaction recreate makenn$ext
my ($h,$l)=printlattice("./mcphas.j",">./results/makenn$ext");
@atoms=();
for($nnn=$nofatoms+1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)
 { # if screeningfile is given - use it and screen charges 
  if($screeningfile){
  system("fillcol 5 c1 results/makenn.a$nnn.pc");
  system("fillcol 6 c1 results/makenn.a$nnn.pc");
  system("mult 8 1 results/makenn.a$nnn.pc 1 2 ".$screeningfile);
  system("mult 8 5 results/makenn.a$nnn.pc 1 3 ".$screeningfile);
  system("mult 8 6 results/makenn.a$nnn.pc 1 4 ".$screeningfile);  
   # run pointc to create derivatives of Blm from pointcharge model for magnetic atoms
  system("pointc -d -o ".$sipf_file[$nnn]." results/makenn.a$nnn.pc 5 6 > results/makenn.a$nnn.sipf");
   } else {
   # run pointc to create derivatives of Blm from pointcharge model for magnetic atoms
  system("pointc -d -o ".$sipf_file[$nnn]." results/makenn.a$nnn.pc > results/makenn.a$nnn.sipf");
          }
  $sipf_file[$nnn]="results/makenn.a$nnn.sipf";
  # set nuclear scattering lengths in this purely magnetic atom zero (for mcdiff)
  setvariable("SCATTERINGLENGTHREAL",0,$sipf_file[$nnn]);
  setvariable("SCATTERINGLENGTHIMAG",0,$sipf_file[$nnn]);

  if($cfph==2){my_rename("results/pointc.dLlm","results/makenn.a$nnn.dLlm");unlink("results/pointc.dBlm");}
  if($cfph==3){my_rename("results/pointc.dBlm","results/makenn.a$nnn.dBlm");unlink("results/pointc.dLlm");}
 }
for($nnn=1;$nnn<=$nofatoms;++$nnn){$nofneighbours[$nnn]=0;}
for($nnn=$nofatoms+1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)
 {my $GG= zeroes(7,28); # to be filled with Gcfph

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
  if($cfph==2){unless (open(Fin,"results/makenn.a$nnn.dLlm")){die "cannot open file results/makenn.a$nnn.dLlm\n";}}
  if($cfph==3){unless (open(Fin,"results/makenn.a$nnn.dBlm")){die "cannot open file results/makenn.a$nnn.dBlm\n";}}
  while(<Fin>){next if /^\s*#/;
                   my @numbers=split(" ",$_); 
                   if($#numbers>=4)
                    {++$nofn;

   $rvec=pdl [$numbers[0],$numbers[1],$numbers[2]]; # relative position
 #  of charged nucleus with respect to magnetic ion in eucliedean coordinates [Angstroem]
   $aabbcc=$rvec x $invrtoijk;$aabbcc=$aabbcc->slice(":,(0)");
   $da=$aabbcc->at(0);
   $db=$aabbcc->at(1);
   $dc=$aabbcc->at(2); # ... transformed to lattice
  
   $abc=pdl [$x[$nnn],$y[$nnn],$z[$nnn]]; # positio vector of magnetic ion in lattice coordinates
   $Rnnn= $dabc x $rtoijk;$Rnnn=$Rnnn->slice(":,(0)"); # ... transformed to euclidean system [Angstroem]
  $Rx=$numbers[0]+$Rnnn->at(0);  # position vector of charged nucleus in euclidean syste [Angstroem]
  $Ry=$numbers[1]+$Rnnn->at(1);
  $Rz=$numbers[2]+$Rnnn->at(2);  # ... needed to calculate Gcfph
# $R=$Rx*$Rx+$Ry*$Ry+$Rz*$Rz; $R=sqrt($R);print $R." ";
  $outstring=sprintf("%+10.6f %+10.6f %+10.6f ",$da, $db ,($dc-0.1/$c));
  $outstringp=sprintf("%+10.6f %+10.6f %+10.6f ",-$da, -$db ,-$dc+0.1/$c);
  # now the off diagonal Elements: Jab Jba Jac Jca Jad Jda .. Jbc Jcb Jbd Jdb ... Jcd Jdc ... etc.
  # a b c d e ... = 11/x 11s/y 10/z  22s 21s 20 21 22  33s 32s ...
  #                 I1    I2    I3    I4  I5
  # for PHONON      ux    uy    uz    dum dum dum ...
  # for so1ion      O11  O11s   O10   O22s O21s ...
  # Hamiltonian Hcfph= - sum_i<j,lmgamma u_gamma(i) (-)dBlm(j)/du_gamma(i)  Olm(j) = 
  #                  =-sumi<j_alphabeta  Ialpha Jalphabeta Ibeta
  #                  =-1/2sumij_alphabeta  Ialpha Jalphabeta Ibeta
  # note: the factor 1/2 comes from the fact, that the pairs are counted twice
  # in the last line. The interaction Jalphabeta is the same as in the lines above.
  # i.e. J12=0 J13=0 J41=- dB22s/dux J14=0 J51=- dB21s/dux J15=0 ...
  # so1ion 5+9+13=27 ... 27x3=81  5+81=86
  # ic1ion 8+9+13=30 ... 27x3=81  8+81=89
  $soff=0;if($cfph==2){$off=3;}
  for($i=6+$soff;$i<=86+$soff;++$i){$outstring.= " ".(-$numbers[$i]);
                                    $outstringp.= " ".(-$numbers[$i]);
                                    }
  for($lm=1;$lm<=27;++$lm){
   $ix=($lm-1)*3+6; # $numbers[$ix]=dBlm/dux
   $iy=($lm-1)*3+6+1; # $numbers[$iy]=dBlm/duy
   $iz=($lm-1)*3+6+2; # $numbers[$iz]=dBlm/duz
   # $i ... dBlm/dux  $i+1  ... dBlm/duy  $i+2 ... dBlm/duz
   
   # voigt component 11 $lm
   $GG->slice("1,$lm")+=-$Rx*$numbers[$ix]/0.5292;  
 # a0=0.5292 Angstroem (needed because $Rx is in Angstroem and 
 # cf parameter derivative is in units of a0 as delivered by pointc program
   # voigt component 22 $lm
   $GG->slice("2,$lm")+=-$Ry*$numbers[$iy]/0.5292;  
   # voigt component 33 $lm
   $GG->slice("3,$lm")+=-$Rz*$numbers[$iz]/0.5292;  
   # voigt component 4=(23) $lm
   $GG->slice("4,$lm")+=-0.5*($Ry*$numbers[$iz]+$Rz*$numbers[$iy])/0.5292;  
   # voigt component 5=(13) $lm
   $GG->slice("5,$lm")+=-0.5*($Rx*$numbers[$iz]+$Rz*$numbers[$ix])/0.5292;  
   # voigt component 6=(12) $lm
   $GG->slice("6,$lm")+=-0.5*($Rx*$numbers[$iy]+$Ry*$numbers[$ix])/0.5292;  

   }

  $outstring.= "\n";
  $outstringp.= "\n";
   $ph[$nnn].=$outstring;
  # push the outstring to the appropriate PHONON ions neighbour list
  $ph[$ni[$nofn]].=$outstringp;
  ++$nofneighbours[$ni[$nofn]];
                    }
                   }
  close Fin;
  if($cfph==2){unlink("results/makenn.a$nnn.dLlm");}
  if($cfph==3){unlink("results/makenn.a$nnn.dBlm");}

  for($j=1;$j<=27;++$j){for($i=1;$i<=6;++$i){
  $Gcfph[$nnn].=" ".$GG->at($i,$j);}}
 }
# here output the phononic atoms with module phonon first
for($nnn=1;$nnn<=$nofatoms;++$nnn)
 {push @atoms, ("#*************************************************************************\n");
  push @atoms, ("#! da=".$x[$nnn]." [a] db=".$y[$nnn]." [b] dc=".$z[$nnn]." nofneighbours=".$nofneighbours[$nnn]." diagonalexchange=2  sipffilename=".$sipf_file[$nnn]."\n");
  push @atoms, ("# crystal field phonon interaction parameters from pointcharge calculation\n");
  push @atoms, ("# da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
  if($cfph==2){
  push @atoms, ("#! symmetricexchange=0 indexexchange= 1,7 2,7 3,7 1,8 2,8 3,8 1,9 2,9 3,9 1,10 2,10 3,10 1,11 2,11 3,11 ");
  #O4m
  push @atoms, (" 1,19 2,19 3,19 1,20 2,20 3,20 1,21 2,21 3,21 1,22 2,22 3,22 1,23 2,23 3,23 1,24 2,24 3,24 1,25 2,25 3,25 1,26 2,26 3,26 1,27 2,27 3,27 ");
  #O6m
  push @atoms, (" 1,39 2,39 3,39 1,40 2,40 3,40 1,41 2,41 3,41 1,42 2,42 3,42 1,43 2,43 3,43 1,44 2,44 3,44 1,45 2,45 3,45 1,46 2,46 3,46 1,47 2,47 3,47 1,48 2,48 3,48 1,49 2,49 3,49 1,50 2,50 3,50 1,51 2,51 3,51 \n");
              }
  if($cfph==3){
  push @atoms, ("#! symmetricexchange=0 indexexchange= 1,4 2,4 3,4 1,5 2,5 3,5 1,6 2,6 3,6 1,7 2,7 3,7 1,8 2,8 3,8 ");
  #O4m
  push @atoms, (" 1,16 2,16 3,16 1,17 2,17 3,17 1,18 2,18 3,18 1,19 2,19 3,19 1,20 2,20 3,20 1,21 2,21 3,21 1,22 2,22 3,22 1,23 2,23 3,23 1,24 2,24 3,24 ");
  #O6m
  push @atoms, (" 1,36 2,36 3,36 1,37 2,37 3,37 1,38 2,38 3,38 1,39 2,39 3,39 1,40 2,40 3,40 1,41 2,41 3,41 1,42 2,42 3,42 1,43 2,43 3,43 1,44 2,44 3,44 1,45 2,45 3,45 1,46 2,46 3,46 1,47 2,47 3,47 1,48 2,48 3,48 \n");
              }
# here we should enter the cf- phonon interactions for the MODULE=phonon oscillators   
#  1) calculate  number of neighbours (only the magnetic atoms) 
#  2) fill values from results/makenn.a$nnn.dBlm  [all done above for magnetic atoms and filled in here]
  push @atoms, $ph[$nnn];
 }

# now output the magnetic ions ...
for($nnn=$nofatoms+1;$nnn<=$nofatoms+$nofmagneticatoms;++$nnn)
 {push @atoms, ("#*************************************************************************\n");
  push @atoms, ("#! da=".$x[$nnn]." [a] db=".$y[$nnn]." [b] dc=".($z[$nnn]+0.1/$c)." nofneighbours=".($nofneighbours[$nnn]+1)." diagonalexchange=2  sipffilename=".$sipf_file[$nnn]."\n");
  push @atoms, ("# crystal field phonon interaction parameters from pointcharge calculation\n");
 
 push @atoms, ("# the mixing terms Gcfph in meV\n");
 push @atoms, ("#! Gindices=");
 my $Gind=pdl[7,8,9,10,11,19,20,21,22,23,24,25,26,27,39,40,41,42,43,44,45,46,47,48,49,50,51];

if($cfph==2){ for($i=1;$i<=27;++$i){
for($j=1;$j<=6;++$j){
push @atoms, sprintf(" %i,%i",$j,$Gind->at($i-1));
}}
}

if($cfph==3){ for($i=1;$i<=27;++$i){
for($j=1;$j<=6;++$j){
push @atoms, sprintf(" %i,%i",$j,$Gind->at($i-1)-3);
}}
}
push @atoms, "\n";
 push @atoms, ("#! G=");
 push @atoms, $Gcfph[$nnn];
push @atoms, "\n";

  push @atoms, ("# da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
    if($cfph==2){
  push @atoms, ("#! symmetricexchange=0 indexexchange= 7,1 7,2 7,3 8,1 8,2 8,3 9,1 9,2 9,3 10,1 10,2 10,3 11,1 11,2 11,3 ");
  # O4m
  push @atoms, (" 19,1 19,2 19,3 20,1 20,2 20,3 21,1 21,2 21,3 22,1 22,2 22,3 23,1 23,2 23,3 24,1 24,2 24,3 25,1 25,2 25,3 26,1 26,2 26,3 27,1 27,2 27,3 ");
  # O6m
  push @atoms, (" 39,1 39,2 39,3 40,1 40,2 40,3 41,1 41,2 41,3 42,1 42,2 42,3 43,1 43,2 43,3 44,1 44,2 44,3 45,1 45,2 45,3 46,1 46,2 46,3 47,1 47,2 47,3 48,1 48,2 48,3 49,1 49,2 49,3 50,1 50,2 50,3 51,1 51,2 51,3 \n"); 
               }
  if($cfph==3){
  push @atoms, ("#! symmetricexchange=0 indexexchange= 4,1 4,2 4,3 5,1 5,2 5,3 6,1 6,2 6,3 7,1 7,2 7,3 8,1 8,2 8,3 ");
  # O4m
  push @atoms, (" 16,1 16,2 16,3 17,1 17,2 17,3 18,1 18,2 18,3 19,1 19,2 19,3 20,1 20,2 20,3 21,1 21,2 21,3 22,1 22,2 22,3 23,1 23,2 23,3 24,1 24,2 24,3 ");
  # O6m
  push @atoms, (" 36,1 36,2 36,3 37,1 37,2 37,3 38,1 38,2 38,3 39,1 39,2 39,3 40,1 40,2 40,3 41,1 41,2 41,3 42,1 42,2 42,3 43,1 43,2 43,3 44,1 44,2 44,3 45,1 45,2 45,3 46,1 46,2 46,3 47,1 47,2 47,3 48,1 48,2 48,3 \n"); 
               }
  # here come the cf-phonon interactions for the MODULE=so1ion magnetic ions 
  # open pointcharge file 
  push @atoms, $ph[$nnn];

# and finally recreate corresponding phononic atoms sipf file to be a MODULE=PHONON file
if($cfphr){store_phonon_file($sipf_file[$nph[$nnn]]);}

}
endprint($h,$l);
} # fi cfph recreate makenn.j

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
} # fi bvk consistency check

print "created files: results/makenn$ext     (interaction parameters)\n";
print "               results/makenn.a*.pc   (pointcharge environment files)\n";
if($rfunc){
print "               results/makenn.func.r  (interaction as function of r)\n";
}
if($bvk){
print "               results/makenn.Cel     (elastic constants)\n";
}
print "********************************************************\n";
print "               end of program makenn\n";
print " Reference: M. Rotter et al. PRB 68 (2003) 144418\n";
print "********************************************************\n";
if($npc){unlink("results/makenn.a*.pc");}
   exit;

#----------------clear MODPAR in sipf files -------------------------------
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
                                                           $Kxyread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR5\s*=\s*[^\s\;\n\t\*]+/$1MODPAR5=$Kxyread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR6\s*=\s*/) {($Kxzread)=($line=~m/^(?:#!|[^#])*\bMODPAR6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxzread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR6\s*=\s*[^\s\;\n\t\*]+/$1MODPAR6=$Kxzread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR7\s*=\s*/) {($Kyzread)=($line=~m/^(?:#!|[^#])*\bMODPAR7\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kyzread=0;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR7\s*=\s*[^\s\;\n\t\*]+/$1MODPAR7=$Kyzread/g;}
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

# ---------------------- add K to to MODPARS in sipf files ----------------------------
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
                if($line=~/^(#!|[^#])*\bMODPAR6\s*=\s*/) {($Kxzread)=($line=~m/^(?:#!|[^#])*\bMODPAR6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kxzread+=$Kxz/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR6\s*=\s*[^\s\;\n\t\*]+/$1MODPAR6=$Kxzread/g;}
                if($line=~/^(#!|[^#])*\bMODPAR7\s*=\s*/) {($Kyzread)=($line=~m/^(?:#!|[^#])*\bMODPAR7\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
                                                           $Kyzread+=$Kyz/2;
                                                          $line=~s/(^(#!|[^#])*?\b)MODPAR7\s*=\s*[^\s\;\n\t\*]+/$1MODPAR7=$Kyzread/g;}
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
   my ($Gmix,$gJthis,$gJ,$sipffilethis,$sipffile,$r,$rx,$ry,$rz)=@_;
   my $n;

 if($bvk==1)
 { # here do the Born von Karman calculation using spring constants
   # $long_spring and $trans_spring are the spring constants in N/m
  $a0 = .5292e-10;#(m)
  $J2meV=1/1.60217646e-22; #1 millielectron volt = 1.60217646 . 10-22 joules
  $jaa=0;$jbb=0;$jcc=0;$jab=0;$jbc=0;$jac=0;
  for($n=1;$n<=$nof_springs;++$n){$a1=$atom_m[$n];$a2=$atom_n[$n];
         if(abs($r-$bondlength[$n])<$SMALLdr
            &&(($sipffilethis=~/(results\/)?$a2/&&$sipffile=~/(results\/)?$a1/)
             ||($sipffilethis=~/(results\/)?$a1/&&$sipffile=~/(results\/)?$a2/)
              )
            )
             {# bond found - do something
  $kL=$long_spring[$n]*$a0*$a0*$J2meV; 
  $kT=$trans_spring[$n]*$a0*$a0*$J2meV; 
  $Kxx= -(( $kL-$kT)* $rx * $rx + $kT* $r *$r)  /$r /$r;$jaa-=$Kxx;
  $Kyy= -(( $kL-$kT)* $ry * $ry + $kT* $r *$r)  /$r /$r;$jbb-=$Kyy;
  $Kzz= -(( $kL-$kT)* $rz * $rz + $kT* $r *$r)  /$r /$r;$jcc-=$Kzz;

  $Kxy= -($kL-$kT) * ( $rx * $ry) /$r /$r ;$jab-=$Kxy;
  $Kyz= -($kL-$kT) * ( $ry * $rz) /$r /$r ;$jbc-=$Kyz;
  $Kxz= -($kL-$kT) * ( $rx * $rz) /$r /$r ;$jac-=$Kxz;

  #if(abs($r*$r-$rx*$rx-$ry*$ry-$rz*$rz)>0.00001){die "inconsistent r\n";}
   # add also something to the Knn and Kmm in $sipffilethis and $sipffile !!! 
   addK($sipffilethis);
   addK($sipffile);
   # increase elastic constant $cel according to this bond
   # convert N/m in  meV/A^2: Nm=Joule  m=10^10A meV=1.60217646 . 10-22 joules
   # 1 N/m= 1 J/m^2 = meV /1.60217646e-2 A 
   $cL=$long_spring[$n]/1.60217646e-2;
   $cT=$trans_spring[$n]/1.60217646e-2;
# factor of 0.5 inserted, because each bond is called twice
# here when running makenn
   $Cel->slice('1,1')+=0.5*(($cL-$cT)*$rx*$rx*$rx*$rx/$r/$r+$cT*$rx*$rx);
   $Cel->slice('2,2')+=0.5*(($cL-$cT)*$ry*$ry*$ry*$ry/$r/$r+$cT*$ry*$ry);
   $Cel->slice('3,3')+=0.5*(($cL-$cT)*$rz*$rz*$rz*$rz/$r/$r+$cT*$rz*$rz);
   $Cel->slice('4,4')+=0.5*(($cL-$cT)*$ry*$rz*$ry*$rz/$r/$r+0.25*$cT*($ry*$ry+$rz*$rz));
   $Cel->slice('5,5')+=0.5*(($cL-$cT)*$rx*$rz*$rx*$rz/$r/$r+0.25*$cT*($rx*$rx+$rz*$rz));
   $Cel->slice('6,6')+=0.5*(($cL-$cT)*$rx*$ry*$rx*$ry/$r/$r+0.25*$cT*($ry*$ry+$rx*$rx));

   $Cel->slice('1,2')+=0.5*(($cL-$cT)*$rx*$rx*$ry*$ry/$r/$r);

   $Cel->slice('1,3')+=0.5*(($cL-$cT)*$rx*$rx*$rz*$rz/$r/$r);
   $Cel->slice('2,3')+=0.5*(($cL-$cT)*$ry*$ry*$rz*$rz/$r/$r);

   $Cel->slice('1,4')+=0.5*(($cL-$cT)*$rx*$rx*$ry*$rz/$r/$r);
   $Cel->slice('2,5')+=0.5*(($cL-$cT)*$ry*$ry*$rx*$rz/$r/$r);
   $Cel->slice('3,6')+=0.5*(($cL-$cT)*$rz*$rz*$rx*$ry/$r/$r);
   
   $Cel->slice('1,5')+=0.5*(($cL-$cT)*$rx*$rx*$rx*$rz/$r/$r+0.5*$cT*$rx*$rz);
   $Cel->slice('1,6')+=0.5*(($cL-$cT)*$rx*$rx*$rx*$ry/$r/$r+0.5*$cT*$rx*$ry);

   $Cel->slice('2,4')+=0.5*(($cL-$cT)*$ry*$ry*$ry*$rz/$r/$r+0.5*$cT*$ry*$rz);
   $Cel->slice('2,6')+=0.5*(($cL-$cT)*$ry*$ry*$rx*$ry/$r/$r+0.5*$cT*$ry*$rx);

   $Cel->slice('3,4')+=0.5*(($cL-$cT)*$rz*$rz*$ry*$rz/$r/$r+0.5*$cT*$rz*$ry);
   $Cel->slice('3,5')+=0.5*(($cL-$cT)*$rz*$rz*$rx*$rz/$r/$r+0.5*$cT*$rz*$rx);

   $Cel->slice('4,5')+=0.5*(($cL-$cT)*$ry*$rz*$rx*$rz/$r/$r+0.25*$cT*$rx*$ry);
   $Cel->slice('4,6')+=0.5*(($cL-$cT)*$ry*$rz*$rx*$ry/$r/$r+0.25*$cT*$rx*$rz);
   $Cel->slice('5,6')+=0.5*(($cL-$cT)*$rx*$rz*$rx*$ry/$r/$r+0.25*$cT*$rz*$ry);
             
 # now deal with the mixing term phonon-strain interaction constants Gmix
 # units should be meV (because in the phonon module the interaction operators 
 # I1=P1/a0 I2=P2/a0 I3=P3/a0 with a0=.5292e-10m=0.5292A
 # $cL and $cT above have units mev/A^2, i.e. we add a factor of 0.5292
   $cL=$long_spring[$n]*0.5292/1.60217646e-2;
   $cT=$trans_spring[$n]*0.5292/1.60217646e-2;
   # Voigt notation in first slot Gmix(11)1 
   $Gmix->slice('1,1')+=0.5*(-2*($cL-$cT)*$rx*$rx*$rx/$r/$r-$cT*($rx+$rx));
   # Gmix(11)2
   $Gmix->slice('1,2')+=0.5*(-2*($cL-$cT)*$rx*$rx*$ry/$r/$r);
   # Gmix(11)3
   $Gmix->slice('1,3')+=0.5*(-2*($cL-$cT)*$rx*$rx*$rz/$r/$r);
   # Gmix(22)1
   $Gmix->slice('2,1')+=0.5*(-2*($cL-$cT)*$ry*$ry*$rx/$r/$r);
   # Gmix(22)2
   $Gmix->slice('2,2')+=0.5*(-2*($cL-$cT)*$ry*$ry*$ry/$r/$r-$cT*($ry+$ry));
   # Gmix(22)3
   $Gmix->slice('2,3')+=0.5*(-2*($cL-$cT)*$ry*$ry*$rz/$r/$r);
   # Gmix(33)1
   $Gmix->slice('3,1')+=0.5*(-2*($cL-$cT)*$rz*$rz*$rx/$r/$r);
   # Gmix(33)2
   $Gmix->slice('3,2')+=0.5*(-2*($cL-$cT)*$rz*$rz*$ry/$r/$r);
   # Gmix(33)3
   $Gmix->slice('3,3')+=0.5*(-2*($cL-$cT)*$rz*$rz*$rz/$r/$r-$cT*($rz+$rz));

   # Gmix(23)1
   $Gmix->slice('4,1')+=0.5*(-2*($cL-$cT)*$ry*$rz*$rx/$r/$r);
   # Gmix(23)2
   $Gmix->slice('4,2')+=0.5*(-2*($cL-$cT)*$ry*$rz*$ry/$r/$r-$cT*$rz);
   # Gmix(23)3
   $Gmix->slice('4,3')+=0.5*(-2*($cL-$cT)*$ry*$rz*$rz/$r/$r-$cT*$ry);

   # Gmix(13)1
   $Gmix->slice('5,1')+=0.5*(-2*($cL-$cT)*$rx*$rz*$rx/$r/$r-$cT*$rz);
   # Gmix(13)2
   $Gmix->slice('5,2')+=0.5*(-2*($cL-$cT)*$rx*$rz*$ry/$r/$r);
   # Gmix(13)3
   $Gmix->slice('5,3')+=0.5*(-2*($cL-$cT)*$rx*$rz*$rz/$r/$r-$cT*$rx);

   # Gmix(12)1
   $Gmix->slice('6,1')+=0.5*(-2*($cL-$cT)*$rx*$ry*$rx/$r/$r-$cT*$ry);
   # Gmix(12)2
   $Gmix->slice('6,2')+=0.5*(-2*($cL-$cT)*$rx*$ry*$ry/$r/$r-$cT*$rx);
   # Gmix(12)3
   $Gmix->slice('6,3')+=0.5*(-2*($cL-$cT)*$rx*$ry*$rz/$r/$r);
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


 return ([$Gmix,$jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc]);  
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

      if ($r1x==0&&/^(#!|[^#])*r1a\s*=\s*/){($r1x)=extract("r1a",$_);}
      if ($r1y==0&&/^(#!|[^#])*r1b\s*=\s*/){($r1y)=extract("r1b",$_);}
      if ($r1z==0&&/^(#!|[^#])*r1c\s*=\s*/){($r1z)=extract("r1c",$_);}
      if ($r2x==0&&/^(#!|[^#])*r2a\s*=\s*/){($r2x)=extract("r2a",$_);}
      if ($r2y==0&&/^(#!|[^#])*r2b\s*=\s*/){($r2y)=extract("r2b",$_);}
      if ($r2z==0&&/^(#!|[^#])*r2c\s*=\s*/){($r2z)=extract("r2c",$_);}
      if ($r3x==0&&/^(#!|[^#])*r3a\s*=\s*/){($r3x)=extract("r3a",$_);}
      if ($r3y==0&&/^(#!|[^#])*r3b\s*=\s*/){($r3y)=extract("r3b",$_);}
      if ($r3z==0&&/^(#!|[^#])*r3c\s*=\s*/){($r3z)=extract("r3c",$_);}

      if ($r1x==0&&/^(#!|[^#])*r1x\s*=\s*/){($r1x)=extract("r1x",$_);}
      if ($r1y==0&&/^(#!|[^#])*r1y\s*=\s*/){($r1y)=extract("r1y",$_);}
      if ($r1z==0&&/^(#!|[^#])*r1z\s*=\s*/){($r1z)=extract("r1z",$_);}
      if ($r2x==0&&/^(#!|[^#])*r2x\s*=\s*/){($r2x)=extract("r2x",$_);}
      if ($r2y==0&&/^(#!|[^#])*r2y\s*=\s*/){($r2y)=extract("r2y",$_);}
      if ($r2z==0&&/^(#!|[^#])*r2z\s*=\s*/){($r2z)=extract("r2z",$_);}
      if ($r3x==0&&/^(#!|[^#])*r3x\s*=\s*/){($r3x)=extract("r3x",$_);}
      if ($r3y==0&&/^(#!|[^#])*r3y\s*=\s*/){($r3y)=extract("r3y",$_);}
      if ($r3z==0&&/^(#!|[^#])*r3z\s*=\s*/){($r3z)=extract("r3z",$_);}

      if ($nofatoms==0){($nofatoms)=extract("nofatoms",$_);}
      if ($nofcomponents==0){($nofcomponents)=extract("nofcomponents",$_);}


      if (/^(#!|[^#])*nofneighbours\s*=\s*/){++$n;

                                   ($nofneighbours[$n])=extract("nofneighbours",$_);
                                   ($diagonalexchange)=extract("diagonalexchange",$_);
                                   if($diagonalexchange>1){print "# Warning program makenn: diagonalexchange=$diagonalexchange not implemented setting diagonalexchange=0\n";$diagonalexchange=0;}
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
                                    if($module[$n]=~/ic1ion/&&$classdip){die "Error makenn in $sipffilename: for ic1ion module the classical dipole interaction is not implemented - you can use rkky or kaneyoshi\n";}
                                    if($module[$n]=~/icf1ion/&&$classdip){die "Error makenn in $sipffilename: for icf1ion module the classical dipole interaction is not implemented - you can use rkky or kaneyoshi\n";}
                                    if($module[$n]=~/phonon/&&$classdip){die "Error makenn in $sipffilename: for phonon module the classical dipole cannot be calculated - use option -bvk\n";}

                                    if($cfph!=0){($magnetic)=extractfromfile("MAGNETIC",$sipffilename); 
                                     if($magnetic!=0){unless($module[$n]=~/so1ion/||$module[$n]=~/ic1ion/||$module[$n]=~/icf1ion/){die "Error makenn:  ion $n in mcphas.j sipf=$sipffilename for option -cfph magnetic MODULES must be so1ion or ic1ion or icf1ion\n"; }
                                                      if($module[$n]=~/ic1ion/||$module[$n]=~/icf1ion/){unless($cfph==1||$cfph==2){die"Error makenn: for option -cfph mixing so1ion and ic1ion/icf1ion modules is not possible\n";}
                                                                                $cfph=2;} # for ic1ion and icf1ion module set it to 2 
                                                      if($module[$n]=~/so1ion/){unless($cfph==1||$cfph==3){die"Error makenn: for option -cfph mixing so1ion and ic1ion modules is not possible\n";}
                                                                                $cfph=3;} # for so1ion module set it to 3 
                                                      ++$nofmagneticatoms;$nph[$nofmagneticatoms+$nofatoms]=$n;
                                                      $sipf_file[$nofmagneticatoms+$nofatoms]=$sipffilename;
                                                      $gJ[$nofmagneticatoms+$nofatoms]=$gJ[$n];
                                                      $x[$nofmagneticatoms+$nofatoms]=$x[$n];
                                                      $y[$nofmagneticatoms+$nofatoms]=$y[$n];
                                                      $z[$nofmagneticatoms+$nofatoms]=$z[$n];
                                                      $charge[$nofmagneticatoms+$nofatoms]=$charge[$n];
                                                     }  }   
                                   if ($bvk){unless($module[$n]=~/phonon/){die "Error makenn in $sipffilename: option -bvk makes only sense for phonon module \n";}
                                             clearMP($sipffilename);
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
    unless (open($h,$filein)) {die "cannot open file $filein\n";}
#get comment lines from input file
    my $aaa=0;
     while(<$h>)
     {#next if /^\s*#/;
      $text=$_;
      if ($aaa==0){($aaa)=extract("a",$_);}
      last if($aaa!=0); # stop printing comments when crystallographic information starts
      print $l ($text);
     }
     print $l "#! a=".$a." b=".$b." c=".$c." alpha=".$alpha." beta=".$beta." gamma=".$gamma."\n";
     print $l "#! r1a=".$r1x." r2a=".$r2x." r3a=".$r3x."\n";
     print $l "#! r1b=".$r1y." r2b=".$r2y." r3b=".$r3y." primitive lattice vectors [a][b][c]\n";
     print $l "#! r1c=".$r1z." r2c=".$r2z." r3c=".$r3z."\n#\n";

     if ($aaa==0){die "ERROR makenn: unable to find 'a=' in file $filein\n";}
     
 return ($h,$l);
}



sub printneighbourlist {
  my ($Gmix,$h,$nofn,$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb)=@_;
     
     push @atoms,"#*************************************************************************\n";
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
      push @atoms,$text;
# the next lines are to advance $h to the end of the numeric table

      while(<$h>)
      {last if ($nn0==0);
       unless (/^\s*#/){--$nn0;}
      }



if ($bvk==1)
{push @atoms, "# it follows the Born von Karman model according to  springs read from file $bvk_file\n";
 push @atoms, "# the mixing terms Gmix in meV m\n";
 push @atoms, "#! Gindices=";for($i=1;$i<=6;++$i){for($j=1;$j<=3;++$j){push @atoms, " ".$i.",".$j;}}
 push @atoms, "\n";
 push @atoms, "#! G=";for($i=1;$i<=6;++$i){for($j=1;$j<=3;++$j){push @atoms, sprintf(" %+10.9e",$Gmix->at($i,$j));}}
 push @atoms, "\n";
}
elsif ($rkky==1)
{push @atoms, ("# it follows output of RKKY interaction according to J(R)=A.cos(2.kf.R)/(2.kf.R)^3 with A=$scale meV and kf=$kf A^-1 generated by makenn\n");}
elsif ($rkky==3)
{push @atoms, ("# it follows output of RKKY interaction according to J(R)=A.[sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 with A=$scale meV and kf=$kf A^-1 generated by makenn\n");}
elsif ($rkky==2)
{push @atoms, ("# kaneyoshi parametrization for the Bethe-Slater curve J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2] for scale A=$scale meV D=$D A alpha=$aa\n");}
elsif ($rkky==4)
{push @atoms, ("# it follows output of RKKY interaction J(R)=A.cos(2.kfR)/(2.kfR)^3 for scale A=$scale meV and\n");
 push @atoms, ("# kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n");}
elsif($rkky==5)
{push @atoms, ("# kaneyoshi parametrization for the Bethe-Slater curve J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2] for scale A=$scale meV\n");
 push @atoms, ("# with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2) with Da=$Da A Db=$Db A Dc=$Dc A and alpha=$aa\n");}
elsif($rkky==6)
 {push @atoms, ("# it follows output of RKKY interaction  J(R)=A [sin(2.kfR)-2.kfR.cos(2.kfR)]/(2.kfR)^4 for scale A=$scale meV\n");
  push @atoms, ("# kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n");}
elsif($readtable>0)
  {push @atoms, ("# it follows output generated from interactions read from table $table_file\n");
  }
else
{push @atoms, ("# it follows output of classical DD interaction generated by makenn\n");}

if ($djdx) {push @atoms, ("# - derivative with respect to x\n");}
if ($djdy) {push @atoms, ("# - derivative with respect to y\n");}
if ($djdz) {push @atoms, ("# - derivative with respect to z\n");}


    if($alpha!=90||$beta!=90||$gamma!=90)
     {push @atoms, ("#da[a]    db[b]     dc[c]       Jii[meV]  Jjj[meV]  Jkk[meV]  Jij[meV]  Jji[meV]  Jik[meV]  Jki[meV]  Jjk[meV]  Jkj[meV] with j||b, k||(a x b) and i normal to k and j\n");}
    else
     {push @atoms, ("#da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");}

my $pcout=">./results/makenn.a".$nnn.".pc";
unless (open($l1,$pcout)){die "cannot open file $pcout\n";}
print $l1 "#-------------------------------------------------------------------------------------\n";
print $l1 "#!  table with neighbors and charges for atom n=$nnn at da=".$x[$nnn]." db=".$y[$nnn]." dc=".$z[$nnn]." sipffilename=".$sipf_file[$nnn]."\n";
print $l1 "# output of program makenn:, Reference: M. Rotter et al. PRB 68 (2003) 144418\n";
print $l1 "# lattice [A]:\n";
# print $l1 "# $rtoijk\n";
print $l1 sprintf("#! ax= %+10.6f  ay=%+10.6f az=%+10.6f\n",$rtoijk->index(0)->at(0),$rtoijk->index(1)->at(0),$rtoijk->index(2)->at(0)); 
print $l1 sprintf("#! bx= %+10.6f  by=%+10.6f bz=%+10.6f\n",$rtoijk->index(0)->at(1),$rtoijk->index(1)->at(1),$rtoijk->index(2)->at(1)); 
print $l1 sprintf("#! cx= %+10.6f  cy=%+10.6f cz=%+10.6f\n",$rtoijk->index(0)->at(2),$rtoijk->index(1)->at(2),$rtoijk->index(2)->at(2)); 
print $l1 "#------------------------------------------------------------------------------------\n";
    if($alpha!=90||$beta!=90||$gamma!=90)
     {print $l1 "#orthonormal coordinate system xyz is defined with respect to abc as y||b, z||(a x b) and x normal to k and j\n#charge[|e|]  dx[A]   dy[A]   dz[A]        da[a]    db[b]    dc[c]   distance[A] atomnr\n";}
     else
     {print $l1 "#1            2         3         4              5          6          7         8             9\n";
      print $l1 "#charge[|e|]  da[A]     db[A]     dc[A]          da[a]      db[b]      dc[c]     distance[A]   atomnr\n";}

 for ($n1=1;$n1<(($rn->dims)[0]);++$n1)
 {next if($rn->index($n)->at($n1)==0);
  # the position xyz is relative position (not absolute coordinate of neighbour)
 push @atoms, sprintf("%+10.6f %+10.6f %+10.6f ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
 $ddd=$an->index($n)->at($n1);
  print $l1 sprintf("%8s   %+10.6f %+10.6f %+10.6f     ",$charge[$ddd],$in->index($n)->at($n1),$jn->index($n)->at($n1),$kn->index($n)->at($n1));
  print $l1 sprintf("%+10.6f %+10.6f %+10.6f ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
  print $l1 sprintf("%+10.6f     %s\n",$rn->index($n)->at($n1),$ddd);

  if (($classdip==1)||$bvk==1||($readtable>0&&$DM>0)||($readtable>0&&$Jp>0)) # here anisotropic interaction comes in
   {      push @atoms, sprintf("%+10.9e %+10.9e %+10.9e ",$Jaa->index($n)->at($n1),$Jbb->index($n)->at($n1),$Jcc->index($n)->at($n1));
          for($i=4;$i<=$nofcomponents;++$i){push @atoms, "0 ";} # add other components
          push @atoms, sprintf("%+10.9e %+10.9e ",$Jab->index($n)->at($n1),$Jba->index($n)->at($n1));
          push @atoms, sprintf("%+10.9e %+10.9e ",$Jac->index($n)->at($n1),$Jca->index($n)->at($n1));
          for($i=4;$i<=$nofcomponents;++$i){push @atoms, "0 0 ";} # add other components
          push @atoms, sprintf("%+10.9e %+10.9e ",$Jbc->index($n)->at($n1),$Jcb->index($n)->at($n1));
          for($i=4;$i<=$nofcomponents;++$i){for($ii=$i+1;$ii<=$nofcomponents;++$ii){push @atoms, "0 0 ";}} # add other components         
   }
   else  #here the isotropic interaction is written
   {push @atoms, sprintf("%+10.9e %+10.9e %+10.9e ",$Jaa->index($n)->at($n1),$Jbb->index($n)->at($n1),$Jcc->index($n)->at($n1));
    for($i=4;$i<=$nofcomponents;++$i){push @atoms, "0 ";} # add other components
    
   }

 if ($calcdist==1) {push @atoms, sprintf("%+10.9e a%s",$rn->index($n)->at($n1),$ddd);}

  push @atoms, "\n";

 }

 close $l1;
}


sub outCel {
my ($l)=@_;
print $l "#! Primitive Unit Cell Volume [A^3]: pVol=".$pVolume."\n";
print $l "# Nonzero Elastic constants [meV per primitive crystal unit cell] \n";
print $l "# in Voigt notation only first index<=second index has to be given\n";
print $l "# because the constants are symmetric Celij=Celji\n";
print $l "# Elastic constants refer to the Euclidean coordinate system ijk defined\n";
print $l "# with respect to abc as j||b, k||(a x b) and i normal to k and j\n";

# here output the elastic constants table
$i1=0;print $l "#! ";
for($i=1;$i<=6;++$i){ if($i1>0){print $l "#! ";}
for($j=$i;$j<=6;++$j){
if(abs($Cel->at($i,$j))>1e-6){++$i1;print $l sprintf(" Cel%i%i=%+10.9g",$i,$j,$Cel->at($i,$j));}
                     }if($i1>0){print $l "\n";}}
if($i1==0){print $l "\n";}

# 1meV= 1.60218e-22 J
# 1 A= 1e-10 m
$confact=0.1*1.60218/$pVolume;
print $l "#! unit conversion:  1 meV/Primitive Unit Cell Volume =".$confact." GPa\n";
}

sub endprint {
  my ($h,$l)=@_;  

     close $h; 

# print $l ("#*************************************************************************\n");
 outCel($l);
if($bvk){my $l1= new FileHandle;open($l1,">results/makenn.Cel"); outCel($l1);close $l1;}
$text="#! nofatoms= $nofatoms  nofcomponents=$nofcomponents  number of atoms in primitive unit cell/number of components of each spin\n";
 if($cfph!=0){$nofatomsnew=$nofatoms+$nofmagneticatoms;}
 if($cfph==2){$text="#! nofatoms= $nofatomsnew  nofcomponents=51  number of atoms in primitive unit cell/number of components of each spin\n";}
 if($cfph==3){$text="#! nofatoms= $nofatomsnew  nofcomponents=48  number of atoms in primitive unit cell/number of components of each spin\n";}
 print $l ($text);

     print $l @atoms;
     close $l;


}
# rename file
sub my_rename {
my ($s,$d)=@_;
unless (rename $s,$d)
          {unless(open (Fout, ">$d"))     
      {die "\n error:unable to write to $d\n";}
      open (Fin, $s);
      while($line=<Fin>){ print Fout $line;}
      close Fin;
      close Fout;
      system "del $s"; 
     }

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
# sets variable in file
# 
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to set this number  to 0.20 in this file just use:         
#
# setvariable("sta",0.20,"data.dat");
# 
# 
#
sub setvariable {
my ($varnam,$value,$file)=@_;
my $line;   
unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   open ( Fout, ">range.out");$ri=0;
  while($line=<Fin>)
     { if ($line=~/^(#!|[^#])*?\b$varnam\s*=/) {
                                            #here write modified parameter set to line
       # $line=~s/(^(#!|[^#])*?\b)$varnam[ \t]*=[ \t]*([^ \;\n\r\t\*]*)([ \;\n\r\t\*])/$1$varnam=$value$4/g;
       $line=~s/$varnam[ \t]*=[ \t]*([^ \;\n\r\t\*]*)([ \;\n\r\t\*])/$varnam=$value$2/g;
                                               }
       print Fout $line;
      }
      close Fin;
      close Fout;

     unless (rename "range.out",$file)
     {unless(open (Fout, ">$file"))     
      {die "\n error:unable to write to $file\n";}
      open (Fin, "range.out");
      while($line=<Fin>){ print Fout $line;}
      close Fin;
      close Fout;
      system "del range.out"; 
     }
}


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
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\EA-Za-z]+)/);}}
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\" to extract variable $variable\n";
             }
             return $value;
            }
# **********************************************************************************************
sub usage ()
{
print STDOUT << "EOF";

 usage: makenn 23.3 [options] 

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


 option -rkky A(meV) kf(1/A) calculates the rkky interaction
              according to J(R)=A.cos(2.kf.R)/(2.kf.R)^3
              scaling A<0, kf should be the Fermi wavevector (usually
              between 0.3-2.5 A^-1 depending on the electrondensity^0.333)
 option -rkky3d A(meV) ka(1/A) kb(1/A) kc(1/A) calculates the rkky interaction
              according to J(R)=A.cos(2.kfR)/(2.kfR)^3
              scaling A<0, kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2)
 option -rkkz A(meV) kf(1/A) calculates the rkky interaction
              according to J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4
              scaling A>0, kf should be the Fermi wavevector
 option -rkkz3d A(meV) ka(1/A) kb(1/A) kc(1/A)  calculates the rkky interaction
              according to J(R)=A [sin(2.kfR)-2.kfR.cos(2.kfR)]/(2.kfR)^4
              scaling A>0, kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2)
 option -kaneyoshi A(meV) D(A) alpha  calculates the kaneyoshi
             parametrization for the Bethe-Slater
              curve: J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2]
              with D corresponding to the orbital radius
              the exponential alpha is conveniently put to  about 1
 option -kaneyoshi3d A(meV) Da(A) Db(A) Dc(A) alpha  calculates the 3d-kaneyoshi
             parametrization for the Bethe-Slater
              curve: J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2]
              with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2)
              the exponential alpha is conveniently put to  about 1
 option -bvk filename
              for phonons: take Born van Karman model with longitudinal and
              transversal spring constants from file - file format, columns:
              #   atom_n_sipf atom_n'_sipf bondlength(A) Clong(N/m) Ctrans(N/m)
              mind: into MODPAR2-6 in *.sipf the Einstein-oscillator paramters 
              are written, too. Omit filename to create a sample file with
              longitudinal springs:Clong=$bvkA*exp(-$bvkalpha*r/A*r/A) N/m
	      Output: file makenn.Cel is created containing just the elastic constants

 option -cfph [screeningfile.r] [-r]
              calculate crystal field phonon interaction: mcphas.j lists 
              magnetic and non magnetic atoms with charges defined in the 
              sipf files by CHARGE= variable. For magnetic atoms the sipf 
              file the variable MAGNETIC=1 has to be set and information 
              about the ion has to be present (IONTYPE etc.). 
              Foreach magnetic ion a new site is created resembling the magnetic 
	      electron charge cloud and this new site is shifted
              0.1 A along c in order to not overlap with the original site.
              For option -r the original magnetic site sipf is replaced automatically
	      to use an sipf file with the MODULE=phonon  similar to all the other
              nonmagnetic sites (which should have a PHONON module). 
	      For the new magnetic site the program
              pointc is used by makenn with option -d to calculate derivatives
              dBlm/du which are inserted as interaction 
              parameters between MODULE=phonon and MODULE=so1ion sites.
              The new magnetic sites single ion property files are named
              results/makenn.a*.sipf and contain crystal field paramters Blm 
              calculated by pointc.
              In order to use the resulting file results/makenn.j a phonon
              model has to be set up, 
              and the phonon model has to be added to makenn.j,  e.g. by
              program addj. Moreover magnetic electron sites sipf files are required, 
              which are created in results/makenn.a*.sipf. 
              Note: a screening file can be used to define distance dependent 
              screening of charges for the pointcharge model calculation
              format: col1 distance r (Angstroem) col 2 screening factor 
              for B2m, col 3 for B4m and col 4 for B6m
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

        -d      puts to the last column the distance of the neighbors (A)
     The neigbours of each atom are stored in separate files
 results\/makenn.a*.pc, which can be used with the program pointc to evaluate
 the pointcharge model and calculate crystal field parameters.
	-npc   do not keep results/makenn.a*.pc files
	-rfunc  create results/makenn.func.r file with distance dependence of
                rkky and kaneyoshi functions (only for options rkky* and kaneyoshi*)
        -djdx
        -djdy
        -djdz   create files makenn.djdx .djdy .djdz instead of makenn.j, respectively
 these contain the derivatives of the interaction paramters with respect to
 displacement of the neighbor in x,y and z direction respectively. These derivatives
 are useful for the calculation of exchange striction effects. 

 xyz refers to a right handed Euclidean coordinate system with 
 y||b, z||(a x b) and x perpendicular to y and z.

EOF
}

# Rounds up a value, using int() function
sub my_ceil
{
    my $t = int($_[0]);
    return ($t != $_[0]) ? ($_[0] < 0) ? $t : $t + 1 : $_[0];
}

# Rounds down a value, using int() function
sub my_floor
{
    return ($_[0] < 0) ? int($_[0]) - 1 : int($_[0]);
}



sub store_phonon_file
{my ($filename)=@_;
 ($IONTYPE)=extractfromfile("IONTYPE",$filename);
 ($CHARGE)=extractfromfile("CHARGE",$filename);
 ($MAGNETIC)=extractfromfile("MAGNETIC",$filename);
 ($SCATTERINGLENGTHREAL)=extractfromfile("SCATTERINGLENGTHREAL",$filename);
 ($SCATTERINGLENGTHIMAG)=extractfromfile("SCATTERINGLENGTHIMAG",$filename);
my $at=$IONTYPE;$at =~ s/_//g; $at =~ s/[0-9]+[A-Z]+//g; $at =~ s/[0-9]+//g; $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g;
 $at = lc $at; $at = ucfirst $at; 
 unless($MASS=${$mass{$at}}[0]){print "WARNING: $at mass not found in internal table for new file $filename - please edit and set mass by hand \n"; }
print "Recreating $filename and making it a sipf file with phonon module with mass $at $MASS\n";
 
 open (FOUT, ">$filename");
 print FOUT "#!MODULE=phonon\n";
    print FOUT $sipfheader;
    print FOUT "IONTYPE=$IONTYPE\n";
    print FOUT "CHARGE=$CHARGE\n";
    print FOUT "MAGNETIC=$MAGNETIC\n";
    print FOUT "\n";
print FOUT << "EOF";
MODPAR1=$MASS  #mass in(m0)
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
            #                6       umax restriction in y and z direction
EOF
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)\n";
    print FOUT "#                      with s=sin(theta)/lambda=Q/4pi\n";
    print FOUT "# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>\n";
    print FOUT "# unit of DWF is [A^2]\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "DWF=0\n";
    print FOUT "\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Neutron Scattering Length (10^-12 cm) (can be complex)\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "SCATTERINGLENGTHREAL=$SCATTERINGLENGTHREAL\n";
    print FOUT "SCATTERINGLENGTHIMAG=$SCATTERINGLENGTHIMAG\n";
    print FOUT "#  ... note: - if an occupancy other than 1.0 is needed, just reduce \n";
    print FOUT "#              the scattering length linear accordingly\n";
close FOUT;
}
 
sub print_time_estimate_until_end #input :ratio = nofpointstodo / nofpointsdone
{my ($ratio)=@_;

 if(!defined $last_call_time){$last_call_time=time();$startcputime=$last_call_time;}
 if ((time()-$last_call_time)>60)# print only if 1 minutes passed since last print
   { # estimate time until end 
    my $time_since_start=time()-$startcputime;
    my $time_until_end=$time_since_start*$ratio; 
       
    my $curtime=time()+$time_until_end;my $loctime=localtime($curtime);
    print STDERR "- end ~at $loctime";
    my $min=int($time_until_end/60);
    if($min<2){
    print STDERR " in $time_until_end s";
    }else
    {my $hours=int($min/60);
     if($hours<2) 
     {print STDERR  " in $min min";
     } else
     {my $days=int($hours/24);
      if($days<2)
      { print STDERR " in $hours h";
      } else
      {my $months=int($days/30);
      if($months<24)
      {print STDERR  " in $months months";}
      else {print STDERR " in ".int($months/12)." years";}
     }}}
    $last_call_time=time();print STDERR "\n";
   }
}
