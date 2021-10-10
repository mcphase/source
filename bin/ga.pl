#!/usr/bin/perl

use FileHandle;

#use Math::Trig;
use PDL;


print "#********************************************************\n";
print "# program  ga 210930 \n";
print "#********************************************************\n";
$PI=3.141592654;

unless ($#ARGV>=0) 

{
print STDOUT << "EOF";

#   use as:
#    
#     ga  rmax atomnr qh qk ql n 
#
#     rmax  .... maximum distance rmax [Angstroem] of neighbour to take into account neighbour
#     atomnr .... crystal field parameters of this atom will be considered (atomnr refers to
#                 the list in mcphas.j)
#     (qh qk ql) ... components of the phonon q-vector in terms of reciprocal lattice
#                    
#     n     ... n frames of the phonon oscillation will be considered
#     
#     The eigenvector of the phonon should be stored in the sipf files
#     with the components 
#     PHONEVRX PHONEVRY PHONEVRZ PHONEVIX PHONEVIY PHONEVIZ
#     where X,Y and Z refers to a cartesian coordinate system XYZ such
#     that Y||b, Z||axb, X perpendicular to Y and Z
#
#  1) reads mpchas.j and gets  charges and eigenvector components (PHONEVR* PHONEVI*) from .sipf files
#      MIND: actually PHONEVRX etc should contain the normalized eigenvector 
#            components divided by the sqrt(M), where M is the mass-number of the nucleus
#  2) calculates for one magnetic ion the pointcharge positions and outputs ga0.pc which
#     resembles makenn.a1 (positions and point charges), but has in addition also
#     additional columns with unit cell vector 
#  3) calculate also ga1.pc, in such a way: take from command line the q-vector (qh qk ql) and calculate
#     from this and eigenvector components the oscillation of all atoms 
#     (for n different phases, n read from the command line), divide charge of each atom by n
#     and put into ga1.pc this charge and relative position to the one magnetic ion, such that
#     running pointc on this file will give the crystal field parameters as modified when the 
#     phonon oscillation is excited.
#  4) calculate a file ga.pc, which is a list taken from the difference of (3) and (2),
#     i.e. using this file with pointc will give the CF-phonon coupling energy for a phonon 
#     with amplitude  corresponding to the eigenvector components taken in Angstroem

#  Theory:
# Sj ... deviation Vector of atom position from equilibrium position 
#       (due to phononic vibration)
#
# Hcf_phon= sum_ijralpha  dBalpha(i)/dSjr Sjr Oalpha(i)      (1) 
#
# from Nowotny Festkoerperphysik p. 138 (r=nui index of atom in unit cell & euclidean coordinate)
# 
# Sjr= 1/sqrt(N) sum_q,s  sqrt[hbar/(2 Mr omega_q,s)] exp(i q.Rj) x
#                        x er_q,s ( a^cross_-qs  + a_qs ) 
#
#  with j primitive lattice vector index, N number of atoms, 
#  s phonon branch, q vector, Mr nuclear mass, omega_qs phonon frequency
#
# inserting in (1) yields
#
# Hcf_phon= 1/sqrt(N) sum_ijralphaqs  dBalpha(i)/dSj sqrt[hbar/(2 Mr omega_q,s)] exp(i q(Rj-Ri)) x
#                        x er_q,s ( a^cross_-qs  + a_qs ) Oalpha(i) exp(i qRi)
#               
# compare with Becker 2021  
#
# Hcf_phon= -1/sqrt(N) sum_iqsalpha  g_alpha(qs) ( a^cross_-qs  + a_qs ) Oalpha(i) exp(i qRi)
#
# we identify
#
#  g_alpha(qs)= - sum_j  dBalpha(i)/dSjr er_q,s  sqrt[hbar/(2 Mr omega_q,s)] exp(i q(Rj-Ri))
#
#
# for the file g_alpha.s we need the absolute value |g_alpha(qs)|
# because only this enters into the phononic cross section
#
# in order to evaluate dBalpha(i)/dSjr the program pointc can be used with the
# option -d, creating the file results/pointc.dBlm with
# position and dBalpha/duj  where duj=dSj/a0 
#
# 
EOF
 exit 0;}

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($rmax) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($atomnr) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($qh) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($qk) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($ql) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($nfi) = eval $ARGV[0];


my ($latt,$p) = getlattice("./mcphas.j"); # gets lattice and atomic positions
my ($a,$b,$c,$alpha,$beta,$gamma,$nofatoms,$nofcomponents) = @{$latt};
 print "# cf- phonon coupling calculated for atom atomnr=".$atomnr."with sipffile ".$sipf_file[$atomnr]."\n";
 print "# a=".$a." b=".$b." c=".$c." alpha=".$alpha." beta=".$beta." gamma=".$gamma."\n";
 print "# primitive lattice[abc]:".$p."\n";

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
if ($rtoijk->at(2,2)<=0){die "ERROR ga: alpha beta and gamma geometrically inconsistent\n";}
$t.=sqrt($t);
#print $t;
#print $rtoijk;
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
 print "primitive lattice[A]:".$p."\n";
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

 
# determine $nmin,$nmax by looking at a cube with side 3rmax
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
  
print "$n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";

     #initialize output file results/makenn.j
print "number of atoms = $nofatoms\n calculating ...\n";
 for ($nnn=$atomnr;$nnn<=$atomnr;++$nnn)    
 { print "atom $nnn ...\n";
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

     my ($inp)=new PDL();
     my ($jnp)=new PDL();
     my ($knp)=new PDL();
      
     my ($NNA)=new PDL();
     my ($NNB)=new PDL();
     my ($NNC)=new PDL();
     
     my ($Nr1)=new PDL();
     my ($Nr2)=new PDL();
     my ($Nr3)=new PDL();

  for ($n1=$n1min;$n1<=$n1max;++$n1){ 
  for ($n2=$n2min;$n2<=$n2max;++$n2){ 
  for ($n3=$n3min;$n3<=$n3max;++$n3){  
  for ($nz=1;$nz<=$nofatoms;++$nz){  
   $dabc=pdl [($x[$nz]-$x[$nnn]),($y[$nz]-$y[$nnn]),($z[$nz]-$z[$nnn])];
   $rvec= $dabc x $rtoijk;$rvec=$rvec->slice(":,(0)");
   $rvec+=$n1*$p->slice(",(0)")+$n2*$p->slice(",(1)")+$n3*$p->slice(",(2)");
   $Rnvec=$n1*$p->slice(",(0)")+$n2*$p->slice(",(1)")+$n3*$p->slice(",(2)");

   $rr=inner($rvec, $rvec);
   $r=sqrt($rr->at());
   $aabbcc=$rvec x $invrtoijk;$aabbcc=$aabbcc->slice(":,(0)");
   $xx=$aabbcc->at(0);
   $yy=$aabbcc->at(1);
   $zz=$aabbcc->at(2);
   $NABC=$Rnvec x $invrtoijk;$NABC=$NABC->slice(":,(0)");
   $NA=$NABC->at(0);
   $NB=$NABC->at(1);
   $NC=$NABC->at(2);

   if ($r<=$rmax && $r>0){#save neighbour j format
          
    $an=$an->append( pdl ([$nz]));
    $rn=$rn->append( pdl ([$r]));
    $xn=$xn->append( pdl ([$xx]));
    $yn=$yn->append( pdl ([$yy]));
    $zn=$zn->append( pdl ([$zz]));
    $in=$in->append( pdl ([$rvec->at(0)]));
    $jn=$jn->append( pdl ([$rvec->at(1)]));
    $kn=$kn->append( pdl ([$rvec->at(2)]));
    $NNA=$NNA->append( pdl ([$NA]));
    $NNB=$NNB->append( pdl ([$NB]));
    $NNC=$NNC->append( pdl ([$NC]));
    $Nr1=$Nr1->append( pdl ([$n1]));
    $Nr2=$Nr2->append( pdl ([$n2]));
    $Nr3=$Nr3->append( pdl ([$n3]));

                    
for ($nn1=0;$nn1<$nfi;++$nn1){  
$fi=$nn1*2*$PI/$nfi;

#dx=real(ex exp(i (Q R)+i fi))=
#  =real[(exr+i exi)(cos(QR+ fi)+i sin(QR+ fi)) ]
#  =exr cos(QR+ fi) - exi sin (QR+ fi)
#
# QR= (qh, qk, ql) . (Ra,Rb,Rc)=2pi (h da + k db + l dc)
    $QR=2*$PI*($qh*$xx+$qk*$yy+$ql*$zz);
    $inp=$inp->append( pdl ([$rvec->at(0)+($evrx[$n]-$evrx[$nnn])*cos($QR+$fi)-($evix[$n]-$evix[$nnn])*sin($QR+$fi)]));
    $jnp=$jnp->append( pdl ([$rvec->at(1)+($evry[$n]-$evry[$nnn])*cos($QR+$fi)-($eviy[$n]-$eviy[$nnn])*sin($QR+$fi)]));
    $knp=$knp->append( pdl ([$rvec->at(2)+($evrz[$n]-$evrz[$nnn])*cos($QR+$fi)-($eviz[$n]-$eviz[$nnn])*sin($QR+$fi)]));
                               }
                                     
                      }
    }}}}  

   $n= qsorti($rn); 
   $nofneighbours[$nnn]=(($rn->dims)[0]-1);
  
   printneighbourlist(">./results/ga0.pc","for equilibrium position",$nofneighbours[$nnn],$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$NNA,$NNB,$NNC,$Nr1,$Nr2,$Nr3);
   printneighbourlistp(">./results/ga1.pc","for excited phonon with q=($qh,$qk,$ql)",$nofneighbours[$nnn],$n,$an,$rn,$inp,$jnp,$knp,$NNA,$NNB,$NNC);
   printneighbourlistp(">./results/ga.pc","for excited phonon with q=($qh,$qk,$ql)",$nofneighbours[$nnn],$n,$an,$rn,$inp,$jnp,$knp,$NNA,$NNB,$NNC);
for ($n1=1;$n1<=$nofatoms;++$n1){$charge[$n1]=-1.0*$charge[$n1];}
   printneighbourlist(">>./results/ga.pc","for coupling parameters to phonon q=($qh,$qk,$ql)",$nofneighbours[$nnn],$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$NNA,$NNB,$NNC,$Nr1,$Nr2,$Nr3);

 }


  

 print "created files: results/ga*.pc (pointcharge environment files)\n";
print "********************************************************\n";
print "               end of program ga\n";
print "********************************************************\n";

   exit;

#-----------------------------------------------------------------------
sub printneighbourlistp {
  my ($file,$comment,$nofn,$n,$an,$rn,$in,$jn,$kn,$NNA,$NNB,$NNC)=@_;
     if ($nofn=="-1"){$nofn="0";}

my $pcout="$file";
unless (open($l1,$pcout)){die "cannot open file $pcout\n";}
print $l1 "#-------------------------------------------------------------------------------------\n";
print $l1 "#  table with neighbors and charges for atom $nnn\n";
print $l1 "# output of program ga $comment\n";
print $l1 "#-------------------------------------------------------------------------------------\n";
    if($alpha!=90||$beta!=90||$gamma!=90)
     {print $l1 "#orthonormal coordinate system ijk is defined with respect to abc as j||b, k||(a x b) and i normal to k and j\n#charge[|e|]  di[A]   dj[A]   dk[A]        atomnr\n";}
     else
     {print $l1 "#charge[|e|]  da[A]     db[A]     dc[A]            atomnr\n";}

 for ($n1=1;$n1<(($rn->dims)[0]);++$n1)
 {next if($rn->index($n)->at($n1)==0);
  $ddd=$an->index($n)->at($n1);
  for($nn1=0;$nn1<$nfi;++$nn1)
{ $np=($n->at($n1)-1)*$nfi+$nn1+1;
 print $l1 sprintf("%8s   %+10.6f %+10.6f %+10.6f     ",$charge[$ddd]/$nfi,$in->at($np),$jn->at($np),$kn->at($np));
 print $l1 sprintf("    %s\n",$ddd);
}
   
 }

 close $l1;

}



sub printneighbourlist {   
  my ($file,$comment,$nofn,$gJ,$n,$an,$rn,$xn,$yn,$zn,$in,$jn,$kn,$NNA,$NNB,$NNC,$Nr1,$Nr2,$Nr3)=@_;
     if ($nofn=="-1"){$nofn="0";}

my $pcout="$file";
unless (open($l1,$pcout)){die "cannot open file $pcout\n";}
print $l1 "#-------------------------------------------------------------------------------------\n";
print $l1 "#  table with neighbors and charges for atom $nnn\n";
print $l1 "# output of program ga $comment\n";
print $l1 "#-------------------------------------------------------------------------------------\n";
    if($alpha!=90||$beta!=90||$gamma!=90)
     {print $l1 "#orthonormal coordinate system ijk is defined with respect to abc as j||b, k||(a x b) and i normal to k and j\n#charge[|e|]  di[A]   dj[A]   dk[A]        da[a]    db[b]    dc[c]   distance[A] atomnr\n";}
     else
     {print $l1 "#charge[|e|]  da[A]     db[A]     dc[A]          da[a]      db[b]      dc[c]     distance[A]   atomnr na nb nc nr1 nr2 nr3\n";}

 for ($n1=1;$n1<(($rn->dims)[0]);++$n1)
 {next if($rn->index($n)->at($n1)==0);
  # the position xyz is relative position (not absolute coordinate of neighbour)
# print $l sprintf("%+10.6f %+10.6f %+10.6f ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
 $ddd=$an->index($n)->at($n1);
  print $l1 sprintf("%8s   %+10.6f %+10.6f %+10.6f     ",$charge[$ddd],$in->index($n)->at($n1),$jn->index($n)->at($n1),$kn->index($n)->at($n1));
  print $l1 sprintf("%+10.6f %+10.6f %+10.6f ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
  print $l1 sprintf("%+10.6f     %s %+10.6f %+10.6f %+10.6f  ",$rn->index($n)->at($n1),$ddd,$NNA->index($n)->at($n1),$NNB->index($n)->at($n1),$NNC->index($n)->at($n1));
  print $l1 sprintf("%+10.6f %+10.6f %+10.6f  \n",$Nr1->index($n)->at($n1),$Nr2->index($n)->at($n1),$Nr3->index($n)->at($n1));
 

   
 }

 close $l1;

}








# Get lattic data, reading it from file 

sub getlattice {
    my ($file) = @_;
    my $h = new FileHandle;
    my $n = 0;
    $nofcomponents=0;
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
                                   ($evrx[$n])=extractfromfile("PHONEVRX",$sipffilename);
                                     ($evix[$n])=extractfromfile("PHONEVIX",$sipffilename);
                                   ($evry[$n])=extractfromfile("PHONEVRY",$sipffilename);
                                     ($eviy[$n])=extractfromfile("PHONEVIY",$sipffilename);
                                   ($evrz[$n])=extractfromfile("PHONEVRZ",$sipffilename);
                                     ($eviz[$n])=extractfromfile("PHONEVIZ",$sipffilename);
                                     ($gJ[$n])=extractfromfile("GJ",$sipffilename);                       
                                     $sipf_file[$n]=$sipffilename;
				  }

     }

     close $h; 

     if ($n!=$nofatoms) {print STDOUT "Failed to read data file \"$file\": wrong number of atoms\n";

                         return undef;}
     if ($alpha<=0) { die "ERROR ga: reading unit cell angle alpha=$alpha <=0\n";}
     if ($beta<=0) { die "ERROR ga: reading unit cell angle beta=$beta <=0\n";}
     if ($gamma<=0) { die "ERROR ga: reading unit cell angle gamma=$gamma <=0\n";}

     if($nofcomponents==0) {$nofcomponents=3;}

     return ([$a,$b,$c,$alpha,$beta,$gamma,$nofatoms,$nofcomponents],pdl [[$r1x,$r1y,$r1z],[$r2x,$r2y,$r2z],[$r3x,$r3y,$r3z]]);

    } else {
	print STDOUT "Warning: failed to read lattice from file \"$file\"\n";
	return undef;
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
