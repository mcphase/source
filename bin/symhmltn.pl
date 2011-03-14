#!/usr/bin/perl
#
# symhmltn.pl - Determines the symmetry allowed multipolar order terms in the Hamiltonian given a point symmetry
#
# This program is part of the McPhase package, licensed under the GNU GPL v2. Please see the COPYING
# file
#
# Tue Dec 21 10:30:03 WEST 2010 - Duc Le and Martin Rotter

# We want to calculate the linear combinations of multipolar operator products (representing the
# exchange interactions) which are invariant under all operations of a particular point group.
# To do this we:
#   1. Take the multipolar operators as transforming according to the group (T^l)x(T^l) where l denotes the
#      multipolar order (1==dipole, 2==quadrupole, etc). and the group T transform in the same way as 
#      spherical harmonics.
#   2. We decompose the product (T^l)x(T^l) into IRREPS of the point group of interest, from the character of
#      the group T^l which is sin((l+1/2)alpha)/sin(alpha/2) where alpha is the angle or rotation for a particular
#      operator of the point group of interest. We take mirror planes to have alpha=pi. 
#   3. Find the number of times the unit IRREP of the point group is included in the decomposition of 
#      (T^l)x(T^l). This is the number of independent terms in the Hamiltonian.
#   4. Determine a projection operator that projects the basis states of (T^l)x(T^l) into a subspace which
#      transforms according to the unit IRREP. This is given by:
#             P^{A1} = (1/order) sum_G { chi(G).T(G) }
#      where G labels the group elements. The operators T(G) are determined in matrix form from the rotation 
#      operators given by Buckmaster, phys. stat. sol., v13, p9, 1972. chi(G) is the character of the element.
#      The projection operator must be hermitian, and idempotent, e.g. P=P' and P^n=P, and is checked in the code.
#   5. Operate this projection operator on the basis states of T^l systematically until the required number of 
#      independent terms is obtained.

$debug=1;
$SMALL=1e-3;

require 'ptgptabs.pl';
require 'rotstev.pl';
use Math::Trig;
use Getopt::Long;

GetOptions("help"=>\$helpflag,
           "buckmaster"=>\$usebuckop);

$usestevop=!$usebuckop;

sub sign { if(@_[0]>=0) { return 1; } else { return -1; } }
sub round { return int(@_[0] + sign(@_)*0.5); }

if ($#ARGV<1 || $helpflag) { 
   print " $0 - program to determine the symmetry allowed multipolar terms in\n";
   print "                a Hamiltonian given a point symmetry.\n\n";
   print " Syntax: $0 <PTGP> <l>\n\n";
   print " where <PTGP> is a (case insensitive) string denoting the point group.\n";
   print "              Both Hermman-Maguin and Schoenflies symbols are allowed.\n";
   print "       <l> is the multipolar order.\n";
   print "              (e.g. 2 for quadrupoles, 3 for octupoles, etc.)\n";
   print " By default the program outputs the allowed multipolar Hamiltonian in\n";
   print " terms of Stevens operators. If a flag -b is given, Buckmaster/Wybourne\n";
   print " operators (transforming like spherical harmonics, rather than tesseral\n";
   print " harmonics like the Stevens operators) are used instead.\n";
   exit(0);
}

# Determines if the input point group is in HM or Schoenflies notation, and converts to Schoenflies
$ptgpsym=$ARGV[0]; $J=$ARGV[1];
if ($ptgpsym=~m/m/ || $ptgpsym!~/[CcDdSsTtOoV]/) { $ptgpsym=hm2schoenflies($ARGV[0]); }
if ($ptgpsym!~/^[VS]/) { $ptgpsym =~ s/c/C/; $ptgpsym =~ s/d/D/; $ptgpsym =~ s/o/O/; $ptgpsym =~ s/t/T/;
$ptgpsym =~ s/H/h/; $ptgpsym =~ s/V/v/; $ptgpsym =~ s/D$/d/; $ptgpsym =~ s/I/i/; $ptgpsym =~ s/S/s/; }

if ($debug) { print $ptgpsym,"\n"; }

# Gets the rotation angles, operators and character table for this group
for $iop (0..($ptgp{$ptgpsym}[0]-1)) { push(@ptops,$ptgp{$ptgpsym}[$iop+2]); push(@angs,$rang{$ptgpsym}[$iop]); }
$order = $ptgp{$ptgpsym}[1]; %mychartab = %{ $chartab{$ptgpsym} };
@irreps = keys %mychartab;

if ($debug) { print "order = $order\n"; }
if ($debug) { print join("\t",@irreps),"\n"; }
if ($debug) { print "\t"; foreach (@angs) { print $_/pi,"\t"; } print "\n"; }
if ($debug) { print "\n\t",join("\t",@ptops),"\n"; }
if ($debug) { foreach(@irreps) { print $_,"\t",join("\t",@{ $mychartab{$_} }),"\n"; } }

# Calculates the character of the Rotation group in 2J+1 dimensions and the class multiplicities
for $ialpha (0..$#angs) {
   $alpha = $angs[$ialpha]; $r_op=$ptops[$ialpha];
   if($r_op=~/[Ei]/) {
      push(@charR, 2*$J+1); }
   elsif($alpha==0) {
      push(@charR, 0); }
   else {
      push(@charR, (sin(($J+0.5)*$alpha)/sin($alpha/2))); }
   if (abs($charR[-1])<$SMALL) { $charR[-1] = 0; }
   $r_op =~ /^(.)/; $cm=$1; $cm =~ s/[CSsEi]/1/; 
   push(@classmult, $cm); 
}

if ($debug) { print "\t",join("\t",@classmult),"\n"; }
if ($debug) { print "R_J\t",join("\t",@charR),"\n\n"; }

# Determines the unit representation.
for $i_irrep (0..$#irreps) { 
   push(@cmpfl,0); $unitfl=1;
   # Checks if the character is complex, we need to include the complex conjugate terms as well...
   for $i_ops (0..$#ptops) { 
      if(Im($mychartab{$irreps[$i_irrep]}[$i_ops])!=0) { $cmpfl[$i_irrep]=1; $unitfl=0; } 
      elsif ($mychartab{$irreps[$i_irrep]}[$i_ops]!=1) { $unitfl=0; }
   } 
   if($unitfl==1) { $iunit=$i_irrep; }
}

if ($debug) { print join("\t",@cmpfl),"\n\n"; }
if ($debug) { print "Unit representation is ",$irreps[$iunit],"\n"; }

# Determines the reduction of the rotation group in terms of IRREPS of this point group
#    m_irrep = 1/order * sum_(operators p) classmult_p . conj(char_irrep_p) . char_rot_p
for $i_irrep (0..$#irreps) {
   $mval=0;
   if ($cmpfl[$i_irrep]==1) {
      for $idr (1..2) { for $i_ops (0..$#ptops) {
        #$mychar = $mychartab{$irreps[$i_irrep]}[$i_ops]; if($idr==2) { $mychar=~$mychar; }
         $mychar = $mychartab{$irreps[$i_irrep]}[$i_ops]; if($idr==2) { $mychar=cplx(Re($mychar),-Im($mychar)); }
         $mval = $mval + $classmult[$i_ops]*$mychar*($charR[$i_ops]**2);
      } }
      $mval=$mval/2; $mychartab{$irreps[$i_irrep]}[0]=2;
   }
  else {
      for $i_ops (0..$#ptops) {
         $mval = $mval + $mychartab{$irreps[$i_irrep]}[$i_ops]*($charR[$i_ops]**2) * $classmult[$i_ops];
      }
   }
   if(abs(Im($mval))<$SMALL) { $mval = Re($mval); }
   if(abs($mval)<$SMALL) { $mval = 0; }
   push(@m_irrep,$mval/$order);
}

if ($debug) { print join("\t",@m_irrep),"\n"; }

foreach (0..$#m_irrep) { 
   if(abs($m_irrep[$_])>$SMALL) { push(@outstr,$m_irrep[$_],"*",$irreps[$_],"(",$mychartab{$irreps[$_]}[0],") + "); } }
   $outstr = join("",@outstr); $outstr =~ s/\+\ $//; print $outstr,"\n";

# Creates the operators T(G) which transforms the basis vector of Stevens operators under the operation of each 
#   element of the point group.
$j2p = 2*$J+1; $thisop = $rops{$ptgpsym}; $rotop=eye($j2p);
for $ii (0..$#{$thisop}) {
   push (@opchain,$rotop);
   for $jj (1..$#{$thisop->[$ii]}) {
      if($thisop->[$ii][$jj] =~ /p/) { 
#        my $tmp=zeros($j2p,$j2p); 
#        for $ix(0..($j2p-1)) { for $iy(0..($j2p-1)) { $tmp->[$ix][$iy]= ~($opchain[-1]->[$iy][$ix]); } }
#        push(@opchain,$tmp);
      }
      else {
         $D = rotk($J,$thisop->[$ii][$jj][1],$thisop->[$ii][$jj][0]); 
         # Check whether we want Stevens or Buckmaster/Wybourne operators
         if($usestevop) {
            my $tmp = mmult($D,$invA[$J-1]); my $tmp2 = mmult($A[$J-1],$tmp); push(@opchain,mmult($tmp2,$opchain[-1])); }
         else { push(@opchain,mmult($D,$opchain[-1])); }
      } 
   } 
   push(@Tgops,$opchain[-1]); @opchain=();
}

if ($debug) { print "\n";
   for $i_op (0..$#Tgops) {
      print "Matrix representation of operator $rops{$ptgpsym}->[$i_op][0]:\n";
      for $ii (0..($#{$Tgops[$i_op]})) { print "["; for $jj (0..($#{$Tgops[$i_op]})) { 
         $elem = $Tgops[$i_op]->[$ii][$jj];
         if(abs(Re($elem))<$SMALL) { $elem=cplx(0,Im($elem)); }
         if(abs(Im($elem))<$SMALL) { $elem=cplx(Re($elem),0); }
         $elem->display_format('format' => '%5.2f');
         print $elem,"\t"; } print "]\n"; 
      }
      $tracem=0; for $ii (0..(2*$J)) { $tracem+=$Tgops[$i_op]->[$ii][$ii]; }; print "Trace=$tracem\n\n";
   }
}

# Constructs the projection operator from the operators T(G)
$j2s = $j2p**2; $pA1s = zeros($j2s,$j2s); $Aj = zeros($j2s,$j2s); $invAj = zeros($j2s,$j2s);
for $i_ops (0..$#Tgops) {  # Pairs: 4,9 6,11 8,5 10,7
#for $j_ops (0..$#Tgops) {
#for $i_ops (8,5,10,7) { 
   for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) {
      $m1i = $ii % $j2p; $m1j = $jj % $j2p; $m2i = int($ii/$j2p); $m2j = int($jj/$j2p);
      $pA1s->[$ii][$jj] += $Tgops[$i_ops]->[$m1i][$m1j] * $Tgops[$i_ops]->[$m2i][$m2j];
   } 
   if ($i_ops==0) { 
      $v1=$m1i-$J; $v2=$m2i-$J; push(@wvec,sprintf("|%2i,%2i> ",$v1,$v2)); 
      my $namestring = sprintf("O%i%i.O%i%i",$J,$v1,$J,$v2); $namestring =~ s/-([0-7])/\1S/g;
      push(@opname,$namestring);
   } } 
#}
}

# Sets near-zero elements to zero. Criteria is the $SMALL variable declared at the top.
for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) { 
   if(abs(Re($pA1s->[$ii][$jj]))<$SMALL) { $pA1s->[$ii][$jj] = cplx(0,Im($pA1s->[$ii][$jj])); }
   if(abs(Im($pA1s->[$ii][$jj]))<$SMALL) { $pA1s->[$ii][$jj] = cplx(Re($pA1s->[$ii][$jj]),0); }
   $pA1s->[$ii][$jj]/=$order;
} }

# Finds zero columns in the projection operator
for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) { if(abs($pA1s->[$ii][$jj])>$SMALL) { push(@nonzero,$ii); last; } } }

# Prints the projection matrix, or if $J is large, prints it removing all zeros
if ($debug) { 
   $mtrace=0;
   if($J<=1) { print "\nThe projection operator for the unit representation, p_{A1}:\n";
      print "        "; for $ii(0..($j2s-1)) { print $wvec[$ii]; } print "\n";
      for $ii(0..($j2s-1)) { print $wvec[$ii],"["; for $jj(0..($j2s-1)) {
         $el=$pA1s->[$ii][$jj]; $el->display_format('format' => '%5.2g'); print $el,"\t";
      } print "];\n"; $mtrace+=$pA1s->[$ii][$ii]; } print "Trace=$mtrace\n"; 
      print "Nonzero columns are: ",join(",",@nonzero),"\n"; 
   }
   else {
      for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) { if(abs($pA1s->[$jj][$ii])>$SMALL) { push(@nnz,$ii); last; } } }
      %seen = (); @uniqj = grep { ! $seen{$_} ++ } @nnz;      # More stolen code: http://docstore.mik.ua/orelly/perl/cookbook/ch04_07.htm
      %seen = (); @uniqi = grep { ! $seen{$_} ++ } @nonzero;
      print "Truncated projection operator matrix is:\n        "; foreach $ii(@uniqi) { print $wvec[$ii]; } print "\n";
      foreach $ii(@uniqi) { print $wvec[$ii],"["; foreach $jj(@uniqj) { 
         $el=$pA1s->[$ii][$jj]; $el->display_format('format' => '%5.2g'); print $el,"\t"; } print "];\n"; 
         $mtrace+=$pA1s->[$ii][$ii]; } print "Trace=$mtrace\n";
   }
}

#while ( my ($key, $value) = each(%seen) ) { print "$key => $value\n"; }

# Checks that projection operator is correct
#$conjsum=0; for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) { $conjsum+=($pA1s->[$ii][$jj]-~($pA1s->[$jj][$ii])); }}
$conjsum=0; for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) { 
  my $ipr=Im($proj[$_]->[$jj][$ii]); my $vconj=Re($proj[$_]->[$jj][$ii]); 
  if(abs($ipr)>$SMALL) { $vconj=cplx(Re($proj[$_]->[$jj][$ii]),-Im($proj[$_]->[$jj][$ii])); } 
  $conjsum+=abs($proj[$_]->[$ii][$jj]-$vconj); 
} }
if (abs($conjsum)<$SMALL) { print "Projection matrix is hermitian. "; } else { print "Error: projection matrix not hermitian. "; }
$pA1s2 = mmult($pA1s,$pA1s); $p2sum=0; for $ii(0..($j2s-1)) { for $jj(0..($j2s-1)) { $p2sum+=($pA1s2->[$ii][$jj]-$pA1s->[$ii][$jj]); }}
if (abs($p2sum)<$SMALL) { print "Projection matrix is idempotent.\n"; } else { print "Error: projection matrix pA1^2!=pA1.\n"; }

# Determines the maximum value of each of the column vectors of the projection operator in order to normalise them for comparing duplicates
for(0..($j2s-1)) { $maxcol=0; for $jj(0..($j2s-1)) { if(abs($pA1s->[$jj][$_])>$maxcol) { $maxcol=abs($pA1s->[$jj][$_]); } } push(@mxcv,abs($maxcol)); }
$pA1n=zeros($j2s,$j2s); for(0..($j2s-1)) { if($mxcv[$_]>$SMALL) { for $jj(0..($j2s-1)) { $pA1n->[$jj][$_] = $pA1s->[$jj][$_]/$mxcv[$_]; } } }

# Checks whether we have the right number of operators - if not check for duplicates
if (abs($#nonzero+1-$m_irrep[$iunit])>$SMALL) { 
#  @duplicate = @nonzero;
#  print join(",",@duplicate),"\n";
   my %idu; if($J>1) { for(0..$#uniqi) { $idu{$uniqi[$_]}=$_; } }
   @duplicate = ();
   foreach $ic(@nonzero) { 
      $doneflag=0; foreach (@duplicate) { if($_==$ic) { $doneflag=1; last; } } if($doneflag) { next; }
      my @ident; my @identP;
      for $ii(0..($j2s-1)) { if($ii!=$ic) { $plussum=0; $minussum=0; for $jj(0..($j2s-1)) { 
         $plussum+=abs($pA1n->[$jj][$ii]+$pA1n->[$jj][$ic]); $minussum+=abs($pA1n->[$jj][$ii]-$pA1n->[$jj][$ic]);
      } if((abs($plussum)<$SMALL) || (abs($minussum)<$SMALL)) { push(@ident,$ii); } } }
      if ($debug && $#ident>-1)  { 
         if($J>1) { foreach(@ident) { push(@identP,$idu{$_}); } print "Columns ",join(",",@identP)," is same as column $idu{$ic}.\n"; } 
	 else     { if($#ident>-1) { print "Columns ",join(",",@ident)," is same as column $ic.\n"; } }
      }
      foreach (@ident) { push(@duplicate,$_); }
#     foreach $id(@ident) { for $idp (0..$#duplicate) { if($duplicate[$idp]==$id) { delete $duplicate[$idp]; } } push(@idents,\@ident); }
   }
#  print join(",",@duplicate),"\n";
#  print @idents,"\n";
#  print "$#idents\n";
   for $ic(0..$#nonzero) { foreach (@duplicate) { if($_==$nonzero[$ic]) { delete($nonzero[$ic]); } } }
   foreach (@nonzero) { if(defined($_)) { push(@newnz,$_); } }
   my @newnzP; if($J>1) { foreach(@newnz) { push(@newnzP,$idu{$_}); } } else { @newnzP=@newnz; }
   if ($debug) { print "Nonzero, non-duplicate columns are: ",join(",",@newnzP),"\n"; }
   if (abs($#newnz+1-$m_irrep[$iunit])>$SMALL) { 
#     # Try to normalise the projection operator matrix to find duplicate columns
#     my $pA1n = zeros($j2s,$j2s); 
#     for (0..($j2s-1)) { 
#        $normv=0; $maxv=0; 
#        for $jj (0..($j2s-1)) {
#           if($pA1s->[$jj][$_]>$maxv) { $maxv=$pA1s->[$jj][$_]; } 
#           $normv += $pA1s->[$jj][$_]*$pA1s->[$jj][$_]; 
#        }
#        if(abs(Re($normv))<$SMALL){ $normv= cplx(0,Im($normv));} if(abs(Im($normv))<$SMALL){ $normv= cplx(Re($normv),0); }
#        if(abs(Re($maxv))<$SMALL) { $maxv = cplx(0,Im($maxv)); } if(abs(Im($maxv))<$SMALL) { $maxv = cplx(Re($maxv),0); }
#        push(@normcolv,$normv); if(abs($normv)<$SMALL) { $normv=1; } push(@maxnormcolv,$maxv/$normv);
#     } 
##    for (0..($j2s-1)) { print $normcolv[$_]," "; } print "\n"; for (0..($j2s-1)) { print $maxnormcolv[$_]," "; } print "\n";
#     if (abs($#newnz+1-$m_irrep[$iunit])>$SMALL) { 
         print "$ptgpsym Error: number of allowed terms from projection operator does not match character analysis. ";
         print "Expected $m_irrep[$iunit], found ", $#newnz+1,".\n"; 
#     }
   }
}
else { @newnz = @nonzero; }

# Construct the expression of the symmetry allowed Hamiltonian.
$inz=0;
foreach (@newnz) {
#  $normcol=0; for $jj(0..($j2s-1)) { $normcol+=$pA1s->[$jj][$_]*$pA1s->[$jj][$_]; }
#  $maxcol=0; for $jj(0..($j2s-1)) { if($pA1s->[$jj][$_]>$maxcol) { $maxcol=$pA1s->[$jj][$_]; } }
   push(@outhmltn,sprintf("K%s|",chr(97+$inz))); $inz++;
#  for $jj(0..($j2s-1)) { if(abs($pA1s->[$jj][$_])>$SMALL) { push(@outhmltn,$pA1s->[$jj][$_]/abs($maxcol)); push(@outhmltn,$opname[$jj]); } }
   for $jj(0..($j2s-1)) { 
      if(abs($pA1s->[$jj][$_])>$SMALL) { my $val=$pA1s->[$jj][$_]/$mxcv[$_]; $val=round($val*1e4)/1e4; push(@outhmltn,$val); push(@outhmltn,$opname[$jj]); } }
   push(@outhmltn,"\n");
}
$outstr = join(":",@outhmltn); $outstr =~ s/\|:/\(/g; $outstr =~ s/-1/-/g; $outstr =~ s/\:([02-9])/ +\1/g; $outstr =~ s/\:1/+/g; 
#$outstr =~ s/:O/\ O/g;
$outstr =~ s/\(1\:/\(/g; $outstr =~ s/://g; $outstr =~ s/\n/)\ +\ /g; $outstr =~ s/\+\ $//; $outstr =~ s/\+\ \-/-/g; $outstr =~ s/\s+$//;
#$outstr =~ s/(K[a-z])\(([0-9]*)/\2\1\(/g;
print "\nThe symmetry allowed Hamiltonian for $ptgpsym is $outstr\n";

# For dipolar interactions, also print out it in terms of the operators Jx,Jy,Jz
if ($J==1 && !$usebuckop) {
  #$outstr =~ s/O11S\.O11S/-Jy.Jy/g; 
   $outstr =~ s/O11S/Jy/g; $outstr =~ s/O11/Jx/g; $outstr =~ s/O10/Jz/g;
   print "\nThe symmetry allowed Hamiltonian for $ptgpsym is $outstr\n";
}

# For quadrupolar interactions, also print out it in terms of the operators Qxy,Qyz,Qzz,Qzx,Qx2y2
if ($J==2 && !$usebuckop) {
   $outstr =~ s/(O[0-2S]*)\.O22S/2\1.Qxy/g; $outstr =~ s/O22S\./2Qxy./g; $outstr =~ s/22Qxy/4Qxy/g;
   $outstr =~ s/O21S/Qyz/g; $outstr =~ s/O20/Qzz/g; $outstr =~ s/O21/Qzx/g; $outstr =~ s/O22/Qx2y2/g; 
   print "\nThe symmetry allowed Hamiltonian for $ptgpsym is $outstr\n";
}
