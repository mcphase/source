#!/usr/bin/perl
#
# cfsplit.pl - Determines the symmetry allowed splitting for a particular J and point group.
#
# This program is part of the McPhase package, licensed under the GNU GPL v2. Please see the COPYING
# file
#
# Thu Dec 16 15:23:01 WEST 2010 - Duc Le and Martin Rotter

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
   print " $0 - program to determine the symmetry allowed CF splitting from\n";
   print "                point group character tables.\n\n";
   print " Syntax: $0 <PTGP> <J>\n\n";
   print " where <PTGP> is a (case insensitive) string denoting the point group.\n";
   print "              Both Hermman-Maguin and Schoenflies symbols are allowed.\n";
   print "       <J> is the total angular momentum number of the ground multiplet\n";
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
if ($debug) { for (keys %{ $chartab{$ptgpsym} }) { print $_,"\t",join("\t",@{ $chartab{$ptgpsym}{$_} }),"\n"; } }

# Calculates the character of the Rotation group in 2J+1 dimensions and the class multiplicities
for $ialpha (0..$#angs) {
  $alpha = $angs[$ialpha]; $r_op=$ptops[$ialpha];
  if($r_op=~/[Ei]/) {
    push(@charR, 2*$J+1); }
  elsif($alpha==0) {
    push(@charR, 0); }
  else {
    push(@charR, (sin(($J+0.5)*$alpha)/sin($alpha/2))); }
  if (abs($charR[-1])<1e-3) { $charR[-1] = 0; }
  $r_op =~ /^(.)/; $cm=$1; $cm =~ s/[CSsEi]/1/; 
  push(@classmult, $cm); 
}

if ($debug) { print "\t",join("\t",@classmult),"\n"; }
if ($debug) { print "R_J\t",join("\t",@charR),"\n\n"; }

# Checks if the character is complex, we need to include the complex conjugate terms as well...
for $i_irrep (0..$#irreps) { push(@cmpfl,0); for $i_ops (0..$#ptops) { 
   if(Im($mychartab{$irreps[$i_irrep]}[$i_ops])!=0) { $cmpfl[$i_irrep]=1; } } }

if ($debug) { print join("\t",@cmpfl),"\n\n"; }

# Determines the reduction of the rotation group in terms of IRREPS of this point group
#    m_irrep = 1/order * sum_(operators p) classmult_p . conj(char_irrep_p) . char_rot_p
for $i_irrep (0..$#irreps) {
  $mval=0; push(@dims,$mychartab{$irreps[$i_irrep]}[0]);
  if ($cmpfl[$i_irrep]==1) {
    for $idr (1..2) { for $i_ops (0..$#ptops) {
      $mychar = $mychartab{$irreps[$i_irrep]}[$i_ops]; if($idr==2) { $mychar=cplx(Re($mychar),-Im($mychar)); }
      $mval = $mval + $classmult[$i_ops]*$mychar*$charR[$i_ops];
    } }
    $mval/=2; $dims[-1]=2; #$mychartab{$irreps[$i_irrep]}[0]=2;
  }
  else {
    for $i_ops (0..$#ptops) {
      $mval = $mval + $classmult[$i_ops]*$mychartab{$irreps[$i_irrep]}[$i_ops]*$charR[$i_ops];
    }
  }
  if(abs(Im($mval))<1e-3) { $mval = Re($mval); }
  push(@m_irrep,$mval/$order);
}

if ($debug) { print join("\t",@m_irrep),"\n"; }

print "\nAllowed CF levels are: (number paretheses is multiplicity)\n";
foreach (0..$#m_irrep) { 
  if(abs($m_irrep[$_])>1e-3) { push(@outstr,$m_irrep[$_],"*",$irreps[$_],"(",$dims[$_],") + "); } }
  $outstr = join("",@outstr); $outstr =~ s/\+\ $//; print $outstr,"\n";

if((2*$J)%2==1) { exit(0); }

# Calculates which Jz levels below to which representations (multiplets) using the projection operator

# Creates the operators T(G) which transforms the basis vector of Stevens operators under the operation of each 
#   element of the point group.
$j2p = 2*$J+1; $thisop = $rops{$ptgpsym}; $rotop=eye($j2p);
for $ii (0..$#{$thisop}) {
   push (@opchain,$rotop);
   for $jj (1..$#{$thisop->[$ii]}) {
      $D = rotk($J,$thisop->[$ii][$jj][1],$thisop->[$ii][$jj][0]); 
      # Check whether we want Stevens or Buckmaster/Wybourne operators
      if($usestevop) {
         my $tmp = mmult($D,$invA[$J-1]); my $tmp2 = mmult($A[$J-1],$tmp); push(@opchain,mmult($tmp2,$opchain[-1])); }
      else { push(@opchain,mmult($D,$opchain[-1])); }
   } 
   push(@Tgops,$opchain[-1]); @opchain=();
}
for $i_op (0..$#Tgops) { for $ix (0..(2*$J)) { for $iy (0..(2*$J)) { $elem = $Tgops[$i_op]->[$ix][$iy];
   if(abs(Re($elem))<$SMALL) { $elem=cplx(0,Im($elem)); } if(abs(Im($elem))<$SMALL) { $elem=cplx(Re($elem),0); }
   $Tgops[$i_op]->[$ix][$iy] = $elem; 
} } }

#if ($debug) { print "\n"; for $i_op (0..$#Tgops) { print "Matrix representation of operator $rops{$ptgpsym}->[$i_op][0]:\n";
#   for $ii (0..($#{$Tgops[$i_op]})) { print "["; for $jj (0..($#{$Tgops[$i_op]})) { 
#      $elem = $Tgops[$i_op]->[$ii][$jj]; $elem->display_format('format' => '%5.2f'); print $elem,"\t"; } print "]\n"; }
#   $tracem=0; for $ii (0..(2*$J)) { $tracem+=$Tgops[$i_op]->[$ii][$ii]; }; print "Trace=$tracem\n\n";
#} }
#
#@irops = @{$IRops{$ptgpsym}};
#for $i_irrep (0..$#irreps) {
#   my $prop=zeros($j2p,$j2p);
#   if ($cmpfl[$i_irrep]==1) {
#      for $idr (1..2) { for $i_ops (0..$#Tgops) {
#        #$mychar = $mychartab{$irreps[$i_irrep]}[$irops[$i_ops]]; if($idr==2) { $mychar=~($mychar); }
#         $mychar = $mychartab{$irreps[$i_irrep]}[$irops[$i_ops]]; if($idr==2) { $mychar=cplx(Re($mychar),-Im($mychar)); }
#         for (0..($j2p-1)) { for $jj(0..($j2p-1)) { $prop->[$_][$jj] += $Tgops[$i_ops]->[$_][$jj]*($mychar); } }
#      } }
#   }
#   else {
#      for $i_ops (0..$#Tgops) { 
#         $mychar = $mychartab{$irreps[$i_irrep]}[$irops[$i_ops]]; #if(Im($mychar)>$SMALL) { $mychar=~$mychar; }
#         for (0..($j2p-1)) { for $jj(0..($j2p-1)) { $prop->[$_][$jj] += $Tgops[$i_ops]->[$_][$jj]*$mychar; } }
#      }
#   }
#   for $ii(0..($j2p-1)) { for $jj(0..($j2p-1)) { $prop->[$ii][$jj]*=($mychartab{$irreps[$i_irrep]}[0]/$order); } }
#   push(@proj,$prop);
#}
#
#if ($debug) { for(0..$#irreps) { print "Projection operator for IRREP $irreps[$_]:\n"; 
#   for $ii(0..($j2p-1)) { print "["; for $jj(0..($j2p-1)) { $elem = $proj[$_]->[$ii][$jj];
#      if(abs(Re($elem))<$SMALL) { $elem=cplx(0,Im($elem)); } if(abs(Im($elem))<$SMALL) { $elem=cplx(Re($elem),0); }
#      $elem->display_format('format' => '%5.2f'); print $elem,"\t"; } print "]\n"; 
#   }
#   # Checks that projection operator is correct
#   $conjsum=0; 
#   for $ii (0..(2*$J)) { for $jj(0..(2*$J)) { $conjsum+=abs($proj[$_]->[$ii][$jj]-cplx(Re($proj[$_]->[$jj][$ii]),-Im($proj[$_]->[$jj][$ii]))); } }
#   if (abs($conjsum)<$SMALL) { print "Projection matrix is hermitian. "; } else { print "Error: projection matrix not hermitian. "; }
#   $proj2 = mmult($proj[$_],$proj[$_]); $p2sum=0; for $ii(0..(2*$J)) { for $jj(0..(2*$J)) { $p2sum+=($proj2->[$ii][$jj]-$proj[$_]->[$ii][$jj]); } }
#   if (abs($p2sum)<$SMALL) { print "Projection matrix is idempotent.\n"; } else { print "Error: projection matrix pA1^2!=pA1.\n"; }
#   $tracem=0; for $ii (0..(2*$J)) { $tracem+=$proj[$_]->[$ii][$ii]; }; 
#   if (Im($tracem)<$SMALL) { $tracem=Re($tracem); } if (abs($tracem)<$SMALL) { $tracem=0; } print "Trace=$tracem\n\n";
#} }

