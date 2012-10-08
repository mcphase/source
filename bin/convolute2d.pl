#!/usr/bin/perl
# File: convolute2d.pl
# Does a 1D convolution (like with convolute) for lines in y for a set of x,y,z data points.
#
# Syntax: perl convolute2d.pl col_X col_Y col_Z input_file ccol_Y ccol_Z resolution_file
#
# Example: To plot a contour of a dispersion with a (constant) resolution function:
#
#          perl convolute2d.pl 5 9 10 results/mcdisp.qei 1 2 resolution.dat > results/mcdisp.clc
#          displaycontour 1 2 3 results/mcdisp.clc

use POSIX qw(ceil floor);

$c1 = $ARGV[0]; $c2 = $ARGV[1]; $c3 = $ARGV[2]; $file1 = $ARGV[3];
$cx = $ARGV[4]; $cy = $ARGV[5];                 $file2 = $ARGV[6];

# determine range of convolution function data
$minr=1e100; $maxr=-1e100; $delta=1e100; $oldnr=-1e100; $ii=0;$i=1;
unless (open (Fin1, $file2)){die "\n error:unable to open $file2\n";}
while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);
        if($numbers1[$cx-1]<$minr) {$minr=$numbers1[$cx-1];}
        if($numbers1[$cx-1]>$maxr) {$maxr=$numbers1[$cx-1];}
        if(abs($numbers1[$cx-1]-$oldnr)<$delta){$delta=abs($numbers1[$cx-1]-$oldnr);}
        $oldnr=$numbers1[$cx-1];
        #store convolution function values
        ++$ii;
        $cxvalues[$ii]=$numbers1[$cx-1];
        $cyvalues[$ii]=$numbers1[$cy-1];
        }
     }

$Xtest = -100000; $first = 1;
$ii1=0; $min=1e100; $max=-1e100; $c1--; $c2--; $c3--;

unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}
while(<Fin1>) {
  if ($_!~/^\s*\#/) {              # Ignores comment lines
    chomp; $_ =~ s/^\s*//;         # Removes trailing newline and leading spaces
    @Line = split(/\s+/);          # Splits line into columns
    if($Line[$c1]!=$Xtest) {       # Check if x has changed, and if so
      $Xtest = $Line[$c1];         #    updates index of where the blocks of x are.
      push @Xdat, $Line[$c1];
      push @Ybegin, $ii1;
      if($first!=1) { 
        push @Yend, $ii1;
      }
      else { $first = 0; }
    }
    $ii1++;
    if($Line[$c2]>$max) { $max = $Line[$c2]; }  # Determines maximum value of y
    if($Line[$c2]<$min) { $min = $Line[$c2]; }  # Determines minimum value of y
    push @Ydat, $Line[$c2];        # Appends to Y data
    push @Zdat, $Line[$c3];        # Appends to Z data
  }
}
push @Ybegin, $ii1;
close(Fin1);
$max = ceil($max); $min = floor($min);

# determine range and step of x column for output [min,max] with step delta
$min+=$minr; $max+=$maxr;

for($iX=0; $iX<=$#Xdat; $iX++) {
   $i=1;
   for($y=$min; $y<=$max; $y+=$delta) {
      $z=0;
      for($ip=$Ybegin[$iX]; $ip<$Yend[$iX]; $ip++) {
        #if($Idat[$ip]>1e-7) {
        #  $z += $Idat[$ip] * exp(-4*log(2)*(($Et-$Edat[$ip])/$fwhm)**2);
        #}
         $dd=$y-$Ydat[$ip];
         if ( ($dd>$minr) && ($dd<$maxr) ) {
            $imax=$ii; $imin=1;
            $zi=$Zdat[$ip];
            while($imax-$imin>1)  {   # intervallschachtelung
               if($dd<$cxvalues[$i]) {
                  $imax=$i; $i=int($imin+($i-$imin)/2); }
               else {
                  $imin=$i; $i=int($i+($imax-$i)/2); }
            }
            $z += $zi*($cyvalues[$imin]+($dd-$cxvalues[$imin])*($cyvalues[$imin+1]-$cyvalues[$imin])/($cxvalues[$imin+1]-$cxvalues[$imin]));
         }
      }
      print sprintf("%10.9e %10.9e %10.9e\n",$Xdat[$iX],$y,$z);
   }
}

