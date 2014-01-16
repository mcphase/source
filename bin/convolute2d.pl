#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}


use POSIX qw(ceil floor);

unless ($#ARGV >11) 
{ print STDOUT << "EOF";
# Program: convolute2d
#
# Does a 2D convolution (like with convolute) 
#
# Syntax: convolute2d col_X col_Y col_Z input_file ccol_X ccol_Y ccol_Z resolution_file minx maxx Nx miny maxy Ny
#
# Example: To plot a contour of a dispersion with a gaussian resolution function of
#          fwhm=0.2 in h and fwhm=0.5meV in energy
#
#          gauss2d 0.2 0.5 0 0.05 -0.5 0.5 0.1 -1.0 +1.0 > resolution.dat
#          convolute2d 5 9 10 results/mcdisp.qei 1 2 3 resolution.dat 1 2 11 0 20 40 > results/mcdisp.clc
#          displaycontour 1 2 3 results/mcdisp.clc
EOF
exit(1);
}
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cx =eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;$cy =eval $ARGV[1]; 
$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;$cz =eval $ARGV[2]; 
$file1 = $ARGV[3];
$ARGV[4]=~s/exp/essp/g;$ARGV[4]=~s/x/*/g;$ARGV[4]=~s/essp/exp/g;$c1 =eval $ARGV[4]; 
$ARGV[5]=~s/exp/essp/g;$ARGV[5]=~s/x/*/g;$ARGV[5]=~s/essp/exp/g;$c2 =eval $ARGV[5]; 
$ARGV[6]=~s/exp/essp/g;$ARGV[6]=~s/x/*/g;$ARGV[6]=~s/essp/exp/g;$c3 =eval $ARGV[6]; 
$file2 = $ARGV[7];
$ARGV[8]=~s/exp/essp/g;$ARGV[8]=~s/x/*/g;$ARGV[8]=~s/essp/exp/g;$lx=eval $ARGV[8];
$ARGV[9]=~s/exp/essp/g;$ARGV[9]=~s/x/*/g;$ARGV[9]=~s/essp/exp/g;$ux=eval $ARGV[9];
$ARGV[10]=~s/exp/essp/g;$ARGV[10]=~s/x/*/g;$ARGV[10]=~s/essp/exp/g;$Nx=eval $ARGV[10];
$ARGV[11]=~s/exp/essp/g;$ARGV[11]=~s/x/*/g;$ARGV[11]=~s/essp/exp/g;$ly=eval $ARGV[11];
$ARGV[12]=~s/exp/essp/g;$ARGV[12]=~s/x/*/g;$ARGV[12]=~s/essp/exp/g;$uy=eval $ARGV[12];
$ARGV[13]=~s/exp/essp/g;$ARGV[13]=~s/x/*/g;$ARGV[13]=~s/essp/exp/g;$Ny=eval $ARGV[13];

# determine range of convolution function data
$minx=1e100; $maxx=-1e100; 
$miny=1e100; $maxy=-1e100; 
 $ic=0;
unless (open (Fin1, $file2)){die "\n error:unable to open $file2\n";}
while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);
        if($numbers1[$c1-1]<$minx) {$minx=$numbers1[$c1-1];}
        if($numbers1[$c1-1]>$maxx) {$maxx=$numbers1[$c1-1];}
        if($numbers1[$c2-1]<$miny) {$miny=$numbers1[$c2-1];}
        if($numbers1[$c2-1]>$maxy) {$maxy=$numbers1[$c2-1];}
        #store  function values
        $c1values[$ic]=$numbers1[$c1-1];
        $c2values[$ic]=$numbers1[$c2-1];
        $c3values[$ic]=$numbers1[$c3-1];
        ++$ic;
        }
     }
close Fin1;

# determine range of  function data
$minrx=1e100; $maxrx=-1e100; 
$minry=1e100; $maxry=-1e100; 
 $ii=0;
unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}
while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);
        if($numbers1[$cx-1]+$maxx>$lx&&
           $numbers1[$cx-1]+$minx<$ux&&
           $numbers1[$cy-1]+$maxy>$ly&&
           $numbers1[$cy-1]+$miny<$uy)
        {if($numbers1[$cx-1]<$minrx) {$minrx=$numbers1[$cx-1];}
         if($numbers1[$cx-1]>$maxrx) {$maxrx=$numbers1[$cx-1];}
         if($numbers1[$cy-1]<$minry) {$minry=$numbers1[$cy-1];}
         if($numbers1[$cy-1]>$maxry) {$maxry=$numbers1[$cy-1];}
         #store  function values
         $cxvalues[$ii]=$numbers1[$cx-1];
         $cyvalues[$ii]=$numbers1[$cy-1];
         $czvalues[$ii]=$numbers1[$cz-1];
         ++$ii;
         }
      }
     }
close Fin1;
# initialize output piddle
@a=($Nx,$Ny);
for($ix=0;$ix<$Nx;++$ix){
for($iy=0;$iy<$Ny;++$iy){$a[$ix][$iy]=0;}}

$dx=($ux-$lx)/($Nx-1);
$dy=($uy-$ly)/($Ny-1);

# calculate convolution 
for($i=0;$i<$ii;++$i) # take data points
 { for($j=0;$j<$ic;++$j) # take convolution function points
   {#determine which piddle it should go into
    $x=$cxvalues[$i]+$c1values[$j];
    $ix=($x-$lx)/$dx;
    $ixf=floor($ix);
    $ixc=ceil($ix);
    $y=$cyvalues[$i]+$c2values[$j];
    $iy=($y-$ly)/$dy;
    $iyf=floor($iy);
    $iyc=ceil($iy);
  if($ixf>=0&&$iyf>=0&&$ixc<$Nx&&$iyc<$Ny)
    {$wxc=$ix-$ixf;#if($wxc<0){print "error: $ix $ixf";exit(1);}
     $wxf=1-$wxc; #if($wxc>1){print "error: $ix $ixf";exit(1);}

     $wyc=$iy-$iyf;
     $wyf=1-$wyc;
     
     $z=$czvalues[$i]*$c3values[$j];
    if($ixc>0&&$iyc>0){
     $a[$ixf][$iyf]+=$wxf*$wyf*$z;
     $a[$ixf][$iyc]+=$wxf*$wyc*$z;
     $a[$ixc][$iyf]+=$wxc*$wyf*$z;
     $a[$ixc][$iyc]+=$wxc*$wyc*$z;}
    }
   } #next $j
 } # next $i

# print results
#print "results:\n";
for($ix=0;$ix<$Nx;++$ix){
for($iy=0;$iy<$Ny;++$iy){
$x=$lx+$ix*$dx;
$y=$ly+$iy*$dy;
$z=$a[$ix][$iy];
print $x." ".$y." ".$z."\n";
}
print "#\n";
}