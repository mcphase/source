#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}


#use POSIX qw(ceil floor);

# Rounds up a value, using int() function
sub my_ceil
{
    my $t = int($_[0]);
    return ($t != $_[0]) ? ($_[0] < 0) ? $t : $t + 1 : $_[0];
}

# Rounds down a value, using int() function
sub my_floor
{   my $t = int($_[0]);
    return ($t != $_[0]) ? ($_[0] < 0) ? $t - 1 : $t : $_[0]-1;
}

# Rounds a value, using int() function
sub my_round
{
    return ($_[0] < 0) ? int($_[0] - 0.5) : int($_[0] + 0.5);
}


# test
#for $i ( qw/3 5.5 -3.5 0.0005 -1.4 4.00005 -4.0005 -0.0005/) {
#     local $,="\t"; local $\="\n";
#      print "ceil:  ", $i, my_ceil($i);
#     print "floor: ", $i, my_floor($i);
#}


 unless ($^O=~/MSWin/){$ds="\"";}

unless ($#ARGV >11) 
{ print STDERR << "EOF";
# Program: convolute2d
#
# Does a 2D convolution (like with convolute) 
#
# Syntax: 1) convolute2d col_X col_Y col_Z input_file ccol_X ccol_Y ccol_Z resolution_file minx maxx Nx miny maxy Ny
#     or  2) convolute2d col_X col_Y col_Z input_file ccol_X ccol_Y ccol_Z resolution_file -d[f] dx dy dz data_file [sx sy]
#
# formula:  f(x,y) = sum_i,j zi cj(x-xi,y-yi)
#
# input file: xi yi zi  
# resolution_file: rxj ryj cj(rxi,ryi)
#
# input_file and resolution_file may be unsorted, equal steps not necessary
# data_file has to be a Nx times Ny points grid with equal x and y spacing
#
# output to stdout 3 columns with
#       column  1(2) x(y) range from minx(y) to maxx(y) with Nx(y) steps and column 3 result of the convolution
#       if syntx 2 is used,  an additional column 4 is added containing the scaled dz data and
#                 standard deviation, volume of curves etc are output to
#                 stdout and to environment variables MCPHASE_STA etc.
#
# Options: -df      instead of -d will not output convolution but only standard deviations
#          sx sy    columns sx,sy of the datafile may contain 
#                         stretching factors in x,y direction for the convolution function in 
#                         order to allow for x,y dependent resolution		
#
# Note on normalisation assuming that the input functions are normalised and equal spaced: i.e.
#                         sum_i zi(xi,yi) dx dy =1  
#                         sum_j cj(rxi,ryi) drx dry =1  
#     
#  stepwidth of function and resolution function are dx,dy and drx,dry respectively
#  stepwidth of output dfx=(maxx-minx)/Nx,dfy=(maxy-miny)/Ny
#
#                        multiply f(x,y) by dx*dy*drx*dry/(dfx*dfy) 
#
#     to get a normalised function obeying sum_k f(xk,yk) dfx dfy = 1
# 
#
# Example: To plot a contour of a dispersion with a gaussian resolution function of
#          fwhm=0.2 in h and fwhm=0.5meV in energy
#
#          gauss2d 0.2 0.5 0 0.05 -0.5 0.5 0.1 -1.0 +1.0 > resolution.dat
#          convolute2d 5 9 10 results/mcdisp.qei 1 2 3 resolution.dat 1 2 11 0 20 40 > results/mcdisp.clc
#          displaycontour 1 2 3 results/mcdisp.clc
EOF
exit(1);
}else{print STDERR "#* $0 *\n";}
if ($^O=~/MSWin/){print "\@echo off\n";}
print "echo $ds# $0 @ARGV$ds\n";
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cx =eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;$cy =eval $ARGV[1]; 
$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;$cz =eval $ARGV[2]; 
$file1 = $ARGV[3];
$ARGV[4]=~s/exp/essp/g;$ARGV[4]=~s/x/*/g;$ARGV[4]=~s/essp/exp/g;$c1 =eval $ARGV[4]; 
$ARGV[5]=~s/exp/essp/g;$ARGV[5]=~s/x/*/g;$ARGV[5]=~s/essp/exp/g;$c2 =eval $ARGV[5]; 
$ARGV[6]=~s/exp/essp/g;$ARGV[6]=~s/x/*/g;$ARGV[6]=~s/essp/exp/g;$c3 =eval $ARGV[6]; 
$file2 = $ARGV[7];
if ($ARGV[8]=~/-d/) # usage 2 with datafile
{
$ARGV[9]=~s/exp/essp/g;$ARGV[9]=~s/x/*/g;$ARGV[9]=~s/essp/exp/g;$dx =eval $ARGV[9]; 
$ARGV[10]=~s/exp/essp/g;$ARGV[10]=~s/x/*/g;$ARGV[10]=~s/essp/exp/g;$dy =eval $ARGV[10]; 
$ARGV[11]=~s/exp/essp/g;$ARGV[11]=~s/x/*/g;$ARGV[11]=~s/essp/exp/g;$dz =eval $ARGV[11]; 
$file3 = $ARGV[12];
if($#ARGV>13)
 {
$ARGV[13]=~s/exp/essp/g;$ARGV[13]=~s/x/*/g;$ARGV[13]=~s/essp/exp/g;$sx =eval $ARGV[13]; 
$ARGV[14]=~s/exp/essp/g;$ARGV[14]=~s/x/*/g;$ARGV[14]=~s/essp/exp/g;$sy =eval $ARGV[14]; 
 }
}
else
{ # usage 1
$ARGV[8]=~s/exp/essp/g;$ARGV[8]=~s/x/*/g;$ARGV[8]=~s/essp/exp/g;$lx=eval $ARGV[8];
$ARGV[9]=~s/exp/essp/g;$ARGV[9]=~s/x/*/g;$ARGV[9]=~s/essp/exp/g;$ux=eval $ARGV[9];
$ARGV[10]=~s/exp/essp/g;$ARGV[10]=~s/x/*/g;$ARGV[10]=~s/essp/exp/g;$Nx=eval $ARGV[10];
$ARGV[11]=~s/exp/essp/g;$ARGV[11]=~s/x/*/g;$ARGV[11]=~s/essp/exp/g;$ly=eval $ARGV[11];
$ARGV[12]=~s/exp/essp/g;$ARGV[12]=~s/x/*/g;$ARGV[12]=~s/essp/exp/g;$uy=eval $ARGV[12];
$ARGV[13]=~s/exp/essp/g;$ARGV[13]=~s/x/*/g;$ARGV[13]=~s/essp/exp/g;$Ny=eval $ARGV[13];
}

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

# determine range of  data_file
$mindx=1e100; $maxdx=-1e100; 
$mindy=1e100; $maxdy=-1e100; 
$di=0;$deltax=1e100;$deltay=1e100;
if ($ARGV[8]=~/-d/) # usage 2 with datafile
{unless (open (Fin1, $file3)){die "\n error:unable to open $file3\n";}
while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);
         if($numbers1[$dx-1]<$mindx) {$mindx=$numbers1[$dx-1];}
         if($numbers1[$dx-1]>$maxdx) {$maxdx=$numbers1[$dx-1];}
         if($numbers1[$dy-1]<$mindy) {$mindy=$numbers1[$dy-1];}
         if($numbers1[$dy-1]>$maxdy) {$maxdy=$numbers1[$dy-1];}
         #store  function values
         $dxvalues[$di]=$numbers1[$dx-1];
         $dyvalues[$di]=$numbers1[$dy-1];
         $dzvalues[$di]=$numbers1[$dz-1];
if($#ARGV>13){$sxvalues[$di]=$numbers1[$sx-1];
              $syvalues[$di]=$numbers1[$sy-1];
             }
        if($di>0){$try=abs($dxvalues[$di]-$dxvalues[$di-1]);if($try>0&&$try<$deltax){$deltax=$try;} # try to find $deltax and $deltay
                  $try=abs($dyvalues[$di]-$dyvalues[$di-1]);if($try>0&&$try<$deltay){$deltay=$try;}
                 }
         ++$di;
         
      }
     }
close Fin1;
$lx=$mindx;$ux=$maxdx;$Nx=1+my_round(($ux-$lx)/$deltax);
$ly=$mindy;$uy=$maxdy;$Ny=1+my_round(($uy-$ly)/$deltay);

}

# determine range of  function in input_file 
$minrx=1e100; $maxrx=-1e100; 
$minry=1e100; $maxry=-1e100; 
 $ii=0;
unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}
while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);
        if($numbers1[$cx-1]+$maxx>=$lx&&
           $numbers1[$cx-1]+$minx<=$ux&&
           $numbers1[$cy-1]+$maxy>=$ly&&
           $numbers1[$cy-1]+$miny<=$uy)
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
$deltax=($ux-$lx)/($Nx-1);
$deltay=($uy-$ly)/($Ny-1);
for($ix=0;$ix<$Nx;++$ix){
for($iy=0;$iy<$Ny;++$iy){$a[$ix][$iy]=0;}}

#-----------------------------------------------
if ($ARGV[8]=~/-d/) # usage 2 with datafile - store data and stretching factors in maps data stretchx stretchy
{for($ix=0;$ix<$Nx;++$ix){
for($iy=0;$iy<$Ny;++$iy){$data[$ix][$iy]=0;}}

if($#ARGV>13){
    for($ix=0;$ix<$Nx;++$ix){
    for($iy=0;$iy<$Ny;++$iy){$stretchx[$ix][$iy]=0;$stretchy[$ix][$iy]=0;}}
              }

for($j=0;$j<$di;++$j) # take data_file points
   {#determine which piddle it should go into
if($a[0][0]!=0){print STDERR "error j=$j a00=".$a[0][0]."\n";exit(1);}

    $x=$dxvalues[$j];
    $ix=($x-$lx)/$deltax;
    $ixf=my_floor($ix);
    $ixc=my_ceil($ix);
    $y=$dyvalues[$j];
    $iy=($y-$ly)/$deltay;
    $iyf=my_floor($iy);
    $iyc=my_ceil($iy);
  if($ixf>=-1&&$iyf>=-1&&$ixc<=$Nx&&$iyc<=$Ny)
    {$wxc=$ix-$ixf;

     $wyc=$iy-$iyf; 
     $wyf=1-$wyc;
     #if($wxc<0){print "error wxc<0: $ix $ixf";exit(1);}
     #if($wxc>1){print "error wxc>1: $ix $ixf";exit(1);}
     #if($wyc<0){print "error wyc<0: $iy $iyf";exit(1);}
     #if($wyc>1){print "error wyc>1: $iy $iyf";exit(1);}
#if($ixc<0){print STDERR "error ixf=$ixf ixc=$ixc ";exit(1);}
#if($iyc<0){print STDERR "error iyf=$iyf iyc=$iyc ";exit(1);}
#if($ixf>$Nx){print STDERR "error ixf=$ixf ixc=$ixc Nx=$Nx";exit(1);}
#if($iyf>$Ny){print STDERR "error iyf=$iyf iyc=$iyc Ny=$Ny";exit(1);}

     $z=$dzvalues[$j];$zx=$sxvalues[$j];$zy=$syvalues[$j];
    # distribute $z to 4 nearest grid points according to distance
    if($ixf>=0){$wxf=1-$wxc;
                if($iyf>=0){
     $data[$ixf][$iyf]+=$wxf*$wyf*$z;if($#ARGV>13){$stretchx[$ixf][$iyf]+=$wxf*$wyf*$zx;$stretchy[$ixf][$iyf]+=$wxf*$wyf*$zy;}
}
               if($iyc<$Ny){
     $data[$ixf][$iyc]+=$wxf*$wyc*$z;if($#ARGV>13){$stretchx[$ixf][$iyc]+=$wxf*$wyc*$zx;$stretchy[$ixf][$iyc]+=$wxf*$wyc*$zy;}
}
                }
    if($ixc<$Nx){if($iyf>=0){
     $data[$ixc][$iyf]+=$wxc*$wyf*$z;if($#ARGV>13){$stretchx[$ixc][$iyf]+=$wxc*$wyf*$zx;$stretchy[$ixc][$iyf]+=$wxc*$wyf*$zy;}
}
                if($iyc<$Ny){
     $data[$ixc][$iyc]+=$wxc*$wyc*$z;if($#ARGV>13){$stretchx[$ixc][$iyc]+=$wxc*$wyc*$zx;$stretchy[$ixc][$iyc]+=$wxc*$wyc*$zy;}
}
                }
    }
   } #next $j

}

#-----------------------------------------------------------------
$strx=1;$stry=1;
# calculate convolution 
for($i=0;$i<$ii;++$i) # take data points
 {    if ($ARGV[8]=~/-d/&&$#ARGV>13){ # take into account stretching factor ... determine it for the data point $i
                                     $x=$cxvalues[$i];
                                     $ix=($x-$lx)/$deltax;
                                     $ixf=my_floor($ix);
                                     $ixc=my_ceil($ix);
                                     $y=$cyvalues[$i];
                                     $iy=($y-$ly)/$deltay;
                                     $iyf=my_floor($iy);
                                     $iyc=my_ceil($iy);
                                     # if($ixf>-1&&$iyf>-1&&$ixc<$Nx&&$iyc<$Ny)
                                      $wxc=$ix-$ixf;$wxf=1-$wxc;if($ixf<0){$ixf=0;$ixc=0;$wxc=1;$wxf=0;}
                                                                if($ixc>$Nx-1){$ixf=$Nx-1;$ixc=$Nx-1;$wxc=1;$wxf=0;}
                                      $wyc=$iy-$iyf;$wyf=1-$wyc;if($ixyf<0){$iyf=0;$iyc=0;$wyc=1;$wyf=0;}
                                                                if($iyc>$Ny-1){$iyf=$Ny-1;$iyc=$Ny-1;$wyc=1;$wyf=0;}
                                      $strx=$wxf*$wyf*$stretchx[$ixf][$iyf]+$wxc*$wyc*$stretchx[$ixc][$iyc]+$wxc*$wyf*$stretchx[$ixc][$iyf]+$wxf*$wyc*$stretchx[$ixf][$iyc];
                                      $stry=$wxf*$wyf*$stretchy[$ixf][$iyf]+$wxc*$wyc*$stretchy[$ixc][$iyc]+$wxc*$wyf*$stretchy[$ixc][$iyf]+$wxf*$wyc*$stretchy[$ixf][$iyc];
                                    }

  for($j=0;$j<$ic;++$j) # take convolution function points
   {#determine which piddle it should go into
    $x=$cxvalues[$i]+$c1values[$j]*$strx;
    $ix=($x-$lx)/$deltax;
    $ixf=my_floor($ix);
    $ixc=my_ceil($ix);
    $y=$cyvalues[$i]+$c2values[$j]*$stry;
    $iy=($y-$ly)/$deltay;
    $iyf=my_floor($iy);
    $iyc=my_ceil($iy);
  if($ixf>=-1&&$iyf>=-1&&$ixc<=$Nx&&$iyc<=$Ny)
    {$wxc=$ix-$ixf;$wxf=1-$wxc;
     $wyc=$iy-$iyf;$wyf=1-$wyc;
     #if($wxc<0){print "error wxc<0: $ix $ixf";exit(1);}
     #if($wxc>1){print "error wxc>1: $ix $ixf";exit(1);}
     #if($wyc<0){print "error wyc<0: $iy $iyf";exit(1);}
     #if($wyc>1){print "error wyc>1: $iy $iyf";exit(1);}
# if($ixf<2&&$iyf==4){print STDERR "error : ix=$ix ixf=$ixf ixc=$ixc z=$z wxf=$wxf wxc=$wxc\n";}
     $z=$czvalues[$i]*$c3values[$j]/$strx/$stry;
       # distribute $z to 4 nearest grid points according to distance
    if($ixf>=0){
                if($iyf>=0){ $a[$ixf][$iyf]+=$wxf*$wyf*$z;}
               if($iyc<$Ny){ $a[$ixf][$iyc]+=$wxf*$wyc*$z;}
                }            
    if($ixc<$Nx){if($iyf>=0){ $a[$ixc][$iyf]+=$wxc*$wyf*$z;}
                if($iyc<$Ny){ $a[$ixc][$iyc]+=$wxc*$wyc*$z;}
                }
    }
   } #next $j
 } # next $i

# print results
$sta=0;$volcalc=0;$voldata=0;
for($ix=0;$ix<$Nx;++$ix){
for($iy=0;$iy<$Ny;++$iy){
$x=$lx+$ix*$deltax;
$y=$ly+$iy*$deltay;
$z=$a[$ix][$iy];
$zd=$data[$ix][$iy];
 $sta+=($z-$zd)*($z-$zd);
$volcalc+=$dx*$dy*$z;
$voldata+=$dx*$dy*$zd;
if ($ARGV[8]=~/-d/){;}else {print "echo $ds".sprintf("%+10.9e %10.9e %10.9e$ds\n",$x,$y,$z);}
}
# print "echo $ds#\n";
}

if ($ARGV[8]=~/-d/)
{

$stanorm=0;$stacalc=0;
if ($voldata==0) {print STDERR "Error reading data points or volume below data points zero\n";}
  $scale=$volcalc/$voldata;
for($ix=0;$ix<$Nx;++$ix){
for($iy=0;$iy<$Ny;++$iy){
$x=$lx+$ix*$deltax;
$y=$ly+$iy*$deltay;
$z=$a[$ix][$iy];
$zd=$data[$ix][$iy];
 $stanorm+=($z-$zd*$scale)*($z-$zd*$scale);
if ($ARGV[8]=~/-df/){;}else {print "echo $ds".sprintf("%+10.9e %10.9e %10.9e %10.9e$ds\n",$x,$y,$z,$zd*$scale);}


}}
print STDOUT << "EOF";
echo $ds# ***************************************************************$ds
echo $ds# result of: convolute2d $c1 $c2 $c3 $file1 $cx $cy $cz $file2 $dx $dy $dz $file3 $sx $sy $ds
echo $ds#convolution of data in $file1 (x column is $c1 y column is $c2 z column is $c3)$ds
echo $ds#with resolution function from $file2 (x column is $cx y column is $cy z column is $cz)$ds
echo $ds#evaluated at data points from $file3 (x column $dx y column $dy $z column $dz)$ds
EOF
if ($#ARGV>13)
{
print STDOUT << "EOF";
echo $ds#using stretching factor for resolution function in columns $sx $sy in $file3$ds
EOF
}
print STDOUT << "EOF";
echo $ds#                 the above output contains data from $file3 as given, however$ds
echo $ds#                 with a scaled column $d2. Two additional columns are added $ds
echo $ds#                 containing the calculated results of the convolution and the$ds 
echo $ds#                 original unscaled data. The following environment variables are set:$ds
EOF
 print sprintf("echo $ds#$ds\necho $ds#!sta=%+10.9e MCPHASE_STA ... sum of squared deviations (data - convolution result)^2$ds\n",$sta);
 print  sprintf("echo $ds#!voldata=%+10.9e MCPHASE_VOLDATA$ds\n",$voldata);
 print sprintf( "echo $ds#!volcalc=%+10.9e MCPHASE_VOLCALC$ds\n",$volcalc);
 print  sprintf("echo $ds#!column %i scaled by$ds\necho $ds#!scale_factor=%+10.9e MCPHASE_SCALEFACTOR ... volcalc/voldata$ds\n",$dz,$scale);
 print  sprintf("echo $ds#!sta_of_normalized_curves=%+10.9e MCPHASE_STA_OF_NORMALIZED_CURVES ... sum of squared deviations (data*scale_factor-convolution result)^2$ds\n",$stanorm);
if ($scale==0){$stacalc=1e100;}else{$stacalc=$stanorm/$scale/$scale;}
 print  sprintf("echo $ds#!sta_of_normalized_calc=%+10.9e MCPHASE_STA_OF_NORMALIZED_CALC ... sum of squared deviations (data-convolution result/scale_factor)^2$ds\n",$stacalc);

if ($^O=~/MSWin/){
print "set MCPHASE_STA=$sta\n";
print "set MCPHASE_VOLDATA=$voldata\n";
print "set MCPHASE_VOLCALC=$volcalc\n";
print "set MCPHASE_SCALEFACTOR=$scale\n";
print "set MCPHASE_STA_OF_NORMALIZED_CURVES=$stanorm\n";
print "set MCPHASE_STA_OF_NORMALIZED_CALC=$stacalc\n";
                  }
                 else
                  {
print "export MCPHASE_STA=$sta\n";
print "export MCPHASE_VOLDATA=$voldata\n";
print "export MCPHASE_VOLCALC=$volcalc\n";
print "export MCPHASE_SCALEFACTOR=$scale\n";
print "export MCPHASE_STA_OF_NORMALIZED_CURVES=$stanorm\n";
print "export MCPHASE_STA_OF_NORMALIZED_CALC=$stacalc\n";
                  }


}



print STDOUT << "EOF";
echo $ds#$ds
echo $ds#                     McPhase Software$ds
echo $ds#$ds
echo $ds#                     please reference$ds
echo $ds#$ds
echo $ds# M. Rotter and A. Boothroyd Phys. Rev. B 79 (2009) 140405R$ds
echo $ds#$ds
echo $ds# ***************************************************************$ds
EOF