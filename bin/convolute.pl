#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}

 unless ($^O=~/MSWin/){$ds="\"";}


unless ($#ARGV >4) 
{print STDERR " program convolute:

            Usage: convolute c1 c2 file cx cy convfuncfile [d1 d2 datafile [d3]]
	        
	           convolutes data given as column c1 vs column c2 in file 
                 (data pairs xi,yi) with the convolution function given in 
		     columnn cx vs cy  of convfuncfile (function c(x))
		     
                 Range and stepwidth of output is determined from range 
                 and step of convfuncfile unless a datafile is given. 
                 Values out of range of convfile are assumed to be zero,
		     convfile has to be sorted according to ascending x.

                 If a  datafile is given, with data column d1 and d2,the result
                 of the convolution is calculated for x-values of data column 
                 d1 and f(x) is compared to data in column d2 - a standard 
                 deviation sta is calculated as a sum of squared deviations.
                 As output the datafile is given, however with a scaled column
                 d2 and 2 additional columns are added containing the calculated
                 results of the convolution and the original unscaled data.
                 standard deviation, area of curves etc are output to
                 stdout and to environment variables MCPHASE_STA etc.
		     
		     Formula: f(x)=sum_i yi c(x-xi) , 
		     
                 output is written to stdout.
		     
                 Option: a column d3 of the data file may contain 
                         a stretching factor for the convolution function in 
                         order to allow for x dependent resolution \n";
 exit 0;}else{print STDERR "#* $0 *\n";}


print "echo $ds# $0 @ARGV$ds\n";
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$c1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$c2=eval $ARGV[0];shift @ARGV;
$file1=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cy=eval $ARGV[0];shift @ARGV;
$file2=$ARGV[0];shift @ARGV;

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$d1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$d2=eval $ARGV[0];shift @ARGV;
$file3=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$d3=eval $ARGV[0];shift @ARGV;

# determine range of convolution function data
$minr=1e100;
$maxr=-1e100;
$delta=1e100;
$oldnr=-1e100;   
$ii=0;$i=1;
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

$ii1=0;
$min=1e100;
$max=-1e100;
      unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}   

while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);

	if($numbers1[$c1-1]<$min) {$min=$numbers1[$c1-1];}
	if($numbers1[$c1-1]>$max) {$max=$numbers1[$c1-1];}
	#store  values
	$xvalues[$ii1]=$numbers1[$c1-1];
	$yvalues[$ii1]=$numbers1[$c2-1];	
	++$ii1;
	}
     }
   close Fin1;     

unless ($file3)
{
# determine range and step of x column for output [min,max] with step delta
$min+=$minr;
$max+=$maxr;

# convolute each x point separately and output   
for($x=$min;$x<$max;$x+=$delta)   
   { $y=0;
# convolute each x point separately and output   
     $y=0;
#     open (Fin1, $file1);
   $i1=0;
   foreach(@xvalues)
     {
#       if ($line1=~/^\s*#/) {;}
#       else{# take one point of data file one
#             $line1=~s/D/E/g;@numbers1=split(" ",$line1);
#	     $xi=$numbers1[$c1-1];
#	     $yi=$numbers1[$c2-1];
             # calculate c(x-xi)
	     $dd=$x-$xvalues[$i1];

	     if (($dd<$maxr)&&($dd>$minr))
	      {$imax=$ii;
	       $imin=1;
               $yi=$yvalues[$i1];
	       # intervallschachtelung
	       while($imax-$imin>1)
	       {
	       if($dd<$cxvalues[$i])
	         {$imax=$i;$i=int($imin+($i-$imin)/2);}else
		 {$imin=$i;$i=int($i+($imax-$i)/2);}
	       }
               $y+=$yi*($cyvalues[$imin]+($dd-$cxvalues[$imin])*($cyvalues[$imin+1]-$cyvalues[$imin])/($cxvalues[$imin+1]-$cxvalues[$imin])); 
	       
	      } 
             ++$i1;
   
     }
#     close Fin1;
        	     
    print "echo $ds".sprintf("%+10.9e %10.9e$ds\n",$x,$y);
   }
}
else
{#here do something if the "data"-file3 is given ....
 $sta=0;$areadata=0;$areacalc=0;$iii=0;
    unless (open (Fin3, $file3)){die "\n error:unable to open $file3\n";}   
 while($line3=<Fin3>)
 {
  unless ($line3=~/^\s*#/)
  {# get x-values from file 3 
             $line3=~s/D/E/g;@numbers3=split(" ",$line3);
	     $x=$numbers3[$d1-1];
	     $yd=$numbers3[$d2-1];
             unless($d3){$stretch=1;}else{$stretch=$numbers3[$d3-1];}
	
# convolute each x point separately and output   
     $y=0;
#     open (Fin1, $file1);
   $i1=0;
   foreach(@xvalues)
     {
#       if ($line1=~/^\s*#/) {;}
#       else{# take one point of data file one
#             $line1=~s/D/E/g;@numbers1=split(" ",$line1);
#	     $xi=$numbers1[$c1-1];
#	     $yi=$numbers1[$c2-1];}
             # calculate c(x-xi)
	     $dd=$x-$xvalues[$i1];
	     if (($dd<$maxr*$stretch)&&($dd>$minr*$stretch))
	      {
               $yi=$yvalues[$i1];
	       $imax=$ii;
	       $imin=1;
	       # intervallschachtelung
	       while($imax-$imin>1)
  	        {
	        if($dd<$cxvalues[$i]*$stretch)
	         {$imax=$i;$i=int($imin+($i-$imin)/2);}else
		 {$imin=$i;$i=int($i+($imax-$i)/2);}
  	        }
#	       print $imin." ".$i." ".$imax." ".$ii."\n";
               $y+=$yi*($cyvalues[$imin]/$stretch+($dd-$cxvalues[$imin]*$stretch)*($cyvalues[$imin+1]/$stretch-$cyvalues[$imin]/$stretch)/($cxvalues[$imin+1]*$stretch-$cxvalues[$imin]*$stretch)); 
       }
             ++$i1;
     }
#     close Fin1;

        if($iii>0)
	  {$areadata+=($yd+$ydold)*($x-$xold)/2;
	   $areacalc+=($y+$yold)*($x-$xold)/2;
	  }
	$xold=$x;
	$ydold=$yd;
	$yold=$y;
        ++$iii;
        $ydata[$iii]=$yd;
	$ycalc[$iii]=$y;
		     
    $sta+=($y-$yd)*($y-$yd);
   }
  } 
 close Fin3;
  $stanorm=0;
  if ($areadata==0) {print STDERR "Error reading data points or area below data points zero\n";}
  $scale=$areacalc/$areadata;
 $ii=0;
 open (Fin3, $file3);
 while($line3=<Fin3>)
 {
  if ($line3=~/^\s*#/) {$line3=~s/\n//;print "echo $ds".$line3."$ds\n";}
  else
  {# get x-values from file 3 
             $line3=~s/D/E/g;@numbers3=split(" ",$line3);
	     $x=$numbers3[$d1-1];
	     $yd=$numbers3[$d2-1];
#             $if($d3){$stretch=$numbers3[$d3-1];}else{$stretch=1;}
   ++$ii;
  $stanorm+=($ydata[$ii]*$scale-$ycalc[$ii])*($ydata[$ii]*$scale-$ycalc[$ii]); 
  $yorig=$numbers3[$d2-1];
   $numbers3[$d2-1]*=$scale;
            	  $i=0;print "echo $ds";
		  foreach (@numbers3)
		  {++$i;print $numbers3[$i-1]." ";}     
    print sprintf(" %+10.9e %10.9e$ds\n",$ycalc[$ii],$yorig);
  }
 } 
print STDOUT << "EOF";
echo $ds# ***************************************************************$ds
echo $ds# result of: convolute $c1 $c2 $file1 $cx $cy $file2 $d1 $d2 $file3 $d3$ds
echo $ds#convolution of data in $file1 (x column is $c1 y column is $c2)$ds
echo $ds#with resolution function from $file2 (x column is $cx y column is $cy)$ds
echo $ds#evaluated at data point from $file3 (x column $d1 y column $d2)$ds
EOF
if ($d3)
{
print STDOUT << "EOF";
echo $ds#using stretching factor for resolution function in column $d3 in $file3$ds
EOF
}
print STDOUT << "EOF";
echo $ds#                 the above output contains data from $file3 as given, however$ds
echo $ds#                 with a scaled column $d2. Two additional columns are added $ds
echo $ds#                 containing the calculated results of the convolution and the$ds 
echo $ds#                 original unscaled data. The following environment variables are set:$ds
EOF
 print sprintf("echo $ds#$ds\necho $ds#!sta=%+10.9e MCPHASE_STA ... sum of squared deviations (data - convolution result)^2$ds\n",$sta);
 print  sprintf("echo $ds#!areadata=%+10.9e MCPHASE_AREADATA$ds\n",$areadata);
 print sprintf( "echo $ds#!areacalc=%+10.9e MCPHASE_AREACALC$ds\n",$areacalc);
 print  sprintf("echo $ds#!column %i scaled by$ds\necho $ds#!scale_factor=%+10.9e MCPHASE_SCALEFACTOR ... areacalc/areadata$ds\n",$d2,$scale);
 print  sprintf("echo $ds#!sta_of_normalized_curves=%+10.9e MCPHASE_STA_OF_NORMALIZED_CURVES ... sum of squared deviations (data*scale_factor-convolution result)^2$ds\n",$stanorm);
if ($scale==0){$stacalc=1e100;}else{$stacalc=$stanorm/$scale/$scale;}
 print  sprintf("echo $ds#!sta_of_normalized_calc=%+10.9e MCPHASE_STA_OF_NORMALIZED_CALC ... sum of squared deviations (data-convolution result/scale_factor)^2$ds\n",$stacalc);
}

# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_STA=$sta\n";
#print Fout "set MCPHASE_AREADATA=$areadata\n";
#print Fout "set MCPHASE_AREACALC=$areacalc\n";
#print Fout "set MCPHASE_SCALEFACTOR=$scale\n";
#print Fout "set MCPHASE_STA_OF_NORMALIZED_CURVES=$stanorm\n";
#print Fout "set MCPHASE_STA_OF_NORMALIZED_CALC=$stacalc\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_STA=$sta\n";
#print Fout "export MCPHASE_AREADATA=$areadata\n";
#print Fout "export MCPHASE_AREACALC=$areacalc\n";
#print Fout "export MCPHASE_SCALEFACTOR=$scale\n";
#print Fout "export MCPHASE_STA_OF_NORMALIZED_CURVES=$stanorm\n";
#print Fout "export MCPHASE_STA_OF_NORMALIZED_CALC=$stacalc\n";
#close Fout;
 if ($^O=~/MSWin/){
print "set MCPHASE_STA=$sta\n";
print "set MCPHASE_AREADATA=$areadata\n";
print "set MCPHASE_AREACALC=$areacalc\n";
print "set MCPHASE_SCALEFACTOR=$scale\n";
print "set MCPHASE_STA_OF_NORMALIZED_CURVES=$stanorm\n";
print "set MCPHASE_STA_OF_NORMALIZED_CALC=$stacalc\n";
                  }
                 else
                  {
print "export MCPHASE_STA=$sta\n";
print "export MCPHASE_AREADATA=$areadata\n";
print "export MCPHASE_AREACALC=$areacalc\n";
print "export MCPHASE_SCALEFACTOR=$scale\n";
print "export MCPHASE_STA_OF_NORMALIZED_CURVES=$stanorm\n";
print "export MCPHASE_STA_OF_NORMALIZED_CALC=$stacalc\n";
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
