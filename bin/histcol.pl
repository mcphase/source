#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

unless ($#ARGV >1) 

{print " program histcol  used to generate a histogram of a column in a data file and writes it\n";
 print " to histcol.out and stdout. Average and sum of squared deviations per point are calculated.\n";
 print " usage: histcol column step  *.*   \n column=column, step=stepwidth of histogram points (or -n 100 for 100 steps)\n *.* .. filenname\n";

 exit 0;}
print "# $0 @ARGV\n";
$ARGV[0]=~s/x/*/g;$col=eval $ARGV[0];shift @ARGV;
$step=$ARGV[0];shift @ARGV;
if ($step=~/-n/){$ARGV[0]=~s/x/*/g;$nofsteps=eval $ARGV[0];shift @ARGV; }
         else {$step=eval $step;}

# get minimum and maximum and calculate average
$min=1e100;$max=-1e100;$average=0;$nn=0;
  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error histcol:unable to open $file\n";}
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);++$nn;$average+=$numbers[$col-1];
		  if ($max<$numbers[$col-1]){$max=$numbers[$col-1];}
    		  if ($min>$numbers[$col-1]){$min=$numbers[$col-1];}
            }
      }
   close Fin;
   }
$average/=$nn;

if ($step=~/-n/){$step=($max-$min)/$nofsteps;}
if($max<=$min){die "Error histcol $file: maximum equal or less than minimum\n";}
if ($step/($max-$min)<1e-3) {die "Error histcol $file: not more than 1000 steps allowed \n";}
# determine histogram
@histo=();
 # histogramm steps (not more than 1000)
for($hx=0;$hx<=int(($max-$min)/$step)+1;++$hx){$histo[$hx]=0;}
$sta=0;
  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error histcol:unable to open $file\n";}
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);$sta+=($average-$numbers[$col-1])*($average-$numbers[$col-1]);
		$hx=int(($numbers[$col-1]-$min)/$step);
                ++$histo[$hx];
            }
      }
   close Fin;
   }
   $sta/=$nn;

   open(Fout,">histcol.out");
   print Fout "#{Histogram of column $col in file(s) @ARGV\n";
   print STDOUT "#{Histogram of column $col in file(s) @ARGV\n";
   for($hx=0;$hx<=int(($max-$min)/$step)+1;++$hx)
   {print Fout (($hx+0.5)*$step+$min)."   ".($histo[$hx])."\n";
    print STDOUT (($hx+0.5)*$step+$min)."   ".($histo[$hx])."\n";
   } 
   print Fout "#! AVERAGE=$average (sum of squared deviations)/(number of points) STAPP=$sta\n";
   print STDOUT "#! AVERAGE=$average (sum of squared deviations)/(number of points) STAPP=$sta\n";
   close Fout;


#\end{verbatim} 