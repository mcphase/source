#!/usr/bin/perl
use File::Copy;

BEGIN{@ARGV=map{glob($_)}@ARGV}
unless ($#ARGV >0)
{
print STDERR << "EOF";

 program to extract the y-value of a function y(x) by averaging
 over an interval xvalue+-dx,
 
 usage: perl getvalue.pl colx coly xvalue dx filename

 note:  colx has to be sorted
        if colx=0 then the x axis is assumed to be the line number (not considering
        comment lines)
        if you do not want any integration at interval boundaries to happen, 
        remember to put dx to zero.

 examples:

 1) from a spectrum in file spectrum.dat with Energy in column 1 and Intensity 
    in column 2 get Intensity at 
    35 meV averaging datapoint over +-1.5meV

. getvalue 1 2 35 1.5 spectrum.dat


 2) in a function stored in column 2 (x values) vs column 5 (yvalues) of file function.dat
    get y value at x=23.4 by linear interpolation

. getvalue 2 5 23.4 0 function.dat

 3) in datafile f.dat get value of 46th datapoint in column 3 

. getvalue 0 3 46 0 f.dat

 

 output: y-value           to stdout and stored in env. variable MCPHASE_YVALUE
         1/y-value         to stdout and stored in MCPHASE_YVALUE_INVERSE
         standarddeviation to stdout and stored in MCPHASE_STA
EOF
# clean bat files
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");close Fout;
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");close Fout;
exit(1);
}else{print STDERR "#* $0 *\n";}
$ARGV[0]=~s/x/*/g;$colx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$coly=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$xvalue=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$dE=eval $ARGV[0];shift @ARGV; 
foreach(@ARGV)
{$filename=$_;
($yvalue,$sta)=getvalue_by_averaging_over_intervaldE($xvalue,$colx,$coly,$dE,$filename);
if (abs($yvalue)>1e-300){$yinv=1/$yvalue;}else{$yinv=" ";}
print "echo \"#! in colx= $colx  coly = $coly of  $filename the xvalue=$xvalue +- dx=$dE corresponds\"\n";
print "echo \"#! to the yvalue=$yvalue  (1/yvalue=$yinv)";if($sta>0){print "deviations sta=$sta\"\n";}else{print"\"\n";}
} 
# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_YVALUE=$yvalue\n";
#print Fout sprintf("set MCPHASE_YVALUE_ROUNDED_INT=%.0f\n",$yvalue);
#print Fout "set MCPHASE_YVALUE_INVERSE=$yinv\n";
#print Fout sprintf("set MCPHASE_YVALUE_INVERSE_ROUNDED_INT=%.0f\n",$yinv);
#print Fout "set MCPHASE_STA=$sta\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_YVALUE=$yvalue\n";
#print Fout sprintf("export MCPHASE_YVALUE_ROUNDED_INT=%.0f\n",$yvalue);
#print Fout "export MCPHASE_YVALUE_INVERSE=$yinv\n";
#print Fout sprintf("export MCPHASE_YVALUE_INVERSE_ROUNDED_INT=%.0f\n",$yinv);
#print Fout "export MCPHASE_STA=$sta\n";
#printclose Fout;

 if ($^O=~/MSWin/){
print  "set MCPHASE_YVALUE=$yvalue\n";
print  sprintf("set MCPHASE_YVALUE_ROUNDED_INT=%.0f\n",$yvalue);
print  "set MCPHASE_YVALUE_INVERSE=$yinv\n";
print  sprintf("set MCPHASE_YVALUE_INVERSE_ROUNDED_INT=%.0f\n",$yinv);
print  "set MCPHASE_STA=$sta\n";
                  }
                 else
                  {
print "export MCPHASE_YVALUE=$yvalue\n";
print sprintf("export MCPHASE_YVALUE_ROUNDED_INT=%.0f\n",$yvalue);
print "export MCPHASE_YVALUE_INVERSE=$yinv\n";
print sprintf("export MCPHASE_YVALUE_INVERSE_ROUNDED_INT=%.0f\n",$yinv);
print "export MCPHASE_STA=$sta\n";
                  }



exit(0);

sub getvalue_by_averaging_over_intervaldE { 
#integrates intensity between $constx+-$dE
my ($constx,$colx,$coly,$dE,$file)=@_;
  unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
  $Iav=0;$j=0;$esum=0;$nofpoints=0;$order=1;
  while($line=<Fin>){
   unless($line=~/^\s*#/||$line=~/^\s*\n/)
   {$line=~s/D/E/g;@numbers=split(" ",$line);
    unshift(@numbers,$j+1); # put into first column the line number 
    if ($j==0){@numbers1=@numbers;}else{if($order==-1&&$numbers[$colx]>$numbers1[$colx]){die "Error getvalue: column $colx not sorted\n";}
             $order=1;if($numbers[$colx]<$numbers1[$colx]){$order=-1;}} # determine if ascending or descending data in colx
            ++$j; # read dataline
   unless(0==($numbers[$colx]-$numbers1[$colx])||
              $constx-$dE>$numbers[$colx]||
              $constx-$dE>$numbers1[$colx]||
              $constx+$dE<$numbers[$colx]||
              $constx+$dE<$numbers1[$colx])  # check if dataline and previous dataline are in range constx +- dE
             {++$nofpoints;
  $Iav+=$order*($numbers[$coly]+$numbers1[$coly])/2*($numbers[$colx]-$numbers1[$colx]);
  $esum+=$order*($numbers[$colx]-$numbers1[$colx]);
       	     } 
# now treat a boundary point correctly (i.e. either dataline or previous datalin is not in range constx+- dE)
unless(0==($numbers[$colx]-$numbers1[$colx]))
   { if($dE==0) # dE=0 means interpolate between points in datafile
     {if(($order*$numbers[$colx]>=$order*$constx)&($order*$constx>=$numbers1[$colx]*$order)) # interpolate if constx in interval
      {$esum=1;$Iav=$numbers1[$coly]+($numbers[$coly]-$numbers1[$coly])*($constx-$numbers1[$colx])/($numbers[$colx]-$numbers1[$colx]);
      }# print $esum.$Iav.$numbers[$colx].$numbers1[$colx].$constx."\n";
     }
     else { 
    if(($order*($constx-$dE)>$order*($numbers[$colx]))&($order*($constx-$dE)<$order*($numbers1[$colx])))
                  {++$nofpoints;# print STDERR "#hello1";
                   $bx=$constx-$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=-($by+$numbers1[$coly])/2*($bx-$numbers1[$colx]);
                   $esum+=-$order*($bx-$numbers1[$colx]);
                  }
                 if(($order*($constx-$dE)>$order*($numbers1[$colx]))&($order*($constx-$dE)<$order*($numbers[$colx])))
                  {++$nofpoints;# print STDERR "#hello2";
                   $bx=$constx-$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=-($by+$numbers[$coly])/2*($bx-$numbers[$colx]);
                   $esum+=-$order*($bx-$numbers[$colx]);
                  }
                 if(($order*($constx+$dE)>$order*($numbers[$colx]))&($order*($constx+$dE)<$order*($numbers1[$colx])))
                  {++$nofpoints;# print STDERR "#hello3";
                   $bx=$constx+$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=($by+$numbers[$coly])/2*($bx-$numbers[$colx]);
                   $esum+=$order*($bx-$numbers[$colx]);
                  }
                 if(($order*($constx+$dE)>$order*($numbers1[$colx]))&($order*($constx+$dE)<$order*($numbers[$colx])))
                  {++$nofpoints;# print STDERR "#hello4";
                   $bx=$constx+$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=($by+$numbers1[$coly])/2*($bx-$numbers1[$colx]);
                   $esum+=$order*($bx-$numbers1[$colx]);
                  }
                 }}
   @numbers1=@numbers;
   }}
  close Fin;
  if($j==1&&$numbers[$colx]==$constx){# only one point in file and specified, special ...
            $esum=1;$Iav=$numbers[$coly];
           }
  if (abs($esum)<1e-10){print STDERR " first xvalue sum on averaging is zero for $file nofpoints=$nofpoints- maybe $constx out of range of x values\n";<stdin>;}
 #print STDERR "# j=$j Iav=$Iav esum=$esum\n";  
 $Iav/=$esum; 
  my $sta=0; # here calculate sta (scattering of data in interval dE)
  unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
  $j=0;$esum=0;
  while($line=<Fin>){
   unless($line=~/^\s*#/||$line=~/^\s*\n/){$line=~s/D/E/g;@numbers=split(" ",$line);unshift(@numbers,$j+1);
         if ($j==0){@numbers1=@numbers;}else{$order=1;if($numbers[$colx]<$numbers1[$colx]){$order=-1;}}
            ++$j;
                 unless(0==($numbers[$colx]-$numbers1[$colx])||$constx-$dE>$numbers[$colx]||$constx-$dE>$numbers1[$colx]
                                                                 ||$constx+$dE<$numbers[$colx]||$constx+$dE<$numbers1[$colx])
                  {
  my $d=($Iav-($numbers[$coly]+$numbers1[$coly])/2);
  $sta+=$d*$d*$order*($numbers[$colx]-$numbers1[$colx]);
  $esum+=$order*($numbers[$colx]-$numbers1[$colx]);
       	          }
#print "$constx $dE ".$numbers[$colx]." ".$numbers1[$colx]."\n";
   @numbers1=@numbers;
   }} 
  if (abs($esum)<1e-300){print "echo \"# getvalue: xvalues variance on averaging is too small ($esum) for $file\"\n";$sta=-1;}
  else{$sta/=$esum;}
  close Fin;
  return ($Iav,$sta);
}

# **********************************************************************************************
# extracts variable from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# ($standarddeviation)=extract("sta","data.dat");
#
# ... it stores 0.24 in the variable $standarddeviation
#
sub extract {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^.*$var\s*=/) {($value)=($line=~m|$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);}                                        }
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
