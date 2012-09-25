#!/usr/bin/perl
use File::Copy;

BEGIN{@ARGV=map{glob($_)}@ARGV}
unless ($#ARGV >0)
{
print STDOUT << "EOF";

 program to get the y-value of a function by averaging
 over an interval xvalue+-dx,
 note: colx has to be sorted

 usage: perl getvalue.pl colx coly xvalue dx filename

 note: if colx=0 then the x axis is assumed to be the line number (not considering
       comment lines)
 output: the y-value is written to stdout and to env. varaibale MCPHASE_YVALUE
         1/y-value is written to stdout MCPHASE_YVALUE_INVERSE
         standarddeviation to stdaout and MCPHASE_STA
EOF
# clean bat files
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");close Fout;
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");close Fout;
exit(1);
}

$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;
$xvalue=$ARGV[0];shift @ARGV;
$dE=$ARGV[0];shift @ARGV; 
foreach(@ARGV)
{$filename=$_;
($yvalue,$sta)=getvalue_by_averaging_over_intervaldE($xvalue,$colx,$coly,$dE,$filename);
if (abs($yvalue)>1e-300){$yinv=1/$yvalue;}else{$yinv=" ";}
print "#! in colx= $colx  coly = $coly of  $filename the xvalue=$xvalue +- dx=$dE corresponds\n";
print "#! to the yvalue=$yvalue  (1/yvalue=$yinv)";if($sta>0){print "deviations sta=$sta\n";}else{print"\n";}
} 
# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_YVALUE=$yvalue\n";
print Fout sprintf("set MCPHASE_YVALUE_ROUNDED_INT=%.0f\n",$yvalue);
print Fout "set MCPHASE_YVALUE_INVERSE=$yinv\n";
print Fout sprintf("set MCPHASE_YVALUE_INVERSE_ROUNDED_INT=%.0f\n",$yinv);
print Fout "set MCPHASE_STA=$sta\n";
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_YVALUE=$yvalue\n";
print Fout sprintf("export MCPHASE_YVALUE_ROUNDED_INT=%.0f\n",$yvalue);
print Fout "export MCPHASE_YVALUE_INVERSE=$yinv\n";
print Fout sprintf("export MCPHASE_YVALUE_INVERSE_ROUNDED_INT=%.0f\n",$yinv);
print Fout "export MCPHASE_STA=$sta\n";
close Fout;

exit(0);

sub getvalue_by_averaging_over_intervaldE { #integrates intensity between $constx+-$dE
my ($constx,$colx,$coly,$dE,$file)=@_;
  unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
  $Iav=0;$j=0;$esum=0;$nofpoints=0;$order=1;
  while($line=<Fin>){
   unless($line=~/^\s*#/||$line=~/^\s*\n/){$line=~s/D/E/g;@numbers=split(" ",$line);unshift(@numbers,$j+1);
         if ($j==0){@numbers1=@numbers;}else{$order=1;if($numbers[$colx]<$numbers1[$colx]){$order=-1;}}
            ++$j;
                 unless(0==($numbers[$colx]-$numbers1[$colx])||$constx-$dE>$numbers[$colx]||$constx-$dE>$numbers1[$colx]
                                                                 ||$constx+$dE<$numbers[$colx]||$constx+$dE<$numbers1[$colx])
                  {++$nofpoints;
  $Iav+=$order*($numbers[$coly]+$numbers1[$coly])/2*($numbers[$colx]-$numbers1[$colx]);
  $esum+=$order*($numbers[$colx]-$numbers1[$colx]);
       	          } # now treat a boundary point correctly
                 if(($order*($constx-$dE)>$order*($numbers[$colx]))&($order*($constx-$dE)<$order*($numbers1[$colx])))
                  {++$nofpoints;#print "hello1";
                   $bx=$constx-$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=-($by+$numbers1[$coly])/2*($bx-$numbers1[$colx]);
                   $esum+=-$order*($bx-$numbers1[$colx]);
                  }
                 if(($order*($constx-$dE)>$order*($numbers1[$colx]))&($order*($constx-$dE)<$order*($numbers[$colx])))
                  {++$nofpoints;#print "hello2";
                   $bx=$constx-$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=-($by+$numbers[$coly])/2*($bx-$numbers[$colx]);
                   $esum+=-$order*($bx-$numbers[$colx]);
                  }
                 if(($order*($constx+$dE)>$order*($numbers[$colx]))&($order*($constx+$dE)<$order*($numbers1[$colx])))
                  {++$nofpoints;#print "hello3";
                   $bx=$constx+$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=($by+$numbers[$coly])/2*($bx-$numbers[$colx]);
                   $esum+=$order*($bx-$numbers[$colx]);
                  }
                 if(($order*($constx+$dE)>$order*($numbers1[$colx]))&($order*($constx+$dE)<$order*($numbers[$colx])))
                  {++$nofpoints;#print "hello4";
                   $bx=$constx+$dE;
                   $by=$numbers[$coly]+($bx-$numbers[$colx])*($numbers1[$coly]-$numbers[$coly])/($numbers1[$colx]-$numbers[$colx]);
                   $Iav+=($by+$numbers1[$coly])/2*($bx-$numbers1[$colx]);
                   $esum+=$order*($bx-$numbers1[$colx]);
                  }
   @numbers1=@numbers;
   }}
  close Fin;
  if (abs($esum)<1e-10){print "\n first xvalue sum on averaging is zero for $file nofpoints=$nofpoints- maybe $constx out of range of x values\n";<stdin>;}
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
  if (abs($esum)<1e-300){print "\n getvalue: xvalues variance on averaging is too small ($esum) for $file\n";$sta=-1;}
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
