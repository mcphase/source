#!/usr/bin/perl
use File::Copy;

BEGIN{@ARGV=map{glob($_)}@ARGV}
unless ($#ARGV >0)
{
print STDOUT << "EOF";

 program to get the y-value of a function by averaging
 over an interval xvalue+-dE,
 note: colx has to be sorted

 usage: perl getvalue.pl colx coly xvalue dE filename

 output: the y-value is written to stdout and environment variable MCPHASE_YVALUE
         1/y-value is written to stdout MCPHASE_YVALUE_INVERSE
         standarddeviation to stdaout and MCPHASE_STA
EOF
exit(1);
}

$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;
$xvalue=$ARGV[0];shift @ARGV;
$dE=$ARGV[0];shift @ARGV;
foreach(@ARGV)
{$filename=$_;
($yvalue,$sta)=getvalue_by_averaging_over_intervaldE($xvalue,$colx,$coly,$dE,$filename);
if (abs($yvalue)>1e-10){$yinv=1/$yvalue;}else{$yinv=" ";}
print "#! in colx= $colx  coly = $coly of  $filename the xvalue=$xvalue +- dE=$dE corresponds\n";
print "#! to the yvalue=$yvalue  (1/yvalue=$yinv) deviations sta=$sta\n";
}
# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_YVALUE=$yvalue\n";
print Fout "set MCPHASE_YVALUE_INVERSE=$yinv\n";
print Fout "set MCPHASE_STA=$sta\n";
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_YVALUE=$yvalue\n";
print Fout "export MCPHASE_YVALUE_INVERSE=$yinv\n";
print Fout "export MCPHASE_STA=$sta\n";
close Fout;

exit(0);

sub getvalue_by_averaging_over_intervaldE { #integrates intensity between $constx+-$dE
my ($constx,$colx,$coly,$dE,$file)=@_;
  unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
  $Iav=0;$j=0;$esum=0;$nofpoints=0;$order=1;
  while($line=<Fin>){
   unless($line=~/^\s*#/){$line=~s/D/E/g;@numbers=split(" ",$line);
         if ($j==0){@numbers1=@numbers;}else{$order=1;if($numbers[$colx-1]<$numbers1[$colx-1]){$order=-1;}}
            ++$j;
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1])||$constx-$dE>$numbers[$colx-1]||$constx-$dE>$numbers1[$colx-1]
                                                                 ||$constx+$dE<$numbers[$colx-1]||$constx+$dE<$numbers1[$colx-1])
                  {++$nofpoints;
  $Iav+=$order*($numbers[$coly-1]+$numbers1[$coly-1])/2*($numbers[$colx-1]-$numbers1[$colx-1]);
  $esum+=$order*($numbers[$colx-1]-$numbers1[$colx-1]);
       	          } # now treat a boundary point correctly
                 if(($order*($constx-$dE)>$order*($numbers[$colx-1]))&($order*($constx-$dE)<$order*($numbers1[$colx-1])))
                  {++$nofpoints;#print "hello1";
                   $bx=$constx-$dE;
                   $by=$numbers[$coly-1]+($bx-$numbers[$colx-1])*($numbers1[$coly-1]-$numbers[$coly-1])/($numbers1[$colx-1]-$numbers[$colx-1]);
                   $Iav+=-($by+$numbers1[$coly-1])/2*($bx-$numbers1[$colx-1]);
                   $esum+=-$order*($bx-$numbers1[$colx-1]);
                  }
                 if(($order*($constx-$dE)>$order*($numbers1[$colx-1]))&($order*($constx-$dE)<$order*($numbers[$colx-1])))
                  {++$nofpoints;#print "hello2";
                   $bx=$constx-$dE;
                   $by=$numbers[$coly-1]+($bx-$numbers[$colx-1])*($numbers1[$coly-1]-$numbers[$coly-1])/($numbers1[$colx-1]-$numbers[$colx-1]);
                   $Iav+=-($by+$numbers[$coly-1])/2*($bx-$numbers[$colx-1]);
                   $esum+=-$order*($bx-$numbers[$colx-1]);
                  }
                 if(($order*($constx+$dE)>$order*($numbers[$colx-1]))&($order*($constx+$dE)<$order*($numbers1[$colx-1])))
                  {++$nofpoints;#print "hello3";
                   $bx=$constx+$dE;
                   $by=$numbers[$coly-1]+($bx-$numbers[$colx-1])*($numbers1[$coly-1]-$numbers[$coly-1])/($numbers1[$colx-1]-$numbers[$colx-1]);
                   $Iav+=($by+$numbers[$coly-1])/2*($bx-$numbers[$colx-1]);
                   $esum+=$order*($bx-$numbers[$colx-1]);
                  }
                 if(($order*($constx+$dE)>$order*($numbers1[$colx-1]))&($order*($constx+$dE)<$order*($numbers[$colx-1])))
                  {++$nofpoints;#print "hello4";
                   $bx=$constx+$dE;
                   $by=$numbers[$coly-1]+($bx-$numbers[$colx-1])*($numbers1[$coly-1]-$numbers[$coly-1])/($numbers1[$colx-1]-$numbers[$colx-1]);
                   $Iav+=($by+$numbers1[$coly-1])/2*($bx-$numbers1[$colx-1]);
                   $esum+=$order*($bx-$numbers1[$colx-1]);
                  }

   @numbers1=@numbers;
   }}
  close Fin;
  if (abs($esum)<1e-10){print "\n firstenergy sum on averaging is zero for $file nofpoints=$nofpoints- maybe $constx out of range of energy values\n";<stdin>;}
  $Iav/=$esum;
  my $sta=0;
  unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
  $j=0;$esum=0;
  while($line=<Fin>){
   unless($line=~/^\s*#/){$line=~s/D/E/g;@numbers=split(" ",$line);
         if ($j==0){@numbers1=@numbers;}else{$order=1;if($numbers[$colx-1]<$numbers1[$colx-1]){$order=-1;}}
            ++$j;
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1])||$constx-$dE>$numbers[$colx-1]||$constx-$dE>$numbers1[$colx-1]
                                                                 ||$constx+$dE<$numbers[$colx-1]||$constx+$dE<$numbers1[$colx-1])
                  {
  my $d=($Iav-($numbers[$coly-1]+$numbers1[$coly-1])/2);
  $sta+=$d*$d*$order*($numbers[$colx-1]-$numbers1[$colx-1]);
  $esum+=$order*($numbers[$colx-1]-$numbers1[$colx-1]);
       	          }
#print "$constx $dE ".$numbers[$colx-1]." ".$numbers1[$colx-1]."\n";
   @numbers1=@numbers;
   }}
  if (abs($esum)<1e-10){print "\n energy sum on averaging is zero for $file\n";<stdin>;}
  close Fin;$sta/=$esum;
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
