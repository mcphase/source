#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >0) 

{print " program sumcol  used to sum columnx, the result goes to stdout\n";
 print " and is written to the environment variables MCPHASE_SUM, \n";
 print " MCPHASE_STA,MCPHASE_STAPPOINT,MCPHASE_SUMABS,MCPHASE_SUMABSPPOINT\n";
 print " usage: sumcol colx *.*   \n colx=columnx \n *.* .. filenname\n";

 exit 0;}

 

$ARGV[0]=~s/x/*/g;$colx=eval $ARGV[0];shift @ARGV;





  foreach (@ARGV)

  {

   $file=$_;

    unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

                  $i=0;$sum=0;$sta=0;$abs=0;

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

                  ++$i;

                  $sum+=$numbers[$colx-1];

                  $sta+=$numbers[$colx-1]*$numbers[$colx-1];

                  $abs+=abs($numbers[$colx-1]);

           }

      }

      close Fin;

   print ">\n";


   print "sum of values sum=".$sum."\n";

   print "number of points n=".$i."\n";
   if($i<1){$i=1;}
   $stappoint=$sta/$i;
   $absppoint=$abs/$i;
   print "standard deviation sum_i (value_i*value_i) sta=".$sta."\n";

   print "standard deviation/n stappoint=".$stappoint."\n";

   print "sum of absolute values sumabs=".$abs."\n";

   print "sum of absolute values/n sumabsppoint=".$absppoint."\n";

   }

# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_SUM=$sum\n";
print Fout "set MCPHASE_STA=$sta\n";
print Fout "set MCPHASE_STAPPOINT=$stappoint\n";
print Fout "set MCPHASE_SUMABS=$abs\n";
print Fout "set MCPHASE_SUMABSPPOINT=$absppoint\n";
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_SUM=$sum\n";
print Fout "export MCPHASE_STA=$sta\n";
print Fout "export MCPHASE_STAPPOINT=$stappoint\n";
print Fout "export MCPHASE_SUMABS=$abs\n";
print Fout "export MCPHASE_SUMABSPPOINT=$absppoint\n";
close Fout;

#\end{verbatim} 

