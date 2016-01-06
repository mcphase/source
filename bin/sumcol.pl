#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >0) 

{print STDERR " program sumcol  used to sum columnx in a file, the result goes to the file and total sum to stdout\n";
 print STDERR " and is written to the environment variables MCPHASE_SUM, \n";
 print STDERR " MCPHASE_STA,MCPHASE_STAPPOINT,MCPHASE_SUMABS,MCPHASE_SUMABSPPOINT\n";
 print STDERR " usage: sumcol colx *.*   \n colx=columnx \n *.* .. filenname\n";

 exit 0;}

 

$ARGV[0]=~s/x/*/g;$colx=eval $ARGV[0];shift @ARGV;





  foreach (@ARGV)

  {

   $file=$_;

    unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "echo \"<".$file;

                  $i=0;$sum=0;$sta=0;$abs=0;
   open (Fout, ">range.out");
   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

                  ++$i;

                  $sum+=$numbers[$colx-1];

                  $sta+=$numbers[$colx-1]*$numbers[$colx-1];

                  $abs+=abs($numbers[$colx-1]);
                  $numbers[$colx-1]=$sum;
                  foreach (@numbers)
		  {print Fout $_." ";}     
                  print Fout "\n";
           }

      }

      close Fin;
close Fout;

       unless (rename "range.out",$file)

          {unless(open (Fout, ">$file"))     

      {die "\n error:unable to write to $file\n";}

      open (Fin, "range.out");

      while($line=<Fin>){ print Fout $line;}

      close Fin;

      close Fout;

      system "del range.out"; 

     }

   print ">\"\n";


   print "echo \"#! sum of values sum=".$sum."\"\n";

   print "echo \"#! number of points n=".$i."\"\n";
   if($i<1){$i=1;}
   $stappoint=$sta/$i;
   $absppoint=$abs/$i;
   print "echo \"#! standard deviation sum_i (value_i*value_i) sta=".$sta."\"\n";

   print "echo \"#! standard deviation\"\necho \"#! stappoint=".$stappoint."\"\n";

   print "echo \"#! sum of absolute values sumabs=".$abs."\"\n";

   print "echo \"#! sum of absolute values\"\necho \"#! sumabsppoint=".$absppoint."\"\n";

   }

# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_SUM=$sum\n";
#print Fout "set MCPHASE_STA=$sta\n";
#print Fout "set MCPHASE_STAPPOINT=$stappoint\n";
#print Fout "set MCPHASE_SUMABS=$abs\n";
#print Fout "set MCPHASE_SUMABSPPOINT=$absppoint\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_SUM=$sum\n";
#print Fout "export MCPHASE_STA=$sta\n";
#print Fout "export MCPHASE_STAPPOINT=$stappoint\n";
#print Fout "export MCPHASE_SUMABS=$abs\n";
#print Fout "export MCPHASE_SUMABSPPOINT=$absppoint\n";
#close Fout;

 if ($^O=~/MSWin/){
print  "set MCPHASE_SUM=$sum\n";
print  "set MCPHASE_STA=$sta\n";
print  "set MCPHASE_STAPPOINT=$stappoint\n";
print  "set MCPHASE_SUMABS=$abs\n";
print  "set MCPHASE_SUMABSPPOINT=$absppoint\n";
                  }
                 else
                  {
print  "export MCPHASE_SUM=$sum\n";
print  "export MCPHASE_STA=$sta\n";
print  "export MCPHASE_STAPPOINT=$stappoint\n";
print  "export MCPHASE_SUMABS=$abs\n";
print  "export MCPHASE_SUMABSPPOINT=$absppoint\n";
                  }


#\end{verbatim} 

