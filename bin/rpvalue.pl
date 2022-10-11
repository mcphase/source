#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print STDERR " program rpvalue  used to calculate the rpvalue from 2 columns in a file\n";
 print STDERR " usage: rpvalue col1 col2  *.*   \n col=columns \n *.* .. filenname\n";
 print STDERR " the rpvalue is defined as\n 100*[sum_{allpoints} abs(col2-col1)]/[sum_{allpoints}abs(col1)]\n";
 print STDERR " the result ist printed to stdout and stored in environment variabele MCPHASE_RP\n";

 exit 0;}else{print STDERR "#* $0 *\n";}

 

$ARGV[0]=~s/x/*/g;$col1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$col2=eval $ARGV[0];shift @ARGV;



  foreach (@ARGV)

  {$nofpoints=0;$rpvalue=0;

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "echo \"<".$file;

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

            $rpvalue+=abs($numbers[$col1-1]-$numbers[$col2-1]);

            $nofpoints+=abs($numbers[$col1-1]);

           }

     } 

      $rpvalue/=$nofpoints/100;

      print " rpvalue=".$rpvalue." ";

      close Fin;

   print ">\"\n";

   }
# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_RP=$rpvalue\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_RP=$rpvalue\n";
#close Fout;

 if ($^O=~/MSWin/){
print "set MCPHASE_RP=$rpvalue\n";
                  }
                 else
                  {
print "export MCPHASE_RP=$rpvalue\n";
                  }



#\end{verbatim} 

