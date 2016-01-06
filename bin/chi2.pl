#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print STDERR "# program chi2  used to calculate the chi-squared from 3 columns in a file\n";
 print STDERR "# if col1, col2 and col3 are calculation, experiment and experimental error, respectively, then\n";
 print STDERR "# chisquared is defined as\n# 1/nofpoints *sum_i{allpoints} (col2 - col1)^2/col3^2\n";
 print STDERR "# for each data point a line sta= (col2 - col1)^2  col3^2 is output to stdout \n";
 print STDERR "# and written to the environment variable MCPHASE_CHI2 \n";
 print STDERR "# a column containing the value of chisquared during each step of the\n";
 print STDERR "# summation is added to the file\n";
 print STDERR "# usage: chi2 col1 col2 col3  *.*   \n# col=columns \n# *.* .. filenname\n";
 exit 0;}

 

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$col1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$col2=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$col3=eval $ARGV[0];shift @ARGV;

$nofpoints=0;$chi2=0;
  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "echo \"<".$file.">\"\n";
   open (Fout, ">range.out");
   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);
            $s=($numbers[$col1-1]-$numbers[$col2-1]);
            $e=$numbers[$col3-1];
            $s2=$s*$s;$e2=$e*$e;
            $chi2+=$s2/$e2;
            print "echo \"#! sta= $s2 $e2\"\n";
            $nofpoints+=1;
            $i=0;
             foreach (@numbers)
             {++$i;
  		   print Fout $numbers[$i-1]." ";}     
              print Fout ($chi2/$nofpoints)."\n";
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



   }
      $chi2/=$nofpoints;

      print "echo \"#! chi2=".$chi2."\"\n";
# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_CHI2=$chi2\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_CHI2=$chi2\n";
#close Fout;

 if ($^O=~/MSWin/){
print "set MCPHASE_CHI2=$chi2\n";
                  }
                 else
                  {
print "export MCPHASE_CHI2=$chi2\n";
                  }

      


#\end{verbatim} 

