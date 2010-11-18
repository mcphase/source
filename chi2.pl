#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program chi2  used to calculate the chi-squared from 3 columns in a file\n";
 print " if col1, col2 and col3 are calculation, experiment and experimental error, respectively, then\n";
 print " chisquared is defined as\n 1/nofpoints *sum_i{allpoints} (col2 - col1)^2/col3^2\n";
 print " usage: chi2 col1 col2 col3  *.*   \n col=columns \n *.* .. filenname\n";
 exit 0;}

 

$col1=$ARGV[0];shift @ARGV;
$col2=$ARGV[0];shift @ARGV;
$col3=$ARGV[0];shift @ARGV;


  foreach (@ARGV)

  {$nofpoints=0;$chi2=0;

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

            $chi2+=($numbers[$col1-1]-$numbers[$col2-1])*($numbers[$col1-1]-$numbers[$col2-1])/$numbers[$col3-1]/$numbers[$col3-1];

            $nofpoints+=1;

           }

     } 

      $chi2/=$nofpoints;

      print " chi2=".$chi2." ";

      close Fin;

   print ">\n";

   }



#\end{verbatim} 

