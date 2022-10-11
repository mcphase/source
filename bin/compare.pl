#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV ==1) 

{print "\n program compare  \n\n used to compare data file1 and data file2\n";

 print " all columns and rows are compared and a standard deviation\n"; 

 print " is calculated accoding to sum_i (file1_i - file2_i)^2\n";

 print " this standard deviation is output to stdout, e.g. as sta=143.3\n\n";

 print " usage: compare file1 file2   \n\n file1, file2 .. filennames\n";

 exit 0;}else{print STDERR "#* $0 *\n";}

 

$file1=$ARGV[0];shift @ARGV;

$file2=$ARGV[0];



   unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}   



   unless (open (Fin2, $file2)){die "\n error:unable to open $file2\n";}   



   $sta=0;

   while($line1=<Fin1>)

     {

       if ($line1=~/^\s*#/) {}

       else{

              while(($line2=<Fin2>)=~/^\s*#/){}

             $line1=~s/D/E/g;$line2=~s/D/E/g;

             @numbers1=split(" ",$line1);

             @numbers2=split(" ",$line2);

		  $i=0;

		  foreach (@numbers2)

		  {++$i;

               if (i<=$#numbers1+1){$dd=$numbers2[$i-1]-$numbers1[$i-1];$sta+=$dd*$dd;}

               }     

           }

      }

      close Fin1;

      close Fin2;



print "sta=$sta\n";

     





#\end{verbatim} 