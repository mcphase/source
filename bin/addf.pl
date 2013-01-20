#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program addf  used to add column1 in data file1 and column 2 in data file2 of same length\n";

 print " the results goes to col2 in datafile 2\n";

 print " usage: addf col1 file1 col2 file2   \n col1=column 1, col2 =column2 \n file1, file2 .. filenname\n";

 exit 0;}

 

$ARGV[0]=~s/x/*/g;$col1=eval $ARGV[0];shift @ARGV;

$file1=$ARGV[0];shift @ARGV;

$ARGV[0]=~s/x/*/g;$col2=eval $ARGV[0];shift @ARGV;

$file2=$ARGV[0];shift @ARGV;



     unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}   
     unless (open (Fin2, $file2)){die "\n error:unable to open $file2\n";}   
   

   open (Fout2, ">swapf2.out");

   while($line1=<Fin1>)

     {

       if ($line1=~/^\s*#/) {}

       else{

              while(($line2=<Fin2>)=~/^\s*#/){print Fout2 $line2;}

             $line1=~s/D/E/g;$line2=~s/D/E/g;

             @numbers1=split(" ",$line1);

             @numbers2=split(" ",$line2);

	     

		       $numbers2[$col2-1]+=$numbers1[$col1-1];

		  $i=0;

		  foreach (@numbers2)

		  {++$i;print Fout2 $numbers2[$i-1]." ";}     

            print Fout2 "\n";

           }

      }

      close Fin1;

      close Fin2;

      close Fout2;



       unless (rename "swapf2.out",$file2)

         {unless(open (Fout, ">$file"))     

      {die "\n error:unable to write to $file\n";}

      open (Fin, "swapf2.out");

      while($line=<Fin>){ print Fout $line;}

      close Fin;

      close Fout;

      system "del swapf2.out"; 

     }





#\end{verbatim} 