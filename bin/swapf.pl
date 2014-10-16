#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program swapf  used to swap column1 in data file1 and column 2 in data file2 of same length\n";

 print " usage: swapc col1 file1 col2 file2   \n col1=column 1, col2 =column2 \n file1, file2 .. filenname\n";

 exit 0;}

 

$ARGV[0]=~s/x/*/g;$col1=eval $ARGV[0];shift @ARGV;
$file1=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$col2=eval $ARGV[0];shift @ARGV;
$file2=$ARGV[0];shift @ARGV;



   unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}   
   unless (open (Fin2, $file2)){die "\n error:unable to open $file2\n";}   

   
   open (Fout1, ">swapf1.out");

   open (Fout2, ">swapf2.out");

   while($line1=<Fin1>)

     {

       if ($line1=~/^\s*#/) {print Fout1 $line1;}

       else{

              while(($line2=<Fin2>)=~/^\s*#/){print Fout2 $line2;}



             $line1=~s/D/E/g;@numbers1=split(" ",$line1);

             $line2=~s/D/E/g;@numbers2=split(" ",$line2);

	     

           	       $store=$numbers1[$col1-1];

		       $numbers1[$col1-1]=$numbers2[$col2-1];

		       $numbers2[$col2-1]=$store;

		  $i=0;

		  foreach (@numbers1)

		  {++$i;print Fout1 $numbers1[$i-1]." ";}     

		  $i=0;

		  foreach (@numbers2)

		  {++$i;print Fout2 $numbers2[$i-1]." ";}     

            print Fout1 "\n";

            print Fout2 "\n";

           }

      }

      close Fin1;

      close Fin2;

      close Fout1;

      close Fout2;

       unless (rename "swapf1.out",$file1)

          {unless(open (Fout, ">$file1"))     

      {die "\n error:unable to write to $file1\n";}

      open (Fin, "swapf1.out");

      while($line=<Fin>){ print Fout $line;}

      close Fin;

      close Fout;

      system "del swapf1.out"; 

     }



       unless (rename "swapf2.out",$file2)

         {unless(open (Fout, ">$file2"))     

      {die "\n error:unable to write to $file2\n";}

      open (Fin, "swapf2.out");

      while($line=<Fin>){ print Fout $line;}

      close Fin;

      close Fout;

      system "del swapf2.out"; 

     }





#\end{verbatim} 