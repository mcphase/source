#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program form used to reformat numbers in columns using a special number format\n\n";

 print " usage: fform col1 col2 format *.*   \n\n col1=1st column,\n col2=last column to be formatted\n format=output format such as 4.4g \n *.* .. filenname\n";

 exit 0;}

 

$col1=$ARGV[0];shift @ARGV;

$col2=$ARGV[0];shift @ARGV;

$format=$ARGV[0];shift @ARGV;



  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

  open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

           	  $i=0;++$j;

		  foreach (@numbers)

		  {++$i;

		  if ($i>=$col1&&$i<=$col2) 

                {print Fout sprintf("%".$format." ",$numbers[$i-1]);}

		  else

                {print Fout $numbers[$i-1]." ";}     

              }

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



   print ">\n";

   }



#\end{verbatim} 