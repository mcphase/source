#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}
use Time::Local;
#\begin{verbatim}



unless ($#ARGV >1) 

{print " program form used to reformat numbers in columns using a special number format

 usage: fform col1 col2 format *.*   

 col1=1st column,
 col2=last column to be formatted
 format=output format such as 4.4g or -t, -t means that col1 to col2 are taken to be
        year month day hour min second  and the column with 
        second is replaced by a timestamp in seconds
   
 *.* .. filenname\n";

 exit 0;}else{print STDERR "#* $0 *";}

 

$ARGV[0]=~s/x/*/g;$col1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$col2=eval $ARGV[0];shift @ARGV;

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

                {unless($format=~/-t/){print Fout sprintf("%".$format." ",$numbers[$i-1]);}
                       else
                      {if($i!=$col1+5){print Fout $numbers[$i-1]." ";}
                       else
                      {$sec=$numbers[$col1+5-1];
                       $min=$numbers[$col1+4-1];
                       $hour=$numbers[$col1+3-1];
                       $mday=$numbers[$col1+2-1];
                       $mon=$numbers[$col1+1-1];
                       $year=$numbers[$col1+0-1];
                       print Fout timelocal($sec,$min,$hour,$mday,$mon-1,$year)." ";}
                      }

                }
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



   print ">";

   }
print "\n";


#\end{verbatim} 