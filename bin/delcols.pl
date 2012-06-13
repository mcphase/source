#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1)

{print " program delcols  used to delete n columns in a data file\n";
 print " usage: delcols col n  *.*   \n col=first column to be deleted\n n=number of columns to be deleted\n *.* .. filenname\n";
 exit 0;}

$column=$ARGV[0];shift @ARGV;
$n=$ARGV[0];shift @ARGV;
  foreach (@ARGV)
  {   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
  open (Fout, ">range.out");
   $j=0;
   while($line=<Fin>)
     {
       if ($line=~/^\s*#/) {print Fout $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
           	  $i=0;++$j;
		  foreach (@numbers)
		  {++$i;
		  if ($i<$column||$i>$column+$n-1){print Fout $numbers[$i-1]." ";}

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