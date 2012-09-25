#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1)

{print " program newcols  used to create new columns by copying contents of a given column n times\n";
 print "  .... instead of column col there will be n columns in the data file after running this program\n";
 print " usage: newcols col n  *.*   \n col=existing column \n n=number of new columns\n *.* .. filenname\n";

 exit 0;}

 

$column=$ARGV[0];shift @ARGV;
$n=$ARGV[0];shift @ARGV;
if ($column<1||$n<1){die "Error newcols: invalid column number $column or number of new columns $n\n";}


  foreach (@ARGV)

  {

   $file=$_;

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
		  print Fout $numbers[$i-1]." ";
		  if ($i==$column) { # insert column $column $n times
                   for($k=1;$k<=$n;++$k){print Fout $numbers[$i-1]." ";}
                                   }
                  }
		  if ($column>=$#numbers+2){print STDERR "Warning: column $column does not exist in dataline $j of file $file- no new column(s) created in this line\n";}
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

