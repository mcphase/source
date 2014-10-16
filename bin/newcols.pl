#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1)

{print " program newcols  used to create new columns by copying contents of a given column n times\n";
 print "  .... instead of column col there will be n columns in the data file after running this program\n";
 print " usage: newcols col n [options] *.*   \n col=existing column \n n=number of new columns\n";
print " *.* .. filenname\noptions: -c 12.3  ... all the new columns will be filled with  the constant 12.3\n";
print "-n ... all the new columns will be filled with the line number\n";
 exit 0;}

 

$ARGV[0]=~s/x/*/g;$column=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$n=eval $ARGV[0];shift @ARGV;
if ($column<1||$n<1){die "Error newcols: invalid column number $column or number of new columns $n\n";}
$dc=0;if ($ARGV[0]=~m/-c/){$dc=1;shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$c=eval $ARGV[0];shift @ARGV;}
$dn=0;if ($ARGV[0]=~m/-n/){$dn=1;shift @ARGV;}

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
                   for($k=1;$k<=$n;++$k){if($dc){print Fout $c." ";}
                                         elsif($dn){print Fout $j." ";}
                                         else{print Fout $numbers[$i-1]." ";}
                                        }
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

