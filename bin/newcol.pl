#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >0) 

{print " program newcol  used to create a new column (containing line numbers)\n";

 print " usage: newcol col [option] *.*   \n col=column \n *.* .. filenname\noptions: -c 12.3  ... the new column will be filled with the constant 12.3\n";

 exit 0;}else{print STDERR "#* $0 *";}

 
$dc=0;
$ARGV[0]=~s/x/*/g;$column=eval $ARGV[0];shift @ARGV;
if ($ARGV[0]=~m/-c/){$dc=1;shift @ARGV;$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$c=eval $ARGV[0];shift @ARGV;}



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

		  if ($i==$column) {if($dc){print Fout $c." ";}
                                    else{print Fout $j." ";}
                                   }

		  print Fout $numbers[$i-1]." ";}

		  if ($column==$#numbers+2){if($dc){print Fout $c." ";}
                                    else{print Fout $j." ";}
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

