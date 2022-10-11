#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program swapc  used to swap columnx to columny in data file\n";

 print " usage: swapc colx coly  *.*   \n colx=columnx, coly=columny \n *.* .. filenname\n";

 exit 0;}else{print STDERR "#* $0 *";}

 

$ARGV[0]=~s/x/*/g;$colx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$coly=eval $ARGV[0];shift @ARGV;



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

           	  $i=0;$store=$numbers[$coly-1];

		       $numbers[$coly-1]=$numbers[$colx-1];

		       $numbers[$colx-1]=$store;

		  foreach (@numbers)

		  {++$i;print Fout $numbers[$i-1]." ";}     

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