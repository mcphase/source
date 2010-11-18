#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >0) 

{print " program sum  used to sum columnx, the result goes to stdout\n";

 print " usage: sum colx *.*   \n colx=columnx \n *.* .. filenname\n";

 exit 0;}

 

$colx=$ARGV[0];shift @ARGV;





  foreach (@ARGV)

  {

   $file=$_;

    unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

                  $i=0;$sum=0;$sta=0;$abs=0;

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

                  ++$i;

                  $sum+=$numbers[$colx-1];

                  $sta+=$numbers[$colx-1]*$numbers[$colx-1];

                  $abs+=abs($numbers[$colx-1]);

           }

      }

      close Fin;

   print ">\n";


   print "sum of values sum=".$sum."\n";

   print "number of points n=".$i."\n";
   if($i<1){$i=1;}
   $stappoint=$sta/$i;
   $absppoint=$abs/$i;
   print "standard deviation sum_i (value_i*value_i) sta=".$sta."\n";

   print "standard deviation/n stappoint=".$stappoint."\n";

   print "sum of absolute values sumabs=".$abs."\n";

   print "sum of absolute values/n sumabsppoint=".$absppoint."\n";

   }



#\end{verbatim} 

