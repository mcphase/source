#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program addcolerr  used to addm columnx to columny and corresponding error columns in data file, the result goes to coly and colyerr\n";

 print " usage: addcol colx colxerr coly colyerr  *.*   \n colx=columnx,colxerr= error corresponding to column x,
 coly=columny, colyerr= error corresponding to column y \n *.* .. filenname\n";
 print " Error addition is done by colyerr=sqrt(colxerr^2+colyerr^2)\n";

 exit 0;}




$colx=$ARGV[0];shift @ARGV;
$colxerr=$ARGV[0];shift @ARGV;

$coly=$ARGV[0];shift @ARGV;
$colyerr=$ARGV[0];shift @ARGV;


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

           	  $i=0;$numbers[$coly-1]+=$numbers[$colx-1];
                       $numbers[$colyerr-1]=sqrt($numbers[$colxerr-1]*$numbers[$colxerr-1]+$numbers[$colyerr-1]*$numbers[$colyerr-1]);
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

   print ">\n";

   }



#\end{verbatim} 

