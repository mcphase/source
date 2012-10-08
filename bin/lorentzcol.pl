#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



use Math::Trig;

unless ($#ARGV >0) 

{print " program lorentzcol  used to calculate a lorentzian from a column\n";
 print " usage: lorentzcol col position fwhm area *.*   \n col=column\n *.* .. filenname\n

the formula for a lorentz curve is:
lorentz(x)=area/3.1415/fwhm/(1.0+(x-position)^2/fwhm^2)\n";

 exit 0;}

 

$column=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;
$pos=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;
$fwhm=eval{$ARGV[0]/2};shift @ARGV;
$ARGV[0]=~s/x/*/g;
$area=eval $ARGV[0];shift @ARGV;


  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
   open (Fout, ">range.out");
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);
           	  $i=0;++$j;
		  foreach (@numbers)
		  {++$i;
		  if ($i==$column) {$x=$numbers[$i-1]-$pos;
                                    $x1=$x*$x/$fwhm/$fwhm;
                                    if(abs($x1)<100){$numbers[$i-1]=$area/3.1415/$fwhm/(1.0+$x1);}
                                             else {$numbers[$i-1]=0;}
                                   }
		  print Fout $numbers[$i-1]." ";}     
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