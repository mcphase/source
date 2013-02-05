#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



use Math::Trig;

unless ($#ARGV >0) 

{print " program gausscol  used to calculate a gaussian from a column\n";
 print " usage: gausscol col position fwhm area *.*   \n col=column\n *.* .. filenname\n

 the formula for a gaussian is:
sigma=fwhm/2/sqrt(2*log(2))
gauss(x)=area/sqrt(2*3.1415)/sigma* exp(-(x-position)^2/2sigma^2)\n";

 exit 0;}

 

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$column=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$pos=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$fwhm=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$area=eval $ARGV[0];shift @ARGV;


  foreach (@ARGV)
  {
   $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
   open (Fout, ">range.out");
    $sigma=$fwhm/2/sqrt(2*log(2));



   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);
           	  $i=0;++$j;
		  foreach (@numbers)
		  {++$i;
		  if ($i==$column) {$x=$numbers[$i-1]-$pos;
                                    $x1=-$x*$x/$sigma/$sigma/2;
                                    if(abs($x1)<100){$numbers[$i-1]=$area/sqrt(2*3.1415)/$sigma* exp(-$x*$x/$sigma/$sigma/2);}
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