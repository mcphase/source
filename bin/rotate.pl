#!/usr/bin/perl

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program rotate  used to rotate coordinate axes\n";
  print STDERR << "EOF";
 
  usage: rotate xcol ycol angle  *.*

 xcol,yxol=columns containing x and y , 
 angle=angle of rotation around z
 *.* .. filenname
 
 the rotation is done using the following formula
 x'= cos(angle)*x+sin(angle)*y
 y'=-sin(angle)*x+cos(angle)*y

EOF
 exit 0;}

 

$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;
$angle=eval{$ARGV[0]*3.1415926535897932384626433832795/180};shift @ARGV;



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

                  $x=$numbers[$colx-1];
                  $y=$numbers[$coly-1];
                  $numbers[$colx-1]= cos($angle)*$x+sin($angle)*$y;
                  $numbers[$coly-1]=-sin($angle)*$x+cos($angle)*$y;
		  foreach (@numbers)

		  {++$i;

		 

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