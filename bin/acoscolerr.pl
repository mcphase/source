#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



use Math::Trig;
unless ($#ARGV >0) 
{print " program acoscolerr  used to calculate arccosine of a column\n";
 print " usage: acoscolerr col colerr *.*   \n col=column\n colerr= corresponding error column\n *.* .. filenname\n";
 exit 0;}

$column=$ARGV[0];shift @ARGV;
$columnerr=$ARGV[0];shift @ARGV;

  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
   open (Fout, ">range.out");
   while($line=<Fin>)
     { if ($line=~/^\s*#/) {print Fout $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
           	  $i=0;++$j;
		  foreach (@numbers)
		  {++$i;
		  if ($i==$column) {$numbers[$i-1]=acos($numbers[$i-1]);}
		  if ($i==$columnerr) {$numbers[$i-1]=abs($numbers[$i-1]/sin(acos($numbers[$column-1])));}
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