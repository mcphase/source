#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1)

{print " program extractmatrix  used to extract a matrix block from a data file\n";
 print " usage: extractmatrix colmin colmax linemin linemax  *.*   \n \
         colmin, colmax, linemin, linemax = matrix columns and lines \n \
          *.* .. filenname\n";
 exit 0;}else{print STDERR "#* $0 *\n";}

$ARGV[0]=~s/x/*/g;$colmin=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$colmax=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$linemin=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$linemax=eval $ARGV[0];shift @ARGV;

  foreach (@ARGV)
  { $file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
    $j=0;$firstline=0;
   while($line=<Fin>)
     {
       if ($line=~/^\s*#/) {print $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
           	  $i=0;++$j;
                  if($linemin<=$j&&$j<=$linemax){
                  if($firstline==0){$firstline=1; 
print "# extracted matrix column $colmin to $colmax and lines $linemin to $linemax from file ".$file."\n";
                                    }
		  foreach (@numbers)
		  {++$i;
		  if ($colmin<=$i&&$i<=$colmax){print $numbers[$i-1]." ";}

		  }     

            print  "\n";}

           }

      }

      close Fin;

  }

#\end{verbatim} 
