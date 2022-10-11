#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program dif  used to differentiate d(columny)/d(columnx) in data file, the result goes to coly\n";

 print " differentiation is done with averaging n points\n";

 print " usage: dif colx coly  n *.*   \n colx=columnx, coly=columny, n=number of points to average \n *.* .. filenname\n";

 exit 0;}else{print STDERR "#* $0 *";}

 

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$colx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$coly=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$n=eval $ARGV[0];shift @ARGV;


  foreach (@ARGV)
  {$file=$_;$ii=-1;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
   print "<".$file;
   open (Fout, ">range.out");
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}
      else{$line=~s/D/E/g;@numbers=split(" ",$line);
                                     # here we put in the new data line into the field used to calculate the output
                                     if($ii<$n){++$ii;}else{shift @field;}
                                     for($k=0;$k<=$#numbers;++$k){$field[$ii][$k]=$numbers[$k];}
                                    #for($i=0;$i<=$ii;++$i){for($k=0;$k<=$#numbers;++$k) {print $field[$i][$k]." ";} print "\n";}
                                    #print "\nhello\n";
                                    # here calculate the derivative !!! and average the rest
                                    @numout=@numbers;
                                    for($i=0;$i<$ii;++$i){for($k=0;$k<=$#numbers;++$k) {$numout[$k]+=$field[$i][$k];}}
                                    for($k=0;$k<=$#numbers;++$k) {$numout[$k]/=$ii+1;}
                                    $xm=$numout[$colx-1];
                                    $ym=$numout[$coly-1];

                                    $xxm = 0;$xym = 0; 
                                    for($i=0;$i<=$ii;++$i){$xym+=$field[$i][$colx-1]*$field[$i][$coly-1];
                                                           $xxm+=$field[$i][$colx-1]*$field[$i][$colx-1];} 
                                    $xxm /= $ii+1; $xym /= $ii+1;

              if(($xxm - $xm * $xm)!=0&&$ii==$n)
                 {             #  here calculate the slope of the column jj% with respect to column ii%
                                   $numout[$coly-1]=($xym - $ym * $xm) / ($xxm - $xm * $xm);
                   $i=0;
		   foreach (@numout)
		   {print Fout $numout[$i]." ";++$i;}     
                    print Fout "\n";
                 } 
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