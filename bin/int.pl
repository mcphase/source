#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program int  used to integrate columnx vs columny=f(x) in data file\n";

 print " integration is done point by point, the result goes to the data file,\n";
 print " the total integral INT=integral f(x)dx is printed to stdout and set to the \n";
 print " environment variable MCPHASE_INT\n";
 print " usage: int [-m] colx coly  *.*   \n colx=columnx, coly=columny \n *.* .. filenname\n";
 print "  option -m: n-th moments are calculated according to mu_1=integral x f(x)dx/INT\n";
 print "             mu_n=integral (x-mu_1)^n f(x)dx/INT\n";
 print "             the results go to stdout and environment variables MCPHASE_INT_MU_1,\n";
 print "             MCPHASE_INT_MU_2,MCPHASE_INT_MU_3 ... are set\n";
 exit 0;}

 

$colx=$ARGV[0];shift @ARGV;
$calcm=0;if($colx=~/-m/){$colx=$ARGV[0];shift @ARGV;$calcm=1;}
$coly=$ARGV[0];shift @ARGV;

  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

  $integral=0;$j=0;$mu[1]=0;
    open (Fout, ">range.out");
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}
      else{$line=~s/D/E/g;@numbers=split(" ",$line);
           if ($j==0){@numbers1=@numbers;}
            ++$j;
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))
		  {@numout=@numbers;
                   $d=0.5*($numbers[$coly-1]+$numbers1[$coly-1])*($numbers[$colx-1]-$numbers1[$colx-1]);
                   $integral+=$d;
                   $mu[1]+=$d*0.5*($numbers[$colx-1]+$numbers1[$colx-1]);
		  } 
	    @numbers1=@numbers;$i=0;
            foreach (@numbers)
		  {++$i;if ($i!=$coly){print Fout $numbers[$i-1]." ";}else{print Fout $integral." ";}
		  }     
            print Fout "\n";
           }
      }
      close Fin;
      close Fout;$nmax=1;
   if($calcm==1){ # here calculate mth moments 1- $nmax
                  $nmax=6;
                 unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
                 $center=$mu[1]/$integral;@mu= (0) x ($nmax+1);            # (0,0,0,0,0)
                 $j=0;
                  while($line=<Fin>)
                  {if ($line=~/^\s*#/) {}
                   else{$line=~s/D/E/g;@numbers=split(" ",$line);
                   if ($j==0){@numbers1=@numbers;}
                   ++$j;
                   unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))
  	   	      {@numout=@numbers;
                       $d=0.5*($numbers[$coly-1]+$numbers1[$coly-1])*($numbers[$colx-1]-$numbers1[$colx-1]);
                       for($k=1;$k<=$nmax;++$k)
                         {$d*=(0.5*($numbers[$colx-1]+$numbers1[$colx-1])-$center);$mu[$k]+=$d;}
		      }
	           @numbers1=@numbers;
                       }
                  }
                 close Fin;$mu[1]=$center*$integral;
                }
   $i=0; foreach(@mu){$mu[$i]/=$integral;++$i}
   print " x=col".$colx." y=f(x)=col".$coly." INT=".$integral;
   for($i=1;$i<=$nmax;++$i){print " mu_$i=".$mu[$i];}
   print ">\n";
       unless (rename "range.out",$file)
      {unless(open (Fout, ">$file"))     
      {die "\n error:unable to write to $file\n";}
      open (Fin, "range.out");
      while($line=<Fin>){ print Fout $line;}
      close Fin;
      close Fout;
      system "del range.out";  }


   }

# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_INT=$integral\n";
for($i=1;$i<=$nmax;++$i){print Fout "set MCPHASE_INT_MU_$i=".$mu[$i]."\n";}
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_INT=$integral\n";
for($i=1;$i<=$nmax;++$i){print Fout "set MCPHASE_INT_MU_$i=".$mu[$i]."\n";}
close Fout;

#\end{verbatim} 