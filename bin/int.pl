#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print STDOUT <<"EOF";


program int  used to integrate columnx vs columny=f(x) in data file

integration is done point by point, the result goes to the data file,
the total integral INT=integral f(x)dx is printed to stdout and set to the 
environment variable MCPHASE_INT.
if an error column is given, errors dy_i are summed as
 err=sqrt(sum_i dy_i^2 (x_i-x_i-1)^2)
and MCPHASE_INT_ERR=err is set.

usage: int [-m] colx coly[ecolyerr]  *.*

input:
  colx, coly, colyerr  columns containing x and y=f(x) and yerror values
  *.*                  filennames

option -m: n-th moments are calculated according to mu_1=integral x f(x)dx/INT
           mu_n=integral (x-mu_1)^n f(x)dx/INT

output:
the results go to stdout and environment variables MCPHASE_INT_MU_1,
MCPHASE_INT_MU_2,MCPHASE_INT_MU_3 ... are set

EOF
 exit 0;}

 

$colx=$ARGV[0];shift @ARGV;
$calcm=0;if($colx=~/-m/){$ARGV[0]=~s/x/*/g;$colx=eval $ARGV[0];shift @ARGV;$calcm=1;}
                  else {$colx=~s/x/*/g;$colx=eval $colx;}
$coly=$ARGV[0];shift @ARGV;
if ($coly=~/e/){$_=$coly;($colyerr)=/e(\d*)/;($coly)=/(\d*)e/;}else{$colyerr=0;}
$coly=~s/x/*/g;$coly=eval $coly;
$colyerr=~s/x/*/g;$colyerr=eval $colyerr;


  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

  $integral=0;$j=0;$mu[1]=0;$muerr2[1]=0;$err2=0;
    open (Fout, ">range.out");
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}
      else{$line=~s/D/E/g;@numbers=split(" ",$line);
           if ($j==0){@numbers1=@numbers;}
            ++$j;
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))
		  {@numout=@numbers;
                   $dx=($numbers[$colx-1]-$numbers1[$colx-1]);
                   $y=0.5*($numbers[$coly-1]+$numbers1[$coly-1]);
                   $integral+=$y*$dx;
                   $mu[1]+=$y*$dx*0.5*($numbers[$colx-1]+$numbers1[$colx-1]);
                   if($colyerr>0)
                   {$dy2=$numbers[$colyerr-1]*$numbers[$colyerr-1]/4+$numbers1[$colyerr-1]*$numbers1[$colyerr-1]/4;
                   $err2+=$dy2*$dx*$dx;
                   }
		  } 
	    @numbers1=@numbers;
            $numbers[$coly-1]=$integral;
            if($colyerr>0){$numbers[$colyerr-1]=sqrt($err2);}
            $i=0;foreach (@numbers){++$i;print Fout $numbers[$i-1]." ";}     
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
  if($colyerr>0){$err=sqrt($err2);print " INT_ERROR=".$err;}
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
if($colyerr>0){print Fout "set MCPHASE_INT_ERR=".$err."\n";}
for($i=1;$i<=$nmax;++$i){print Fout "set MCPHASE_INT_MU_$i=".$mu[$i]."\n";}
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_INT=$integral\n";
if($colyerr>0){print Fout "export MCPHASE_INT_ERR=".$err."\n";}
for($i=1;$i<=$nmax;++$i){print Fout "set MCPHASE_INT_MU_$i=".$mu[$i]."\n";}
close Fout;

#\end{verbatim} 