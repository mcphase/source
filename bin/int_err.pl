#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program int_err  used to integrate columnx vs columny in data file\n";
 print " integration is done point by point - taking into account summation of error\n";
 print " errors dy_i are summed as: err=sqrt(sum_i dy_i^2 (x_i-x_i-1)^2)\n";
 print " usage: int colx coly colerr *.*   \n colx=columnx, coly=columny, colerr= error column \n *.* .. filenname\n";
 print " results is output to stdout and environment variables MCPHASE_INT and MCPHASE_ERR\n";
exit 0;}

 

$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;
$colerr=$ARGV[0];shift @ARGV;


  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;
  $integral=0;$j=0;$err2=0;
    open (Fout, ">range.out");
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
            if ($j==0){@numbers1=@numbers;}
            ++$j;
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))
		  {@numout=@numbers;
                   $dx=($numbers[$colx-1]-$numbers1[$colx-1]);
                   $y=($numbers[$coly-1]+$numbers1[$coly-1])/2;
                   $dy2=$numbers[$colerr-1]*$numbers[$colerr-1]/4+$numbers1[$colerr-1]*$numbers1[$colerr-1]/4;
                   $integral+=$y*$dx;
                   $err2+=$dy2*$dx*$dx;

		  } 
	    @numbers1=@numbers;
            $numbers[$coly-1]=$integral;
            $numbers[$colerr-1]=sqrt($err2);
            $i=0;
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
      system "del range.out";  }
   $err=sqrt($err2);
   print " x=col".$colx." y=f(x)=col".$coly." integral=".$integral." $err=$err>\n";
   }

# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_INT=$integral\n";
print Fout "set MCPHASE_ERR=$err\n";
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_INT=$integral\n";
print Fout "export MCPHASE_ERR=$err\n";
close Fout;

#\end{verbatim} 