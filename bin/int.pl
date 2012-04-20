#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program int  used to integrate columnx vs columny in data file\n";

 print " integration is done point by point, result goes to stdout and \n";
 print " environment variable MCPHASE_INT\n";
 print " usage: int colx coly  *.*   \n colx=columnx, coly=columny \n *.* .. filenname\n";

 exit 0;}

 

$colx=$ARGV[0];shift @ARGV;

$coly=$ARGV[0];shift @ARGV;



  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

  $integral=0;$j=0;
    open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {print Fout $line;}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

            if ($j==0){@numbers1=@numbers;}
            ++$j;
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))

		  {@numout=@numbers;

                   $integral+=($numbers[$coly-1]+$numbers1[$coly-1])/2*($numbers[$colx-1]-$numbers1[$colx-1]);

		  } 

	    @numbers1=@numbers;
            $i=0;
            foreach (@numbers)
		  {++$i;
		  if ($i!=$coly){print Fout $numbers[$i-1]." ";}else{print Fout $integral." ";}
		  }     
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
 
   print " x=col".$colx." y=f(x)=col".$coly." integral=".$integral.">\n";

   }

# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_INT=$integral\n";
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_INT=$integral\n";
close Fout;

#\end{verbatim} 