#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}



if ($#ARGV<1)

{print "program fermicol: calculates Fermi function column in file

    use as fermicol col T filename

          col     ...... column containing Energy values (eV) relative to EF
          T       ...... Temperatur in K
	  filename... filename
\n";
exit(0);
}

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$col=eval $ARGV[0]; shift @ARGV;--$col;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$T=eval $ARGV[0]; shift @ARGV;

 foreach (@ARGV)
  {$file=$_; open (Fout, ">range.out");
  if(open(Fin,$file))
  {      while($line=<Fin>)
     {if ($line=~/^\s*#/){print Fout $line;}
      else
      {$line=~s/D/E/g;@numbers=split(" ",$line);$i=0;
		  foreach (@numbers)
		  {
                   if ($i==$col) {# calculate Fermi function
                                    $x=$numbers[$i]/8.617e-5/$T;
                                    if ($x>5){$fd=0;}
                                    else {if ($x<-5){$fd=1;}
                                          else {$fd=1.0/(exp($x)+1.0);}
                                         }
                                   print Fout $fd." ";} # print out results
                    else {print Fout $numbers[$i]." ";}
                   ++$i;
 		  }
            print Fout "\n";++$i;
      }
     }
   close Fin;
  }
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

