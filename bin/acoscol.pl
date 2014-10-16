#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use Math::Trig;
unless ($#ARGV >0)
{print " program acoscol  used to calculate arccosine of a column\n";
 print " usage: acoscol col[ecolerr] *.*   \n col=column\n optional colerr= corresponding error column, e.g. 3e4 means argumentvalues are in column 3 and error in column 4\n *.* .. filenname\n";
 exit 0;}
$ARGV[0]=~s/x/*/g;
$column=eval $ARGV[0];shift @ARGV;
if ($column=~/e/){$_=$column;($columnerr)=/e(\d*)/;($column)=/(\d*)e/;}else{$columnerr=0;}
#print $column." ".$columnerr."\n";exit;

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
