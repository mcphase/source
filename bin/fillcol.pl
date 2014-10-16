#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}
use Math::Trig;
unless ($#ARGV >1) 
{print stdout << "EOF";
 program fillcol  used to fill column with numbers in data file
 usage: fillcol col expression *.*  

 col       ....   column
 expression ...   e.g. 'tan(c1x7.52)+c2'  
                  here c1,c2,.. refer to column 1,2,...
                  operations are multiplication (x), division(/)
                  addition (+), subtraction(-), power (xx)
                  trigonometric functions tan,cos,sin, asin,atan,acos
                  exp,(natural) log                  
 *.*       ....   filenname
 
EOF
 exit 0;}

$ARGV[0]=~s/x/*/g;$col=$ARGV[0];shift @ARGV;
$col=eval $col;

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$exp=$ARGV[0];shift @ARGV;
$exp =~s/c(\d+)/\$numbers[\1-1]/g; print $exp."\n";

  foreach (@ARGV)
  {$file=$_;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
   print "<".$file;
   open (Fout, ">range.out");
   while($line=<Fin>)
   {
      if ($line=~/^\s*#/) {print Fout $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);                       
                       $numbers[$col-1]=eval $exp;
           	  $i=0;foreach (@numbers)
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
      system "del range.out";
     }
   print ">\n";
   }



