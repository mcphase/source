#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >6)
{print STDOUT <<"EOF";

  adderr: program to adderr functions y1(x1) with y1err and y2(x2) with y2err
       taken from data file1 and data file2

  usage: adderr colx1 coly1 coly1err file1 colx2 coly2 coly2err file2

  input:
  file1, file2         filennames
  colx, coly, colyerr  columns containing x and y=f(x) and yerror values

  output:
  file1            contains in coly1=coly1+f2(colx1)
                   and in   coly1err=sqrt[coly1err^2+f2err(colx1)^2]
                   f2(colx1) is calculated by linear
                   interpolation
                   f2(colx1)=
                   (colx1-colx2(n)*(coly2(n+1)-coly2(n))/(colx2(n+1)-colx2(n))
                   f2err(colx1)=
                   (colx1-colx2(n)*(coly2err(n+1)-coly2err(n))/(colx2(n+1)-colx2(n))

  note:            colx2 has to be sorted in file2

EOF
 exit 0;}

 
$colx1=$ARGV[0];shift @ARGV;
$coly1=$ARGV[0];shift @ARGV;
$coly1err=$ARGV[0];shift @ARGV;
$file1=$ARGV[0];shift @ARGV;
$colx2=$ARGV[0];shift @ARGV;
$coly2=$ARGV[0];shift @ARGV;
$coly2err=$ARGV[0];shift @ARGV;
$file2=$ARGV[0];shift @ARGV;

#read file2
$nn=-1;
   unless (open (Fin2, $file2)){die "\n error:unable to open $file2\n";}
   while($line=<Fin2>)
     {if ($line=~/^\s*#/){}
      else{++$nn;($x[$nn],$y[$nn],$yerr[$nn],@numbers)=getxy($line,$colx2,$coly2,$coly2err);
          }
     }
   close Fin2;
if ($nn<1) {die "ERROR add: too few points in file $file2\n";}
$order=+1;if($x[1]-$x[0]<0){$order=-1;}

#do addition
   unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}   
$n=0;
   
open (Fout, ">range.out");
   while($line1=<Fin1>)
     {if ($line1=~/^\s*#/) {print Fout $line1;}
      else{ ($x1,$y1,$y1err,@numout)=getxy($line1,$colx1,$coly1,$coly1err);
# find $n
while($x1<=$order*$x[$n]&$n>0){--$n;}
while($x1>=$order*$x[$n+1]&$n<$nn-1){++$n;}
unless($n==0||$n==$nn-1) # do not extrapolate
           { # treat equal values correctly
              while($x[$n]==$x[$n+1]&$n<$nn-1){++$n;}
              while($x[$n]==$x[$n+1]&$n>0){--$n;}
             #print $n;
              # do addition using linear interpolation
              $numout[$coly1-1]=$y1+$y[$n]+($x1-$x[$n])*($y[$n+1]-$y[$n])/($x[$n+1]-$x[$n]);
                  $y2err=$yerr[$n]+($x1-$x[$n])*($yerr[$n+1]-$yerr[$n])/($x[$n+1]-$x[$n]);
              $numout[$coly1err-1]=sqrt($y1err*$y1err+$y2err*$y2err);

		 $i=0;
		   foreach (@numout)
		   {print Fout $numout[$i]." ";++$i;}
                    print Fout "\n";
                 
           }
          }
      }
      close Fin1;
  
       unless (rename "range.out",$file1)
     {unless(open (Fout, ">$file1"))
      {die "\n error:unable to write to $file1\n";}
      open (Fin, "range.out");
      while($line=<Fin>){ print Fout $line;}
      close Fin;
      close Fout;
      system "del range.out";
     }

sub getxy { my ($line,$colx,$coly,$colyerr)=@_;
           $line=~s/D/E/g;my @numbers=split(" ",$line);
           return ($numbers[$colx-1],$numbers[$coly-1],$numbers[$colyerr-1],@numbers);
           }

# **********************************************************************************************
# extracts variable from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# ($standarddeviation)=extract("sta","data.dat");
#
# ... it stores 0.24 in the variable $standarddeviation
#
sub extract {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";$value="";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);}}
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************

