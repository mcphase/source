#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >4)
{print STDOUT <<"EOF";

  add: program to add functions y1(x1) (with optional y1err) and y2(x2) (with optional y2err)
       taken from data file1 and data file2

  usage: add colx1 coly1[ecoly1err] file1 colx2 coly2[ecoly2err] file2

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
 exit 0;}else{print"#* add *\n";}


$ARGV[0]=~s/x/*/g;$colx1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$coly1=$ARGV[0];shift @ARGV;
if ($coly1=~/e/){$_=$coly1;($coly1err)=/e(\d*)/;($coly1)=/(\d*)e/;}else{$coly1err=0;}
$coly1= eval $coly1; $coly1err= eval $coly1err;
$file1=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$colx2=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$coly2=$ARGV[0];shift @ARGV;
if ($coly2=~/e/){$_=$coly2;($coly2err)=/e(\d*)/;($coly2)=/(\d*)e/;}else{$coly2err=0;}
$coly2= eval $coly2; $coly2err= eval $coly2err;
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
#print $order."\n";
#do addition
   unless (open (Fin1, $file1)){die "\n error:unable to open $file1\n";}
$n=0;

open (Fout, ">range.out");
   while($line1=<Fin1>)
     {if ($line1=~/^\s*#/) {print Fout $line1;}
      else{ ($x1,$y1,$y1err,@numout)=getxy($line1,$colx1,$coly1,$coly1err);
# find $n
while($order*$x1<$order*$x[$n]&$n>=0){--$n;}
while($order*$x1>$order*$x[$n+1]&$n<=$nn-1){++$n;}
unless($n<0||$n>$nn-1) # do not extrapolate
           { # treat equal values correctly
              while($x[$n]==$x[$n+1]&$n<=$nn-1){++$n;}
              while($x[$n]==$x[$n+1]&$n>=0){--$n;}
           unless($n<0||$n>$nn-1) # do not extrapolate
             { #print $n;
              # do addition using linear interpolation
              $numout[$coly1-1]=$y1+$y[$n]+($x1-$x[$n])*($y[$n+1]-$y[$n])/($x[$n+1]-$x[$n]);
              if($coly1err>0&&$coly2err>0)
              {     $y2err=$yerr[$n]+($x1-$x[$n])*($yerr[$n+1]-$yerr[$n])/($x[$n+1]-$x[$n]);
              $numout[$coly1err-1]=sqrt($y1err*$y1err+$y2err*$y2err);
              }
		 $i=0;
		   foreach (@numout)
		   {print Fout $numout[$i]." ";++$i;}
                    print Fout "\n";

              }
             else
             {$n=0;}
           }
      else
         {$n=0;}
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
           if ($colyerr>0){$yerr=$numbers[$colyerr-1];}else{$yerr=0;}
           return ($numbers[$colx-1],$numbers[$coly-1],$yerr,@numbers);
           }

