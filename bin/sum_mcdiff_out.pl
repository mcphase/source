#!/usr/bin/perl
use Getopt::Long;
#use Math::ElFun;
use PDL;
#use PDL::Slatec;
#use Switch;
use File::Copy;
$fileout="./results/mcdiffsummed.out";
@AARGV=@ARGV;
GetOptions("help"=>\$helpflag,
           "average"=>\$avflag);
usage() if $helpflag;
if ($avflag) {$text="averaged";}else{$text="summed";}
$n=0;
$numbersold[6]=0;
$numbersold[7]=0;
$numbersold[8]=0;
$numbersold[11]=0;
unless (open(Fin,"./results/mcdiff.out")){ die "Error program sum_mcdiff_out: file  ./results/mcdiff.out not found\n (type 'sum_mcdiff_out -h' to get help )\n";}
open(Fout,">".$fileout);
while($line=<Fin>)
{if ($line=~/^\s*#/){print Fout $line;}
 else
 { $line=~s/D/E/g; @numbers=split(" ",$line);
if($numbers[12]<0){++$n;
$numbersold[6]+=$numbers[6];
$numbersold[7]+=$numbers[7];
$numbersold[8]+=$numbers[8];
$numbersold[11]+=$numbers[11];
print Fout sprintf("# $text by sum_mcdiff_out.pl %6.3f %6.3f %6.3f \n",
                    $numbers[0],
                    $numbers[1],
                    $numbers[2]);
}
else
{$numbers[6]+=$numbersold[6];++$n;
$numbers[7]+=$numbersold[7];
$numbers[8]+=$numbersold[8];
$numbers[11]+=$numbersold[11];
unless ($avflag) {$n=1;}
print Fout sprintf("%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %5.4E %8.4f %8.4f %8.4f %5.4E %8.4f %8.4f\n",
                    $numbers[0],
                    $numbers[1],
                    $numbers[2],
                    $numbers[3],
                    $numbers[4],
                    $numbers[5],
                    $numbers[6]/$n,
                    $numbers[7]/$n,
                    $numbers[8]/$n,
                    $numbers[9],
                    $numbers[10],
                    $numbers[11]/$n,
                    $numbers[12],
                    $numbers[13]);
$numbersold[6]=0;$n=0;
$numbersold[7]=0;
$numbersold[8]=0;
$numbersold[11]=0;

}
 }
}

print STDERR << "EOF";
 all calculated reflections with Iobs=-1 have been $text

 and stored in output file $fileout

EOF
exit;
sub usage() {

  print STDERR << "EOF";

    sum_mcdiff_out: program to sum all calculated intensities (columns 7 8 9 12)
                    in mcdiff.out which have been summed from experiment

    usage: sum_mcdiff_out -h
           sum_mcdiff_out -av
           sum_mcdiff_out

     -h           : this (help) message
     -av          : average instead of summing all calculated intensities 

    output files:

    $fileout

EOF

  exit;

}
