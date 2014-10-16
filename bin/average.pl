#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}
use Getopt::Long;
sub usage()
{print STDOUT << "EOF";

program average  used  to AVERAGE data by deleting close data points

use as:

 average [-h] [-help]
         [-middle]       middle point is taken (default)
         [-av]           points are averaged
         [-sum]           points are added
         [-sumcol=12,13,14] but points in column 12,13,14 are added
	 [-median]       median of points is calculated and kept
         [-dmin=0.4]    takes instead of 15 lines a variable number
                         of lines determined by the condition that
                         data in column 15 is closer than dmin 
         15 filenames    takes sets of 15 lines in data file and averages data

EOF
 exit 0;}

  usage() if ($#ARGV<1);

@p = @ARGV;
@pars = @ARGV;
if (join('',@pars) =~/\-/) {
  # see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions
  GetOptions("help"=>\$helpflag,
             "av"=>\$avflag,
             "median"=>\$medianflag,
             "sum"=>\$sumflag,
             "sumcol=s"=>\$sumcol,
             "dmin=s"=>\$dmin);
  usage() if $helpflag;
}
if($dmin){$dd=$dmin;$dd=~s/exp/essp/g;$dd=~s/x/*/g;$dd=~s/essp/exp/g;$dmin=eval $dd;}
@ARGV=@p;while (join('',@pars) =~ /\-/){shift @ARGV;@pars = @ARGV;}
# if ($dmin) {shift @ARGV;} # removed because it squeezes out the column MR 1.10.2011
$ARGV[0]=~s/exp/essssssssssssssp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essssssssssssssp/exp/g;$n=eval $ARGV[0];shift @ARGV;

if ($avflag) {print "averaging ...\n";}
elsif ($sumflag) {print "summing  ...\n";}
elsif ($medianflag) {print "taking median ...\n";}
else {print "taking middle point ...\n";}
if ($sumcol){print "...but summing column $sumcol ...\n";@scol=split(",",$sumcol);}
if ($dmin) {print "taking blocks of lines with data in column $n closer than $dmin\n";}
else {print "taking $n lines\n";}

--$n; # necessary because arrays indices start at 0 and not at 1
  foreach (@ARGV)
  {$file=$_;$ii=-1;
   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
   print "<".$file;
   open (Fout, ">range.out");
   while($line=<Fin>)
     {if ($line=~/^\s*#/) {print Fout $line;}
      else{$line=~s/D/E/g;@numbers=split(" ",$line);

                       if($dmin){
                                      
                                 if($ii==-1||(abs($field[0][$n]-$numbers[$n])<$dmin)){++$ii;
                                      for($k=0;$k<=$#numbers;++$k){$field[$ii][$k]=$numbers[$k];}
                                      # for($i=0;$i<=$ii;++$i){for($k=0;$k<=$#numbers;++$k) {print $field[$i][$k]." ";} print "\n";}
                                      # print "\nhello\n";
                                      }
                                 else {emptyblock();}
                                }
                            else{
                                    # here we put in the new data line into the field used to calculate the output
                                    if($ii<$n){++$ii;}
                                    for($k=0;$k<=$#numbers;++$k){$field[$ii][$k]=$numbers[$k];}
                                   
                                    emptyblock() if($ii==$n) ;
                                 }
           }
      }
      if($ii>-1){emptyblock();}
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

sub emptyblock()
{@numout=();
 for($k=0;$k<=$#numbers;++$k) {$avfl=$avflag;$sumfl=$sumflag;$medianfl=$medianflag;
                  if ($sumcol){foreach(@scol){if($_==$k+1){$sumfl=1;}}}
                  if ($avfl||$sumfl) {# here  average
                                                  for($i=0;$i<=$ii;++$i){$numout[$k]+=$field[$i][$k];}
                                                  if($avflag){$numout[$k]/=($ii+1);}
                                          }
                                    elsif ($medianfl) {
                                                  $i=$ii;
                                                  while($i>1){$max=-1e10;$min=1e10;
                                                          for($kk=0;$kk<=$ii;++$kk)
                                                                {if($field[$kk][$k]<$min&&$field[$kk][$k]>-1e10){$min=$field[$kk][$k];$kkmin=$kk;}
                                                                 if($field[$kk][$k]>$max&&$field[$kk][$k]<+1e10){$max=$field[$kk][$k];$kkmax=$kk;}
                                                                }
                                                              $field[$kkmin][$k]=-1e10;$field[$kkmax][$k]=+1e10;
                                                              $i-=2;
                                                             }
                                                   for($kk=0;$kk<=$ii;++$kk)
                                                          {if($field[$kk][$k]>-1e10&&$field[$kk][$k]<+1e10){$numout[$k]=$field[$kk][$k];}
                                                          }
                                                 }
                                    else         {# here take middle point
                                                  $i=$ii/2;
                                                 $numout[$k]=$field[$i][$k];
                                                 }
                                 }
                  $ii=0;$i=0;for($k=0;$k<=$#numbers;++$k){$field[$ii][$k]=$numbers[$k];} # put new line to field
		   foreach (@numout)
		   {print Fout $numout[$i]." ";++$i;}
                    print Fout "\n";
                 }