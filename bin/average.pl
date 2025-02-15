#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}
use Getopt::Long;
sub usage()
{print STDOUT << "EOF";

program average  used  to reduce data by deleting close data points

use as:

 average [option] 15 filenames    takes sets of 15 lines in data file and averages data

 options:    
         [-h] [-help]
         [-middle]       middle point is taken (default)
         [-first]        first point is taken 
         [-last]         last point is taken 
         [-av]           points are averaged
         [-sum]          points are added
         [-sumcol=12,13,14] for column 12,13,14 datapoints are added up and the sum is output
	 [-median]       median of points is calculated and kept
         [-dmin=0.4]     takes instead of 15 lines a variable number
                         of lines determined by the condition that
                         data in column 15 is closer than dmin 
         
EOF
 exit 0;}

  usage() if ($#ARGV<1);

print"#* average *";

@p = @ARGV;
@pars = @ARGV;
if (join('',@pars) =~/\-/) {
  # see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions
  GetOptions("help"=>\$helpflag,
             "av"=>\$avflag,
             "median"=>\$medianflag,
             "last"=>\$lastflag,
             "first"=>\$firstflag,
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
elsif ($lastflag) {print "taking last point ...\n";}
elsif ($firstflag) {print "taking first point ...\n";}
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
              # store line 
                       if($dmin){
                                    # check if block needs to be emptied ( printed  ) before new line can be stored
                                 if($ii>-1&&(abs($field[0][$n+1]-$numbers[$n])>$dmin)){emptyblock();}
              ++$ii;$field[$ii][0]=$#numbers;
              for($k=0;$k<=$#numbers;++$k){$field[$ii][$k+1]=$numbers[$k];}
                                }
                            else{
              ++$ii;$field[$ii][0]=$#numbers;
              for($k=0;$k<=$#numbers;++$k){$field[$ii][$k+1]=$numbers[$k];}
                                    # here we put in the new data line into the field used to calculate the output
                                    if($ii>=$n){emptyblock();}
                                 } 
                                   
           }
      }
      if($ii>-1){ print "hh $ii\n";emptyblock();}
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
 for($k=0;$k<=$field[0][0];++$k) {$avfl=$avflag;$sumfl=$sumflag;$medianfl=$medianflag;$firstfl=$firstflag;$lastfl=$lastflag;
                  if ($sumcol){foreach(@scol){if($_==$k+1){$sumfl=1;}}}
                  if ($avfl||$sumfl) {# here  average
                                                  for($i=0;$i<=$ii;++$i){$numout[$k]+=$field[$i][$k+1];}
                                                  if($avflag){$numout[$k]/=($ii+1);}
                                          }
                                    elsif ($medianfl) {
                                                  $i=$ii;
                                                  while($i>1){$max=-1e10;$min=1e10;
                                                          for($kk=0;$kk<=$ii;++$kk)
                                                                {if($field[$kk][$k+1]<$min&&$field[$kk][$k+1]>-1e10){$min=$field[$kk][$k+1];$kkmin=$kk;}
                                                                 if($field[$kk][$k+1]>$max&&$field[$kk][$k+1]<+1e10){$max=$field[$kk][$k+1];$kkmax=$kk;}
                                                                }
                                                              $field[$kkmin][$k+1]=-1e10;$field[$kkmax][$k+1]=+1e10;
                                                              $i-=2;
                                                             }
                                                   for($kk=0;$kk<=$ii;++$kk)
                                                          {if($field[$kk][$k+1]>-1e10&&$field[$kk][$k+1]<+1e10){$numout[$k]=$field[$kk][$k+1];}
                                                          }
                                                 }
                                    elsif ($firstfl)        {# here take first point
                                                  $i=0;
                                                 $numout[$k]=$field[$i][$k+1];
                                                 }
                                    elsif ($lastfl)        {# here take last point
                                                  $i=$ii;
                                                 $numout[$k]=$field[$i][$k+1];
                                                 }
                                    else         {# here take middle point
                                                  $i=int($ii/2); 
                                                 $numout[$k]=$field[$i][$k+1];
                                                 }
                                 }
                  $ii=-1;$i=0;foreach (@numout)
		   {print Fout $numout[$i]." ";++$i;}
                    print Fout "\n";
}
