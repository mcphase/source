#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use Getopt::Long;

unless ($#ARGV >3) 
{print STDERR " program posls (position least squares computation):

            Usage: posls [options] cp datafile cpos cI calcfile
	        
                computes a measure of goodness of a list of peaks
                given  in calcfile (in column cpos with intensities in column cI)
                by comparing it to experimental data in datafile
                (given as a list of peak positions in column cp)

                The following quantities are computed:

1 sta                               = 1/N sum_i |weight(i)|*[Pexp(i) - nearestPcalc]^[2*sign(weight(i))]
2 sta_without_antipeaks             = 1/N sum_i_with_weight(i)>0  weight(i)*[Pexp(i) - nearestPcalc]^2
3 sta_without_weights               = 1/N sum_i [Pexp(i) - nearestPcalc]^[2*sign(weight(i))]
4 sta_without_antipeaks_weights     = 1/N sum_i_with_weight(i)>0 [Pexp(i) - nearestPcalc]^2

Pexp(i) ........ an experimental peak position taken from column cp in datafile
weight(i) ...... corresponding weight, by default set to 1, or set by option -w 
nearestPcalc ... indicates the calculated energy which is next to Pexp(i)
N............... number of contributions to the sum

Options:

-w cw           read weight from column cw, a negative value indicates an antipeak - i.e.
                a position where no peak is observed 
-i 0.001        take into account calculated energy only if it is larger than Imin=0.001 (default Imin=0)
-s 1              output to stdout sta=sta (default)
-s 2              output to stdout sta=sta_without_antipeaks
-s 3              output to stdout sta=sta_without_weights
-s 4              output to stdout sta=sta_without_antipeaks_weights
-r              reverse sums, i.e. sum over calculated peaks meaning calculate
             1 sta                               = 1/N sum_i |weight(nearestexp)|*[Pcalc(i) - nearestPexp]^[2*sign(weight(nearestexp))]
             2 sta_without_antipeaks             = 1/N sum_i weight(nearestexp)*[Pcalc(i) - nearestPexp_with_weight>0]^2
             3 sta_without_weights               = 1/N sum_i [Pcalc(i) - nearestPexp]^[2*sign(weight(i))]
             4 sta_without_antipeaks_weights     = 1/N sum_i [Pcalc(i) - nearestPexp_with_weight>0]^2

	          \n";
 exit 0;}else{print STDERR "#* $0 *\n";}

$Imin=0;
# Parses command line options
die "exiting" unless GetOptions("w=s"=>\$cw,
                                "i=s"=>\$Imin,
                                "s=s"=>\$s,
                                "r"=>\$r
                               );

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cp=eval $ARGV[0];shift @ARGV;
$datafile=$ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cpos=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$cI=eval $ARGV[0];shift @ARGV;
$calcfile=$ARGV[0];shift @ARGV;

$ii=0;$i=1;
      unless (open (Fin1, $calcfile)){die "\n error:unable to open $calcfile\n";}   
while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);
	if($#numbers1>=$cp){
      #store calc pos int values
     if($numbers1[$cI-1]>$Imin)
	{++$ii;
	$PC[$ii]=$numbers1[$cpos-1];
	$IC[$ii]=$numbers1[$cI-1];	
        }
	                   }    
           }

     }


$ii1=0;
      unless (open (Fin1, $datafile)){die "\n error:unable to open $datafile\n";}   

while($line1=<Fin1>)
     {
       if ($line1=~/^\s*#/) {;}
       else{$line1=~s/D/E/g;@numbers1=split(" ",$line1);

	#store  data values
	if($#numbers1>=$cp){
	if(defined $cw){unless(($s==2||$s==4)&&$numbers1[$cw-1]<0)
                         {# only keep exp peak if needed
                            $PE[$ii1]=$numbers1[$cp-1];
                            $w[$ii1]=$numbers1[$cw-1];++$ii1;
                         } 
                       }	
        else
         {$PE[$ii1]=$numbers1[$cp-1];$w[$ii1]=1;++$ii1;}
	
	}
            }
     }
   close Fin1;     

$sta=0;$ii1=0;$N=0;
unless($r)
{foreach(@PE)
 {# search nearestPcalc
  $dmin=1e100;$pe=$_;
  if($s==3||$s==4){if($w[$ii1]>0){$w[$ii1]=1;}else{$w[$ii1]=-1;}}
   foreach(@PC){$d=abs($pe-$_);if($d<$dmin){$dmin=$d;}}
   $ss=$dmin*$dmin;
 #1 sta  = sum_i |weight(i)|*[Pexp(i) - nearestPcalc]^[2*sign(weight(i))]
 #2 sta_without_antipeaks             = sum_i_with_weight(i)>0  weight(i)*[Pexp(i) - nearestPcalc]^2
 #3 sta_without_weights               = sum_i [Pexp(i) - nearestPcalc]^[2*sign(weight(i))]
 #4 sta_without_antipeaks_weights     = sum_i_with_weight(i)>0 [Pexp(i) - nearestPcalc]^2
  if($w[$ii1]>0){ $sta+=$w[$ii1]*$ss;++$N;}
  if($w[$ii1]<0){ $sta-=$w[$ii1]/($ss+0.0001);++$N;}
 ++$ii1;
 }
}
else
{foreach(@PE){if($s==3||$s==4){if($w[$ii1]>0){$w[$ii1]=1;}else{$w[$ii1]=-1;}}++$ii1;}
 foreach(@PC)
 {# search nearestPexp
  $dmin=1e100;$pc=$_;
    $ii=0;foreach(@PE){$d=abs($pc-$_);if($d<$dmin){$dmin=$d;$ii1=$ii;}++$ii;}
   $ss=$dmin*$dmin;
  if($w[$ii1]>0){ $sta+=$w[$ii1]*$ss;++$N;}
  if($w[$ii1]<0){ $sta-=$w[$ii1]/($ss+0.0001);++$N;}
  
 ++$ii1;
 }
}
if($N>0){$sta/=$N;}
print "N=$N\nsta=$sta\n";
