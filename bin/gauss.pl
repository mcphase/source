#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program used to calculate a gaussian \n";

 print " usage: gauss fwhm stp min max

 the formula for a gaussian is:
sigma=fwhm/2/sqrt(2*log(2))
gauss(x)=1.0/sqrt(2*3.14159265359)/sigma* exp(-x^2/2sigma^2)\n";



 exit 0;}

print "# $0 @ARGV\n";
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$fwhm=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$stp=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$min=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$max=eval $ARGV[0];shift @ARGV;



for($i=$min;$i<$max;$i+=$stp)

{$sigma=$fwhm/2/sqrt(2*log(2));

$gauss=1.0/sqrt(2*3.14159265359)/$sigma* exp(-$i*$i/$sigma/$sigma/2);

print $i." ".$gauss."\n";

}

