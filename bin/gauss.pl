#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program used to calculate a gaussian \n";

 print " usage: gauss fwhm stp min max\n";

 exit 0;}

print "# $0 @ARGV\n";

$fwhm=$ARGV[0];shift @ARGV;

$stp=$ARGV[0];shift @ARGV;

$min=$ARGV[0];shift @ARGV;

$max=$ARGV[0];shift @ARGV;



for($i=$min;$i<$max;$i+=$stp)

{$sigma=$fwhm/2/sqrt(2*log(2));

$gauss=1.0/sqrt(2*3.1415)/$sigma* exp(-$i*$i/$sigma/$sigma/2);

print $i." ".$gauss."\n";

}

