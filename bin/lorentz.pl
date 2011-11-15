#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program used to calculate a lorentzian \n";

 print " usage: lorentz fwhm stp min max

the formula for a lorentz curve is:
lorentz(x)=1.0/3.1415/fwhm/(1.0+x^2/fwhm^2)
\n";

 exit 0;}

print "# $0 @ARGV\n";

$fwhm=$ARGV[0]/2;shift @ARGV;

$stp=$ARGV[0];shift @ARGV;

$min=$ARGV[0];shift @ARGV;

$max=$ARGV[0];shift @ARGV;



for($i=$min;$i<$max;$i+=$stp)

{

$lorentz=1.0/3.1415/$fwhm/(1.0+$i*$i/$fwhm/$fwhm);

print $i." ".$lorentz."\n";

}

