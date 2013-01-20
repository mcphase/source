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
$ARGV[0]=~s/x/*/g;$fwhm=eval {$ARGV[0]/2};shift @ARGV;
$ARGV[0]=~s/x/*/g;$stp=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$min=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/x/*/g;$max=eval $ARGV[0];shift @ARGV;



for($i=$min;$i<$max;$i+=$stp)

{

$lorentz=1.0/3.1415/$fwhm/(1.0+$i*$i/$fwhm/$fwhm);

print $i." ".$lorentz."\n";

}

