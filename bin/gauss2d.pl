#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program used to calculate a gaussian \n";

 print " usage: gauss2d fwhm1 fwhm2 rotangle stpx minx maxx stpy miny maxy

 the formula for a gaussian is:
sigma=fwhm/2/sqrt(2*log(2))
gauss(x)=1.0/sqrt(2*3.14159265359)/sigma* exp(-x^2/2sigma^2)\n";



 exit 0;}

print "# $0 @ARGV\n";
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$fwhm1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$fwhm2=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$rotangle=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$stpx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$minx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$maxx=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$stpy=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$miny=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$maxy=eval $ARGV[0];shift @ARGV;

$rotangle*=3.141592654/180;

for($x=$minx;$x<$maxx;$x+=$stpx)
{
for($y=$miny;$y<$maxy;$y+=$stpy)
{
$sigma1=$fwhm1/2/sqrt(2*log(2));
$sigma2=$fwhm2/2/sqrt(2*log(2));

$u=$x*cos($rotangle)-$y*sin($rotangle);
$v=$x*sin($rotangle)+$y*cos($rotangle);

$gauss=1.0/sqrt(2*3.14159265359)/$sigma1* exp(-$u*$u/$sigma1/$sigma1/2);
$gauss*=1.0/sqrt(2*3.14159265359)/$sigma2* exp(-$v*$v/$sigma2/$sigma2/2);

print $x." ".$y." ".$gauss."\n";
} print "#\n";
}

