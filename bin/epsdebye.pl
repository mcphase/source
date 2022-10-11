#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#\begin{verbatim}



unless ($#ARGV >2) 

{print " program epsdebye:

            Usage: epsdebye Tmax dT Tdebye scale Tnorm [d1 d2 datafile]

	        

		    calculates the strain epsilon using the debye model
		    according to the following formula:

		    epsilon=scale*T*D(Tdebye/T)+const

		    with

		    D(z)=3/z^3 integral_0^z [x^3/(e^x-1)]dx

                 Range is from zero to Tmax in stepwidths dT
		 unless a datafile is given. The const is chosen such that
		 the result epsilon is zero at T=Tnorm. 

                 If a  datafile is given, with data column d1 and d2,the strain
                 is calculated for T-values of data column d1 and epsilon
		  is compared to data in column d2 - a standard 
                 deviation sta is calculated as a sum of squared deviations.
                 As output the datafile is given, an additional is column added 
		 containing the calculated strain epsilon. The datafile has to
		 be sorted according to descending T values !!!

                 output is written to stdout.\n";

 exit 0;}else{print STDERR "#* $0 *\n";}

print "# $0 @ARGV\n";
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Tmax=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$dT=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Tdebye=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$scale=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Tnorm=eval $ARGV[0];shift @ARGV;


$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$d1=eval $ARGV[0];shift @ARGV;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$d2=eval $ARGV[0];shift @ARGV;

$file1=$ARGV[0];shift @ARGV;



      $dx=0.0001; # dx for approximating debye integral



      $int=0;

      $zold=$dx;

#calculate normalization constant

      $z=$Tdebye/$Tnorm;



      for($x=$zold;$x<$z;$x+=$dx)

        {$int+=$x*$x*$x*$dx/(exp($x)-1);

        }	     

	$zold=$z;

	

      $D=3*$int/$z/$z/$z;

      $const=-$scale*$Tnorm*$D;



unless ($file1)

{

      $int=0;

      $zold=$dx;

# calculate strain

for($T=$Tmax;$T>$dT;$T-=$dT)   

    { $z=$Tdebye/$T;

      

      for($x=$zold;$x<$z;$x+=$dx)

        {$int+=$x*$x*$x*$dx/(exp($x)-1);

        }	     

	$zold=$z;

	

      $D=3*$int/$z/$z/$z;

      $eps=$scale*$T*$D+$const;

    print  $T." ".$eps."\n";

   }

}

else

{#here do something if the "data"-file1 is given ....

 $sta=0;$iii=0;

      $int=0;

      $zold=$dx;



 open (Fin1, $file1);

 while($line1=<Fin1>)

 {

  if ($line1=~/^\s*#/) {print $line1;}

  else

  {# get x-values from file 3 

             $line1=~s/D/E/g;@numbers1=split(" ",$line1);

	     $T=$numbers1[$d1-1];

	     $yd=$numbers1[$d2-1];

	

# calculate strain

     $z=$Tdebye/$T;

     if ($z<$zold) {print "error file $file1: temperature in column $d1 not sorted descendingly";  exit(1);}

      

      for($x=$zold;$x<$z;$x+=$dx)

        {$int+=$x*$x*$x*$dx/(exp($x)-1);

        }	     

	$zold=$z;

	

      $D=3*$int/$z/$z/$z;

      $eps=$scale*$T*$D+$const;

    

               	  $i=0;

		  foreach (@numbers1)

		  {++$i;print $numbers1[$i-1]." ";}     

    print  $eps."\n";		     

    $sta+=($eps-$yd)*($eps-$yd);

   }

  } 

 close Fin1;

 print "#! sta=$sta\n";

}



#\end{verbatim} 