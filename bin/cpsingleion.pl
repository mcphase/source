#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}


unless($#ARGV>=3)

  {print "Program cpsingleion - calculates specific heat c=du/dT from output file *.levels.cef of singleion program \n";

   print "             use as :  cpsingleion tmin tmax deltat file.levels.cef [-option]\n";

   print "             alternatively \n";

   print "             use as: cpso1ion col1 col2 datafile file.levels.cef [-option]\n";

   print "                     (take cp-data from datafile and calculate\n";

   print "                      standard deviation)\n";

   print "     output is written to stdout, energies have to be given in meV,\n";

   print "     temperatures tmin tmax deltat in Kelvin\n";

   print "     Options: -s   .... calculate entropy  s=integral c dT/T   (J/molK) instead of cp\n";

   print "              -f   .... calculate free energy f=-kT ln(z) (J/mol) instead of cp\n";

   print "              -u   .... calculate magnetic energy u= sum_i Ei exp(-Ei/kT)/z (J/mol) instead of cp\n";

   print "              -z   .... calculate partition sum z= sum_i exp(-Ei/kT) instead of cp\n";

   exit(1);

  }else{print STDERR "#* $0 *\n";}

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Tmin=eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;$Tmax=eval $ARGV[1];
$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;$deltaT=eval $ARGV[2];
$filename=$ARGV[3];

      $cptext="c=du/dT(J/molK)";

      if($ARGV[4]=~/-f/){$cptext="f=-kT ln(z)(J/mol)";}

       if($ARGV[4]=~/-z/){$cptext="z=sum_i exp(-Ei/kT)";}

       if($ARGV[4]=~/-u/){$cptext="u=sum_i Ei exp(-Ei/kT)/z (J/mol)";}

       if($ARGV[4]=~/-s/){$cptext="s=integral c dT/T(J/molK)";}



unless (open (Fin,$filename)) {unless (open (Fin,"so1ion.out")){print "ERROR cpso1ion: file $filename not found\n";exit(1);}else{print "#reading so1ion.out\n"; }}

else {print "#reading $filename\n";}



$noflevels=1; # initialize noflevels

# read energies from so1ion.out

  while($line=<Fin>)

  {if($line=~/^.*\QEigenvalues\E/){
       ($n)=($line=~m/^(?:#!|[^#])*\bEigenvalues\s*=\s*(.*)/);
       @E=split(" ",$n);
       $noflevels=$#E;
        $energyshift=$E[0];
       # ($noflevels)=extract("noflevels",$line);
      # ($energyshift)=extract("Eshift",$line);#($line=~m|\QEnergy shift\E\s*\Q(Eshift)\E\s*:\s*([\d.eEdD\Q-\E\Q+\E]+)|);
      # ($E[$i])=($line=~m|\QE( \E$i\Q)\E\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);
      # ($deg[$i])=extract("Degeneracy",$line);#($line=~m|.*\QDegeneracy\E\s*:\s*([\d.eEdD]+)|);
      for($i=0;$i<=$noflevels;++$i){$E[$i]-=$energyshift; $deg[$i]=1;
      print "#! E($i)=".$E[$i]."meV degeneracy ".$deg[$i]."-fold\n"; 
         }
  }     
   }
  close Fin;



$dT=0.001; # this is fixed value for calculation of derivative of energy



if (open (Fin,$deltaT))

{

$colT=$Tmin;

$colcp=$Tmax;

$ii=0;$sta=0;

# read temperatures and calculate cp + comp to experiment

  print  "# T(K)  $cptext  exp-$cpctext\n";

  while($line=<Fin>)

  {

       if ($line=~/^\s*#/) {}

       else{$line=~s/D/E/g;@numbers=split(" ",$line);

           	$T0=$numbers[$colT-1];

		$cpexp=$numbers[$colcp-1];

            ($cpclc)=cp(); # CALCULATE cp !!!!

	 ++$ii;

	 $sta+=($cpclc-$cpexp)*($cpclc-$cpexp);

         print   "$T $cpclc $cpexp\n";

        }

  }

  $sta/=$ii;

  print  "#sta=$sta\n";       

  close Fin;

}

else

{print  "# T(K)  $cptext\n";

 for($T0=$Tmin;$T0<$Tmax;$T0+=$deltaT)

 {($cpclc)=cp();

  print "$T0  $cpclc\n";

 }

}

exit(0);



sub cp {

$T=$T0+$dT/2;

            $Zp=0;$Up=0;

            for ($i=0;$i<=$noflevels;++$i)

                {$x=$E[$i]/$T/0.0862;

                 $Zp+=$deg[$i]*exp(-$x);

                 $Up+=$deg[$i]*$E[$i]*exp(-$x);

                 }

            $Up/=$Zp*0.0862; # magnetic energy per ion in K

            $Up+=$energyshift/0.0862; # shift energy 

            $Up*=1.38066*6.023; # magnetic energy in J per mol

	    

	    $T=$T0-$dT/2;
            $Zm=0;$Um=0;

            for ($i=0;$i<=$noflevels;++$i)

                {$x=$E[$i]/$T/0.0862;

                 $Zm+=$deg[$i]*exp(-$x);

                 $Um+=$deg[$i]*$E[$i]*exp(-$x);

                 }

            $Um/=$Zm*0.0862; # magnetic energy per ion in K

            $Um+=$energyshift/0.0862; # shift energy 

            $Um*=1.38066*6.023; # 0.0862meV/K*1.602e-22J/meV*6.023e23ion/mol 

                                # magnetic energy in J per mol

         

	 $cpc=($Up-$Um)/$dT; # specific heat 	    

       if($ARGV[4]=~/-f/){$cpc=($energyshift/0.0862-$T0*log(0.5*($Zm+$Zp)))*1.38066*6.023;} # helmholtz function

       if($ARGV[4]=~/-z/){$cpc=0.5*($Zm+$Zp)*exp(-$energyshift/$T0/0.0862);} # partition sum

       if($ARGV[4]=~/-u/){$cpc=0.5*($Um+$Up);} # magnetic energy

       if($ARGV[4]=~/-s/){$cpc=0.5*($Um+$Up)/$T0-($energyshift/$T0/0.0862-log(0.5*($Zm+$Zp)))*1.38066*6.023;} # entropy



return $cpc;

}

# **********************************************************************************************
# extracts number from string
#
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
#
# ... it stores 0.3 in the variable $standarddeviation
#
sub extract {
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************

# **********************************************************************************************
# extracts string from string
#
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
#
# ... it stores 0.3 in the variable $standarddeviation
#
sub extractstring {
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*\b([^\n\s]+)[\s\n]/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************


# **********************************************************************************************
# extracts number from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# ($standarddeviation)=extractfromfile("sta","data.dat");
#
# ... it stores 0.24 in the variable $standarddeviation
#
sub extractfromfile {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";$value="";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);}}
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
