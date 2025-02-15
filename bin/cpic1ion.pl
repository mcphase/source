#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}



unless($#ARGV>=2)

  {print "Program cpic1ion - calculates specific heat from output file results/ic1ion.out\n";

   print "             use as :  cpic1ion tmin tmax deltat [-option]\n"; 

   print "             alternatively \n";

   print "             use as: cpic1ion col1 col2 datafile [-option]\n";

   print "                     (take cp-data from datafile and calculate\n";

   print "                      standard deviation)\n";

   print "     output is written to stdout, energies have to be given in meV,\n";

   print "     temperatures tmin tmax deltat in Kelvin\n";

   print "     Options: -s   .... calculate entropy  (J/molK) instead of cp\n";

   print "              -f   .... calculate free energy (J/mol) instead of cp\n";

   print "              -u   .... calculate magnetic energy (J/mol) instead of cp\n";

   print "              -z   .... calculate partition sum instead of cp\n";

   exit(1);

  }else{print STDERR "#* $0 *\n";}

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$Tmin=eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;$Tmax=eval $ARGV[1];
$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;$deltaT=eval $ARGV[2];

      $cptext="cp(J/molK)";

      if($ARGV[3]=~/-f/){$cptext="f(J/mol)";}

       if($ARGV[3]=~/-z/){$cptext="z";}

       if($ARGV[3]=~/-u/){$cptext="u(J/mol)";}

       if($ARGV[3]=~/-s/){$cptext="s(J/molK)";}



unless (open (Fin,"results/ic1ion.out")) {unless (open (Fin,"ic1ion.out")){print "ERROR cpic1ion: file results/ic1ion.out not found\n";exit(1);}else{print "#reading ic1ion.out\n"; }}

else {print "#reading results/ic1ion.out\n";}



$noflevels=0; # initialize noflevels

# read energies from ic1ion.out

  while($line=<Fin>)

  {if($line=~/^.*#/)
   {if($line=~/^.*\QEnergy(\E/){unless ($line=~/^.*\QEnergy(meV)\E/)
                                      {print "ERROR cpic1ion: energies in ic1ion.out must be in meV ! - ".$line."\n";exit(1);}
                               }
    if($line=~/^.*\QEnergy offset\E.*E0/) # read energyshift
      {($energyshift)=($line=~m|\QE0\E\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);
       print "#Energy shift: $energyshift meV\n";
      } 
      else
      {print $line;}
   }
   else
   {if($line=~/^[\d.eEdD\Q-\E\Q+\E]/)
    {++$noflevels;
    @numbers=split(" ",$line);
    $E[$noflevels]=$numbers[0];$deg[$noflevels]=1;
    print "#! E($noflevels)=".$E[$noflevels]."meV degeneracy ".$deg[$noflevels]."-fold\n";
    }
   }
  }

  close Fin;



$dT=0.1; # this is fixed value for calculation of derivative of energy



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

            for ($i=1;$i<=$noflevels;++$i)

                {$x=$E[$i]/$T/0.0862;

                 $Zp+=$deg[$i]*exp(-$x);

                 $Up+=$deg[$i]*$E[$i]*exp(-$x);

                 }

            $Up/=$Zp*0.0862; # magnetic energy per ion in K

            $Up+=$energyshift/0.0862; # shift energy 

            $Up*=1.38066*6.023; # magnetic energy in J per mol

	    

	    $T=$T0-$dT/2;

            $Zm=0;$Um=0;

            for ($i=1;$i<=$noflevels;++$i)

                {$x=$E[$i]/$T/0.0862;

                 $Zm+=$deg[$i]*exp(-$x);

                 $Um+=$deg[$i]*$E[$i]*exp(-$x);

                 }

            $Um/=$Zm*0.0862; # magnetic energy per ion in K

            $Um+=$energyshift/0.0862; # shift energy 

            $Um*=1.38066*6.023; # 0.0862meV/K*1.602e-22J/meV*6.023e23ion/mol 

                                # magnetic energy in J per mol

         

	 $cpc=($Up-$Um)/$dT; # specific heat 	    

       if($ARGV[3]=~/-f/){$cpc=($energyshift/0.0862-$T0*log(0.5*($Zm+$Zp)))*1.38066*6.023;} # helmholtz function

       if($ARGV[3]=~/-z/){$cpc=0.5*($Zm+$Zp)*exp(-$energyshift/$T0/0.0862);} # partition sum

       if($ARGV[3]=~/-u/){$cpc=0.5*($Um+$Up);} # magnetic energy

       if($ARGV[3]=~/-s/){$cpc=0.5*($Um+$Up)/$T0-($energyshift/$T0/0.0862-log(0.5*($Zm+$Zp)))*1.38066*6.023;} # entropy



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
