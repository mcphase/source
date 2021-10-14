#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}
use PDL;
# use PDL::Slatec;
use File::Copy;

# copy parameter files to specific locations
# remember a few useful perl commands 
# copy("alter_name","neuer_name");
# rename("alter_name","neuer_name");
# chdir('../.');

unless ($#ARGV >0) 
{print "***************************************************\n";
 print " program searchspace used to search parameter space\n\n";
 print " usage: searchspace i [-jpglog 3.5e-1 file.jpg] * \n\n";
 print " i ... level of parameter search (i=0 for first scan, for i>0 a\n";
 print " file searchspace.i with parameter sets must exist)\n\n";
 print " * .. filename(s) of paramter file(s)\n";
 print " \n... there must exist a program calcsta, the output of which\n";
 print " contains 'sta = 249' - this program is called from this output\n";
 print " at every iteration step and  the standard sta deviation is \n";
 print " recorded in filesearchspace.i  i=1,2,3, ....\n\n";
 print " localminima are also recorded (searchspace.i.localminima)\n\n";
 print STDOUT << "EOF";
 option -jpglog: if sta is less then 3.5e-1, then the image file.jpg 
 (which should be created by calcsta) is copied to results/parsetnr.jpg


EOF
 print " Note: at each level searchspace covers parameterspace\n";
 print "       with a more dense net of points. The minimum distance\n";
 print "       of these points is taken from the parameter-stepwidths\n\n";
 print "***************************************************\n";
 print " alternative usage:  searchspace -27.03 results/searchspace.3 * \n\n";
 print " If used with option -n, the program reads a parameter set \n";
 print " nr 27.03 from the searchspace file (e.g. results/searchspace.3)\n";
 print " and puts this parameter set to the parameter files *\n";
 print " and *.forfit\n\n";
 print "***************************************************\n";
 print " <Press enter to close>";$in=<STDIN>;
 exit 0;}

#in dos version
# mv -> rename
# cp -> copy
# .\calcsta -> calcsta.bat
# rm -> del


#********************************************************************************
 #load from files * parameters, stepwidths, interval and copy original files to *.fit
#format STDOUT_TOP =
#parameter   [ value,     min,      max,     err,       stp     ]   
#.
format STDOUT =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf("%s [%+e,%+e,%+e,%+e,%+e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i])
.
format Fout =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf ("%s [%+e,%+e,%+e,%+e,%+e]",$parnam[$ii],$par[$ii],$parmin[$ii],$parmax[$ii],$parerr[$ii],$parstp[$ii])
.

 @parnam=();@par=();@parmin=();@parmax=();@parerr=();@parstp=();@parav=();
 $ARGV[0]=~s/x/*/g; $searchlevel=eval $ARGV[0]; shift @ARGV;
  if($searchlevel<0){# in this case only read line number -$searchlevel from input file and update parameter files
                     $inputfile=$ARGV[0]; shift @ARGV;
                    }
# treat option -jpglog
$jpglog=0;
 if($ARGV[0]=~/-jpglog/)
{shift @ARGV;$ARGV[0]=~s/x/*/g; $jpglog=eval $ARGV[0];shift @ARGV;
 $jpgimagefile=$ARGV[0];shift @ARGV;
}


 while(!open(Fout,">results/searchspace.status")){print "Error opening file results/searchspace.status\n";<STDIN>;}
  print Fout "parameter[value,      min,           max,           (not used)   ,minimum meshwidth]\n";
 foreach (@ARGV)
 {$file=$_; if(mycopy ($file,$file.".bak")){print "\n error copying $file \n";print " <Press enter to close>";$in=<STDIN>;exit 1;} 
   unless (open (Fin, $file.".forfit")){print "\n error:unable to open $file\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}   
   while($line=<Fin>)
{while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/) {++$#par;#load another parameter
				 ($parname)=($line=~m/(?:#!|[^#])*?\b(par\w+)\s*\Q[\E/);
                                 foreach(@parnam){if (/$parname/){print "ERROR searchspace: parameter $parname occurs more than one time in input files\n"; print " <Press enter to close>";$in=<STDIN>;exit 1;}}
                                 $parnam[$#par]=$parname;
				 ($par[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*([^,]+)/);
				 ($parmin[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*([^,]+)/);
				 ($parmax[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parerr[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parstp[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^\Q]\E]+)/);
                                 $i=$#par;write STDOUT;$ii=$i;write Fout;
				  #check if parmin<=parmax
                          if ($parmin[$#par]>$parmax[$#par]) 
                            {print "ERROR searchspace reading parameterrange: parmin > parmax\n";
                              print " <Press enter to close>";$in=<STDIN>;exit 1;
                             }

                         $line=~s/(?:#!|[^#])*?\bpar\w+\s*\Q[\E//;                     
				 }
     } close Fin;
 }  
    if ($#par<0) {print "Error searchspace: no parameters found in input files @ARGV\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
    print ($#par+1);print " parameters found\n"; close Fout;

  if($searchlevel<0){# in this case only read parameter set -$searchlevel from input file and update parameter files
                    print "reading points from file $inputfile\n";
                    open (Fin1,$inputfile);
                    while($line=<Fin1>)
                    {if ($line=~/^\s*#/) {;}
                     else{$line=~s/D/E/g;@numbers=split(" ",$line);
                          if($numbers[0]==-$searchlevel){
                                                $i=0;foreach(@par){$par[$i]=$numbers[$i+1];++$i}
                                                $sta=$numbers[$i+1]; #sta at input parameter position
                                                writefiles();
                                                # output parameter information to status file and to screen
                                                open(Fout,">results/searchspace.status");
                                                print Fout ".......... searchspace stopped...............\n"; 
                                                print Fout ".............................................\n";
                                                print Fout "parameter set nr ".(-$searchlevel)." read from $inputfile\n";
                                                print Fout ".............................................\n";
                                                print Fout "parameter[value,      min,           max,           (not used)   ,minimum meshwidth]\n";
                                                $ii=0;foreach(@par){write Fout;++$ii;}
                                                print Fout "..............................sta=$sta.......\n";
                                                print Fout "these values are stored in parameter files * and *.forfit \n";
                                                print Fout ".............................................\n";
                                                close Fout;
                                                print STDOUT ".......... searchspace stopped...............\n"; 
                                                print STDOUT ".............................................\n";
                                                print STDOUT "parameter set nr ".(-$searchlevel)." read from $inputfile\n";
                                                print STDOUT ".............................................\n";
                                                print STDOUT "parameter[value,      min,           max,           (not used)   ,minimum meshwidth]\n";
                                                $i=0;foreach(@par){write STDOUT;++$i;}
                                                print STDOUT ".............................................\n";
                                                print STDOUT "these values are stored in parameter files * and *.forfit \n";
                                                print STDOUT ".............................sta=$sta........\n";
                                                exit 0;
                                                }
                         }
                    }
                   print "Error program searchspace, loading parameter set ".(-$searchlevel)." failed !\n"; 
                   print " <Press enter to close>";$in=<STDIN>;exit 1;
                   }


#*******************************************************************************
print" starting search of level $searchlevel ...\n";
$minnumber=1;$pointcounter=0;
$starttime=time;$stamin=1e10;
# open and initialize output file
print "storing points in file searchspace.$searchlevel\n";
open (Foutlevel,">results/searchspace.$searchlevel");
print Foutlevel "#ParSetNr ";foreach(@parnam){print Foutlevel $_." ";}print Foutlevel "sta variance chisquared\n";
open(Foutlocalmin,">results/searchspace.$searchlevel.localminima");
print Foutlocalmin "#ParSetNr ";foreach(@parnam){print Foutlocalmin $_." ";}print Foutlocalmin "sta\n";

if ($searchlevel==0)
{$pointstocalculate=($#par+1)*2;
            $i=0;foreach(@par){$par[$i]=$parmin[$i]+($parmax[$i]-$parmin[$i])/2;++$i}
            @parorigin=@par;$pointcounterorigin=$pointcounter+1;
            ($sta)=sta(); #sta at input parameter position
       $staorigin=$sta;$i=0;$minimum=1;
       foreach(@par)
        {$dpar=($parmax[$i]-$parmin[$i])/2**($searchlevel+2);
         if (abs($dpar)>abs($parstp[$i]))
         {
          $par[$i]+=$dpar;
          ($sta)=sta();
           read_write_statusfile();
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; # foreach (@ARGV)
                                          # {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                          # mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
                                          $pointcountermin=$pointcounter+1;
			}

          if ($sta<=$staorigin){$minimum=0;}
           write_set();++$pointcounter; 
          $par[$i]-=2*$dpar;
          ($sta)=sta();
          read_write_statusfile();
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; #foreach (@ARGV)
                                       #{$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                       # mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
                                          $pointcountermin=$pointcounter+1;
			}
          if ($sta<=$staorigin){$minimum=0;}
           write_set();++$pointcounter;
          $par[$i]+=$dpar;
         }
         last if ($sta==0);
        ++$i;
       }
       if ($pointcounter!=0)
       { if ($minimum==1&&$sta!=0)
        {#the input parameters are probably a (local) minimum - so save to file results/searchspace.$searchlevel.localminima
          ($sta)=sta();
#        last if ($sta==0);
       if ($sta<$stamin){$stamin=$sta; #foreach (@ARGV)
                                       #{$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                       # mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
                                          $pointcountermin=$pointcounter+1;
			}
#          foreach (@ARGV){$file=$_; mycopy ($file.".forfit",$file.".$searchlevel.$minnumber");} # do not save these lots of files, this is confusing
          print Foutlocalmin $pointcounterorigin." ";
          $ii=0;foreach(@parorigin){$dd=sprintf("%e ",$parorigin[$ii]);print Foutlocalmin $dd;++$ii} print Foutlocalmin $staorigin."\n";
	  ++$minnumber;
       }
        }
}
else
{print "reading points from file searchspace.".($searchlevel-1)."\n";
 $pointstocalculate=0;open (Fin1,"results/searchspace.".($searchlevel-1));
 while($line=<Fin1>){if ($line=~/^\s*#/) {;} else {++$pointstocalculate;}}close Fin1;
 $pointstocalculate*=2*($#par+1);
 open (Fin1,"results/searchspace.".($searchlevel-1));
 while($line=<Fin1>)
  {
   if ($line=~/^\s*#/) {;}
   else{$line=~s/D/E/g;@numbers=split(" ",$line);
            $pointcounterorigin=$numbers[0];$i=0;foreach(@par){$par[$i]=$numbers[$i+1];++$i}
            $sta=$numbers[$i+1]; #sta at input parameter position
       $staorigin=$sta;$i=0;$minimum=1;@parorigin=@par;
       foreach(@par)
        { $dpar=($parmax[$i]-$parmin[$i])/2**($searchlevel+2);
         if (abs($dpar)>abs($parstp[$i]))
         {$par[$i]+=$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;close Fin;
          if ($line=~/exiting searchspace/){$sta=0;}
          last if ($sta==0);
          read_write_statusfile();
          if ($sta<$stamin){$stamin=$sta;# foreach (@ARGV)
                                       #{$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                       # mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
                                          $pointcountermin=$pointcounter+1;
			}

          if ($sta<=$staorigin){$minimum=0;}
           write_set();++$pointcounter; 
          $par[$i]-=2*$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;close Fin;
          if ($line=~/exiting searchspace/){$sta=0;}
          last if ($sta==0);
          read_write_statusfile();
          if ($sta<$stamin){$stamin=$sta; #foreach (@ARGV)
                                       #{$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                       # mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
                                          $pointcountermin=$pointcounter+1;
			}
          if ($sta<=$staorigin){$minimum=0;}
           write_set();++$pointcounter; 
          $par[$i]+=$dpar;
          }
         last if ($sta==0);
        ++$i;
       }
         last if ($sta==0);
         last if ($pointcounter==0);
       if ($minimum==1&&$sta!=0)
        {#the input parameters are probably a (local) minimum - so save a file
          ($sta)=sta();if ($sta>$staorigin*1.00000001){$warning="#WARNING - calculation of sta gave another result at some point than stored in file searchspace.".($searchlevel-1)."\n";}
       if ($sta<$stamin&&$sta!=0){# $stamin=$sta;# foreach (@ARGV)
                                       #{$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                       # mycopy ($file.".forfit",$file.".forfit.min.$searchlevel");}
                                       #   $pointcountermin=$pointcounter+1;
			}
#          foreach (@ARGV){$file=$_; mycopy ($file.".forfit",$file.".$searchlevel.$minnumber");}
          print Foutlocalmin $pointcounterorigin." ";
          $ii=0;foreach(@parorigin){$dd=sprintf("%e ",$parorigin[$ii]);print Foutlocalmin $dd;++$ii} print Foutlocalmin $staorigin."\n";
	  ++$minnumber;
        }
       }
    
  }
  close Fin1;
}
$hours=(time-$starttime)/3600;
$estimate=sprintf("%6.2f",$hours*2*($#par+1));
print "nofpts=$pointcounter points calculated in  $hours h.\n Time estimate for next level ".($searchlevel+1).": $estimate h\n";
print Foutlevel "#$pointcounter points calculated in  $hours h.\n# Time estimate for next level ".($searchlevel+1).": $estimate h\n";
print Foutlevel "#".($#ssta+1)." contributions to sta found in output of calcsta ...\n";
   if($chisquared){print Foutlevel "#sta=chisquared(=1/".($#ssta+1)."sum deviations^2/experrors^2)\n";}
               else{print Foutlevel "#sta=variance s^2(=1/".($#ssta+1)."sum deviations^2)\n";}
print Foutlevel "#Mininimal found sta minsta=$stamin for parameter set nr ".($pointcountermin+$searchlevel/100)." in file  results/searchspace.$searchlevel \n";
print Foutlevel "#noflocalminima=".($minnumber-1)." local minima found (see results/searchspace.$searchlevel.localminima)\n";
print Foutlocalmin "#noflocalminima=".($minnumber-1)." local minima found (see results/searchspace.$searchlevel.localminima)\n";
close Foutlevel;close Foutlocalmin;
# foreach (@ARGV) {$file=$_; mycopy($file.".min.$searchlevel",$file);}
     open(Fout,">results/searchspace.status");
print Fout ".......... searchspace stopped...............\n"; 
print Fout ".............................................\n";
print Fout "nofpts=$pointcounter points calculated in  $hours h.\n";
print Fout ".............................................\n";
print Fout ($#ssta+1)." contributions to sta found in output of calcsta ...\n";
if($chisquared){print Fout "sta=chisquared(=1/".($#ssta+1)."sum deviations^2/experrors^2)\n";}
               else{print Fout "sta=variance s^2(=1/".($#ssta+1)."sum deviations^2)\n";}
print Fout "Mininimal found sta=$stamin for parameter set nr ".($pointcountermin+$searchlevel/100)." in file  results/searchspace.$searchlevel \n";
print Fout "noflocalminima=".($minnumber-1)." local minima found (see results/searchspace.$searchlevel.localminima)\n";
print Fout ".............................................\n";
print Fout "Time estimate for next level ".($searchlevel+1).": $estimate h\n";
print Fout ".............................................\n";
close Fout;
exit 0;

# END OF MAIN PROGRAM
#****************************************************************************** 


sub sta {local $SIG{INT}='IGNORE';
 writefiles();
 # print "#call routine calcsta to calculate standard deviation\n";
 if ($^O=~/MSWin/){
                   if(system ("calcsta.bat 1e10 > results\\searchspace.sta")){print "\n error executing calcsta.bat\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
                  }
 else
                  {
                   if(system ("./calcsta 1e10 > results/searchspace.sta")){print "\n error executing calcsta.bat\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
                  }

 open (Fin,"./results/searchspace.sta"); $i6=0;$errc=1;
 while($line=<Fin>){
           if($line=~/^(#!|[^#])*?\bsta\s*=/) {($staline)=($line=~m/(?:#!|[^#])*?\bsta\s*=\s*([\d.eEdD\Q-\E\Q+\E\s]+)/);
                                               $staline=~s/D/E/g;my @ss=split(" ",$staline);
                                               $ssta[$i6]=$ss[0];
                                               if($errc==1){if ($#ss>0){$eerr[$i6]=$ss[1];}else{$errc=0;}}
                                               ++$i6;
                                              }
                   }
 close Fin;

 mydel("./results/searchspace.sta");
 $delta= sqrt PDL->new(@ssta);
 $c=PDL->new(@par);
# print $delta;
 $s2=inner($delta,$delta)/($#ssta+1); # this is s^2
 $sta=$s2;
 if($errc>0) #if errors are given we can minimize chisquared and calculate covariance matrix
 {my  $err=   sqrt PDL->new(@eerr);
  $delta=$delta/$err;
  $chisquared=inner($delta,$delta)/($#ssta+1); # this is chisquared
  # if we have errors present we rather minimize chi2
  $sta=$chisquared;
 }
 return $sta;
}
  
sub writefiles {
 #print "#write modified parameterset to files *\n";
 foreach (@ARGV)
 {$file=$_; open (Fin, $file.".forfit");open (Fout1, ">".$file);open (Fout2,">results/searchspace.par");
   while($line=<Fin>)
     {$modline=$line;
      if ($line=~/^(#!|[^#])*?\bpar/) {#here write modified parameter set to line
                            while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/)
			          {my $i=0;#insert a parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     if ($line=~/^(#!|[^#])*?\b$pnam\s*\Q[\E/)
                                        {$line=~s|$pnam\s*\Q[\E[^\Q]\E]*\Q]\E|$par[$i] |;
                                   $dd=sprintf("%s [%e,%e,%e,%e,%e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i]);
                                   $modline=~s|$pnam\s*\Q[\E[^\Q]\E]+\Q]\E|$dd|;
					}
				     ++$i;
				    }
				  }

                             while ($line=~/^(?:#!|[^#])*?\bfunction\s*\Q[\E(.*?)\Q]\E/)
			          {$expression=" ".$1." ";
				   #substitute into the expression the paramter values
				   my $i=0;#insert a function parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     while ($expression=~/.*\W$pnam\W/)
                                        {$expression=~s|$pnam|($par[$i])|;
					}
				     ++$i;
				    }
				    # calculate the expression by a little perl program
				    #open (Foutcc, ">./results/ccccccc.ccc");
				    #printf Foutcc "#!/usr/bin/perl\nprint ".$expression.";\n";
				    #close Foutcc;$systemcall="perl ./results/ccccccc.ccc > ./results/cccccc1.ccc";
                                    #if ($^O=~/MSWin/){$systemcall=~s|\/|\\|g;}
				    #if(system $systemcall){print "error evaluating expression in results/ccccc*";print " <Press enter to close>";$in=<STDIN>;exit 1;}
				    #open (Fincc,"results/cccccc1.ccc");
				    #$data=<Fincc>; close Fincc;
				    #mydel("./results/ccccccc.ccc");mydel("./results/cccccc1.ccc");
				    $data= eval $expression;
                                    # $data contains now the result of the mathematical expression
                                    $line=~s|function\s*\Q[\E.*?\Q]\E|$data |;
				   }



                            } print Fout1 $line;print Fout2 $modline;
     } close Fin;close Fout1;close Fout2;
     if (mycopy("results/searchspace.par",$file.".forfit"))
     {print "\n error copying results/searchspace.par to  $file.forfit\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
 }
}

sub write_set()
{       my $dd=sprintf("%.2f ",($pointcounter+1)+$searchlevel/100);print Foutlevel $dd;
        my $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta." ".$s2." ".$chisquared."\n";
        if($jpglog>$sta){$dd=sprintf("%g.jpg",($pointcounter+1)+$searchlevel/100);
                         print Foutlevel '#<img src="'.$dd.'">'."\n";
                         mycopy($jpgimagefile,"./results/".$dd);
                        }
}

sub mycopy { my ($file1,$file2)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;$file2=~s|\/|\\|g;
                               return system("copy ".$file1." ".$file2);
                              }
                 else
                              {return system("cp -f ".$file1." ".$file2);
                              }

           }

sub mydel  { my ($file1)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;
                               return system("del ".$file1);
                              }
                 else
                              {return system("rm ".$file1);
                              }

           }
    
sub read_write_statusfile {
          open(Fin,"./results/searchspace.status");$line=<Fin>;close Fin;
          if ($line=~/exiting searchspace/){$sta=0;}
          else
     {open(Fout,">./results/searchspace.status");$ii=0;
     $est=sprintf("%6.2f",($pointstocalculate-$pointcounter)*(time-$starttime)/3600/($pointcounter+1));
     $passed=sprintf("%6.2f",(time-$starttime)/3600);
     print Fout ($#ssta+1)." contributions to sta found in output of calcsta ...\n";
     print Fout "Searchlevel $searchlevel parameter set nr ".($pointcounter+1)."  Current sta=$sta \n";
     print Fout "Minimum sta=$stamin for parameter set nr ".($pointcountermin+$searchlevel/100)."\n";
    if($chisquared){print Fout "sta=chisquared(=1/".($#ssta+1)."sum deviations^2/experrors^2)\n";}
               else{print Fout "sta=variance s^2(=1/".($#ssta+1)."sum deviations^2)\n";}
     print Fout "--------------------------------------------------------------------\n";
     print Fout "Time since Start: $passed h Estimated time to complete: $est h\n";
     print Fout "--------------------------------------------------------------------\n";
     print Fout "parameter[value,      min,           max,           (not used)   ,minimum meshwidth]\n";
     foreach(@par){write Fout;++$ii;}
     print Fourt $warning;
     close Fout;}

                          }