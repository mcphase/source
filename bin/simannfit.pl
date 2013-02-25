#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use PDL;

use PDL::Slatec;

unless ($#ARGV >0) 
{print " program simannfit used to perform simulating annealing\n";
 print " usage: simannfit 10 [options] * \n  10 .. initial temperature \n * .. filename(s) of paramter file(s)\n";
 print " ... there must exist a program calcsta.bat (windows) / calcsta (linux), the output of which\n";
 print " contains 'sta = 249' - this program is called from this output\n";
 print " at every iteration step and  the standard sta deviation is minimzed\n";
 print " simannfit gives as a parameter a maximum number for sta - if during the\n";
 print " calculation of sta in calcsta.bat this number is exceeded calcsta.bat can exit (saves time)\n";
 print " option -t sets time limit until program end (in seconds)\n";
 print " option -s gives maximal number of iteration steps to be done\n";
 print " <Press enter to close>";$in=<STDIN>;
 exit 0;}


# extra version for dos because operating system commands are present:
#system("alias 'copy'='cp -f'");
#system("alias 'rename'='mv'");
#system("alias 'del'='rm'");
#system("alias 'calcsta'='./calcsta'");

#********************************************************************************
 #load from files * parameters, stepwidths, interval and copy original files to *.fit
#format STDOUT_TOP =
#parameter   [ value,     min,      max,     err,       stp     ]   
#.
format STDOUT =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf ("%s [%+e,%+e,%+e,%+e,%+e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i])
.
format Fout =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf ("%s [%+e,%+e,%+e,%+e,%+e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i])
.

 @parnam=();@par=();@parmin=();@parmax=();@parerr=();@parstp=();@parav=();@thisparstp=();
				 @parhisto=();@parhistostp=();@perlhistostart=();$hh=0;
  $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$stattemp=eval $ARGV[0]; shift @ARGV;
  $starttime=time;$maxtim=1e10;$maxstep=1e24;
  if ($ARGV[0]=~"-t") {shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$maxtim=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0]=~"-s") {shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$maxstep=eval $ARGV[0]; shift @ARGV;}

 while(!open(Fout,">results/simannfit.status")){print "Error opening file results/simannfit.status\n";<STDIN>;}
   print Fout "parameter[value,      min,           max,           variation,     stepwidth]\n";
  foreach (@ARGV)
 {$file=$_; if(mycopy ($file,$file.".bak")){print "\n warning copying $file not possible - press enter to continue\n";<stdin>;}
   unless (open (Fin, $file.".forfit")){die "\n error:unable to open $file.forfit\n";}   
    while($line=<Fin>)
 {while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/) {++$#par;#load another parameter
				 ($parname)=($line=~m/(?:#!|[^#])*?\b(par\w+)\s*\Q[\E/);
                                 foreach(@parnam){if (/$parname/){print "ERROR simannfit: parameter $parname occurs more than one time in input files\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}}
                                 $parnam[$#par]=$parname;
				 ($par[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*([^,]+)/);
				 ($parmin[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*([^,]+)/);
				 ($parmax[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parerr[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parstp[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^\Q]\E]+)/);
                         $parstp[$#par]=abs($parstp[$#par]);
                                 $i=$#par;write STDOUT;write Fout;
				 #check if parmin<=parmax
                          if ($parmin[$#par]>$parmax[$#par]) 
                            {print "ERROR simannfit reading parameterrange: parmin > parmax\n";
                             print " <Press enter to close>";$in=<STDIN>;exit 1;
                             }
                          if($parstp[$i]<=0)
                            {print "ERROR simannfit reading parameterstepwidth: parstp <=0 \n";
                             print " <Press enter to close>";$in=<STDIN>;exit 1;
                             }

                         $line=~s/(?:#!|[^#])*?\bpar\w+\s*\Q[\E//;                     
				 $perlhistostart[$i]=$hh;# histogramm steps (not more than 1000)
                                 $parhistostp[$i]=$parstp[$i];if(($parmax[$i]-$parmin[$i])/$parstp[$i]>1000){$parhistostp[$i]=($parmax[$i]-$parmin[$i])/1000;}
				  for($hx=0;$hx<=int(($parmax[$i]-$parmin[$i])/$parhistostp[$i])+1;++$hx)
                                   {$parhisto[($hx+$perlhistostart[$i])]=0;}
                                  $hh+=int(($parmax[$i]-$parmin[$i])/$parhistostp[$i])+2;

                                 }
     } close Fin;
 }  
    if ($#par<0) {print "Error simannfit: no parameters found in input files @ARGV\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
   close Fout;
print "initialize parameter storage\n";
$parstore = zeroes $#par+3,$#par+1;
$deltastore =  PDL->nullcreate(0);
$nof_calcsta_calls=0;
#*******************************************************************************
# fitting loop
 print ($#par+1);print " parameters found - testing calculation of sta\n";
$rnd=1;$stasave=1e20;
 ($sta)=sta();$stps=1;$noofupdates=0;$stepnumber=0;

if($sta>0)
{print "starting fit\n";
 $SIG{INT} = \&catch_zap;
 while($sta>0)
 {  $stasave=$sta;
 # modify parameters
 print "\n ...next fitting loop ...\n";
 @parsav=@par;$i=0;
 foreach(@par){$rnd=rand;$thisparstp[$i]=($rnd-0.5)*$parstp[$i]*$stps;$par[$i]+=$thisparstp[$i];
               if ($par[$i]<$parmin[$i]) {$thisparstp[$i]=$parmin[$i]-$par[$i];$par[$i]=$parmin[$i];}
	       if ($par[$i]>$parmax[$i]) {$thisparstp[$i]=$parmax[$i]-$par[$i];$par[$i]=$parmax[$i];}
               write STDOUT; 
	       ++$i;}
 $rnd=rand;
   ($sta)=sta(); # CALCULATE sta !!!!
   ++$stepnumber;
   print " ...  current sta=$sta, statistical T=$stattemp, step ratio=$stps\nsta of stored parameters=$stasave\n";
   open(Fin,"./results/simannfit.status");$line=<Fin>;
    if ($line=~/exiting simannfit/){$sta=0;close Fin;}
    else
    {close Fin;
     read_write_statusfile();
    }

 if (time-$starttime>$maxtim){$sta=0;print "\n maximum time for fitting reached - stopping fit\n";}
 if ($stepnumber>$maxstep){$sta=0;print "\n maximum step number for fitting reached - stopping fit\n";}
 if ($sta==0) {#recover ol pars
      @par=@parsav; 
      }
 last if ($sta==0);
 if ($sta>$stasave)
  {if($rnd>exp(-($sta-$stasave)/$stattemp))
     {#recover ol pars and adapt parstep to step not so big in this direction
      @par=@parsav;$i=0;foreach(@parstp){if($parstp[$i]>$parhistostp[$i]/1000){$parstp[$i]-=0.1*abs($thisparstp[$i]);}++$i;}
      $stps*=0.999;$sta=$stasave; 
      if ($stps<0.01){$stps=10;}# if stepwidth decreased too much make large steps again to get out of side minimum !!!
     }
   else
     {$stattemp=$stattemp*0.999;
     }     
  }  
 else
  {$stattemp=$stattemp*0.995;
   #update errors
   $i=0;
   print "STA DECREASED in the last LOOP:\n";
   foreach(@par){$p=$_;
   $parav[$i]=($parav[$i]*$noofupdates + $p)/($noofupdates+1);
   $parerr[$i]=sqrt($parerr[$i]*$parerr[$i]*$noofupdates+
                   ($p-$parav[$i])*($p-$parav[$i]))/($noofupdates+1);     
   $parstp[$i]+=0.1*abs($thisparstp[$i]); # adapt parstp to be more bold in this direction
   $hx=int(($p-$parmin[$i])/$parhistostp[$i]);
   ++$parhisto[($hx+$perlhistostart[$i])];
   open(Fout,">./results/".$parnam[$i].".hst");
   print Fout "#{Histogram of parameter ".$parnam[$i]."\n# value vs. number of  occurrences in good solutions (sta decreased)}\n";
   for($hx=0;$hx<=int(($parmax[$i]-$parmin[$i])/$parhistostp[$i])+1;++$hx)
   {print Fout (($hx+0.5)*$parhistostp[$i]+$parmin[$i])."   ".($parhisto[($hx+$perlhistostart[$i])])."\n";
   } close Fout;

   ++$i;} ++$noofupdates;
   #printout current parameters
   $i=0;foreach(@par){write STDOUT;++$i;}
   print "   sta of stored parameters=$sta\n";
  }	       
 print "*";
 }	       

   print "best fit:\n";
   $i=0;foreach(@par){write STDOUT;++$i;}
   ($sta)=sta(); # CALCULATE sta !!!!
}
else
{print "sta=0 already - not fit required !?\n";exit;}
# calculate covariance matrix
# for($i6=0;$i6<=$#par;++$i6){set $parstore,0,$i6,$par[$i6];}
#$parstore
#$deltastore

print $parstore;
print $store_counter."\n";
# rotate back the storage to make last set in column 0 ...
for($i6=1;$i6<=$store_counter;++$i6){
$deltastore=rotate $deltastore,-1;
$parstore= rotate $parstore,-1;
}

print $parstore;
$b=$parstore->slice(0)->copy;
$parstore.=$parstore-$b;
print $parstore;
$V=$parstore->slice('1:-1');
$c=$deltastore->slice(0)->copy;
$deltastore.=$deltastore-$c;
$delta=$deltastore->slice('1:-1');
#print $V; # here we have calculated V
#print $delta;
# now we need to cut out the zero vector (if present) !!!
$i6=$#par+3;
while(sum($V->slice(0)*$V->slice(0))>1e-10&&$i6>0){
$V=rotate $V,1;$delta=rotate $delta,1;--$i6;
}
$delta=$delta->slice('1:-1');
$V=$V->slice('1:-1');
print $V;
print $delta;

#{$cov="calculation of covariance matrix not successfull because last n steps of simulated annealing were not orthogonal in parameter space - restart simannfit and try again ...\n";}
#print $cov;
#print $Fij;

   print "best fit:\n";
   $i=0;foreach(@par){write STDOUT;++$i;}
  if($chisquared){print "      sta=chi2=$sta (=sum deviations^2/experrors^2)\n    variance s2=$s2 (=sum deviations^2)\n";
     print  "Covariance matrix( may be not successfull because last n steps of simulated annealing may be not necessarily\n orthogonal in parameter space - if this happens restart and try again):\n";
     }
          else {  print "      sta=variance=s2=$s2 (=sum deviations^2)\n";}


# move files
# foreach (@ARGV)
# {$file=$_; #if(mydel("$file.fit")){die "\n error deleting $file.fit \n";}
#            #if(mycopy ($file,$file.".fit")){die "\n error copying file $file\n";}
#            #if(mydel($file)){die "\n error deleting $file \n";}
#            #if(mycopy ($file.".par ",$file)){die "\n error copying $file.par \n";}
# }
     open(Fout,">results/simannfit.status");print Fout " ... simannfit stopped\n";
     print Fout ($#ssta+1)." contributions to sta found in output of calcsta ...\n";
     print Fout "best fit:\n";
  if($chisquared){print Fout "      sta=chi2=$sta (=sum deviations^2/(".($#ssta+1)."*experrors^2))\n    variance s2=$s2 (=sum deviations^2/".($#ssta+1).")\n";}
          else {  print Fout "      sta=variance=s2=$s2 (=sum deviations^2/".($#ssta+1).")\n";}
     print Fout "----------------------------------------------------------------------------------------\n";
     print Fout "Statistical Temperature=$stattemp      Step Ratio=$stps\n";
     print Fout "----------------------------------------------------------------------------------------\n";
     $est=sprintf("%6.2f",(time-$starttime)/3600);$maxtimest=sprintf("%6.2f",($maxtim)/3600);
     print Fout "Time since start of simannfit: $est hours (limit:$maxtimest), $stepnumber steps (limit:$maxstep)\n";
     print Fout "----------------------------------------------------------------------------------------\n";
     print Fout "parameter[value,      min,           max,           variation,     stepwidth]\n";
     $i=0;     foreach(@par){$parcent=int(10*($par[$i]-$parmin[$i])/(1e-10+$parmax[$i]-$parmin[$i]));
                        print Fout "|";for($jsw=0;$jsw<=9;++$jsw){
                                       if ($jsw==$parcent){print Fout "*";}else{print Fout "-";}
                                                                 }
                        print Fout "|";print Fout sprintf ("%s [%+e,%+e,%+e,%+e,%+e]\n",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i]);
                   ++$i;}
  if($chisquared){
     print Fout "Covariance matrix( may be not successfull because last n steps of simulated annealing may be not necessarily\n orthogonal in parameter space - if this happens restart and try again):\n";
$Fij=$delta x matinv($V);
 $FtF=$Fij->xchg(0,1) x $Fij;
 $cov=$sta*matinv($FtF);  # multiply by sta=chi2 in order to get covariance matrix
     print $cov; print Fout $cov;
                $i=0;foreach(@par){print $parnam[$i]." error=".(sqrt($cov->at($i,$i)))."\n";
                                   print Fout $parnam[$i]." error=".(sqrt($cov->at($i,$i)))."\n";
                                   ++$i;}
     }
     close Fout;
print " <Press enter to close>";$in=<STDIN>;exit 0;
# END OF MAIN PROGRAM
#****************************************************************************** 

sub catch_zap {
 my $signame = shift;
# foreach (@ARGV)
# {$file=$_; mycopy($file,$file.".fit");
#            mycopy($file.".par",$file);}
 die "SIG$signame stopping fit\n";
}
#****************************************************************************** 


sub sta {#local $SIG{INT}='IGNORE';
 #print "#write modified parameterset to files *\n";
 foreach (@ARGV)
 {$file=$_; open (Fin, $file.".forfit");open (Fout1, ">".$file);open (Fout2,">./results/simannfit.par");
   while($line=<Fin>)
     {$modline=$line;
      if ($line=~/^(#!|[^#])*?\bpar/) {#here write modified parameter set to line
                            while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/) 
			          {$i=0;#insert a parameter
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
				   $i=0;#insert a function parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     while ($expression=~/.*\W$pnam\W/)
                                        {$expression=~s|$pnam|($par[$i])|;
					}
				     ++$i;
				    }
				    # calculate the expression by a little perl program
				    #unless(open (Foutcc, ">./results/ccccccc.ccc")){die "error simannfit: could not open temporary file results/ccccc.ccc/n";}
				    #printf Foutcc "#!/usr/bin/perl\nprint ".$expression.";\n";
				    #close Foutcc;
				    #close Foutcc;$systemcall="perl ./results/ccccccc.ccc > ./results/cccccc1.ccc";
                                    #if ($^O=~/MSWin/){$systemcall=~s|\/|\\|g;}
				    #if(system $systemcall){die "error evaluating expression in results/ccccc*";}
				    #unless(open (Fincc,"./results/cccccc1.ccc")){die "error simannfit: could not open temporary file results/ccccc1.ccc/n";}
				    #$data=<Fincc>; close Fincc;
				    #mydel ("./results/ccccccc.ccc");
                                    #mydel ("./results/cccccc1.ccc");
				    $data=eval $expression;
                                    # $data contains now the result of the mathematical expression
                                    $line=~s|function\s*\Q[\E.*?\Q]\E|$data |;
				   }



                            } print Fout1 $line;print Fout2 $modline;
     } close Fin;close Fout1;close Fout2;
     if (mycopy("./results/simannfit.par",$file.".forfit"))
     {die "\n error copying results/simannfit.par to  $file.forfit\n";}
 }

# print "#call routine calcsta.bat to calculate standard deviation\n";
$staboundary=$stasave-log($rnd+1e-10)*$stattemp;
 if ($^O=~/MSWin/){
                   if(system ("calcsta.bat $staboundary > results\\simannfit.sta")){die "\n error executing calcsta.bat\n";}
                  }
 else
                  {
                   if(system ("./calcsta $staboundary > ./results/simannfit.sta")){die "\n error executing calcsta\n";}
                  }	
 open (Fin,"./results/simannfit.sta");  $i6=0;$errc=1;
 while($line=<Fin>){
           if($line=~/^(#!|[^#])*?\bsta\s*=/) {($staline)=($line=~m/(?:#!|[^#])*?\bsta\s*=\s*([\d.eEdD\Q-\E\Q+\E\s]+)/);
                                               $staline=~s/D/E/g;my @ssn=split(" ",$staline);
                                               $ssta[$i6]=$ssn[0];
                                               if($errc==1){if ($#ssn>0){$eerr[$i6]=$ssn[1];}else{$errc=0;}}
                                               ++$i6;
                                              }
                   }
 close Fin;
 mydel ("./results/simannfit.sta"); 
# print @ssta;
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

# $sta= ... sum of @ssta
# $deltastore= ... @ssta
# $parstore= ....@par
if($nof_calcsta_calls<$#par+3)
{# extend storage of delta
 $deltastore=$deltastore->append($delta->dummy(0));
 $store_counter=$nof_calcsta_calls;
}
else
{#rotate
# $deltastore=rotate $deltastore,1; this would be good but does not work pdl bug
 if ($store_counter>$#par+2){$store_counter=0;}
 for($i6=0;$i6<=$#ssta;++$i6){set $deltastore,$store_counter,$i6,$delta->at($i6);}
}
# $parstore= rotate $parstore,1;  this would be good but does not work pdl bug
 for($i6=0;$i6<=$#par;++$i6){set $parstore,$store_counter,$i6,$par[$i6];}
 ++$nof_calcsta_calls;++$store_counter; print $store_counter." ".$#par."\n";
 return $sta;
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
     open(Fout,">./results/simannfit.status");$i=0;
     print Fout ($#ssta+1)." contributions to sta found in output of calcsta ...\n";
     if($chisquared){print Fout "Current sta=chi2=$sta (=sum deviations^2/(".($#ssta+1)."*experrors^2)) sta of stored parameters=$stasave\n";}
              else {print Fout " Current     sta=variance=s2=$s2 (=sum deviations^2/".($#ssta+1).")   sta of stored parameters=$stasave\n";}
     print Fout "----------------------------------------------------------------------------------------\n";
     print Fout "Statistical Temperature=$stattemp      Step Ratio=$stps\n";
     print Fout "----------------------------------------------------------------------------------------\n";
     $est=sprintf("%6.2f",(time-$starttime)/3600);$maxtimest=sprintf("%6.2f",($maxtim)/3600);
     print Fout "Time since start of simannfit: $est hours (limit:$maxtimest), $stepnumber steps (limit:$maxstep)\n";
     print Fout "----------------------------------------------------------------------------------------\n";
     print Fout "parameter[value,      min,           max,           variation,     stepwidth]\n";
     foreach(@par){$parcent=int(10*($par[$i]-$parmin[$i])/(1e-10+$parmax[$i]-$parmin[$i]));
                        print Fout "|";for($jsw=0;$jsw<=9;++$jsw){
                                       if ($jsw==$parcent){print Fout "*";}else{print Fout "-";}
                                                                 }
                        print Fout "|";print Fout sprintf ("%s [%+e,%+e,%+e,%+e,%+e]\n",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i]);
                   ++$i;}
     close Fout;

                          }