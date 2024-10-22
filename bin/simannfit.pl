#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use PDL;
use Cwd;

# use PDL::Slatec;

unless ($#ARGV >0) 
{print STDOUT << "EOF";
  program simannfit used to perform simulating annealing fit / paramter search

  usage: simannfit 10 [options] * 
  
  10 .. initial temperature 
  * .. filename(s) of paramter file(s)
 
  ... there must exist a program calcsta.bat (windows) / calcsta (linux), the output of which
  contains 'sta = 249' - this program is called from this output
  at every iteration step and  the standard sta deviation is minimzed
  simannfit gives as a parameter a maximum number for sta - if during the
  calculation of sta in calcsta.bat this number is exceeded calcsta.bat can exit (saves time)
  
  
  OPTIONS for ...
  
  ... STEPWIDTHS
    -w 1.4  ... before starting simannfit, multiply all stepwidths by factor 1.4 
    -r 0.2  ... before starting simannfit, set all stepwidths to parameter 
            range=max-min times 0.2, 
    -f 0.2  ... before starting simannfit, set all stepwidths to parameter value 
            times 0.2, however never smaller than parameterrange/1000 
    -p 20  ... probing parameter space: stepwidths are not decreased during
	    fitting, new set of pars parn are generated and for all sets par in 
	    results/simannfit.p and results/simannfit.n are scanned and
	    distance d=sum_i^N (par_i-parn_i)^2/stepwidht_i^2 is calculated. If d>N=nofparameters
	    for all sets par, then sta is calculated. if  sta<=sta of initial parameters
	    then parn is appended to the list in results/simannfit.p
	    up to 20 parameters are appended to this list, after that the program stops.
	    the program also stops if step ratio reaches a value > 11 (step ratio 
	    is a factor applied to all (maximum) stepwidths before generating a step in the random walk
	    It is decreased each successful step in the random walk and increased to find an
	    acceptable step (distance criteria) to try
	    thus, to probe a large region of parameter space do not set statistical
	    temperature too low to enable efficient stepping. On the other hand to
	    probe only the vicinity of a optimum in every direction of parameter
	    space a small statistical inital temperature is necessary
	    - in this way a series of equally good solutions can be explored.

 ... PARALLEL PROCESSING
    -s0 filename ... instead of results/simannfit.n use filename to store/read 
    -s1 filename ... instead of results/simannfit.p use filename to store/read 
    -i 23 ... for storing in results/simannfit.* use index .23 instead of .0 as suffix 

 ... LOGGING 
     -n 50  ... specifies that every 50 steps the parameters should be
            stored in file results/simannfit.n (appending existing file)
     -jpglog 3.5e-1 file.jpg    ... works only if option -n is present:
             if sta is less then 3.5e-1, then the image file.jpg 
            (which should be created by calcsta) is copied to results/parsetnr.jpg
            a html tag is added to results/simannfit.n or results/simannfit.p if 
	    stored
    -log 1.3 batchfile.bat   ... works only if option -n is present:
             if sta is less then 1.3, then execute the file 
             batchfile, with sta as argument, which can be adressed in the batch
	     file with \$1 (linux) or  \%1 (windows)
    -h       Histograms are stored in results/par*.hst for review of the variation 
             of parameters during the run of simannfit. 
    -d       store percentage of different contributions to sta in file results/simannfit.dst 

... TERMINATION 
   -l 0.2 ... sets limit - program will end if sta < 0.2
   -t 100 ... sets time limit until program end to 100 seconds
   -s 132 ... gives maximal number of iteration steps 132 to be done
   -c     ... continue at end of program - do not ask for pressing enter 

EOF
#print " <Press enter to close>";$in=<STDIN>;
 exit 0;}else{print STDERR "#* $0 *\n";}


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

 @parnam=();@par=();@parmin=();@parminini=();@parmax=();@parmaxini=();@parerr=();@parstp=();@parav=();@thisparstp=();
				 @parhisto=();@parhistostp=();@perlhistostart=();$hh=0;
  $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$stattemp=eval $ARGV[0]; shift @ARGV;
  $starttime=time;$maxtim=1e10;$maxstep=1e24;$tablestep=0;$stepset=0;$limsta=-1e100;$tableoffset=0;
  $options=1;$probe=0;$stepfact=1;$cont=0;$jpglog=0;$index=0;$log=0;$hist=0;$dist=0;
  $s0file="results/simannfit.n";
  $s1file="results/simannfit.p";
  while($options==1)
  {$options=0;
  if ($ARGV[0] eq '-l') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$limsta=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-t') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$maxtim=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-s') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$maxstep=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-n') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$tablestep=eval $ARGV[0]; shift @ARGV;}  
  if ($ARGV[0] eq '-p') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$probe=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-w') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$stepfact=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-r') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$rangfact=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-f') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$stepset=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-i') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$index=eval $ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-s0') {$options=1;shift @ARGV;$s0file=$ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-s1') {$options=1;shift @ARGV;$s1file=$ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-jpglog') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$jpglog=eval $ARGV[0]; shift @ARGV;$jpgimagefile=$ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-log') {$options=1;shift @ARGV; $ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$log=eval $ARGV[0]; shift @ARGV;$logbatchfile=$ARGV[0]; shift @ARGV;}
  if ($ARGV[0] eq '-c') {$options=1;shift @ARGV; $cont=1;}
  if ($ARGV[0] eq '-h') {$options=1;shift @ARGV; $hist=1;}
  if ($ARGV[0] eq '-d') {$options=1;shift @ARGV; $dist=1;}
  }

 while(!open(Fout,">results/simannfit.status")){print "Error opening file results/simannfit.status\n";<STDIN>;}
   print Fout "simannfit running in ".cwd()."\n";
   print Fout "parameter[value,      min,           max,           variation,     stepwidth]\n";
  foreach (@ARGV)
 {$file=$_; if(mycopy ($file,$file.".bak")){print "\n warning copying $file not possible - press enter to continue\n";<stdin>;}
   unless (open (Fin, $file.".forfit")){die "\n error:unable to open $file.forfit\n";}   
    while($line=<Fin>)
 {while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/) {++$#par;#load another parameter
				 ($parname)=($line=~m/(?:#!|[^#])*?\b(par\w+)\s*\Q[\E/);
                                 foreach(@parnam){if ($_ eq $parname){print "ERROR simannfit: parameter $parname occurs more than one time in input files\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}}
                                 $parnam[$#par]=$parname;
				 ($par[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*([^,]+)/);
				 ($parmin[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*([^,]+)/);$parminini[$#par]=$parmin[$#par];
				 ($parmax[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);$parmaxini[$#par]=$parmax[$#par];
				 ($parerr[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parstp[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^\Q]\E]+)/);
                         if($stepset>0){$parstp[$#par]=abs($par[$#par])*$stepset;
                                        if($parstp[$#par]<($parmax[$#par]-$parmin[$#par])/1000)
                                        {$parstp[$#par]=($parmax[$#par]-$parmin[$#par])/1000;}
                                        }
                         if($rangfact){$parstp[$#par]=$rangfact*($parmax[$#par]-$parmin[$#par]);}
                         $parstp[$#par]=$stepfact*abs($parstp[$#par]);
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

  if($tablestep!=0){writeini($s0file);}
  if($probe!=0){writeini($s1file);}


print "initialize parameter storage\n";
$parstore = zeroes $#par+3,$#par+1;
$deltastore =  PDL->nullcreate(0);
$nof_calcsta_calls=0;
#*******************************************************************************
# fitting loop
 print ($#par+1);print " parameters found - testing calculation of sta\n";
$rnd=1;$stasave=1e20;
 ($sta)=sta();$stps=1;$noofupdates=0;$stepnumber=0;$stastart=$sta;
  if($maxstep==1){open(Fout,">results/simannfit.status");
                  print Fout " maxstep=1 ... simannfit stopped after initial run of calcsta\n";
                  print Fout ($#ssta+1)." contributions to sta found in output of calcsta ...\n";
                  print Fout "sta=$sta\n";
                  close Fout;
                  print " maxstep=1 ... simannfit stopped after initial run of calcsta\n";
                  print "sta=$sta\n";
                  print " <Press enter to close>";if($cont==0){$in=<STDIN>;}
                 exit 0;}
if($tablestep!=0){ write_set(">>$s0file");}
if($sta>0)
{print "starting fit\n";
 $SIG{INT} = \&catch_zap;
 while($sta>0)
 {  $stasave=$sta;
 # modify parameters
 print "\n ...next fitting loop ...\n";
 @parsav=@par;$dmin=0;
 if($probe>0){$stps=1;} # if probing is desired initialize step to 1 and make it larger if no set is found ...
 while($dmin<$#par+1&&$stps<11.0)
 {$i=0;
  foreach(@par){$rnd=rand;$thisparstp[$i]=($rnd-0.5)*$parstp[$i]*$stps;$par[$i]+=$thisparstp[$i];
               if ($par[$i]<$parmin[$i]) {$thisparstp[$i]=$parmin[$i]-$par[$i];$par[$i]=$parmin[$i];}
	       if ($par[$i]>$parmax[$i]) {$thisparstp[$i]=$parmax[$i]-$par[$i];$par[$i]=$parmax[$i];}
	       ++$i;}
   $dmin=1e10; 
   if($probe>0)# check distance to other parameter sets
               {if(open(Fin1,$s1file))
                {while(($line=<Fin1>)&&($dmin>$#par+1))
                    {unless ($line=~/^\s*#/)
                        {$line=~s/D/E/g;@numbers=split(" ",$line);
                          $d=0;$i=0;foreach(@par){$dd=($par[$i]-$numbers[$i+1])/$parstp[$i];$d+=$dd*$dd; ++$i;}
                          if($d<$dmin){$dmin=$d;}
                         }
                     }   
                close Fin1;
                } 
               if(open(Fin1,$s0file))
                {while(($line=<Fin1>)&&($dmin>$#par+1))
                    {unless ($line=~/^\s*#/)
                        {$line=~s/D/E/g;@numbers=split(" ",$line);
                          $d=0;$i=0;foreach(@par){$dd=($par[$i]-$numbers[$i+1])/$parstp[$i];$d+=$dd*$dd; ++$i;}
                          if($d<$dmin){$dmin=$d;}
                         }
                     }   
                close Fin1;
                } 
               }
   if($dmin<$#par+1){$stps*=1.001;if($stps<2){@par=@parsav;}}
  }
if($stps<11){
 print " .. calculating sta ..\n";

 $rnd=rand;
   ($sta)=sta(); # CALCULATE sta !!!!
   ++$stepnumber;
   
   if($tablestep!=0&&$stepnumber%$tablestep==0){ write_set(">>$s0file");}
   if($probe>0&&$sta<=$stastart){--$probe;print "#dmin=$dmin>Npar=".($#par+1)." parset stored in $s1file - $probe other sets to be found, continuing ...\n";write_set(">>$s1file");                      
                last if ($probe==0);
                                }
   print " ...  current sta=$sta, statistical T=$stattemp, step ratio=$stps\nsta of stored parameters=$stasave\n";
   open(Fin,"./results/simannfit.status");$line=<Fin>;
    if($sta==0){@parsav=@par;} # if sta really got zero... store parameters as parsav 
    if ($line=~/exiting simannfit/){$sta=0;close Fin;}
    else
    {close Fin;
     read_write_statusfile();
    }
          } else {$sta=0;print "\n stepsize = $stps > 11.0 - stopping fit\n";}
 if (time-$starttime>$maxtim){$sta=0;print "\n maximum time for fitting reached - stopping fit\n";}
 if ($stepnumber>$maxstep){$sta=0;print "\n maximum step number for fitting reached - stopping fit\n";}
 if ($sta==0) {#recover old pars
      @par=@parsav; 
      }
 last if ($sta<$limsta);
 if ($sta>$stasave)
  {if($rnd>exp(-($sta-$stasave)/$stattemp))
     {#recover old pars and adapt parstep to step not so big in this direction
      @par=@parsav;$i=0;if($probe==0){foreach(@parstp){if($parstp[$i]>$parhistostp[$i]/1000){$parstp[$i]-=0.1*abs($thisparstp[$i]);}++$i;}}
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
   if($probe==0){if($parstp[$i]<($parmax[$i]-$parmin[$i])/2){$parstp[$i]+=0.1*abs($thisparstp[$i]);}
                 else{$parstp[$i]=($parmax[$i]-$parmin[$i])/2;}
                } # adapt parstp to be more bold in this direction
   $hx=int(($p-$parminin[$i])/$parhistostp[$i]);
   ++$parhisto[($hx+$perlhistostart[$i])];
   if($hist>0){open(Fout,">./results/".$parnam[$i].".hst");
              print Fout "#{Histogram of parameter ".$parnam[$i]."\n# value vs. number of  occurrences in good solutions (sta decreased)}\n";
              for($hx=0;$hx<=int(($parmaxini[$i]-$parminini[$i])/$parhistostp[$i])+1;++$hx)
              {print Fout (($hx+0.5)*$parhistostp[$i]+$parminini[$i])."   ".($parhisto[($hx+$perlhistostart[$i])])."\n";
              } close Fout;
              }

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
     open(Fout,">results/simannfit.status");print Fout cwd()." ... simannfit stopped\n";
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
     print Fout "Covariance matrix( may be not successfull because last n steps of simulated annealing may be not necessarily\n orthogonal in parameter space - if this happens  (you will get error=0) restart and try again):\n";
$Fij=$delta x inv($V);
 $FtF=$Fij->xchg(0,1) x $Fij;
 $cov=$sta*inv($FtF);  # multiply by sta=chi2 in order to get covariance matrix
     print $cov; print Fout $cov;
                $i=0;foreach(@par){if($cov->at($i,$i)>0){print $parnam[$i]." error=".(sqrt($cov->at($i,$i)))."\n";
                                   print Fout $parnam[$i]." error=".(sqrt($cov->at($i,$i)))."\n";}
                                   else {print $parnam[$i]." covariance matrix not successful error=0\n";
                                         print Fout $parnam[$i]." covariance matrix not successful error=0\n";
                                         }
                                   ++$i;}
     }
     close Fout;if($tablestep!=0){++$stepnumber;write_set(">>$s0file");}
print " <Press enter to close>";if($cont==0){$in=<STDIN>;}
exit 0;
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
                                        {# check if parmin or parmax have been changed by user and update
                                 ($parminnew)=($line=~m/(?:#!|[^#])*?\b$pnam\s*\Q[\E\s*[^,]+\s*,\s*([^,]+)/);
                                 ($parmaxnew)=($line=~m/(?:#!|[^#])*?\b$pnam\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
                                 unless($parminnew>=$parmaxnew||$par[$i]>$parmaxnew||$par[$i]<$parminnew){  
                                  if($parminnew!=$parmin[$i]){print "parmin of parameter $pnam changed from ".$parmin[$i]." to $parminnew \n";$parmin[$i]=$parminnew;}
				  if($parmaxnew!=$parmax[$i]){print "parmax of parameter $pnam changed from ".$parmax[$i]." to $parmaxnew \n";$parmax[$i]=$parmaxnew;}
								}
				          $line=~s|$pnam\s*\Q[\E[^\Q]\E]*\Q]\E|$par[$i]|;
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
                                    $line=~s|function\s*\Q[\E.*?\Q]\E|$data|;
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
 if($dist>0){open(Fout1, ">results/simannfit.dst");
            print Fout1 "# ".($#ssta+1). " contributions to sta\n";
            print Fout1 "# number  vs percentage vs sta \n";$i6=0;
            foreach(@ssta){++$i6;
                           print Fout1 sprintf("%i %6.2f %g\n",$i6,(100*$_/(($#ssta+1)*$s2)),$_);
                          }
            close Fout1;
            }
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
     print Fout "#! working directory: ".cwd()."\n#! NS=".($#ssta+1)." contributions to sta found in output of calcsta ...\n";
     if($chisquared){print Fout "#! Current sta=chi2=$sta (=sum deviations^2/(".($#ssta+1)."*experrors^2))\n#! sta_of_stored_parameters=$stasave  initial_sta=$stastart\n";}
              else {print Fout "#! Current     sta=variance=s2=$s2 (=sum deviations^2/".($#ssta+1).")  \n#! sta_of_stored_parameters=$stasave  initial_sta=$stastart\n";}
     print Fout "#----------------------------------------------------------------------------------------\n";
     print Fout "#! Statistical_Temperature=$stattemp      Step_Ratio=$stps\n";
     print Fout "#----------------------------------------------------------------------------------------\n";
     $est=sprintf("%6.2f",(time-$starttime)/3600);$maxtimest=sprintf("%6.2f",($maxtim)/3600);
     print Fout "#! Time_since_start_of_simannfit=$est hours (tlimit=$maxtimest), N=$stepnumber steps (slimit=$maxstep)\n";
     if($probe>0){print Fout "#! P=$probe sets of different parameters to be found for $s1file\n";}
     print Fout "#----------------------------------------------------------------------------------------\n";
     print Fout "# parameter[value,      min,           max,           variation,     stepwidth]\n";
     foreach(@par){$parcent=int(10*($par[$i]-$parmin[$i])/(1e-10+$parmax[$i]-$parmin[$i]));
                        print Fout "#|";for($jsw=0;$jsw<=9;++$jsw){
                                       if ($jsw==$parcent){print Fout "*";}else{print Fout "-";}
                                                                 }
                        print Fout "|";print Fout sprintf ("%s [%+e,%+e,%+e,%+e,%+e]\n",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i]);
                   ++$i;}
     close Fout;

                          }

sub writeini()
   {my ($filename)=@_;
    #local *FH; removed 28.1.2017 because FH was not returned correctly 
    if(open(FH,$filename)){print "Appending table to $filename ...\n";
       while($line=<FH>){unless ($line=~/^\s*#/)
                                {$line=~s/D/E/g;@numbers=split(/ |\./,$line);
                                 if($index==$numbers[1]){
                                 if($numbers[0]>=$tableoffset){$tableoffset=$numbers[0]+1;}
                                                   }
                                }
                      
                               }
    close FH;} 
    while(!open(FH,">>$filename")){print "Error opening file $filename\n";<STDIN>;}
   print "storing points in file $filename\n";
   print FH "#IterationNr ";foreach(@parnam){print FH $_." ";}print FH "sta variance chisquared\n";
   close FH;
   }

sub write_set()
{ my ($filename)=@_;
  unless(open(FH,$filename)){die "Error openening $filename\n";}
     my $dd=sprintf("%i.%i ",$stepnumber+$tableoffset,$index);print FH $dd;
        my $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print FH $dd;++$ii} 
        print FH $sta." ".$s2." ".$chisquared."\n";
        if($sta<$jpglog){$dfile=sprintf("%i.%i.jpg",$stepnumber+$tableoffset,$index);                         
                         mycopy($jpgimagefile,"./results/".$dfile);
                         print FH '#<img src="'.$dfile.'">'."\n";
                         }
  close FH;
        if($sta<$log){system("$logbatchfile $sta");
                         }
}
