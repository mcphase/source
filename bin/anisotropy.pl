#!/usr/bin/perl
#use Getopt::Long;
#use Math::Trig;
use PDL;
#use PDL::Slatec;
#use Switch;
use File::Copy;

usage() if $#ARGV<2;
$PI=3.141592654;
print STDOUT << "EOF";
#*******************************************************
#               program anisotropy
#        M Rotter 2013 - calculating magnetic anistropy 
#*******************************************************
EOF

system("anisotropyit @ARGV");
print STDOUT << "EOF";
#*******************************************************
#you may now view an anisotropy plot by command
#
#    display 8 13 results/anisotropy.out
#
#              end of program anisotropy
# Reference: M. Rotter et al. PRB 68 (2003) 144418
#*******************************************************
EOF
system("display 8 13 results/anisotropy.out");
exit;

sub usage() {

  print STDERR << "EOF";

    anistropy: program to calculate the magnetic anisotropy by doing a mcphas or
               single ion calculation for different external magnetic field
               directions and evaluating the expectation value of the magnetic 
               moment.

    usage: anisotropy -h
           anisotropy  T H xn yn zn nofsteps [-r sipffilename Hxc1 Hxc2 ... Hxcnofcomponents]
           anisotropy  T H -p nofthetasteps [-r sipffilename Hxc1 Hxc2 ... Hxcnofcomponents]

     -h           : this (help) message
      T           : temperature in Kelvin
      H           : absolute value of the external magnetic field (T)
      xn,yn,zn    : direction normal to plane, in which the anisotropy
                    should be calculated ... e.g. if you want to
                    calculate the anisotropy in the xy plane, then
                    enter xn yn zn = 0 0 1
      nofsteps    : number of steps to be calculated 
     -p           : calculate polycrystal average
      nofthetastepssteps    : number of theta steps to be calculated for polycrystal average

    option:
    -r sipffilename: filename of single ion parameter file
                      Hxc1,Hxc2,... are the exchange field components (meV)
                     (exchange field is kept constant, external magnetic
                     field is rotated in the anisotropy calculation)

    output files:

    ./results/anisotropy.out  contains anisotropy information

EOF

  exit;

}
# **********************************************************************************************
sub calc()
{$AARGV[0]=~s/exp/essp/g;$AARGV[0]=~s/x/*/g;$AARGV[0]=~s/essp/exp/g;$T=eval $AARGV[0];shift @AARGV;$Tm=$T+0.1;
 $AARGV[0]=~s/exp/essp/g;$AARGV[0]=~s/x/*/g;$AARGV[0]=~s/essp/exp/g;$H=eval $AARGV[0];shift @AARGV;$Hm=$H;
 $AARGV[0]=~s/exp/essp/g;$AARGV[0]=~s/x/*/g;$AARGV[0]=~s/essp/exp/g;
 $AARGV[1]=~s/exp/essp/g;$AARGV[1]=~s/x/*/g;$AARGV[1]=~s/essp/exp/g;
 $AARGV[2]=~s/exp/essp/g;$AARGV[2]=~s/x/*/g;$AARGV[2]=~s/essp/exp/g;
 $direction= pdl[eval $AARGV[0],eval $AARGV[1],eval $AARGV[2]];shift @AARGV;shift @AARGV;shift @AARGV;
 $AARGV[0]=~s/exp/essp/g;$AARGV[0]=~s/x/*/g;$AARGV[0]=~s/essp/exp/g;$nstp=eval $AARGV[0];shift @AARGV;
  if ($nstp<1) {die "Error program anisotropy: nofsteps <1\n";}
 $module=$AARGV[0];shift @AARGV;
 for ($module)
  {  if    (/ic1ion/)  {$sipffile=$AARGV[0];shift @AARGV;}
     elsif (/icf1ion/) {$sipffile=$AARGV[0];shift @AARGV;}
     elsif (/so1ion/)  {shift @AARGV;$sipffile=$AARGV[0];shift @AARGV;$cfpartype=$AARGV[0]; shift @AARGV;}
     else {die "ERROR program anisotropy: module $module not implemented\n";}
  }

 $direction/=sqrt(inner($direction,$direction));
# now get $r1 and $r2 which are the basis vectors of the plane of rotation
 $x=pdl [1,0,0];
if (abs(inner($direction,$x))>0.95){$x=pdl[0,1,0];}
 $r1=$x - $direction * inner($direction,$x);
 $r1/=sqrt(inner($r1,$r1));
 $r2=crossp($direction,$r1);
 $r2/=sqrt(inner($r2,$r2));


  $h1x=$r1->at(0);
  $h1y=$r1->at(1);
  $h1z=$r1->at(2);
      unless (open(Fout,">".$fileout)){die "cannot open file $fileout\n";}
    print Fout << "EOF";
# output file of program: anisotropy @command
#! displayxtext=azimuth(deg) in plane perpendicular to $direction
#! displayytext=M||H(mub)
#! displaytitle= Anisotropy plot az=0 corresponds to [$h1x $h1y $h1z]
# phi(deg) theta(deg) T[K] |H|[T] Hx[T] Hy[T] Hz[T] azimuth(deg) |M|[mb] Mx[mb] My[mb] Mz[mb] MparallelH[mb]
EOF
close Fout;
# here comes the loop for different field directions ....
 for($az=0;$az<2*$PI-0.00001;$az+=2*$PI/$nstp)
 {$Hv=$H*(cos($az)*$r1+sin($az)*$r2);
  $Hx=$Hv->at(0);
  $Hy=$Hv->at(1);
  $Hz=$Hv->at(2);
  $phi=$PI/2;if($Hy<0){$phi=3*$PI/2;}
  if($Hx>0.001&$Hy>=0){$phi=atan($Hy/$Hx);}
  if($Hx>0.001&$Hy<0){$phi=atan($Hy/$Hx)+2*$PI;}
  if($Hx<-0.001){$phi=atan($Hy/$Hx)+$PI;}
  $theta=acos($Hz/$H);

  for ($module)
  {  if   (/ic1ion/)  {ic1ion();}
     elsif(/icf1ion/)  {ic1ion();}
     elsif(/so1ion/)  {so1ion();}
     else {die "ERROR program anisotropy: module $module not implemented\n";}
  }
  unless (open(Fout,">>".$fileout)){die "cannot open file $fileout\n";}
  print Fout  sprintf("%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.3f   %6.3f   %6.3f %6.3f %6.3f %6.3f\n",$phi*180/$PI,$theta*180/$PI,$T,$H,$Hx,$Hy,$Hz,$az*180/$PI,$M,$Mx,$My,$Mz,$Mp);
  close Fout;
 }
    unlink("anisotropy.sipf");
}
#***************************************************************************************
sub ic1ion()
{ #copy($sipffile,"anisotropy.sipf");
  open(F,">anisotropy.sipf");
  open (Fin,$sipffile);
  while($line=<Fin>) {next if ($line=~/calcmag/);
                      next if ($line=~/xT/);
                      next if ($line=~/xHa/);
                      next if ($line=~/xHb/);
                      next if ($line=~/xHc/);
                      next if ($line=~/xmin/);
                      next if ($line=~/xstep/);
                      next if ($line=~/xmax/);
                      next if ($line=~/yT/);
                      next if ($line=~/yHa/);
                      next if ($line=~/yHb/);
                      next if ($line=~/yHc/);
                      next if ($line=~/ymin/);
                      next if ($line=~/ystep/);
                      next if ($line=~/ymax/);
                      print F $line;}
  print F << "EOF";

calcmag
 xT   = 1
 xHa  = 0
 xHb  = 0
 xHc  = 0
 xstep= 1
 yT   = 0
 ystep= 1
EOF
print F sprintf(" xmin = %8.3f\n",$T);
print F sprintf(" xmax = %8.3f\n",$Tm);
print F sprintf(" yHa = %8.3f\n",$Hx);
print F sprintf(" yHb = %8.3f\n",$Hy);
print F sprintf(" yHc = %8.3f\n",$Hz);
print F sprintf(" ymin = %8.3f\n",$H);
print F sprintf(" ymax = %8.3f\n",$Hm);
close F;
  system("$module anisotropy.sipf");
  ($M)=extractnumber(1,5,"./results/ic1ion.mag");
  ($Mx)=extractnumber(1,6,"./results/ic1ion.mag");
  ($My)=extractnumber(1,7,"./results/ic1ion.mag");
  ($Mz)=extractnumber(1,8,"./results/ic1ion.mag");
  ($Mp)=extractnumber(1,9,"./results/ic1ion.mag");

}
# **********************************************************************************************
sub so1ion()
{ #copy($sipffile,"anisotropy.sipf");
  open(F,">anisotropy.sipf");
  open (Fin,$sipffile);
  while($line=<Fin>) {next if ($line=~/TEMP\s*=/);
                      next if ($line=~/T\s*=/);
                      next if ($line=~/Bx\s*=/);
                      next if ($line=~/By\s*=/);
                      next if ($line=~/Bz\s*=/);
                      print F $line;}
print F sprintf(" TEMP=%8.3f\n",$T);
print F sprintf(" Bx = %8.3f\n",$Hx);
print F sprintf(" By = %8.3f\n",$Hy);
print F sprintf(" Bz = %8.3f\n",$Hz);
close F;
  system("so1ion -r anisotropy.sipf $cfpartype");
  ($Mx)=extract("mx","./results/so1ion.out");
  ($My)=extract("my","./results/so1ion.out");
  ($Mz)=extract("mz","./results/so1ion.out");
   $M=sqrt($Mx*$Mx+$My*$My+$Mz*$Mz);
   $Mp=($Mx*$Hx+$My*$Hy+$Mz*$Hz)/$H;
}

# **********************************************************************************************
# extracts variable from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# ($standarddeviation)=extract("sta","data.dat");
#
# ... it stores 0.24 in the variable $standarddeviation
#
sub extract {
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^.*$var\s*=/) {($value)=($line=~m|$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);}                                        }
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
# **********************************************************************************************
# extracts number in col from line file
#
# for example somewhere in a file data.dat has 3 columns
#  # some header
#  0 2    3
#  2 0.24 3
#  4 2.3  2
#
# ($number)=extractnumber(3,2,"data.dat");
#
# ... it stores 2.3 in the variable $number
#
sub extractnumber {
             my ($line,$col,$filename)=@_;
             if(open (Fin,$filename))
             {my $i=0;while(($read=<Fin>)&&$i<$line){
                                                if($read=~/^\s*#/) {;}
                                                else{++$i;
                                         if($i==$line){@numbers=split(" ",$read);close Fin;return $numbers[$col-1];}
                                                    }

                                               }
               close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";exit(1);
             }
            }
# **********************************************************************************************