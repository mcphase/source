#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}


if ($#ARGV<1)

{print "program fitfermi: calculates Fermi function fit to data in file 

    use as fitfermi T EF fwhm min max filename

          T       ...... Temperature in K
          EF      ...... initial value of Fermi Energy (eV)
          fwhm    ...... initial value for resolution  (eV)
          min,max ...... energy range of fit
	  filename... filename (column 1 is energy in eV and col 2 intensity)

    fermifunction is defined as
          f(E)=b+l(E-EF)+[d+k(E-EF)]/(exp((E-EF)/kT)+1)
          f(E) is convoluted with a Gaussian function of given fwhm
               and the  result is compared to exp data

    output: - files can be found in directiry results
            - filename.fit is created with fitted function and parameter values


\n";
exit(0);
}

$T=$ARGV[0]; shift @ARGV;
$EF=$ARGV[0]; shift @ARGV;
$fwhm=$ARGV[0]; shift @ARGV;
$min=$ARGV[0]; shift @ARGV;
$max=$ARGV[0]; shift @ARGV;
system("mkdir results");
 foreach (@ARGV)
  {$file=$_;
    if ($^O=~/MSWin/){
   open (Fout, ">calcsta.bat.forfit");
print Fout << "EOF";
copy $file exp.dat

call range -d 1 function[$min-4*parfwhm] function[$max+4*parfwhm] exp.dat
rem shift energy by EF
call shiftc 1 function[-parEF] exp.dat
call newcols 1 4 exp.dat

rem col 2 - const bkg
call factcol 2 0 exp.dat
call shiftcol 2 1 exp.dat
rem col 3 - linear bkg - ok
rem col 4 fermi function
rem ...
call fermicol 4 300 exp.dat
rem col 5 (E-EF)*Fermifunction
call multcol 4 5 exp.dat

rem create resolution function from fwhm and convolute col 2-5 with that
call gauss parfwhm [$fwhm,1.000000e-002,1.000000e+000,2.828650e-003,$fwhm/10] function[parfwhm/10] function[-parfwhm*1.6] function[parfwhm*1.6] > results\\gauss.dat
EOF

$gausspoints=33; # has to be consistent with the number of points in gauss function above

print Fout << "EOF";
rem here we convolute storing the result of the convolution
rem in columns 7 8 9 10
call delcols 7 100 exp.dat
call convolute 1 2 exp.dat 1 2 results/gauss.dat 1 2 exp.dat  > exp2.cvt
call delcol 8 exp2.cvt
call convolute 1 3 exp2.cvt 1 2 results/gauss.dat 1 3 exp2.cvt  > exp.dat
call delcol 9 exp.dat
call convolute 1 4 exp.dat 1 2 results/gauss.dat 1 4 exp.dat  > exp2.cvt
call delcol 10 exp2.cvt
call convolute 1 5 exp2.cvt 1 2 results/gauss.dat 1 5 exp2.cvt  > exp.dat
call delcol 11 exp.dat
del exp2.cvt

rem  before linreg we define the range of fit min-EF to max-EF
call range -d 1 function[$min-parEF] function[$max-parEF]  exp.dat

call linreg  6 4 exp.dat

rem now column 11 should contain the calculated curve

call shiftcol 1 parEF [$EF,$min,$max,2.270382e-003,0.1] exp.dat
copy exp.dat $file.fit
EOF
close Fout;}
else # for linux system
{  open (Fout, ">calcsta.forfit");
print Fout << "EOF";
cp $file exp.dat

range -d 1 function[$min-4*parfwhm] function[$max+4*parfwhm] exp.dat
rem shift energy by EF
shiftc 1 function[-parEF] exp.dat
newcols 1 4 exp.dat

rem col 2 - const bkg
 factcol 2 0 exp.dat
 shiftcol 2 1 exp.dat
rem col 3 - linear bkg - ok
rem col 4 fermi function
rem ...
 fermicol 4 300 exp.dat
rem col 5 (E-EF)*Fermifunction
 multcol 4 5 exp.dat

rem create resolution function from fwhm and convolute col 2-5 with that
gauss parfwhm [$fwhm,1.000000e-002,1.000000e+000,2.828650e-003,$fwhm/10] function[parfwhm/10] function[-parfwhm*1.6] function[parfwhm*1.6] > results/gauss.dat
EOF

$gausspoints=33; # has to be consistent with the number of points in gauss function above

print Fout << "EOF";
rem here we convolute storing the result of the convolution
rem in columns 7 8 9 10
delcols 7 100 exp.dat
convolute 1 2 exp.dat 1 2 results/gauss.dat 1 2 exp.dat  > exp2.cvt
delcol 8 exp2.cvt
convolute 1 3 exp2.cvt 1 2 results/gauss.dat 1 3 exp2.cvt  > exp.dat
delcol 9 exp.dat
convolute 1 4 exp.dat 1 2 results/gauss.dat 1 4 exp.dat  > exp2.cvt
delcol 10 exp2.cvt
convolute 1 5 exp2.cvt 1 2 results/gauss.dat 1 5 exp2.cvt  > exp.dat
delcol 11 exp.dat
rm exp2.cvt

rem  before linreg we define the range of fit min-EF to max-EF
range -d 1 function[$min-parEF] function[$max-parEF]  exp.dat

 linreg  6 4 exp.dat

rem now column 11 should contain the calculated curve

shiftcol 1 parEF [$EF,$min,$max,2.270382e-003,0.1] exp.dat
cp exp.dat $file.fit
EOF
close Fout;
}
   if ($^O=~/MSWin/){
   system ("echo e | perl %MCPHASE_DIR%\\bin\\simannfit.pl 1e-9 -s 1 calcsta.bat ");
   system ("display 1 6 $file.fit 1 11 $file.fit");
   system ("start /B java simannfitstatus");
   system ("echo e | perl %MCPHASE_DIR%\\bin\\simannfit.pl 1e-9 calcsta.bat ");
   system("factcol 2 $gausspoints exp.dat");
   system("factcol 3 $gausspoints exp.dat");
   system("factcol 4 $gausspoints exp.dat");
   system("factcol 5 $gausspoints exp.dat"); 
   system ("calcsta.bat > results\\calcsta.out");
   system('substitute "# kinetic Energy(eV) vs Intensity_normalized_normalized_to_frames_nofscans vs Error_normalized_to_frames_nofscans vs Intensity vs Error" "# E(eV) b l(E-EF) d/(exp((E-EF)/kT)+1) k(E-EF)]/(exp((E-EF)/kT)+1) Intensity_normalized_normalized_to_frames_nofscans b_conv l(E-EF)_conv d/(exp((E-EF)/kT)+1)_conv k(E-EF)]/(exp((E-EF)/kT)+1)_conv fitted_Intensity" '.$file.'.fit');
   system("comment 1 100 results\\simannfit.status");
   system("appendfile $file.fit results\\simannfit.status");
 
   ($a7)=extract("a7","results/calcsta.out");
   ($a8)=extract("a8","results/calcsta.out");
   ($a9)=extract("a9","results/calcsta.out");
   ($a10)=extract("a10","results/calcsta.out");
   $a7*=$gausspoints;
   $a8*=$gausspoints;
   $a9*=$gausspoints;
   $a10*=$gausspoints;
   system("copy calcsta.bat results");system("del calcsta.bat");
   system("del calcsta.bat.forfit");
   system("del calcsta.bat.bak");
#   system("del calcsta.bat.forfit.bak");
   
#   system("del exp.dat");
   print "\n#! fitted parameters for file=$file\n";
   system("type results\\simannfit.status");
   } else
   {
   system ("echo e | perl $MCPHASE_DIR/bin/simannfit.pl 1e-9 -s 1 calcsta ");
   system ("display 1 6 $file.fit 1 11 $file.fit");
   system ("java simannfitstatus &");
   system ("echo e | perl $MCPHASE_DIR/bin/simannfit.pl 1e-9 calcsta ");
   system("factcol 2 $gausspoints exp.dat");
   system("factcol 3 $gausspoints exp.dat");
   system("factcol 4 $gausspoints exp.dat");
   system("factcol 5 $gausspoints exp.dat");
   system ("./calcsta > results/calcsta.out");
   system('substitute "# kinetic Energy(eV) vs Intensity_normalized_normalized_to_frames_nofscans vs Error_normalized_to_frames_nofscans vs Intensity vs Error" "# E(eV) b l(E-EF) d/(exp((E-EF)/kT)+1) k(E-EF)]/(exp((E-EF)/kT)+1) Intensity_normalized_normalized_to_frames_nofscans b_conv l(E-EF)_conv d/(exp((E-EF)/kT)+1)_conv k(E-EF)]/(exp((E-EF)/kT)+1)_conv fitted_Intensity" '.$file.'.fit');
   system("comment 1 100 results/simannfit.status");
   system("appendfile $file.fit results/simannfit.status");
 
   ($a7)=extract("a7","results/calcsta.out");
   ($a8)=extract("a8","results/calcsta.out");
   ($a9)=extract("a9","results/calcsta.out");
   ($a10)=extract("a10","results/calcsta.out");
   $a7*=$gausspoints;
   $a8*=$gausspoints;
   $a9*=$gausspoints;
   $a10*=$gausspoints;
   system("copy calcsta results");system("rm calcsta");
   system("del calcsta.forfit");
   system("del calcsta.bak");
#   system("del calcsta.forfit.bak");
   
#   system("del exp.dat");
   print "\n#! fitted parameters for file=$file\n";
   system("cat results/simannfit.status");
}
   
      
   
   print "#! constant background b=$a7\n";
   print "#! linear background l=$a8 eV^-1\n";
   print "#! DOS at EF d=$a9 \n";
   print "#! slope of DOS at EF k=$a10 eV^-1\n";
open (Fout,">>$file.fit");
print Fout << "EOF";
#! constant background b=$a7
#! linear background l=$a8 eV^-1
#! DOS at EF d=$a9 
#! slope of DOS at EF k=$a10 eV^-1
#   fermifunction is defined as
#          f(E)=b+l(E-EF)+[d+k(E-EF)]/(exp((E-EF)/kT)+1)
#          f(E) is convoluted with a Gaussian function of given fwhm
#               and the  result is compared to exp data
EOF

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


