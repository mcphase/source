#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use File::Copy;
if ($#ARGV<1)

{print "program fitcol: simple function fit to data

    use as fitcol coldata prog col parprog1 parprog2 ... []  in filename
               []=[and prog col parprog1 parprog2 ... [and ...]]

          coldata ...... column number of data column to be fitted
          prog    ...... program name, e.g. shiftcol
          col     ...... column to which prog should be applied
          parprog1 ..... parameter of the program prog, which should be fitted
          filename... filename

    in order to fit data in column coldata, prog is run many times on
    the column col with varying parameter set parprog1 ..., the result
    is scaled to fit best the experimental data. If several programs
    are combined with option 'and' then the best linear combination of the
    results is calculated by linear regression to fit coldata.

    Starting values for the parameters parprog are taken from the
    commandline. initial Stepwidths are chosen 10percent of parameter value, or may
    be given by adding them to the parameter with an 's', e.g. 100.3s0.1
    If a parameter should not be fitted and kept fix, add an 'f', e.g. 100.3f

    output: - files can be found in directory results
            - filename.fit is created with fitted function and parameter values

    examples:
    1) to fit a gaussian to column 2 in datafile expdat (with xvalues in column 1)
    with starting values 132.3, 0.5 and 10 for position, fwhm and area, respectively:

    fitcol 2 gausscol 1 132.3 0.5 10 in exp.dat


    2) to do the same fit but with a background create a column 3, fill it
       with constant values and use echo as a fake column manipulation program
      doing nothing.

   newcol 3 exp.dat
   factcol 3 0 exp.dat
   shiftcol 3 1 exp.dat
   fitcol 2 gausscol 1 132.3 0.5 10 and rem 3 in exp.dat

   3) to do the same with a linear background, put into a column 4 the x values

   newcol 4 exp.dat
   factcol 4 0 exp.dat
   addcol 1 4 exp.dat
   fitcol 2 gausscol 1 132.3 0.5 10 and rem 3 and rem 4 in exp.dat

   4) to fit two gaussians with fixed fwhm and stepping in position initially
    only with 0.1

   fitcol 2 gausscol 1 132.3s0.1 0.5f 10 and  gausscol 1 100.3s0.1 0.5f 10 and rem 3 and rem 4 in exp.dat

\n";
exit(0);
}else{print STDERR "#* $0 *\n";}

mkdir ("./results");
$ARGV[0]=~s/x/*/g;$coldata=eval $ARGV[0]; shift @ARGV;

# slice out the programs
@index{@ARGV} = (0..$#ARGV);@progs=@ARGV[0..$index{"in"}-1];$file=$ARGV[$index{"in"}+1];
#slice out different program strings:
@progs=split(/and/,join(' ',@progs));

# first create an intermediate file with additional columns
# for each program to be used
copy($file,"./range0.out");
# create columns corresponding to the number of programs
# and insert them after data column, all filled with zeroes
system("newcol ".($coldata+1)." range0.out");
system("factcol ".($coldata+1)." 0 range0.out");
if($#progs>0){system("newcols ".($coldata+1)." ".($#progs)." range0.out");}
# fill these columns with values from progs operating columns
$i=$coldata;
foreach(@progs)
{@prog=split;
if($prog[1]>$coldata){$prog[1]+=$#progs+1;}
$i++;system("addcol ".$prog[1]." $i range0.out");
 }
# now range0.out is prepared

# for linux prepare calcsta instead of calcsta.bat
open (Fout, ">calcsta");
print Fout << "EOF";
shopt -s expand_aliases
alias del='rm'
alias copy='cp'
alias pause='read'
alias cls='kill'
alias call=''
alias start='. ./do'
alias taskkill='#'
alias wait='sleep'
EOF
print Fout "./calcsta.bat \$*";
close Fout;

#prepare calcsta.bat.forfit
open (Fout, ">calcsta.bat.forfit");
print Fout << "EOF";
copy range0.out range1.out
EOF
# some statements to execute the programs in the fit
$i=1;
foreach(@progs)
{@prog=split;
 $program=$prog[0];shift @prog;shift @prog;
 # here create parameter table
 $j=0;
 foreach(@prog)
  {$old=$prog[$j];
   # check if stepwidth is given
   $step=$old/10;
   if (/s/){($old)=/(.*)s/;$step=/.*s(.*)/;}
   # check if parameter needs to be varied
   if (/f/){($prog[$j])=/(.*)f/;}
      else {$prog[$j]="par".$i.$program.$j."[".$old.",-1e20,+1e20,0,$step]";}
  ++$j;
  }
 $pars=join(' ',@prog);
 print Fout "call ".$program." ".($i+$coldata)." ".$pars." range1.out\n";
 ++$i;
}

print Fout "call linreg $coldata ".($#progs+1)." range1.out\n";
print Fout << "EOF";
copy range1.out $file.fit
EOF
close Fout;

   system ("echo e | perl %MCPHASE_DIR%/bin/simannfit.pl 1e-9 -s 1 calcsta.bat ");
$i=0;
foreach(@progs)
{@prog=split;
if($prog[1]>$coldata){$prog[1]+=$#progs+2;}
$i++;system ("display ".$prog[1]." $coldata $file.fit ".$prog[1]." ".($coldata+$#progs+2)." $file.fit");
 }
   system ("start /B java simannfitstatus");
   system ("echo e | perl %MCPHASE_DIR%/bin//simannfit.pl 1e-9 calcsta.bat ");
   system ("calcsta.bat > calcsta.out");

open (Fout,">>results/simannfit.status");
$i=$coldata;
foreach(@progs)
{@prog=split;
if($prog[1]>$coldata){$prog[1]+=$#progs+1;}
$i++;($a)=extract("a$i","calcsta.out");
print Fout "#! to sum fitcurve from data columns use for column $i in $file.fit coefficient a$i=$a\n";
 }
close Fout;

   system("appendfile $file.fit results/simannfit.status");
 
   copy ("calcsta.bat","results/calcsta.bat");mydel("calcsta.bat");
   copy ("calcsta","results/calcsta");mydel("calcsta");
   copy ("calcsta.bat.forfit","results/calcsta.bat.forfit");mydel("calcsta.bat.forfit");
   copy ("calcsta.bat.bak","results/calcsta.bat.bak");mydel("calcsta.bat.bak");
   copy ("calcsta.out","results/calcsta.out");mydel("calcsta.out");
   copy ("range0.out","results/range0.out");mydel("range0.out");mydel("range1.out");

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


sub mydel  { my ($file1)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;
                               return system("del ".$file1);
                              }
                 else
                              {return system("rm ".$file1);
                              }

           }

