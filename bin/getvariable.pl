#!/usr/bin/perl
use File::Copy;

BEGIN{@ARGV=map{glob($_)}@ARGV}
unless ($#ARGV >0)
{
print STDOUT << "EOF";

 program to get the value of a variable from a file 
  (e.g. somewhere in a file there is 
   a statement T=4.3 and you want to get out this 4.3)

 usage: perl getvariable variablename filename

 output: the variable value is written to stdout and environment variable 
         MCPHASE_GETVARIABLE_VALUE, the name is stored in 
         MCPHASE_GETVARIABLE_NAME
        mind lines starting with # are ignored (unless these start with #!)

EOF
# clean bat files
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");close Fout;
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");close Fout;
exit(1);
}

$name=$ARGV[0];shift @ARGV;
foreach(@ARGV)
{$filename=$_;
($value)=extract($name,$filename);
print "#! getvariable: in file $filename  $name=$value\n";
}
# for setting environment variables
open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
print Fout "set MCPHASE_GETVARIABLE_NAME=$name\n";
print Fout "set MCPHASE_GETVARIABLE_VALUE=$value\n";
close Fout;

open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
print Fout "export MCPHASE_GETVARIABLE_NAME=$name\n";
print Fout "export MCPHASE_GETVARIABLE_VALUE=$value\n";
close Fout;

exit(0);


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
sub extract {$value="not found";
             my ($variable,$filename)=@_;
             $var="\Q$variable\E";
             if(open (Fin,$filename))
             {while($line=<Fin>){
                if($line=~/^.*$var\s*=/) {($value)=($line=~m|$var\s*=\s*([^\s]+)|);}                                        }
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
