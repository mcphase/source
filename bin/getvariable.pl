#!/usr/bin/perl
use File::Copy;
use Getopt::Long;

BEGIN{@ARGV=map{glob($_)}@ARGV}
unless ($#ARGV >0)
{
print STDERR << "EOF";
#
# program to get the value of a variable from a file 
#  (e.g. somewhere in a file there is 
#   a statement T=4.3 and you want to get out this 4.3)
#
# usage: perl getvariable [options] variablename filename
#
# options: -c 24.13   compare the value with 24.13+-0.01 (error corresponds to last 
#                    digit, and exit with failure message if extracted y-value is
#                    not corresponding
#
# output: the variable value is written to stdout and environment variable 
#         MCPHASE_GETVARIABLE_VALUE, the name is stored in 
#         MCPHASE_GETVARIABLE_NAME
#        mind lines starting with # are ignored (unless these start with #!)
#
EOF
# clean bat files
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");close Fout;
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");close Fout;
exit(1);
}else{print STDERR "#* $0 *\n";}

GetOptions("c=s"=>\$compare);


$name=$ARGV[0];shift @ARGV;
foreach(@ARGV)
{$filename=$_;
($value)=extract($name,$filename);
print STDERR "#! getvariable: in file $filename  $name=$value\n";
}
# for setting environment variables
#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat.bat");
#print Fout "set MCPHASE_GETVARIABLE_NAME=$name\n";
#print Fout "set MCPHASE_GETVARIABLE_VALUE=$value\n";
#close Fout;

#open (Fout,">$ENV{'MCPHASE_DIR'}/bin/bat");
#print Fout "export MCPHASE_GETVARIABLE_NAME=$name\n";
#print Fout "export MCPHASE_GETVARIABLE_VALUE=$value\n";
#close Fout;

 if ($^O=~/MSWin/){
print "set MCPHASE_GETVARIABLE_NAME=$name\n";
print "set MCPHASE_GETVARIABLE_VALUE=$value\n";
                  }
                 else
                  {
print "export MCPHASE_GETVARIABLE_NAME=$name\n";
print "export MCPHASE_GETVARIABLE_VALUE=$value\n";
                  }

if(defined $compare)
{@d=split("e|E|d|D",$compare);
 $d[0]=~s/\d(?=[\d\.]*?\d)/X/g;
 $d[0]=~s/\d/1/g;
 $d[0]=~s/X/0/g;
 $err=join("e",@d);
print STDERR "accuracy=".$err."\n";
 if(abs($compare-$value)>abs($err))
{die "Error getvalue comparing  $compare to extracted value $name=$value from file $filename \n";}
}

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
                if($line=~/^.*$var\s*=/) {($value)=($line=~m|$var\s*=\s*([^\s^>^<^=]+)|);}                                        }
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
