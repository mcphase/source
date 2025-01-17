#!/usr/bin/perl
use File::Copy;
use Getopt::Long;
use Cwd qw(getcwd);

BEGIN{@ARGV=map{glob($_)}@ARGV}
unless ($#ARGV >=0)
{
print STDERR << "EOF";

program to test a script
 
 usage: perl test.pl [options] test.bat

will execute the script test.bat line by line and stop with an error message
if an error occurs.

For testing  goto directory demo and look out for testing batch files
named test*.bat

- use such a file to test mcphase, e.g.:
mctest test_ndcu2_mcphasit.bat

- to test a mcphase installation with all tests type
mctest test*.bat

- to test a specific command occuring in several test-batch files, e.g. makenn, try
to run all these batch files. First, check by searching all test-batch files with
grep makenn test*.bat
 ... and run them by:
grep -l makenn test*.bat | xargs mctest

EOF
exit(1);
}else{print STDERR "#* $0 *\n";}

foreach(@ARGV)
{$file=$_;

unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
$linenr=0;   
while($line=<Fin>)

     {++$linenr;
      # do some substitutions to allow windows and linux batches 
if ($^O=~/MSWin/){$line=~s/^\s*rem/#/;
                  $line=~s/^\s*REM/#/;
                  $line=~s/^\s*rm\s/del /;
                  $line=~s/^\s*cp\s/copy /;


                 } else
                 {$line=~s/^\s*call\s//;
                  $line=~s/^\s*del\s/rm /;
                  $line=~s/^\s*copy\s/cp /;
                  $line=~s/^\s*rem/#/;
                  $line=~s/^\s*REM/#/;
                  
                 }


     unless ($line=~/^\s*[#\n]/){
               if($line=~/^\s*cd/){$line=~s/^\s*cd//;$line=~s/\n//;$line=~s/\s*//;
                                    unless(chdir($line)){die "\nError  executing command \ncd  $line in $file line number $linenr\n";}
                                  } # cd
                                  else
                                  {
                                   if(system($line)){die "\nError  executing command \n $line in $file line number $linenr \n";}
                                  }
                                } # no comment line
   } # next line
close Fin;
print "#************************\n";
print "# mctest <".$file. "> OK\n";
print "#************************\n";

} # next batch file
foreach(@ARGV)
{$file=$_;
print "# mctest <".$file. "> OK\n";
}
print "#**** END program mctest ******\n";


