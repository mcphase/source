#!/usr/bin/perl
use Getopt::Long;

GetOptions("help"=>\$helpflag);
usage() if $helpflag;
print STDERR "#* $0 *\n";
print STDOUT << "EOF";
*******************************************************
setting up exchange paramter fit 
*******************************************************
reading mcphas.j
EOF
open(Fin,"mcphas.j");

print STDOUT << "EOF";
writing mcphas.j.forfit
EOF


$lowboundary=-1e0;
$highboundary=1e0;
$stepwidth=1e-6;

open(Fout,">mcphas.j.forfit");


  while($line=<Fin>)
     { if ($line=~/^\s*#/) {print Fout $line;}
      else{$line=~s/D/E/g;@numbers=split(" ",$line);
           for($i=0;$i<=2;++$i) {print Fout $numbers[$i]." ";}
           for($i=3;$i<=$#numbers;++$i)
           {if(abs($numbers[$i])==0){print Fout "0.0 ";}
            else
             {$ii=-1;
              foreach(@values)
               {++$ii;
               if(abs($_-$numbers[$i])==0){print Fout "function\[$parname[$ii]\] ";$ii=-2;last;}
               }
              if($ii==$#values)   # introduce new parameter 
                 {$values[$ii+1]=$numbers[$i];
                  $parname[$ii+1]="par".$numbers[0].$numbers[1].$numbers[2]."_".($i-2);
                  $parname[$ii+1]=~s/\./_/g;$parname[$ii+1]=~s/-/m/g;$parname[$ii+1]=~s/\+/p/g;
                  print Fout $parname[$ii+1]."\[".$numbers[$i].",".$lowboundary.",".$highboundary.",0,".$stepwidth."\] ";
                 }
              }
           }
           print Fout "\n";
          }
     }
close Fout;
close Fin;

print STDOUT << "EOF";
      now you must setup a batch program calcsta 
      to calculate the standard deviation 
      
      then you can start
      a fit by simannfit or searchspace
EOF
sub usage() {

  print STDERR << "EOF";

    setup_mcphasjforfit: program to setup a fit of exchange parameters by
                         creating mcphas.j.forfit from mcphas.j
                    
    usage: setup_mcphasjforfit [-h]

     -h          : this (help) message

    required input files:

    mcphas.j (+ single ion parameter files)
                 :  structural information including all magnetic atoms

    output files:

    mcphas.j.forfit  : all interaction parameters values are substituted
                       with parJxxx[0.0,-1e0,1e0,0,1e-6]

    - all values which are zero are left untouched
    - all values  which are exactly equal are identified to be the same paremeter
 
    - after running this program you must setup a file calcsta 
      to calculate the standard deviation and then you can start
      a fit by simannfit or searchspace
EOF

  exit;

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
