#!/usr/bin/perl
use Getopt::Long;

GetOptions("help"=>\$helpflag,
           "prefix|p=s"=>\$prefix);
usage() if $helpflag||$#ARGV<1;
print STDERR "#* $0 *\n";
$ARGV[0]=~s/x/*/g;$ARGV[0]=eval $ARGV[0];
$ARGV[1]=~s/x/*/g;$ARGV[1]=eval $ARGV[1];
if ($#ARGV>2) { 
$ARGV[2]=~s/x/*/g;$ARGV[2]=eval $ARGV[2];
$ARGV[3]=~s/x/*/g;$ARGV[3]=eval $ARGV[3];
            }
print STDOUT << "EOF";
******************************************************************
* setup_mcdisp_mf 221011 setting up mcdisp.mf to be used by mcdisp
EOF
if ($#ARGV>2) { 
print STDOUT "T=$ARGV[0] K Ha=$ARGV[1] T Hb=$ARGV[2] T Hc=$ARGV[3] T\n";
$err=system ("spins -f results/".$prefix."mcphas.mf $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]  > mcdisp.mf");
             }
else
            {
print STDOUT "x=$ARGV[0]  y=$ARGV[1] \n";
$err=system ("spins -f results/".$prefix."mcphas.mf $ARGV[0] $ARGV[1]  > mcdisp.mf");
             }
if($err){unlink "mcdisp.mf";exit(1);}



print STDOUT << "EOF";
*******************************************************
reading results/mcphas.mf
....writing mcdisp.mf
reading results/mcphas.sps
....writing results/spins.*
EOF

system ("spins $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]");


print STDOUT << "EOF";

    You can start now mcdisp - required input files
    are mcdisp.mf (just created), mcdisp.par, mcphas.j
    and corresponding single ion property files *sipf.

    You can view the magnetic structure in postscriptfiles
    results/spins*.eps, by fp_studio results/spins.fst and
    by javaview results/spins.jvx
*********** end setup_mcdisp_mf ***************************
EOF
exit;

sub usage() {

  print STDERR << "EOF";

    setup_mcdisp_mf: program to setup mcdisp.mf with information on meanfields
                    to be used by program mcdisp. Note, you must
                    have done a mcphas calculation to stabilise
                    a magnetic structure at the desired Temperature/Field.
                    setup_mcdisp_mf reads the results of this calculation
                    from results/mcphas.mf and puts the meanfields into
                    mcdisp.mf

    usage: setup_mcdisp_mf [options] T Ha Hb Hc
             or
           setup_mcdisp_mf [options] x y

     -h          : this (help) message
      T          : Temperature (K)
      Ha,Hb,Hc   : Magnetic Field (T)
      x,y        : x,y values of point in the phase diagram

    required input files:

    results/mcphas.mf
                 :  result of a mcphas calculation

    output files:

    mcdisp.mf    :  required input file for mcdisp

    - after running this program you can start mcdisp to do the calculation
      of dispersion of excitations or diffuse scattering

    options:

    -prefix 001  :  instead of results/mcphas.mf read results/001mcphas.mf
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
