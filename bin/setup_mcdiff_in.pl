#!/usr/bin/perl
use Getopt::Long;

GetOptions("help"=>\$helpflag,
           "prefix|p=s"=>\$prefix);
usage() if $helpflag||$#ARGV<1;
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;$ARGV[0]=eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;$ARGV[1]=eval $ARGV[1];
if ($#ARGV>2) { 
$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;$ARGV[2]=eval $ARGV[2];
$ARGV[3]=~s/exp/essp/g;$ARGV[3]=~s/x/*/g;$ARGV[3]=~s/essp/exp/g;$ARGV[3]=eval $ARGV[3];
             }
print STDOUT << "EOF";
*******************************************************
setting up mcdisp.mf to be used by mcdisp
EOF
if ($prefix){$pp="-prefix ".$prefix;}
if ($#ARGV>2) { 
print STDOUT "T=$ARGV[0] K Ha=$ARGV[1] T Hb=$ARGV[2] T Hc=$ARGV[3] T\n";
system ("spins $pp $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]");
             }
else
            {
print STDOUT "x=$ARGV[0]  y=$ARGV[1] \n";
system ("spins $pp $ARGV[0] $ARGV[1]");
             }
print STDOUT << "EOF";
*******************************************************
reading results/mcphas.mf
.... trying to calculate results/spins.*
reading results/mcphas.mf by program spins ...
writing results/spins.* ...
EOF


print "\n\nSetting up mcdiff.in ...\n";

open (Fout,">mcdiff.in");
{open (Fin,"results/spins.out");
while(<Fin>) {print Fout $_;}
close Fout,Fin;
}

print STDOUT << "EOF";

    mcdiff.in generated: you can start now mcdiff
    However, please remember to set wavelength, lorentzfactor etc.
    in input file mcdiff.in

    You can view the magnetic structure in postscriptfiles
    results/spins*.eps, by fp_studio results/spins.fst and
    by javaview results/spins.jvx

EOF
exit;

sub usage() {

  print STDERR << "EOF";

    setup_mcdiff_in: program to setup mcdiff_in with information on spinconfiguration
                    to be used by program mcdiff. Note, you must
                    have done a mcphas calculation to stabilise
                    a magnetic structure at the desired Temperature/Field.
                    setup_mcdiff_in reads the results of this calculation
                    from results/mcphas.sps and results/mcphas.mf generates an
                    input file mcdiff.in

    usage: setup_mcdiff_in [option] T Ha Hb Hc
              or
           setup_mcdiff_in [option] x y

     -h          : this (help) message
      T          : Temperature (K)
      Ha,Hb,Hc   : Magnetic Field (T)
      x,y        : x,y values of point in the phase diagram

    required input files:

    results/mcphas.mf
                 :  result of a mcphas calculation

    output files:

    mcdiff.in    :  required input file for mcdiff

    - after running this program you can start mcdiff to do the calculation
      magnetic diffraction pattern

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
