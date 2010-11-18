#!/usr/bin/perl
use Getopt::Long;

GetOptions("help"=>\$helpflag);
usage() if $helpflag||$#ARGV<2;
print STDOUT << "EOF";
*******************************************************
setting up exchange paramter fit for propagation vector
tau=($ARGV[0] $ARGV[1] $ARGV[2])
*******************************************************
reading mcphas.j
EOF

$nofatoms=extract("nofatoms","mcphas.j");
$nofcomponents=extract("nofcomponents","mcphas.j");

print STDOUT << "EOF";
writing mcdisp.par
EOF

open(Fout,">mcdisp.par");
print Fout << "EOF";
# Parameter file  mcdisp.par - created by jqfit_setup.pl
#<!--mcdisp.mcdisp.par>
#*********************************************************************
# mcdisp - program to calculate the dispersion of magnetic excitations
# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
#*********************************************************************
#
# mcdisp calculates the neutron scattering cross section dsigma/dOmegadE' [barn/sr/meV/f.u.]
#           f.u.=crystallogrpaphic unit cell (r1xr2xr3) for inelastic and diffuse scattering
#
# depending on what is kept constant it follows either kf or ki (1/A)
#!kf=10
#
# for full calculation of the dynamical susceptibility (option "-r", inversion of the MF-RPA equation
# for each point in Q-omega space) the minimum and maximum energy has to be given (energy stepwidth is
# equal to the parameter epsilon given in the command line after "-r")
#
#!emin=0.1
#!emax=10
#
# optional parameter is extended_eigenvector_dimension
# which is used to define, how many components of the
# eigenvector should be in the ouput to file mcdisp.qee
# (important for charge density movies)
#!extended_eigenvector_dimension=3
#
# It follows either
#
# (i) a Q vector mesh to be mapped in the calculation
#hmin=0 hmax=0.5 deltah=0.1
#kmin=0 kmax=0.5 deltak=0.1
#lmin=0 lmax=0.5 deltal=0.1
#
# or (if no mesh is given the program takes)
# (ii) a list of Q vectors with (optional) energies of observed excitations to be fitted
# h k l [E1(meV) E2(meV) ...]
# THE FIRST Q VECTOR IS THE PROPAGATION VECTOR - J(Q) SHOULD BE MAXIMAL FOR THIS VECTOR
$ARGV[0] $ARGV[1] $ARGV[2]
# now follows a list of other propagation vectors to be tested and compared
EOF
for($h=0;$h<=0.5;$h+=0.1){
for($k=0;$k<=0.5;$k+=0.1){
for($l=0;$l<=0.5;$l+=0.1){
if(abs($h-$ARGV[0]) +abs($k-$ARGV[1]) +abs($l-$ARGV[2])>0.02 )
{print Fout sprintf("%2.1f  %2.1f %2.1f\n",$h,$k,$l);}
}}}
close Fout;

print STDOUT << "EOF";
writing mcdisp.mf
EOF

open (Fout, ">mcdisp.mf");
print Fout << "EOF";
# Parameter file  mcdisp.mf - created by jqfit_setup.pl
#<!--mcdisp.mcdisp.mf>
#*********************************************************************
# mcdisp - program to calculate the dispersion of magnetic excitations
# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
#*********************************************************************
#'T'             temperature T
#'Ha' 'Hb' 'Hc'  magnetic field
#'n'             number of atoms in magnetic unit cell
#'nofatoms'      number of atoms in primitive crystal unit cell
#'nofcomponents' dimension of moment vector of a magnetic atoms
T=1.4 Ha=5 Hb=0 Hc=0 n=1 nofatoms=$nofatoms nofcomponents=$nofcomponents
 1.0
EOF
for($i=2;$i<=$nofatoms*$nofcomponents;++$i) {print Fout "0.0000\n";}
close Fout;
print STDOUT << "EOF";
writing calcsta
EOF

open (Fout,">calcsta");
print Fout << "EOF";
 perl calcsta.pl
EOF
close Fout;

open (Fout,">calcsta.bat");
print Fout << "EOF";
 perl calcsta.pl
EOF
close Fout;

print STDOUT << "EOF";
writing calcsta.pl.forfit
EOF
open(Fout,">calcsta.pl.forfit");
print Fout << "EOF";
#!\\usr\\bin\\perl
use File::Copy;

# **************************************************
# fit parameter file created by jqfit_setup.pl
# **************************************************

#-------------------------------
# uncomment one of the following suggestions for calculating
# interaction constants
#-------------------------------

# OPTION I
# Bethe Slater curve with a fixed value  of -0.015 at 3.012 A
\$makenn="makenn 14 -kaneyoshi function[-0.015*exp(paralpha*3.012*3.012/(parD*parD))/(-3.012*3.012/(parD*parD)+3.012*3.012*3.012*3.012/(parD*parD)/(parD*parD))] parD [2.101336e+00,1.500000e+00,7.500000e+00,1.091110e-02,3.041289e-01]  paralpha [5.000000e-01,5.000000e-01,4.500000e+00,2.269634e-03,4.963513e-01] -d";

# OPTION II
# RKKY function with a fixed value  of 0.015 at 3.012 A
#\$makenn="makenn 20 -rkky function[-0.015*8*parkf*parkf*parkf*3.012*3.012*3.012/cos(2*parkf*3.012)] parkf [4.965124e-01,3.000000e-01,2.500000e+00,5.412158e-02,8.487453e-01] -d";

# OPTION III
# Bethe Slater curve with a fixed value
# of 0.015 at 3.012 A, with anisotropy (requires less than cubic symmetry)
#\$makenn="makenn 14 -kaneyoshi3d function[-0.015*exp(paralpha*3.012*3.012/(parDa*parDa+parDc*parDc))/(-3.012*3.012/(parDa*parDa+parDc*parDc)+3.012*3.012*3.012*3.012/(parDa*parDa+parDc*parDc)/(parDa*parDa+parDc*parDc))] parDa [2.000000e+00,2.000000e+00,7.500000e+00,3.187591e-02,2.769258e+00] function[parDa] parDc [7.500000e+00,2.000000e+00,7.500000e+00,8.656280e-02,5.956705e+00] paralpha [3.193363e+00,1.000000e+00,4.000000e+00,9.762790e-03,3.384279e+00] -d";

# OPTION IV
# RKKY curve with a fixed value
# of 0.015 at 3.012 A, with anisotropy (requires less than cubic symmetry)
# \$makenn="makenn 10 -rkky3d function[-0.015*8*parka*parka*parka*3.012*3.012*3.012/cos(2*parka*3.012)] parka [3.331250e+000,1.000000e-001,4.500000e+000,6.670118e-004,4.534191e-002] function[parka] parkc [2.300000e+000,1.000000e-001,4.500000e+000,5.194754e-004,8.930458e-002] -d";


#-------------------------------
system("\$makenn");     # start makenn to calculate the interactions
copy("results/makenn.j","mcphas.j"); # copy interactions to file mcphas.j
system("delcol 7 mcphas.j"); # remove redundant columns
system("delcol 7 mcphas.j");

system("mcdispit -jq "); # start mcdisp to calculate J(Q)

(\$sta)=extract("sta","results/mcdisp.jq");
# --- let sta be greater than one if tau is not maximum ...
# so a value smaller than one indicates that J(Q) with maximum at tau has
# been found ...
print "sta=".(\$sta+1)."\\n";
exit;

# **********************************************************************************************
# extracts variable from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# (\$standarddeviation)=extract("sta","data.dat");
#
# ... it stores 0.24 in the variable \$standarddeviation
#
sub extract {
             my (\$variable,\$filename)=\@_;
             \$var="\\Q\$variable\\E";
             if(open (Fin,\$filename))
             {while(\$line=<Fin>){
                if(\$line=~/^.*\$var\\s*=/) {(\$value)=(\$line=~m|\$var\\s*=\\s*([\\d.eEdD\\Q-\\E\\Q+\\E]+)|);}                                        }
              close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \\"\$filename\\"\\n";
             }
             return \$value;
            }
# **********************************************************************************************
EOF
close Fout;

print STDOUT << "EOF";
writing fit.bat
EOF

open(Fout,">fit.bat");
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

call display 7 4 results/makenn.j &
call simannfit 1e-3 calcsta.pl
EOF

print STDOUT << "EOF";

    You can start now immediately a fit of exchange
    parameters. Edit calcsta.pl.forfit and fit.bat
    to fine tune the fit according to your needs.
    Then type:    fit.bat

EOF
exit;

sub usage() {

  print STDERR << "EOF";

    $0: program to setup a fit of exchange parameters in order
                    to reproduce an experimental propagation vector

    usage: $0 [-h] [--help] h k l

     -h          : this (help) message
      hkl        : Miller indices of propagation vector

    required input files:

    mcphas.j (+ single ion paramter files)
                 :  structural information including all magnetic atoms

    output files:

    mcdisp.par   :  contains propagation vector and list of other hkl to
                    be probed
    mcdisp.mf    :  required input file for mcdisp
    calcsta      :  required input file for simannfit and searchspace
    calcsta.pl.forfit: file with fitparameters for Bethe slater, RKKY fits
    fit.bat      :  batch to start the fit

    - after running this program you can start immediately a fit of exchange
    parameters. Edit calcsta.pl.forfit and fit.bat to fine tune the fit 
    according to your needs.
    During fitting a value of sta < 1 indicates, that the maximum of J(Q) is
    at the propagation vector tau. How much it is below one depends on the
    magnitude of J(Q) for the competing wavevectors in the list in mcdisp.par.
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
