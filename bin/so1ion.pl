#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use File::Copy;
# copy parameter files to specific locations
# remember a few useful perl commands
# copy("alter_name","neuer_name");
# rename("alter_name","neuer_name");
# chdir('../.');
print STDOUT << "EOF";

Progam so1ion - calculation of single ion problems in LS coupling
schemes (Hee >> Hso >>  Hcef ~ Hze)

EOF

system "so1ionit @ARGV";

# check if command line is reading file for neutron intensity calculation (option -r)
unless ($ARGV[0]=~m/-r/&&$ARGV[2]=~m/-[B,L,b,l]/) {exit(0);}
# if yes, check if file exists and if it starts with #!MODULE=so1ion
unless (open(Fin,$ARGV[1])){print "error: file $ARGV[1] does not exist\n";exit(0);}
$line=<Fin>;close Fin;
unless($line=~m/\Q#!MODULE=so1ion\E/){print "error: file $ARGV[1] does not start with #!MODULE=so1ion\n";exit(0);}
# if yes then do also a mcdisp calculation to create so1ion.trs ...

($T)=extract("T","$ARGV[1]");
($gJ)=extract("GJ","$ARGV[1]");
($Bx)=extract("Bx","$ARGV[1]");
($By)=extract("By","$ARGV[1]");
($Bz)=extract("Bz","$ARGV[1]");



# set up mcdisp.mf
if (open(Fin,"mcdisp.mf")) {@mcdispmf=<Fin>;close Fin;}
open (Fout, ">mcdisp.mf");
print Fout << "EOF";
#****************************************************
# spins - display spinconfiguration at given htpoint
# Author: Martin Rotter mcphas version 4.0
#****************************************************
#! T=$T Ha=$Bx Hb=$By Hc=$Bz n=1 spins nofatoms=1 in primitive basis nofcomponents=3 - momentum configuration <J(i)>
EOF
$MB=5.788378E-02;
print Fout ($gJ*$Bx*$MB)."\n";
print Fout ($gJ*$By*$MB)."\n";
print Fout ($gJ*$Bz*$MB)."\n";
close Fout;

# set up mcdisp.par
if (open(Fin,"mcdisp.par")) {@mcdisppar=<Fin>;close Fin;}
open (Fout, ">mcdisp.par");
print Fout << "EOF";
# autocreated Parameter file  mcdisp.par
#<!--mcdisp.mcdisp.par>
#!kf=100000
#!extended_eigenvector_dimension=3
0.1 0 0
0 0.1 0
0 0 0.1
EOF
close Fout;

# set up mcphas.j
if (open(Fin,"mcphas.j")) {@mcphasj=<Fin>;close Fin;}
open (Fout, ">mcphas.j");
print Fout << "EOF";
# autocreated file for so1ion
#<!--mcphase.mcphas.j-->
#***************************************************************
# Lattice Constants (A)
#! a=1 b=1 c=1 alpha=  90 beta=  90 gamma=  90
#! r1a= 1 r2a=   0 r3a=   0
#! r1b= 0 r2b=   1 r3b=   0   primitive lattice vectors [a][b][c]
#! r1c= 0 r2c=   0 r3c=   1
#! nofatoms=1  nofcomponents=3  number of atoms in primitive unit cell/number of components of each spin
# ****************************************************************************
# ****************************************************************************
#! da=   0 [a] db=   0 [b] dc=   0 [c]   nofneighbours=0 diagonalexchange=1 gJ=$gJ cffilename=$ARGV[1]
# da[a]      db[b]      dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]
EOF
close Fout;
if (open(Fin,"results/mcdisp.trs")) {@mcdisptrs=<Fin>;close Fin;}

#
# start mcdispit -c
if(system "mcdispit -c > range.out ") {print "problem starting mcdispit";}
else
{print "Transitions and single ion neutron intensities stored in results/so1ion.trs";}

mydel("range.out");
# remove files
if (@mcdispmf){open (Fout, ">mcdisp.mf"); print Fout @mcdispmf;close Fout;}
else {mydel("mcdisp.mf");}
if (@mcdisppar){open (Fout, ">mcdisp.par"); print Fout @mcdisppar;close Fout;}
else {mydel("mcdisp.par");}
if (@mcphasj){open (Fout, ">mcphas.j"); print Fout @mcphasj;close Fout;}
else {mydel("mcphas.j");}


# rename mcdisp.trs into so1ion.trs
 chdir('./results');
 rename("mcdisp.trs","so1ion.trs");
 chdir('./..');

if (@mcdisptrs){open (Fout, ">results/mcdisp.trs"); print Fout @mcdisptrs;close Fout;}



sub mydel  { my ($file1)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;
                               return system("del ".$file1);
                              }
                 else
                              {return system("rm ".$file1);
                              }

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
