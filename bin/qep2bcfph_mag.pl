#!/usr/bin/perl
# 
unless ($#ARGV>3) {die "
# perl program to convert output file mcdisp.qep into
# results/g_alpha.s to be used for bcfph 
# model calculations of crystal field spectra with CF-Phonon Interaction 
#
# usage: qep2bcfph_mag $rmax $anr  qepfilename nofatoms MASS1 MASS2   ... MASS_nofatoms
#
#  qepfilename .. filename of output file of mcdisp, e.g. results/mcdisp.qep
#                 all q-vectors stored in this file will be used and averaged
#                 in the calculation
#  nofatoms ..... number of atoms in primitive basis
#  MASS* ........ masses of atoms\n";}
#
# read paramters from command line

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($rmax) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($anr) = eval $ARGV[0];
shift @ARGV; 
my ($filename)=$ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($nofa) = eval $ARGV[0];

for($i=1;$i<=$nofa;++$i){
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$MASS[$i] = eval $ARGV[0];
}

unless (open (Fin, $filename)){die "\n error:unable to open $filename\n";}
$j=1;$N=0;$s=0;open(Foutev,">phonon.ev");close Foutev;

while($line=<Fin>)
     { if ($line=~/^\s*#/) {}
 else{$line=~s/D/E/g;@col=split(" ",$line); $ev=$col[8]." ";$ol[0]=$col[8];
      $line=<Fin>; $line=<Fin>;  
      while($line=~/^\s*[^#^\n]/){@cl=split(" ",$line);$ol[$j]=$cl[0];++$j;$line=<Fin>;}
      $line=<Fin>;$line=<Fin>;
      while($line=~/^\s*[^#^\n]/){@cl=split(" ",$line);$ol[$j]=$cl[0];++$j;$line=<Fin>;}
      $j=1;unless($#ol==$nofa*6){die "ERROR qep2galphas - $nofa nofatoms x 6 not equal to eigenvector dimension $#ol\n";}
      if($col[8]>=0)
{++$s;


# mcdisp outputs as eigenvector an oscillation amplitude ... to get a conventional
 # phonon eigenvector we have to multiply these by sqrt(mass)
 for($i=1;$i<=$nofa;++$i){
       $ol[3*($i-1)+1]*=sqrt($MASS[$i]);$ol[3*$nofa+3*($i-1)+1]*=sqrt($MASS[$i]);
       $ol[3*($i-1)+2]*=sqrt($MASS[$i]);$ol[3*$nofa+3*($i-1)+2]*=sqrt($MASS[$i]);
       $ol[3*($i-1)+3]*=sqrt($MASS[$i]);$ol[3*$nofa+3*($i-1)+3]*=sqrt($MASS[$i]);
                         }
# calculate norm
$norm=-$ol[0]*$ol[0];foreach(@ol){ $norm+=$_*$_;} 
# normalize eigenvector
        for($i=1;$i<=$#ol;++$i){$ol[$i]/=sqrt($norm);$ev=$ev." ".$ol[$i];}
# ... now we have from mcdisp output mcdisp.qep calculated a phonon polarisation eigenvector
# with similar properties as the polarisation eigenvectors output in parlinsly d46 files
#$norm=-$ol[0]*$ol[0];foreach(@ol){ $norm+=$_*$_;} 
open(Foutev,">>phonon.ev");

print Foutev $ev."\n"; 
# print "$col[8] meV Inuc=$col[11] |ev|^2=$norm  2dim=$#ol\n";
close Foutev;
if($s==$nofa*3){$s=0;++$N;
print "N=$N ( $col[12]  $col[13] $col[14] ) 2dim=$#ol\n";

if($N==1){system("makegalphas $rmax $anr $col[12]  $col[13] $col[14] $nofa @MASS ");
        }
else     {system("makegalphas 0     $anr $col[12]  $col[13] $col[14] $nofa @MASS ");} 
        open(Foutev,">phonon.ev");close Foutev;
        open(Foutev,">>results/g_alpha.s");print Foutev "# N=$N \n";close Foutev;
        }}

     }} 
      
close Fin;
