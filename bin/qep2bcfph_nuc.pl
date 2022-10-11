#!/usr/bin/perl
# 
unless ($#ARGV>4) {die "
# perl program to convert output file mcdisp.qep into
# phonon.ev and phonon.int to be used for bcfph
# model calculations of coherent nuclear spectra with CF-Phonon Interaction 
#
# usage: qep2bcfph_nuc h k l qepfilename nofatoms MASS1 MASS2   ... MASS_nofatoms
#
#  (h k l) ...... q-Vector in  nonprimitive basis, to be retrieved from qepfilename
#  qepfilename .. filename of output file of mcdisp, e.g. results/mcdisp.qep
#  nofatoms ..... number of atoms in primitive basis
#  MASS* ........ masses of atoms in sequence as in mcphas.j\n";}else{print STDERR "#* $0 *\n";}
#
# read paramters from command line

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($H) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($K) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($L) = eval $ARGV[0];
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
open(Foutint,">phonon.int");
open(Foutev,">phonon.ev");
$j=1;
while($line=<Fin>)
     { if ($line=~/^\s*#/) {}
 else{$line=~s/D/E/g;@col=split(" ",$line); $ev=$col[8]." ";$ol[0]=$col[8];
      $line=<Fin>; $line=<Fin>;  
      while($line=~/^\s*[^#^\n]/){@cl=split(" ",$line);$ol[$j]=$cl[0];++$j;$line=<Fin>;}
      $line=<Fin>;$line=<Fin>;
      while($line=~/^\s*[^#^\n]/){@cl=split(" ",$line);$ol[$j]=$cl[0];++$j;$line=<Fin>;}
      $j=1;unless($#ol==$nofa*6){die "ERROR qep2bcfph - $nofa nofatoms x 6 not equal to eigenvector dimension $#ol\n";}
      if($col[8]>=0&&$H==$col[4]&&$K==$col[5]&&$L==$col[6])
{

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

print Foutev $ev."\n"; print Foutint $col[8]." ".$col[11]."\n";
                  print "$col[8] meV Inuc=$col[11] |ev|^2=$norm  2dim=$#ol\n";

}
     }} 
      
close Fin;close Foutint;close Foutev;
