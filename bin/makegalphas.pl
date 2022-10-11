#!/usr/bin/perl
# 
unless ($#ARGV>5) {die "
# perl program to determine the crystal field phonon coupling from the
# pointcharge model and store the coupling constants in results/g_alpha.s
#
# usage: makegalphas rmax anr Qh Qk Ql nofatoms mass1 mass2 mass3 ... mass_nofatoms
#                                        
#  (Qh Qk Ql) q-Vector in  primitive basis
#  rmax .... maximum distance [Angstroem] to be considered in pointcharge model
#            rmax=0 triggers the program to take pointcharge calculation from previous run
#                   (mind: all other parameters except Qh Qk Ql and phonon.ev have to be the
#                          same in this case, results are appended to results/g_alpha.s)
#  anr  .... index of magnetic atom 
#  nofatoms ... total number of atoms in primitive unit cell followed by
#  mass*  ...  mass number of these atoms
#
# required input files:
#  phonon.ev  ... required input file with Parlinsky d46 file type Phonon eigenvectors
#  mcphas.j & *.sipf  ... crystal structure data, in .sipf files CHARGE= 
#                  indicates the point charge to be used in the calculation
#
# required programs: ga, pointc, newcols, swapf
#
# output:
# 
#  results/g_alpha.s   file with Crystal Field Phonon coupling constants
#  results/makegalpha.sipf  sipf file of magnetic ion with crystal field parameters
#              calculated according to the pointcharge model
#
# formula:
#
# g_alpha(qs)= - sum_j  dBalpha(i)/dSj ej_q,s  sqrt[hbar/(2 Mj omega_q,s)] exp(i q(Rj-Ri))
#
# Balpha(i) .... Stevens CF Parameters alpha=lm at crystal site i according to 
#                pointcharge model calculation
# Sj ........... phonon displacement at site j
# hbar ......... Planck Constant h/2pi
# Mj ........... mass of nucleus at site j
# omega_q,s .... phonon frequency at wave vector q, branch s
# Rj ........... position vector of site j\n";}else{print STDERR "#* $0 *\n";}

# read paramters from command line

$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($rmax) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($anr) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($Qh) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($Qk) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($Ql) = eval $ARGV[0];
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($nofa) = eval $ARGV[0];

for($i=1;$i<=$nofa;++$i){
shift @ARGV; 
$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
$MASS[$i] = eval $ARGV[0];
}

$PI=3.141592654;
$file=">./results/g_alpha.s";
unless($rmax==0){

# initialize  output file 
unless (open (Fout, $file)){die "\n error:unable to open $file\n";}
print Fout "# alpha s |g_alpha(s)| Re(g_alpha(s)) Im(g_alpha(s)) [meV] (h k l)  omegaqs[meV] \n"; close Fout;

# the ga programm creates neighbours from mcphas.j structural information and 
# from this a list of relative positions to Ce3+ is created and stored in 
# ga0.pc  ...  without oscillation, i.e. normal pointcharge model 

       system("ga $rmax $anr $Qh $Qk $Ql 1"); # we set n=1 i.e. no oscillation simulated, Qh Qk Ql irrelevant

#    with pointc the derivatives dBalpha/du are created and stored in file results\pointc.dBlm
       system("pointc -d Ce3+ results/ga0.pc > results/makegalpha.sipf");

# create from dBalpha/du  cf-phonon interaction according to
#  g_alpha(qs)= - sum_j  dBalpha(i)/dSj ej_q,s  sqrt[hbar/(2 Mj omega_q,s)] exp(i q(Rj-Ri))
#  DBalpa/du are in file pointc.dBlm 
#  ejqs/Mj eigenvector of phonons are in col[1-15]+icol[16-30]
#  omegaqs=col[0] in THz

# enter self-ion into last line of ga0.pc
system ("echo 0 0 0 0 0 0 0 0 $anr 0 0 0 0 0 0 0 0 0 >> results/ga0.pc");
# swap from ga0.pc information about  Nr1 Nr2 Nr3 atomnr to pointc.dBlm to have everything needed in one file ...
system(" newcols 1 4 results/pointc.dBlm");
system(" swapf 1 results/pointc.dBlm 9 results/ga0.pc");
system(" swapf 2 results/pointc.dBlm 13 results/ga0.pc");
system(" swapf 3 results/pointc.dBlm 14 results/ga0.pc");
system(" swapf 4 results/pointc.dBlm 15 results/ga0.pc");
system(" swapf 5 results/pointc.dBlm 16 results/ga0.pc");
system(" swapf 6 results/pointc.dBlm 17 results/ga0.pc");
system(" swapf 7 results/pointc.dBlm 18 results/ga0.pc");
# now atomnumnber is in column 1 and integer unit cell numbers with respect to primitive lattice
#  Nr1 Nr2 Nr3 are in columns 2-4 of file pointc.dBlm
# atomic positions dr1 dr2 dr3 with respect to primitive lattice are in columns 5-7 of file pointc.dBlm
# using $Qh $Qk $Ql we can evaluate Q.R

}

unless (open (Fout, ">".$file)){die "\n error:unable to open $file\n";}

$fileev="phonon.ev"; unless (open (Fin, $fileev)){die "\n error:unable to open $fileev\n";}
$j=0;
while($line=<Fin>)
     { if ($line=~/^\s*#/) {}
       else{$line=~s/D/E/g;@col=split(" ",$line);
           	  ++$j;  
       if($#col>1){

#  g_alpha(qs)= - sum_j  dBalpha(i)/dSj ej_q,s  sqrt[hbar/(2 Mj omega_q,s)] exp(i q(Rj-Ri))
# we  need to take into account the factor 1/sqrt(Mj) in the above formula and multiply
# phonon eigenvector components accordingly

#now we have to divide by sqrt(mass) to get an oscillation amplitude

for($i=1;$i<=$nofa;++$i){
       $col[3*($i-1)+1]/=sqrt($MASS[$i]);$col[3*$nofa+3*($i-1)+1]/=sqrt($MASS[$i]);
       $col[3*($i-1)+2]/=sqrt($MASS[$i]);$col[3*$nofa+3*($i-1)+2]/=sqrt($MASS[$i]);
       $col[3*($i-1)+3]/=sqrt($MASS[$i]);$col[3*$nofa+3*($i-1)+3]/=sqrt($MASS[$i]);
                       }
      

# now open results/pointc.dBlm and sum up galphas
$fil="results/pointc.dBlm";
unless (open (FD, $fil)){die "\n error:unable to open $fil\n";}
my @gas=(),@gasr=(),@gasi=();
for($ilm=1;$ilm<=27;++$ilm){$gasr[$ilm]=0;$gasi[$ilm]=0;} # clear array  for summation
$cn=0;
while($line=<FD>)
     { if ($line=~/^\s*#/) {}
       else{$line=~s/D/E/g;@cl=split(" ",$line);           	  
       if($#cl>1&&$cl[0]>0){ ++$cn; # mind in pointc.dBlm last line is for the self momvement of the magnetic ion
#  g_alpha(qs)= - sum_j  dBalpha(i)/dSj ej_q,s  sqrt[hbar/(2 Mj omega_q,s)] exp(i q(Rj-Ri))
#     Sj=uj*a0  
#  g_alpha(qs)= - sum_j  dBalpha(i)/duj ej_q,s  sqrt[hbar^2/(2 Mj a0 a0 hbar omega_q,s)] exp(i q(Rj-Ri))
# 
# divide by sqrt of Energy (omegaqs)
$factor=-1/sqrt($col[0]); #  MIND:here is still missing  sqrt(hbar hbar/(a0 a0 2 m0 meV)) !!!!!!!
#    a0=0.5292 A= 0.5292e-10 m
#    m0=1.672e-27 kg
#     hbar=1.055e-34 Js
#    1 meV=1.602e-22 J
#  hbar/a0 a0 2 m0 meV= (1.055e-34 Js)^2 / ( 2 0.5292 0.5292e-20 m^2 1.672e-27 kg 1.602e-22 J ) =
#  = 1.055 1.055e+1 / (2 0.5292 0.5292 1.672 1.602 ) = 7.419
$factor*=sqrt(7.419);  # so  units should be correct now 

# evaluate exp(iQRj-Ri)=$expr + i $expi  Rj and Ri lattice vectors
$QR=2*$PI*($Qh*$cl[1]+$Qk*$cl[2]+$Ql*$cl[3]);
$expr=cos($QR);      
$expi=sin($QR);      

# evaluate exp(iQ(rmu)) .... rmu  atomic positions in primitive unit cell
$rmu1=$cl[4]-$cl[1];
$rmu2=$cl[5]-$cl[2];
$rmu3=$cl[6]-$cl[3];

$QRp=2*$PI*($Qh*$rmu1+$Qk*$rmu2+$Ql*$rmu3);
$exprp=cos($QRp);   
$expip=sin($QRp);      


$atomnr=$cl[0];
 ($nofatoms)=extract("nofatoms","mcphas.j");
 unless($nofatoms==$nofa){die "ERROR makgalphas.pl: nofatoms in commandline not equal to nofatoms in mcphas.j\n";}
 #if ($j==3&&$cn<2){print "\n phonon nr $j \n neighbournr atomnr xyz  polarization/srqt(Mmu) conventional : parlinsky \n";}
for($xyz=1;$xyz<=3;++$xyz){ # loop euclidean coordinates of eigenvectors
$ejrparlinsky=$col[($atomnr-1)*3+$xyz]; # eigenvector real and imag parts
$ejiparlinsky=$col[($atomnr-1)*3+$xyz+3*$nofatoms]; # $nofatoms is the number of atoms in the unit cell !!! 

# transform parlinsky vectors to standard vectors by phase transformation ... see manual parlinsky:
# "Polarization vectors e(k, j; μ), calculated by Phonon, differ from the conventional
#   once e(kτ , j; μ), which usually are defined in the τ Brillouin zone they belong to.
#   The Phonon’s polarization vector e(k,j;μ) is defined for the wave vector k pointing out 
#   from the absolute origin of the reciprocal space. The conventional polarization vector
#    eτ (kτ , j) is defined for the wave vector kτ pointing out from the center of a given 
#   Brillouin zone labelled by reciprocal lattice vector τ. The relation between the two type
#   of polarization vectors is easy to derive. Inserting the relation k = τ + kτ to the 
#   equation Eq.(34), one finds e(k, j, μ) = eτ (kτ , j; μ)exp[−2πτ · rμ] "
#  
$ejr=$ejrparlinsky*$exprp-$ejiparlinsky*$expip;
$eji=$ejrparlinsky*$expip+$ejiparlinsky*$exprp;
#if ($j==3&&$cn<20){print "$cn $atomnr $xyz  $ejr + i $eji : $ejrparlinsky + i $ejiparlinsky \n ";}

for($ilm=1;$ilm<=27;++$ilm){
# dB22S/dux ... $cl[10] ... gas[1] usw.
# $gasr[1]+=$cl[9+$xyz]*($ejr*$expr-$eji*$expi)*$factor;
# $gasi[1]+=$cl[9+$xyz]*($ejr*$expi+$eji*$expr)*$factor;
# dB21S/dux ... $cl[13] ... gas[2] usw.
# $gasr[2]+=$cl[12+$xyz]*($ejr*$expr-$eji*$expi)*$factor;
# $gasi[2]+=$cl[12+$xyz]*($ejr*$expi+$eji*$expr)*$factor;
# ...
$gasr[$ilm]+=$cl[9+($ilm-1)*3+$xyz]*($ejr*$expr-$eji*$expi)*$factor;
$gasi[$ilm]+=$cl[9+($ilm-1)*3+$xyz]*($ejr*$expi+$eji*$expr)*$factor;
 
    }
                          }
       }}}
close FD;


for($ilm=1;$ilm<=27;++$ilm){$gas[$ilm]=sqrt($gasr[$ilm]*$gasr[$ilm]+$gasi[$ilm]*$gasi[$ilm]);
          if($ilm<6){$alpha=$ilm+3;}                 #dB2m
          else{if($ilm<15){$alpha=$ilm+10;}                 #dB4m
               else{$alpha=$ilm+21;}                 #dB6m
               }
print Fout "$alpha $j $gas[$ilm] $gasr[$ilm] $gasi[$ilm] $Qh $Qk $Ql $col[0]\n";

$gasr[$ilm]=0;$gasi[$ilm]=0;
}

           }}     }
close Fin;

close Fout;

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
             my ($variable,$filename)=@_; my $line="";
             $var="\Q$variable\E";$value="";
             if(open (Fin1,$filename))
             {while($line=<Fin1>){
                if($line=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($line=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);}}
              close Fin1;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";
             }
             return $value;
            }
# **********************************************************************************************
