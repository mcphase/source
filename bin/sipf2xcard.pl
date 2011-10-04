#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

use Getopt::Long;
use Math::Trig;
#use PDL;
#use PDL::Slatec;
use Switch;
use File::Copy;


@command=@ARGV;
@AARGV=@ARGV;
GetOptions("help"=>\$helpflag);
usage() if $helpflag||$#AARGV<0;
$file=$AARGV[0];shift @AARGV;

$PI=3.141592654;
print STDOUT << "EOF";
//*******************************************************
//                 program sipf2xcard
//         converts sipf file $file to xcard file
//*******************************************************
EOF
unless (open(Fin,$file)){die "cannot open file $file\n";}
$line=<Fin>;

($module)=($line=~m|\#\!(.*)|);
if ($module=~/MODULE\=/){($module)=($line=~m|\#\!MODULE\=(.*)|);}
if ($module=~/.*ic1ion/){$module="ic1ion";}
if ($module=~/.*icf1ion/){$module="icf1ion";}
if ($module=~/.*so1ion/){$module="so1ion";}
# print $module;exit;
# initialise CEF parameters
$L20=0.0;
$L21s=0.0;
$L22s=0.0;
$L21=0.0;
$L22=0.0;
$L40=0.0;
$L41s=0.0;
$L42s=0.0;
$L43s=0.0;
$L44s=0.0;
$L41=0.0;
$L42=0.0;
$L43=0.0;
$L44=0.0;
$L60=0.0;
$L61=0.0;
$L62=0.0;
$L63=0.0;
$L64=0.0;
$L65=0.0;
$L66=0.0;
$L61s=0.0;
$L62s=0.0;
$L63s=0.0;
$L64s=0.0;
$L65s=0.0;
$L66s=0.0;
switch ($module)
{  case "ic1ion"  {ic1ion()}
   case "icf1ion" {icf1ion()}
   case "so1ion"  {so1ion()}
  else {die "ERROR program sipf2xcard: module $module not implemented\n";}
}

# transform CEF parameters to eV
$L20*=0.001;
$L21s*=0.001;
$L22s*=0.001;
$L21*=0.001;
$L22*=0.001;
$L40*=0.001;
$L41s*=0.001;
$L42s*=0.001;
$L43s*=0.001;
$L44s*=0.001;
$L41*=0.001;
$L42*=0.001;
$L43*=0.001;
$L44*=0.001;
$L60*=0.001;
$L61*=0.001;
$L62*=0.001;
$L63*=0.001;
$L64*=0.001;
$L65*=0.001;
$L66*=0.001;
$L61s*=0.001;
$L62s*=0.001;
$L63s*=0.001;
$L64s*=0.001;
$L65s*=0.001;
$L66s*=0.001;

$F2*=0.001;
$F4*=0.001;
$F6*=0.001;
$xi*=0.001;

$mb=0.05788e-3; # bohr magneton in eV/T
$Bx*=$mb;$By*=$mb;$Bz*=$mb;$Bdirx=1;$Bdiry=0;$Bdirz=0;
$Babs=$Bx*$Bx+$By*$By+$Bz*$Bz; # set absolute value of applied magnetic field in eV
if ($Babs>0){$Babs=sqrt($Babs);
             $Bdirx=$Bx/$Babs;
             $Bdiry=$By/$Babs;
             $Bdirz=$Bz/$Babs;
            }
open (FOUT,">T.in");
print FOUT "$T\n";
close FOUT;
print "file T.in created\n";
open (FOUT,">$file.xcardx");
print FOUT << "EOF";
 XCRD:
 //XCARD for $file
 (


   // all parameters in eV
   Bapplied={$Babs,$Bdirx,$Bdiry,$Bdirz};

    // CRYSTAL FIELD
    //10Dq(#i1 $conf )=1.0;
    // Ak are the coefficients of <<spherical harmonics *sqrt(4pi/(2l+1))>>, i.e.
    // Llm are the coefficients of <<tesseral harmonics (small z)>>
    // Ll0=Al0,
    // Alm=(-1)^m*[Llm - i * Ll-m]    for m>0
    // Al-m=[Llm+ i Ll-m]
    // i.e. for L40=5.74755   A40=5.74755 eV
    //      for L44=3.43482   A44=3.43482,   A4-4=3.43482,
    //Ak(#i1 $conf)={4,0,5.74755,
    //                4,4,3.43482,
    //               4,-4,3.43482};
EOF


                                       print FOUT "    Hcf={2,0,$L20,0\n";
if (abs($L21*$L21+$L21s*$L21s)>1e-10) {print FOUT "              ,2,1,".(-$L21).",".( $L21s).",2,-1,$L21,$L21s\n";}
if (abs($L22*$L22+$L22s*$L22s)>1e-10) {print FOUT "              ,2,2,".( $L22).",".(-$L22s).",2,-2,$L22,$L22s\n";}
if (abs($L40)>1e-10)                  {print FOUT "              ,4,0,$L40,0\n";}
if (abs($L41*$L41+$L41s*$L41s)>1e-10) {print FOUT "              ,4,1,".(-$L41).",".( $L41s).",4,-1,$L41,$L41s\n";}
if (abs($L42*$L42+$L42s*$L42s)>1e-10) {print FOUT "              ,4,2,".( $L42).",".(-$L42s).",4,-2,$L42,$L42s\n";}
if (abs($L43*$L43+$L43s*$L43s)>1e-10) {print FOUT "              ,4,3,".(-$L43).",".( $L43s).",4,-3,$L43,$L43s\n";}
if (abs($L44*$L44+$L44s*$L44s)>1e-10) {print FOUT "              ,4,4,".( $L44).",".(-$L44s).",4,-4,$L44,$L44s\n";}
if (abs($L60)>1e-10)                  {print FOUT "              ,6,0,$L60,0\n";}
if (abs($L61*$L61+$L61s*$L61s)>1e-10) {print FOUT "              ,6,1,".(-$L61).",".( $L61s).",6,-1,$L61,$L61s\n";}
if (abs($L62*$L62+$L62s*$L62s)>1e-10) {print FOUT "              ,6,2,".( $L62).",".(-$L62s).",6,-2,$L62,$L62s\n";}
if (abs($L63*$L63+$L63s*$L63s)>1e-10) {print FOUT "              ,6,3,".(-$L63).",".( $L63s).",6,-3,$L63,$L63s\n";}
if (abs($L64*$L64+$L64s*$L64s)>1e-10) {print FOUT "              ,6,4,".( $L64).",".(-$L64s).",6,-4,$L64,$L64s\n";}
if (abs($L65*$L65+$L65s*$L65s)>1e-10) {print FOUT "              ,6,5,".(-$L65).",".( $L65s).",6,-5,$L65,$L65s\n";}
if (abs($L66*$L66+$L66s*$L66s)>1e-10) {print FOUT "              ,6,6,".( $L66).",".(-$L66s).",6,-6,$L66,$L66s\n";}
print FOUT "          };\n";

print FOUT << "EOF";
 )
 CNFG:
         $conf
   #i1   $nof_electrons
   #f1   $nof_electrons


 PARA:
    // if RED is not given, RED is equal to 0.8 and Rk parameters are scaled according to
    // RED !!!
    RED=1.0;

 EXEC:

   Mode =xas;
   Mag={$conf};
   Ninit=14;
   Range={2000,0,0.4,0.001};
   // this is to define that we want to calculate expectation value of Y40 and Y44
   //Goprt="A40;A44;A4m4";

 OPRT:
EOF
if($spdf=~/f/){print FOUT "      Rk(#i1 $conf $conf)={$F2,$F4,$F6};\n";}
else {print FOUT "      Rk(#i1 $conf $conf)={$F2,$F4};\n";}
print FOUT << "EOF";
   CAk(#i1 $conf)=Hcf;
    Zta(#i1 $conf)=$xi;

    // magnetic field in eV as muB*H
    Ba(#i1 $conf)=Bapplied;

 OPRT:
EOF
if($spdf=~/f/){print FOUT "      Rk(#i1 $conf $conf)={$F2,$F4,$F6};\n";}
else {print FOUT "      Rk(#i1 $conf $conf)={$F2,$F4};\n";}
print FOUT << "EOF";
   CAk(#i1 $conf)=Hcf;
    Zta(#i1 $conf)=$xi;

    // magnetic field in eV as muB*H
    Ba(#i1 $conf)=Bapplied;

    // here we say we want output of neutron spectra in dipole approximation for powder
    // thus we have to do 3 calculations: for 1) <i|Mx|f><f|Mx|i>, for 2) <i|My|f><f|My|i>
    // and for 3) <i|Mz|f><f|Mz|i>   and subsequently take  (sum of  the spectra) * 0.6667!!
 OPRT:

EOF
close FOUT;
copy("$file.xcardx","$file.xcardy");
copy("$file.xcardx","$file.xcardz");
open (FOUT,">>$file.xcardx");
print FOUT << "EOF";
    Ba(#i1 #f1 $conf $conf)={1,1,0,0};
 XEND:
 STOP:
EOF
open (FOUT,">>$file.xcardy");
print FOUT << "EOF";
    Ba(#i1 #f1 $conf $conf)={1,0,1,0};
 XEND:
 STOP:
EOF
open (FOUT,">>$file.xcardz");
print FOUT << "EOF";
    Ba(#i1 #f1 $conf $conf)={1,0,0,1};
 XEND:
 STOP:
EOF


print STDOUT << "EOF";
//*******************************************************
// xcard $file.xcard created ...
// you may now start an xtls calculation by
//        x900 $file.xcard
// or
//        x900G $file.xcard
//
// .... output is in xobj
//
// mind: - the directories xobjs and xwrk must exist
//       - midnight commander mc is not compatible with xmtl
//*******************************************************
EOF

exit;

sub usage() {

  print STDERR << "EOF";

    sipf2xcard: program to convert sipf file to xcard (for xtls)

    usage: sipf2xcard -h
           sipf2xcard sipffilename

     -h           : this (help) message
      sipffilename: filename of single ion parameter file

    output:

    sipffilename.xcard    xcard for XTLS

EOF

  exit;

}
# **********************************************************************************************
sub ic1ion()
{close Fin;
 unless (open(Fin,$file)){die "cannot open $file\n";}
 while($line=<Fin>){
                if($line=~/^\s*T\s*=/) {($T)=($line=~m|T\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                   }
 close Fin;

 system("ic1ionit -h $file > ic1ion.log");
 chdir "results";
 unless (open(Fin,"ic1ion.out")){die "cannot open file ic1ion.out\n";}
 while($line=<Fin>){
                if($line=~/^# Free ion configuration:/) {($spdf)=($line=~m|# Free ion configuration:\s*([spdf])|);
                           ($nof_electrons)=($line=~m|# Free ion configuration:\s*[spdf]\^(\d*)|);
                   print STDERR "Give the main quantum number n of the n$spdf^$nof_electrons configuration\n";
                   $n= <STDIN>;$n=~s/\n//;$conf=$n.$spdf;}
                if($line=~/^.*Bx\s*=/) {($Bx)=($line=~m|Bx\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*By\s*=/) {($By)=($line=~m|By\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*Bz\s*=/) {($Bz)=($line=~m|Bz\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*F\^2\s*=/) {($F2)=($line=~m|F\^2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*F\^4\s*=/) {($F4)=($line=~m|F\^4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*xi\s*=/)  {($xi)=($line=~m|xi\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*XI\s*=/)  {($xi)=($line=~m|XI\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*zeta\s*=/)  {($xi)=($line=~m|zeta\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*ZETA\s*=/)  {($xi)=($line=~m|ZETA\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L20\s*=/) {($L20)=($line=~m|L20\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L21\s*=/) {($L21)=($line=~m|L21\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L22\s*=/) {($L22)=($line=~m|L22\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L2-1\s*=/) {($L21s)=($line=~m|L2-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L2-2\s*=/) {($L22s)=($line=~m|L2-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L40\s*=/) {($L40)=($line=~m|L40\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L41\s*=/) {($L41)=($line=~m|L41\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L42\s*=/) {($L42)=($line=~m|L42\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L43\s*=/) {($L43)=($line=~m|L43\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L44\s*=/) {($L44)=($line=~m|L44\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-1\s*=/) {($L41s)=($line=~m|L4-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-2\s*=/) {($L42s)=($line=~m|L4-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-3\s*=/) {($L43s)=($line=~m|L4-3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-4\s*=/) {($L44s)=($line=~m|L4-4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
if($spdf=~/f/){
                if($line=~/^.*F\^6\s*=/) {($F6)=($line=~m|F\^6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L60\s*=/) {($L60)=($line=~m|L60\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L61\s*=/) {($L61)=($line=~m|L61\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L62\s*=/) {($L62)=($line=~m|L62\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L63\s*=/) {($L63)=($line=~m|L63\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L64\s*=/) {($L64)=($line=~m|L64\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L65\s*=/) {($L65)=($line=~m|L65\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L66\s*=/) {($L66)=($line=~m|L66\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-1\s*=/) {($L61s)=($line=~m|L6-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-2\s*=/) {($L62s)=($line=~m|L6-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-3\s*=/) {($L63s)=($line=~m|L6-3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-4\s*=/) {($L64s)=($line=~m|L6-4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-5\s*=/) {($L65s)=($line=~m|L6-5\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-6\s*=/) {($L66s)=($line=~m|L6-6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
              }
                    }
              close Fin;

 chdir "..";
 print STDERR "configuration is $conf\n";
 print STDERR "F^2=$F2 meV F^4=$F4 meV F^6=$F6 meV ZETA=$xi meV\n";
 print STDERR "Bx=$Bx T By=$By T Bz=$Bz T    T=$T K\n";

}
# **********************************************************************************************
sub icf1ion()
{close Fin;
 unless (open(Fin,$file)){die "cannot open $file\n";}
 while($line=<Fin>){
                if($line=~/^.*T\s*=/) {($T)=($line=~m|T\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                   }
 close Fin;

 system("icf1ionit $file > icf1ion.log");
 chdir "results";
 unless (open(Fin,"icf1ion.out")){die "cannot open file icf1ion.out\n";}
 while($line=<Fin>){
                if($line=~/^# Free ion configuration:/) {($spdf)=($line=~m|# Free ion configuration:\s*([spdf])|);
                           ($nof_electrons)=($line=~m|# Free ion configuration:\s*[spdf]\^(\d)|);
                   print STDERR "Give the main quantum number n of the n$spdf^$nof_electrons configuration\n";
                   $n= <STDIN>;$n=~s/\n//;$conf=$n.$spdf;}
                if($line=~/^.*Bx\s*=/) {($Bx)=($line=~m|Bx\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*By\s*=/) {($By)=($line=~m|By\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*Bz\s*=/) {($Bz)=($line=~m|Bz\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*xi\s*=/)  {($xi)=($line=~m|xi\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*XI\s*=/)  {($xi)=($line=~m|XI\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*zeta\s*=/)  {($xi)=($line=~m|zeta\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*ZETA\s*=/)  {($xi)=($line=~m|ZETA\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L20\s*=/) {($L20)=($line=~m|L20\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L21\s*=/) {($L21)=($line=~m|L21\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L22\s*=/) {($L22)=($line=~m|L22\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L2-1\s*=/) {($L21s)=($line=~m|L2-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L2-2\s*=/) {($L22s)=($line=~m|L2-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L40\s*=/) {($L40)=($line=~m|L40\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L41\s*=/) {($L41)=($line=~m|L41\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L42\s*=/) {($L42)=($line=~m|L42\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L43\s*=/) {($L43)=($line=~m|L43\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L44\s*=/) {($L44)=($line=~m|L44\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-1\s*=/) {($L41s)=($line=~m|L4-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-2\s*=/) {($L42s)=($line=~m|L4-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-3\s*=/) {($L43s)=($line=~m|L4-3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-4\s*=/) {($L44s)=($line=~m|L4-4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
if($spdf=~/f/){
                if($line=~/^.*L60\s*=/) {($L60)=($line=~m|L60\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L61\s*=/) {($L61)=($line=~m|L61\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L62\s*=/) {($L62)=($line=~m|L62\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L63\s*=/) {($L63)=($line=~m|L63\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L64\s*=/) {($L64)=($line=~m|L64\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L65\s*=/) {($L65)=($line=~m|L65\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L66\s*=/) {($L66)=($line=~m|L66\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-1\s*=/) {($L61s)=($line=~m|L6-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-2\s*=/) {($L62s)=($line=~m|L6-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-3\s*=/) {($L63s)=($line=~m|L6-3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-4\s*=/) {($L64s)=($line=~m|L6-4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-5\s*=/) {($L65s)=($line=~m|L6-5\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-6\s*=/) {($L66s)=($line=~m|L6-6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
              }
                    }
              close Fin;

 chdir "..";
 print STDERR "configuration is $conf\n";
$F2=1e10;# 1st and 2nd hunds rule -> strong Hee in comparison to Hso
$F4=1e10;
$F6=1e10;
 print STDERR "F^2=$F2 meV F^4=$F4 meV F^6=$F6 meV ZETA=$xi meV\n";
 print STDERR "Bx=$Bx T By=$By T Bz=$Bz T    T=$T K\n";

}
# **********************************************************************************************
sub so1ion()
{close Fin;
$spdf="f";$n=4;$conf=$n.$spdf;
 system("ic1ion $file > ic1ion.log");
 chdir "results";
 unless (open(Fin,"ic1ion.out")){die "cannot open file ic1ion.out\n";}
 while($line=<Fin>){
                if($line=~/^.*Ne\s*=/) {($nof_electrons)=($line=~m|Ne\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*Bx\s*=/) {($Bx)=($line=~m|Bx\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*By\s*=/) {($By)=($line=~m|By\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*Bz\s*=/) {($Bz)=($line=~m|Bz\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*T\s*=/) {($T)=($line=~m|T\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L20\s*=/) {($L20)=($line=~m|L20\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L21\s*=/) {($L21)=($line=~m|L21\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L22\s*=/) {($L22)=($line=~m|L22\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L2-1\s*=/) {($L21s)=($line=~m|L2-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L2-2\s*=/) {($L22s)=($line=~m|L2-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L40\s*=/) {($L40)=($line=~m|L40\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L41\s*=/) {($L41)=($line=~m|L41\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L42\s*=/) {($L42)=($line=~m|L42\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L43\s*=/) {($L43)=($line=~m|L43\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L44\s*=/) {($L44)=($line=~m|L44\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-1\s*=/) {($L41s)=($line=~m|L4-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-2\s*=/) {($L42s)=($line=~m|L4-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-3\s*=/) {($L43s)=($line=~m|L4-3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L4-4\s*=/) {($L44s)=($line=~m|L4-4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
if($spdf=~/f/){
                if($line=~/^.*L60\s*=/) {($L60)=($line=~m|L60\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L61\s*=/) {($L61)=($line=~m|L61\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L62\s*=/) {($L62)=($line=~m|L62\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L63\s*=/) {($L63)=($line=~m|L63\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L64\s*=/) {($L64)=($line=~m|L64\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L65\s*=/) {($L65)=($line=~m|L65\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L66\s*=/) {($L66)=($line=~m|L66\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-1\s*=/) {($L61s)=($line=~m|L6-1\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-2\s*=/) {($L62s)=($line=~m|L6-2\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-3\s*=/) {($L63s)=($line=~m|L6-3\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-4\s*=/) {($L64s)=($line=~m|L6-4\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-5\s*=/) {($L65s)=($line=~m|L6-5\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
                if($line=~/^.*L6-6\s*=/) {($L66s)=($line=~m|L6-6\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|); }
              }
                    }
              close Fin;

 chdir "..";
 print STDERR "configuration is $conf\n";
$F2=1e10;# 1st and 2nd hunds rule -> strong Hee in comparison to Hso
$F4=1e10;
$F6=1e10;
$xi=1e7; # Hso
 print STDERR "F^2=$F2 meV F^4=$F4 meV F^6=$F6 meV ZETA=$xi meV\n";
 print STDERR "Bx=$Bx T By=$By T Bz=$Bz T    T=$T K\n";

}


# **********************************************************************************************
# extracts variable from file
#
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:
#
# ($standarddeviation)=fileextract("sta","data.dat");
#
# ... it stores 0.24 in the variable $standarddeviation
#
sub fileextract {
             my ($variable,$filename)=@_;
             my $var="\Q$variable\E";
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
# **********************************************************************************************
# extracts number in col from line file
#
# for example somewhere in a file data.dat has 3 columns
#  # some header
#  0 2    3
#  2 0.24 3
#  4 2.3  2
#
# ($number)=extractnumber(3,2,"data.dat");
#
# ... it stores 2.3 in the variable $number
#
sub extractnumber {
             my ($line,$col,$filename)=@_;
             if(open (Fin,$filename))
             {my $i=0;while(($read=<Fin>)&&$i<$line){
                                                if($read=~/^\s*#/) {;}
                                                else{++$i;
                                         if($i==$line){@numbers=split(" ",$read);close Fin;return $numbers[$col-1];}
                                                    }

                                               }
               close Fin;
       	     }
             else
             {
             print STDERR "Warning: failed to read data file \"$filename\"\n";exit(1);
             }
            }
# **********************************************************************************************