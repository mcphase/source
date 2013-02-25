#! /usr/bin/perl

use FileHandle;

use PDL;

use PDL::Slatec;



$pi=3.1415;



unless ($#ARGV>=2) 

{print " program to create mcdisp.ini from mcphas.j for the calculation of neutron powder spectra\n\n";

print " usage: powdermagnon 0.3 2 0.1 10 0.5 30\n\n";

print " meaning take mcphas.j, generate a reflection list and put it to mcdisp.par \n";

print " 0.3 ....qmin   [1/A] minimal q vector\n";

print " 2   ....qmax   [1/A] maximal q vector\n";

print " 0.1 ....deltaq [1/A] stepwidth in q\n";

print " 10  ....number of steps in polar coordinate theta\n";

print " 0.5 ....Emin   [meV] minimal energy\n";

print " 30  ....Emax   [meV] maximal energy\n";

print "        (steps in fi fixed by dfi=dtheta*pi/(4*sin(theta))\n";



print "\n";

print "... then run mcdisp\n";

print "\n";

print ".... then restart powdermagnon by: powdermagnon -r results/mcdisp.qei 0.5 30 0.2 \n";

print " 0.5 ....Emin   [meV] minimal energy\n";

print " 30  ....Emax   [meV] maximal energy\n";

print " 0.5 ....deltaE [meV] energy stepwidth\n";

print " The powder average is generated - output is printed to the console (STDOUT)\n";

 exit 0;}


$ARGV[0]=~s/exp/essp/g;$ARGV[0]=~s/x/*/g;$ARGV[0]=~s/essp/exp/g;
my ($qmin) = eval $ARGV[0];
$ARGV[1]=~s/exp/essp/g;$ARGV[1]=~s/x/*/g;$ARGV[1]=~s/essp/exp/g;
my ($qmax) =eval $ARGV[1];



if($qmin=="-r")

{#read output file mcdisp.qei and create powder average 

 print "#Emin=".$Emin." meV   Emax=".$Emax." meV deltaE=".$deltaE." meV\n";

 print "#Ha[T] Hb[T] Hc[T] T[K] Q[A^-1] energy[meV] powderintensity [barn/sr/f.u.]   f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n";

    my $h = new FileHandle;

    open($h,$qmax);


$ARGV[2]=~s/exp/essp/g;$ARGV[2]=~s/x/*/g;$ARGV[2]=~s/essp/exp/g;
my ($Emin) = eval $ARGV[2];
$ARGV[3]=~s/exp/essp/g;$ARGV[3]=~s/x/*/g;$ARGV[3]=~s/essp/exp/g;
my ($Emax) = eval $ARGV[3];
$ARGV[4]=~s/exp/essp/g;$ARGV[4]=~s/x/*/g;$ARGV[4]=~s/essp/exp/g;
my ($deltaE) = eval $ARGV[4];

$n=int(($Emax-$Emin)/$deltaE);

my (@ints)=();$#ints=$n;
my (@intsbey)=();$#intsbey=$n;



$qold=0;$counter=1;

$qhold=0;$qkold=0;$qlold=0;

while(<$h>)

 {next if /^\s*#/;

  $line=$_;

  $line=~s/D/E/g;@numbers=split(" ",$line);

  $q=$numbers[7];

  if($q!=$qold&&$qold!=0){

     for($i=0;$i<=$n;++$i){

     $E=$Emin+($i-0.5)*$deltaE;$ints[$i]/=$counter;$intsbey[$i]/=$counter;

     print $numbers1[0]." ".$numbers1[1]." ".$numbers1[2]." ".$numbers1[3]." ";

     print $qold." ".$E." ".$ints[$i]." ".$intsbey[$i]."\n";}

     $counter=0;@ints=0;@intsbey=0;

               }

  if($numbers[8]<=$Emax&&$numbers[8]>=$Emin)
  {  $i=int(($numbers[8]-$Emin)/$deltaE);
     $ints[$i]+=$numbers[9];
     $intsbey[$i]+=$numbers[10];
  }
  $qold=$q;

  $qh=$numbers[4];

  $qk=$numbers[5];

  $ql=$numbers[6];

  if($qh!=$qhold||$qk!=$qkold||$ql!=$qlold){++$counter;}

  $qhold=$qh;

  $qkold=$qk;

  $qlold=$ql;
  @numbers1=@numbers;
 }

close $h;

exit;

}



my ($deltaq) = $ARGV[2];

my ($nn)=$ARGV[3];

my ($Emin) = $ARGV[4];

my ($Emax) = $ARGV[5];


$dtheta=$pi/$nn+0.000001;



my ($latt,$p) = getlattice("./mcphas.j");

my ($a,$b,$c,$nofatoms,$nofcomponents) = @{$latt};

 print "lattice: a=".$a." b=".$b." c=".$c."\n";

 print "creating file mcdisp.par\n";

    my $h = new FileHandle;

    open($h,">./mcdisp.par");

print $h "# file mcdisp.par created by program powdermagnon\n";
print $h "#!<--mcdisp.mcdisp.ini>\n";
print $h "#qmin=".$qmin." A^-1  qmax=".$qmax." A^-1 deltaq=".$deltaq." A^-1 no of thetasteps=".$nn."\n";
print $h "# lattice: a=".$a." b=".$b." c=".$c."\n";
print $h "#*********************************************************************\n";
print $h "# mcdisp - program to calculate the dispersion of magnetic excitations\n";
print $h "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n";
print $h "#*********************************************************************\n";
print $h "#\n";
print $h "# mcdisp calculates the neutron scattering cross section dsigma/dOmegadE' [barn/sr/meV/f.u.]\n";
print $h "#           f.u.=crystallogrpaphic unit cell (r1xr2xr3) for inelastic and diffuse scattering\n";
print $h "#\n";
print $h "# depending on what is kept constant it follows either kf or ki (1/A)\n";
print $h "#!kf=2\n";
print $h "# \n";
print $h "# emin and emax define the energy range in which neutron intensities are calculated\n";
print $h "# for full calculation of the dynamical susceptibility (option -r, inversion of the MF-RPA equation \n";
print $h "# for each point in Q-omega space) the minimum and maximum energy has to be given (energy stepwidth is \n";
print $h "# equal to the parameter epsilon given in the command line after -r)\n";
print $h "#\n";
print $h "#!emin=$Emin\n";
print $h "#!emax=$Emax\n";
print $h "#\n";
print $h "# optional parameter is extended_eigenvector_dimension\n";
print $h "# which is used to define, how many components of the\n";
print $h "# eigenvector should be in the ouput to file mcdisp.qee\n";
print $h "# (important for charge density movies)\n";
print $h "#! extended_eigenvector_dimension=3\n";
print $h "#\n";
print $h "#  a list of Q vectors with - generated by powdermagnon\n";
print $h "# h k l\n";
print $h "#\n";


#loop q

for ($q=$qmin;$q<=$qmax;$q+=$deltaq){ 



#loop sphere

 $h1=0;

 $h2=0;

 $h3=$c*$q/2/$pi;



 print $h "$h1 $h2 $h3\n";



 for ($theta=$dtheta;$theta<$pi;$theta+=$dtheta){ 

  $dfi=$dtheta*$pi/4/sin($theta);

 for ($fi=0;$fi<2*$pi;$fi+=$dfi){ 

 

 $h1=$a*$q*sin($theta)*cos($fi)/2/$pi;

 $h2=$b*$q*sin($theta)*sin($fi)/2/$pi;

 $h3=$c*$q*cos($theta)/2/$pi;



 print $h "$h1 $h2 $h3\n";

 }}

 $h1=0;

 $h2=0;

 $h3=-$c*$q/2/$pi;



 print $h "$h1 $h2 $h3\n";



}

    close $h; 

 print "please edit mcdisp.par, mcdisp.mf and run mcdisp\n";

 print "then restart powdermagnon by: powdermagnon -r results\mcdisp.qei Emin Emax deltaE\n";

   exit;

#-----------------------------------------------------------------------


# Get lattic data, reading it from file 

sub getlattice {

    my ($file) = @_;

    my $h = new FileHandle;

    my $n = 0;

    $nofcomponents=0;

#     my @xlist = ();

  # input data int piddle

  if(open($h,$file))

  {      while(<$h>)

     {#next if /^\s*#/;

      # detect a= b= c= ...

         
      if ($a==0){($a)=extract("a",$_);}
      if ($b==0){($b)=extract("b",$_);}
      if ($c==0){($c)=extract("c",$_);}

      if ($r1x==0){($r1x)=extract("r1a",$_);}
      if ($r1y==0){($r1y)=extract("r1b",$_);}
      if ($r1z==0){($r1z)=extract("r1c",$_);}
      if ($r2x==0){($r2x)=extract("r2a",$_);}
      if ($r2y==0){($r2y)=extract("r2b",$_);}
      if ($r2z==0){($r2z)=extract("r2c",$_);}
      if ($r3x==0){($r3x)=extract("r3a",$_);}
      if ($r3y==0){($r3y)=extract("r3b",$_);}
      if ($r3z==0){($r3z)=extract("r3c",$_);}

      if ($r1x==0){($r1x)=extract("r1x",$_);}
      if ($r1y==0){($r1y)=extract("r1y",$_);}
      if ($r1z==0){($r1z)=extract("r1z",$_);}
      if ($r2x==0){($r2x)=extract("r2x",$_);}
      if ($r2y==0){($r2y)=extract("r2y",$_);}
      if ($r2z==0){($r2z)=extract("r2z",$_);}
      if ($r3x==0){($r3x)=extract("r3x",$_);}
      if ($r3y==0){($r3y)=extract("r3y",$_);}
      if ($r3z==0){($r3z)=extract("r3z",$_);}

      if ($nofatoms==0){($nofatoms)=extract("nofatoms",$_);}
      if ($nofcomponents==0){($nofcomponents)=extract("nofcomponents",$_);}

      if (/^(#!|[^#])*nofneighbours\s*=\s*/){++$n;

                                   ($nofneighbours[$n])=extract("nofneighbours",$_);
                                   ($x[$n])=extract("da",$_);
                                   ($y[$n])=extract("db",$_);
                                   ($z[$n])=extract("dc",$_);
                                   if (/^.*\Qx=\E/){($x[$n])=extract("x",$_);}
				   if (/^.*\Qy=\E/){($y[$n])=extract("y",$_);}
				   if (/^.*\Qz=\E/){($z[$n])=extract("z",$_);}
				     ($sipffilename)=extractstring("sipffilename",$_);
                                     ($charge[$n])=extractfromfile("CHARGE",$sipffilename);
                                     if($charge[$n]==""){$charge[$n]=$sipffilename;}                                 
                                             #               print "$sipffilename  charge=".$charge[$n]."\n";
                                                           

				  }

     }

     close $h; 

     if ($n!=$nofatoms) {print STDOUT "Failed to read data file \"$file\": wrong number of atoms\n";

                         return undef;}

     if($nofcomponents==0) {$nofcomponents=3;}

     return ([$a,$b,$c,$nofatoms,$nofcomponents],pdl [[$r1x,$r1y,$r1z],[$r2x,$r2y,$r2z],[$r3x,$r3y,$r3z]]);

    } else {

	print STDOUT "Warning: failed to read data file \"$file\"\n";

	return undef;

    }

}





# **********************************************************************************************
# extracts number from string
# 
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
# 
# ... it stores 0.3 in the variable $standarddeviation
#
sub extract { 
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************

# **********************************************************************************************
# extracts string from string
# 
# ($standarddeviation)=extract("sta","sta=0.3");
# ($standarddeviation)=extract("sta","#!sta=0.3 # sta=0.2");  # i.e. comments are ignored unless followed by !
# 
# ... it stores 0.3 in the variable $standarddeviation
#
sub extractstring { 
             my ($variable,$string)=@_;
             $var="\Q$variable\E";
             $value="";
             if($string=~/^(#!|[^#])*\b$var\s*=\s*/) {($value)=($string=~m/^(?:#!|[^#])*\b$var\s*=\s*\b([^\n\s]+)[\s\n]/);
#             ($value)=($string=~m|^?:\#!|[^\#])*($var\s*=\s*([^\s]*))|);
             return $value;}
            }
# **********************************************************************************************




# **********************************************************************************************
# extracts number from file
# 
# for example somewhere in a file data.dat is written the text "sta=0.24"
# to extract this number 0.24 just use:         
#
# ($standarddeviation)=extractfromfile("sta","data.dat");
# 
# ... it stores 0.24 in the variable $standarddeviation
#
sub extractfromfile {
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
