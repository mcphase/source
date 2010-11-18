#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

# gcalc.pl

# Calculates the gap size of a capacitance dilatometer from its capacitance and physical

# dimensions, taking into account the thermal expansion of the silver casing and sapphire

# screw bearing.

#

# By Martin Rotter (GCALC.BAS). Ported to perl by Duc Le.

#

# Reference: Rotter et al. Rev. Sci. Instr. vol 69. pp 2742



use constant PI    => 4 * atan2(1, 1);

use Getopt::Long;

use Math::Interpolate qw(derivatives constant_interpolate

                          linear_interpolate robust_interpolate);



# ------------------------------ Parses input variables ------------------------------#



sub usage() {

  print STDERR << "EOF";



    $0:

    Calculates  the  gap  size  of  a  capacitance  dilatometer  given  its

    capacitance  and physical dimensions,  taking into  account the thermal

    expansion of its silver construction and saphire bearings.

    

    Based on GCALC.BAS by Martin Rotter. 

    Ref: Rotter et al. Rev. Sci. Instr. vol 69. pp 2742



    usage: $0 [-h] [--help] 

              [--c0 initial_capacitance_in_pF]

	      [--b centre_pivot_distance_in_mm] 

	      [--ri inner_radius_in_mm]

	      [--ra outer_radius_in_mm]

	      [--input input_datafile]

	      [--temp temperature]

	      [--silver silver_expansion_filename]



    The data file must have the  capacitance  in the third (3) column,  and

    the temperature in the second  (2)  column.  Alternatively, if data was

    taken at constant temperature,  the option --temp may be specified with

    this temperature (e.g. for magnetostriction).



    If the options  are omitted,  and  instead  a  string  is  given on the

    command,  it will be assumed that the parameters are given in the above

    order. E.g.



    $0 c0 b ri ra input_datafile silver_expansion_file



    Otherwise if the parameters are omitted entirely, the default values of

    c0 = 4pF ; b = 9.794mm ; ri = 2.33 ; ra = 6.25 ; will be assumed (which

    corresponds to the W cell). Additionly the recommended values of silver

    expansion from Touloukian et al. will be used. E.g.



    $0 input_datafile



EOF

  exit;

}



# Initial parameters defaults

$c0 = 4;    # Initial capacitance (pF)

$b = 9.794; # Distance between centre of capacitor to pivot (mm)

$ri = 2.33; # Inner plate radius (mm)

$ra = 6.25; # Outer plate radius (mm)



@pars = @ARGV;

if (join('',@pars) =~ /\-[\-cbri]/) {

  # see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions

  GetOptions("help"=>\$helpflag,

  	     "c0=f"=>\$c0,

	     "b=f"=>\$b,

	     "ri=f"=>\$ri,

	     "ra=f"=>\$ra,

	     "input=s"=>\$inputfile,

	     "temp=f"=>\$temp,

	     "field=f"=>\$field,

	     "silver=f"=>\$silverfile);

	   

  usage() if $helpflag;



  if (!$inputfile && !$ARGV[0]) { 

    print STDERR "$0: requires at least input file on commandline.\n";

    usage();

    exit;

  }



  if ($ARGV[0]) { $inputfile = $ARGV[0]; }

  #print "c0=$c0\tb=$b\tri=$ri\tra=$ra\tinput=$inputfile\n";

}

else {

  if (!$ARGV[0]) { 

    print STDERR "$0: no c0\n"; usage(); exit;

  }

  else { 

    if ($ARGV[0] =~ /\d+/) { 

      $c0 = $ARGV[0]; 

      if (!$ARGV[1]) { $inputfile = $ARGV[0]; $c0 = 4; break }  else { $b  = $ARGV[1]; }

      if (!$ARGV[2]) { print STDERR "$0: no ri\n";usage();exit;} else { $ri = $ARGV[2]; }

      if (!$ARGV[3]) { print STDERR "$0: no ra\n";usage();exit;} else { $ra = $ARGV[3]; }

      if (!$ARGV[4]) { print STDERR "$0: no input\n";usage();exit;} else { $inputfile = $ARGV[4]; } 

      if ($ARGV[5]) { $silverfile = $ARGV[5]; } 

    }

    else { 

      $inputfile = $ARGV[0]; 

    }

  }



  #print "c0=$c0\tb=$b\tri=$ri\tra=$ra\tinput=$inputfile\tsilver=$silverfile\n";

}



# -------------------------------- Inputs data and lit -------------------------------#



if ($silverfile) {

  open(FHSilver,$silverfile);

  #open(FHSilver,'silver.touloukian');

  while(<FHSilver>) {                      # Reads in the Rc Cernox calibratin file

    chomp;                                 # Cuts out newlines

    if ($_ =~/^[\.0-9\s+]+/) {             # Ignores non-digit (non-data) values

      $_ =~ s/^\s+//;                      # removes whitespace at start of string.

      $_ =~ s/\s+$//;                      # removes whitespace at end of string.

      @Line = split(/\s+/);

      push @AgT, $Line[0];

      push @Agdll, $Line[1];

      #print "$AgT[-1] $Agdll[-1]\n";

    } 

  }

  close(FHSilver);

}

else {

# Data from pp298, Thermophysical properties of matter, vol. 12: Metallic Elements and

#   Alloys, YS Touloukian, RK Kirby, RE Taylor and PD Desai, Plenum Press 1981.

# error is +/- 3% over all temperature range

# approximated by: \Delta L/L_0 = -0.515 + 1.647e-3*T + 3.739e-7*T**2 + 6.283e-11*T**3

# T	\Delta L/L_0	alpha_x_10^6

  @AgT = (5,25,50,100,200,293,400,500,600,700,800,900,1000,1200);

  @Agdll = (-0.00410,-0.00409,-0.00396,-0.00336,-0.00171,0,0.00207,0.00409,0.00620,0.00842,0.01072,0.01312,0.01568,0.02110);

}



#open(FHSaphire,'saphire.lit');

#while(<FHSaphire>) {                     # Reads in the Rc Cernox calibratin file

#  chomp;                                 # Cuts out newlines

#  if ($_ =~/^[\.0-9\s+]+/) {             # Ignores non-digit (non-data) values

#    $_ =~ s/^\s+//;                      # removes whitespace at start of string.

#    $_ =~ s/\s+$//;                      # removes whitespace at end of string.

#    @Line = split(/\s+/);

#    push @sapT, $Line[0];

#    push @sapdll, $Line[1];

#    #print "$sapT[-1] $sapdll[-1]\n";

#  } 

#}

#close(FHSaphire);



# Value of the thermal expansion of saphire from literature

# Source: G.K. White, Thermochim. Acta. vol. 218, pp. 83 (1993) 

#@sapT = (20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 293, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000);

#@sapdll = (-635E-6, -635E-6, -635E-6, -634E-6, -633E-6, -631E-6, -627E-6, -622E-6, -615E-6, -595E-6, -565E-6, -524E-6, -473E-6, -410E-6, -212E-6, 0, 37E-6, 324E-6, 642E-6, 986E-6, 1350E-6, 1730E-6, 2125E-6, 2530E-6, 2940E-6, 3790E-6, 4665E-6, 5560E-6, 6485E-6, 7430E-6, 8400E-6, 9400E-6, 10420E-6, 11460E-6, 12520E-6, 13600E-6, 14690E-6, 15810E-6);



open(FHDatafile,$inputfile);

while(<FHDatafile>) {

  if ($_ =~/^[\.0-9\s+]+/) {             # Ignores non-digit (non-data) values

    if ($_ =~/^\s+$/) { next; }          # Ignores if line is empty

    chomp;

    @Line = split(/\s+/);

    if ($field) { push @H, $field }

    else { push @H, shift(@Line); }

    if ($temp) { push @T, $temp }

    else { push @T, shift(@Line); }

    push @C, shift(@Line);

    push @otherdat, join("\t",@Line);

    #print "H=$H[-1]\tT=$T[-1]\tC=$C[-1]\t$otherdat[-1]\n";

  }

  else {

    print "#$0:Oldheader: $_";

  }

}

close(FHDatafile);



print "#Field(T)\tTemperature(K)\tGap(mm)\tother_data(see_above_header)\n";



# -------------------------------- Sets Initial values -------------------------------#



# Physical constants. Taken from NIST Reference on Constants, Units, and 

# Uncertainty, http://physics.nist.gov/cuu/Constants/

#$eps0  = 8.854187817e-12;       # F/m - permitvity of free space

$eps0  = 8.854187817e-3;         # pF/mm - permitvity of free space



# Calculates the areas of the plates

$Aa0 = PI*$ra**2;                # Outer

$Ai0 = PI*$ri**2;                # Inner



$k0 = $eps0 * ($Aa0-$Ai0)/$c0;   # Initial pivot distance k(T_0)

$ds = 0.8;                       # thickness of saphire washer in mm.



for $iloop (0..$#H) {

  # Interpolate to find dl/l of silver at present temperature

  $dllAg = robust_interpolate($T[$iloop], \@AgT, \@Agdll);

  #$dllSa = robust_interpolate($T[$iloop], \@sapT, \@sapdll);



  #$kT = k0 * (1 + dllAg)**2;     # pivot distance k(T)

  #$kT = $k0 + 2*$ds*($dllAg-$dllSa);

  $kT = $k0;

  #$AaT = $Aa0 * (1 + $dllAg)**2; # outer plate area A_a(T)

  #$AiT = $Ai0 * (1 + $dllAg)**2; # inner plate area A_i(T)

  $AaT = $Aa0 * (1 + $dllAg);    # outer plate area A_a(T)

  $AiT = $Ai0 * (1 + $dllAg);    # inner plate area A_i(T)



  $ccalc = 1000;                 # convergence iteration limit



  # Initialises gap size d

  $d = ($AaT-$AiT) * $eps0 / $C[$iloop]; 

  $dstep = 0.0001;



  if ( $d < ($kT*$ra/($b+$ra)) ) {

    $d = 1.1*$kT*$ra / ($b+$ra);

  }



  $oldsign = 1;

  while ( abs($C[$iloop]-$ccalc)>1e-9 & ($dstep/$d)>1e-15 ) {

    $gamma_a = ($ra/$b) * ($kT-$d)/$d;

    $gamma_i = ($ri/$b) * ($kT-$d)/$d;

    $ccalc = (2*$eps0/$d) * ( ($AaT*(1-sqrt(1-$gamma_a**2))/$gamma_a**2) );

    $ccalc+= (2*$eps0/$d) * (-($AiT*(1-sqrt(1-$gamma_i**2))/$gamma_i**2) );

    # The <=> operator returns -1, 0, 1 depending on whether left is less, 

    # equal to or more than the right hand side.

    $newsign = ($ccalc-$C[$iloop]) <=> 0;

    $d = $d + $newsign*$dstep;

    if ($oldsign!=$newsign || $newsign==0) {

      $oldsign = $newsign;

      $dstep = $dstep/10;

    }

  }

  

  $gap[$iloop] = $d;



  print "$H[$iloop]\t$T[$iloop]\t$gap[$iloop]\t$otherdat[$iloop]\n";

}

