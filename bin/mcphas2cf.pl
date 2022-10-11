#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

# rotmat - program to rotate Crystal Field parameters to the McPhase convention:

#   x||a, y||b, z||c --> x||c. y||a, z||b



use Getopt::Long;



# The following subroutines were shamelessly stolen from the Perl Cookbook

# by O'Reily (Recipe 2.14), and used for the matrix multiplications.



sub mmult {

    my ($m1,$m2) = @_;

    my ($m1rows,$m1cols) = matdim($m1);

    my ($m2rows,$m2cols) = matdim($m2);



    unless ($m1cols == $m2rows) {       # raise exception

        die "IndexError: matrices don't match: $m1cols != $m2rows";

    }



    my $result = [];

    my ($i, $j, $k);



    for $i (range($m1rows)) {

        for $j (range($m2cols)) {

            for $k (range($m1cols)) {

                $result->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];

            }

        }

    }

    return $result;

}



sub range { 0 .. ($_[0] - 1) }



sub veclen {

    my $ary_ref = $_[0];

    my $type = ref $ary_ref;

    if ($type ne "ARRAY") { die "$type is bad array ref for $ary_ref" }

    return scalar(@$ary_ref);

}



sub matdim {

    my $matrix = $_[0];

    my $rows = veclen($matrix);

    my $cols = veclen($matrix->[0]);

    return ($rows, $cols);

}



# End of stolen code!



sub usage() {

  print STDERR << "EOF";



    $0:

    Calculates a set of  crystal field  parameters for  Stevens equivalent

    operators in which Jx || a, Jy || b and Jz || c, from a set of crystal

    field parameters with Jx || c, Jy || a, and Jz || b. It is for use  by

    the program  mcphase to  convert parameters in its  convention back to

    the standard convention.



    The calculations are  done by means of matrix  multiplication based on

    the method of Buckmaster (phys. stat. sol. a, vol 13,  pp 9, 1972) and

    Rudowicz (J. Phys: Solid State Phys., vol 18, pp 1415, 1985).   



    usage: $0 [-h] [--help] 

              [-i input_file] [--input input_file]

	      [-o output_file] [--output output_file]

              [-v] [--verbose]



     -h          : this (help) message

     -i in_file  : input CF parameters file in cfield or mcphase formats

     -o out_file : output CF parameters file in mcphase format

     -v          : verbose mode. Will print out input parameters as read.



    if -i is omitted, the program will  assume the input CF parameters are

          given on the command line in the format: Bkq=x.xx,Bkq=x.xx, etc.

          e.g. $0 B20=0.21,B40=0.0005,B60=0.051,B66=0.626

    negative q parameters such as B_2^{-2},  are specified as:  B22S, with 

          an 'S' at the end, as per the McPhase convention.

    you may also  specify the ion type by a dding another  parameter after

          the CF parameters: e.g. $0 B20=0.21,B40=0.5 Pr3+



    if -o is omitted, the program prints the parameters to standard output. 



EOF

  exit;

}



# see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions



GetOptions("help"=>\$helpflag,

	   "input=s"=>\$input,

	   "output=s"=>\$output,

           "verbose"=>\$verbose);

	   

usage() if $helpflag;
print STDERR "#* $0 *\n";


if (!$input && !$ARGV[0]) { 

  print STDERR "$0: requires at least input file or CF parameters on command line.\n";

  usage();

  exit;

}



# Defines the inverse rotation matrices.



my $MR2 = [

            [    0, -1/2,    0,    0,    0],

            [    0,    0,    0,    1,    0],

            [    0,    0, -1/2,    0,  1/2],

            [    2,    0,    0,    0,    0],

            [    0,    0, -3/2,    0, -1/2]

          ];



my $MR4 = [

            [    0,  1/8,    0,  7/8,    0,    0,    0,    0,    0],

            [    0,    0,    0,    0,    0, -7/4,    0, -3/4,    0],

            [    0, -1/4,    0,  1/4,    0,    0,    0,    0,    0],

            [    0,    0,    0,    0,    0, -3/4,    0,  1/4,    0],

            [    0,    0,    0,    0,  3/8,    0, -1/8,    0,  1/8],

            [    1,    0, -1/2,    0,    0,    0,    0,    0,    0],

            [    0,    0,    0,    0,  5/2,    0, -1/2,    0, -1/2],

            [   -1,    0, -7/2,    0,    0,    0,    0,    0,    0],

            [    0,    0,    0,    0, 35/8,    0,  7/8,    0,  1/8],

          ];



my $MR6 = [

  [      0,  -1/32,      0, -11/32,      0, -33/32,      0,      0,      0,      0,      0,      0,      0],

  [      0,      0,      0,      0,      0,      0,      0,  33/16,      0,  33/16,      0,   5/16,      0],

  [      0,    1/8,      0,    3/8,      0,   -3/8,      0,      0,      0,      0,      0,      0,      0],

  [      0,      0,      0,      0,      0,      0,      0,  15/16,      0,  -1/16,      0,  -5/16,      0],

  [      0,  -5/32,      0,   9/32,      0,  -5/32,      0,      0,      0,      0,      0,      0,      0],

  [      0,      0,      0,      0,      0,      0,      0,    5/8,      0,   -3/8,      0,    1/8,      0],

  [      0,      0,      0,      0,      0,      0,  -5/16,      0,   1/16,      0,  -1/16,      0,   1/16],

  [    3/4,      0,   -1/2,      0,    1/4,      0,      0,      0,      0,      0,      0,      0,      0],

  [      0,      0,      0,      0,      0,      0,-105/32,      0,  17/32,      0,  -5/32,      0, -15/32],

  [   -5/8,      0,   -5/4,      0,    9/8,      0,      0,      0,      0,      0,      0,      0,      0],

  [      0,      0,      0,      0,      0,      0, -63/16,      0,   3/16,      0,  13/16,      0,   3/16],

  [    3/8,      0,   11/4,      0,   33/8,      0,      0,      0,      0,      0,      0,      0,      0],

  [      0,      0,      0,      0,      0,      0,-231/32,      0, -33/32,      0, -11/32,      0,  -1/32],

];



my %B=();



if ($input) {

  open (input_file, $input) or die "$0: cannot open $input for input CF parameters";



  while(<input_file>) {                                   # Selects out lines with crystal field parameters

    if ($_ =~ s/(B[0-9\ CcSs]*)\s*[=:]\s*([-\.\de]*)// ) { # () are groups which may be access with $1, $2 etc.

                                                          # * means match previous char any number of times.

      $ky = $1; $vl = $2;                                 # Parameters are of form Bkq = x.xx or Bkq : x.xx

      if ($vl != "") {                                    # \s matches whitespace characters.

        $ky =~ s/[cC ]//g;

        $B{$ky}=$vl;                                      # Assigns values of CF parameters to a hash.

      }



      if ($_ =~ s/(B[0-9\ CcSs]*)\s*[=:]\s*([-\.\de]*)// ) {

        $ky = $1; $vl = $2;                               # Second loop to get sine (-q) params in cfield

        if ($vl != "") {                                  #   input files.

          $ky =~ s/[ ]//g; $ky =~ s/s/S/g;

          $B{$ky}=$vl;                                    # Assigns values of CF parameters to a hash.

        }

      }  

    }

    if ($_ =~ /I[Oo][Nn].*[=:]([ \w]\w*\+)/ ) { $Ion = $1; }

  }

  close (input_file);

  print "\n";

}

else {

  @par = split (/,/, $ARGV[0]);

  foreach (@par) { 

    /(B[0-9S]*)=([-\.\d]*)/;                           # Parameters are in the form: Bkq=x.xx,Bkq=-x.xx,etc.

    $B{$1} = $2;

  }

  if ($ARGV[1]) { $Ion = $ARGV[1]; }

}



my $B2 = [ [$B{"B22S"}],[$B{"B21S"}],[$B{"B20"}],[$B{"B21"}],[$B{"B22"}] ];

my $R2 = mmult($MR2, $B2);

my $B4 = [ [$B{"B44S"}],[$B{"B43S"}],[$B{"B42S"}],[$B{"B41S"}],[$B{"B40"}],[$B{"B41"}],[$B{"B42"}],[$B{"B43"}],[$B{"B44"}] ];

my $R4 = mmult($MR4, $B4);

my $B6 = [ [$B{"B66S"}],[$B{"B65S"}],[$B{"B64S"}],[$B{"B63S"}],[$B{"B62S"}],[$B{"B61S"}],[$B{"B60"}],[$B{"B61"}],[$B{"B62"}],[$B{"B63"}],[$B{"B64"}],[$B{"B65"}],[$B{"B66"}] ];

my $R6 = mmult($MR6, $B6);



my %Bo2 = (0, "B22S", 1, "B21S", 2, "B20 ", 3, "B21 ", 4, "B22 ");

my %Bo4 = (0, "B44S", 1, "B43S", 2, "B42S", 3, "B41S", 4, "B40 ", 5, "B41 ", 6, "B42 ", 7, "B43 ", 8, "B44 ");

my %Bo6 = (0, "B66S", 1, "B65S", 2, "B64S", 3, "B63S", 4, "B62S", 5, "B61S", 6, "B60 ", 7, "B61 ", 8, "B62 ", 9, "B63 ", 10, "B64 ", 11, "B65 ", 12, "B66 ");



if (!$output) {

  if ($verbose) {

    print "Input rotated crystal field parameters (x||c, y||a, z||b):\n";

    while ( ($key, $value) = each %B) {

      print "$key = $value\n";

    }

    if ($Ion) { print "\nIontype is ".$Ion."\n"; }

    print "\nOutput crystal field parameters (x||a, y||b, z||c):\n";

  }

  for $i (0 .. 5) { if ($_ = $R2->[$i][0]) { print "$Bo2{$i} = $_\n"; } }

  for $i (0 .. 9) { if ($_ = $R4->[$i][0]) { print "$Bo4{$i} = $_\n"; } }

  for $i (0 .. 13) { if ($_ = $R6->[$i][0]) { print "$Bo6{$i} = $_\n"; } }

}

else {

  open (outfile, ">$output") or die "$0: cannot open $output for output.";

  

  print outfile << "EOF";

#!cfield

#<!--mcphase.sipf-->

EOF

  

  $Ion =~ s/ //;

  print outfile " IONTYPE=$Ion\n";

  for $i (0 .. 5) { if ($_ = $R2->[$i][0]) { print outfile " $Bo2{$i} = $_\n"; } }

  for $i (0 .. 9) { if ($_ = $R4->[$i][0]) { print outfile " $Bo4{$i} = $_\n"; } }

  for $i (0 .. 13) { if ($_ = $R6->[$i][0]) { print outfile " $Bo6{$i} = $_\n"; } }

  print outfile << "EOF";





#   D = 2 * pi / Q

#   s = 1 / 2 / D: sintheta = lambda * s



#Nd3+ magnetische formfaktoren nach international tables

#     ff(i,1...7)= <j0(kr)>-terms A,a,B,b,C,c,D

FFj0A=0.2953 FFj0a= 17.6846 FFj0B=0.2923 FFj0b=6.7329 FFj0C=0.4313 FFj0c=5.3827 FFj0D=-0.0194

#     ff(i,8..14)= <j2(kr)>-terms A,a,B,b,C,c,D

FFj2A=0.9809 FFj2a=18.0630 FFj2B = 1.8413 FFj2b= 7.7688 FFj2C=0.9905  FFj2c=2.8452 FFj2D=0.0120



#  ....

#   j0 = ff(1) * EXP(-ff(2) * s * s) + ff(3) * EXP(-ff(4) * s * s)

#   j0 = j0 + ff(5) * EXP(-ff(6) * s * s) + ff(7)

#   j2 = ff(8) * s * s * EXP(-ff(9) * s * s) + ff(10) * s * s * EXP(-ff(11) * s * s)

#   j2 = j2 + ff(12) * s * s * EXP(-ff(13) * s * s) + s * s * ff(14)

#   FQ = (j0 + j2 * (2 / gJ - 1))





# Debey-Waller Factor

DWF=0

# debeywallerfactor = EXP(-2 * DWF *s*s)

EOF



  close (outfile);

}

