#!/usr/bin/perl
# rotmat - program to rotate Crystal Field parameters about an azimuthal angle fi and polar angle theta

BEGIN{@ARGV=map{glob($_)}@ARGV}

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
    Rotates a set of  crystal field  parameters for  Stevens equivalent
    operators by an azimuthal angle fi about the original z axis and
    a polar angle theta about the new y axis. A right hand axis system is 
    assumed and a positive rotation is one which advances a right-hand 
    screw in a positive direction along the axis.
    The calculations are  done by means of matrix  multiplication based on
    the method of Buckmaster (phys. stat. sol. a, vol 13,  pp 9, 1972) and
    Rudowicz (J. Phys: Solid State Phys., vol 18, pp 1415, 1985).   
    usage: $0 [-h] [--help] 
              [-i input_file] [--input input_file]
	      [-o output_file] [--output output_file]
              [-v] [--verbose] [-th theta] [-fi fi] [CF parameters]
     -h          : this (help) message
     -i in_file  : input CF parameters file in so1ion formats
     -o out_file : output CF parameters file in mcphase format
     -v          : verbose mode. Will print out input parameters as read.
     -th	 : polar angle theta in degrees
     -fi         : azimuthal angle fi in degrees
    if -i is omitted, the program will  assume the input CF parameters are
          given on the command line in the format: Bkq=x.xx,Bkq=x.xx, etc.
          e.g. $0 "B20=0.21,B40=0.0005,B60=0.051,B66=0.626"
          negative q parameters such as B_2^{-2},  
          are specified as:  B22S, with an 'S' at the end, as per the 
          McPhase convention. you may also  specify the ion type by adding 
          another parameter after the CF parameters: e.g. 
	  $0 "B20=0.21,B40=0.5" Pr3+
	  Note that the quotation marks are required in Windows, but not
	  in Linux.
    if -o is omitted, the program prints the parameters to standard output. 
EOF
  exit;

}
# see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions

$theta=0; $fi=0;
GetOptions("help"=>\$helpflag,
	   "input=s"=>\$input,
	   "output=s"=>\$output,
           "verbose"=>\$verbose,
	   "th=s"=>\$theta,
	   "fi=s"=>\$fi);
$theta=~s/x/*/g;$theta= eval $theta;
$fi=~s/x/*/g;$fi=eval $fi;
usage() if $helpflag;

if (!$input && !$ARGV[0]) { 
  print STDERR "$0: requires at least input file or CF parameters on command line.\n";
  usage();
  exit;
}

my %B=();
if ($input) {
  open (input_file, $input) or die "$0: cannot open $input for input CF parameters";
  while(<input_file>) {                                   # Selects out lines with crystal field parameters
  if ($_=~/^\s*#/) {}
  else{  
   if ($_ =~ s/(B[0-9\ CcSs]+)\s*[=:]\s*([-\.\de]*)// ) { # () are groups which may be access with $1, $2 etc.
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

print "theta=".$theta." deg  fi=".$fi." deg\n\n";
$PI=3.14159265;

# Defines the rotation matrices 
$CF=cos($fi/180*$PI);
$SF=sin($fi/180*$PI);
$CF2=cos(2*$fi/180*$PI);
$SF2=sin(2*$fi/180*$PI);
$CF3=cos(3*$fi/180*$PI);
$SF3=sin(3*$fi/180*$PI);
$CF4=cos(4*$fi/180*$PI);
$SF4=sin(4*$fi/180*$PI);
$CF5=cos(5*$fi/180*$PI);
$SF5=sin(5*$fi/180*$PI);
$CF6=cos(6*$fi/180*$PI);
$SF6=sin(6*$fi/180*$PI);
$CT=cos($theta/180*$PI);
$ST=sin($theta/180*$PI);
$CT2=cos(2*$theta/180*$PI);
$ST2=sin(2*$theta/180*$PI);
# in Rudowicz (J. Phys: Solid State Phys., vol 18, pp 1415, 1985)
# there is givien Slmm' such that
# Olm(J)=sum_m' Slmm' Olm'(J')
# ... it follows
# Blm'_new = sum_m ( Blm_orig * Slmm')
# i.e. in Slmm' transformation matrices for the operators the m and m' has to be exchanged
#      to get the transformation of the crystal field parameters
my $MR2 = [
            [   $CT*$CF2                  , -$ST*$CF / 2  ,0                     ,$ST*$SF / 2  ,  -$CT*$SF2              ],  
            [    2*$ST*$CF2               ,$CT*$CF        ,0                     ,-$CT*$SF     ,  -2*$ST*$SF2            ],  
            [   $ST**2*$SF2 / 2           ,$ST2*$SF / 4   ,( 3*$CT**2 - 1 ) / 2  ,$ST2*$CF / 4 ,$ST**2*$CF2 / 2          ],  
            [   $ST2*$SF2                 ,$CT2*$SF       , -3*$ST2              ,$CT2*$CF     ,$ST2*$CF2                ],  
            [   $SF2 * ( $CT**2 + 1 ) / 2 ,-$ST2*$SF / 4  , 3*$ST**2 / 2         ,-$ST2*$CF / 4,$CF2 *( $CT**2 + 1 ) / 2 ]
          ];
my $MR4 = [
            [  $CF4*$CT *( $CT**2 + 1 ) / 2        ,-$ST*$CF3 * ( 3*$CT**2 + 1 ) / 8      ,  7*$ST**2*$CF2*$CT / 4               ,-7*$ST**3*$CF / 8                     ,0                              ,7*$ST**3*$SF / 8,- 7*$ST**2*$SF2*$CT / 4,$ST*$SF3 * ( 3*$CT**2 + 1 ) / 8,  - $SF4*$CT *( $CT**2 + 1 ) / 2    ],
            [  $ST*$CF4 * ( 3*$CT**2 + 1 )        , $CF3*$CT * ( 9*$CT**2 - 5 ) / 4      ,-7*$ST*$CF2 * ( 3*$CT**2 - 1 ) / 2    , 21*$ST**2*$CF*$CT / 4                 ,0                              ,- 21*$ST**2*$SF*$CT / 4,7*$ST*$SF2 * ( 3*$CT**2 - 1 ) / 2,- $SF3*$CT * ( 9*$CT**2 - 5 ) / 4,  - $ST*$SF4 * ( 3*$CT**2 + 1 )    ],
            [  $ST**2*$CF4*$CT                    , $ST*$CF3 * ( 3*$CT**2 - 1 ) / 4      ,  $CF2*$CT * ( 7*$CT**2 - 5 ) / 2     ,-$ST*$CF * ( 7*$CT**2 - 1 ) / 4        ,0                              ,$ST*$SF * ( 7*$CT**2 - 1 ) / 4,- $SF2*$CT * ( 7*$CT**2 - 5 ) / 2,- $ST*$SF3 * ( 3*$CT**2 - 1 ) / 4,  - $ST**2*$SF4*$CT],
            [  $ST**3*$CF4                       , 3*$ST**2*$CF3*$CT / 4                ,  $ST*$CF2 *( 7*$CT**2 - 1 ) / 2       , $CF*$CT *( 7*$CT**2 - 3 ) / 4          ,0                              ,- $SF*$CT *( 7*$CT**2 - 3 ) / 4,- $ST*$SF2* ( 7*$CT**2 - 1 ) / 2,- 3*$ST**2*$SF3*$CT / 4,  - $ST**3*$SF4],
            [$ST**4*$SF4 / 8                     ,$ST**3*$SF3*$CT / 8                   ,$ST**2*$SF2 * ( 7*$CT**2 - 1 ) / 8    ,$ST2*$SF *( 7*$CT**2 - 3 ) / 16         ,( 35*$CT**4 - 30*$CT**2 + 3 ) / 8,$ST2*$CF *( 7*$CT**2 - 3 ) / 16,$ST**2*$CF2 * ( 7*$CT**2 - 1 ) / 8,$ST**3*$CF3*$CT / 8,  $ST**4*$CF4 / 8],
            [$ST**3*$SF4*$CT                      ,$ST**2*$SF3 *( 4*$CT**2 - 1 ) / 4      ,$ST2*$SF2 * ( 7*$CT**2 - 4 ) / 2      ,$SF * ( 28*$CT**4 - 27*$CT**2 + 3 ) / 4,- 5*$ST2 *( 7*$CT**2 - 3 ) / 2    ,$CF * ( 28*$CT**4 - 27*$CT**2 + 3 ) / 4,$ST2*$CF2 * ( 7*$CT**2 - 4 ) / 2,$ST**2*$CF3 *( 4*$CT**2 - 1 ) / 4,  $ST**3*$CF4*$CT],
            [$ST**2*$SF4 * ( $CT**2 + 1 ) / 2     ,$ST*$SF3*$CT**3 / 2                   ,$SF2 * ( 7*$CT**4 - 6*$CT**2 + 1 ) / 2,- $ST2*$SF * ( 7*$CT**2 - 4 ) / 4      ,5*$ST**2 * ( 7*$CT**2 - 1 ) / 2  ,- $ST2*$CF * ( 7*$CT**2 - 4 ) / 4,$CF2 * ( 7*$CT**4 - 6*$CT**2 + 1 ) / 2,$ST*$CF3*$CT**3 / 2,  $ST**2*$CF4 * ( $CT**2 + 1 ) / 2],
            [$ST2*$SF4 * ( $CT**2 + 3 ) / 2       ,$SF3 * ( 4*$CT**4 + 3*$CT**2 - 3 ) / 4,- 7*$ST*$SF2*$CT**3                   ,7*$ST**2*$SF * ( 4*$CT**2 - 1 ) / 4    ,- 35*$ST**3*$CT                  ,7*$ST**2*$CF * ( 4*$CT**2 - 1 ) / 4,- 7*$ST*$CF2*$CT**3,$CF3 * ( 4*$CT**4 + 3*$CT**2 - 3 ) / 4,  $ST2*$CF4 * ( $CT**2 + 3 ) / 2],
            [$SF4 *  ( $CT**4 + 6*$CT**2 + 1 ) / 8,- $ST2*$SF3 *( $CT**2 + 3 ) / 16       ,7*$ST**2*$SF2 * ( $CT**2 + 1 ) / 8    ,- 7*$ST**3*$SF*$CT / 8                 ,35*$ST**4 / 8                   ,- 7*$ST**3*$CF*$CT / 8, 7*$ST**2*$CF2 * ( $CT**2 + 1 ) / 8,- $ST2*$CF3 *( $CT**2 + 3 ) / 16,  $CF4 *  ( $CT**4 + 6*$CT**2 + 1 ) / 8]
          ];
my $MR6 = [
  [ $CT*$CF6 * ( 3*$CT**4 + 10*$CT**2 + 3 ) / 16  ,-$ST*$CF5* ( 5*$CT**4 + 10*$CT**2 + 1 ) / 32   , 11*$ST**2*$CT*$CF4 *($CT**2 + 1 ) / 8     ,-11*$ST**3*$CF3* ( 3*$CT**2 + 1 ) / 32         , 33*$ST**4*$CT*$CF2 / 16                         ,-33*$ST**5*$CF / 32                                 ,      0                                           ,33*$ST**5*$SF / 32                        ,- 33*$ST**4*$CT*$SF2 / 16                    ,11*$ST**3*$SF3* ( 3*$CT**2 + 1 ) / 32            ,- 11*$ST**2*$CT*$SF4* ($CT**2 + 1 ) / 8           ,$ST*$SF5* ( 5*$CT**4 + 10*$CT**2 + 1 ) / 32,- $CT*$SF6 * ( 3*$CT**4 + 10*$CT**2 + 3 ) / 16],
  [ 3*$ST*$CF6* ( 5*$CT**4 + 10*$CT**2 + 1 ) / 8  , $CT*$CF5 *( 25*$CT**4 + 10*$CT**2 - 19 ) / 16 ,-11*$ST*$CF4 * ( 5*$CT**4 - 1 ) / 4        , 33*$ST**2*$CT*$CF3 * ( 5*$CT**2 - 1 ) / 16    ,-33*$ST**3*$CF2 *( 5*$CT**2 - 1 ) / 8            , 165*$ST**4*$CT*$CF / 16                            ,      0                                           ,- 165*$ST**4*$CT*$SF / 16                 ,33*$ST**3*$SF2 *( 5*$CT**2 - 1 ) / 8         , - 33*$ST**2*$CT*$SF3 * ( 5*$CT**2 - 1 ) / 16    ,11*$ST*$SF4 * ( 5*$CT**4 - 1 ) / 4                ,- $CT*$SF5 *( 25*$CT**4 + 10*$CT**2 - 19 ) / 16,- 3*$ST*$SF6* ( 5*$CT**4 + 10*$CT**2 + 1 ) / 8],
  [ 3*$ST**2*$CT*$CF6 *( $CT**2 + 1 ) / 4         , $ST*$CF5 * ( 5*$CT**4 - 1 ) / 8               , $CT*$CF4* ( 11*$CT**4 - 10*$CT**2 + 1 )/2 ,-3*$ST*$CF3 * ( 11*$CT**4 - 8*$CT**2 + 1 ) / 8 , 3*$ST**2*$CT*$CF2 * ( 11*$CT**2 - 5 ) / 4       ,-3*$ST**3*$CF * ( 11*$CT**2 - 1 ) / 8               ,      0                                           ,3*$ST**3*$SF * ( 11*$CT**2 - 1 ) / 8      , - 3*$ST**2*$CT*$SF2 * ( 11*$CT**2 - 5 ) / 4 ,3*$ST*$SF3 * ( 11*$CT**4 - 8*$CT**2 + 1 ) / 8    ,- $CT*$SF4* ( 11*$CT**4 - 10*$CT**2 + 1 ) / 2     ,- $ST*$SF5 * ( 5*$CT**4 - 1 ) / 8,- 3*$ST**2*$CT*$SF6 *( $CT**2 + 1 ) / 4 ],
  [ 5*$ST**3*$CF6 * ( 3*$CT**2 + 1 ) / 8          ,  5*$ST**2*$CT*$CF5* ( 5*$CT**2 - 1 ) / 16     , 5*$ST*$CF4 *( 11*$CT**4 - 8*$CT**2 + 1)/4 , $CT*$CF3 * ( 165*$CT**4 - 206*$CT**2 +57)/16  ,-3*$ST*$CF2* ( 55*$CT**4 - 42*$CT**2 + 3 ) / 8   , 15*$ST**2*$CT*$CF * ( 11*$CT**2 - 3 ) / 16         ,      0                                           ,- 15*$ST**2*$CT*$SF *(11*$CT**2-3)/16     ,3*$ST*$SF2* ( 55*$CT**4 - 42*$CT**2 + 3 )/8  ,- $CT*$SF3 * ( 165*$CT**4 - 206*$CT**2 + 57 ) /16,- 5*$ST*$SF4 *( 11*$CT**4 - 8*$CT**2 + 1 ) / 4    , - 5*$ST**2*$CT*$SF5* ( 5*$CT**2 - 1 ) / 16,- 5*$ST**3*$SF6 * ( 3*$CT**2 + 1 ) / 8],
  [ 15*$ST**4*$CT*$CF6 / 16                       , 5*$ST**3*$CF5 *( 5*$CT**2 - 1 ) / 32          , 5*$ST**2*$CT*$CF4 * ( 11*$CT**2 - 5 ) / 8 , 3*$ST*$CF3* ( 55*$CT**4 - 42*$CT**2 + 3)/ 32  , $CT*$CF2 *( 165*$CT**4 - 186*$CT**2 + 37 ) / 16 ,-5*$ST*$CF *( 33*$CT**4 - 18*$CT**2 + 1 ) / 32      ,      0                                           ,5*$ST*$SF *(33*$CT**4-18*$CT**2+1)/32     ,- $CT*$SF2 *( 165*$CT**4 - 186*$CT**2+37)/16 ,- 3*$ST*$SF3* ( 55*$CT**4 - 42*$CT**2 + 3 ) / 32 ,- 5*$ST**2*$CT*$SF4 * ( 11*$CT**2 - 5 ) / 8       ,- 5*$ST**3*$SF5 *( 5*$CT**2 - 1 ) / 32,- 15*$ST**4*$CT*$SF6 / 16],
  [ 3*$ST**5*$CF6 / 4                             , 5*$ST**4*$CT*$CF5 / 8                         , $ST**3*$CF4* ( 11*$CT**2 - 1 ) / 2        , 3*$ST**2*$CT*$CF3 * ( 11*$CT**2 - 3 ) / 8     , $ST*$CF2 *( 33*$CT**4 - 18*$CT**2 + 1 ) / 4     ,  $CT*$CF * ( 33*$CT**4 - 30*$CT**2 + 5 ) / 8       ,      0                                           ,- $CT*$SF * (33*$CT**4-30*$CT**2+5)/8     ,- $ST*$SF2 *( 33*$CT**4 - 18*$CT**2 +1)/ 4   ,- 3*$ST**2*$CT*$SF3 * ( 11*$CT**2 - 3 ) / 8      ,- $ST**3*$SF4* ( 11*$CT**2 - 1 ) / 2              ,- 5*$ST**4*$CT*$SF5 / 8,- 3*$ST**5*$SF6 / 4 ],
  [$ST**6*$SF6 / 16                               ,$ST**5*$CT*$SF5 / 16                           ,$ST**4*$SF4 * ( 11*$CT**2 - 1 ) / 16       ,$ST**3*$CT*$SF3* ( 11*$CT**2 - 3 ) / 16        ,$ST**2*$SF2 * ( 33*$CT**4 - 18*$CT**2 + 1 ) / 16 ,$ST2*$SF * ( 33*$CT**4 - 30*$CT**2 + 5 ) / 32       ,( 231*$CT**6 - 315*$CT**4 + 105*$CT**2 - 5 ) / 16 ,$ST2*$CF*(33*$CT**4-30*$CT**2+5) / 32     ,$ST**2*$CF2 * ( 33*$CT**4 - 18*$CT**2 +1)/16 ,$ST**3*$CT*$CF3* ( 11*$CT**2 - 3 ) / 16          ,$ST**4*$CF4 * ( 11*$CT**2 - 1 ) / 16              ,$ST**5*$CT*$CF5 / 16 ,$ST**6*$CF6 / 16],
  [3*$ST**5*$CT*$SF6 / 4                          ,$ST**4*$SF5 * ( 6*$CT**2 - 1 ) / 8             ,$ST**3*$CT*$SF4 * ( 33*$CT**2 - 13 ) / 4   ,3*$ST**2*$SF3 * ( 22*$CT**4 - 15*$CT**2 + 1)/8 ,$ST2*$SF2 * ( 99*$CT**4 - 102*$CT**2 + 19 ) / 8  ,$SF*( 198*$CT**6 - 285*$CT**4 + 100*$CT**2 - 5 ) / 8,- 21*$ST2* ( 33*$CT**4 - 30*$CT**2 + 5 ) / 8      ,$CF*(198*$CT**6-285*$CT**4+100*$CT**2-5)/8,$ST2*$CF2 * ( 99*$CT**4 - 102*$CT**2 +19)/8  ,3*$ST**2*$CF3 * ( 22*$CT**4 - 15*$CT**2 + 1 ) / 8,$ST**3*$CT*$CF4 * ( 33*$CT**2 - 13 ) / 4          ,$ST**4*$CF5 * ( 6*$CT**2 - 1 ) / 8,3*$ST**5*$CT*$CF6 / 4],
  [15*$ST**4*$SF6 * ( $CT**2 + 1 ) / 32           ,5*$ST**3*$CT*$SF5* ( 3*$CT**2 + 1 ) / 32       ,5*$ST**2*$SF4 * ( 33*$CT**4-10*$CT**2+1)/32,3*$ST2*$SF3 * ( 55*$CT**4 - 50*$CT**2 + 11)/64 ,$SF2*(495*$CT**6 - 735*$CT**4 + 289*$CT**2-17)/32,- 5*$ST2*$SF* ( 99*$CT**4 - 102*$CT**2 + 19 ) / 64  ,105*$ST**2 * ( 33*$CT**4 - 18*$CT**2 + 1 ) / 32   , - 5*$ST2*$CF*(99*$CT**4-102*$CT**2+19)/64,$CF2*(495*$CT**6-735*$CT**4+289*$CT**2-17)/32,3*$ST2*$CF3 * ( 55*$CT**4 - 50*$CT**2 + 11 ) / 64,5*$ST**2*$CF4 * ( 33*$CT**4 - 10*$CT**2 + 1 ) / 32,5*$ST**3*$CT*$CF5* ( 3*$CT**2 + 1 ) / 32,15*$ST**4*$CF6 * ( $CT**2 + 1 ) / 32],
  [5*$ST**3*$CT*$SF6 * ( $CT**2 + 3 ) / 8         ,5*$ST**2*$SF5 *( 2*$CT**4 + 3*$CT**2 - 1 ) / 16,5*$ST2*$SF4* ( 11*$CT**4 + 2*$CT**2 - 5)/16,$SF3*( 110*$CT**6 - 105*$CT**4+12*$CT**2-1)/16 ,- 3*$ST2*$SF2 * ( 55*$CT**4 - 50*$CT**2 + 11 )/16,15*$ST**2*$SF * ( 22*$CT**4 - 15*$CT**2 + 1 ) / 16  ,- 105*$ST**3*$CT *( 11*$CT**2 - 3 ) / 8           ,15*$ST**2*$CF*(22*$CT**4-15*$CT**2+1) / 16, - 3*$ST2*$CF2 * ( 55*$CT**4-50*$CT**2+11)/16,$CF3*( 110*$CT**6 - 105*$CT**4 + 12*$CT**2 -1)/16,5*$ST2*$CF4*(11*$CT**4 + 2*$CT**2 - 5 ) / 16      ,5*$ST**2*$CF5 *( 2*$CT**4 + 3*$CT**2 - 1 ) / 16,5*$ST**3*$CT*$CF6 * ( $CT**2 + 3 ) / 8 ],
  [3*$ST**2*$SF6 *( $CT**4 + 6*$CT**2 + 1 ) / 16  ,$ST2*$SF5 * ( 3*$CT**4 + 10*$CT**2 - 5 ) / 32  ,$SF4*( 33*$CT**6+35*$CT**4-65*$CT**2+13)/16,- 3*$ST2*$SF3 *( 11*$CT**4 + 2*$CT**2 - 5 ) /32,3*$ST**2*$SF2* ( 33*$CT**4 - 10*$CT**2 + 1 ) / 16,- 3*$ST**3*$CT*$SF * ( 33*$CT**2 - 13 ) / 16        ,63*$ST**4 * ( 11*$CT**2 - 1 ) / 16                , - 3*$ST**3*$CT*$CF*(33*$CT**2-13 ) / 16  ,3*$ST**2*$CF2* ( 33*$CT**4 - 10*$CT**2+1)/ 16,- 3*$ST2*$CF3 *( 11*$CT**4 + 2*$CT**2 - 5 ) / 32 ,$CF4*(33*$CT**6+35*$CT**4-65*$CT**2 + 13 ) / 16   ,$ST2*$CF5 * ( 3*$CT**4 + 10*$CT**2 - 5 ) / 32 ,3*$ST**2*$CF6 *( $CT**4 + 6*$CT**2 + 1 ) / 16],
  [3*$ST2*$SF6 * ( $CT**4 + 10*$CT**2 + 5 ) / 16  ,$SF5 * ( 6*$CT**6 + 35*$CT**4 - 20*$CT**2-5)/16,- 11*$ST2*$SF4 * ( 3*$CT**4+10*$CT**2-5)/16,33*$ST**2*$SF3 * ( 2*$CT**4 + 3*$CT**2 -1)/16  ,- 33*$ST**3*$CT*$SF2 *( 3*$CT**2 + 1 ) / 8       ,33*$ST**4*$SF * ( 6*$CT**2 - 1 ) / 16               , - 693*$ST**5*$CT / 8                             ,33*$ST**4*$CF * ( 6*$CT**2 - 1 ) / 16     , - 33*$ST**3*$CT*$CF2 *( 3*$CT**2 + 1 ) / 8  ,33*$ST**2*$CF3 * ( 2*$CT**4 + 3*$CT**2 - 1 ) / 16,- 11*$ST2*$CF4 * ( 3*$CT**4 + 10*$CT**2 - 5 ) / 16,$CF5 * ( 6*$CT**6 + 35*$CT**4 - 20*$CT**2 - 5 ) / 16,3*$ST2*$CF6 * ( $CT**4 + 10*$CT**2 + 5 ) / 16],
  [$SF6 *( $CT**6 + 15*$CT**4 + 15*$CT**2 + 1)/32 ,- $ST2*$SF5 * ( $CT**4 + 10*$CT**2 + 5 ) / 64  ,11*$ST**2*$SF4 *( $CT**4 + 6*$CT**2+1) / 32,- 11*$ST**3*$CT*$SF3* ( $CT**2 + 3 ) / 32      ,33*$ST**4*$SF2 *( $CT**2 + 1 ) / 32              ,- 33*$ST**5*$CT*$SF / 32                            ,231*$ST**6 / 32                                   ,- 33*$ST**5*$CT*$CF / 32                  ,33*$ST**4*$CF2 *( $CT**2 + 1 ) / 32          ,- 11*$ST**3*$CT*$CF3* ( $CT**2 + 3 ) / 32        ,11*$ST**2*$CF4 *( $CT**4 + 6*$CT**2 + 1 ) / 32    ,- $ST2*$CF5 * ( $CT**4 + 10*$CT**2 + 5 ) / 64,$CF6 *( $CT**6 + 15*$CT**4 + 15*$CT**2 + 1 ) / 32]
];

my $B2 = [ [$B{"B22S"}],[-$B{"B21S"}],[$B{"B20"}],[-$B{"B21"}],[$B{"B22"}] ];
my $R2 = mmult($MR2, $B2);
my $B4 = [ [$B{"B44S"}],[-$B{"B43S"}],[$B{"B42S"}],[-$B{"B41S"}],[$B{"B40"}],[-$B{"B41"}],[$B{"B42"}],[-$B{"B43"}],[$B{"B44"}] ];
my $R4 = mmult($MR4, $B4);
my $B6 = [ [$B{"B66S"}],[-$B{"B65S"}],[$B{"B64S"}],[-$B{"B63S"}],[$B{"B62S"}],[-$B{"B61S"}],[$B{"B60"}],[-$B{"B61"}],[$B{"B62"}],[-$B{"B63"}],[$B{"B64"}],[-$B{"B65"}],[$B{"B66"}] ];
my $R6 = mmult($MR6, $B6);
my %Bo2 = (0, "B22S", 1, "B21S", 2, "B20 ", 3, "B21 ", 4, "B22 ");
my %Bo4 = (0, "B44S", 1, "B43S", 2, "B42S", 3, "B41S", 4, "B40 ", 5, "B41 ", 6, "B42 ", 7, "B43 ", 8, "B44 ");
my %Bo6 = (0, "B66S", 1, "B65S", 2, "B64S", 3, "B63S", 4, "B62S", 5, "B61S", 6, "B60 ", 7, "B61 ", 8, "B62 ", 9, "B63 ", 10, "B64 ", 11, "B65 ", 12, "B66 ");

if (!$output) {
  if ($verbose) {
    print "Crystal field parameters input:\n";
    while ( ($key, $value) = each %B) {
      print "$key = $value\n";
    }
    if ($Ion) { print "\nIontype is ".$Ion."\n"; }
    print "\nRotated crystal field parameters:\n";
  }
  for $i (0 .. 5) { if ($_ = $R2->[$i][0]) { if(abs($_)>1e-8) {if($i % 2){$_=-$_}; print "$Bo2{$i} = $_\n"; } } }
  for $i (0 .. 9) { if ($_ = $R4->[$i][0]) { if(abs($_)>1e-8) {if($i % 2){$_=-$_}; print "$Bo4{$i} = $_\n"; } } }
  for $i (0 .. 13) { if ($_ = $R6->[$i][0]) { if(abs($_)>1e-8){if($i % 2){$_=-$_}; print "$Bo6{$i} = $_\n"; } } }
}
else {
  open (outfile, ">$output") or die "$0: cannot open $output for output.";
 
if ($input) {
  open (input_file, $input) or die "$0: cannot open $input for input CF parameters";
  while(<input_file>) {                                   # Selects out lines with crystal field parameters
   unless ($_ =~ s/(B[0-9\ CcSs]+)\s*[=:]\s*([-\.\de]*)// ) { # () are groups which may be access with $1, $2 etc.
    print outfile $_;                                         # * means match previous char any number of times.       
    }}
  close (input_file);}
else 
 { print outfile << "EOF";
#!MODULE=so1ion
#<!--mcphase.sipf-->
EOF
  
  $Ion =~ s/ //;
  print outfile " IONTYPE=$Ion\n";
  print outfile << "EOF";
#   D = 2 * pi / Q
#   s = 1 / 2 / D: sintheta = lambda * s
# Debey-Waller Factor
DWF=0
# debeywallerfactor = EXP(-2 * DWF *s*s)
EOF
}
 for $i (0 .. 5) { if ($_ = $R2->[$i][0]) { if(abs($_)>1e-8) {if($i % 2){$_=-$_}; print outfile "$Bo2{$i} = $_\n"; } } }
  for $i (0 .. 9) { if ($_ = $R4->[$i][0]) { if(abs($_)>1e-8) {if($i % 2){$_=-$_}; print  outfile "$Bo4{$i} = $_\n"; } } }
  for $i (0 .. 13) { if ($_ = $R6->[$i][0]) { if(abs($_)>1e-8){if($i % 2){$_=-$_}; print  outfile "$Bo6{$i} = $_\n"; } } }
  close (outfile);
}

