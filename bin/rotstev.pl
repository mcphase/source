#!/usr/bin/perl
# rotstev.pl 
#
# Generates rotation matrices to rotate Stevens' Operator equivalents by
# phi about y-axis and theta about z-axis. Based on the method of Buckmaster
# (Physica Status Solidi A vol 13, pp 9, 1972) and Rudowicz (J. Phys: C: 
# Solid State Physics, vol 18, pp 1415, 1985).
#
# NB: phi is denoted by f and theta by t in this program
#
# by Duc Le 2006-2010 - duc.le@helmholtz-berlin.de
#
# Initial version 2006
# Updated Thu Jan 20 18:33:07 WEST 2011 - mdl - Added k=3,5 (octupolar/5-pol)

use Math::Complex;

# The following subroutines were shamelessly stolen from the Perl Cookbook
# by O'Reily, and used for the matrix multiplications.

sub mmult {
  my ($m1,$m2) = @_;
  my ($m1rows,$m1cols) = matdim($m1);
  my ($m2rows,$m2cols) = matdim($m2);

  unless ($m1cols == $m2rows) {       # raise exception
    die "IndexError: matrices don't match: $m1cols != $m2rows";
  }

  my $result = [];
  my ($i, $j, $k, $m1el, $m2el, $prod_el);

  for $i (range($m1rows)) {
    for $j (range($m2cols)) {
      $prod_el = 0;
      for $k (range($m1cols)) {
        $m1el = $m1->[$i][$k];
        $m2el = $m2->[$k][$j];
        $prod_el += $m1el * $m2el;
      }
      $result->[$i][$j] = $prod_el;
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

sub zeros {
    my ($rows,$cols) = @_;
    my $ii,$jj;
    my @result;
    for $ii(0..($rows-1)) { my @blank; for $jj(0..($cols-1)) { push(@blank,0); } push(@result,\@blank); }
    return \@result;
}

sub eye {
    my $result = zeros($_[0],$_[0]);
    my $ii; for $ii(0..(@_[0]-1)) { $result->[$ii][$ii]=1; }
    return $result;
}

# Routines to calculate the rotation matrices
# --------------------------------------- O1 ------------------------------------- #
sub rot1 { my ($t,$f) = @_; my $D;

my $D10p1 = (sqrt(2)/2)*sin($t);
my $D10p0 = cos($t);
my $D1p1m1 = cplxe(1, $f) *  (1/2)*(1-cos($t));
my $D1p1p1 = cplxe(1, $f) *  (1/2)*(1+cos($t));
my $D1p1p0 = cplxe(1, $f) * -(sqrt(2)/2)*sin($t);
my $D1m1p0 = cplxe(1,-$f) *  (sqrt(2)/2)*sin($t);
my $D1m1m1 = cplxe(1,-$f) *  (1/2)*(1+cos($t));
my $D1m1p1 = cplxe(1,-$f) *  (1/2)*(1-cos($t));

$D = [ 
           [$D1m1m1, $D1m1p0, $D1m1p1],
           [-$D10p1,  $D10p0,  $D10p1],
           [$D1p1m1, $D1p1p0, $D1p1p1],
         ];
return $D;
}
# --------------------------------------- O2 ------------------------------------- #
sub rot2 { my ($t,$f) = @_; my $D;

my $D20p2 = (sqrt(6)/4)*sin($t)**2;
my $D20p1 = (sqrt(6)/2)*sin($t)*cos($t);
my $D20p0 = (1/2)*(3*cos($t)**2 - 1);

my $D2m1m2 = cplxe(1,-$f) * -(1/2)*sin($t)*(1+cos($t));
my $D2m1p2 = cplxe(1,-$f) *  (1/2)*sin($t)*(1-cos($t));
my $D2m1m1 = cplxe(1,-$f) *  (1/2)*(1+cos($t))*(2*cos($t)-1);
my $D2m1p1 = cplxe(1,-$f) *  (1/2)*(1-cos($t))*(2*cos($t)+1);
my $D2m1p0 = cplxe(1,-$f) *  (sqrt(6)/2)*sin($t)*cos($t);

my $D2p1m2 = cplxe(1, $f) * -(1/2)*sin($t)*(1-cos($t));
my $D2p1p2 = cplxe(1, $f) *  (1/2)*sin($t)*(1+cos($t));
my $D2p1m1 = cplxe(1, $f) *  (1/2)*(1-cos($t))*(2*cos($t)+1);
my $D2p1p1 = cplxe(1, $f) *  (1/2)*(1+cos($t))*(2*cos($t)-1);
my $D2p1p0 = cplxe(1, $f) * -(sqrt(6)/2)*sin($t)*cos($t);

my $D2m2m2 = cplxe(1,-2*$f) *  (1/4)*(1+cos($t))**2;
my $D2m2p2 = cplxe(1,-2*$f) *  (1/4)*(1-cos($t))**2;
my $D2m2m1 = cplxe(1,-2*$f) *  (1/2)*sin($t)*(1+cos($t));
my $D2m2p1 = cplxe(1,-2*$f) *  (1/2)*sin($t)*(1-cos($t));
my $D2m2p0 = cplxe(1,-2*$f) *  (sqrt(6)/4)*sin($t)**2;

my $D2p2m2 = cplxe(1, 2*$f) *  (1/4)*(1-cos($t))**2;
my $D2p2p2 = cplxe(1, 2*$f) *  (1/4)*(1+cos($t))**2;
my $D2p2m1 = cplxe(1, 2*$f) * -(1/2)*sin($t)*(1-cos($t));
my $D2p2p1 = cplxe(1, 2*$f) * -(1/2)*sin($t)*(1+cos($t));
my $D2p2p0 = cplxe(1, 2*$f) *  (sqrt(6)/4)*sin($t)**2;

$D = [
	   [$D2m2m2, $D2m2m1, $D2m2p0, $D2m2p1, $D2m2p2],
           [$D2m1m2, $D2m1m1, $D2m1p0, $D2m1p1, $D2m1p2],
           [$D20p2, -$D20p1,  $D20p0,  $D20p1,  $D20p2 ],
           [$D2p1m2, $D2p1m1, $D2p1p0, $D2p1p1, $D2p1p2],
           [$D2p2m2, $D2p2m1, $D2p2p0, $D2p2p1, $D2p2p2]
         ];
return $D;
}

# --------------------------------------- O3 ------------------------------------- #
sub rot3 { my ($t,$f) = @_; my $D;

my $D30p0 = (1/2)*cos($t)*(5*cos($t)**2-3);
my $D30p1 = (sqrt(3)/4)*sin($t)*(5*cos($t)**2-1);
my $D30p2 = (sqrt(30)/4)*(sin($t)**2)*cos($t);
my $D30p3 = (sqrt(5)/4)*(sin($t)**3);

my $D3m1p0 = cplxe(1,-$f) *  (sqrt(3)/4)*sin($t)*(5*cos($t)**2-1);
my $D3m1p1 = cplxe(1,-$f) *  (1/8)*(1-cos($t))*(15*cos($t)**2+10*cos($t)-1);
my $D3m1m1 = cplxe(1,-$f) *  (1/8)*(1+cos($t))*(15*cos($t)**2-10*cos($t)-1);
my $D3m1p2 = cplxe(1,-$f) *  (sqrt(10)/8)*sin($t)*(1-cos($t))*(3*cos($t)+1);
my $D3m1m2 = cplxe(1,-$f) * -(sqrt(10)/8)*sin($t)*(1+cos($t))*(3*cos($t)-1);
my $D3m1p3 = cplxe(1,-$f) *  (sqrt(15)/8)*(sin($t)**2)*(1-cos($t));
my $D3m1m3 = cplxe(1,-$f) *  (sqrt(15)/8)*(sin($t)**2)*(1+cos($t));

my $D3p1p0 = cplxe(1, $f) * -(sqrt(3)/4)*sin($t)*(5*cos($t)**2-1);
my $D3p1p1 = cplxe(1, $f) *  (1/8)*(1+cos($t))*(15*cos($t)**2-10*cos($t)-1);
my $D3p1m1 = cplxe(1, $f) *  (1/8)*(1-cos($t))*(15*cos($t)**2+10*cos($t)-1);
my $D3p1p2 = cplxe(1, $f) *  (sqrt(10)/8)*sin($t)*(1+cos($t))*(3*cos($t)-1);
my $D3p1m2 = cplxe(1, $f) * -(sqrt(10)/8)*sin($t)*(1-cos($t))*(3*cos($t)+1);
my $D3p1p3 = cplxe(1, $f) *  (sqrt(15)/8)*(sin($t)**2)*(1+cos($t));
my $D3p1m3 = cplxe(1, $f) *  (sqrt(15)/8)*(sin($t)**2)*(1-cos($t));

my $D3m2p0 = cplxe(1,-2*$f) *  (sqrt(30)/4)*(sin($t)**2)*cos($t);
my $D3m2p1 = cplxe(1,-2*$f) *  (sqrt(10)/8)*sin($t)*(1-cos($t))*(3*cos($t)+1);
my $D3m2m1 = cplxe(1,-2*$f) *  (sqrt(10)/8)*sin($t)*(1+cos($t))*(3*cos($t)-1);
my $D3m2p2 = cplxe(1,-2*$f) *  (1/4)*((1-cos($t))**2)*(3*cos($t)+2);
my $D3m2m2 = cplxe(1,-2*$f) *  (1/4)*((1+cos($t))**2)*(3*cos($t)-2);
my $D3m2p3 = cplxe(1,-2*$f) *  (sqrt(6)/8)*sin($t)*(1-cos($t))**2;
my $D3m2m3 = cplxe(1,-2*$f) * -(sqrt(6)/8)*sin($t)*(1+cos($t))**2;

my $D3p2p0 = cplxe(1, 2*$f) *  (sqrt(30)/4)*(sin($t)**2)*cos($t);
my $D3p2p1 = cplxe(1, 2*$f) * -(sqrt(10)/8)*sin($t)*(1+cos($t))*(3*cos($t)-1);
my $D3p2m1 = cplxe(1, 2*$f) * -(sqrt(10)/8)*sin($t)*(1-cos($t))*(3*cos($t)+1);
my $D3p2p2 = cplxe(1, 2*$f) *  (1/4)*((1+cos($t))**2)*(3*cos($t)-2);
my $D3p2m2 = cplxe(1, 2*$f) *  (1/4)*((1-cos($t))**2)*(3*cos($t)+2);
my $D3p2p3 = cplxe(1, 2*$f) *  (sqrt(6)/8)*sin($t)*(1+cos($t))**2;
my $D3p2m3 = cplxe(1, 2*$f) * -(sqrt(6)/8)*sin($t)*(1-cos($t))**2;

my $D3m3p0 = cplxe(1,-3*$f) *  (sqrt(5)/4)*(sin($t)**3);
my $D3m3p1 = cplxe(1,-3*$f) *  (sqrt(15)/8)*(sin($t)**2)*(1-cos($t));
my $D3m3m1 = cplxe(1,-3*$f) *  (sqrt(15)/8)*(sin($t)**2)*(1+cos($t));
my $D3m3p2 = cplxe(1,-3*$f) *  (sqrt(6)/8)*sin($t)*(1-cos($t))**2;
my $D3m3m2 = cplxe(1,-3*$f) *  (sqrt(6)/8)*sin($t)*(1+cos($t))**2;
my $D3m3p3 = cplxe(1,-3*$f) *  (1/8)*(1-cos($t))**3;
my $D3m3m3 = cplxe(1,-3*$f) *  (1/8)*(1+cos($t))**3;

my $D3p3p0 = cplxe(1, 3*$f) * -(sqrt(5)/4)*(sin($t)**3);
my $D3p3p1 = cplxe(1, 3*$f) *  (sqrt(15)/8)*(sin($t)**2)*(1+cos($t));
my $D3p3m1 = cplxe(1, 3*$f) *  (sqrt(15)/8)*(sin($t)**2)*(1-cos($t));
my $D3p3p2 = cplxe(1, 3*$f) * -(sqrt(6)/8)*sin($t)*(1+cos($t))**2;
my $D3p3m2 = cplxe(1, 3*$f) * -(sqrt(6)/8)*sin($t)*(1-cos($t))**2;
my $D3p3p3 = cplxe(1, 3*$f) *  (1/8)*(1+cos($t))**3;
my $D3p3m3 = cplxe(1, 3*$f) *  (1/8)*(1-cos($t))**3;

$D = [
	   [$D3m3m3, $D3m3m2, $D3m3m1, $D3m3p0, $D3m3p1, $D3m3p2, $D3m3p3],
	   [$D3m2m3, $D3m2m2, $D3m2m1, $D3m2p0, $D3m2p1, $D3m2p2, $D3m2p3],
           [$D3m1m3, $D3m1m2, $D3m1m1, $D3m1p0, $D3m1p1, $D3m1p2, $D3m1p3],
           [-$D30p3, $D30p2, -$D30p1,  $D30p0,  $D30p1,  $D30p2,  $D30p3],
           [$D3p1m3, $D3p1m2, $D3p1m1, $D3p1p0, $D3p1p1, $D3p1p2, $D3p1p3],
           [$D3p2m3, $D3p2m2, $D3p2m1, $D3p2p0, $D3p2p1, $D3p2p2, $D3p2p3],
	   [$D3p3m3, $D3p3m2, $D3p3m1, $D3p3p0, $D3p3p1, $D3p3p2, $D3p3p3],
         ];
return $D;
}

# --------------------------------------- O4 ------------------------------------- #
sub rot4 { my ($t,$f) = @_; my $D;

my $D4p0p4 =  sqrt(70)/16  * sin($t)**4;
my $D4p0p3 =  sqrt(35)/4   * sin($t)**3 * cos($t);
my $D4p0p2 =  sqrt(10)/8   * sin($t)**2 * (7*cos($t)**2 - 1);
my $D4p0p1 =  sqrt(5)/4    * sin($t)    * cos($t) * (7*cos($t)**2 - 3);
my $D4p0p0 =  (1/8)        * (35*cos($t)**4 - 30*cos($t)**2 + 3);

my $D4p1p0 = -cplxe(1, $f) * sqrt(5)/4     * sin($t)     * cos($t)     * (7*cos($t)**2-3);
my $D4p1p1 =  cplxe(1, $f) * (1/8)         * (1+cos($t)) * (28*cos($t)**3 - 21*cos($t)**2 - 6*cos($t) + 3);
my $D4p1m1 =  cplxe(1, $f) * (1/8)         * (1-cos($t)) * (28*cos($t)**3 + 21*cos($t)**2 - 6*cos($t) - 3);
my $D4p1p2 =  cplxe(1, $f) * sqrt(2)/8     * sin($t)     * (1+cos($t)) * (14*cos($t)**2 - 7*cos($t) - 1);
my $D4p1m2 = -cplxe(1, $f) * sqrt(2)/8     * sin($t)     * (1-cos($t)) * (14*cos($t)**2 + 7*cos($t) - 1);
my $D4p1p3 =  cplxe(1, $f) * sqrt(7)/8     * sin($t)**2  * (1+cos($t)) * ( 4*cos($t) - 1);
my $D4p1m3 =  cplxe(1, $f) * sqrt(7)/8     * sin($t)**2  * (1-cos($t)) * ( 4*cos($t) + 1);
my $D4p1p4 =  cplxe(1, $f) * sqrt(14)/8    * sin($t)**3  * (1+cos($t));
my $D4p1m4 = -cplxe(1, $f) * sqrt(14)/8    * sin($t)**3  * (1-cos($t));

my $D4m1p0 =  cplxe(1,-$f) * sqrt(5)/4     * sin($t)     * cos($t)     * (7*cos($t)**2-3);
my $D4m1p1 =  cplxe(1,-$f) * (1/8)         * (1-cos($t)) * (28*cos($t)**3 + 21*cos($t)**2 - 6*cos($t) - 3);
my $D4m1m1 =  cplxe(1,-$f) * (1/8)         * (1+cos($t)) * (28*cos($t)**3 - 21*cos($t)**2 - 6*cos($t) + 3);
my $D4m1p2 =  cplxe(1,-$f) * sqrt(2)/8     * sin($t)     * (1-cos($t)) * (14*cos($t)**2 + 7*cos($t) - 1);
my $D4m1m2 = -cplxe(1,-$f) * sqrt(2)/8     * sin($t)     * (1+cos($t)) * (14*cos($t)**2 - 7*cos($t) - 1);
my $D4m1p3 =  cplxe(1,-$f) * sqrt(7)/8     * sin($t)**2  * (1-cos($t)) * ( 4*cos($t) + 1);
my $D4m1m3 =  cplxe(1,-$f) * sqrt(7)/8     * sin($t)**2  * (1+cos($t)) * ( 4*cos($t) - 1);
my $D4m1p4 =  cplxe(1,-$f) * sqrt(14)/8    * sin($t)**3  * (1-cos($t));
my $D4m1m4 = -cplxe(1,-$f) * sqrt(14)/8    * sin($t)**3  * (1+cos($t));

my $D4p2p0 =  cplxe(1, 2*$f) * sqrt(10)/8  * sin($t)**2 * (7*cos($t)**2 - 1);
my $D4p2p1 = -cplxe(1, 2*$f) * sqrt(2)/8   * sin($t)    * (1+cos($t)) * (14*cos($t)**2 - 7*cos($t) - 1);
my $D4p2m1 = -cplxe(1, 2*$f) * sqrt(2)/8   * sin($t)    * (1-cos($t)) * (14*cos($t)**2 + 7*cos($t) - 1);
my $D4p2p2 =  cplxe(1, 2*$f) * (1/4)       * (1+cos($t))**2 * (7*cos($t)**2 - 7*cos($t) + 1);
my $D4p2m2 =  cplxe(1, 2*$f) * (1/4)       * (1-cos($t))**2 * (7*cos($t)**2 + 7*cos($t) + 1);
my $D4p2p3 =  cplxe(1, 2*$f) * sqrt(14)/8  * sin($t)    * (1+cos($t))**2 * (2*cos($t) - 1);
my $D4p2m3 = -cplxe(1, 2*$f) * sqrt(14)/8  * sin($t)    * (1-cos($t))**2 * (2*cos($t) + 1);
my $D4p2p4 =  cplxe(1, 2*$f) * sqrt(7)/8   * sin($t)**2 * (1+cos($t))**2;
my $D4p2m4 =  cplxe(1, 2*$f) * sqrt(7)/8   * sin($t)**2 * (1-cos($t))**2;

my $D4m2p0 =  cplxe(1,-2*$f) * sqrt(10)/8  * sin($t)**2 * (7*cos($t)**2 - 1);
my $D4m2p1 =  cplxe(1,-2*$f) * sqrt(2)/8   * sin($t)    * (1-cos($t)) * (14*cos($t)**2 + 7*cos($t) - 1);
my $D4m2m1 =  cplxe(1,-2*$f) * sqrt(2)/8   * sin($t)    * (1+cos($t)) * (14*cos($t)**2 - 7*cos($t) - 1);
my $D4m2p2 =  cplxe(1,-2*$f) * (1/4)       * (1-cos($t))**2 * (7*cos($t)**2 + 7*cos($t) + 1);
my $D4m2m2 =  cplxe(1,-2*$f) * (1/4)       * (1+cos($t))**2 * (7*cos($t)**2 - 7*cos($t) + 1);
my $D4m2p3 =  cplxe(1,-2*$f) * sqrt(14)/8  * sin($t)    * (1-cos($t))**2 * (2*cos($t) + 1);
my $D4m2m3 = -cplxe(1,-2*$f) * sqrt(14)/8  * sin($t)    * (1+cos($t))**2 * (2*cos($t) - 1);
my $D4m2p4 =  cplxe(1,-2*$f) * sqrt(7)/8   * sin($t)**2 * (1-cos($t))**2;
my $D4m2m4 =  cplxe(1,-2*$f) * sqrt(7)/8   * sin($t)**2 * (1+cos($t))**2;

my $D4p3p0 = -cplxe(1, 3*$f) * sqrt(35)/4  * sin($t)**3 * cos($t);
my $D4p3p1 =  cplxe(1, 3*$f) * sqrt(7)/8   * sin($t)**2 * (1+cos($t)) * (4*cos($t) - 1);
my $D4p3m1 =  cplxe(1, 3*$f) * sqrt(7)/8   * sin($t)**2 * (1-cos($t)) * (4*cos($t) + 1);
my $D4p3p2 = -cplxe(1, 3*$f) * sqrt(14)/8  * sin($t)    * (1+cos($t))**2 * (2*cos($t) - 1);
my $D4p3m2 = -cplxe(1, 3*$f) * sqrt(14)/8  * sin($t)    * (1-cos($t))**2 * (2*cos($t) + 1);
my $D4p3p3 =  cplxe(1, 3*$f) * (1/8)       * (1+cos($t))**3 * (4*cos($t) - 3);
my $D4p3m3 =  cplxe(1, 3*$f) * (1/8)       * (1-cos($t))**3 * (4*cos($t) + 3);
my $D4p3p4 =  cplxe(1, 3*$f) * sqrt(2)/8   * sin($t)    * (1+cos($t))**3;
my $D4p3m4 = -cplxe(1, 3*$f) * sqrt(2)/8   * sin($t)    * (1-cos($t))**3;

my $D4m3p0 =  cplxe(1,-3*$f) * sqrt(35)/4  * sin($t)**3 * cos($t);
my $D4m3p1 =  cplxe(1,-3*$f) * sqrt(7)/8   * sin($t)**2 * (1-cos($t)) * (4*cos($t) + 1);
my $D4m3m1 =  cplxe(1,-3*$f) * sqrt(7)/8   * sin($t)**2 * (1+cos($t)) * (4*cos($t) - 1);
my $D4m3p2 =  cplxe(1,-3*$f) * sqrt(14)/8  * sin($t)    * (1-cos($t))**2 * (2*cos($t) + 1);
my $D4m3m2 =  cplxe(1,-3*$f) * sqrt(14)/8  * sin($t)    * (1+cos($t))**2 * (2*cos($t) - 1);
my $D4m3p3 =  cplxe(1,-3*$f) * (1/8)       * (1-cos($t))**3 * (4*cos($t) + 3);
my $D4m3m3 =  cplxe(1,-3*$f) * (1/8)       * (1+cos($t))**3 * (4*cos($t) - 3);
my $D4m3p4 =  cplxe(1,-3*$f) * sqrt(2)/8   * sin($t)    * (1-cos($t))**3;
my $D4m3m4 = -cplxe(1,-3*$f) * sqrt(2)/8   * sin($t)    * (1+cos($t))**3;

my $D4p4p0 =  cplxe(1, 4*$f) * sqrt(70)/16 * sin($t)**4;
my $D4p4p1 = -cplxe(1, 4*$f) * sqrt(14)/8  * sin($t)**3 * (1+cos($t));
my $D4p4m1 = -cplxe(1, 4*$f) * sqrt(14)/8  * sin($t)**3 * (1-cos($t));
my $D4p4p2 =  cplxe(1, 4*$f) * sqrt(7)/8   * sin($t)**2 * (1+cos($t))**2;
my $D4p4m2 =  cplxe(1, 4*$f) * sqrt(7)/8   * sin($t)**2 * (1-cos($t))**2;
my $D4p4p3 = -cplxe(1, 4*$f) * sqrt(2)/8   * sin($t)    * (1+cos($t))**3;
my $D4p4m3 = -cplxe(1, 4*$f) * sqrt(2)/8   * sin($t)    * (1-cos($t))**3;
my $D4p4p4 =  cplxe(1, 4*$f) * (1/16)      * (1+cos($t))**4;
my $D4p4m4 =  cplxe(1, 4*$f) * (1/16)      * (1-cos($t))**4;

my $D4m4p0 =  cplxe(1,-4*$f) * sqrt(70)/16 * sin($t)**4;
my $D4m4p1 =  cplxe(1,-4*$f) * sqrt(14)/8  * sin($t)**3 * (1-cos($t));
my $D4m4m1 =  cplxe(1,-4*$f) * sqrt(14)/8  * sin($t)**3 * (1+cos($t));
my $D4m4p2 =  cplxe(1,-4*$f) * sqrt(7)/8   * sin($t)**2 * (1-cos($t))**2;
my $D4m4m2 =  cplxe(1,-4*$f) * sqrt(7)/8   * sin($t)**2 * (1+cos($t))**2;
my $D4m4p3 =  cplxe(1,-4*$f) * sqrt(2)/8   * sin($t)    * (1-cos($t))**3;
my $D4m4m3 =  cplxe(1,-4*$f) * sqrt(2)/8   * sin($t)    * (1+cos($t))**3;
my $D4m4p4 =  cplxe(1,-4*$f) * (1/16)      * (1-cos($t))**4;
my $D4m4m4 =  cplxe(1,-4*$f) * (1/16)      * (1+cos($t))**4;

$D = [
             [ $D4m4m4,  $D4m4m3,  $D4m4m2,  $D4m4m1,  $D4m4p0,  $D4m4p1,  $D4m4p2,  $D4m4p3,  $D4m4p4], 
             [ $D4m3m4,  $D4m3m3,  $D4m3m2,  $D4m3m1,  $D4m3p0,  $D4m3p1,  $D4m3p2,  $D4m3p3,  $D4m3p4], 
             [ $D4m2m4,  $D4m2m3,  $D4m2m2,  $D4m2m1,  $D4m2p0,  $D4m2p1,  $D4m2p2,  $D4m2p3,  $D4m2p4], 
             [ $D4m1m4,  $D4m1m3,  $D4m1m2,  $D4m1m1,  $D4m1p0,  $D4m1p1,  $D4m1p2,  $D4m1p3,  $D4m1p4], 
             [ $D4p0p4, -$D4p0p3,  $D4p0p2, -$D4p0p1,  $D4p0p0,  $D4p0p1,  $D4p0p2,  $D4p0p3,  $D4p0p4], 
             [ $D4p1m4,  $D4p1m3,  $D4p1m2,  $D4p1m1,  $D4p1p0,  $D4p1p1,  $D4p1p2,  $D4p1p3,  $D4p1p4], 
             [ $D4p2m4,  $D4p2m3,  $D4p2m2,  $D4p2m1,  $D4p2p0,  $D4p2p1,  $D4p2p2,  $D4p2p3,  $D4p2p4], 
             [ $D4p3m4,  $D4p3m3,  $D4p3m2,  $D4p3m1,  $D4p3p0,  $D4p3p1,  $D4p3p2,  $D4p3p3,  $D4p3p4], 
             [ $D4p4m4,  $D4p4m3,  $D4p4m2,  $D4p4m1,  $D4p4p0,  $D4p4p1,  $D4p4p2,  $D4p4p3,  $D4p4p4]
         ];
return $D;
}

# --------------------------------------- O5 ------------------------------------- #
sub rot5 { my ($t,$f) = @_; my $D;

my $D50p0  = (1/8)*cos($t)*(63*cos($t)**4-70*cos($t)**2+15);
my $D50p1  = (sqrt(30)/16)*sin($t)*(21*cos($t)**4-14*cos($t)**2+1);
my $D50p2  = (sqrt(210)/8)*(sin($t)**2)*cos($t)*(3*cos($t)**2-1);
my $D50p3  = (sqrt(35)/16)*(sin($t)**3)*(9*cos($t)**2-1);
my $D50p4  = (3*sqrt(70)/16)*(sin($t)**4)*cos($t);
my $D50p5  = (3*sqrt(7)/16)*(sin($t)**5);

my $D5p1p0 = cplxe(1,   $f) * -(sqrt(30)/16)*sin($t)*(21*cos(4)**4-14*cos($t)**2+1);
my $D5m1p0 = cplxe(1,  -$f) *  (sqrt(30)/16)*sin($t)*(21*cos(4)**4-14*cos($t)**2+1);
my $D5p1p1 = cplxe(1,   $f) *  (1/16)*(1+cos($t))*(105*cos($t)**4-84*cos($t)**3-42*cos($t)**2+28*cos($t)+1);
my $D5p1m1 = cplxe(1,   $f) *  (1/16)*(1-cos($t))*(105*cos($t)**4+84*cos($t)**3-42*cos($t)**2-28*cos($t)+1);
my $D5m1p1 = cplxe(1,  -$f) *  (1/16)*(1-cos($t))*(105*cos($t)**4+84*cos($t)**3-42*cos($t)**2-28*cos($t)+1);
my $D5m1m1 = cplxe(1,  -$f) *  (1/16)*(1+cos($t))*(105*cos($t)**4-84*cos($t)**3-42*cos($t)**2+28*cos($t)+1);
my $D5p1p2 = cplxe(1,   $f) *  (sqrt(7)/8)*sin($t)*(1+cos($t))*(15*cos($t)**3-9*cos($t)**2-3*cos($t)+1);
my $D5p1m2 = cplxe(1,   $f) *  (sqrt(7)/8)*sin($t)*(1-cos($t))*(15*cos($t)**3+9*cos($t)**2-3*cos($t)-1);
my $D5m1p2 = cplxe(1,  -$f) *  (sqrt(7)/8)*sin($t)*(1-cos($t))*(15*cos($t)**3+9*cos($t)**2-3*cos($t)-1);
my $D5m1m2 = cplxe(1,  -$f) *  (sqrt(7)/8)*sin($t)*(1+cos($t))*(15*cos($t)**3-9*cos($t)**2-3*cos($t)+1);
my $D5p1p3 = cplxe(1,   $f) *  (sqrt(42)/32)*(sin($t)**2)*(1+cos($t))*(15*cos($t)**2-6*cos($t)-1);
my $D5p1m3 = cplxe(1,   $f) *  (sqrt(42)/32)*(sin($t)**2)*(1-cos($t))*(15*cos($t)**2+6*cos($t)-1);
my $D5m1p3 = cplxe(1,  -$f) *  (sqrt(42)/32)*(sin($t)**2)*(1-cos($t))*(15*cos($t)**2+6*cos($t)-1);
my $D5m1m3 = cplxe(1,  -$f) *  (sqrt(42)/32)*(sin($t)**2)*(1+cos($t))*(15*cos($t)**2-6*cos($t)-1);
my $D5p1p4 = cplxe(1,   $f) *  (sqrt(21)/16)*(sin($t)**3)*(1+cos($t))*(5*cos($t)-1);
my $D5p1m4 = cplxe(1,   $f) *  (sqrt(21)/16)*(sin($t)**3)*(1-cos($t))*(5*cos($t)+1);
my $D5m1p4 = cplxe(1,  -$f) *  (sqrt(21)/16)*(sin($t)**3)*(1-cos($t))*(5*cos($t)+1);
my $D5m1m4 = cplxe(1,  -$f) *  (sqrt(21)/16)*(sin($t)**3)*(1+cos($t))*(5*cos($t)-1);
my $D5p1p5 = cplxe(1,   $f) *  (sqrt(210)/32)*(sin($t)**4)*(1+cos($t));
my $D5p1m5 = cplxe(1,   $f) *  (sqrt(210)/32)*(sin($t)**4)*(1-cos($t));
my $D5m1p5 = cplxe(1,  -$f) *  (sqrt(210)/32)*(sin($t)**4)*(1-cos($t));
my $D5m1m5 = cplxe(1,  -$f) *  (sqrt(210)/32)*(sin($t)**4)*(1+cos($t));

my $D5p2p0 = cplxe(1, 2*$f) *  (sqrt(210)/8)*(sin($t)**2)*cos($t)*(3*cos($t)**2-1);
my $D5m2p0 = cplxe(1,-2*$f) *  (sqrt(210)/8)*(sin($t)**2)*cos($t)*(3*cos($t)**2-1);
my $D5p2p1 = cplxe(1, 2*$f) * -(sqrt(7)/8)*sin($t)*(1+cos($t))*(15*cos($t)**3-9*cos($t)**2-1);
my $D5p2m1 = cplxe(1, 2*$f) * -(sqrt(7)/8)*sin($t)*(1-cos($t))*(15*cos($t)**3+9*cos($t)**2+1);   # Something fishy here... see Buckmaster
my $D5m2p1 = cplxe(1,-2*$f) *  (sqrt(7)/8)*sin($t)*(1-cos($t))*(15*cos($t)**3+9*cos($t)**2+1);
my $D5m2m1 = cplxe(1,-2*$f) *  (sqrt(7)/8)*sin($t)*(1+cos($t))*(15*cos($t)**3-9*cos($t)**2-1);   # Something fishy here... see Buckmaster
my $D5p2p2 = cplxe(1, 2*$f) *  (1/4)*((1+cos($t))**2)*(15*cos($t)**3-18*cos($t)**2+3*cos($t)+1);
my $D5p2m2 = cplxe(1, 2*$f) *  (1/4)*((1-cos($t))**2)*(15*cos($t)**3+18*cos($t)**2+3*cos($t)-1);
my $D5m2p2 = cplxe(1,-2*$f) *  (1/4)*((1-cos($t))**2)*(15*cos($t)**3+18*cos($t)**2+3*cos($t)-1);
my $D5m2m2 = cplxe(1,-2*$f) *  (1/4)*((1+cos($t))**2)*(15*cos($t)**3-18*cos($t)**2+3*cos($t)+1);
my $D5p2p3 = cplxe(1, 2*$f) *  (sqrt(6)/16)*sin($t)*((1+cos($t))**2)*(15*cos($t)**2-12*cos($t)+1);
my $D5p2m3 = cplxe(1, 2*$f) *  (sqrt(6)/16)*sin($t)*((1-cos($t))**2)*(15*cos($t)**2+12*cos($t)+1);
my $D5m2p3 = cplxe(1,-2*$f) *  (sqrt(6)/16)*sin($t)*((1-cos($t))**2)*(15*cos($t)**2+12*cos($t)+1);
my $D5m2m3 = cplxe(1,-2*$f) *  (sqrt(6)/16)*sin($t)*((1+cos($t))**2)*(15*cos($t)**2-12*cos($t)+1);
my $D5p2p4 = cplxe(1, 2*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1+cos($t))**2)*(5*cos($t)**2-2);
my $D5p2m4 = cplxe(1, 2*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1-cos($t))**2)*(5*cos($t)**2+2);
my $D5m2p4 = cplxe(1,-2*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1-cos($t))**2)*(5*cos($t)**2+2);
my $D5m2m4 = cplxe(1,-2*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1+cos($t))**2)*(5*cos($t)**2-2);
my $D5p2p5 = cplxe(1, 2*$f) *  (sqrt(30)/16)*(sin($t)**3)*((1+cos($t))**2);
my $D5p2m5 = cplxe(1, 2*$f) *  (sqrt(30)/16)*(sin($t)**3)*((1-cos($t))**2);
my $D5m2p5 = cplxe(1,-2*$f) *  (sqrt(30)/16)*(sin($t)**3)*((1-cos($t))**2);
my $D5m2m5 = cplxe(1,-2*$f) *  (sqrt(30)/16)*(sin($t)**3)*((1+cos($t))**2);

my $D5p3p0 = cplxe(1, 3*$f) * -(sqrt(35)/16)*(sin($t)**3)*(9*cos($t)**2-1);
my $D5m3p0 = cplxe(1,-3*$f) *  (sqrt(35)/16)*(sin($t)**3)*(9*cos($t)**2-1);
my $D5p3p1 = cplxe(1, 3*$f) *  (sqrt(42)/32)*(sin($t)**2)*(1+cos($t))*(15*cos($t)**2-6*cos($t)-1);
my $D5p3m1 = cplxe(1, 3*$f) *  (sqrt(42)/32)*(sin($t)**2)*(1-cos($t))*(15*cos($t)**2+6*cos($t)-1);
my $D5m3p1 = cplxe(1,-3*$f) *  (sqrt(42)/32)*(sin($t)**2)*(1-cos($t))*(15*cos($t)**2+6*cos($t)-1);
my $D5m3m1 = cplxe(1,-3*$f) *  (sqrt(42)/32)*(sin($t)**2)*(1+cos($t))*(15*cos($t)**2-6*cos($t)-1);
my $D5p3p2 = cplxe(1, 3*$f) * -(sqrt(6)/16)*sin($t)*((1+cos($t))**2)*(15*cos($t)**2-12*cos($t)+1);
my $D5p3m2 = cplxe(1, 3*$f) * -(sqrt(6)/16)*sin($t)*((1-cos($t))**2)*(15*cos($t)**2+12*cos($t)+1);
my $D5m3p2 = cplxe(1,-3*$f) *  (sqrt(6)/16)*sin($t)*((1-cos($t))**2)*(15*cos($t)**2+12*cos($t)+1);
my $D5m3m2 = cplxe(1,-3*$f) *  (sqrt(6)/16)*sin($t)*((1+cos($t))**2)*(15*cos($t)**2-12*cos($t)+1);
my $D5p3p3 = cplxe(1, 3*$f) *  (1/32)*((1+cos($t))**3)*(45*cos($t)**2-54*cos($t)+13);
my $D5p3m3 = cplxe(1, 3*$f) *  (1/32)*((1-cos($t))**3)*(45*cos($t)**2+54*cos($t)+13);
my $D5m3p3 = cplxe(1,-3*$f) *  (1/32)*((1-cos($t))**3)*(45*cos($t)**2+54*cos($t)+13);
my $D5m3m3 = cplxe(1,-3*$f) *  (1/32)*((1+cos($t))**3)*(45*cos($t)**2-54*cos($t)+13);
my $D5p3p4 = cplxe(1, 3*$f) *  (3*sqrt(2)/32)*sin($t)*((1+cos($t))**3)*(5*cos($t)-3);
my $D5p3m4 = cplxe(1, 3*$f) *  (3*sqrt(2)/32)*sin($t)*((1-cos($t))**3)*(5*cos($t)+3);
my $D5m3p4 = cplxe(1,-3*$f) *  (3*sqrt(2)/32)*sin($t)*((1-cos($t))**3)*(5*cos($t)+3);
my $D5m3m4 = cplxe(1,-3*$f) *  (3*sqrt(2)/32)*sin($t)*((1+cos($t))**3)*(5*cos($t)-3);
my $D5p3p5 = cplxe(1, 3*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1+cos($t))**3);
my $D5p3m5 = cplxe(1, 3*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1-cos($t))**3);
my $D5m3p5 = cplxe(1,-3*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1-cos($t))**3);
my $D5m3m5 = cplxe(1,-3*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1+cos($t))**3);

my $D5p4p0 = cplxe(1, 4*$f) *  (3*sqrt(70)/16)*(sin($t)**4)*cos($t);
my $D5m4p0 = cplxe(1,-4*$f) *  (3*sqrt(70)/16)*(sin($t)**4)*cos($t);
my $D5p4p1 = cplxe(1, 4*$f) * -(sqrt(21)/16)*(sin($t)**3)*(1+cos($t))*(5*cos($t)-1);
my $D5p4m1 = cplxe(1, 4*$f) * -(sqrt(21)/16)*(sin($t)**3)*(1-cos($t))*(5*cos($t)+1);
my $D5m4p1 = cplxe(1,-4*$f) *  (sqrt(21)/16)*(sin($t)**3)*(1-cos($t))*(5*cos($t)+1);
my $D5m4m1 = cplxe(1,-4*$f) *  (sqrt(21)/16)*(sin($t)**3)*(1+cos($t))*(5*cos($t)-1);
my $D5p4p2 = cplxe(1, 4*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1+cos($t))**2)*(5*cos($t)-2);
my $D5p4m2 = cplxe(1, 4*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1-cos($t))**2)*(5*cos($t)+2);
my $D5m4p2 = cplxe(1,-4*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1-cos($t))**2)*(5*cos($t)+2);
my $D5m4m2 = cplxe(1,-4*$f) *  (sqrt(3)/8)*(sin($t)**2)*((1+cos($t))**2)*(5*cos($t)-2);
my $D5p4p3 = cplxe(1, 4*$f) * -(3*sqrt(2)/32)*sin($t)*((1+cos($t))**3)*(5*cos($t)-3);
my $D5p4m3 = cplxe(1, 4*$f) * -(3*sqrt(2)/32)*sin($t)*((1-cos($t))**3)*(5*cos($t)+3);
my $D5m4p3 = cplxe(1,-4*$f) *  (3*sqrt(2)/32)*sin($t)*((1-cos($t))**3)*(5*cos($t)+3);
my $D5m4m3 = cplxe(1,-4*$f) *  (3*sqrt(2)/32)*sin($t)*((1+cos($t))**3)*(5*cos($t)-3);
my $D5p4p4 = cplxe(1, 4*$f) *  (1/16)*((1+cos($t))**4)*(5*cos($t)-4);
my $D5p4m4 = cplxe(1, 4*$f) *  (1/16)*((1-cos($t))**4)*(5*cos($t)+4);
my $D5m4p4 = cplxe(1,-4*$f) *  (1/16)*((1-cos($t))**4)*(5*cos($t)+4);
my $D5m4m4 = cplxe(1,-4*$f) *  (1/16)*((1+cos($t))**4)*(5*cos($t)-4);
my $D5p4p5 = cplxe(1, 4*$f) *  (sqrt(10)/32)*sin($t)*((1+cos($t))**4);
my $D5p4m5 = cplxe(1, 4*$f) *  (sqrt(10)/32)*sin($t)*((1-cos($t))**4);
my $D5m4p5 = cplxe(1,-4*$f) *  (sqrt(10)/32)*sin($t)*((1-cos($t))**4);
my $D5m4m5 = cplxe(1,-4*$f) *  (sqrt(10)/32)*sin($t)*((1+cos($t))**4);

my $D5p5p0 = cplxe(1, 5*$f) * -(3*sqrt(7)/16)*sin($t)**5;
my $D5m5p0 = cplxe(1,-5*$f) *  (3*sqrt(7)/16)*sin($t)**5;
my $D5p5p1 = cplxe(1, 5*$f) *  (sqrt(210)/32)*(sin($t)**4)*(1+cos($t));
my $D5p5m1 = cplxe(1, 5*$f) *  (sqrt(210)/32)*(sin($t)**4)*(1-cos($t));
my $D5m5p1 = cplxe(1,-5*$f) *  (sqrt(210)/32)*(sin($t)**4)*(1-cos($t));
my $D5m5m1 = cplxe(1,-5*$f) *  (sqrt(210)/32)*(sin($t)**4)*(1+cos($t));
my $D5p5p2 = cplxe(1, 5*$f) * -(sqrt(30)/16)*(sin($t)**3)*((1+cos($t))**2);
my $D5p5m2 = cplxe(1, 5*$f) * -(sqrt(30)/16)*(sin($t)**3)*((1-cos($t))**2);
my $D5m5p2 = cplxe(1,-5*$f) *  (sqrt(30)/16)*(sin($t)**3)*((1-cos($t))**2);
my $D5m5m2 = cplxe(1,-5*$f) *  (sqrt(30)/16)*(sin($t)**3)*((1+cos($t))**2);
my $D5p5p3 = cplxe(1, 5*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1+cos($t))**3);
my $D5p5m3 = cplxe(1, 5*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1-cos($t))**3);
my $D5m5p3 = cplxe(1,-5*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1-cos($t))**3);
my $D5m5m3 = cplxe(1,-5*$f) *  (3*sqrt(5)/32)*(sin($t)**2)*((1+cos($t))**3);
my $D5p5p4 = cplxe(1, 5*$f) *  (sqrt(10)/32)*sin($t)*((1+cos($t))**4);
my $D5p5m4 = cplxe(1, 5*$f) *  (sqrt(10)/32)*sin($t)*((1-cos($t))**4);
my $D5m5p4 = cplxe(1,-5*$f) *  (sqrt(10)/32)*sin($t)*((1-cos($t))**4);
my $D5m5m4 = cplxe(1,-5*$f) *  (sqrt(10)/32)*sin($t)*((1+cos($t))**4);
my $D5p5p5 = cplxe(1, 5*$f) *  (1/32)*((1+cos($t))**5);
my $D5p5m5 = cplxe(1, 5*$f) *  (1/32)*((1-cos($t))**5);
my $D5m5p5 = cplxe(1,-5*$f) *  (1/32)*((1-cos($t))**5);
my $D5m5m5 = cplxe(1,-5*$f) *  (1/32)*((1+cos($t))**5);

$D = [
             [ $D5m5m5,  $D5m5m4,  $D5m5m3,  $D5m5m2,  $D5m5m1,  $D5m5p0,  $D5m5p1,  $D5m5p2,  $D5m5p3,  $D5m5p4,  $D5m5p5], 
             [ $D5m4m5,  $D5m4m4,  $D5m4m3,  $D5m4m2,  $D5m4m1,  $D5m4p0,  $D5m4p1,  $D5m4p2,  $D5m4p3,  $D5m4p4,  $D5m4p5], 
             [ $D5m3m5,  $D5m3m4,  $D5m3m3,  $D5m3m2,  $D5m3m1,  $D5m3p0,  $D5m3p1,  $D5m3p2,  $D5m3p3,  $D5m3p4,  $D5m3p5], 
             [ $D5m2m5,  $D5m2m4,  $D5m2m3,  $D5m2m2,  $D5m2m1,  $D5m2p0,  $D5m2p1,  $D5m2p2,  $D5m2p3,  $D5m2p4,  $D5m2p5], 
             [ $D5m1m5,  $D5m1m4,  $D5m1m3,  $D5m1m2,  $D5m1m1,  $D5m1p0,  $D5m1p1,  $D5m1p2,  $D5m1p3,  $D5m1p4,  $D5m1p5], 
             [-$D50p5 ,  $D50p4 , -$D50p3 ,  $D50p2 , -$D50p1 ,  $D50p0 ,  $D50p1 ,  $D50p2 ,  $D50p3 ,  $D50p4 ,  $D50p5 ], 
             [ $D5p1m5,  $D5p1m4,  $D5p1m3,  $D5p1m2,  $D5p1m1,  $D5p1p0,  $D5p1p1,  $D5p1p2,  $D5p1p3,  $D5p1p4,  $D5p1p5], 
             [ $D5p2m5,  $D5p2m4,  $D5p2m3,  $D5p2m2,  $D5p2m1,  $D5p2p0,  $D5p2p1,  $D5p2p2,  $D5p2p3,  $D5p2p4,  $D5p2p5], 
             [ $D5p3m5,  $D5p3m4,  $D5p3m3,  $D5p3m2,  $D5p3m1,  $D5p3p0,  $D5p3p1,  $D5p3p2,  $D5p3p3,  $D5p3p4,  $D5p3p5], 
             [ $D5p4m5,  $D5p4m4,  $D5p4m3,  $D5p4m2,  $D5p4m1,  $D5p4p0, -$D5p4p1,  $D5p4p2,  $D5p4p3,  $D5p4p4,  $D5p4p5], 
             [ $D5p5m5,  $D5p5m4,  $D5p5m3,  $D5p5m2,  $D5p5m1, -$D5p5p0,  $D5p5p1,  $D5p5p2,  $D5p5p3,  $D5p5p4,  $D5p5p5]
         ];
return $D;
}

# --------------------------------------- O6 ------------------------------------- #
sub rot6 { my ($t,$f) = @_; my $D;

my $D60p0 =  (1/16)          * (231*cos($t)**6 - 315*cos($t)**4 + 105*cos($t)**2 - 5);
my $D60p1 =  sqrt(6*7)/16    * sin($t)    * cos($t) * (33*cos($t)**4 - 30*cos($t)**2 + 5);
my $D60p2 =  sqrt(3*5*7)/32  * sin($t)**2 * (33*cos($t)**4 - 18*cos($t)**2 + 1);
my $D60p3 =  sqrt(3*5*7)/16  * sin($t)**3 * cos($t) * (11*cos($t)**2 - 3);
my $D60p4 =  3*sqrt(2*7)/32  * sin($t)**4 * (11*cos($t)**2 - 1);
my $D60p5 =  3*sqrt(7*11)/16 * sin($t)**5 * cos($t);
my $D60p6 =  sqrt(3*7*11)/32 * sin($t)**6;

my $D6p1p0 = -cplxe(1, $f) * sqrt(6*7)/16    * sin($t)    * cos($t)     * (33*cos($t)**4 - 30*cos($t)**2 + 5);
my $D6p1p1 =  cplxe(1, $f) * (1/16) * (1+cos($t)) * (198*cos($t)**5 - 165*cos($t)**4 - 120*cos($t)**3 + 90*cos($t)**2 + 10*cos($t) - 5);
my $D6p1m1 =  cplxe(1, $f) * (1/16) * (1-cos($t)) * (198*cos($t)**5 + 165*cos($t)**4 - 120*cos($t)**3 - 90*cos($t)**2 + 10*cos($t) + 5);
my $D6p1p2 =  cplxe(1, $f) * sqrt(2*5)/32    * sin($t)    * (1+cos($t)) * (99*cos($t)**4 - 66*cos($t)**3 - 36*cos($t)**2 + 18*cos($t) + 1);
my $D6p1m2 = -cplxe(1, $f) * sqrt(2*5)/32    * sin($t)    * (1-cos($t)) * (99*cos($t)**4 + 66*cos($t)**3 - 36*cos($t)**2 - 18*cos($t) + 1);
my $D6p1p3 =  cplxe(1, $f) * 3*sqrt(2*5)/32  * sin($t)**2 * (1+cos($t)) * (22*cos($t)**3 - 11*cos($t)**2 - 4*cos($t) + 1);
my $D6p1m3 =  cplxe(1, $f) * 3*sqrt(2*5)/32  * sin($t)**2 * (1-cos($t)) * (22*cos($t)**3 + 11*cos($t)**2 - 4*cos($t) - 1);
my $D6p1p4 =  cplxe(1, $f) * sqrt(3)/16      * sin($t)**3 * (1+cos($t)) * (33*cos($t)**2 - 11*cos($t) - 2);
my $D6p1m4 = -cplxe(1, $f) * sqrt(3)/16      * sin($t)**3 * (1-cos($t)) * (33*cos($t)**2 + 11*cos($t) - 2);
my $D6p1p5 =  cplxe(1, $f) * sqrt(6*11)/32   * sin($t)**4 * (1+cos($t)) * (6*cos($t) - 1);
my $D6p1m5 =  cplxe(1, $f) * sqrt(6*11)/32   * sin($t)**4 * (1-cos($t)) * (6*cos($t) + 1);
my $D6p1p6 =  cplxe(1, $f) * 3*sqrt(2*11)/32 * sin($t)**5 * (1+cos($t));
my $D6p1m6 = -cplxe(1, $f) * 3*sqrt(2*11)/32 * sin($t)**5 * (1-cos($t));

my $D6m1p0 =  cplxe(1,-$f) * sqrt(6*7)/16    * sin($t)    * cos($t)     * (33*cos($t)**4 - 30*cos($t)**2 + 5);
my $D6m1p1 =  cplxe(1,-$f) * (1/16) * (1-cos($t)) * (198*cos($t)**5 + 165*cos($t)**4 - 120*cos($t)**3 - 90*cos($t)**2 + 10*cos($t) + 5);
my $D6m1m1 =  cplxe(1,-$f) * (1/16) * (1+cos($t)) * (198*cos($t)**5 - 165*cos($t)**4 - 120*cos($t)**3 + 90*cos($t)**2 + 10*cos($t) - 5);
my $D6m1p2 =  cplxe(1,-$f) * sqrt(2*5)/32    * sin($t)    * (1-cos($t)) * (99*cos($t)**4 + 66*cos($t)**3 - 36*cos($t)**2 - 18*cos($t) + 1);
my $D6m1m2 = -cplxe(1,-$f) * sqrt(2*5)/32    * sin($t)    * (1+cos($t)) * (99*cos($t)**4 - 66*cos($t)**3 - 36*cos($t)**2 + 18*cos($t) + 1);
my $D6m1p3 =  cplxe(1,-$f) * 3*sqrt(2*5)/32  * sin($t)**2 * (1-cos($t)) * (22*cos($t)**3 + 11*cos($t)**2 - 4*cos($t) - 1);
my $D6m1m3 =  cplxe(1,-$f) * 3*sqrt(2*5)/32  * sin($t)**2 * (1+cos($t)) * (22*cos($t)**3 - 11*cos($t)**2 - 4*cos($t) + 1);
my $D6m1p4 =  cplxe(1,-$f) * sqrt(3)/16      * sin($t)**3 * (1-cos($t)) * (33*cos($t)**2 + 11*cos($t) - 2);
my $D6m1m4 = -cplxe(1,-$f) * sqrt(3)/16      * sin($t)**3 * (1+cos($t)) * (33*cos($t)**2 - 11*cos($t) - 2);
my $D6m1p5 =  cplxe(1,-$f) * sqrt(6*11)/32   * sin($t)**4 * (1-cos($t)) * (6*cos($t) + 1);
my $D6m1m5 =  cplxe(1,-$f) * sqrt(6*11)/32   * sin($t)**4 * (1+cos($t)) * (6*cos($t) - 1);
my $D6m1p6 =  cplxe(1,-$f) * 3*sqrt(2*11)/32 * sin($t)**5 * (1-cos($t));
my $D6m1m6 = -cplxe(1,-$f) * 3*sqrt(2*11)/32 * sin($t)**5 * (1+cos($t));

my $D6p2p0 =  cplxe(1, 2*$f) * sqrt(3*5*7)/32   * sin($t)**2 * (33*cos($t)**4 - 18*cos($t)**2 + 1);
my $D6p2p1 = -cplxe(1, 2*$f) * sqrt(2*5)/32     * sin($t)    * (1+cos($t))    * (99*cos($t)**4 - 66*cos($t)**3 - 36*cos($t)**2 + 18*cos($t) + 1); 
my $D6p2m1 = -cplxe(1, 2*$f) * sqrt(2*5)/32     * sin($t)    * (1-cos($t))    * (99*cos($t)**4 + 66*cos($t)**3 - 36*cos($t)**2 - 18*cos($t) + 1); 
my $D6p2p2 =  cplxe(1, 2*$f) * (1/64) *  (1+cos($t))**2 * (495*cos($t)**4 - 660*cos($t)**3 + 90*cos($t)**2 + 108*cos($t) - 17);
my $D6p2m2 =  cplxe(1, 2*$f) * (1/64) *  (1-cos($t))**2 * (495*cos($t)**4 + 660*cos($t)**3 + 90*cos($t)**2 - 108*cos($t) - 17);
my $D6p2p3 =  cplxe(1, 2*$f) * (3/32)           * sin($t)    * (1+cos($t))**2 * (55*cos($t)**3 - 55*cos($t)**2 + 5*cos($t) + 3);
my $D6p2m3 = -cplxe(1, 2*$f) * (3/32)           * sin($t)    * (1-cos($t))**2 * (55*cos($t)**3 + 55*cos($t)**2 + 5*cos($t) - 3);
my $D6p2p4 =  cplxe(1, 2*$f) * sqrt(6*5)/64     * sin($t)**2 * (1+cos($t))**2 * (33*cos($t)**2 - 22*cos($t) + 1);
my $D6p2m4 =  cplxe(1, 2*$f) * sqrt(6*5)/64     * sin($t)**2 * (1-cos($t))**2 * (33*cos($t)**2 + 22*cos($t) + 1);
my $D6p2p5 =  cplxe(1, 2*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1+cos($t))**2 * (3*cos($t) - 1);
my $D6p2m5 = -cplxe(1, 2*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1-cos($t))**2 * (3*cos($t) + 1);
my $D6p2p6 =  cplxe(1, 2*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1+cos($t))**2;
my $D6p2m6 =  cplxe(1, 2*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1-cos($t))**2;

my $D6m2p0 =  cplxe(1,-2*$f) * sqrt(3*5*7)/32   * sin($t)**2 * (33*cos($t)**4 - 18*cos($t)**2 + 1);
my $D6m2p1 =  cplxe(1,-2*$f) * sqrt(2*5)/32     * sin($t)    * (1-cos($t))    * (99*cos($t)**4 + 66*cos($t)**3 - 36*cos($t)**2 - 18*cos($t) + 1); 
my $D6m2m1 =  cplxe(1,-2*$f) * sqrt(2*5)/32     * sin($t)    * (1+cos($t))    * (99*cos($t)**4 - 66*cos($t)**3 - 36*cos($t)**2 + 18*cos($t) + 1); 
my $D6m2p2 =  cplxe(1,-2*$f) * (1/64) *  (1-cos($t))**2 * (495*cos($t)**4 + 660*cos($t)**3 + 90*cos($t)**2 - 108*cos($t) - 17);
my $D6m2m2 =  cplxe(1,-2*$f) * (1/64) *  (1+cos($t))**2 * (495*cos($t)**4 - 660*cos($t)**3 + 90*cos($t)**2 + 108*cos($t) - 17);
my $D6m2p3 =  cplxe(1,-2*$f) * (3/32)           * sin($t)    * (1-cos($t))**2 * (55*cos($t)**3 + 55*cos($t)**2 + 5*cos($t) - 3);
my $D6m2m3 = -cplxe(1,-2*$f) * (3/32)           * sin($t)    * (1+cos($t))**2 * (55*cos($t)**3 - 55*cos($t)**2 + 5*cos($t) + 3);
my $D6m2p4 =  cplxe(1,-2*$f) * sqrt(6*5)/64     * sin($t)**2 * (1-cos($t))**2 * (33*cos($t)**2 + 22*cos($t) + 1);
my $D6m2m4 =  cplxe(1,-2*$f) * sqrt(6*5)/64     * sin($t)**2 * (1+cos($t))**2 * (33*cos($t)**2 - 22*cos($t) + 1);
my $D6m2p5 =  cplxe(1,-2*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1-cos($t))**2 * (3*cos($t) + 1);
my $D6m2m5 = -cplxe(1,-2*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1+cos($t))**2 * (3*cos($t) - 1);
my $D6m2p6 =  cplxe(1,-2*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1-cos($t))**2;
my $D6m2m6 =  cplxe(1,-2*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1+cos($t))**2;

my $D6p3p0 = -cplxe(1, 3*$f) * sqrt(3*5*11)/16  * sin($t)**3 * cos($t)        * (11*cos($t)**2 - 3);
my $D6p3p1 =  cplxe(1, 3*$f) * 3*sqrt(2*5)/32   * sin($t)**2 * (1+cos($t))    * (22*cos($t)**3 - 11*cos($t)**2 - 4*cos($t) + 1);
my $D6p3m1 =  cplxe(1, 3*$f) * 3*sqrt(2*5)/32   * sin($t)**2 * (1-cos($t))    * (22*cos($t)**3 + 11*cos($t)**2 - 4*cos($t) - 1);
my $D6p3p2 = -cplxe(1, 3*$f) * (3/32)           * sin($t)    * (1+cos($t))**2 * (55*cos($t)**3 - 55*cos($t)**2 + 5*cos($t) + 3);
my $D6p3m2 = -cplxe(1, 3*$f) * (3/32)           * sin($t)    * (1-cos($t))**2 * (55*cos($t)**3 + 55*cos($t)**2 + 5*cos($t) - 3);
my $D6p3p3 =  cplxe(1, 3*$f) * (1/32)           *           (1+cos($t))**3 * (110*cos($t)**3 - 165*cos($t)**2 + 60*cos($t) - 1);
my $D6p3m3 =  cplxe(1, 3*$f) * (1/32)           *           (1-cos($t))**3 * (110*cos($t)**3 + 165*cos($t)**2 + 60*cos($t) + 1);
my $D6p3p4 =  cplxe(1, 3*$f) * sqrt(5*6)/32     * sin($t)    * (1+cos($t))**3 * (11*cos($t)**2 - 11*cos($t) + 2);
my $D6p3m4 = -cplxe(1, 3*$f) * sqrt(5*6)/32     * sin($t)    * (1-cos($t))**3 * (11*cos($t)**2 + 11*cos($t) + 2);
my $D6p3p5 =  cplxe(1, 3*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1+cos($t))**3 * (2*cos($t) - 1);
my $D6p3m5 =  cplxe(1, 3*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1-cos($t))**3 * (2*cos($t) + 1);
my $D6p3p6 =  cplxe(1, 3*$f) * sqrt(5*11)/32    * sin($t)**3 * (1+cos($t))**3;
my $D6p3m6 = -cplxe(1, 3*$f) * sqrt(5*11)/32    * sin($t)**3 * (1-cos($t))**3;

my $D6m3p0 =  cplxe(1,-3*$f) * sqrt(3*5*11)/16  * sin($t)**3 * cos($t)        * (11*cos($t)**2 - 3);
my $D6m3p1 =  cplxe(1,-3*$f) * 3*sqrt(2*5)/32   * sin($t)**2 * (1-cos($t))    * (22*cos($t)**3 + 11*cos($t)**2 - 4*cos($t) - 1);
my $D6m3m1 =  cplxe(1,-3*$f) * 3*sqrt(2*5)/32   * sin($t)**2 * (1+cos($t))    * (22*cos($t)**3 - 11*cos($t)**2 - 4*cos($t) + 1);
my $D6m3p2 =  cplxe(1,-3*$f) * (3/32)           * sin($t)    * (1-cos($t))**2 * (55*cos($t)**3 + 55*cos($t)**2 + 5*cos($t) - 3);
my $D6m3m2 =  cplxe(1,-3*$f) * (3/32)           * sin($t)    * (1+cos($t))**2 * (55*cos($t)**3 - 55*cos($t)**2 + 5*cos($t) + 3);
my $D6m3p3 =  cplxe(1,-3*$f) * (1/32)           *           (1-cos($t))**3 * (110*cos($t)**3 + 165*cos($t)**2 + 60*cos($t) + 1);
my $D6m3m3 =  cplxe(1,-3*$f) * (1/32)           *           (1+cos($t))**3 * (110*cos($t)**3 - 165*cos($t)**2 + 60*cos($t) - 1);
my $D6m3p4 =  cplxe(1,-3*$f) * sqrt(5*6)/32     * sin($t)    * (1-cos($t))**3 * (11*cos($t)**2 + 11*cos($t) + 2);
my $D6m3m4 = -cplxe(1,-3*$f) * sqrt(5*6)/32     * sin($t)    * (1+cos($t))**3 * (11*cos($t)**2 - 11*cos($t) + 2);
my $D6m3p5 =  cplxe(1,-3*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1-cos($t))**3 * (2*cos($t) + 1);
my $D6m3m5 =  cplxe(1,-3*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1+cos($t))**3 * (2*cos($t) - 1);
my $D6m3p6 =  cplxe(1,-3*$f) * sqrt(5*11)/32    * sin($t)**3 * (1-cos($t))**3;
my $D6m3m6 = -cplxe(1,-3*$f) * sqrt(5*11)/32    * sin($t)**3 * (1+cos($t))**3;
 
my $D6p4p0 =  cplxe(1, 4*$f) * 3*sqrt(2*7)/32   * sin($t)**4 * (11*cos($t)**2 - 1);
my $D6p4p1 = -cplxe(1, 4*$f) * sqrt(3)/16       * sin($t)**3 * (1+cos($t))    * (33*cos($t)**2 - 11*cos($t) - 2);
my $D6p4m1 = -cplxe(1, 4*$f) * sqrt(3)/16       * sin($t)**3 * (1-cos($t))    * (33*cos($t)**2 + 11*cos($t) - 2);
my $D6p4p2 =  cplxe(1, 4*$f) * sqrt(6*5)/64     * sin($t)**2 * (1+cos($t))**2 * (33*cos($t)**2 - 22*cos($t) + 1);
my $D6p4m2 =  cplxe(1, 4*$f) * sqrt(6*5)/64     * sin($t)**2 * (1-cos($t))**2 * (33*cos($t)**2 + 22*cos($t) + 1);
my $D6p4p3 = -cplxe(1, 4*$f) * sqrt(6*5)/32     * sin($t)    * (1+cos($t))**3 * (11*cos($t)**2 - 11*cos($t) + 2);
my $D6p4m3 = -cplxe(1, 4*$f) * sqrt(6*5)/32     * sin($t)    * (1-cos($t))**3 * (11*cos($t)**2 + 11*cos($t) + 2);
my $D6p4p4 =  cplxe(1, 4*$f) * (1/32)           *              (1+cos($t))**4 * (33*cos($t)**2 - 44*cos($t) + 13);
my $D6p4m4 =  cplxe(1, 4*$f) * (1/32)           *              (1-cos($t))**4 * (33*cos($t)**2 + 44*cos($t) + 13);
my $D6p4p5 =  cplxe(1, 4*$f) * sqrt(2*11)/32    * sin($t)    * (1+cos($t))**4 * (3*cos($t) - 2);
my $D6p4m5 = -cplxe(1, 4*$f) * sqrt(2*11)/32    * sin($t)    * (1-cos($t))**4 * (3*cos($t) + 2);
my $D6p4p6 =  cplxe(1, 4*$f) * sqrt(6*11)/64    * sin($t)**2 * (1+cos($t))**4;
my $D6p4m6 =  cplxe(1, 4*$f) * sqrt(6*11)/64    * sin($t)**2 * (1-cos($t))**4;

my $D6m4p0 =  cplxe(1,-4*$f) * 3*sqrt(2*7)/32   * sin($t)**4 * (11*cos($t)**2 - 1);
my $D6m4p1 =  cplxe(1,-4*$f) * sqrt(3)/16       * sin($t)**3 * (1-cos($t))    * (33*cos($t)**2 + 11*cos($t) - 2);
my $D6m4m1 =  cplxe(1,-4*$f) * sqrt(3)/16       * sin($t)**3 * (1+cos($t))    * (33*cos($t)**2 - 11*cos($t) - 2);
my $D6m4p2 =  cplxe(1,-4*$f) * sqrt(6*5)/64     * sin($t)**2 * (1-cos($t))**2 * (33*cos($t)**2 + 22*cos($t) + 1);
my $D6m4m2 =  cplxe(1,-4*$f) * sqrt(6*5)/64     * sin($t)**2 * (1+cos($t))**2 * (33*cos($t)**2 - 22*cos($t) + 1);
my $D6m4p3 =  cplxe(1,-4*$f) * sqrt(6*5)/32     * sin($t)    * (1-cos($t))**3 * (11*cos($t)**2 + 11*cos($t) + 2);
my $D6m4m3 =  cplxe(1,-4*$f) * sqrt(6*5)/32     * sin($t)    * (1+cos($t))**3 * (11*cos($t)**2 - 11*cos($t) + 2);
my $D6m4p4 =  cplxe(1,-4*$f) * (1/32)           *              (1-cos($t))**4 * (33*cos($t)**2 + 44*cos($t) + 13);
my $D6m4m4 =  cplxe(1,-4*$f) * (1/32)           *              (1+cos($t))**4 * (33*cos($t)**2 - 44*cos($t) + 13);
my $D6m4p5 =  cplxe(1,-4*$f) * sqrt(2*11)/32    * sin($t)    * (1-cos($t))**4 * (3*cos($t) + 2);
my $D6m4m5 = -cplxe(1,-4*$f) * sqrt(2*11)/32    * sin($t)    * (1+cos($t))**4 * (3*cos($t) - 2);
my $D6m4p6 =  cplxe(1,-4*$f) * sqrt(6*11)/64    * sin($t)**2 * (1-cos($t))**4;
my $D6m4m6 =  cplxe(1,-4*$f) * sqrt(6*11)/64    * sin($t)**2 * (1+cos($t))**4;

my $D6p5p0 = -cplxe(1, 5*$f) * 3*sqrt(7*11)/16  * sin($t)**5 * cos($t);
my $D6p5p1 =  cplxe(1, 5*$f) * sqrt(6*11)/32    * sin($t)**4 * (1+cos($t))    * (6*cos($t) - 1);
my $D6p5m1 =  cplxe(1, 5*$f) * sqrt(6*11)/32    * sin($t)**4 * (1-cos($t))    * (6*cos($t) + 1);
my $D6p5p2 = -cplxe(1, 5*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1+cos($t))**2 * (3*cos($t) - 1);
my $D6p5m2 = -cplxe(1, 5*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1-cos($t))**2 * (3*cos($t) + 1);
my $D6p5p3 =  cplxe(1, 5*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1+cos($t))**3 * (2*cos($t) - 1);
my $D6p5m3 =  cplxe(1, 5*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1-cos($t))**3 * (2*cos($t) + 1);
my $D6p5p4 = -cplxe(1, 5*$f) * sqrt(2*11)/32    * sin($t)    * (1+cos($t))**4 * (3*cos($t) - 2);
my $D6p5m4 = -cplxe(1, 5*$f) * sqrt(2*11)/32    * sin($t)    * (1-cos($t))**4 * (3*cos($t) + 2);
my $D6p5p5 =  cplxe(1, 5*$f) * (1/32)           *              (1+cos($t))**5 * (6*cos($t) - 5);
my $D6p5m5 =  cplxe(1, 5*$f) * (1/32)           *              (1-cos($t))**5 * (6*cos($t) + 5);
my $D6p5p6 =  cplxe(1, 5*$f) * sqrt(3)/32       * sin($t)    * (1+cos($t))**5;
my $D6p5m6 = -cplxe(1, 5*$f) * sqrt(3)/32       * sin($t)    * (1-cos($t))**5;

my $D6m5p0 =  cplxe(1,-5*$f) * 3*sqrt(7*11)/16  * sin($t)**5 * cos($t);
my $D6m5p1 =  cplxe(1,-5*$f) * sqrt(6*11)/32    * sin($t)**4 * (1-cos($t))    * (6*cos($t) + 1);
my $D6m5m1 =  cplxe(1,-5*$f) * sqrt(6*11)/32    * sin($t)**4 * (1+cos($t))    * (6*cos($t) - 1);
my $D6m5p2 =  cplxe(1,-5*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1-cos($t))**2 * (3*cos($t) + 1);
my $D6m5m2 =  cplxe(1,-5*$f) * sqrt(3*5*11)/32  * sin($t)**3 * (1+cos($t))**2 * (3*cos($t) - 1);
my $D6m5p3 =  cplxe(1,-5*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1-cos($t))**3 * (2*cos($t) + 1);
my $D6m5m3 =  cplxe(1,-5*$f) * sqrt(3*5*11)/32  * sin($t)**2 * (1+cos($t))**3 * (2*cos($t) - 1);
my $D6m5p4 =  cplxe(1,-5*$f) * sqrt(2*11)/32    * sin($t)    * (1-cos($t))**4 * (3*cos($t) + 2);
my $D6m5m4 =  cplxe(1,-5*$f) * sqrt(2*11)/32    * sin($t)    * (1+cos($t))**4 * (3*cos($t) - 2);
my $D6m5p5 =  cplxe(1,-5*$f) * (1/32)           *              (1-cos($t))**5 * (6*cos($t) + 5);
my $D6m5m5 =  cplxe(1,-5*$f) * (1/32)           *              (1+cos($t))**5 * (6*cos($t) - 5);
my $D6m5p6 =  cplxe(1,-5*$f) * sqrt(3)/32       * sin($t)    * (1-cos($t))**5;
my $D6m5m6 = -cplxe(1,-5*$f) * sqrt(3)/32       * sin($t)    * (1+cos($t))**5;

my $D6p6p0 =  cplxe(1, 6*$f) * sqrt(3*7*11)/32  * sin($t)**6;
my $D6p6p1 = -cplxe(1, 6*$f) * 3*sqrt(2*11)/32  * sin($t)**5 * (1+cos($t));
my $D6p6m1 = -cplxe(1, 6*$f) * 3*sqrt(2*11)/32  * sin($t)**5 * (1-cos($t));
my $D6p6p2 =  cplxe(1, 6*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1+cos($t))**2;
my $D6p6m2 =  cplxe(1, 6*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1-cos($t))**2;
my $D6p6p3 = -cplxe(1, 6*$f) * sqrt(5*11)/32    * sin($t)**3 * (1+cos($t))**3;
my $D6p6m3 = -cplxe(1, 6*$f) * sqrt(5*11)/32    * sin($t)**3 * (1-cos($t))**3;
my $D6p6p4 =  cplxe(1, 6*$f) * sqrt(6*11)/64    * sin($t)**2 * (1+cos($t))**4;
my $D6p6m4 =  cplxe(1, 6*$f) * sqrt(6*11)/64    * sin($t)**2 * (1-cos($t))**4;
my $D6p6p5 = -cplxe(1, 6*$f) * sqrt(3)/32       * sin($t)    * (1+cos($t))**5;
my $D6p6m5 = -cplxe(1, 6*$f) * sqrt(3)/32       * sin($t)    * (1-cos($t))**5;
my $D6p6p6 =  cplxe(1, 6*$f) * (1/64)           *              (1+cos($t))**6;
my $D6p6m6 =  cplxe(1, 6*$f) * (1/64)           *              (1-cos($t))**6;

my $D6m6p0 =  cplxe(1,-6*$f) * sqrt(3*7*11)/32  * sin($t)**6;
my $D6m6p1 =  cplxe(1,-6*$f) * 3*sqrt(2*11)/32  * sin($t)**5 * (1-cos($t));
my $D6m6m1 =  cplxe(1,-6*$f) * 3*sqrt(2*11)/32  * sin($t)**5 * (1+cos($t));
my $D6m6p2 =  cplxe(1,-6*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1-cos($t))**2;
my $D6m6m2 =  cplxe(1,-6*$f) * 3*sqrt(5*11)/64  * sin($t)**4 * (1+cos($t))**2;
my $D6m6p3 =  cplxe(1,-6*$f) * sqrt(5*11)/32    * sin($t)**3 * (1-cos($t))**3;
my $D6m6m3 =  cplxe(1,-6*$f) * sqrt(5*11)/32    * sin($t)**3 * (1+cos($t))**3;
my $D6m6p4 =  cplxe(1,-6*$f) * sqrt(6*11)/64    * sin($t)**2 * (1-cos($t))**4;
my $D6m6m4 =  cplxe(1,-6*$f) * sqrt(6*11)/64    * sin($t)**2 * (1+cos($t))**4;
my $D6m6p5 =  cplxe(1,-6*$f) * sqrt(3)/32       * sin($t)    * (1-cos($t))**5;
my $D6m6m5 =  cplxe(1,-6*$f) * sqrt(3)/32       * sin($t)    * (1+cos($t))**5;
my $D6m6p6 =  cplxe(1,-6*$f) * (1/64)           *              (1-cos($t))**6;
my $D6m6m6 =  cplxe(1,-6*$f) * (1/64)           *              (1+cos($t))**6;

$D = [
     [$D6m6m6, $D6m6m5, $D6m6m4, $D6m6m3, $D6m6m2, $D6m6m1, $D6m6m0, $D6m6p1, $D6m6p2, $D6m6p3, $D6m6p4, $D6m6p5, $D6m6p6],
     [$D6m5m6, $D6m5m5, $D6m5m4, $D6m5m3, $D6m5m2, $D6m5m1, $D6m5m0, $D6m5p1, $D6m5p2, $D6m5p3, $D6m5p4, $D6m5p5, $D6m5p6],
     [$D6m4m6, $D6m4m5, $D6m4m4, $D6m4m3, $D6m4m2, $D6m4m1, $D6m4m0, $D6m4p1, $D6m4p2, $D6m4p3, $D6m4p4, $D6m4p5, $D6m4p6],
     [$D6m3m6, $D6m3m5, $D6m3m4, $D6m3m3, $D6m3m2, $D6m3m1, $D6m3m0, $D6m3p1, $D6m3p2, $D6m3p3, $D6m3p4, $D6m3p5, $D6m3p6],
     [$D6m2m6, $D6m2m5, $D6m2m4, $D6m2m3, $D6m2m2, $D6m2m1, $D6m2m0, $D6m2p1, $D6m2p2, $D6m2p3, $D6m2p4, $D6m2p5, $D6m2p6],
     [$D6m1m6, $D6m1m5, $D6m1m4, $D6m1m3, $D6m1m2, $D6m1m1, $D6m1m0, $D6m1p1, $D6m1p2, $D6m1p3, $D6m1p4, $D6m1p5, $D6m1p6],
     [$D60p6 ,-$D60p5,  $D60p4, -$D60p3,  $D60p2, -$D60p1,  $D60p0,  $D60p1 , $D60p2,  $D60p3,  $D60p4,  $D60p5,  $D60p6 ],
     [$D6p1m6, $D6p1m5, $D6p1m4, $D6p1m3, $D6p1m2, $D6p1m1, $D6p1p0, $D6p1p1, $D6p1p2, $D6p1p3, $D6p1p4, $D6p1p5, $D6p1p6],
     [$D6p2m6, $D6p2m5, $D6p2m4, $D6p2m3, $D6p2m2, $D6p2m1, $D6p2p0, $D6p2p1, $D6p2p2, $D6p2p3, $D6p2p4, $D6p2p5, $D6p2p6],
     [$D6p3m6, $D6p3m5, $D6p3m4, $D6p3m3, $D6p3m2, $D6p3m1, $D6p3p0, $D6p3p1, $D6p3p2, $D6p3p3, $D6p3p4, $D6p3p5, $D6p3p6],
     [$D6p4m6, $D6p4m5, $D6p4m4, $D6p4m3, $D6p4m2, $D6p4m1, $D6p4p0, $D6p4p1, $D6p4p2, $D6p4p3, $D6p4p4, $D6p4p5, $D6p4p6],
     [$D6p5m6, $D6p5m5, $D6p5m4, $D6p5m3, $D6p5m2, $D6p5m1, $D6p5p0, $D6p5p1, $D6p5p2, $D6p5p3, $D6p5p4, $D6p5p5, $D6p5p6],
     [$D6p6m6, $D6p6m5, $D6p6m4, $D6p6m3, $D6p6m2, $D6p6m1, $D6p6p0, $D6p6p1, $D6p6p2, $D6p6p3, $D6p6p4, $D6p6p5, $D6p6p6]
          ];
return $D;
}
# ------------------------------------ Finished !! ------------------------------- #

# Transformation matrix between Stevens and Buckmaster operator equivalents
$j=cplx(0,1);
$A1 = [
         [ $j/sqrt(2),	0, 	$j/sqrt(2) ],
         [ 0,		1,		0  ],
         [ 1/sqrt(2),	0,	-1/sqrt(2) ]
      ];
$invA1 = [
         [ -$j/sqrt(2),	0, 	 1/sqrt(2) ],
         [ 0,		1,		0  ],
         [ -$j/sqrt(2),	0,	-1/sqrt(2) ]
      ];
#$invA1 = $A1; 

$A2 = [
	[$j,	0,	0,	   0,	  -$j],
	[0,	$j/2,	0,	  $j/2,    0],
	[0,	0,	sqrt(6),  0,	   0],
	[0,	1/2,	0,	  -1/2,	   0],
	[1,	0,	0,	   0,	   1]
      ];
$invA2 = [
	[     -$j/2,         0,         0,         0,       1/2],
	[         0,       -$j,         0,         1,         0],
	[         0,         0, sqrt(1/6),        0,         0],
	[         0,       -$j,         0,        -1,         0],
	[      $j/2,         0,         0,         0,       1/2]
      ];

$A3 = [
        [$j*sqrt(2),0,   0,     0,     0,         0, $j*sqrt(2)],
	[0, $j/sqrt(3),  0,     0,     0,   -$j/sqrt(3),  0],
	[0,   0, $j*sqrt(10/3), 0, $j*sqrt(10/3), 0,      0],
	[0,   0,         0, sqrt(10), 0,          0,      0],
	[0,   0,    sqrt(10/3), 0,   -sqrt(10/3), 0,      0],
	[0,  1/sqrt(3),  0,     0,     0,     1/sqrt(3),  0],
        [sqrt(2),0,      0,     0,     0,         0,   -sqrt(2)]
      ];
$invA3 = [
        [-$j/sqrt(8), 0,  0,   0,     0,       0, 1/sqrt(8)],
	[0, -$j*sqrt(3/4),0,   0,     0,   sqrt(3/4),  0],
	[0,   0,-$j*sqrt(3/40),0, sqrt(3/40),  0,      0],
	[0,   0,      0,   1/sqrt(10),0,       0,      0],
	[0,   0,-$j*sqrt(3/40),0,-sqrt(3/40),  0,      0],
	[0, $j*sqrt(3/4), 0,   0,     0,   sqrt(3/4),  0],
        [-$j/sqrt(8), 0,  0,   0,     0,       0,-1/sqrt(8)]
      ];

$A4 = [
	[2*$j,0,        0,      0,        0,      0,         0,     0,	 -2*$j],
	[0,$j/sqrt(2),  0,      0,        0,      0,         0,	$j/sqrt(2),  0],
	[0,   0, $j*sqrt(7),    0,        0,      0, -$j*sqrt(7),   0,	     0],
	[0,   0,        0,$j*sqrt(7/2),   0,$j*sqrt(7/2),    0,     0,	     0],
	[0,   0,        0,      0,    2*sqrt(70), 0,         0,     0,	     0],
	[0,   0,        0,   sqrt(7/2),   0,  -sqrt(7/2),    0,     0,	     0],
	[0,   0,    sqrt(7),    0,	  0,      0,     sqrt(7),   0,	     0],
	[0, 1/sqrt(2),  0,      0,	  0,      0,         0,	-1/sqrt(2),  0],
	[2,   0,        0,      0,	  0,      0,         0,     0,	     2]
     ];
$invA4 = [
	[-$j/4,    0,       0,       0,        0,      0,         0,     0,      1/4],
	[   0,-$j/sqrt(2),  0,       0,        0,      0,         0,  1/sqrt(2),  0],
	[   0,     0, -$j/sqrt(28),  0,        0,      0,  0.5/sqrt(7),  0,       0],
	[   0,     0,       0,-$j/sqrt(2*7),   0,  1/sqrt(2*7),   0,     0,       0],
	[   0,     0,       0,       0,  0.5/sqrt(70), 0,         0,     0,       0],
	[   0,     0,       0,-$j/sqrt(2*7),   0, -1/sqrt(2*7),   0,     0,       0],
	[   0,     0,  $j/sqrt(28),  0,        0,      0,  0.5/sqrt(7),  0,       0],
	[   0,-$j/sqrt(2),  0,       0,        0,      0,         0, -1/sqrt(2),  0],
	[ $j/4,    0,       0,       0,        0,      0,         0,     0,      1/4]
      ];

$A5 = [
        [$j*sqrt(8), 0, 0,      0,      0,      0,      0,      0,      0,      0,$j*sqrt(8)],
        [0,2*$j/sqrt(5),0,      0,      0,      0,      0,      0,      0,-2*$j/sqrt(5),0],
        [0,     0,$j*sqrt(72/5),0,      0,      0,      0,      0,$j*sqrt(72/5),0,      0],
        [0,     0,      0,$j*sqrt(3/5), 0,      0,      0,-$j*sqrt(3/5),0,      0,      0],
        [0,     0,      0,      0,$j*sqrt(84/5),0,$j*sqrt(84/5),0,      0,      0,      0],
        [0,     0,      0,      0,      0,  6*sqrt(14), 0,      0,      0,      0,      0],
        [0,     0,      0,      0,   sqrt(84/5),0,  -sqrt(84/5),0,      0,      0,      0],
        [0,     0,      0,   sqrt(3/5), 0,      0,      0,    sqrt(3/5),0,      0,      0],
        [0,     0,   sqrt(72/5),0,      0,      0,      0,      0,  -sqrt(72/5),0,      0],
        [0,  2/sqrt(5), 0,      0,      0,      0,      0,      0,      0,  2/sqrt(5),  0],
        [sqrt(8), 0,    0,      0,      0,      0,      0,      0,      0,      0,-sqrt(8)],
      ];
$invA5 = [
        [-$j/sqrt(32),0,  0,      0,        0,      0,      0,     0,      0,     0, 1/sqrt(32)],
        [0,-$j*sqrt(5/16),0,      0,        0,      0,      0,     0,      0, sqrt(5/16), 0],
        [0,     0,-$j*sqrt(5/288),0,        0,      0,      0,     0, sqrt(5/288),0,      0],
        [0,     0,        0,-$j*sqrt(15/36),0,      0,      0, sqrt(15/36),0,     0,      0],
        [0,     0,        0,      0,-$j*sqrt(5/336),0, sqrt(5/336),0,      0,     0,      0],
        [0,     0,        0,      0,        0, 1/sqrt(504), 0,     0,      0,     0,      0],
        [0,     0,        0,      0,-$j*sqrt(5/336),0,-sqrt(5/336),0,      0,     0,      0],
        [0,     0,        0, $j*sqrt(15/36),0,      0,      0, sqrt(15/36),0,     0,      0],
        [0,     0,-$j*sqrt(5/288),0,        0,      0,      0,     0,-sqrt(5/288),0,      0],
        [0, $j*sqrt(5/16),0,      0,        0,      0,      0,     0,      0, sqrt(5/16), 0],
        [-$j/sqrt(32), 0, 0,      0,        0,      0,      0,     0,      0,     0,-1/sqrt(32)],
      ];

$A6 = [
	[4*$j,  0,      0,        0,      0,       0,    0,      0,       0,       0,     0,       0, -4*$j],
	[0,2*$j/sqrt(3),0,        0,      0,       0,    0,      0,       0,       0,     0, 2*$j/sqrt(3),0],
	[0,     0,4*$j*sqrt(11/6),0,      0,       0,    0,      0,       0,       0,-4*$j*sqrt(11/6),0,  0],
	[0,     0,      0,2*$j*sqrt(11/5),0,       0,    0,      0,       0,$j*sqrt(44/5),0,       0,     0],
	[0,     0,      0,        0,$j*sqrt(176/5),0,    0,      0,-$j*sqrt(176/5),0,     0,       0,     0],
	[0,     0,      0,        0,      0,$j*sqrt(22), 0,$j*sqrt(2*11), 0,       0,     0,       0,     0],
        [0,     0,      0,        0,      0,       0,4*sqrt(231),0,       0,       0,     0,       0,     0],
	[0,     0,      0,        0,      0,   sqrt(22), 0,  -sqrt(2*11), 0,       0,     0,       0,     0],
	[0,     0,      0,        0,  sqrt(176/5), 0,    0,      0, sqrt(176/5),   0,     0,       0,     0],
	[0,     0,      0,  2*sqrt(11/5), 0,       0,    0,      0,       0,-2*sqrt(11/5),0,       0,     0],
	[0,     0,  4*sqrt(11/6), 0,      0,       0,    0,      0,       0,       0, 4*sqrt(11/6),0,     0],
	[0,  2/sqrt(3), 0,        0,      0,       0,    0,      0,       0,       0,     0,-2/sqrt(3),   0],
	[4,     0,      0,        0,      0,       0,    0,      0,       0,       0,     0,       0,     4]
      ];
$invA6 = [
	[ -$j/8,                0,0,0,0,0,0,0,0,0,0,0,      1/8           ],
	[0, -$j*sqrt(3)/4,        0,0,0,0,0,0,0,0,0,    sqrt(3)/4,      0],
	[0,0, -$j*sqrt(3/2/11)/4,   0,0,0,0,0,0,0,   sqrt(3/2/11)/4, 0,0],
	[0,0,0, -$j*sqrt(5/11)/4,     0,0,0,0,0,   sqrt(5/11)/4,    0,0,0],
	[0,0,0,0, -$j*sqrt(5/11)/8,     0,0,0,    sqrt(5/11)/8,    0,0,0,0],
	[0,0,0,0,0, -$j*sqrt(1/2/11)/2,   0,    sqrt(1/2/11)/2,  0,0,0,0,0],
	[0,0,0,0,0,0,         sqrt(1/3/7/11)/4,             0,0,0,0,0,0],
	[0,0,0,0,0, -$j*sqrt(1/2/11)/2,  0,   -sqrt(1/2/11)/2,  0,0,0,0,0],
	[0,0,0,0,  $j*sqrt(5/11)/8,    0,0,0,   sqrt(5/11)/8,     0,0,0,0],
	[0,0,0, -$j*sqrt(5/11)/4,    0,0,0,0,0,  -sqrt(5/11)/4,     0,0,0],
	[0,0,  $j*sqrt(3/2/11)/4,  0,0,0,0,0,0,0,   sqrt(3/2/11)/4,  0,0],
	[0, -$j*sqrt(3)/4,       0,0,0,0,0,0,0,0,0,  -sqrt(3)/4,        0],
	[   $j/8,              0,0,0,0,0,0,0,0,0,0,0,       1/8            ]
      ];

@A = ($A1, $A2, $A3, $A4, $A5, $A6);
@invA = ($invA1, $invA2, $invA3, $invA4, $invA5, $invA6);

sub rotk {
  my ($k,$t,$f) = @_;
  my $D;
  if (!$f) { $f = 0; }
  if (!$t) { $t = 0; }
  unless ($k>0 && $k<7) {
    die "Sorry only operators up to rank 6 supported\n"; }

     if($k==1) { $D = rot1($t,$f); }
  elsif($k==2) { $D = rot2($t,$f); }
  elsif($k==3) { $D = rot3($t,$f); }
  elsif($k==4) { $D = rot4($t,$f); }
  elsif($k==5) { $D = rot5($t,$f); }
  elsif($k==6) { $D = rot6($t,$f); }

# for $ii (0..(2*$J)) { print "["; for $jj (0..(2*$J)) {
#    print $D->[$ii][$jj],"\t"; } print "]\n"; }
}

# The operator in the rotated frame. Rx(y) is equivalent to {Oxy} in Rudowicz notation
# and Ox_y, Oxmy is equivalent to [Oxy], [OxyM] in Rudowicz

#my $B2 = [ [$O2m2],[$O2m1],[$O2p0],[$O2p1],[$O2p2] ];

#my $tmp1 = mmult($invA2,$B2);
#my $tmp2 = mmult($D2, $tmp1);
#my $R2 = mmult($A2, $tmp2);

#my $S2 = [];
#for $l (0 .. 4) { 
#  $S2->[$l] = $R2->[$l][0]->sub(f=>$phi_val, t=>$theta_val);
#  $S2->[$l] = $S2->[$l]->sub(six=>6); 
#  $M2->[0][$l] = $S2->[$l]->sub(Otmt=>1,Otmo=>0,Otpz=>0,Otpo=>0,Otpt=>0);
#  $M2->[1][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>1,Otpz=>0,Otpo=>0,Otpt=>0);
#  $M2->[2][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>0,Otpz=>1,Otpo=>0,Otpt=>0);
#  $M2->[3][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>0,Otpz=>0,Otpo=>1,Otpt=>0);
#  $M2->[4][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>0,Otpz=>0,Otpo=>0,Otpt=>1);
#}

#print "Test\n";
