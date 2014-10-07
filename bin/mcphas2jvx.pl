#!/usr/bin/perl
#
# mcphas2jvx.pl - script to convert mcphas.j and other input files into a Javaview .jvx for viewing
#
#

use PDL;
use PDL::Slatec;
use Getopt::Long;

# Loads the tables of element information (ionic size and colours)
push @INC, $ENV{'MCPHASE_DIR'}.'/bin/';
require 'elements.pl';

# Set to 1 to print extra information, and to not create any files, but to pipe everything to STDOUT
$debug = 0;

@parn = ("a", "b", "c",
         "alpha", "beta", "gamma",
         "nofatoms", "nofcomponents");
@atpn = ("da", "db", "dc", 
         "nofneighbours",
         "diagonalexchange", 
         "sipffilename");
$outputfile = "results/mcphas.jvx";

# Parses command line options
GetOptions("help"=>\$helpflag,
           "debug"=>\$debug,
           "output"=>\$outputfile);

if ($#ARGV<0 || $helpflag) {
   print " $0:\n";
   print "   - script to convert mcphas.j and other input files into a Javaview .jvx for viewing\n\n";
   print " Syntax: $0 [mcphas.j]\n\n";
   print " where [mcphas.j] is any of the mcphas.j, mcphas_magnetic_atoms.j or mcphas_all_atoms.j input files.\n\n";
   print " Options include:\n";
   print "    --help          - prints this message\n";
   print "    --output [FILE] - outputs JVX data to [FILE] instead of default results/mcphas.jvx\n\n";
   exit(0);
}

$jfile = $ARGV[0];

# ------------------------------------------------------------------------------------------------------------------------ #
# Reads in the mcphas.j file, and determines the lattice, atomic positions and exchanges. 
# ------------------------------------------------------------------------------------------------------------------------ #
$maxsumj=0; $atidx=0;
while(<>) {
  $_ =~ s/\R//g;            # safe chomp
  if($_ =~ /^#!/) {         # extracts the structure/lattice parameters
    foreach $id (@parn) {
      if($_ =~ /\s+$id\s*=/) {
        ($pars{$id}) = ($_ =~ m/\s+$id\s*=\s*([\d.eEdD\-\+]+)/);
      }
    }
    foreach $id (@atpn) {
      if($_ =~ /\s+$id\s*=/) {
        push @{$atoms{$id}}, ($_ =~ m/\s+$id\s*=\s*([\d.eEdD\-\+\.\w]+)/);
      }
    }
  } 
  if($_ =~ /^#/) { 
    if(@{$neighbours[$atidx]} != "") { $atidx++; }
    @{$neighbours[$atidx]}=();
  } else {
    @ll = split; $sumj=0; for $ii(3..$#ll) { $sumj+=abs($ll[$ii]); }
    if(abs($sumj)>1e-6) { 
      push @{$neighbours[$atidx]}, $_." ".(sprintf "%5.2f",$sumj); 
      if($sumj>$maxsumj) { $maxsumj = $sumj; }
    }
  }
}

# Loops over the list of nearest neighbours and determines the number of unique ones (from sumj)
for $ii (0..$#{$atoms{"da"}}) {
  %uniqtmp = ();
  for $jj (0..$#{$neighbours[$ii]}) {
    @pos = split(" ",${$neighbours[$ii]}[$jj]);
    $uniqtmp{$pos[-1]} = 1;
  }
  for (keys %uniqtmp) { push @{$uniqneig[$ii]}, $_; } 
  @{$uniqneig[$ii]} = sort {$b cmp $a} @{$uniqneig[$ii]};
}

if($debug) {
  for (keys %pars) {
    print STDERR "$_:$pars{$_}\n";
  }

  print STDERR "nofatoms=".$pars{"nofatoms"}."\tlastIndexRead=".$#{$atoms{"da"}}."\n";

  for $ii (0..$#{$atoms{"da"}}) {
    for (@atpn) {
      print STDERR "$_=${$atoms{$_}}[$ii]\t";
    }
    print STDERR "\n";
    print STDERR "lastneighbour_index=".$#{$neighbours[$ii]}."\n";
    for $jj (0..$#{$neighbours[$ii]}) {
      print STDERR ${$neighbours[$ii]}[$jj]."\n";
    }
    print STDERR "unique exchanges=";
    foreach (@{$uniqneig[$ii]}) { print STDERR "$_ "; }
    print STDERR "\n";
  }
}

# Determine the transformation matrix to convert from fractional to Cartesian coordinates and back 
#   (using McPhase convention, y||b, z_|_ab, a perp to both y,z)
$PI=3.141592654;
$a = $pars{"a"};
$b = $pars{"b"};
$c = $pars{"c"};
$alpha = $pars{"alpha"}*$PI/180;
$beta  = $pars{"beta"}*$PI/180;
$gamma = $pars{"gamma"}*$PI/180;
$sc = sin($gamma); $ca = cos($alpha); $cb = cos($beta); $cc = cos($gamma);
$v = ($cb-$ca*$cc); $r13 = $c*$v/$sc; $r23 = $c*$ca; $r33 = $c*$c-$r13*$r13-$r23*$r23;
if(abs($r13)>$c || $r33<=0) { die "Error: alpha, beta, gamma are geometrically inconsistent\n"; }
$rtoijk = pdl [ [ $a*$sc, $a*$cc, 0 ],
                [      0,     $b, 0 ],
                [   $r13,   $r23, sqrt($r33) ] ];
$invrtoijk = matinv($rtoijk);

if($debug==1) {
  print STDERR $rtoijk;
  print STDERR $invrtoijk;
}

# Loops through the list of atoms, and look into the sipf to determine element/ion name
for $ii (0..$#{$atoms{"da"}}) {
  $sipf = ${$atoms{"sipffilename"}}[$ii];
  if(!(defined $sipfseen{$sipf})) {
    $sipfseen{$sipf} = 1;
    open(FSIPF,$sipf);
    while(<FSIPF>) {
      $_ =~ s/\R//g;            # safe chomp
      if($_ =~ /#!MODULE/) {
        ($module{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w]+)/); }
      if($_ =~ /IONTYPE/) {
        ($iontyp{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w\=]+)/); }
    }
    close FSIPF;
    $attab = $element{$iontyp{$sipf}};
    if(!$attab) {
      $at = $sipf; $at =~ s/\.sipf//;
      $at =~ s/_//g; $at =~ s/[0-9]+[A-Za-z\-\+]+//g; $at =~ s/[0-9]+//g; 
      $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g;
      $at = lc $at; $at = ucfirst $at;
      $attab = $element{$at};
      if(!$attab) {
        print STDERR "Warning: unknown element: ".$iontyp{$sipf}." - assuming it's carbon!\n";
        $attab = $element{"C"};
      }
    }
    $atcol{$sipf} = ${$attab}[5];
    $r_ion{$sipf} = ${$attab}[4];
    if($debug) { print STDERR "|$sipf|:\t$module{$sipf}\t$iontyp{$sipf}\t$r_ion{$sipf}\t@{$atcol{$sipf}}\n"; }
  }
}

# ------------------------------------------------------------------------------------------------------------------------ #
# Outputs the JVX file for JavaView
# ------------------------------------------------------------------------------------------------------------------------ #
if($debug) { *FOUT = *STDOUT; } else { open (FOUT, ">$outputfile"); }

# Draws the crystallographic unit cell
print FOUT "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n";
print FOUT "<!DOCTYPE jvx-model SYSTEM \"http://www.javaview.de/rsrc/jvx.dtd\">\n";
print FOUT "<jvx-model>\n";
print FOUT "  <title>$title</title>\n";
print FOUT "  <geometries>\n";
print FOUT "    <geometry name=\"crystallographic unit cell\">\n";
print FOUT "      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n";
print FOUT "        <points>\n";
printf FOUT "          <p>          % 10.5f% 10.5f% 10.5f </p>\n", 0,0,0;
printf FOUT "          <p name=\"a\"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,0),$rtoijk->at(1,0),$rtoijk->at(2,0);
printf FOUT "          <p name=\" \"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,0)+$rtoijk->at(0,1),$rtoijk->at(1,0)+$rtoijk->at(1,1),$rtoijk->at(2,0)+$rtoijk->at(2,1);
printf FOUT "          <p name=\"b\"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,1),$rtoijk->at(1,1),$rtoijk->at(2,1);
printf FOUT "          <p name=\"c\"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,2),$rtoijk->at(1,2),$rtoijk->at(2,2);
printf FOUT "          <p name=\" \"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,0)+$rtoijk->at(0,2),$rtoijk->at(1,0)+$rtoijk->at(1,2),$rtoijk->at(2,0)+$rtoijk->at(2,2);
printf FOUT "          <p name=\" \"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,0)+$rtoijk->at(0,1)+$rtoijk->at(0,2),$rtoijk->at(1,0)+$rtoijk->at(1,1)+$rtoijk->at(1,2),$rtoijk->at(2,0)+$rtoijk->at(2,1)+$rtoijk->at(2,2);
printf FOUT "          <p name=\" \"> % 10.5f% 10.5f% 10.5f </p>\n",$rtoijk->at(0,1)+$rtoijk->at(0,2),$rtoijk->at(1,1)+$rtoijk->at(1,2),$rtoijk->at(2,1)+$rtoijk->at(2,2);
print FOUT "          <thickness>0.0</thickness>\n";
print FOUT "          <colorTag type=\"rgb\">255 0 0</colorTag>\n";
print FOUT "                  <labelAtt horAlign=\"head\" visible=\"show\" font=\"fixed\" verAlign=\"top\">\n";
print FOUT "                          <xOffset>0</xOffset>\n";
print FOUT "                          <yOffset>0</yOffset>\n";
print FOUT "                  </labelAtt>\n";
print FOUT "        </points>\n";
print FOUT "      </pointSet>\n";
print FOUT "      <lineSet  arrow=\"hide\" line=\"show\">\n";
print FOUT "        <lines>\n";
print FOUT "          <l>0 1</l>\n";
print FOUT "          <l>1 2</l>\n";
print FOUT "          <l>2 3</l>\n";
print FOUT "          <l>3 0</l>\n";
print FOUT "          <l>0 4</l>\n";
print FOUT "          <l>1 5</l>\n";
print FOUT "          <l>2 6</l>\n";
print FOUT "          <l>3 7</l>\n";
print FOUT "          <l>4 5</l>\n";
print FOUT "          <l>5 6</l>\n";
print FOUT "          <l>6 7</l>\n";
print FOUT "          <l>7 4</l>\n";
print FOUT "          <thickness>1.0</thickness>\n";
printf FOUT "        <color type=\"rgb\">%3d %3d %3d</color>\n",255,0,0;
print FOUT "        </lines>\n";
print FOUT "      </lineSet>\n";
print FOUT "    </geometry>\n";

# Draws each type of *sipf as an individual geometry - so user and make each visible or not separately
for $sipf (keys %sipfseen) {
  @colours = ();
  printf FOUT "    <geometry name=\"%s\">\n", $sipf;
  print FOUT "      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n";
  print FOUT "        <points>\n";
  for $i1(0..1) { for $i2(0..1) { for $i3(0..1) { 
    for $ii (0..$#{$atoms{"da"}}) {
      if(!($sipf eq ${$atoms{"sipffilename"}}[$ii])) { next; }
      $p1 = ${$atoms{"da"}}[$ii]+$i1; $p2 = ${$atoms{"db"}}[$ii]+$i2; $p3 = ${$atoms{"dc"}}[$ii]+$i3;
      if($p1>=0 && $p1<=1 && $p2>=0 && $p2<=1 && $p3>=0 && $p3<=1) {
        $fpos = pdl [ $p1, $p2, $p3 ];
        $cpos = $fpos x $rtoijk;
        printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0);
        push @colours, join(":",@{$atcol{${$atoms{"sipffilename"}}[$ii]}});
      }
    }
  }}}
  printf FOUT "          <thickness>% 10.5f</thickness>\n",$r_ion{$sipf}*10;
  print FOUT "        </points>\n";
  print FOUT "        <colors type=\"rgb\">\n";
  for (@colours) {
    @cl = split(":");
    printf FOUT "          <c> %3d %3d %3d </c>\n", @cl;
  }
  print FOUT "        </colors>\n";
  print FOUT "      </pointSet>\n";
  print FOUT "    </geometry>\n";
}

# Draws the exchange interactions for each atom as lines
for $ii (0..$#{$atoms{"da"}}) {
  if($#{$neighbours[$ii]} == -1) { next; }
  for $kk (0..$#{$uniqneig[$ii]}) {
    $nneigh = 0;
    $fpos = pdl [ ${$atoms{"da"}}[$ii], ${$atoms{"db"}}[$ii], ${$atoms{"dc"}}[$ii] ];
    $cpos = $fpos x $rtoijk;
    printf FOUT "    <geometry name=\"atom %d (type %s) exchange set %d\">\n", $ii+1, ${$atoms{"sipffilename"}}[$ii], $kk+1;
    print FOUT "      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n";
    print FOUT "        <points>\n";
    printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n", $cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0);
    for $jj (0..$#{$neighbours[$ii]}) {
      @pos = split(" ",${$neighbours[$ii]}[$jj]);
      if(abs($pos[-1]-${$uniqneig[$ii]}[$kk])>0.01) { next; }
      $fpos = pdl [ ${$atoms{"da"}}[$ii]+$pos[0], ${$atoms{"db"}}[$ii]+$pos[1], ${$atoms{"dc"}}[$ii]+$pos[2] ];
      $cpos = $fpos x $rtoijk;
      printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n", $cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0);
      $nneigh++;
    }
    print FOUT "          <thickness>0.0</thickness>\n";
    print FOUT "        </points>\n";
    print FOUT "      </pointSet>\n";
    print FOUT "      <lineSet  arrow=\"hide\" line=\"show\">\n";
    print FOUT "        <lines>\n";
    for $jj (0..$nneigh-1) {
      printf FOUT "          <l>0 %d</l>\n",$jj+1;
    }
    printf FOUT "          <thickness>%5.3f</thickness>\n",5*(log(${$uniqneig[$ii]}[$kk]/$maxsumj*9+1)/log(10));
    printf FOUT "        <color type=\"rgb\">%3d %3d %3d</color>\n",@{$atcol{${$atoms{"sipffilename"}}[$ii]}};
    print FOUT "        </lines>\n";
    print FOUT "      </lineSet>\n";
    print FOUT "    </geometry>\n";
  }
}

print FOUT "  </geometries>\n";
print FOUT "</jvx-model>\n";

if($debug==0) { 
  close FOUT; 
  print "$0 completed. You can view the output by typing: javaview $outputfile\n";
}
