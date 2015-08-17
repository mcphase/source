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
         "sipffilename",
         "cffilename");
@rmat = ("r1a","r2a","r3a",
         "r1b","r2b","r3b",
         "r1c","r2c","r3c");
$outputfile = "results/mcphas.jvx";

# Parses command line options
GetOptions("help"=>\$helpflag,
           "debug"=>\$debug,
           "individual"=>\$clusterflag,
           "primitive"=>\$primflag,
           "output=s"=>\$outputfile);

if ($#ARGV<0 || $helpflag) {
   print " $0:\n";
   print "   - script to convert mcphas.j and other input files into a Javaview .jvx for viewing\n\n";
   print " Syntax: $0 [mcphas.j]\n\n";
   print " where [mcphas.j] is any of the mcphas.j, mcphas_magnetic_atoms.j or mcphas_all_atoms.j input files.\n\n";
   print " Options include:\n";
   print "    --help          - prints this message\n";
   print "    --output [FILE] - outputs JVX data to [FILE] instead of default results/mcphas.jvx\n\n";
   print "    --individual    - for input files for the cluster module, plot individual atomic positions\n";
   print "                      and intra-cluster exchange as well as cluster centres and mean-field exchange.\n\n";
   print "    --primitive     - plots only atoms explicitly specified in the .j file (without adding atoms at\n";
   print "                      the edge of the supercell\n";
   print " Note that UNIX switch conventions apply so you can use \"short\" switches (-h, -o -i) as well.\n";
   print " E.g. to plot individual atoms and output to a different file:\n";
   print "   $0 mcphas.j -i -o results/cluster.jvx\n";
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
    foreach $id (@rmat) {
      if($_ =~ /\s+$id\s*=/) {
        push @{$rmath{$id}}, ($_ =~ m/\s+$id\s*=\s*([\d.eEdD\-\+\.\w]+)/);
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
$rmat = pdl [ [ $rmath{"r1a"}[0], $rmath{"r2a"}[0], $rmath{"r3a"}[0] ],
              [ $rmath{"r1b"}[0], $rmath{"r2b"}[0], $rmath{"r3b"}[0] ],
              [ $rmath{"r1c"}[0], $rmath{"r2c"}[0], $rmath{"r3c"}[0] ] ];
$primcell = $rtoijk x $rmat;
$invprim  = matinv($primcell);

if($debug==1) {
  print STDERR $rtoijk;
  print STDERR $invrtoijk;
  print STDERR $rmat;
  print STDERR $primcell;
}

# calc_minmax_scale_relabc() function from spincf_out.cpp
@abc=($a,$b,$c);
for $ii (0..2) {
  $ddd = pdl [ $rmat->at(0,ii), $rmat->at(1,ii), $rmat->at(2,ii), 
           $rmat->at(0,ii)+$rmat->at(1,ii),  
           $rmat->at(0,ii)+$rmat->at(2,ii),  
           $rmat->at(1,ii)+$rmat->at(2,ii), 0, sum($rmat->slice(":,0")) ]; 
  $t = min($ddd)/$abc[$ii]; if (abs($t-int($t))>0.0001) { push @minv, (int($t)-1.)*$abc[$ii]; } else { push @minv, min($ddd); }
  $t = max($ddd)/$abc[$ii]; if (abs($t-int($t))>0.0001) { push @maxv, (int($t)+1.)*$abc[$ii]; } else { push @maxv, max($ddd); }
}
for $ii (0..2) {
  @ddd = ();
  $dd0 = pdl @minv;                        $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @minv; $dd0->set(0,$maxv[0]); $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @minv; $dd0->set(0,$maxv[2]); $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @minv; $dd0->set(0,$maxv[3]); $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @maxv;                        $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @maxv; $dd0->set(0,$minv[0]); $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @maxv; $dd0->set(0,$minv[2]); $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  $dd0 = pdl @maxv; $dd0->set(0,$minv[3]); $dd = $invprim x transpose($dd0); push @ddd, $dd->at(0,$ii);
  push @ijkmin, min($ddd);
  push @ijkmax, max($ddd);
}
for $ii (0..2) { $minv[$ii]/=$abc[$ii]; $maxv[$ii]/=$abc[$ii]; }

# Loops through the list of atoms, and look into the sipf to determine element/ion name
for $ii (0..$#{$atoms{"da"}}) {
  $sipf = ${$atoms{"sipffilename"}}[$ii];
  if(!$sipf) { $sipf = ${$atoms{"cffilename"}}[$ii]; }
  if(!(defined $sipfseen{$sipf})) {
    $sipfseen{$sipf} = 1;
    open(FSIPF,$sipf);
    while(<FSIPF>) {
      $_ =~ s/\R//g;            # safe chomp
      if($_ =~ /#!MODULE/) {
        ($module{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w]+)/); }
      if($_ =~ /IONTYPE/) {
        ($iontyp{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w\=]+)/); }
      if($_ =~ /structurefile/) {
        ($cluste{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w\=]+)/); }
    }
    close FSIPF;
    if($module{$sipf} =~ /cluster/ && $clusterflag) {
      $atidx=0;
      open(FCLUS,$cluste{$sipf}); # looks inside cluster.j file for relative coordinates of cluster atoms and exchanges
      while(<FCLUS>) {
        $_ =~ s/\R//g;            # safe chomp
        if($_ =~ /^#!/) {         # extracts the structure/lattice parameters
          foreach $id (@atpn) {
            if($_ =~ /\s+$id\s*=/) {
              if($id eq "da") { $atidx++; @{$cluster_neighbours{$ii."|".$atidx}}=(); }
              ($cluster_atoms{$ii."|".$atidx."|".$id}) = ($_ =~ m/\s+$id\s*=\s*([\d.eEdD\-\+\.\w]+)/);
            }
          }
        }
        if($_ !~ /^#/) {
          @ll = split; $sumj=0; for $ii(3..$#ll) { $sumj+=abs($ll[$ii]); }
          if(abs($sumj)>1e-6) {
            push @{$cluster_neighbours{$ii."|".$atidx}}, $_." ".(sprintf "%5.2f",$sumj);
            if($sumj>$maxsumj) { $maxsumj = $sumj; }
          }
        }
        $cluster_natom[$ii] = $atidx;
      }
      close FCLUS;
      if($debug) {
        print STDERR "Cluster centre:\t\t\t";
        for (@atpn) { print STDERR "$_=${$atoms{$_}}[$ii]\t"; }
        print STDERR "\n";
        for $atidx (1..$cluster_natom[$ii]) {
          print STDERR "Atom $atidx, relative coordinates:\t";
          for (@atpn) { print STDERR "$_=".$cluster_atoms{$ii."|".$atidx."|".$_}."\t"; }
          print STDERR "\n";
        }
      }
      # Opens the real single ion parameter files for the atoms within the cluster
      for $atidx (1..$cluster_natom[$ii]) {
        $sipf = $cluster_atoms{$ii."|".$atidx."|sipffilename"};
        if(!(defined $cluster_sipfseen{$sipf})) {
          $cluster_sipfseen{$sipf} = 1;
          open(FSIPF,$sipf);
          while(<FSIPF>) {
            $_ =~ s/\R//g;            # safe chomp
            if($_ =~ /#!MODULE/) {
              ($module{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w]+)/); }
            if($_ =~ /IONTYPE/) {
              ($iontyp{$sipf}) = ($_ =~ m/\s*=\s*([\d.eEdD\-\+\.\w\=]+)/); }
          }
          close FSIPF;
          $at = $iontyp{$sipf}; $at =~ s/\.sipf//;
          $at =~ s/_//g; $at =~ s/[0-9]+[A-Za-z\-\+]+//g; $at =~ s/[0-9]+//g;
          $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\/\|\[\]\{\}\$]//g;
          $at = lc $at; $at = ucfirst $at;
          $attab = $element{$at};
          if(!$attab) {
            $at = $sipf; $at =~ s/\.sipf//;
            $at =~ s/_//g; $at =~ s/[0-9]+[A-Za-z\-\+]+//g; $at =~ s/[0-9]+//g;
            $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\/\|\[\]\{\}\$]//g;
            $at = lc $at; $at = ucfirst $at;
            $attab = $element{$at};
            if(!$attab) {
              print STDERR "Warning: unknown element: ".$iontyp{$sipf}." - assuming it's carbon!\n";
              $attab = $element{"C"};
            }
          }
          if($debug) { print STDERR "at=$at\n"; }
          $atcol{$sipf} = ${$attab}[5];
          $r_ion{$sipf} = ${$attab}[4];
          if($debug) { print STDERR "|$sipf|:\t$module{$sipf}\t$iontyp{$sipf}\t$r_ion{$sipf}\t@{$atcol{$sipf}}\n"; }
        }
      }
    } else {
      $at = $iontyp{$sipf}; $at =~ s/\.sipf//;
      $at =~ s/_//g; $at =~ s/[0-9]+[A-Za-z\-\+]+//g; $at =~ s/[0-9]+//g; 
      $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\/\|\[\]\{\}\$]//g;
      $at = lc $at; $at = ucfirst $at;
      $attab = $element{$at};
      if(!$attab) {
        $at = $sipf; $at =~ s/\.sipf//;
        $at =~ s/_//g; $at =~ s/[0-9]+[A-Za-z\-\+]+//g; $at =~ s/[0-9]+//g; 
        $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\/\|\[\]\{\}\$]//g;
        $at = lc $at; $at = ucfirst $at;
        $attab = $element{$at};
        if(!$attab) {
          print STDERR "Warning: unknown element: ".$iontyp{$sipf}." - assuming it's carbon!\n";
          $attab = $element{"C"};
        }
      }
      if($debug) { print STDERR "at=$at\n"; }
      $atcol{$sipf} = ${$attab}[5];
      $r_ion{$sipf} = ${$attab}[4];
      if($debug) { print STDERR "|$sipf|:\t$module{$sipf}\t$iontyp{$sipf}\t$r_ion{$sipf}\t@{$atcol{$sipf}}\n"; }
    }
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

# Determines the extent of the supercell
if ($primflag) {
  @r1 = (0);
  @r2 = (0);
  @r3 = (0);
} else {
  @r1 = (min($rmat->slice("0,:")), max($rmat->slice("0,:")));
  @r2 = (min($rmat->slice("1,:")), max($rmat->slice("1,:")));
  @r3 = (min($rmat->slice("2,:")), max($rmat->slice("2,:")));
}

# Draws each type of *sipf as an individual geometry - so user and make each visible or not separately
for $sipf (keys %sipfseen) {
  if($module{$sipf} =~ /cluster/ && $clusterflag) {
    # Clusters - put a cube at the centre of the cluster, and plots the cluster atoms as spheres connected by lines to centre.
    for $i1(@r1) { for $i2(@r2) { for $i3(@r3) {
      for $ii (0..$#{$atoms{"da"}}) {
        if(!($sipf eq ${$atoms{"sipffilename"}}[$ii]) && !($sipf eq ${$atoms{"cffilename"}}[$ii])) { next; }
        @colours = (); @thicknesses = ();
        $p1 = ${$atoms{"da"}}[$ii]+$i1; $p2 = ${$atoms{"db"}}[$ii]+$i2; $p3 = ${$atoms{"dc"}}[$ii]+$i3;
        if($p1>=$ijkmin[0] && $p1<=$ijkmax[0] && $p2>=$ijkmin[1] && $p2<=$ijkmax[1] && $p3>=$ijkmin[2] && $p3<=$ijkmax[2]) {
          $fpos = pdl [ $p1, $p2, $p3 ];
          $cpos = $fpos x $rtoijk;
          # Plots a cube at the centre of the cluster
          printf FOUT "    <geometry name=\"cluster type %s id %d\">\n", $sipf, $ii;
          print FOUT "      <pointSet dim=\"3\" point=\"show\" color=\"show\" thicknesses=\"show\">\n";
          print FOUT "        <points>\n";
          # Points are: centre of cluster, 8 points defining cube, then atomic centres and these shifted slightly to get lines
          #  as an ultra-thin face because we're not allowed a pointset+lineset+faceset in single geometry
          printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0);
          printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0)+0.001;
          push @colours, "0:0:0"; push @colours, "0:0:0";
          push @thicknesses, 0;   push @thicknesses, 0;
          for $xx (-0.5,0.5) { for $yy (-0.5,0.5) { for $zz (-0.5,0.5) { 
            printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$cpos->at(0,0)+$xx,$cpos->at(1,0)+$yy,$cpos->at(2,0)+$zz;
            push @colours, "0:0:0";
            push @thicknesses, 0;
          }}}
          for $atidx (1..$cluster_natom[$ii]) {
            $q1 = $cluster_atoms{$ii."|".$atidx."|da"}+$p1;
            $q2 = $cluster_atoms{$ii."|".$atidx."|db"}+$p2;
            $q3 = $cluster_atoms{$ii."|".$atidx."|dc"}+$p3;
            $fpos = pdl [ $q1, $q2, $q3 ];
            $cpos = $fpos x $rtoijk;
            printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0);
            printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$cpos->at(0,0),$cpos->at(1,0),$cpos->at(2,0)+0.001;
            push @colours, join(":",@{$atcol{$cluster_atoms{$ii."|".$atidx."|sipffilename"}}}); push @colours, "0:0:0";
            push @thicknesses, $r_ion{$cluster_atoms{$ii."|".$atidx."|sipffilename"}}*10;       push @thicknesses, 0;
          }
          print FOUT "        </points>\n";
          print FOUT "        <colors type=\"rgb\">\n";
          for (@colours) {
            printf FOUT "          <c> %3d %3d %3d </c>\n", split(":",$_);
          }
          print FOUT "        </colors>\n";
          print FOUT "        <thicknesses>\n";
          for (@thicknesses) {
            printf FOUT "          <th> %f </th>\n", $_;
          }
          print FOUT "        </thicknesses>\n";
          print FOUT "      </pointSet>\n";
          print FOUT "      <faceSet face=\"show\" edge=\"show\">\n";
          print FOUT "        <faces>\n";             # - - - 1    z
          print FOUT "          <f> 2 3 5 4 </f>\n";  # - - + 2   / \ (2,6)----(4,8)   (1,2,4,3)
          print FOUT "          <f> 2 3 7 6 </f>\n";  # - + - 3    |    |        |     (1,2,6,5)
          print FOUT "          <f> 3 5 9 7 </f>\n";  # - + + 4    |    |        |     (2,4,8,6)
          print FOUT "          <f> 5 4 8 9 </f>\n";  # + - - 5         |        |     (4,3,7,8)
          print FOUT "          <f> 4 2 6 8 </f>\n";  # + - + 6       (1,5)----(3,7)   (3,1,5,7)
          print FOUT "          <f> 6 7 9 8 </f>\n";  # + + - 7                        (5,6,8,7)
          for $atidx (1..$cluster_natom[$ii]) {       # + + + 8    ------>  y
            printf FOUT "          <f> 0 %d %d 1 </f>\n", $atidx*2+8, $atidx*2+9;
          }
          print FOUT "          <color type=\"rgb\">0 97 255 </color>\n";
          print FOUT "          <colorTag type=\"rgb\">255 0 255</colorTag>\n";
          print FOUT "        </faces>\n";
          print FOUT "      </faceSet>\n";
          print FOUT "    </geometry>\n";
          # Plots the exchange interactions within the cluster as rectangular "lines"
          printf FOUT "    <geometry name=\"cluster type %s id %d exchanges\">\n", $sipf, $ii;
          print FOUT "      <pointSet dim=\"3\" point=\"hide\">\n";
          print FOUT "        <points>\n";
          $ptid = 0; @exln = (); @exth = (); @colours = (); @exfc = ();                  # (0,1) (2,3)
          @cub = ( [0,2,4,6], [0,2,3,1], [2,3,5,4], [4,5,7,6], [6,7,1,0], [1,3,5,7] );   # (6,7) (4,5)
          for $atidx (1..$cluster_natom[$ii]) {
            $q1 = $cluster_atoms{$ii."|".$atidx."|da"}+$p1;
            $q2 = $cluster_atoms{$ii."|".$atidx."|db"}+$p2;
            $q3 = $cluster_atoms{$ii."|".$atidx."|dc"}+$p3;
            $fpos0 = pdl [ $q1, $q2, $q3 ];
            $cpos0 = $fpos0 x $rtoijk;
            for (@{$cluster_neighbours{$ii."|".$atidx}}) {
              @pos = split;
              $fpos = pdl [ $q1+$pos[0], $q2+$pos[1], $q3+$pos[2] ];
              $cpos = $fpos x $rtoijk;
              $axis = norm($cpos0 - $cpos); $perp = norm( crossp $axis, $fpos0 ); 
              $ux=$axis->at(0,0); $uy=$axis->at(1,0); $uz=$axis->at(2,0); 
              for $th(-$PI,-$PI/2,0,$PI/2) {
                $R = pdl [ [ cos($th)+$ux*$ux*(1-cos($th)), $ux*$uy*(1-cos($th))-$uz*sin($th), $ux*$uz*(1-cos($th))+$uy*sin($th) ],
                           [ $uy*$ux*(1-cos($th))+$uz*sin($th), cos($th)+$uy*$uy*(1-cos($th)), $uy*$uz*(1-cos($th))-$ux*sin($th) ],
                           [ $uz*$ux*(1-cos($th))-$uy*sin($th), $uz*$uy*(1-cos($th))+$ux*sin($th), cos($th)+$uz*$uz*(1-cos($th)) ] ];
                $vpos = ($perp x $R)*(0.05*(log($pos[-1]/$maxsumj*9+1)/log(10))) + $cpos;
                printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$vpos->at(0,0),$vpos->at(1,0),$vpos->at(2,0);
                $vpos = ($perp x $R)*(0.05*(log($pos[-1]/$maxsumj*9+1)/log(10))) + $cpos0;
                printf FOUT "          <p> % 10.5f% 10.5f% 10.5f </p>\n",$vpos->at(0,0),$vpos->at(1,0),$vpos->at(2,0);
              }
              for $nn(0..5) { 
                $lin=""; for $mm(0..3) { $lin.=$cub[$nn][$mm]+$ptid.":"; } push @exfc, $lin; 
                push @colours, join(":",@{$atcol{$cluster_atoms{$ii."|".$atidx."|sipffilename"}}});
              } 
              $ptid += 8;
            }
          }
          print FOUT "        </points>\n";
          print FOUT "      </pointSet>\n";
          print FOUT "      <faceSet face=\"show\" edge=\"show\" color=\"show\">\n";
          print FOUT "        <faces>\n";
          for (@exfc) {
            printf FOUT "          <f> %d %d %d %d </f>\n", split(":"); 
          }
          print FOUT "        </faces>\n";
          print FOUT "        <colors>\n";
          for (@colours) {
            printf FOUT "          <c> %3d %3d %3d </c>\n", split(":");
          }
          print FOUT "        </colors>\n";
          print FOUT "      </faceSet>\n";
          print FOUT "    </geometry>\n";
        }
      }
    }}}
  } else {
    @colours = ();
    printf FOUT "    <geometry name=\"%s\">\n", $sipf;
    print FOUT "      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n";
    print FOUT "        <points>\n";
    for $i1(@r1) { for $i2(@r2) { for $i3(@r3) {
      for $ii (0..$#{$atoms{"da"}}) {
        if(!($sipf eq ${$atoms{"sipffilename"}}[$ii]) && !($sipf eq ${$atoms{"cffilename"}}[$ii])) { next; }
        $p1 = ${$atoms{"da"}}[$ii]+$i1; $p2 = ${$atoms{"db"}}[$ii]+$i2; $p3 = ${$atoms{"dc"}}[$ii]+$i3;
        if($p1>=$ijkmin[0] && $p1<=$ijkmax[0] && $p2>=$ijkmin[1] && $p2<=$ijkmax[1] && $p3>=$ijkmin[2] && $p3<=$ijkmax[2]) {
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
