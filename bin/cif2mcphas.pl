#!/usr/bin/perl
#
# cif2mcphas.pl
#
# Perl script to convert a CIF into a mcphas.j input file, similarly to powdercell2j.pl
# Invoke as:
#
# cif2mcphas.pl [CIFFILE]
#
# This program is part of the McPhase package, licensed under the GNU GPL v2. Please see the COPYING file
#
# Fri Sep 26 12:00:53 KST 2014 - Duc Le - mducle@snu.ac.kr

use PDL;
#use PDL::Slatec;
use File::Copy;
use Getopt::Long;
use IPC::Open3;

# Loads the tables of symmetry equivalent positions, and atomic/ionic information
push @INC, $ENV{'MCPHASE_DIR'}.'/bin/';
require 'itc_syms.pl';
require 'elements.pl';

# Set to 1 to print extra information, and to not create any files, but to pipe everything to STDOUT
$debug = 0;

@cpname = ("_cell_length_a",
           "_cell_length_b",
           "_cell_length_c",
           "_cell_angle_alpha",
           "_cell_angle_beta",
           "_cell_angle_gamma");
@datcol = ("_atom_site_label",
           "_atom_type_oxidation_number",
           "_atom_site_fract_x",
           "_atom_site_fract_y",
           "_atom_site_fract_z",
           "_atom_site_Cartn_x",
           "_atom_site_Cartn_y",
           "_atom_site_Cartn_z",
           "_atom_site_type_symbol",
           "_atom_type_symbol");
$isloop = 0;
$loophd = 0;
$lpfdct = 0;
$symlp  = 0;
$poslp  = 0;
$nofatom= 0;

# Parses command line options
GetOptions("help"=>\$helpflag,
           "interactive|i"=>\$interact,
           "debug"=>\$debug,
           "create|c"=>\$create,
           "supercell|s=s"=>\$supersize,
           "pointcharge|pc=f"=>\$pointcharge,
           "screenfile|sf=s"=>\$screeningfile,
           "savepcfile|sp"=>\$savepcfile,
           "savecharges|sc"=>\$savecharges,
           "readpcfile|rp"=>\$readpcfile,
           "charges|ch=s"=>\$inputcharges,
           "nonmagnetic|nm=s"=>\$nonmagnetic,
           "magnetic|m=s"=>\$magnetic,
           "so1ion"=>\$so1ion,
           "ic1ion"=>\$ic1ion,
           "phonon"=>\$phonon,
           "outpos"=>\$checkpos);

if (!$create && ($#ARGV<0 || $helpflag)) {
   print " $0:\n";
   print "   - script to convert a CIF into McPhase input files, similarly to powdercell2j\n\n";
   print " Syntax: $0 [CIFNAME] \n\n";
   print " where [CIFNAME] is a the name of the Crystallographic Information File.\n\n";
   print " Options include:\n";
   print "    --help        or -h  : prints this message\n";
   print "    --create      or -c  : creates a blank CIF file without and structure information.\n";
   print "    --supercell   or -s  : creates a super-cell of size #x#x#. (e.g. -s 2x3x4).\n";
   print "    --interactive or -i  : prompts user for information such as valence states\n";
   print "    --pointcharge or -pc : calculates the crystal field parameters from point charges\n";
   print "                           up to # Angstrom away from magnetic ions (e.g. -p 3.5)\n";
   print "    --savepcfile  or -sp : keep *.pc coordinate files in results folder\n";
   print "    --savecharges or -sc : store pointcharge coordinates in sipf file\n";
   print "    --readpcfile  or -rp : read from *.pc coordinate files in results folder\n";
   print "    --charges     or -ch : input charges to override defaults\n";
   print "    --nonmagnetic or -nm : set nonmagnetic ion overrideing default, e.g. -nm Cu\n";
   print "    --magnetic    or -m : set magnetic ion overrideing default, e.g. -m Cu\n";
   print "    --so1ion      or -so : force use of so1ion for single ion modules.\n";
   print "    --ic1ion      or -ic : force use of ic1ion for single ion modules.\n";
   print "    --phonon      or -ph : force use of phonon for single ion modules.\n";
   print "    --screenfile  or -sf : apply screening function to charges for calc of CF pars\n";
   print "                           e.g. -sf sf.r reads sf.r col 1 r(A) col2 sf(B2m) col3 sf(B4m) col4 sf(B6m)\n";
   print "\n";
   print " By default, this script is automatic, so if the oxidation (valence) states are not\n";
   print "   given in the CIF, it will be guessed at based on the element, as is whether the ion\n";
   print "   is magnetic or not.\n\n";
   print " Also by default, $0 will overwrite all output files (*.j, *.sipf)!\n\n";
   print " If you don't have a CIF of the structure you're studying, you can create a blank CIF\n";
   print "   using the -c option, then fill in the required information (lattice parameters and\n";
   print "   inequivalent site positions, and then rerun $0 on this CIF.\n\n";
   print " For the -pc pointcharge option, the -sp option can be used to keep generated files\n";
   print "   containing table of neighbouring charges named results/<sipfname>.pc\n";
   print "   You can then edit the values of the charges in this and rerun $0 with the -rp\n";
   print "   option to re-read this table with new charges to generate CF parameters in the sipf\n";
   print " You can also input the charges for each element using the -ch option with the syntax:\n";
   print "   $0 -pc <n> -ch O=-2,Ce=+3,Pd=+2\n\n";
   print " Finally, by default $0 will use so1ion for f-electron ions and ic1ion for d-electron\n";
   print "   ions, and phonon for others but this can be overridden using the -so or -ic pr -ph flags.\n";
   exit(0);
} else { print "* cif2mcphas 240328 *\n";}

if ($#ARGV<0) { $cif = "new.cif"; } else { $cif = $ARGV[0]; }

if ($create) {
   if (-e $cif) { 
     print "Warning: File $cif already exists. Do you want to overwrite this? (y/n)";
     $ans = <>; $ans =~ s/[\R\n\r]*//g;
     if($ans !~ /[Yy]/) { print "Exiting without creating file.\n"; exit(0); }
   }
   open (FOUT, ">$cif");
   print FOUT <<EOF;
; Please replace text in square brackets [] with your data
_cell_length_a                     [a_in_Angstrom]
_cell_length_b                     [b_in_Angstrom]
_cell_length_c                     [c_in_Angstrom]
_cell_angle_alpha                  [alpha_in_degrees]
_cell_angle_beta                   [beta_in_degrees]
_cell_angle_gamma                  [gamma_in_degrees]
; The H-M symbol generally needed by cif2mcphas but some spacegroups
; which have only one setting can be identified purely by the number.
; Most orthorhombic or low symmetry spacegroups have multiple centring,
; origin choices or unique axes which will require an H-M symbol.
; If this is the case, and you give just the spacegroup number,
; cif2mcphas may return with an error.
; If you don't use one of these fields, you must comment it out, by
; putting a ";" or "#" at the start of the line.
_symmetry_space_group_name_H-M     [Hermann-Maguin_symbol_with_spaces]
_symmetry_Int_Tables_number        [SpacegroupNumber]
; At the end of the file, please list the non-equivalent
; sites' fractional coordinates of this structure.
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
[AtomSite] [x_coord] [y_coord] [z_coord]
EOF
   print "Empty CIF $cif created. Please edit it and rerun \"$0 $cif\"\n";
   exit(0);
}

# ------------------------------------------------------------------------------------------------------------------------ #
# Subroutines
# ------------------------------------------------------------------------------------------------------------------------ #
sub nearest {
  my ($p, $n) = @_;
  return substr( $n + ( '0.' . '0' x $p . '5' ), 0, $p + length(int($n)) + 1 );
}

# Works out largest common factor of an array.  http://www.perlmonks.org/?node_id=56906
sub gcf {
  my ($x, $y) = @_;
  ($x, $y) = ($y, $x % $y) while $y;
  return $x;
}

sub multigcf {
  my $x = shift;
  $x = gcf($x, shift) while @_;
  return $x;
}

sub trim {
  my $s = shift;
  $s =~ s/^\s+|\s+$//g;
  return $s;
}

# Lifted from makenn.pl - finds nearest neighbours of an ion within some supercell defined by limits in nmx
sub getneighbours {
  my ($i,$posref,$nmxref,$rtoijk,$inputcharges) = @_;
  my @retval = ();
  @pos = @{$posref};
  @nmx = @{$nmxref};
  @p0 = split(":",$pos[$i]);
  $rn = new PDL();
  @dlist = ();
  @rlist = ();
  @alist = ();
  @charg2 = ();
  @charg4 = ();
  @charg6 = ();
  if ($inputcharges) {
    %chhash = split(/,|=/, $inputcharges);
  }
  if ($debug) { print "----------------\n".join("|",@p0)."\n----------------\n"; }
  for $i1 ($nmx[0]..$nmx[3]) {
    for $i2 ($nmx[1]..$nmx[4]) {
      for $i3 ($nmx[2]..$nmx[5]) {
        for (0..$#pos) {
          @ps = split(":",$pos[$_]);
          $dabc = pdl [ ($ps[1]-$p0[1]+$i1), ($ps[2]-$p0[2]+$i2), ($ps[3]-$p0[3]+$i3) ];
          #$rvec = transpose( $rtoijk x transpose($dabc) );  # For OpenBabel convention
          $rvec = $dabc x $rtoijk;                           # For McPhase convention
          $r = sqrt( inner($rvec, $rvec)->at(0) );
          if($r>0 && $r<=$pointcharge) {
            if ($debug) { print join("|",@ps)."\tr=$r\t$dabc\n"; }
            push @dlist, $dabc;              # Relative coordinates
            push @rlist, $rvec;              # Cartesian coordinates
            push @alist, $_;                 # Atom number
            if ($inputcharges) {
              $ps[0] =~ s/^\s+|\s+$//g;      # Trim whitespaces
                 if ($chhash{$ps[7]} ne "") { # if file index is given in charge list - use this changed MR 10.10.2023
              push @charg2, $chhash{$ps[7]};  # Uses user input charges
              push @charg4, $chhash{$ps[7]};  # Uses user input charges
              push @charg6, $chhash{$ps[7]};  # Uses user input charges
                                            } else { 
                 if ($chhash{$ps[8]} ne "") { # ... ifnot - > if atom label is given in charege list - use it
              push @charg2, $chhash{$ps[8]};  # Uses user input charges
              push @charg4, $chhash{$ps[8]};  # Uses user input charges
              push @charg6, $chhash{$ps[8]};  # Uses user input charges
                                            } else {
                if ($chhash{$ps[0]} ne "") { # ... ifnot - > if atom type is given in charge list - use it
              push @charg2, $chhash{$ps[0]};  # Uses user input charges
              push @charg4, $chhash{$ps[0]};  # Uses user input charges
              push @charg6, $chhash{$ps[0]};  # Uses user input charges
                                            } else { # no charge is given in charge list - > use oxydation state
              push @charg2, $ps[6];           # Oxidation state (valence / charge)
              push @charg4, $ps[6];           # Oxidation state (valence / charge)
              push @charg6, $ps[6];           # Oxidation state (valence / charge)
                                                   }
 }}
            } else {
              push @charg2, $ps[6];           # Oxidation state (valence / charge)
              push @charg4, $ps[6];           # Oxidation state (valence / charge)
              push @charg6, $ps[6];           # Oxidation state (valence / charge)
            }
            $rn = $rn->append( pdl([$r]) );  # Distance (Angstrom)
          }
        }
      }
    }
  }
  $n = qsorti($rn);
    if (! -d "results") { mkdir "results"; }
    print "atom ".($i+1)." ... results/$p0[7].pc running pointc ... \n";
    if(!$debug) { open (FOUT, ">results/$p0[7].pc"); } else { *FOUT = *STDOUT; }
    print FOUT "#-------------------------------------------------------------------------------------\n";
    print FOUT "#  table with neighbors and charges for atom ".($i+1)."\n";
    print FOUT "# output of program cif2mcphas, Reference: M. Rotter et al. PRB 68 (2003) 144418\n";
    print FOUT sprintf("#! ax= %+10.6f  ay=%+10.6f az=%+10.6f\n",$rtoijk->index(0)->at(0),$rtoijk->index(1)->at(0),$rtoijk->index(2)->at(0)); 
    print FOUT sprintf("#! bx= %+10.6f  by=%+10.6f bz=%+10.6f\n",$rtoijk->index(0)->at(1),$rtoijk->index(1)->at(1),$rtoijk->index(2)->at(1)); 
    print FOUT sprintf("#! cx= %+10.6f  cy=%+10.6f cz=%+10.6f\n",$rtoijk->index(0)->at(2),$rtoijk->index(1)->at(2),$rtoijk->index(2)->at(2)); 
    print FOUT "#-------------------------------------------------------------------------------------\n";
    if ( (abs($rtoijk->at(1,0))+abs($rtoijk->at(2,0))+abs($rtoijk->at(2,1))) > 0.001 ) {
      print FOUT "#orthonormal coordinate system xyz is defined with respect to abc as y||b, z||(a x b) and x normal to y and z\n";
      print FOUT "#c2[|e|]  dx[A]   dy[A]   dz[A]  c4[|e|] c6[|e|]      da[a]    db[b]    dc[c]   distance[A] atomnr\n";
    } else {
      print FOUT "#c2[|e|]  da[A]     db[A]     dc[A] c4[|e|] c6[|e|]          da[a]      db[b]      dc[c]     distance[A]   atomnr\n"; }
  
  for (0..$#rlist) {
    $id = $n->at($_+1)-1;
    $da = $dlist[$id];
    $xn = $rlist[$id];
    $outstr = sprintf("%8.4g   %+10.6f %+10.6f %+10.6f   %8.4g %8.4g      %+10.6f %+10.6f %+10.6f %+10.6f     %i\n",
      $charg2[$id],$xn->at(0,0),$xn->at(1,0),$xn->at(2,0),$charg4[$id],$charg6[$id],$da->at(0),$da->at(1),$da->at(2),$rn->index($n)->at($_+1),$alist[$id]+1);
     print FOUT $outstr; 
    push @retval, join(":",split(" ",$outstr));
  } 
  if(!$debug) { close FOUT; }
  return \@retval;
}

# ------------------------------------------------------------------------------------------------------------------------ #
# Reads CIF and parses it to get the structure parameters and atomic coordinates.
# ------------------------------------------------------------------------------------------------------------------------ #
while (<>) {
  $_ =~ s/^\s+|\s+$//g;     # safe chomp
  if ($_ =~ /#/) { next; }  # Ignore comments (for ICSD Karlsruhe data)
  if ($_ =~ /^;/) { next; } # Ignore comments (generally)
  if ($_ =~ /_cell/) {
    @line = split; 
    $cellpar{@line[0]} = @line[1];
  }
  if ($_ =~ /loop_/) {      # Check if start of loop structure. We're interested in two loops:
    $isloop = 1;            #   The symmpetry equivalent position loop
    $loophd = 1;            #   and the data loop.
    $symlp = 0;
    $poslp = 0;
    $othlp = 0;
    $lpfdct = 0;
    @poscol = ();
  } elsif($_ =~ /_symmetry_space_group_name_H-M/) { 
    $spagrp = $_;
    $spagrp =~ s/_symmetry_space_group_name_H-M//g;
    $spagrp =~ s/HR/H/; $spagrp =~ s/://g; $spagrp =~ s/['"]//g; $spagrp =~ s/^\s+//g;
  } elsif($_ =~ /_symmetry_Int_Tables_number/) {
    @line = split; $spanum = @line[1];
  } else {                  # First we parse the loop column headers in this else{} block
    if($isloop == 1) {
      if($_ =~ /^\s*_/) {
        if($loophd == 1) {
          if($_ =~ /_symmetry_equiv_pos_as_xyz/ || $_ =~ /_space_group_symop_operation_xyz/)  { 
            $symlp = 1; $symcol = $lpfdct; 
          }
          for $ic (0..$#datcol) {
            $_ =~ s/\s*//g;
            if($_ eq $datcol[$ic]) {
              if($ic>1 && $ic<8) { $poslp = 1; $coordseen[$ic] = 1;} else { $othlp = 1; }
              $poscol[$ic] = $lpfdct;
            }
          }
          $lpfdct++;
        } else {            # see a line with a _parameter -> end of loop
          $isloop = 0;
          $symlp = 0;
          $poslp = 0;
          $othlp = 0;
          $lpfdct = 0;
         #@poscol = ();
        }
      } else {              # Finished reading loop header, but still in the loop. 
        $loophd = 0;
        if($symlp==1) {     # We're in a symmetry equiv pos loop.
          push @sympos, $_;
        } 
        elsif($poslp>0 && $_!~/^\s*$/) {   # We're in the data loop.
          @line = split; 
          if($brokenline) {
            $brokenline = 0; 
            @line = (@line0,@line);
          } elsif($#line < $poscol[2]) {   # data continues on next line
            $brokenline = 1;
            @line0 = @line;
            next;
          }
          for $ic (0..$#datcol) {
            if(!($poscol[$ic]eq"")) { $sub= $line[$poscol[$ic]];
				    # if there are nonnumeric symbols refering to other atomic positions
				    # try to substitute those with numbers
				    if($ic>1&&$ic<8){for $j (0..$nofatom-1) {
                                      if($sub=~/.*$dat[$j][0]/){ #print "found $dat[$j][0] in $sub ... ";
                                                                 $sub=~s/$dat[$j][0]/\($dat[$j][$ic]\)/g; # substitute symbol with atomic coordinate
                                                                 #print "substituting to $sub\n";
                                                                              }
                                      } 
					# if there are errors given 0.4564(12) remove (12) else evaluate expression
                                     if($sub=~/\d+\(\d+\)/){$sub=sprintf "%.5f",$sub; }else{$sub=eval($sub);}
                                                    }
                                     $dat[$nofatom][$ic]=$sub;
                                    }
          }
          $nofatom++;
        }
        elsif($othlp>0) {   # We're in some other loop which has the fields we want
          @line = split;    #  e.g. ICSD (Karlsruhe) puts oxidation states in separate loops
          if(!($poscol[9]eq"")) {          # Use _atom_type_symbol to index these field
            for $ic (0..$#datcol) {
              if(!($poscol[$ic]eq"") && $ic!=9) {
                $otherdat{$line[$poscol[9]].":".$datcol[$ic]} = $line[$poscol[$ic]];
              }
            }
          }
        }
      }
    }
  }
}

# Check if we have a subsidiary table, and if so merge its data with the main table.
for (keys %otherdat) {
  @keydat = split(":");
  for $j(0..$nofatom-1) {
    if ($dat[$j][8] eq $keydat[0]) { 
      for $ic (0..$#datcol) {
        if($keydat[1] =~ /$datcol[$ic]/) { $dat[$j][$ic] = $otherdat{$_}; }
      }
    }
  }
}

# Got the data, now check it is ok
$errct=0;
foreach $cpn (@cpname) {
  if ($cellpar{$cpn} eq "") { print STDERR "Error: cell parameter $cpn not found.\n"; $errct++; }
}
if($errct>0) { die "Sorry, $0 is stupid and needs all lattice constants and angles. (Even for cubic structures)\n"
                  ."  Please edit $cif and make sure each of the above mention parameters is present.\n"; }

@vrn = ("","","x","y","z");
for $ic (2..4) {
  if(!$coordseen[$ic] && !$coordseen[$ic+3]) {
    die "Error: position data is missing the $vrn[$ic] coordinate. $0 requires either the fractional or Cartesian coordinates.\n"
       ."  Please Make sure the data _loop has a column called _atom_site_fract_$vrn[$ic] or _atom_site_Cartn_$vrn[$ic].\n"; }
}

# ------------------------------------------------------------------------------------------------------------------------ #
# Calculates the Cartesian to relative coordinate transformation and converts coordinates if necessary.
# ------------------------------------------------------------------------------------------------------------------------ #

# Determine the transformation matrix to convert from fractional to Cartesian coordinates and back.
#   this is not using the McPhase convention (y||b), but rather using the OpenBabel convention, with x||a
$PI=3.141592654;
$a = $cellpar{"_cell_length_a"};
$b = $cellpar{"_cell_length_b"};
$c = $cellpar{"_cell_length_c"};
# Removes uncertainties in cell parameters
$a =~ s/\([0-9]*\)//;
$b =~ s/\([0-9]*\)//;
$c =~ s/\([0-9]*\)//;
$alpha = $cellpar{"_cell_angle_alpha"}*$PI/180;
$beta  = $cellpar{"_cell_angle_beta"}*$PI/180;
$gamma = $cellpar{"_cell_angle_gamma"}*$PI/180;
$ca = cos($alpha); $cb = cos($beta); $cc = cos($gamma);
$v = 1-$ca*$ca-$cb*$cb-$cc*$cc+2*$ca*$cb*$cc;
if($debug==1) { print STDERR "v=$v\n"; }
if($v<0) { die "Unit cell parameters: alpha=$alpha,beta=$beta,gamma=$gamma are not geometrically consistent.\n"; }
$v=sqrt($v);
$rtoijk = pdl [ [ $a, $b*cos($gamma), $c*cos($beta) ],
                [  0, $b*sin($gamma), $c*(cos($alpha)-cos($beta)*cos($gamma))/sin($gamma) ],
                [  0,              0, $c*$v/sin($gamma) ] ];
#$invrtoijk = matinv($rtoijk);
$invrtoijk = inv($rtoijk);

if($debug==1) {
  print STDERR $rtoijk;
  print STDERR $invrtoijk;
}

# If the data is in Cartesians, convert it back to fractional before applying the symmetry equivalent positions
if($coordseen[2] eq "") {
  for $j (0..$nofatom-1) {
    $cpos = pdl [ ($dat[$j][5]), ($dat[$j][6]), ($dat[$j][7]) ];
    $fpos = $invrtoijk x transpose($cpos);
    $dat[$j][2] = sprintf "%.5f",$fpos->at(0,0); $dat[$j][2]*=1;
    $dat[$j][3] = sprintf "%.5f",$fpos->at(0,1); $dat[$j][3]*=1;
    $dat[$j][4] = sprintf "%.5f",$fpos->at(0,2); $dat[$j][4]*=1;
  }
} else {  # Otherwise, still round the data to 5 d.p. only, and convert to Cartesians for mcdiff.in
  for $j (0..$nofatom-1) { 
    $dat[$j][2] = sprintf "%.5f",$dat[$j][2]; $dat[$j][2]*=1;
    $dat[$j][3] = sprintf "%.5f",$dat[$j][3]; $dat[$j][3]*=1;
    $dat[$j][4] = sprintf "%.5f",$dat[$j][4]; $dat[$j][4]*=1;
  }
}

# ------------------------------------------------------------------------------------------------------------------------ #
# Determines the symmetry equivalent positions, oxidation states and other single ion parameters
# ------------------------------------------------------------------------------------------------------------------------ #

# If the spacegroup symbol is not given but the number is, use it instead if it is unique
if($spagrp eq "" && $spanum>0 && $spanum<231) {
  print STDERR "No Hermann-Maguin symbol detected in CIF, using ITC spacegroup number instead.\n";
  $posspagrp = $spaidx[$spanum-1];
  $possrhom = $rhoms{$spanum};
  if($#{$posspagrp}==0) {
    $spagrp = ${$posspagrp}[0];
    print STDERR "H-M symbol \"$spagrp\" assigned based on ITC spacegroup number $spanum.\n"; 
  } elsif(!($possrhom eq "") && 
       ( ((abs($a-$b)+abs($b-$c))<0.001 && (abs($alpha-$beta)+abs($beta-$gamma))<0.001) || 
         ($cellpar{"_cell_angle_gamma"}==120 && abs($a-$b)<0.001) )) { 
    # It's a rhombohedral group - try to work out if using rhombohedral or hexagonal cell.
    if(abs($b-$c)<0.001) { $spagrp = $possrhom." R"; } else { $spagrp = $possrhom." H"; }
    print STDERR "H-M symbol \"$spagrp\" assigned based on ITC spacegroup number $spanum and lattice parameters.\n"; 
  } else {
    die "Error: spacegroup number $spanum has multiple settings (centring, origin or unique axes) which cannot be determined solely from the ITC number. 
         Please edit the CIF and give the Hermann-Maguin symbol, or use Babel, SPACEGROUP or Vesta to generate another CIF.\n";
  }
}
# If there is no symmetry equivalent positions in the file, we look it up from a table.
if($symcol eq "" || join(",",@sympos)!~/[xyzXYZ]/) { 
  if(!($spagrp eq "")) {
    $symops = $symop{$spagrp};
    if($symops eq "") {
      $tstalias = $alias{$spagrp}; 
      if($tstalias eq "") { die "Error: Hermann-Maugin symbol '$spagrp' is not recognised.\n"; }
      if($tstalias =~ /S$/ && $spagrp !~ /1$/) { 
        print "--------------------------------------------------------------------------------\n";
        print "WARNING: Spacegroup '$spagrp' has 2 origin choices, but you have not specified\n";
        print "   which one to use. Assuming origin choice 1 (see International Tables A).\n";
        print "   This can lead to wrong atom positions (e.g. between spinel/pyrochlore lattices).\n";
        print "   To specify choice 1, use '$spagrp 1' or '$spagrp S' in future.\n";
        print "   To specify choice 2, use '$spagrp 2' or '$spagrp Z' in future.\n";
        print "--------------------------------------------------------------------------------\n"; }
      $symops = $symop{$tstalias};
    }
    if($spagrp =~ /^\s+H/) { 
      print STDERR "Warning: spacegroup $spagrp is an H-centred trigonal cell - this has not been tested.\n"; 
      print STDERR "Please check the generated coordinates in the mcphas.j file carefully,\n";
      print STDERR "or use mcphas2jvx to plot the positions to see they are what you expect.\n"
    }
    @sympos = split(";",$symops);
  } else {
    die "Error:   no _symmetry_equiv_pos_as_xyz or _symmetry_space_group_name_H-M field - I cannot determine the symmetry equivalent positions.
         Please edit the cif file to add these fields, or use a program like Babel or SPACEGROUP to generate another cif.\n"
  } 
} elsif($spagrp eq "" || ($symop{$spagrp} eq "" && $symop{$alias{$spagrp}} eq "")) {
# If the spacegroup string is undefined, attempt to work out what it is from the symmetry equivalent positions
  if(!($spagrp eq "")) { $spagrp0=$spagrp; }
  $maxmatch = 0;
  if($spanum>0 && $spanum<231) { 
    foreach (@{$spaidx[$spanum-1]}) { push @posspg, $_; }
  } else {
    for (keys %symop) { if($_!~/^H/) { push @posspg, $_; } }   # Ignore the H-centred groups for now...
  }
  foreach(@posspg) {
    @testsyms = split(";",$symop{$_});
    if($#testsyms==$#sympos) {
      if($debug) { print STDERR "$_  ---->  "; }
      $matchsum = 0;
      for $tstsymstr(@testsyms) { 
        $x = rand; $y = rand; $z = rand;
        $op = lc $tstsymstr; 
        $op =~ s/^[\-0-9\.]+\s//g; $op =~ s/^\s+[\-a-wA-W0-9]+\s+//g; $op =~ s/'\s+[0-9]+$//; $op =~ s/'//g; $op =~ s/([xyz])/\$$1/g; 
        $op =~ s/^/\$xx=/; $op =~ s/$/;/; $op =~ s/,/;\$yy=/; $op =~ s/,/;\$zz=/;
        eval $op;
        $match = 0;
        for $symposstr0(@sympos) {
          $symposstr=$symposstr0; $symposstr =~ s/^[a-wA-W0-9\s]*'//g; $symposstr =~ s/'//g; $symposstr =~ s/\s+//g;
          if($symposstr eq $tstsymstr) { $match=1; last; }
          else {
            $op = lc $symposstr0; 
            $op =~ s/^[\-0-9\.]+\s//g; $op =~ s/^\s+[\-a-wA-W0-9]+\s+//g; $op =~ s/'\s+[0-9]+$//; $op =~ s/'//g; $op =~ s/([xyz])/\$$1/g; 
            $op =~ s/^/\$x1=/; $op =~ s/$/;/; $op =~ s/,/;\$y1=/; $op =~ s/,/;\$z1=/;
            eval $op;
            if(abs($xx-$x1)<0.001 && abs($yy-$y1)<0.001 && abs($zz-$z1)<0.001) {
              $match=0.5; last;
            }
          }
        }
        if($match) { $matchsum+=$match; }
      }
      if($debug) { print STDERR "$matchsum\n"; }
      if($matchsum==$#sympos) { 
        $spagrp = $_; undef $maybematch; last;
      } elsif($matchsum>$maxmatch) {
        $spagrp = $_; $maybematch=1; $maxmatch = $matchsum; 
      }
    }
  }
  if(defined $spagrp0 && $maxmatch>0) { 
    print STDERR "WARNING: Hermann-Maugin symbol '$spagrp0' is not recognised. Guessing it should be '$spagrp'.\n"; 
  }
  if($maxmatch==0) { 
    print STDERR "WARNING: Hermann-Maugin symbol '$spagrp' ";
    if($spanum>0 && $spanum<231) { print STDERR "with ITC table number $spanum "; }
    print STDERR "is not recognised.\n";
  }
}

# Checks if oxidation state is given in the CIF file (either in the _oxidation_number or _type_symbol fields)
for $j(0..$nofatom-1) {
  push @oxyguess, 0;
  if(defined $dat[$j][1]) { 
    push @oxy, sprintf "%.0f", $dat[$j][1];  # Rounds the number - to give formal oxidation (valence) state
  } elsif(defined $dat[$j][8]) {
    $pp = $dat[$j][8]; $pp =~ s/[A-Za-z]*//;
    if($pp =~ /\+/) { $pp =~ s/\+//; push @oxy, $pp; } elsif($pp =~ /\-/) { $pp =~ s/\-//; push @oxy, (-1)*$pp; } else { push @oxy, ""; }
  } else {
    push @oxy, "";
    $oxyguess[-1]=1;
  }
  if($oxy[$j]==0) { $oxyguess[$j]=1; }
}
# If the oxidation state is undefined, look it up in the table
for $j(0..$nofatom-1) {
  $at = $dat[$j][0]; $at =~ s/_//g; $at =~ s/[0-9]+[A-Z]+//g; $at =~ s/[0-9]+//g; $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g;
  $at = lc $at; $at = ucfirst $at; 
  $attab = $element{$at}; 
  $atmas = $mass{$at}; 

  if(!$attab) { die "Error, unknown element: ".$dat[$j][0]."\n"; }
  if($oxy[$j] eq "" || $oxy[$j]==0) {
    $oxyref = ${$attab}[3];
    if($oxyref =~ /,/) { 
      @a_oxyref = split(",",$oxyref); $oxy[$j] = $a_oxyref[0];
    } else { $oxy[$j] = $oxyref; }
  }
  push @atomicmass, ${$atmas}[0]; 
  push @ismag, ${$attab}[2]; 
  push @realb, ${$attab}[0]/10;  # convert from femtometres to 10^-12cm
  push @imagb, ${$attab}[1]/10;  
  if($ismag[$j]) {
    $eltab = $magions{$at.$oxy[$j]."+"}; 
    if(${$eltab}[0] eq "") { $ismag[$j] = 0; }
    if($nonmagnetic eq $at){ $ismag[$j] = 0; }
    if($nonmagnetic eq $dat[$j][0]){ $ismag[$j] = 0; }
      }else
  {  if($magnetic eq $at){ $ismag[$j] = 1; }
     if($magnetic eq $dat[$j][0]){ $ismag[$j] = 1; }
  }
}

# If in interactive mode, asks the user if a site is magnetic or not, and (if needed) what valence
#   Note interactive mode is disabled if debug mode is active.
if($interact) {
  print "Please enter \"1\", \"yes\", \"y\", etc. if the following is true\n";
  print "\"0\", \"no\", \"n\", if false and a number value for the valence\n";
  print "if applicable (the \"+\" sign is not needed).\n";
  print "The default option is shown in brackets. If you accept the default,\n";
  print "just press \"enter\"\n";
  for $j(0..$nofatom-1) {
    $at = $dat[$j][0]; $at =~ s/_//g; $at =~ s/[0-9]+[A-Z]+//g; $at =~ s/[0-9]+//g; $at =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g;
    $at = lc $at; $at = ucfirst $at;
   #if($oxyguess[$j]) {
      print "What is the valence (oxidation state) of ions on site $dat[$j][0]? (".($oxy[$j]>0?"+":"-")."$oxy[$j])\n";
      if(!$debug) {
        $ans = <>; $ans =~ s/[\R\n\r]*//g;
        if(!($ans eq "")) { 
          if($ans !~ /[0-9\-]/) { 
            print "Sorry, I don't understand the valence state $ans. Assuming default valence: ".($oxy[$j]>0?"+":"-")."$oxy[$j]\n"; 
          } else {
            $ans =~ s/\+//g;
            $oxy[$j] = $ans;
            $eltab = $magions{$at.$oxy[$j]."+"};
            if(${$eltab}[0] eq "") { $ismag[$j] = 0; } else { $ismag[$j] = 1; }
print $dat[$j][0]."\n";
          }}
            if($nonmagnetic eq $at){ $ismag[$j] = 0; }
            if($magnetic eq $at){ $ismag[$j] = 1; }
            if($nonmagnetic eq $dat[$j][0]){ $ismag[$j] = 0; }
            if($magnetic eq $dat[$j][0]){ $ismag[$j] = 1; }
           
        
      }
   #}
    print "Are the ions on site $dat[$j][0] ($at".abs($oxy[$j]).($oxy[$j]>0?"+":"-").") magnetic? (".($ismag[$j]?"true":"false").")\n";
    if(!$debug) { 
      $ans = <>; $ans =~ s/[\R\n\r]*//g;
      if(!($ans eq "")) {
        if($ans=~/[yYtT1]/) { 
          $ismag[$j] = 1;
        } else { 
          $ismag[$j] = 0;
        }
      }
    }
  }
}

# Prints the extracted / converted data
if($debug==1) {
  print STDERR "Spacegroup:|$spagrp|\n";
  foreach $cpn (@cpname) {
    print STDERR "cellpar{$cpn}=$cellpar{$cpn}\n";
  }
  print STDERR "symcol = $symcol\n";
  foreach $ic (0..$#datcol) {
    print STDERR "poscol[$ic]=$poscol[$ic]|$datcol[$ic]\n";
  }
  foreach $so (@sympos) {
    print STDERR "$so\n";
  }
  print STDERR "nofatom=$nofatom\n";
  for $j(0..$nofatom-1) {
    for $ic(0..$#datcol) {
      print STDERR "$dat[$j][$ic]\t";
    }
    print STDERR "\n";
  }
  print STDERR "        mass\toxidation\tismag\tRe[b]\tIm[b]:\n";
  for $j(0..$nofatom-1) {
    print STDERR "$dat[$j][0]\t$atomicmass[$j]\t$oxy[$j]\t$ismag[$j]\t$realb[$j]\t$imagb[$j]\n"; 
  }
}

# Checks that each site position is unique
@same = (); @isame = ();
foreach $ii(0..$nofatom-1) {
  $atom = $dat[$ii][0]; $atom =~ s/_//g; $atom =~ s/[0-9]+[A-Z]+//g;
  $atom =~ s/[0-9]+//g; $atom =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g;
  $atom = lc $atom; $atom1 = ucfirst $atom;
  foreach $jj($ii+1..$nofatom-1) {
    $atom = $dat[$jj][0]; $atom =~ s/_//g; $atom =~ s/[0-9]+[A-Z]+//g;
    $atom =~ s/[0-9]+//g; $atom =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g;
    $atom = lc $atom; $atom2 = ucfirst $atom;
    if($atom1 eq $atom2 && !($atom1 eq -1) && !($atom2 eq -1)) {
      $df = abs($dat[$ii][2]-$dat[$jj][2])+abs($dat[$ii][3]-$dat[$jj][3])+abs($dat[$ii][4]-$dat[$jj][4]);
      if($df<1e-3) { push @same,$jj; push @isame, $ii; }
    }
  }
}
foreach $ii(0..$#same) { 
  print STDERR "Warning: site $dat[$same[$ii]][0] has the same coordinates as site $dat[$isame[$ii]][0].". 
               " $dat[$same[$ii]][0] will be ignored in further calculations.\n";
  $dat[$same[$ii]][0] = -1;
}

if ($inputcharges) {
  %chhash = split(/,|=/, $inputcharges);
}

# Loops over each inequivalent atom and apply each symmetry equivalent operator
for $j(0..$nofatom-1) {
  if($dat[$j][0] eq -1) { push @atoms, -1; next; }
  @epos = ();
  foreach $sym (@sympos) {
    $x = $dat[$j][2];
    $y = $dat[$j][3];
    $z = $dat[$j][4];
    $op = lc $sym; 
    $op =~ s/^[\-0-9\.]+\s//g;
    $op =~ s/^\s+[\-a-wA-W0-9]+\s+//g;
    $op =~ s/'\s+[0-9]+$//;
    $op =~ s/'//g; $op =~ s/([xyz])/\$$1/g; 
    $op =~ s/^/\$xx=/; $op =~ s/$/;/; 
    $op =~ s/,/;\$yy=/; 
    $op =~ s/,/;\$zz=/;
    if($debug) { print STDERR "x=$x; y=$y; z=$z; sym=$op\t---->\t"; }
    eval $op;
    if($debug) { print STDERR "xx=$xx; yy=$yy; zz=$zz;\n"; }
    push @epos, "$xx;$yy;$zz";
  }
  # uniq from: http://perlmaven.com/unique-values-in-an-array-in-perl
  %seen = (); @uepos = grep { ! $seen{$_}++ } @epos; @epos = ();
  foreach $eps (@uepos) {
    @seps=split(";",$eps); 
    for (0..2) { if($seps[$_]<0)  { $seps[$_]+=ceil(abs($seps[$_])); } }
    for (0..2) { if($seps[$_]>=1) { $seps[$_]-=floor(abs($seps[$_])); } }
    $st = "$seps[0];$seps[1];$seps[2]"; 
# $st =~ s/6667/6666/g; $st =~ s/3334/3333/g;  # removed by MR 5.2.2023
    push @epos, $st;
  }
  %seen = (); @uepos = grep { ! $seen{$_}++ } @epos; @epos = ();

  # Loops through and checks there are no duplicates
  @same = ();
  foreach $ii(0..$#uepos) {
    @seps=split(";",@uepos[$ii]); 
    foreach $jj($ii+1..$#uepos) {
      @scps=split(";",@uepos[$jj]); 
      $df = abs($seps[0]-$scps[0])+abs($seps[1]-$scps[1])+abs($seps[2]-$scps[2]);
      if($df<1e-3) { push @same,$jj; }
    }
  }
  foreach (@same) { $uepos[$_]=undef; }
  @uepos = grep defined, @uepos;

  $atom = $dat[$j][0]; $atom =~ s/_//g; $atom =~ s/[0-9]+[A-Z]+//g; 
  $atom =~ s/[0-9]+//g; $atom =~ s/['"\*()\?\+\-\~\^\,\.\%\\\>\=\/\|\[\]\{\}\$]//g; $atom = sprintf "%-4s",$atom;
  $atom = lc $atom; $atom = ucfirst $atom;
  push @atoms, $atom;
  if($ismag[$j]) { $ion = $atom.$oxy[$j]."p"; $ion=~s/\s+//g; } else { $ion = $atom; }

  # Populates the arrays of atom types and positions
  $ioncount = 1;
  foreach $eps (@uepos) {
    @seps=split(";",$eps); 
    $label = $dat[$j][0]; $label =~ s/\s*//g;$lab=$label;
    if ($pointcharge) {
      $label .= "_".$ioncount; $ioncount++;
    }
# increased position accuracy to 12 digits MR 5.2.2023
    push @pos, sprintf "%-4s:% 17.12f:% 17.12f:% 17.12f:%-5s:%d:% 10.5f:%s:%s",$atom,@seps,$ion,$j,$oxy[$j],$label,$lab;
    $seen=0; foreach(@atp) { if($_=~/$atom/) { $seen=1; } } 
    if(!$seen || !defined($atp[-1])) { push @atp, $atom; } $htp{$atom}++;
   #printf $label."% 14.5f% 14.5f% 14.5f\n",@seps;
   #$cpos = pdl [ $seps[0], $seps[1], $seps[2] ];
   #$fpos = $rtoijk x transpose($cpos);
   #$ceps[0] = sprintf "% 14.5f",$fpos->at(0,0); #$ceps[0]*=1;
   #$ceps[1] = sprintf "% 14.5f",$fpos->at(0,1); #$ceps[1]*=1;
   #$ceps[2] = sprintf "% 14.5f",$fpos->at(0,2); #$ceps[2]*=1;
   #print $atom.join(" ",@ceps)."\n";
    if($checkpos) { printf "%-4s% 10.5f% 10.5f% 10.5f\n",$atom,@seps; }
   $mults[$j]++;
  }
  if ($inputcharges) {    # Replace formal valence by user inputted charges
    $atom =~ s/^\s+|\s+$//g;
    if ($chhash{$atom} ne "") {
      $oxy[$j] = $chhash{$atom};
    }
# inserted MR 9.10.23 to be able to set charges for different atom_labels
    if ($chhash{$lab} ne "") {
      $oxy[$j] = $chhash{$lab};
    }

  }
}
#print "$_ $htp{$_}\n" for (keys %htp);
if($debug==1) { foreach (@pos) { print STDERR "$_\n"; } }

# Determine common factors for chemical formula 
for (keys %htp) { push @nat, $htp{$_}; } $denom = multigcf(@nat);
if($debug==1) { print STDERR "denom=$denom\n"; }

if($checkpos) { exit(0); }

# ------------------------------------------------------------------------------------------------------------------------ #
# Use pointc to calculate crystal field parameters from point charges up to R away (R specified on commandline)
# ------------------------------------------------------------------------------------------------------------------------ #
if ($pointcharge) {
  if (!$readpcfile) {
    @nmx = (0,0,0,0,0,0);
    # First determine maximum distance of a basis atom to origin of unit cell
    $distmax = 0;
    for (0..$#pos) {
      @ps = split(":",$pos[$_]);
    
     # check if atom is in unit cell
     # if not translate it into primitive unit cell
     $moved=0;
    while($ps[1]>=1){--$ps[1];$moved=1;}
    while($ps[2]>=1){--$ps[2];$moved=1;}
    while($ps[3]>=1){--$ps[3];$moved=1;}
    while($ps[1]<0){++$ps[1];$moved=1;}
    while($ps[2]<0){++$ps[2];$moved=1;}
    while($ps[3]<0){++$ps[3];$moved=1;}
    if($moved==1){print STDERR "Warning: atom at da=$ps[1] db=$ps[2] dc=$ps[3] was not in primitive unit cell - translated\n"; 
                 }
  # OLD limits nmin nmax
      $dabc = pdl [ ($ps[1]), ($ps[2]), ($ps[3]) ];
      $rvec = transpose( $rtoijk x transpose($dabc) ); 
      $r = sqrt( inner($rvec, $rvec)->at(0) );
      if ($r>$distmax) { 
        $distmax = $r; }
     
    }
    
    if ($debug) { print "distmax = $distmax \n"; }
    $distmax += $pointcharge;
    # Determine $nmin,$nmax by looking at a cube with side 3distmax
    for $i1 (-1,1) {
      for $i2 (-1,1) {
        for $i3 (-1,1) {
          $n = inner($invrtoijk, pdl[$i1*$distmax*1.5,$i2*$distmax*1.5,$i3*$distmax*1.5]);
          for $im (0..2) {
            if (($n->at($im))<$nmx[$im]) { 
              $nmx[$im] = int($n->at($im))-1; }
            if (($n->at($im))>$nmx[$im+3]) { 
              $nmx[$im+3] = int($n->at($im))+1; }
          }
        }
      }
    }
   # NEW 8.4.2024 FASTER: new tighter $nmin$nmax box using cross products
# determine $nmin,$nmax by looking at condition that
# plane of parallelepiped must be more distant than radius of sphere given in command line
# e.g.
$rvec=crossp($rtoijk->slice(":,(1)"),$rtoijk->slice(":,(2)")); # $rvec = $r2 x $r3
$rr=inner($rvec, $rvec);
$rvec/=sqrt($rr);
# (n1max-1) * $r1 . $rvec )> $distmax  ... i.e. n1max > 1.0 + $distmax / ($r1 . $rvec)
# (n1min+1) * $r1. $rvec) < -$distmax           n1min < -1.0 - $distmax / ($r1 . $rvec)
$nmx[3]=my_ceil(1.0+$distmax  / inner($rtoijk->slice(":,(0)"),$rvec));
$nmx[0]=my_floor(-1.0-$distmax  / inner($rtoijk->slice(":,(0)"),$rvec));
# similar ...
$rvec=crossp($rtoijk->slice(":,(2)"),$rtoijk->slice(":,(0)")); # $rvec = $r3 x $r1
$rr=inner($rvec, $rvec);
$rvec/=sqrt($rr);
$nmx[4]=my_ceil(1.0+$distmax  / inner($rtoijk->slice(":,(1)"),$rvec));
$nmx[1]=my_floor(-1.0-$distmax  / inner($rtoijk->slice(":,(1)"),$rvec));
$rvec=crossp($rtoijk->slice(":,(0)"),$rtoijk->slice(":,(1)")); # $rvec = $r1 x $r2
$rr=inner($rvec, $rvec);
$rvec/=sqrt($rr);
$nmx[5]=my_ceil(1.0+$distmax  / inner($rtoijk->slice(":,(2)"),$rvec));
$nmx[2]=my_floor(-1.0-$distmax  / inner($rtoijk->slice(":,(2)"),$rvec));



    if ($debug) {    print "$nmx[0] to $nmx[3], $nmx[1] to $nmx[4], $nmx[2] to $nmx[5]\n"; }
    # Redefine B transformation matrix to follow McPhase convention (y||b)
    $ca = cos($alpha); $cb = cos($beta); $cc = cos($gamma);
    $v = sqrt(1-$ca*$ca-$cb*$cb-$cc*$cc+2*$ca*$cb*$cc);
    $rtoijkmc = pdl [ [ $a*sin($gamma),                                      $a*cos($gamma), 0 ],
                      [  0,                                                  $b,             0 ],             
                      [ $c*(cos($beta)-cos($alpha)*cos($gamma))/sin($gamma), $c*cos($alpha), $c*$v/sin($gamma) ] ];
    if ($debug) { print $rtoijkmc; }
  } 

  for $j (0..$#pos) {
    my @ps = split(":",$pos[$j]);
    if($ismag[$ps[5]]) { 
      $ionname = $ps[4]; $ionname=~s/\s*//g; $ionname =~ s/p/\+/g; $eltab = $magions{$ionname};   # Looks up <r^k> values
      $nofelectrons = ${$eltab}[0]; $nofelectrons =~ s/^[0-9][a-z]//;                             # Looks up nof_electrons
      $so1ionname = ${$eltab}[3]; if ($so1ionname eq "") { $so1ionname = "J=1"; }                 # And ionname for pointc
      $ionkey = $ps[7];
      if (!$readpcfile) {getneighbours($j,\@pos,\@nmx,$rtoijkmc,$inputcharges);
if(defined $screeningfile){
  print "# Applying screening of pointcharges from file $screeningfile.\n";
  system("mult 10 1 results/$ps[7].pc 1 2 ".$screeningfile);
  system("mult 10 5 results/$ps[7].pc 1 3 ".$screeningfile);
  system("mult 10 6 results/$ps[7].pc 1 4 ".$screeningfile);  
                  }
           }
        @neighbours = ();
        open(FIN, "results/$ps[7].pc");
        while (<FIN>) {
          $_ =~ s/\R//g; if ($_ =~ /#/) { next; } 
                                                # chomp and ignore comments
 
          @line = split; push @neighbours, join(":",@line);
        } 
        close(FIN);
     if(!$savepcfile&&!$readpcfile){unlink "results/$ps[7].pc";}
      $pointcstr = join("\n",@neighbours);
      $pointcstr = "$so1ionname\nnof_electrons=$nofelectrons\nR2=${$eltab}[5]\nR4=${$eltab}[6]\nR6=${$eltab}[7]\n".$pointcstr;
      $pointcstr =~ s/:/\ /g;
      if($debug) { print $pointcstr."\n"; }
      # Now pipe list of neighbours to pointc and get its reply.
      unless($savecharges){$omit="-o";}
      $pcpid = open3(\*PCIN,\*PCOUT,\*PCERR,'pointc '.$omit.' -b') || die "Cannot run pointc.\n";
      print PCIN $pointcstr."\n\n";
      # Read pointc results from stdout.
      @pc_head = ();
      @blm_par = ();		
      @llm_par = ();
      $part_id = 0;
      while (<PCOUT>) { 
        if($_ =~ /^#\-/) { $part_id+=0.5; }
        $_ =~ s/([LB]00)/#$1/;
           if($part_id==0)  { push @pc_head, $_; }
        elsif($part_id<1.5) { push @blm_par, $_; }
        elsif($part_id<2.5) { push @llm_par, $_; }
          }
      while(<PCERR>){print STDERR $_;}
      # Wait for pointc to finish so we don't leave zombie processes.
      waitpid($pcpid,0);
      @pc_head[0]="";  # Remove first blank line.
      if($debug) { print "----- pointc output -----".join("",@pc_head,@blm_par,@llm_par)."-------------------------\n"; }
      $Blm{$ionkey} = join("",@pc_head,@blm_par);
      $Llm{$ionkey} = join("",@pc_head,@llm_par);
    }
  }
}

# Convert back from radians
$alpha *= 180/$PI; $beta *= 180/$PI; $gamma *= 180/$PI;

# ------------------------------------------------------------------------------------------------------------------------ #
# Prints out information for user
# ------------------------------------------------------------------------------------------------------------------------ #
print "This is the structure and chemical data extracted from the CIF and/or guessed\n";
print "a=$a b=$b c=$c  alpha=$alpha beta=$beta gamma=$gamma\n";
if($spagrp eq "") {
  print "No spacegroup symbol found in CIF, but the following symmetry equivalent positions will be used:\n";
  foreach $so (@sympos) { print "$so\n"; }
} else { print "spacegroup is $spagrp".(defined $maybematch?" (possibly)":"")."\n"; }
print "--------------------------------------------------------------------------------\n";
print "Label\tElement\tMass\tValence\tMult.\tMagnetic?\tFract_x\tFract_y\tFract_z\n";
print "--------------------------------------------------------------------------------\n";
@magornot = ( "NonMagnetic", "Magnetic" );
for $j(0..$nofatom-1) {
  if(!($dat[$j][0] eq -1)) {
    print "$dat[$j][0]\t$atoms[$j]\t$atomicmass[$j]\t$oxy[$j]\t$mults[$j]\t$magornot[$ismag[$j]]\t$dat[$j][2]\t$dat[$j][3]\t$dat[$j][4]\n";
  }
}
print "--------------------------------------------------------------------------------\n";
print "Mult. is the multiplicity and is calculated by applying all the symmetry\n";
print "   equivalent positions to the coordinates found in the CIF. Note this\n";
print "   assumes full occupation of all sites. If this is wrong, please check\n";
print "   the CIF and file a bug.\n";
print "If the ionic valence or whether the ion is magnetic or not is incorrect,\n";
print "   please use the interactive option to specify it.\n";
print "If positions or structure parameters are not correct, please check the CIF file,\n";
print "   and file a bug.\n";

# Do we want a supercell? If so change a,b,c, lattice parameters to reflect bigger cell
if(defined $supersize) {
  @sss = split("x",$supersize); $sa = $sss[0]; $sb = $sss[1]; $sc = $sss[2];
  if($#sss>2) { $wantindex = 1; } else { $wantindex = 0; }
  print "\n";
  print "You requested a supercell of size $sa x $sb x $sc. Please check this is correct.\n";
  $a0=$a; $b0=$b; $c0=$c;
  $a *= $sa; $b *= $sb; $c *= $sc; 
  $nall = ($#pos+1) * $sa*$sb*$sc;
} else { $sa = 1; $sb = 1; $sc = 1; $nall = ($#pos+1); }

# ------------------------------------------------------------------------------------------------------------------------ #
# Outputs the data to mcphas.j and sipf files.
# ------------------------------------------------------------------------------------------------------------------------ #

# Prints out the mcphas_all_atoms.j without neighbours (use makenn.pl)
if($debug==0) { open (FOUT, ">mcphas_all_atoms.j"); } else { *FOUT = *STDOUT; }
$nnat=$#pos+1;
print FOUT "#<!--mcphase.mcphas.j-->\n";
print FOUT "#***************************************************************\n";
print FOUT "# Lattice and Exchange Parameter file for\n";
print FOUT "# mcphas version 5.2\n";
print FOUT "# - program to calculate static magnetic properties\n";
print FOUT "# reference: M. Rotter JMMM 272-276 (2004) 481\n";
print FOUT "# mcdisp version 5.2\n";
print FOUT "# - program to calculate the dispersion of magnetic excitations\n";
print FOUT "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n";
print FOUT "#***************************************************************\n";
print FOUT "#\n";
print FOUT "# "; foreach (@atp) { $at=$_; $fml=$htp{$_}/$denom; $at=~s/\s+//g; print FOUT "$at ($fml) "; } print FOUT "\n";
print FOUT "#\n";
print FOUT "# Lattice Constants (A)    Spacegroup=$spagrp".(defined $maybematch?" (possibly)":"")."\n";
print FOUT "#\n";
print FOUT "#! a= $a b= $b c= $c  alpha= $alpha beta= $beta gamma= $gamma\n";
print FOUT "#\n";
print FOUT "#! r1a=   1 r2a= 0 r3a=  0\n";
print FOUT "#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]\n";
print FOUT "#! r1c=   0 r2c= 0 r3c=  1\n";
print FOUT "#\n";
print FOUT "#! nofatoms= $nall  nofcomponents=3  number of atoms in primitive unit cell/number of components of each spin\n";
$aid = 1;$nmag=0;
for $si (0..$sa-1) { 
  for $sj (0..$sb-1) { 
    for $sk (0..$sc-1) {
      for (0..$#pos) {
        @ps = split(":",$pos[$_]); $ntype = $htp{$ps[0]}; $ps[0]=~s/\s*//g; $ps[4]=~s/\s*//g; 
        $ions{$ps[7]} = $ps[5]; $ionnames{$ps[7]} = $ps[4];
        if($pointcharge) { $ntype = 1; } else { $ntype = $mults[$ps[5]]; }
        $pa = ($ps[1]+$si)/$sa; $pb = ($ps[2]+$sj)/$sb; $pc = ($ps[3]+$sk)/$sc;
        print FOUT "#********************************************************************* \n";
        print FOUT "#ATOM TYPE $ps[0] ; number of the atom in the UNIT CELL = $aid ; number of the atom within this type = $ntype\n";
        print FOUT "#! da= $pa [a] db= $pb [b] dc= $pc [c] nofneighbours=0 diagonalexchange=1 sipffilename= $ps[7].sipf\n";
        $aid++;
        if($ismag[$ions{$ps[7]}]) { $nmag++; }
      }
    } 
  } 
}
print FOUT "#********************************************************************* \n";
if($debug==0) { close FOUT; }

# Prints out the mcphas_magnetic_atoms.j without neighbours (use makenn.pl)
if($debug==0) { open (FOUT, ">mcphas_magnetic_atoms.j"); }
print FOUT "#<!--mcphase.mcphas.j-->\n";
print FOUT "#***************************************************************\n";
print FOUT "# Lattice and Exchange Parameter file for\n";
print FOUT "# mcphas version 5.2\n";
print FOUT "# - program to calculate static magnetic properties\n";
print FOUT "# reference: M. Rotter JMMM 272-276 (2004) 481\n";
print FOUT "# mcdisp version 5.2\n";
print FOUT "# - program to calculate the dispersion of magnetic excitations\n";
print FOUT "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n";
print FOUT "#***************************************************************\n";
print FOUT "#\n";
print FOUT "# "; foreach (@atp) { $at=$_; $fml=$htp{$_}/$denom; $at=~s/\s+//g; print FOUT "$at ($fml) "; } print FOUT "\n";
print FOUT "#\n";
print FOUT "# Lattice Constants (A)    Spacegroup=$spagrp".(defined $maybematch?" (possibly)":"")."\n";
print FOUT "#\n";
print FOUT "#! a= $a b= $b c= $c  alpha= $alpha beta= $beta gamma= $gamma\n";
print FOUT "#\n";
print FOUT "#! r1a=   1 r2a= 0 r3a=  0\n";
print FOUT "#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]\n";
print FOUT "#! r1c=   0 r2c= 0 r3c=  1\n";
print FOUT "#\n";
print FOUT "#! nofatoms= $nmag  nofcomponents=3  number of atoms in primitive unit cell/number of components of each spin\n";
$aid = 1;$nnonmag = 0;
for $si (0..$sa-1) { 
  for $sj (0..$sb-1) { 
    for $sk (0..$sc-1) {
      for (0..$#pos) {
        @ps = split(":",$pos[$_]); $ntype = $htp{$ps[0]}; $ps[0]=~s/\s*//g; $ps[4]=~s/\s*//g; 
        $ions{$ps[7]} = $ps[5]; $ionnames{$ps[7]} = $ps[4];
        if($pointcharge) { $ntype = 1; } else { $ntype = $mults[$ps[5]]; }
        if($ismag[$ions{$ps[7]}]) {
          $pa = ($ps[1]+$si)/$sa; $pb = ($ps[2]+$sj)/$sb; $pc = ($ps[3]+$sk)/$sc;
          print FOUT "#********************************************************************* \n";
          if($wantindex==1) { print FOUT "# atomindex = [$aid](".($_+1)."|$si,$sj,$sk)\n"; }
          print FOUT "#ATOM TYPE $ps[0] ; number of the atom in the UNIT CELL = $aid ; number of the atom within this type = $ntype\n";
          print FOUT "#! da= $pa [a] db= $pb [b] dc= $pc [c] nofneighbours=0 diagonalexchange=1 sipffilename= $ps[7].sipf\n";
          $aid++;
        } else { $nnonmag++; }
      }
    } 
  } 
}
print FOUT "#********************************************************************* \n";
if($debug==0) { close FOUT; }
copy("mcphas_magnetic_atoms.j","mcphas.j");

# Generates the associated single ion parameter files (sipf). 
$sipfheader = 
  "#<!--mcphase.sipf-->\n".
  "#***************************************************************\n".
  "# Single Ion Parameter File \n".
  "# mcphas  - program to calculate static magnetic properties\n".
  "# reference: M. Rotter JMMM 272-276 (2004) 481\n".
  "# mcdisp  - program to calculate the dispersion of magnetic excitations\n".
  "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n".
  "# mcdiff - program to calculate neutron and magnetic xray diffraction\n".
  "# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n".
  "#***************************************************************\n";
for (keys %ions) {
  if($debug==0) { open (FOUT, ">$_.sipf"); }
  if($ismag[$ions{$_}]&&!$phonon) {
    $ionname = $ionnames{$_}; 
    $ionname =~ s/p/\+/g;
    $eltab = $magions{$ionname};             # Looks up information about the magnetic ions
    $fftab = $magff{$ionname};               #   and form factor
    $nofelectrons = ${$eltab}[0];
    if ($nofelectrons =~ /f/) {
      $so1ionname = ${$eltab}[3];
    } else { 
      $so1ionname = "S=".${$eltab}[1]; }     # Assumes always high spin
    if($so1ion) {
      $modulename = "so1ion";
      $iontype = $so1ionname;
    } elsif ($ic1ion) {
      $modulename = "ic1ion";
      $iontype = $ionname;
    } 
    else {
      if ($nofelectrons =~ /f/) {
        $modulename = "so1ion";
        $iontype = $so1ionname;
      } else {
        $modulename = "ic1ion";
        $iontype = $ionname;
      }
    }
    print FOUT "#!MODULE=$modulename\n";
    $nofelectrons =~ s/^[0-9][a-z]//;
    print FOUT $sipfheader;
    print FOUT "IONTYPE=$iontype\n";
    print FOUT "CHARGE=$oxy[$ions{$_}]\n";
    print FOUT "MAGNETIC=$ismag[$ions{$_}]\n";
    print FOUT "nof_electrons=$nofelectrons\n";
    print FOUT "\n";
    print FOUT "#----------------\n";
    print FOUT "# Lande factor gJ\n";
    print FOUT "#----------------\n";
    print FOUT "GJ=${$eltab}[4]\n";
    print FOUT "\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Radial integrals for point charge calculations\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "R2=${$eltab}[5]\n";
    print FOUT "R4=${$eltab}[6]\n";
    print FOUT "R6=${$eltab}[7]\n";
    print FOUT "\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)\n";
    print FOUT "#                      with s=sin(theta)/lambda=Q/4pi\n";
    print FOUT "# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>\n";
    print FOUT "# unit of DWF is [A^2]\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "DWF=0\n";
    print FOUT "\n";
    print FOUT "#Neutron magnetic form factor coefficients\n";
    print FOUT ${$fftab}[1];
    print FOUT "\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Neutron Scattering Length (10^-12 cm) (can be complex)\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "SCATTERINGLENGTHREAL=$realb[$ions{$_}]\n";
    print FOUT "SCATTERINGLENGTHIMAG=$imagb[$ions{$_}]\n";
    print FOUT "#  ... note: - if an occupancy other than 1.0 is needed, just reduce \n";
    print FOUT "#              the scattering length linear accordingly\n";
    $zk = $zktab{$ionname};
    if(!($zk eq "")) {
      print FOUT "\n";
      print FOUT "#----------------------------------------------------------------------\n";
      print FOUT "# coefficients of Z(K') according to Lovesey (Neutron Scattering) vol.2\n";
      print FOUT "# chapter 11.6.1 page 233: Z(K)= ZKcK-1 * <jK-1(Q)> + ZKcK+1 * <jK+1(Q)>\n";
      print FOUT "#  ... these coefficients are needed to go beyond dipolar approx.\n";
      print FOUT "#      for the neutron magnetic formfactor in rare earth ions\n";
      print FOUT "#----------------------------------------------------------------------\n";
      print FOUT $zk;
    }
    if ($modulename eq "ic1ion") {
      print FOUT "\n";
      print FOUT "#-----------------------------------------------------------------------\n";
      print FOUT "# Free-ion parameters for on-site spin-orbit and Coulomb interactions\n";
      print FOUT "# Obtained from optical spectroscopy data by Carnal et al., J. Chem Phys\n";
      print FOUT "# v90, p3443 (1989) and v96, p8713 (1992).\n";
      print FOUT "#-----------------------------------------------------------------------\n";
      print FOUT "units=meV\n";
      $conf=${$eltab}[0]; $conf=~ s/^[0-9]//;
      print FOUT "conf=$conf\n";
      print FOUT "zeta=${$eltab}[8]\n";
      print FOUT "F2=${$eltab}[9]\n";
      print FOUT "F4=${$eltab}[10]\n";
      if (${$eltab}[0] =~ /f/) {
        print FOUT "F6=${$eltab}[11]\n";
      }
    }
    if ($pointcharge) {
      print FOUT "\n";
      print FOUT "#------------------------------------------------------------------------\n";
      print FOUT "# Crystal field parameters (in meV) from point charges up to $pointcharge A\n";
      print FOUT "#------------------------------------------------------------------------\n";
      if ($modulename eq "so1ion") {
        print FOUT $Blm{$_};
      } else {
        print FOUT $Llm{$_};
      }
    }
  }
  else {
    print FOUT "#!MODULE=phonon\n";
    print FOUT $sipfheader;
    print FOUT "IONTYPE=$_".abs($oxy[$ions{$_}]).($oxy[$ions{$_}]>0?"+":"-")."\n";
    print FOUT "CHARGE=$oxy[$ions{$_}]\n";
    print FOUT "MAGNETIC=$ismag[$ions{$_}]\n";
    print FOUT "\n";
print FOUT << "EOF";
MODPAR1=$atomicmass[$ions{$_}]  #mass in(m0)
MODPAR2=0   # Kxx
MODPAR3=0   # Kyy
MODPAR4=0   # Kzz
MODPAR5=0  # Kxy  in (meV)
MODPAR6=0  # Kxz
MODPAR7=0  # Kyz
MODPAR8=0.4 # umax          maximum (cutoff) for displacement [a0=0.5219 A]
MODPAR9=0   #                0       umax restriction in all directions
            #                1,2,3   umax restriction in x y z direction only
            #                4       umax restriction in x and y direction
            #                5       umax restriction in x and z direction
            #                6       umax restriction in y and z direction
EOF
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)\n";
    print FOUT "#                      with s=sin(theta)/lambda=Q/4pi\n";
    print FOUT "# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>\n";
    print FOUT "# unit of DWF is [A^2]\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "DWF=0\n";
    print FOUT "\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "# Neutron Scattering Length (10^-12 cm) (can be complex)\n";
    print FOUT "#-------------------------------------------------------\n";
    print FOUT "SCATTERINGLENGTHREAL=$realb[$ions{$_}]\n";
    print FOUT "SCATTERINGLENGTHIMAG=$imagb[$ions{$_}]\n";
    print FOUT "#  ... note: - if an occupancy other than 1.0 is needed, just reduce \n";
    print FOUT "#              the scattering length linear accordingly\n";
  }
  if($debug==0) { close FOUT; }
}

# Prints out the mcdiff.in file
if($debug==0) { open (FOUT, ">mcdiff.in"); }
print FOUT << "EOF";
# this file is the input file created by program cif2mcphas 
#<!--mcdiff.mcdiff.in>
#***************************************************************
#      mcdiff is a program for the calculation of elastic
#   neutron diffraction and resonant magnetic Xray scattering
#  reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R
#***************************************************************
# this input file contains 4 sections corresponding to different
# groups of parameters
#
# - all lines have to start with a # sign with the  exception of
#   the lines containing atomic positional parameters
# - the other parameters have to be defined in the corresponding
#   section by statements such as parameter=value
# - the sequence of the parameters within a section is arbitrary
#
#
# %SECTION 1%  OVERALL PARAMETERS
#
#! lambda   = 2.4  wavelength (A)
#
#! thetamax = 60   maximum bragg angle (deg)
#
#! ovalltemp= 0  overall temperature factor (A^2)
#           ...I ~ EXP(-2 * ovalltemp * sintheta^2 / lambda^2)
#                  relation to other notations:
#                  ovalltemp = Biso = 8 pi^2 Uiso^2
#
#! lorentz=0  type of lorentzfactor to be used
#            0.....no lorentzfactor
#            1.....neutron powder flat sample
#            2.....neutron powder cylindrical sample
#            3.....neutron single crystal
#            4.....neutron TOF powder cyl. sample - d-pattern log scaled
#            5.....neutron TOF powder cyl. sample - d-pattern normal scaled
#
#! out10=1    type of desired output in column 10 and 11 of mcdiff.out
#! out11=0    (optional) default is NSF in column 10 and LF in column 11
#            0....LF
#            1....|NSF|[b]
#            2....Re(NSF)[b]
#            3....Im(NSF)[b]
#            4....|MSF|
#            5....|MSF.P|
#            6....Re(MSF.P)
#            7....Im(MSF.P)
#            8....|MSFdip|
#            9....|MSFdip.P|
#            10....Re(MSFdip.P)
#            11....Im(MSFdip.P)
#            12....angl(Q,P)[°]
#            13....i(MSFxMSF*).P
#            14....I+
#            15....I-
#            16....I+/I-
#            17....i(MSFxMSF*)dip.P
#            18....Idip+
#            19....Idip-
#            20....Idip+/Idip-
#            21....2*|MSF.P|/sin^2(angl(Q,P)
#            22....2*|MSFdip.P|/sin^2(angl(Q,P)
#            23....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
#            24....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
#            25....2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
#            26....2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
#
#
#           In the above the intensities I+ and I- are the intensities in a polarised neutron scattering experiment
#           with incoming polarisation parallel (I+) and antiparallel (I-) to P:
#            I+-=LF exp(-OTF Q^2/8pi^2)
#                    [ |NSF/NB|^2 + 3.65/4pi (|MSF|^2-i(MSF x MSF*).P)/NB^2
#                        +-  sqrt(3.65/4pi)/NB^2 (NSF (MSF*.P) + NSF* (MSF.P))]
#
#
#             For some of the above options we need the
#! Pa=  0.0000   Components of Polarisation Vector in terms of lattice vectors P=(Pa * a + Pb * b + Pc *c)
#! Pb=  0.0000   Note: the length of P, i.e. |P| indicates the degree of beam polarisation (|P|<=1)
#! Pc=  1.0000
#
#
#
#
# %SECTION 2% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL
#
#
#! natcryst=$nnonmag      number of nonmagnetic atoms in primitive crystalographic unit cell
#
# it follows a list of nat lines with nonmagnetic atoms
# ... notes: - if an occupancy other than 1.0 is needed, just reduce
#              the scattering length linear accordingly
#            - Debye Waller Factor notation: sqr(Intensity) ~ structure factor ~
#              ~sum_n ()n exp(-2 DWFn sin^2(theta) / lambda^2)=EXP (-Wn),
#              relation to other notations: 2*DWF = B = 8 pi^2 <u^2>, units DWF (A^2)
#
#! use_dadbdc=1
#            - 0 means: da db and dc are not used by the program (unless you enter a line #! use_dadbdc=1),
#               dr1,dr2 and dr3 refer to the primitive lattice given below
#
# Real Imag[scattering length(10^-12cm)]   da(a)    db(b)    dc(c)    dr1(r1)  dr2(r2)  dr3(r3)  DWF(A^2)
EOF
for $si (0..$sa-1) { 
  for $sj (0..$sb-1) { 
    for $sk (0..$sc-1) {
      for (0..$#pos) {
        @ps = split(":",$pos[$_]); $ntype = $htp{$ps[0]}; $ps[0]=~s/\s*//g; $ps[4]=~s/\s*//g; $ions{$ps[7]} = $ps[5];
        if($ismag[$ions{$ps[7]}]==0) {
          $pa = ($ps[1]+$si)/$sa; $pb = ($ps[2]+$sj)/$sb; $pc = ($ps[3]+$sk)/$sc;
          $fpos = pdl [ ($ps[1]+$si), ($ps[2]+$sj), ($ps[3]+$sk) ];
          printf FOUT "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f 0  # %s.sipf\n",
              $realb[$ions{$ps[7]}], $imagb[$ions{$ps[4]}], $pa, $pb, $pc, $fpos->at(0), $fpos->at(1), $fpos->at(2), $ps[7];
        }
      }
    } 
  } 
}
print FOUT << "EOF";
#
#
# %SECTION 3% DESCRIPTION OF THE LATTICE
#
#
# Note: what follows here may directly be taken from the output of program spins
#       (file spins.out) or charges (file charges.out)
# -----------------------------------------------------------------------------
EOF
if($nmag>0){
$na=1;$nb=1;$nc=1;$h=0;$k=0;$l=0;
if($interact){
 print STDOUT "Generating mcdiff.in file with possible magnetic structure -\n";
 print " Please enter the magnetic supercell dimension na nb nc and the propagation vector (h,k,l):\n";
 print " na (>=1) ? (1)\n";$ans = <>; $ans =~ s/[\R\n\r]*//g;if(!($ans eq "")) { $na=$ans;}
 print " nb (>=1) ? (1)\n";$ans = <>; $ans =~ s/[\R\n\r]*//g;if(!($ans eq "")) { $nb=$ans;}
 print " nc (>=1) ? (1)\n";$ans = <>; $ans =~ s/[\R\n\r]*//g;if(!($ans eq "")) { $nc=$ans;}
 print " h ? (0)\n";$ans = <>; $ans =~ s/[\R\n\r]*//g;if(!($ans eq "")) { $h=$ans;}
 print " k ? (0)\n";$ans = <>; $ans =~ s/[\R\n\r]*//g;if(!($ans eq "")) { $k=$ans;}
 print " l ? (0)\n";$ans = <>; $ans =~ s/[\R\n\r]*//g;if(!($ans eq "")) { $l=$ans;}
 }
 copy("mcphas_magnetic_atoms.j","mcphas.j");
 system ("spinsfromq $na $nb $nc $h $k $l > results/mcphas.tst");
 system ("spins -f results/mcphas.tst 1");

unless(open(FIN, "results/spins.out")){die "Error cif2mcphas reading results/spins.out\n";}
while (<FIN>) { last if ($_ =~/%SECTION 3/)   
              }
while(<FIN>){print FOUT $_;}
close FIN;
}else
{
print FOUT << "EOF";
# Lattice Constants (A)
#! a= $a b= $b c= $c  alpha= $alpha beta= $beta gamma= $gamma
#! r1a=   1 r2a=   0 r3a=   0
#! r1b=   0 r2b=   1 r3b=   0   primitive lattice vectors [a][b][c]
#! r1c=   0 r2c=   0 r3c=   1
#*****************
#!T=1 K Ha=0 T Hb= 0 T Hc= 1 T: nr1=1 nr2=1 nr3=1 nat=0 atoms in primitive magnetic unit cell:
EOF
}
if($debug==0) { close FOUT; }

# Creates a results folder for users
if (! -d "results") { mkdir "results"; }

#print "END OF PROGRAM $0\n";
print "\n";
print "Created files:\n";
print "\n";
print " mcphas_all.j                         ... contains all atoms (magnetic and not)\n";
print " mcphas_magnetic_atoms.j == mcphas.j  ... contains only magnetic atoms\n";
for (keys %ions) {
  $atm = $atoms[$ions{$_}]; $atm =~ s/\s+//g;
  printf " %-37s... contains parameters for %s ion\n", "$_.sipf", $atm.abs($oxy[$ions{$_}]).($oxy[$ions{$_}]>0?"+":"-");
}
print "\n";
print "   running \"makenn R\" command  will create from mcphas.j a number of files\n";
print "   \"makenn.aN.pc\" equal with the number of the atoms in the unit cell.\n";
print "   in this file you find the neighbours of the N atom in the unit cell.\n";
print "   It is important to know what number corresponds to each atom type!\n";
print "\n";
print " mcdiff.in                              ... input file for mcdiff\n";
print "   running javaview results/spins.jvx displays the corresponding\n";
print "   spinstructure (collinear), !! mind to edit formfactor information\n";
print "   in the magnetic ions .sipf files before running mcdiff !!\n";
print "\n";


# Rounds up a value, using int() function
sub my_ceil
{
    my $t = int($_[0]);
    return ($t != $_[0]) ? ($_[0] < 0) ? $t : $t + 1 : $_[0];
}

# Rounds down a value, using int() function
sub my_floor
{
    return ($_[0] < 0) ? int($_[0]) - 1 : int($_[0]);
}
