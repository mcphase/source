#!/usr/bin/perl
# 
# icsdread.pl
#
# Perl script to generate a CIF from a record of the Korean ICSD mirror (freely available).
# Search for the structure you want at http://icsd.kisti.re.kr/ and note the ICSD ID number.
# Then invoke this script with:
# 
# icsdread.pl [NUMBER] > outfile.cif
#
# Without the redirection, this script will print to standard output (the command line).
# The script will also read an html file of the ICSD record if that was previously downloaded.
#
# This program is part of the McPhase package, licensed under the GNU GPL v2. Please see the COPYING file
#
# Fri Sep 26 11:58:53 KST 2014 - Duc Le - mducle@snu.ac.kr

use LWP::Simple;
use Getopt::Long;

$debug=0;

# Loads the table of symmetry equivalent positions
push @INC, $ENV{'MCPHASE_DIR'}.'/bin/';
require 'itc_syms.pl';

# Parses command line options
GetOptions("help"=>\$helpflag,
           "debug"=>\$debug);

if ($#ARGV<0 || $helpflag) {
   print " $0:\n";
   print "   - script to generate a CIF from a record of the Korean ICSD mirror (free)\n\n";
   print "     Search for the structure you want at http://icsd.kisti.re.kr/ and note the ICSD ID number.\n";
   print "     Then invoke this script with:\n";
   print "         $0 [NUMBER] > outfile.cif\n\n";
   print " Options include:\n";
   print "    --help         - prints this message\n";
   print "\n";
   print " By default, this script outputs to STDOUT, so the redirction \">\" is needed to pipe to a file.\n";
   print " This script will also open a previously downloaded html file of the ICSD record if given that\n";
   print " instead of an ICSD ID number. (E.g. if the input is a valid file name and not a number)\n";
   exit(0);
}

# Checks whether to read from a file or download from website
if(-e $ARGV[0]) {
  open(FILE,$ARGV[0]) or die "Can't read file 'filename' [$!]\n"; $content = do { local $/; <FILE> }; close(FILE);
} elsif($ARGV[0] =~ /^\d+$/) {
  $content = get "http://icsd.kisti.re.kr/icsd/icsd_view1.jsp?num=$ARGV[0]";
} else {
  die "Error: $ARGV[0] is neither a file nor a ICSD ID number.\n";
}
@infile = split("\n",$content);
 
@fieldname = ("ICSD ID",
              "Chem. Name",
              "Structured formula",
              "Sum formula",
              "Reference",
              "Author",
              "Unit Cell",
              "Volume",
              "Z",
              "Space Group",
              "SG Number",
              "Crystal System",
              "PEARSON");

$tablestart=0;
$headerend=0;
$atomline1=0;
@field=();
@lines=();

foreach (@infile) {
  $_ =~ s/\R//g; 
  if($_=~/ICSD ID/) { $tablestart=1; }
  if($_=~/<tr>/) { 
    $isfield = 1;
    @field=();
    @lines=();
  }
  if($_=~/<\/tr>/) { $isfield = 0; }
  if($isfield==1) { 
    if($_!~/img src/) { 
      if($_=~/<td.*<\/td>/) {                  # Cell already on one line
        push @field, $_; 
      } else {                                 # Puts each table cell in one line <td> .. </td>
        $_=~ s/^\s+|\s+$//g;
        push @lines, $_;
        if($_=~/<\/td>/) { push @field, join(" ",@lines); @lines=(); }
      }
    } 
  }
  if($tablestart==1) {
    if($isfield==0 && $#field>0) {
      if($debug==1) { foreach $fieldline (@field) { print $fieldline."\n"; } print "------------------------------------\n"; }
      if($headerend==0) {
        foreach $fn (@fieldname) {
          if($field[0]=~/$fn/) {
            $val = $field[1]; 
            $val =~ s/<td[\sa-zA-Z0-9=\#\"]*>//g; 
            $val =~ s/^\s+|\s+$|<\/td>//g; 
            $val =~ s/<[a-zA-Z]*>//g;
            $val =~ s/\&nbsp/ /g;
            $fieldval{$fn} = $val;
          }
        }
      } else {                                 # Come to second table with atomic positions
        if($atomline1==0) {
          $atomline1=1;
          foreach $hd (@field) {
            $hd =~ s/<td[\sa-zA-Z0-9=\#\"]*>|^\s+|\s+$|<\/td>|<[\/a-zA-Z]*>//g; 
            push @atomhead,$hd;
          }
        } else {
          @atom=();
          foreach $hd (@field) {
            $hd =~ s/<td[\sa-zA-Z0-9=\#\"]*>|^\s+|\s+$|<\/td>|<[\/a-zA-Z]*>//g; 
            push @atom,$hd;
          }
          push @atoms, [ @atom ];
        }
      }
      @field = ();
    }
    if($_=~/<\/table>/) { 
      if($headerend==1) { $tablestart=0; } else { $headerend=1; }
    }
  }
}

if($debug==1) { 
  foreach $fn (@fieldname) {
    print "$fn=$fieldval{$fn}\n"; 
  }
  print join("\t|",@atomhead)."\n";
  foreach $atomref (@atoms) { 
    print join("\t ",@$atomref)."\n"; 
  }
}

$icsdid = $fieldval{"ICSD ID"};
$chemnm = $fieldval{"Chem. Name"};
$formul = $fieldval{"Structured formula"};
$smform = $fieldval{"Sum formula"};
$refern = $fieldval{"Reference"};
$author = $fieldval{"Author"}; $author =~ s/;/'\n    '/g; $author="    '$author'";
$uncell = $fieldval{"Unit Cell"};
  @cell=split(/\s+/,$uncell); $aa=$cell[1]; $bb=$cell[2]; $cc=$cell[3]; $al=$cell[4]; $be=$cell[5]; $ga=$cell[6];
$vv     = $fieldval{"Volume"};
$Z      = $fieldval{"Z"};
$spagrp = $fieldval{"Space Group"}; $spagrp =~ s/HR/H/;
$sgnum  = $fieldval{"SG Number"};
$crysys = $fieldval{"Crystal System"};
$pearsn = $fieldval{"PEARSON"};

$symops = $symop{$spagrp};
if($symops eq "") {
  $tstalias = $alias{$spagrp}; 
  if($tstalias eq "") { die "Error: Hermann-Maugin symbol '$spagrp' is not recognised.\n"; }
  $symops = $symop{$tstalias};
}
$sympos = "    '".$symops; $sympos =~ s/;/'\n    '/g; $sympos.="'";

$title = $formul; $title =~ s/\s+//g; 

print "data_$icsdid-ICSD\n";
print "_database_code_ICSD                $icsdid\n";
print "_chemical_name_systematic\n";
print "\'$chemnm\'\n";
print "_chemical_formula_structural\n";
print "\'$formul\'\n";
print "_chemical_formula_sum\n";
print "\'$smform'\n";
print "_citation\n";
print "\'$refern'\n";
print "loop_\n";
print "    _publ_author_name\n";
print "$author\n";
print "_cell_length_a                     $aa\n";
print "_cell_length_b                     $bb\n";
print "_cell_length_c                     $cc\n";
print "_cell_angle_alpha                  $al\n";
print "_cell_angle_beta                   $be\n";
print "_cell_angle_gamma                  $ga\n";
print "_cell_volume                       $vv\n";
print "_cell_formula_units_Z              $Z\n";
print "_symmetry_space_group_name_H-M     $spagrp\n";
print "_symmetry_Int_Tables_number        $sgnum\n";
print "loop_\n";
print "    _symmetry_equiv_pos_as_xyz\n";
print $sympos."\n";
print "loop_\n";
print "    _atom_site_label\n";
print "    _atom_type_oxidation_number\n";
print "    _atom_site_symmetry_multiplicity\n";
print "    _atom_site_Wyckoff_symbol\n";
print "    _atom_site_fract_x\n";
print "    _atom_site_fract_y\n";
print "    _atom_site_fract_z\n";
print "    _atom_site_occupancy\n";
print "    _atom_site_B_iso_or_equiv\n";
foreach $atomref (@atoms) { 
 #print "    ".join("",@$atomref[0..1])."\t@$atomref[2]\t".join("\t",split("",@$atomref[3]))."\t".join("\t ",@$atomref[4..$#$atomref])."\n"; 
  $ml = @$atomref[3]; $ml =~ s/[a-z]//g;
  $wy = @$atomref[3]; $wy =~ s/[0-9]//g;
  $np = join("\t ",@$atomref[4..$#$atomref]); $np =~ s/\([0-9]*\)//g; 
  print "    ".join("",@$atomref[0..1])."\t@$atomref[2]\t$ml\t$wy\t$np\n"; 
}
