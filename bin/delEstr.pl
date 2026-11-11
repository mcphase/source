#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

use File::Basename qw(basename);
use File::Glob qw(bsd_glob GLOB_NOCHECK);

my $verbose = 0;
GetOptions('verbose|v' => \$verbose)
    or die "Usage: delEstr.pl [-v|--verbose] NAME=value[,value...] ...\n";

# Print usage and exit if no arguments were supplied.
if (!@ARGV) {
    print STDERR <<'USAGE';

Program delEstr to create a string for the reduce_unitcell -delatoms option
to remove bonding charges from makenn.j file after running makenn -cfph 
It reads results/makenn.j and results/makenn.a*.pc and sipf files.

Bonding Charges are distributed on nearest and next nearest and third nearest
neighbours and so on (sipf file names of these neighbours must not start with E and have the 
phonon module). The coefficients of distributing the bonding charge onto neighbours
are given as command line parameters:

Usage Windows:
 delEstr  [-v] ENAME=value[,value...] [ENAME=value[,value...] ...]

Usage Linux:
.  delEstr  [-v] ENAME=value[,value...] [ENAME=value[,value...] ...]
source delEstr  [-v] ENAME=value[,value...] [ENAME=value[,value...] ...]

Environmental variable MCPHASE_DELESTR is set to contain the distribution string to be
passed to reduce_unitcell, e.g. in the form of
  reduce_unitcell -delatoms $MCPHASE_DELESTR results/makenn.j > mnn.j

Example:
  delEstr.pl ETbO1=0.3,0.7

NAME matches a  bonding charge sipf filename ; for example, ETbO* matches ETbO1.sipf.
USAGE
    exit 1;
}

# Parse NAME=value,value,... arguments.
my %weights;
my @expanded_args;

for my $arg (@ARGV) {
    push @expanded_args, bsd_glob($arg, GLOB_NOCHECK);
}

for my $arg (@expanded_args) {
    for my $item (split /\s+/, $arg) {
        next if $item eq '';

        my ($name, $list) =
           $item =~ /\A([A-Za-z][A-Za-z0-9_.*?]*)=(.+)\z/
            or die "Invalid expression: '$item'\n";

        my @values = split /,/, $list, -1;
        for my $value (@values) {
            die "Invalid numeric value '$value' in '$item'\n"
                unless $value =~
                    /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/;

            die "Value '$value' in '$item' must be greater than 0 and less than 1\n"
                unless $value > 0 && $value <= 1;
        }

        die "Duplicate name '$name'\n" if exists $weights{$name};
        my $sum = 0;
$sum += $_ for @values;

if (sprintf('%.2f', $sum) ne '1.00') {
    die "Cannot renormalize '$name': sum is zero\n" if $sum == 0;

    $_ /= $sum for @values;

    my $renormalized = join(',', map { sprintf('%.2g', $_) } @values);
    warn "Warning: values for '$name' summed to "
        . sprintf('%.10g', $sum)
        . "; renormalized to $renormalized\n";
}

$weights{$name} = \@values;

    }
}

# Find a SIPF file in results/ or in the current directory.
sub find_sipf_file {
    my ($filename) = @_;

    for my $path ("results/$filename", $filename) {
        return $path if -f $path;
    }

    return;
}

# Return true if the first line of a SIPF file declares the phonon module.
# It accepts a file only if its first line starts
#  with #!MODULE=phonon and its filename does not start with E followed by an uppercase letter
sub is_phonon_sipf {
    my ($path) = @_;

    my $filename = basename($path);

    # Exclude names beginning with E followed by an uppercase letter,
    # for example EROC.sipf. Names such as ErOC.sipf are allowed.
    return 0 if $filename =~ /\AE[A-Z]/;

    open my $fh, '<', $path
        or die "Cannot open SIPF file '$path': $!\n";

    my $first_line = <$fh> // '';
    close $fh;

    return $first_line =~ /\A#!MODULE=phonon/ ? 1 : 0;
}


my $makenn_file = 'results/makenn.j';
open my $makenn_fh, '<', $makenn_file
    or die "Cannot open '$makenn_file': $!\n";

# Each sipffilename occurrence gets the next index, starting at 1.
my @sipfs;

while (my $line = <$makenn_fh>) {
    next unless $line =~ /^\s*#!/;

    while ($line =~ /\bsipffilename\s*=\s*([^\s#]+)/g) {
        my $filename = $1;
        my $path = find_sipf_file($filename);

        die "Cannot find SIPF file '$filename' (checked results/ and current directory)\n"
            unless defined $path;

        my $next_index = scalar @sipfs + 1;
my $phonon = is_phonon_sipf($path);

push @sipfs, {
    index    => $next_index,
    filename => $filename,
    path     => $path,
    phonon   => $phonon,
};



    }
}
close $makenn_fh;

# Map each makenn index to its phonon status.
# Expand command-line name patterns against the SIPF filename stems.
my %resolved_weights;
for my $pattern (keys %weights) {
    my $regex = quotemeta($pattern);
    $regex =~ s/\\\*/.*/g;
    $regex =~ s/\\\?/./g;

    my @matches;
    for my $sipf (@sipfs) {
        my $stem = basename($sipf->{filename});
        $stem =~ s/\.sipf\z//i;
        push @matches, $stem if $stem =~ /\A$regex\z/;
    }

    die "No SIPF file in '$makenn_file' matches command-line name '$pattern'\n"
        unless @matches;

    for my $stem (@matches) {
        die "More than one command-line pattern matches SIPF name '$stem'\n"
            if exists $resolved_weights{$stem};

        $resolved_weights{$stem} = $weights{$pattern};
    }
}

%weights = %resolved_weights;

my %phonon_by_index = map { $_->{index} => $_->{phonon} } @sipfs;
my %filename_by_index = map { $_->{index} => $_->{filename} } @sipfs;
my %matched_names;

my @output;

for my $source (@sipfs) {
    my $stem = basename($source->{filename});
    $stem =~ s/\.sipf\z//i;

    # Only process SIPF files that have a corresponding command-line name.
    next unless exists $weights{$stem};
    $matched_names{$stem} = 1;

    my $pc_file = "results/makenn.a$source->{index}.pc";
    open my $pc_fh, '<', $pc_file
        or die "Cannot open '$pc_file': $!\n";

    my @neighbours;
    while (my $line = <$pc_fh>) {
        next if $line =~ /^\s*#/;
        next if $line =~ /^\s*$/;

        my @columns = split /\s+/, $line;
        shift @columns if @columns && $columns[0] eq '';

        # Columns are 1-based: distance is column 8; neighbour index is column 9.
        next unless @columns >= 9;
        next unless $columns[7] =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/;
        next unless $columns[8] =~ /\A\d+\z/;

        my ($distance, $index) = (0 + $columns[7], 0 + $columns[8]);

        # Keep only neighbours whose SIPF file is phonon-module.
        next unless $phonon_by_index{$index};

        push @neighbours, {
            distance => $distance,
            index    => $index,
        };
    }
    close $pc_fh;

    # Sort by distance, then group equal distances into neighbour shells.
    @neighbours = sort {
        $a->{distance} <=> $b->{distance}
            || $a->{index} <=> $b->{index}
    } @neighbours;

    my @shells;
    for my $neighbour (@neighbours) {
        if (!@shells || $neighbour->{distance} != $shells[-1]{distance}) {
            push @shells, {
                distance   => $neighbour->{distance},
                neighbours => [],
            };
        }
        push @{ $shells[-1]{neighbours} }, $neighbour->{index};
    }

    my @parts = ($source->{index});
    my $values = $weights{$stem};
    my @verbose_parts;

    # Each input value applies to the corresponding neighbour shell.
    my $shell_count = @shells < @$values ? scalar(@shells) : scalar(@$values);
    for my $shell_number (0 .. $shell_count - 1) {
        my $shell = $shells[$shell_number];
        my $share = $values->[$shell_number] / @{ $shell->{neighbours} };

        for my $index (@{ $shell->{neighbours} }) {
    my $formatted = sprintf('%.2g', $share);
    $formatted =~ s/^0(?=\.)//;

    push @parts, "$index+--$formatted";

    if ($verbose) {
        my $filename = basename($filename_by_index{$index});
        push @verbose_parts, "$formatted * $filename";
        }
     }

    }

    push @output, join(':', @parts) if @parts > 1;
if ($verbose && @parts > 1) {
    print STDERR basename($source->{filename}), "  --> ",
                 join(' ', @verbose_parts), "\n";
}


}
my @unmatched = grep { !$matched_names{$_} } keys %weights;
die "No SIPF file in '$makenn_file' matches command-line name(s): "
    . join(', ', sort @unmatched) . "\n"
    if @unmatched;

print join(',', @output), "\n";
