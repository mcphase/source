#!/usr/bin/env perl

use strict;
use warnings;

use constant MEV_TO_KELVIN => 11.6045250061657;

my $number_re = qr{
    [+-]?
    (?:
        (?:\d+(?:\.\d*)?)
        |
        (?:\.\d+)
    )
    (?:[eEdD][+-]?\d+)?
}x;

my %theta_name_for = (
    2 => 'ALPHA',
    4 => 'BETA',
    6 => 'GAMMA',
);

# discover_b_parameters($filename) -> @names
#
# Return each active SIPF B_lm assignment once, in input order.
# Dies if the file cannot be read.
sub discover_b_parameters {
    my ($filename) = @_;

    open my $input, '<', $filename
        or die "$filename: cannot open for reading: $!\n";

    my (@names, %seen);
    while (my $line = <$input>) {
        $line = clean_parameter_line($line);
        next unless defined $line;

        while ($line =~ m{
            (?<![A-Za-z0-9_])
            (
                B(?:
                    \d+(?:-\d+|S)?
                )
            )
            (?=\s*=)
        }gx) {
            push @names, $1 unless $seen{$1}++;
        }
    }

    close $input
        or die "$filename: cannot close after reading: $!\n";

    return @names;
}

# parse_b_parameter_name($name) -> \%component
#
# Decode a SIPF name such as B43 or B43S into its rank, signed order, sine
# suffix, and output name. The S suffix denotes the sine/negative-order
# tesseral component. Dies when the name is ambiguous or |m| > l.
sub parse_b_parameter_name {
    my ($name) = @_;

    my ($rank, $signed_order, $sine_suffix);
    if ($name =~ /\AB(\d)(\d)(S?)\z/) {
        ($rank, $signed_order, $sine_suffix) = ($1, $2, $3);
    }
    elsif ($name =~ /\AB(\d)-(\d+)\z/) {
        ($rank, $signed_order, $sine_suffix) = ($1, -$2, '');
    }
    else {
        die "Cannot interpret SIPF parameter '$name' as B_lm\n";
    }

    my $order = abs($signed_order);
    die "Invalid crystal-field parameter '$name': |m|=$order exceeds l=$rank\n"
        if $order > $rank;
    die "Invalid crystal-field parameter '$name': do not combine a negative "
      . "m with an S suffix\n"
        if $signed_order < 0 && $sine_suffix;

    my $a_name;
    $a_name = 'A' . substr($name, 1);

    return {
        rank   => 0 + $rank,
        order  => $signed_order,
        sine   => $sine_suffix,
        a_name => $a_name,
    };
}

# clean_parameter_line($line) -> $line | undef
#
# Apply SIPF comment rules consistently in both file passes: ordinary '#'
# lines are ignored, while '#!' lines contain active assignments. Also removes
# line endings and trailing comments.
sub clean_parameter_line {
    my ($line) = @_;

    chomp $line;
    $line =~ s/\r\z//;
    return if $line =~ /^\s*#(?!\!)/;
    $line =~ s/^\s*#!\s*//;
    $line =~ s/#.*//;

    return $line;
}

# read_parameter_values($filename, @names) -> %value_for
#
# Read only the requested SIPF assignments. Scalar values and fitted
# "parName [current, min, max, ...]" values are accepted; the latter evaluate
# to current. Fortran D exponents are normalized. Dies on malformed or
# duplicate requested assignments, or when the file cannot be read.
sub read_parameter_values {
    my ($filename, @names) = @_;

    open my $input, '<', $filename
        or die "$filename: cannot open for reading: $!\n";

    my %value_for;
    my $line_number = 0;

    while (my $line = <$input>) {
        ++$line_number;
        $line = clean_parameter_line($line);
        next unless defined $line;

        for my $name (@names) {
            my $assignment = qr{
                (?<![A-Za-z0-9_])
                \Q$name\E
                (?![A-Za-z0-9_])
                \s*=\s*
            }x;

            next unless $line =~ /$assignment(.*)\z/;
            my $right_hand_side = $1;
            my $raw_value;

            # A fitted parameter has the form
            # parName [current, minimum, maximum, variation, step].
            # Its effective value is the first item ("current").
            if ($right_hand_side =~
                    /\Apar[A-Za-z0-9_]*\s*\[\s*($number_re)(?=\s*(?:,|\]))/) {
                $raw_value = $1;
            }
            elsif ($right_hand_side =~
                    /\A($number_re)(?=\s*(?:;|\z|[A-Za-z_]\w*\s*=))/) {
                $raw_value = $1;
            }
            else {
                die "$filename:$line_number: $name does not have a "
                  . "numeric or fitted-parameter value\n";
            }

            die "$filename:$line_number: $name is defined more than once\n"
                if exists $value_for{$name};

            $raw_value =~ tr/dD/eE/;
            $value_for{$name} = 0 + $raw_value;
        }
    }

    close $input
        or die "$filename: cannot close after reading: $!\n";

    return %value_for;
}

# usage() -> $text
#
# Return command-line help as a string so main can print it for --help or use
# it as the diagnostic for an invalid argument count.
sub usage {
    return <<"USAGE";
Usage: $0 SIPF_FILE

Find the crystal-field parameters B_lm in SIPF_FILE and convert each one
from meV to A_lm in K/a_0^l using B_lm = A_lm <r^l> theta_l.

Scalar assignments and fitted parameters of the form
  B20=parB20 [current, minimum, maximum, variation, step]
are supported.

The SIPF fields ALPHA, BETA, and GAMMA provide theta_l for ranks 2, 4,
and 6; the corresponding radial integrals are R2, R4, and R6. Every B_lm
field present at those ranks is converted, including odd-m and S components.
B00 is ignored because it is a rank-zero energy offset.
USAGE
}

# main(@arguments) -> $exit_status
#
# Validate the command line, discover and read the required SIPF values, then
# emit one A_lm value for every supported B_lm found. Rank zero is skipped
# because it is a common energy offset; unsupported ranks and missing factors
# are fatal.
sub main {
    my (@arguments) = @_;

    if (@arguments == 1
            && ($arguments[0] eq '-h' || $arguments[0] eq '--help')) {
        print usage();
        return 0;
    }

    die usage() unless @arguments == 1;
    my $input_file = shift @arguments;

    my @b_names = discover_b_parameters($input_file);
    die "$input_file: no crystal-field parameters B_lm found\n"
        unless @b_names;

    my (%component_for, @ranks, %rank_seen);
    for my $b_name (@b_names) {
        my $component = parse_b_parameter_name($b_name);
        $component_for{$b_name} = $component;

        my $rank = $component->{rank};
        next if $rank == 0;
        die "$input_file: cannot convert $b_name: SIPF defines Stevens "
          . "factors only for ranks 2, 4, and 6 "
          . "(ALPHA, BETA, GAMMA)\n"
            unless exists $theta_name_for{$rank};
        push @ranks, $rank unless $rank_seen{$rank}++;
    }

    my @names_to_read = @b_names;
    for my $rank (@ranks) {
        push @names_to_read, $theta_name_for{$rank}, "R$rank";
    }
    my %name_seen;
    @names_to_read = grep { !$name_seen{$_}++ } @names_to_read;

    my %parameter = read_parameter_values($input_file, @names_to_read);

    for my $b_name (@b_names) {
        my $component = $component_for{$b_name};
        my $rank      = $component->{rank};

        # B00 shifts every level equally and has no corresponding SIPF
        # ALPHA/BETA/GAMMA or radial-integral field.
        if ($rank == 0) {
            print STDERR
                "# alm_from_blm: skipping B00 rank-zero energy offset\n";
            next;
        }

        my $theta_name  = $theta_name_for{$rank};
        my $radial_name = "R$rank";
        die "$input_file: missing required SIPF parameter $theta_name\n"
            unless exists $parameter{$theta_name};
        die "$input_file: missing required SIPF parameter $radial_name\n"
            unless exists $parameter{$radial_name};

        my $theta  = $parameter{$theta_name};
        my $radial = $parameter{$radial_name};

        die "$input_file: cannot convert rank $rank: $theta_name is zero\n"
            if $theta == 0;
        die "$input_file: cannot convert rank $rank: $radial_name is zero\n"
            if $radial == 0;

        my $a = $parameter{$b_name} * MEV_TO_KELVIN
              / ($radial * $theta);

        print  "# Alm in units of K/a_0^$rank\n";
        printf "%s = %.17g\n", $component->{a_name}, $a;
    }

    return 0;
}

exit main(@ARGV);
