#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $target;
my $help;
my $output_only = 0;

GetOptions(
    'win'    => sub { set_target('win', \$target) },
    'lin'    => sub { set_target('lin', \$target) },
    'help|h' => \$help,
) or usage();

usage() if $help || !@ARGV;

my @files = @ARGV;
my $host = ($^O =~ /MSWin/i) ? 'win' : 'lin';
$target = $host unless defined $target;

my $out = "#!/bin/bash\nset -e\n";
$out .= bash_helper_function() if $target eq 'lin';

my @nodes;
for my $file (@files) {
    open(my $in, '<', $file)
        or die "mctest: cannot open '$file': $!\n";
    my @lines = <$in>;
    close($in) or die "mctest: cannot close '$file': $!\n";

    my $i = 0;
    my $parsed = parse_block($file, \@lines, \$i, 0);
    push @nodes, @$parsed;
}

my $script = render_script(\@nodes, $target, !$output_only);

# -win/-lin select the output format. If the selected format differs from
# the host, print the translated script without executing it.
if ($output_only) {
    print $script;
    exit 0;
}

my $start_time = time();
my $status = run_script($script, $target);

if ($status == 0) {
    my $elapsed = int(time() - $start_time);
    print "#************************************************\n";
    print "# mctest <$files[0]> in ${elapsed}s OK\n";
    print "#************************************************\n";
}

exit $status;

sub set_target {
    my ($new, $ref) = @_;
    die "mctest: choose only one of -win and -lin\n"
        if defined $$ref && $$ref ne $new;
    $$ref = $new;
    $output_only = 1;
}



sub parse_block {
    my ($file, $lines, $index_ref, $inside_loop) = @_;
    my @nodes;

    while ($$index_ref < @$lines) {
        my $line_no = $$index_ref + 1;
        my $line = $lines->[$$index_ref];
        ++$$index_ref;

        my $trim = $line;
        $trim =~ s/^\s+//;
        $trim =~ s/\s+$//;

        if ($trim =~ /^done\b/i) {
            die "mctest: unexpected 'done' at $file line $line_no\n"
                unless $inside_loop;
            return \@nodes;
        }

        if ($trim =~ /^for\s+([A-Za-z_]\w*)\s+in\s+(.*?)\s*;\s*do\s*$/i) {
            my ($var, $list) = ($1, $2);
            my $body = parse_block($file, $lines, $index_ref, 1);
            push @nodes, {
                type => 'loop',
                file => $file,
                line => $line_no,
                var  => $var,
                list => $list,
                body => $body,
            };
            next;
        }

        if ($trim =~ /^for\s+([A-Za-z_]\w*)\s+in\s+(.*?)\s+do\s*$/i) {
            my ($var, $list) = ($1, $2);
            my $body = parse_block($file, $lines, $index_ref, 1);
            push @nodes, {
                type => 'loop',
                file => $file,
                line => $line_no,
                var  => $var,
                list => $list,
                body => $body,
            };
            next;
        }

        push @nodes, {
            type => 'line',
            file => $file,
            line => $line_no,
            text => $line,
        };
    }

    die "mctest: missing 'done' for loop in $file\n" if $inside_loop;
    return \@nodes;
}

sub bash_helper_function {
    return <<'BASH';
mctest_helper() {
    local helper="$1"
    shift
    local output rc line name value

    # The extensionless legacy files are shell wrappers. Call their Perl
    # implementations directly instead of passing the wrappers to Perl.
    case "$helper" in
        getvalue|getvariable|rpvalue|chi2)
            helper="${helper}.pl"
            ;;
    esac

    output=$("${MCTEST_PERL:-perl}" "$MCPHASE_DIR/bin/$helper" "$@")
    rc=$?
    if [ "$rc" -ne 0 ]; then
        return "$rc"
    fi

    while IFS= read -r line; do
        case "$line" in
            "export "*)
                line=${line#export }
                ;;
            "set "*)
                line=${line#set }
                ;;
            *=*)
                ;;
            *)
                continue
                ;;
        esac

        name=${line%%=*}
        value=${line#*=}
        name=$(printf '%s' "$name" | tr -d '[:space:]')
        if [ -n "$name" ]; then
            export "$name=$value"
        fi
    done <<EOF
$output
EOF
}
BASH
}

sub helper_perl {
    return $ENV{MCTEST_PERL} if defined $ENV{MCTEST_PERL}
        && $ENV{MCTEST_PERL} ne '';
    return 'perl';
}


sub render_script {
    my ($nodes, $target, $diagnostics) = @_;

    my $out;
    if ($target eq 'lin') {
        $out = "#!/bin/bash\nset -e\n";
    } else {
        $out = "\@echo off\n";
    }
    $out .= bash_helper_function() if $target eq 'lin';
    $out .= render_nodes($nodes, $target, undef, $diagnostics);
    return $out;
}

sub render_nodes {
    my ($nodes, $target, $loop_var, $diagnostics) = @_;
    my $out = '';

    for my $node (@$nodes) {
        if ($node->{type} eq 'loop') {
            $out .= render_loop($node, $target, $diagnostics);
            next;
        }

        my $text = $node->{text};
        $text =~ s/\r?\n$//;
        my $translated = translate_line($text, $target, $loop_var);

        # Comments and blank lines are emitted without diagnostics.
        if ($text =~ /^\s*(?:#|rem\b|$)/i || !$diagnostics) {
            $out .= $translated . "\n";
            next;
        }

        my $where = $node->{file} . ':' . $node->{line};
        my $message = 'mctest: failed at ' . $where . ': ' . $text;

        if ($target eq 'lin') {
            $out .= "if { $translated; }; then :; else\n";
            $out .= "  rc=\$?\n";
            $out .= "  printf '%s\\n' " . shell_quote($message) . " >&2\n";
            $out .= "  exit \$rc\n";
            $out .= "fi\n";
        } else {
            $out .= $translated . "\n";
            $out .= "if errorlevel 1 (\n";
            $out .= "  echo $message\n";
            $out .= "  exit /b 1\n";
            $out .= ")\n";
        }
    }

    return $out;
}

sub render_loop {
    my ($node, $target, $diagnostics) = @_;
    my $var = $node->{var};
    my $list = $node->{list};
    my $out = '';

    if ($target eq 'lin') {
        $out .= "for $var in $list; do\n";
        $out .= render_nodes($node->{body}, $target, $var, $diagnostics);
        $out .= "done\n";
        return $out;
    }

    my @values = expand_batch_list($list, $node->{file}, $node->{line});
    my $batch_var = uc(substr($var, 0, 1));
    $out .= "for %%$batch_var in (" . join(' ', @values) . ") do (\n";
    $out .= render_nodes($node->{body}, $target, $batch_var, $diagnostics);
    $out .= ")\n";
    return $out;
}


sub expand_batch_list {
    my ($list, $file, $line) = @_;

    if ($list =~ /^\$\(seq\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\)$/) {
        my ($start, $step, $end) = ($1, $2, $3);
        die "mctest: zero seq step at $file line $line\n" if $step == 0;

        my @values;
        my $value = $start;
        my $guard = 0;
        while (($step > 0 && $value <= $end + 1e-9) ||
               ($step < 0 && $value >= $end - 1e-9)) {
            push @values, format_number($value);
            $value += $step;
            die "mctest: seq range too large at $file line $line\n"
                if ++$guard > 100000;
        }
        return @values;
    }

    # Ordinary lists and globs are passed through for cmd.exe to expand.
    $list =~ s/^\s+//;
    $list =~ s/\s+$//;
    return split(/\s+/, $list);
}

sub format_number {
    my ($n) = @_;
    my $s = sprintf('%.10f', $n);
    $s =~ s/0+$//;
    $s =~ s/\.$//;
    return $s;
}

sub translate_line {
    my ($line, $target, $loop_var) = @_;
    my $out = $line;

    if ($target eq 'lin') {
        $out =~ s/^\s*rem\b/#/i;
        $out =~ s/^\s*call\s+//i;
        $out =~ s/^\s*copy\s+/cp /i;
        $out =~ s/^\s*del\s+/rm /i;
        $out =~ s/^\s*set\s+([A-Za-z_]\w*)\s*=\s*(.*)$/export $1="$2"/i;
        $out =~ s/%([A-Za-z_]\w*)%/\$$1/g;

        # Convert legacy helper invocations to the common assignment protocol.
        if ($out =~ /^\s*(?:\.\s+)?(getvalue|getvariable|rpvalue|chi2)(\.pl)?\b(.*)$/i) {
            my ($helper, $suffix, $args) = ($1, $2, $3);
            $helper .= '.pl' unless defined $suffix;
            return "mctest_helper $helper$args";
        }
        if ($out =~ /^\s*\.\s+delEstr\b(.*)$/i) {
            return "MCPHASE_DELESTR=\$(perl \"\$MCPHASE_DIR/bin/delEstr.pl\"$1); " .
                "export MCPHASE_DELESTR";
        }
    } else {
        $out =~ s/^\s*export\s+([A-Za-z_]\w*)\s*=\s*(.*)$/set "$1=$2"/i;
        $out =~ s/^\s*cp\s+/copy /i;
        $out =~ s/^\s*rm\s+/del /i;
        $out =~ s/\$([A-Za-z_]\w*)/%$1%/g;
        $out =~ s/^\s*#/rem /;

        if ($loop_var) {
            $out =~ s/%\Q$loop_var\E%/%%$loop_var/g;
            $out =~ s/\$\Q$loop_var\E/%%$loop_var/g;
        }

        if ($out =~ /^\s*(?:\.\s+)?(getvalue|getvariable|rpvalue|chi2)(\.pl)?\b(.*)$/i) {
            my ($helper, $suffix, $args) = ($1, $2, $3);
            $helper .= '.pl' unless defined $suffix;
            return "call mctest_helper.bat $helper$args";
        }
    }

    return $out;
}

sub shell_quote {
    my ($s) = @_;
    $s =~ s/'/'\\''/g;
    return "'" . $s . "'";
}


sub run_script {
    my ($script, $target) = @_;
    my ($fh, $tmp) = tempfile('mctestXXXXXX',
        SUFFIX => ($target eq 'lin' ? '.sh' : '.bat'),
        UNLINK => 1);

    print {$fh} $script
        or die "mctest: cannot write temporary script: $!\n";
    close($fh)
        or die "mctest: cannot close temporary script: $!\n";

    if ($target eq 'lin') {
        chmod 0700, $tmp;
        system('bash', $tmp);
    } else {
        system($ENV{COMSPEC} || 'cmd.exe', '/c', $tmp);
    }

    return($? == -1 ? 127 : ($? >> 8));

}

sub usage {
    print STDERR "Usage: mctest.pl [-win|-lin] script [script ...]\n";
    exit 2;
}
