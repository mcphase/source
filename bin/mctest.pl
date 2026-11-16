#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $target;
my $help;
my $output_only = 0;
my $while_id = 0;

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
    my $parsed = parse_block($file, \@lines, \$i, undef);

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
    my ($file, $lines, $index_ref, $terminator) = @_;
    my @nodes;

    while ($$index_ref < @$lines) {
        my $line_no = $$index_ref + 1;
        my $line = $lines->[$$index_ref];
        ++$$index_ref;

        my $trim = $line;
        $trim =~ s/^\s+//;
        $trim =~ s/\s+$//;

        # Close the current block.
        if ($trim =~ /^(done|fi)\b/i) {
            my $found = lc($1);
            die "mctest: unexpected '$found' at $file line $line_no\n"
                unless defined $terminator && $found eq $terminator;
            return \@nodes;
        }

        # Parse: if CONDITION; then   or   if CONDITION
        if ($trim =~ /^if\s+(.+?)\s*(?:;\s*then)?\s*$/i) {
            my $condition = $1;
            $condition =~ s/\s*;\s*then\s*$//i;

            # If "then" wasn't on the if line, require it on the next line.
            if ($trim !~ /;\s*then\s*$/i) {
                my $then_line = $lines->[$$index_ref] // '';
                $then_line =~ s/^\s+//;
                $then_line =~ s/\s+$//;
                die "mctest: expected 'then' after if at $file line $line_no\n"
                    unless $then_line =~ /^then$/i;
                ++$$index_ref;
            }

            my $body = parse_block($file, $lines, $index_ref, 'fi');
            push @nodes, {
                type => 'if',
                file => $file,
                line => $line_no,
                condition => $condition,
                body => $body,
            };
            next;
        }

        # Parse the two while forms used in the script:
        # while [ ... ]
        # do
        if ($trim =~ /^while\s+(.+?)\s*$/i) {
            my $condition = $1;
            $condition =~ s/\s*;\s*do\s*$//i;

            if ($trim !~ /;\s*do\s*$/i) {
                my $do_line = $lines->[$$index_ref] // '';
                $do_line =~ s/^\s+//;
                $do_line =~ s/\s+$//;
                die "mctest: expected 'do' after while at $file line $line_no\n"
                    unless $do_line =~ /^do$/i;
                ++$$index_ref;
            }

            my $body = parse_block($file, $lines, $index_ref, 'done');
            push @nodes, {
                type => 'while',
                file => $file,
                line => $line_no,
                condition => $condition,
                body => $body,
            };
            next;
        }

        # Parse for loops.
        if ($trim =~ /^for\s+([A-Za-z_]\w*)\s+in\s+(.*?)\s*;\s*do\s*$/i
            || $trim =~ /^for\s+([A-Za-z_]\w*)\s+in\s+(.*?)\s+do\s*$/i) {
            my ($var, $list) = ($1, $2);
            my $body = parse_block($file, $lines, $index_ref, 'done');
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

    die "mctest: missing '$terminator' in $file\n" if defined $terminator;
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
        $out = "\@echo off\nsetlocal EnableDelayedExpansion\n";
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
        if ($node->{type} eq 'if') {
            if ($target eq 'lin') {
                $out .= "if $node->{condition}; then\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics
                );
                $out .= "fi\n";
            } else {
                my $batch_condition =
                    batch_file_condition($node->{condition},
                                         $node->{file}, $node->{line});
                $out .= "if $batch_condition (\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics
                );
                $out .= ")\n";
            }
            next;
        }

        if ($node->{type} eq 'while') {
            if ($target eq 'lin') {
                $out .= "while $node->{condition}; do\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics
                );
                $out .= "done\n";
            } else {
                my ($var, $value) =
                    $node->{condition} =~
                    /^\[\s*\$([A-Za-z_]\w*)\s+-eq\s+(-?\d+)\s*\]$/;

                die "mctest: unsupported -win while condition at "
                  . "$node->{file}:$node->{line}: $node->{condition}\n"
                    unless defined $var;

                my $id = ++$while_id;
                my $start = "MCTEST_WHILE_${id}";
                my $end   = "MCTEST_ENDWHILE_${id}";

                $out .= ":$start\n";
                $out .= "if not \"%$var%\"==\"$value\" goto $end\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics
                );
                $out .= "goto $start\n";
                $out .= ":$end\n";
            }
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
        $out =~ s/\$([A-Za-z_]\w*)/!$1!/g;
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
         # Bash positional arguments.
        $out =~ s/\$\*/%*/g;
        $out =~ tr{/}{\\};

        # rm -> delete files. This handles file paths and wildcards, not directories.
        $out =~ s/^\s*rm\s+/del \/q /i;

        # cp -> overwrite destination files.
        $out =~ s/^\s*cp\s+/copy \/Y /i;

        # Convert trailing Bash background marker to cmd.exe's start command.
        if ($out =~ s/\s*&\s*$//) {
            $out = 'start "" /b ' . $out;
        }
        # Replace Bash kill commands with an editable taskkill filter.
        if ($out =~ /^\s*kill\b/) {
            return 'taskkill /FI "WINDOWTITLE eq PLEASE ENTER WINDOW TITLE" /IM *';
        }


    return $out;
}

sub shell_quote {
    my ($s) = @_;
    $s =~ s/'/'\\''/g;
    return "'" . $s . "'";
}

sub batch_file_condition {
    my ($condition, $file, $line) = @_;

    $condition =~ s/^\s+|\s+$//g;

    my ($negated, $path);

    if ($condition =~ /^test\s+-f\s+(.+)$/i) {
        ($negated, $path) = (0, $1);
    } elsif ($condition =~ /^test\s+!\s+-f\s+(.+)$/i) {
        ($negated, $path) = (1, $1);
    } elsif ($condition =~ /^\[\s*!?\s*-f\s+(.+?)\s*\]$/i) {
        $path = $1;
        $negated = ($condition =~ /^\[\s*!/i) ? 1 : 0;
    } else {
        die "mctest: unsupported -win if condition at "
          . "$file:$line: $condition\n";
    }

    $path =~ s/^\s+|\s+$//g;
    $path =~ s/^(['"])(.*)\1$/$2/;
    $path =~ tr{/}{\\};

    return ($negated ? 'not exist' : 'exist') . " \"$path\"";
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
    
print STDERR << "EOF";

program to test a simple script (windows and linux syntax are acceptable)
 
 usage: mctest.pl [-win|-lin] script [script ...]

will execute the script test.bat line by line and stop with an error message
if an error occurs.

Options: 
 -win / -lin     generate from the script  output for windows /linux , which can be
                pasted in  a  .bat / script file, respectively

For testing  goto directory demo and look out for testing batch files
named test*.bat

- use such a file to test mcphase, e.g.:
mctest test_ndcu2_mcphasit.bat

- to test a mcphase installation with all tests type
mctest test*.bat

- to test a specific command occuring in several test-batch files, e.g. makenn, try
to run all these batch files. First, check by searching all test-batch files with
grep makenn test*.bat
 ... and run them by:
grep -l makenn test*.bat | xargs mctest

- generate a list of test by 
ls -1 test* 
  and paste and copy some files into dd using
cat > dd
 terminate by Ctrl^C, then execute these tests by
cat dd | xargs mctest

Note on limitations: mctest will not process correctly some commands
EOF
    exit 2;
}
