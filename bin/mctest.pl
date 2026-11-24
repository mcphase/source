#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use File::DosGlob ();
use File::Spec;

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

my @files;

if ($^O =~ /MSWin/i) {
    for my $arg (@ARGV) {
        my @matches = File::DosGlob::glob($arg);

        if (@matches) {
            push @files, @matches;
        } else {
            # Keep the original argument so the later open reports
            # a useful error if the pattern matched nothing.
            push @files, $arg;
        }
    }
} else {
    @files = @ARGV;
}
warn "mctest: matched files: " . join(" ", @files) . "\n";

my $host = ($^O =~ /MSWin/i) ? 'win' : 'lin';
$target = $host unless defined $target;

my $marker;

if (!$output_only) {
    my ($marker_fh, $marker_path) = tempfile(
        'mctest-progressXXXXXX',
        SUFFIX => '.txt',
        UNLINK => 0
    );

    close($marker_fh)
        or die "mctest: cannot close progress marker: $!\n";

    $marker = File::Spec->rel2abs($marker_path);
}

my $start_time = time();

# Translation-only mode: render each file independently, but print all
# rendered scripts to stdout. Do not add execution diagnostics or markers.
if ($output_only) {
    for my $file (@files) {
        open(my $in, '<', $file)
            or die "mctest: cannot open '$file': $!\n";

        my @lines = <$in>;

        close($in)
            or die "mctest: cannot close '$file': $!\n";

        my $i = 0;
        my $nodes = parse_block($file, \@lines, \$i, undef);

        my @file_nodes = (
            {
                type => 'file_marker',
                file => $file,
            },
            @$nodes,
        );

        my ($script, $command_info) =
            render_script(\@file_nodes, $target, 0, undef);

        print $script;
    }

    exit 0;
}

# Normal execution: each input is parsed, rendered, and executed separately.
my $status = 0;

for my $file (@files) {
    open(my $in, '<', $file)
        or do {
            warn "mctest: cannot open '$file': $!\n";
            $status = 1;
            last;
        };

    my @lines = <$in>;

    close($in)
        or do {
            warn "mctest: cannot close '$file': $!\n";
            $status = 1;
            last;
        };

    my $i = 0;
    my $nodes = parse_block($file, \@lines, \$i, undef);

    my @file_nodes = (
        {
            type => 'file_marker',
            file => $file,
        },
        @$nodes,
    );

    my ($script, $command_info) =
        render_script(\@file_nodes, $target, 1, $marker);


    $status = run_script($script, $target, $marker, $command_info);

    # Stop immediately on the first failed or interrupted input file.
    last if $status != 0;
}

unlink $marker if defined $marker;

if ($status == 0) {
    my $elapsed = int(time() - $start_time);
    print "#************************************************\n";
    print "# mctest <" . join(' ', @files) . ">\n in ${elapsed}s OK\n";
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


sub is_batch_command {
    my ($line) = @_;

    $line =~ s/^\s+//;
    return 0 if $line eq '' || $line =~ /^(?:call|start|rem|echo|set|if|for|goto|cd|copy|del|exit|taskkill)\b/i;

    my $command;
    if ($line =~ /^"([^"]+)"/) {
        $command = $1;
    } elsif ($line =~ /^(\S+)/) {
        $command = $1;
    } else {
        return 0;
    }

    # If a command has an explicit extension, check that exact file.
    if ($command =~ /\.(?:bat|cmd)$/i) {
        return -f $command;
    }
    return 0 if $command =~ /\./;

    # Check the current directory and PATH for command.bat / command.cmd.
    my @dirs = (File::Spec->curdir(), File::Spec->path());
    for my $dir (@dirs) {
        for my $ext (qw(.bat .cmd)) {
            return 1 if -f File::Spec->catfile($dir, "$command$ext");
        }
    }

    return 0;
}

sub parse_block {
    my ($file, $lines, $index_ref, $terminator, $closed_by_ref) = @_;
    my @nodes;

    while ($$index_ref < @$lines) {
        my $line_no = $$index_ref + 1;
        my $line = $lines->[$$index_ref];
        ++$$index_ref;

        my $trim = $line;
        $trim =~ s/^\s+//;
        $trim =~ s/\s+$//;

                # Close the current block, or return the first branch at "else".
        if ($trim =~ /^(done|fi|else)\b/i) {
            my $found = lc($1);

            if ($found eq 'else'
                && defined($terminator)
                && $terminator eq 'fi') {
                $$closed_by_ref = 'else' if defined $closed_by_ref;
                return \@nodes;
            }

            die "mctest: unexpected '$found' at $file line $line_no\n"
                unless defined $terminator && $found eq $terminator;

            $$closed_by_ref = $found if defined $closed_by_ref;
            return \@nodes;
        }


                # Parse: if CONDITION; then  (or "if CONDITION" followed by "then")
        if ($trim =~ /^if\s+(.+?)\s*(?:;\s*then)?\s*$/i) {
            my $condition = $1;
            $condition =~ s/\s*;\s*then\s*$//i;

            if ($trim !~ /;\s*then\s*$/i) {
                my $then_line = defined($lines->[$$index_ref])
                    ? $lines->[$$index_ref] : '';
                $then_line =~ s/^\s+//;
                $then_line =~ s/\s+$//;

                die "mctest: expected 'then' after if at $file line $line_no\n"
                    unless $then_line =~ /^then$/i;
                ++$$index_ref;
            }

            my $closed_by;
            my $then_body = parse_block(
                $file, $lines, $index_ref, 'fi', \$closed_by
            );

            my $else_body = [];
            if (defined $closed_by && $closed_by eq 'else') {
                my $else_closed_by;
                $else_body = parse_block(
                    $file, $lines, $index_ref, 'fi', \$else_closed_by
                );

                die "mctest: missing 'fi' after else at $file line $line_no\n"
                    unless defined $else_closed_by
                        && $else_closed_by eq 'fi';
            }

            push @nodes, {
                type => 'if',
                file => $file,
                line => $line_no,
                condition => $condition,
                body => $then_body,
                else_body => $else_body,
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
                my $do_line = defined($lines->[$$index_ref])
    ? $lines->[$$index_ref] : '';
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
    my ($nodes, $target, $diagnostics, $marker) = @_;
    my %command_info = (_next_id => 0);

    my $out;
    if ($target eq 'lin') {
        $out = "#!/bin/bash\nset -e\n";
    } else {
        $out = "\@echo off\nsetlocal EnableDelayedExpansion\n";
    }

    $out .= bash_helper_function() if $target eq 'lin';

    $out .= render_nodes(
        $nodes, $target, undef, $diagnostics, $marker, \%command_info
    );

    return ($out, \%command_info);
}

sub render_nodes {
    my ($nodes, $target, $loop_var, $diagnostics, $marker, $command_info) = @_;
    my $out = '';

    for my $node (@$nodes) {
     

        if ($node->{type} eq 'file_marker') {
            my $file = $node->{file};
            $file =~ s/[\r\n]//g;
            $out .= "echo MCTEST_FILE: $file\n";
            next;
        }


        if ($node->{type} eq 'loop') {
            $out .= render_loop($node, $target, $diagnostics,$marker, $command_info);
            next;
        }
        if ($node->{type} eq 'if') {
            if ($target eq 'lin') {
                $out .= "if $node->{condition}; then\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics,$marker, $command_info
                );

                if (@{$node->{else_body} || []}) {
                    $out .= "else\n";
                    $out .= render_nodes(
                        $node->{else_body}, $target, $loop_var, $diagnostics,$marker, $command_info
                    );
                }

                $out .= "fi\n";
            } else {
                my $id = ++$while_id;
                my $else_label = "MCTEST_ELSE_${id}";
                my $end_label  = "MCTEST_ENDIF_${id}";
                my $condition = batch_file_condition(
                    $node->{condition}, $node->{file}, $node->{line}
                );

                $out .= "if not $condition goto $else_label\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics,$marker, $command_info
                );

                if (@{$node->{else_body} || []}) {
                    $out .= "goto $end_label\n";
                    $out .= ":$else_label\n";
                    $out .= render_nodes(
                        $node->{else_body}, $target, $loop_var, $diagnostics,$marker, $command_info
                    );
                    $out .= ":$end_label\n";
                } else {
                    $out .= ":$else_label\n";
                }
            }
            next;
        }


        if ($node->{type} eq 'while') {
            if ($target eq 'lin') {
                $out .= "while $node->{condition}; do\n";
                $out .= render_nodes(
                    $node->{body}, $target, $loop_var, $diagnostics,$marker, $command_info
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

my $batch_condition;
if ($node->{condition} =~ /^\[\s*\$([A-Za-z_]\w*)\s+-eq\s+(-?\d+)\s*\]$/) {
    $batch_condition = "!$1! EQU $2";
} else {
    $batch_condition = batch_file_condition(
        $node->{condition}, $node->{file}, $node->{line}
    );
}

$out .= ":$start\n";
$out .= "if not $batch_condition goto $end\n";
$out .= render_nodes($node->{body}, $target, $loop_var, $diagnostics,$marker, $command_info);
$out .= "goto $start\n";
$out .= ":$end\n";

            }
            next;
        }

      
        my $text = $node->{text};
        $text =~ s/\r?\n$//;
        my $translated = translate_line($text, $target, $loop_var);

        # Comments and blank lines are emitted without diagnostics.
        if (
    $text =~ /^\s*(?:#|rem\b|$)/i
    || ($target eq 'win' && $text =~ /^\s*chmod\b/i)
    || !$diagnostics
) {
            $out .= $translated . "\n";
            next;
        }
if (defined $marker && $diagnostics) {
    my $id = ++$command_info->{_next_id};

    $command_info->{$id} = {
    file    => $node->{file},
    line    => $node->{line},
    command => $translated,
};

    if ($target eq 'lin') {
        $out .= qq{printf '%s\\n' '$id' > "\$MCTEST_MARKER"\n};
    } else {
        $out .= qq{>"%MCTEST_MARKER%" echo $id\n};
    }
}


my $where = $node->{file} . ':' . $node->{line};
my $message;
if ($target eq 'lin') {
    $message = 'mctest: failed at ' . $where . ': ' . $text;

    $out .= "if { $translated; }; then :; else\n";
    $out .= "  rc=\$?\n";
    $out .= "  printf '%s\\n' " . shell_quote($message) . " >&2\n";
    $out .= "  exit \$rc\n";
    $out .= "fi\n";
} else {
    $message = 'mctest: failed at ' . $where . ': ' . $translated;

    # Parentheses in an echo inside a parenthesized batch block must be escaped.
    $message =~ s/([()])/^$1/g;

    $out .= "ver >nul\n";          # clear a stale ERRORLEVEL
$out .= $translated . "\n";
$out .= "set \"MCTEST_RC=!ERRORLEVEL!\"\n";
$out .= "if not \"!MCTEST_RC!\"==\"0\" (\n";
$out .= "  echo $message\n";
$out .= "  exit /b !MCTEST_RC!\n";
$out .= ")\n";

}

   }     

    return $out;
}

sub render_loop {
    my ($node, $target, $diagnostics, $marker, $command_info) = @_;
    my $var = $node->{var};
    my $list = $node->{list};
    my $out = '';

    if ($target eq 'lin') {
        $out .= "for $var in $list; do\n";
        $out .= render_nodes($node->{body}, $target, $var, $diagnostics,$marker, $command_info);
        $out .= "done\n";
        return $out;
    }

    my @values = expand_batch_list($list, $node->{file}, $node->{line});
    my $batch_var = uc(substr($var, 0, 1));
    $out .= "for %%$batch_var in (" . join(' ', @values) . ") do (\n";
    $out .= render_nodes($node->{body}, $target, $batch_var, $diagnostics,$marker, $command_info);
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
       $out =~ s{\\}{/}g;

    } else {
# chmod is a Unix command; preserve it as a comment in Windows output.
if ($out =~ /^(\s*)chmod\b/i) {
    $out = $1 . 'rem ' . $out;
    return $out;
}

        $out =~ s/^\s*export\s+([A-Za-z_]\w*)\s*=\s*(.*)$/set "$1=$2"/i;
        $out =~ s/^\s*cp\s+/copy /i;
        $out =~ s/^\s*rm\s+/del /i;
        $out =~ s/\$([A-Za-z_]\w*)/!$1!/g;
        $out =~ s/^\s*#/rem /;

        if ($loop_var) {
            $out =~ s/%\Q$loop_var\E%/%%$loop_var/g;
            $out =~ s/\$\Q$loop_var\E/%%$loop_var/g;
        }

                # Bash source/. command: call the corresponding batch file.
        if ($out =~ /^\s*(?:source|\.)\s+("?)([^\s"']+)\1(.*)$/i) {
            my ($file, $args) = ($2, $3);

            # Use the .bat counterpart for a sourced shell script;
            # add .bat when the command has no extension.
            $file =~ s/\.sh$/.bat/i;
            $file .= '.bat' unless $file =~ /\.(?:bat|cmd)$/i;

            # Single quotes don't quote arguments in cmd.exe.
            $args =~ s/'([^']*)'/"$1"/g;

            $file =~ tr{/}{\\};
            $args =~ s/\$([A-Za-z_]\w*)/!$1!/g;

             return qq{call $file$args};
        }


        if ($out =~ /^\s*(?:\.\s+)?(getvalue|getvariable|rpvalue|chi2)(\.pl)?\b(.*)$/i) {
            my ($helper, $suffix, $args) = ($1, $2, $3);
            my $helper_name = defined $suffix ? $helper . $suffix : $helper;
            return "call mctest_helper.bat $helper_name$args";
        }

        $out =~ s/^\s*([A-Za-z_]\w*)=(.*)$/set "$1=$2"/;

         # Bash positional arguments.
        $out =~ s/\$\*/%*/g;

        # Convert common relative/path arguments, but leave non-path expressions
# such as 0.0863x(1.13+0.06)/(-6) unchanged.
$out =~ s{
    (?<![\w.])
    (
        (?:\.\.?/|[A-Za-z0-9_.-]+/)
        [A-Za-z0-9_.-]+
        (?:/[A-Za-z0-9_.-]+)*
    )
}{
    my $path = $1;
    $path =~ tr{/}{\\};
    $path;
}gex;


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

# cmd.exe does not treat single quotes as argument quotes.
$out =~ s/'([^']*)'/"$1"/g;

        if (is_batch_command($out)) {
            $out = "call $out";
        }

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
    my ($script, $target, $marker, $command_info) = @_;
    my ($fh, $tmp) = tempfile('mctestXXXXXX',
        SUFFIX => ($target eq 'lin' ? '.sh' : '.bat'),
        UNLINK => 1);

    print {$fh} $script
        or die "mctest: cannot write temporary script: $!\n";
    close($fh)
        or die "mctest: cannot close temporary script: $!\n";

    if ($target eq 'lin') {
    chmod 0700, $tmp;
    local $ENV{MCTEST_MARKER} = $marker;
    system('bash', $tmp);

    my $raw_status = $?;

    if ($raw_status == -1) {
        return 127;
    }

    my $signal = $raw_status & 127;
    my $exit_code = $raw_status >> 8;

    # Ctrl-C may appear either as SIGINT or as Bash exit status 130.
    if ($signal == 2 || (!$signal && $exit_code == 130)) {
        my $id;

        if (defined $marker && open my $mf, '<', $marker) {
            $id = <$mf>;
            close $mf;
            chomp $id if defined $id;
        }

        my $info = defined $id ? $command_info->{$id} : undef;

        if ($info) {
            print STDERR "mctest: stopped at line $info->{line} "
                       . "of file $info->{file} in command: "
                       . "$info->{command}\n";
        } else {
            print STDERR "mctest: stopped by Ctrl-C; "
                       . "no command marker was found\n";
        }

        return 130;
    }

    return $signal ? 128 + $signal : $exit_code;
}
 else {
#        system($ENV{COMSPEC} || 'cmd.exe', '/c', $tmp);
        require Win32::Process;

        # Restore the previous value automatically when this block exits.
        local $ENV{MCTEST_MARKER} = $marker;

        my $comspec = $ENV{COMSPEC} || 'cmd.exe';
        my $command = qq{"$comspec" /d /c call "$tmp"};
        my $process;

        Win32::Process::Create(
            $process,
            $comspec,
            $command,
            0,          # do not inherit handles
            0x00000010, # CREATE_NEW_CONSOLE
            '.'
        ) or die "mctest: cannot start cmd.exe: $!\n";

        $process->Wait(Win32::Process::INFINITE());

        my $exit_code;
        $process->GetExitCode($exit_code)
            or die "mctest: cannot get child exit code\n";

        # This is the commonly reported Windows Ctrl-C status.
        if ($exit_code == 3221225786) {
            my $id;

            if (defined $marker && open my $mf, '<', $marker) {
                $id = <$mf>;
                close $mf;
                chomp $id if defined $id;
            }

            my $info = defined $id ? $command_info->{$id} : undef;

            if ($info) {
                print STDERR "mctest: stopped at line $info->{line} "
                           . "of file $info->{file} in command: "
                           . "$info->{command}\n";
            } else {
                print STDERR "mctest: stopped by Ctrl-C; "
                           . "no command marker was found\n";
            }

            return 130;
        }


if ($exit_code != 0) {
    my $id;

    if (defined $marker && open my $mf, '<', $marker) {
        $id = <$mf>;
        close $mf;
        chomp $id if defined $id;
    }

    my $info = defined $id ? $command_info->{$id} : undef;

    if ($info) {
        print STDERR "mctest: failed at line $info->{line} "
                   . "of file $info->{file} in command: "
                   . "$info->{command}\n";
    } else {
        print STDERR "mctest: Windows script failed with exit code "
                   . "$exit_code\n";
    }
}

        return $exit_code;
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
