#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

print "sleeping ".$ARGV[0]."s\n";
$ARGV[0]=~s/x/*/g;sleep eval $ARGV[0];



