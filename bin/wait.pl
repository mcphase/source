#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

print "sleeping ".$ARGV[0]."s\n";

sleep $ARGV[0];



