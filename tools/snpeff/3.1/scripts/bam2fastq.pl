#!/usr/bin/perl

while( $l = <STDIN> ) {
	chomp $l;
	@t = split /\t/, $l;
	print "\@$t[0]\n$t[9]\n+\n$t[10]\n";
}

