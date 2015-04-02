#!/usr/bin/perl

#------------------------------------------------------------------------------
# Split a fasta file (create one file per sequence)
#
#
#------------------------------------------------------------------------------

use strict;

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

my($seq, $name) = ('', '');
my($lineNum, $l, $newName);
#---
# Read fasta file
#---
for($lineNum=0 ; $l = <STDIN> ; $lineNum++ ) {
	chomp $l;
	if( $l =~/^>\s*(.*)\s*$/ ) {
		$newName = $1;
		if( $seq ne "" ) { print "$name\t$seq\n"; } 
		# New sequence
		$name = $newName;
		$seq = "";
	} else { $seq .= $l; }
}

if( $seq ne "" ) { print "$name\t$seq\n"; }

