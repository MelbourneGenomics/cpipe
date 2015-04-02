#!/usr/bin/perl
#-------------------------------------------------------------------------------
#
# Plot a QQ plot (using R) 
# Data is feed as a 1 column of numbers 
#
# Note: Any line that does not match a numeric regular expression, is filtered out).
#
#														Pablo Cingolani
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Parse command line option (file base name)
$base = 'QQ-plot';
if( $ARGV[0] ne '' )	{ $base = $ARGV[0]; }

$pngFile = "$base.png";
$txtFile = "$base.txt";

# Read STDIN and create an R vector
open TXT, "> $txtFile" or die "Cannot open output file '$txtFile'\n";
print TXT "x\n";
for( $ln = 0 ; $l = <STDIN> ; ) {
	chomp $l;

	# Does the string contain exactly one number? (can be float)
	if( $l =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ ) { print TXT "$l\n"; }
}
close TXT;

#---
# Create an R program, save QQ-plot as PNG image
#---

open R, "| R --vanilla --slave " or die "Cannot open R program\n";
print R <<EOF;

qqplot <- function( x, title ) {
	keep <- (x > 0) & (x <= 1) & ( ! is.na(x) );
	x <- x[keep]
	s <- sort(x);
	ly <- -log10(s);

	n <- length(s);
	lx <- -log10( (1:n) / (n+1) )
	
	# Show auto range
	#par( mfrow=c(2,1) );
	#plot( lx, ly, main=title, xlab="-Log[ rank / (N+1) ]", ylab="-Log[ p ]" );
	#abline( 0 , 1 , col='red');

	# Show full range in both plots
	range <- c(0 , max(lx, ly) );
	plot( lx, ly, xlim=range, ylim=range, main=title, xlab="-Log[ rank / (N+1) ]", ylab="-Log[ p ]" );
	abline( 0 , 1 , col='red');
}

png('$pngFile', width = 1024, height = 1024);

data <- read.csv("$txtFile", sep='\t', header = TRUE);
qqplot( data\$x, "$base" );

dev.off();
quit( save='no' )
EOF

close R;

#---
# Show figure
#---

$os = `uname`;
$show = "eog"; 
if( $os =~ "Darwin" )	{ $show = "open"; }
`$show $pngFile &`;

