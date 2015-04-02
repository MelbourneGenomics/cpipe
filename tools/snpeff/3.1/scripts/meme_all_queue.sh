#!/bin/sh

meme=$HOME/tools/meme_4.8.1/bin/meme

for fa in *.fa
do

	dir=`basename $fa .fa`
	size=`cat $fa | wc -c`
	seqs=`cat $fa | grep -v "^>" | wc -l`

	echo $meme $fa -dna  -oc $dir  -mod zoops  -nmotifs 1  -minw 6 -maxw 20 
done
