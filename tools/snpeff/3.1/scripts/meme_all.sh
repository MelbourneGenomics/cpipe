#!/bin/sh

meme=$HOME/tools/meme_4.8.1/bin/meme

for fa in *.fa
do

	dir=`basename $fa .fa`
	size=`cat $fa | wc -c`
	seqs=`cat $fa | grep -v "^>" | wc -l`

	if [ ! -e $dir/meme.html ]
	then
		echo 
		echo "#---"
		echo "# $fa : $seqs sequences"
		echo "#---"
		echo 

		$meme $fa \
			-dna \
			-oc $dir \
			-mod zoops \
			-nmotifs 2 \
			-minw 5 \
			-maxw 15 \ 
			#-maxsize $size
	fi

done
