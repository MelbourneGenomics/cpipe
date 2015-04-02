#!/bin/sh

# Remove some (boring) predictions
grep -v DOWNSTREAM \
	| grep -v UPSTREAM \
	| grep -v INTRON \
	| grep -v UTR_5_PRIME \
	| grep -v UTR_3_PRIME \
	| grep -v INTERGENIC \
	| grep -v "	SYNONYMOUS_CODING" \
	| grep -v WITHIN_NON_CODING_GENE 
