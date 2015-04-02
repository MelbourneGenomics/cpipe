#!/bin/sh

for testNum in 1 2 3 4 
do
	./scripts/snpEffM.sh test $testNum test.1M.vcf > test.1M.out.vcf
	echo
done

