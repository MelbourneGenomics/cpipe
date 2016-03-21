#!/usr/bin/env bash
# parameters
# 1: input
# 2: output
# 3: vep

VARIANTS=`grep -c -v '^#' < $1`
echo "$VARIANTS variant(s) found in $1"
if [ $VARIANTS -eq 0 ];
then
  cp $1 $2
else
  PERL5LIB="$3" perl $3/filter_vep.pl \
    --input_file $1 \
    --filter "Consequence not matches stream" \
    --filter "BIOTYPE match protein_coding" \
    --filter "Feature" \
    --force_overwrite \
    --format vcf \
    -o $2 \
    --only_matched
fi
