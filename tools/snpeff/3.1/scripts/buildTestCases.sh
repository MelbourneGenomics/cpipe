#!/bin/sh

SNPEFF="java -Xmx2G -jar snpEff.jar"

# Test cases hg37
$SNPEFF build -noLog -txt testCase

# Test cases hg37.61
$SNPEFF build -noLog -gtf22 testHg3761Chr15
$SNPEFF build -noLog -gtf22 testHg3761Chr16 

# Test cases hg37.63
$SNPEFF build -noLog -gtf22 testHg3763Chr1
$SNPEFF build -noLog -gtf22 testHg3763Chr20 
$SNPEFF build -noLog -gtf22 testHg3763ChrY 

# Test cases hg37.65
$SNPEFF build -noLog -gtf22 testHg3765Chr22

# Test cases hg37.66
$SNPEFF build -noLog -gtf22 testHg3766Chr1

# Test cases hg37.67
$SNPEFF build -noLog -gtf22 testHg3767Chr21Mt
