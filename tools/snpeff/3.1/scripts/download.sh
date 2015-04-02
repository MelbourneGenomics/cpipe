#!/bin/sh -e

RELEASE=68

#mkdir download
cd download

#---
# Download from ENSEMBL
#---

# # Download GTF files (annotations)
# wget -r -A "*gtf.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/gtf/"
# 
# # Download FASTA files (reference genomes)
# wget -r -A "*toplevel.fa.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"
# 
# # Download CDS sequences
# wget -r -A "*cdna.all.fa.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"
# 
# # Download PROTEIN sequences
# wget -r -A "*.pep.all.fa.gz" "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"


#---
# Download from NCBI (Bacterial genoms)
#---

# wget -r -np -nc -A "gbk,html" "http://ftp.ncbi.nih.gov/genomes/Bacteria/"

#---
# Create directory structure
#---

# # Move all downloaded file to this directory
# mv `find ftp.ensembl.org -type f` .
# 
# # Gene annotations files
# for gtf in *.gtf.gz
# do
# 	short=`../scripts/file2GenomeName.pl $gtf | cut -f 5`
# 	echo ANNOTATIONS: $short
# 
# 	mkdir -p data/$short
# 	cp $gtf data/$short/genes.gtf.gz
# done
#  
# # Reference genomes files
# mkdir -p data/genomes
# for fasta in *.dna.toplevel.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 5`
# 	echo REFERENCE: $genome
# 
# 	cp $fasta data/genomes/$genome.fa.gz
# done
# 
# # CDS genomes files
# for fasta in *.cdna.all.fa.gz
# do
# 	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 5`
# 	echo CDS: $genome
# 
# 	cp $fasta data/$genome/cds.fa.gz
# done
# 
# # Protein seuqence files
# for pep in *.pep.all.fa.gz
# do
# 	short=`../scripts/file2GenomeName.pl $pep | cut -f 5`
# 	echo PROTEIN: $short
# 
# 	mkdir -p data/$short
# 	cp $pep data/$short/protein.fa.gz
# done

#---
# Config file entries
#---

for fasta in *.cdna.all.fa.gz
do
	genome=`../scripts/file2GenomeName.pl $fasta | cut -f 4`
	short=`../scripts/file2GenomeName.pl $fasta | cut -f 5`

	# Individual genome entry
	echo -e "$short.genome : $genome"
	echo -e "$short.reference : ftp://ftp.ensembl.org/pub/release-$RELEASE/gtf/"
	echo
done

# Back to parent dir
cd - > /dev/null

#---
# Create build queue entries
#---

rm -vf queue_build.txt

# Build from refseq files
echo "java -Xmx10G -jar snpEff.jar build -v -refseq hg19"

# Build from TXT files
for genes in `ls data/*/genes.txt* | grep -v hg19`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -txt $genomeName"
done | sort >> queue_build.txt

# Build from GFF2 files
echo "java -Xmx10G -jar snpEff.jar build -v -gff2 amel2" >> queue_build.txt

# Build from GFF3 files
for genes in `ls data/*/genes.gff* | grep -v amel2`
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -gff3 $genomeName"
done | sort >> queue_build.txt

# Build from GTF22 files
for genes in data/*/genes.gtf*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -gtf22 $genomeName"
done | sort >> queue_build.txt

# Build from GenBank files
for genes in data/*/genes.gb*
do
	dir=`dirname $genes`
	genomeName=`basename $dir`
	echo "java -Xmx10G -jar snpEff.jar build -v -genbank $genomeName"
done | sort >> queue_build.txt
