#!/bin/sh -e

queueFile=queue_build_ncbi.txt

#mkdir download/ncbi
cd download/ncbi

#---
# Download from NCBI (Bacterial genoms)
#---

# wget -r -np -nc -A "gbk,html" "http://ftp.ncbi.nih.gov/genomes/Bacteria/"

#---
# Create directory structure
#---

# # Move all downloaded file to this directory
# mv ftp.ncbi.nih.gov/genomes/Bacteria/* .
# rmdir ftp.ncbi.nih.gov/genomes/Bacteria
# rmdir ftp.ncbi.nih.gov/genomes
# rmdir ftp.ncbi.nih.gov

# # Remove old queu file
# touch $queueFile
# rm $queueFile
# 
# for dir in `find . -mindepth 1 -maxdepth 1 -type d `
# do
# 	# Config file entries
# 	gen=`basename $dir`
# 	echo -e "$gen.genome : $gen\n" | tee -a ncbi_append.snpEff.config
# 
# 	# Collapse all fine into one 
# 	cd $dir
# 	cat *.gbk > genes.gb
# 	cd - > /dev/null
# 
# 	# Build queue file
# 	echo "java -Xmx4G -jar snpEff.jar build -v -genbank $gen" >> $queueFile
# done

#---
# Move to data directory
#---

echo Move to DATA dir
for dir in `find . -mindepth 1 -maxdepth 1 -type d `
do
	echo "    $dir"
	mv $dir ../../data/
done

cd -
