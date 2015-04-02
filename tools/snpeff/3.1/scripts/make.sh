#!/bin/sh

VERSION_SNPEFF=3.1
VERSION_SNPSIFT=1.7
VERSION_SNPSQL=0.2

#---
# Build SnpEff
#---

cd $HOME/workspace/SnpEff/

mvn assembly:assembly
cp target/snpEff-$VERSION_SNPEFF-jar-with-dependencies.jar $HOME/snpEff/snpEff.jar


# Install JAR file in local Maven repo
mvn install:install-file \
	-Dfile=target/snpEff-$VERSION_SNPEFF.jar \
	-DgroupId=ca.mcgill.mcb.pcingola \
	-DartifactId=snpEff \
	-Dversion=$VERSION_SNPEFF \
	-Dpackaging=jar \
	-DgeneratePom=true

cd - 

#---
# Build SnpSift
#---
cd $HOME/workspace/SnpSift/

mvn assembly:assembly
cp target/snpSift-$VERSION_SNPSIFT-jar-with-dependencies.jar $HOME/snpEff/SnpSift.jar

# Install JAR file in local Maven repo
mvn install:install-file \
	-Dfile=target/snpSift-$VERSION_SNPSIFT.jar \
	-DgroupId=ca.mcgill.mcb.pcingola \
	-DartifactId=snpSift \
	-Dversion=$VERSION_SNPSIFT \
	-Dpackaging=jar \
	-DgeneratePom=true

cd - 

#---
# Build SnpSql
#---
cd $HOME/workspace/SnpSql/

mvn assembly:assembly
cp target/snpSql-$VERSION_SNPSQL-jar-with-dependencies.jar $HOME/snpEff/SnpSql.jar

# Install JAR file in local Maven repo
mvn install:install-file \
	-Dfile=target/snpSql-$VERSION_SNPSQL.jar \
	-DgroupId=ca.mcgill.mcb.pcingola \
	-DartifactId=snpSql \
	-Dversion=$VERSION_SNPSQL \
	-Dpackaging=jar \
	-DgeneratePom=true

cd - 
