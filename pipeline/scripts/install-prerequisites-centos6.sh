#!/bin/bash
############################################################
#
# Installs Centos 6 dependencies
#
############################################################

# software
sudo yum update
sudo yum install git make gcc poppler-utils
sudo yum install java-1.8.0-openjdk zlib-devel ncurses-devel gcc-c++ R patch

# java home
echo "export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk.x86_64/bin/java" >> ~/.bash_profile
. ~/.bash_profile

# vep and perl
sudo yum install perl-core perl-DBD-MySQL
sudo perl -MCPAN -e "install HTTP::Tiny"
sudo perl -MCPAN -e "install LWP::Simple"

# perl too old?
curl -L https://install.perlbrew.pl | bash
echo "source ~/perl5/perlbrew/etc/bashrc" >> ~/.bash_profile
. ~/.bash_profile
perlbrew install-patchperl
perlbrew install perl-5.18.1
perlbrew switch perl-5.18.1
