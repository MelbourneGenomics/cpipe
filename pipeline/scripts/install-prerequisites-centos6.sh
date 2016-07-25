#!/bin/bash
############################################################
#
# Installs Centos 6 dependencies
#
# this is a set of recommended commands to run in preparation
# for Cpipe.
# evaluate for yourself if you want to run them on your system.
############################################################

# software
sudo yum update
sudo yum install git make gcc poppler-utils
sudo yum install java-1.8.0-openjdk zlib-devel ncurses-devel gcc-c++ R patch
sudo yum install openssl-devel mysql
sudo yum install perl-devel

# java home
echo "export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk.x86_64/bin/java" >> ~/.bash_profile
. ~/.bash_profile

# perl too old?
curl -L https://install.perlbrew.pl | bash
echo "source ~/perl5/perlbrew/etc/bashrc" >> ~/.bash_profile
. ~/.bash_profile
perlbrew install-patchperl
perlbrew install perl-5.18.1
perlbrew switch perl-5.18.1

# vep and perl
#sudo yum install perl-core perl-DBD-MySQL

# don't need sudo if using perlbrew
cpan "install HTTP::Tiny"
perl -MCPAN -e "install DBI"
perl -MCPAN -e "install LWP::Simple"
perl -MCPAN -e "install LWP::Protocol::https"
cpan Archive::Zip
