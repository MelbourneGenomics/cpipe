#!/usr/bin/env bash

# Installs all dependencies for the installer and pipeline
# Note: Only for Ubuntu 16.04

JAVA_VERSON=8u91-b14-0ubuntu4~16.04.1

# Install python using pyenv
curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash
export PATH="/home/ubuntu/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)" >> ~/.bash_profile
pyenv install $PYTHON_VERSION

# Install perl using perlbrew
curl -L https://install.perlbrew.pl | bash
echo "source ~/perl5/perlbrew/etc/bashrc" >> ~/.bash_profile
perlbrew install perl-$PERL_VERSION

#Fix docker issues
sudo mkdir /etc/systemd/system/docker.service.d/
'[Service]
ExecStart=
ExecStart=/usr/bin/docker daemon -H fd:// -s overlay' >> /etc/systemd/system/docker.service.d/overlay.conf
sudo systemctl daemon-reload

#Install R using CRAN repository
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo add-apt-repository 'deb http://cran.ms.unimelb.edu.au/bin/linux/ubuntu xenial/'
sudo apt-get update
sudo apt-get install -y r-base=$R_VERSION\*

# Install Java
sudo apt-get install -y openjdk-8-jre=$JAVA_VERSION\*

# Install apt-getable things
sudo apt-get install -y git make gcc poppler-utils zlib1g-dev ncurses-dev g++ patch libssl-dev unzip maven
