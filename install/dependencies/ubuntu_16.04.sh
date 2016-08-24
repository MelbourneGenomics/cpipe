#!/usr/bin/env bash

# Installs all dependencies for the installer and pipeline on Ubuntu

#Fix for docker apt-get
apt-get update
apt-get install -y software-properties-common

#Fix docker issues
mkdir -p /etc/systemd/system/docker.service.d/
echo '[Service]
ExecStart=
ExecStart=/usr/bin/docker daemon -H fd:// -s overlay' >> /etc/systemd/system/docker.service.d/overlay.conf
systemctl daemon-reload

# Install gradle repos
add-apt-repository -y ppa:cwchien/gradle
apt-get update

#Set ubuntu-specific variables
export JAVA_HOME=/usr

# Install apt-getable things
apt-get install -y git make poppler-utils zlib1g-dev ncurses-dev gcc g++ gfortran patch libssl-dev unzip maven\
 gradle libcurl4-openssl-dev texinfo openjdk-8-jdk python-pip mysql-client xorg-dev libreadline-dev libbz2-dev liblzma-dev\
 libpcre3-dev libsqlite3-dev
