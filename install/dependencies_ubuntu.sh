#!/usr/bin/env bash

# Installs all dependencies for the installer and pipeline on Ubuntu

#Fix docker issues
sudo -s
mkdir /etc/systemd/system/docker.service.d/
echo '[Service]
ExecStart=
ExecStart=/usr/bin/docker daemon -H fd:// -s overlay' >> /etc/systemd/system/docker.service.d/overlay.conf
systemctl daemon-reload

# Install gradle repos
add-apt-repository ppa:cwchien/gradle
apt-get update

# Install apt-getable things
apt-get install -y git make gcc poppler-utils zlib1g-dev ncurses-dev g++ patch libssl-dev unzip maven\
 gradle libcurl4-openssl-dev texinfo openjdk-8-jdk
