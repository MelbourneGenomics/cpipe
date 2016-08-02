#!/usr/bin/env bash

# Installs all dependencies for the installer and pipeline on Ubuntu

#Fix docker issues
sudo mkdir /etc/systemd/system/docker.service.d/
'[Service]
ExecStart=
ExecStart=/usr/bin/docker daemon -H fd:// -s overlay' >> /etc/systemd/system/docker.service.d/overlay.conf
sudo systemctl daemon-reload

# Install gradle repos
sudo add-apt-repository ppa:cwchien/gradle
sudo apt-get update

# Install apt-getable things
sudo apt-get install -y git make gcc poppler-utils zlib1g-dev ncurses-dev g++ patch libssl-dev unzip maven\
 gradle libcurl4-openssl-dev texinfo openjdk-8-jdk
