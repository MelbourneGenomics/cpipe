#!/usr/bin/env bash
CONFIG_FILE=$(readlink -f $(dirname $BASH_SOURCE)/../config.groovy)
CONFIG=`sed 's/\/\/.*$//' $CONFIG_FILE`
eval "$CONFIG"
