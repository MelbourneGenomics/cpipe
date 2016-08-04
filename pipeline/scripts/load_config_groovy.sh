#!/usr/bin/env bash
FILENAME=$1
CONFIG=`sed 's/\/\/.*$//' $FILENAME`
eval "$CONFIG"
