#!/bin/bash
# Start a new subshell, set variables in that, and then give control back to the user. They can then quit using Ctrl+D
ROOT=`dirname $0`
bash --init-file $ROOT/_env
