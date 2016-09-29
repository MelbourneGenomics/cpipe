#!/usr/bin/env bash

    CONFIG_FILE=$(readlink -f $(dirname $BASH_SOURCE)/../config.groovy)

function load_config {
    CONFIG=`sed 's/\/\/.*$//' $CONFIG_FILE`
    eval "$CONFIG"
}

function set_config_variable {
    NAME="$1"
    VALUE="$2"
    cp "$CONFIG_FILE" "$CONFIG_FILE.tmp"
    sed 's,'^[\s]*$NAME'=\("\?\).*$,'$NAME'=\1'$VALUE'\1,g' ${CONFIG_FILE}.tmp > "$CONFIG_FILE"
    rm "${CONFIG_FILE}.tmp"
    load_config
}