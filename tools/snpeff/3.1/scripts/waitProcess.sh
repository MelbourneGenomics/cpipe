#!/bin/sh

PID=$1
while ps -p $PID > /dev/null ; 
do
	sleep 1; 
done
