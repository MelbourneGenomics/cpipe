#!/bin/sh

# Keep interesting predictions
java -jar $HOME/tools/VcfEtc.jar filter \
       "( EFF =~ 'NON_SYN' ) | ( EFF =~ 'CODON') | ( EFF =~ 'SPLICE') | ( EFF =~ 'STOP') | ( EFF =~ 'START') | ( EFF =~ 'FRAME') | ( EFF =~ 'LOST') | ( EFF =~ 'DELETED' )" \

