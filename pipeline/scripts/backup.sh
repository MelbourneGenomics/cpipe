#!/bin/bash
# vim: ts=4:expandtab:cindent:sw=4
####################################################################
#
# Melbourne Genomics Pipeline Backup Script
#
# This script executes an rsync command designed to run nightly to
# backup the Melbourne Genomics production pipeline, batches of data
# and the root Git repository that stores the source code.
#
# Notes:
#
#  -  The target address is a backup server in WEHI.
#
#  -  The backup uses rsync with public key authentication.
#     the key pair is located in the base directory/.ssh folder.
#  
#  -  The backup excludes the "external" folder as that is expected
#     to carry data of users external to Melbourne Genomics for whom
#     MG is not responsible (nor authorized) to backup data for.
#
# Author: Simon Sadedin, simon.sadedin@mcri.edu.au
# Date  : 7th March 2014
#
####################################################################
#
# NOTE: these SLURM settings are included here for reference, but
#       backup using a node does not work because the node ip addresses
#       are not whitelisted on the WEHI end.
#
#SBATCH --job-name=mgha_nightly_backup
#SBATCH --account VR0002
#SBATCH --mem=16384
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH -p terri

# Base directory to which backup should be mirrored
BASE=/vlsci/VR0320/shared

# IP Address of target server
TARGET_IP_ADRESS=128.250.252.218
TARGET_USER=vlsci-syncuser

# Backup failure email recipients
EMAILS="bakker.t@wehi.edu.au ssadedin@gmail.com"

# The source folders that will be backed up
SOURCES="repo production"

# Log file to save output to
LOGFILE="$BASE"/logs/backup_`date +'%y_%m_%d'`.log

####################################################################
#
# Error handler
#
####################################################################
function err() {
    echo ===============================================
    echo "ERROR: $1"
    echo ===============================================
    echo 
    echo "Exiting at "`date`
    echo

    printf "Error Message : $1\n\nTail of log file:\n\n" > $BASE/logs/error.log 
    tail -n 30 $LOGFILE >> $BASE/logs/error.log
    
    mail -s "WARNING: Melbourne Genomics Backup Failure" $EMAILS  < $BASE/logs/error.log

    exit 1
}

####################################################################
#
# Main Backup Routine
#
####################################################################
function backup() {
    echo ===============================================
    echo "Executing backup "`date`
    echo ===============================================
    echo
    cd $BASE || err "Unable to change directory to base dir: $BASE"

    #echo "rsync -r -e \"ssh -i $BASE/.ssh/id_rsa\" repo production $TARGET_USER@$TARGET_IP_ADRESS:"
    rsync -a -v -r -e "ssh -i $BASE/.ssh/id_rsa" $SOURCES $TARGET_USER@$TARGET_IP_ADRESS: || err "Rsync returned failure exit code"
    echo
    echo "Done at "`date`
}

####################################################################
#
# Main Entry Point
#
####################################################################

# Run the actual backup in subshell to capture output
( backup ) > $LOGFILE 2>&1

