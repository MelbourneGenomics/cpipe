#!/usr/bin/env python3
import argparse
from subprocess import check_output
from os import environ
from cpipe import arg_validation
from cpipe.scripts.manage_batch import create_parser
def run(cmd):
    check_output(cmd, shell=True, cwd=environ.get('CPIPE_ROOT'))

def check_java(args):
    # Run the java check if necessary
    if args.check_java:
        run('doit check_java')



# function print_general_usage {
#     echo "${normal}Usage: ./cpipe <CPIPE OPTIONS> COMMAND <COMMAND OPTIONS>"
#     echo "${bold}Commands (type --help after any command for more details):"
#     echo "  ${bold}run${normal}: Runs the analysis pipeline"
#     echo "  ${bold}test${normal}: Runs the pipeline tests"
#     echo "  ${bold}batch${normal}: Creates and modifies analysis batches"
#     echo "  ${bold}genelist${normal}: Creates and modifies genelists"
#     echo "  ${bold}metadata${normal}: Creates and modifies sample metadata files"
#     echo "${bold}Cpipe Options"
#     echo "  ${bold}-b, --batch <batch name>"
#     echo "    ${normal}Specify a batch (a subdirectory inside batches) to use for the run and bpipe commands. Defaults to a batch named 'batch'"
#     echo "  ${bold}--help, --usage"
#     echo "    ${normal}Prints this help page"
#     echo "  ${bold}-j, --no-java-check"
#     echo "    ${normal}Disables the java version check. Only do this if you know what you're doing"
# }
#
# function print_metadata_usage {
#     echo "${normal}Usage: ./cpipe <CPIPE OPTIONS> metadata <SUBCOMMAND>"
#     echo "${bold}metadata subcommands (use --help after each for usage options):"
#     echo "  ${bold}check${normal}: Check an existing metadata file"
#     echo "  ${bold}update${normal}: Update an existing metadata file"
# }
#
# function print_run_usage {
#     echo "${normal}Usage: ./cpipe run <COMMAND OPTIONS>"
#     echo "${bold}Run Options"
#     echo "  ${bold}--help, --usage"
#     echo "    ${normal}Prints this help page"
#     echo "  ${bold}-b, --batch <batch name>"
#     echo "    ${normal}Specify a batch (a subdirectory inside batches) to use for the run and bpipe commands. Defaults to a batch named 'batch'. Setting the batch here is the same as using the flag before the 'run' command."
#     echo "  ${bold}-p, --bpipe-options <options>"
#     echo "    ${normal}Specify options to pass to bpipe. Refer to http://docs.bpipe.org/Commands/run/ for reference"
# }
#
# function check_batch {
#     if [[ -z $BATCH ]]; then
#         >&2 echo 'You must specify a batch using `--batch <batch name>`'
#         exit 1
#     fi
# }

parser = argparse.ArgumentParser(description='Runs core Cpipe functions')
subparsers = parser.add_subparsers(dest='command')

# Run command
run_parser = subparsers.add_parser('run', help='Runs the analysis pipeline')
run_parser.add_argument('batch', type=arg_validation.existing_batch)
run_parser.add_argument('--no-java-check', '-j', type=arg_validation.existing_batch)

# Run command
run_parser = subparsers.add_parser('batch', help='Manage Cpipe batches and metadata files')
run_parser.add_argument('batch', type=arg_validation.existing_batch)

#Parse command line arguments
# JAVA_CHECK=1
# ARGS=$(env POSIXLY_CORRECT=1 getopt -o jb: --long "no-java-check,help,usage,batch:" -n $(basename $BASH_SOURCE) -- "$@")
# eval set -- "$ARGS"
# while true ; do
#     case "$1" in
#         -b|--batch)
#             BATCH=$2
#             shift 2
#         ;;
#         -j|--no-java-check)
#             JAVA_CHECK=0
#             shift 1
#         ;;
#         --help|--usage)
#             print_general_usage
#             exit 0
#         ;;
#         --)
#             shift
#             break
#         ;;
#     esac
# done
#
# case "$1" in
#     batch)
#         #  e.g. docker run cpipe batch add_batch --batch batch_identifier --profile profile_name
#         shift 1
#         cd pipeline/scripts
#         python -m manage_batch "$@" < /dev/stdin
#     ;;
#     genelist)
#         # e.g. docker run cpipe genelist show_bed --profile profile_name
#         shift 1
#         python pipeline/scripts/manage_genelists.py "$@" < /dev/stdin
#     ;;
#     metadata)
#         shift 1
#         case "$1" in
#             --help|--usage)
#                 print_metadata_usage
#                 exit 0
#             ;;
#             check)
#                 shift 1
#                 #e.g. docker run cpipe metadata check < ./batches/batch_identifier/samples.txt
#                 python pipeline/scripts/check_metadata.py "$@" < /dev/stdin
#             ;;
#             update)
#                 shift 1
#                 #e.g. docker run cpipe metadata update --sample sample_name --name prioritised_genes --value "4:ABC1,ABC2" --target ./batches/batch_identifier/samples.txt
#                 python pipeline/scripts/update_metadata.py "$@" < /dev/stdin
#             ;;
#             *)
#                 print_metadata_usage
#                 exit 1
#             ;;
#
#         esac
#     ;;
#     bpipe)
#         check_batch
#         shift 1
#         cd ${ROOT}/batches/${BATCH}/analysis
#         ${BASE}/bpipe $@
#     ;;
#     run)
#         # Parse args
#         shift 1
#         ARGS=$(getopt -o p:b: --long "batch:,bpipe-options:,help,usage" -n $(basename $BASH_SOURCE) -- "$@")
#         eval set -- "$ARGS"
#
#         # Process args - they can specify a batch directory to replace the default 'batch', and they can specify bpipe options manually
#         while true ; do
#             case "$1" in
#                 -b|--batch)
#                     BATCH=$2
#                     shift 2
#                 ;;
#                 -p|--bpipe-options)
#                     BPIPE_OPTIONS=$2
#                     shift 2
#                 ;;
#                 --usage|--help)
#                     print_run_usage
#                     exit 0
#                 ;;
#                 --)
#                     shift
#                     break
#                 ;;
#             esac
#         done
#
#         check_batch
#         check_java
#         mkdir -p batches/${BATCH}/analysis
#         cd ${ROOT}/batches/${BATCH}/analysis
#         ${BASE}/bpipe run ${BPIPE_OPTIONS} ${BASE}/pipeline/pipeline.groovy ${BASE}/batches/${BATCH}/samples.txt < /dev/stdin
#     ;;
#     test)
#         shift 1
#         check_java
#         pipeline/scripts/run_unit_tests.sh && pipeline/scripts/run_tests.sh detect_mutations_test
#     ;;
#     *)
#         echo "Invalid cpipe command!"
#         print_general_usage
#         exit 1
#     ;;
# esac
