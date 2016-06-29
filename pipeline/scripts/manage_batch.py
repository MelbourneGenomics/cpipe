#!/usr/bin/env python
'''
###########################################################################
#
# This file is part of Cpipe.
#
# Cpipe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, under version 3 of the License, subject
# to additional terms compatible with the GNU General Public License version 3,
# specified in the LICENSE file that is part of the Cpipe distribution.
#
# Cpipe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################
'''

import argparse
import collections
import datetime
import glob
import os
import os.path
import sys

DEFAULT_PRIORITY = '1'

def write_log(log, msg):
    '''
        write timestamped log message
    '''
    log.write('%s: %s\n' % (datetime.datetime.now().strftime('%y%m%d-%H%M%S'), msg))

def run_command(cmd, log):
    '''
        execute a command on the shell
    '''
    write_log(log, "executing: {0}...".format(cmd))
    os.system(cmd)
    write_log(log, "executing: {0}: done".format(cmd))

def show_batches(out):
    '''
        look for directories in the batches directory
    '''
    for batch in glob.glob('./batches/*'):
        if os.path.isdir(batch):
            out.write('{0}\n'.format(os.path.basename(batch)))

def add_batch(batch_name, profile_name, exome_name, data_files, force, log):
    '''
        create a new batch
    '''
    batch_directory = './batches/{0}'.format(batch_name)
    # check running from correct directory
    if not os.path.isdir('./batches'):
        write_log(log, "ERROR: batches directory not found. Are you running this in the Cpipe root folder?")
        sys.exit(1)

    # make batch directory
    if not os.path.isdir(batch_directory):
        write_log(log, "Creating directory: {0}...")
        os.mkdir(batch_directory)
        write_log(log, "Creating directory: {0}: done")
    else:
        write_log(log, "Batch directory {0} exists")
        if os.path.exists(os.path.join(batch_directory, "samples.txt")):
            if force:
                write_log(log, "WARNING: samples.txt will be overwritten")
            else:
                write_log(log, "ERROR: samples.txt exists. Remove or specify --force to overwrite")
                sys.exit(1)

    # ensure profile exists
    profile_directory = './designs/{0}'.format(profile_name)
    profile = './designs/{0}/{0}.genes.txt'.format(profile_name)
    profile_bed = './designs/{0}/{0}.bed'.format(profile_name)
    if not os.path.isfile(profile):
        write_log(log, 'ERROR: profile "{0}" not found. Use manage_genelists.py to list or add new profiles'.format(profile_name))
        sys.exit(1)

    # configure exome target if specified
    # create_batch.sh behaviour:
    # - if no exome specified, look for target.bed
    # - look in designs for specified exome
    # - if not there, make one from the design gene list and write to initial-design.bed
    # - write exome_target to target_regions.txt

    if exome_name is None:
        exome_name = '{0}.bed'.format(profile_name)

    # look in designs
    if os.path.exists(os.path.join(profile_directory, profile_name)):
        target_region = os.path.abspath(os.path.join(profile_directory, profile_name))
    # look at abspath
    elif os.path.exists(exome_name):
        target_region = os.path.abspath(exome_name)
    # generate the target regions
    else:
        target_region = os.path.join(profile_directory, "generated_target_regions.bed")
        write_log(log, 'Generating target region from design details and writing to {0}...'.format(target_region))
        run_command("python ./pipeline/scripts/combine_target_regions.py --genefiles {0} --bedfiles {1} --exons {2} > {3}".format(profile, profile_bed, "./designs/genelists/exons.bed", target_region), log)
        write_log(log, 'Generating target region from design details and writing to {0}: done'.format(target_region))
        target_region = os.path.abspath(target_region)

    # write this to the batch
    target_file = os.path.join(batch_directory, "target_regions.txt")
    write_log(log, 'writing {0}...'.format(target_file))
    with open(target_file, 'w') as target_fh:
        target_fh.write('EXOME_TARGET="{0}"\n'.format(target_region))
    write_log(log, 'writing {0}: done'.format(target_file))

    # generate samples.txt
    target_file = os.path.join(batch_directory, "samples.txt")
    if data_files is None:
        data_files = "data/*.fastq.gz"
    else:
        data_files = ' '.join(data_files)
    write_log(log, 'writing {0}...'.format(target_file))
    # for now we outsource sample file creation to the groovy
    run_command("(eval `sed \"s/\\/\\/.*//\" ./pipeline/config.groovy`; cd ./batches/{0}; $GROOVY -cp \"$BASE/tools/groovy-hts-sample-info/v1.1/groovy-hts-sample-info.jar:$BASE/tools/groovy-ngs-utils/1.0.2/groovy-ngs-utils.jar\" $BASE/pipeline/scripts/files_to_sample_info.groovy -batch {0} -disease {1} {2} > samples.txt; cd -)".format(batch_name, profile_name, data_files), log)
    write_log(log, 'writing {0}: done'.format(target_file))

    analysis_directory = os.path.join(batch_directory, "analysis")
    if not os.path.isdir(analysis_directory):
        write_log(log, "Creating directory: {0}...".format(analysis_directory))
        os.mkdir(analysis_directory)
        write_log(log, "Creating directory: {0}: done".format(analysis_directory))
    write_log(log, "To run an analysis: cd {0}; ../../../bpipe run ../../../pipeline/pipeline.groovy ../samples.txt".format(analysis_directory))

def add_sample(batch_name, profile_name, data_files, log):
    '''
        add to an existing sample
    '''
    batch_directory = './batches/{0}'.format(batch_name)
    # generate samples.txt
    target_file = os.path.join(batch_directory, "samples.txt")
    if data_files is None:
        data_files = "data/*.fastq.gz"
    else:
        data_files = ' '.join(data_files)

    write_log(log, 'writing {0}...'.format(target_file))
    run_command("(eval `sed \"s/\\/\\/.*//\" ./pipeline/config.groovy`; cd ./batches/{0}; $GROOVY -cp \"$BASE/tools/groovy-hts-sample-info/v1.1/groovy-hts-sample-info.jar:$BASE/tools/groovy-ngs-utils/1.0.2/groovy-ngs-utils.jar\" $BASE/pipeline/scripts/files_to_sample_info.groovy -batch {0} -noheader -disease {1} {2} >> samples.txt; cd -)".format(batch_name, profile_name, data_files), log)
    write_log(log, 'writing {0}: done'.format(target_file))

def show_batch(batch_name, out):
    '''
        show details of a single batch
    '''
    # batch exists?
    batch_directory = './batches/{0}'.format(batch_name)
    if not os.path.isdir(batch_directory):
        out.write('Batch {0}: not found\n'.format(batch_name))
        sys.exit(0)
    # samples.txt
    samples = os.path.join(batch_directory, 'samples.txt')
    if os.path.exists(samples):
        lines = open(samples, 'r').readlines()
        out.write('BATCH {0} contains {1} sample(s)\n'.format(batch_name, len(lines)-1))
        header = lines[0].strip('\n').split('\t')
        justify = max([len(x) for x in header])
        aggregate = collections.defaultdict(set)
        for line in lines[1:]:
            fields = line.strip('\n').split('\t')
            for idx, field in enumerate(fields):
                if len(field) > 0:
                    aggregate[header[idx]].add(field)
        # display aggregate of sample:
        for header in sorted(aggregate.keys()):
            out.write('{0}: {1}\n'.format(header.rjust(justify), ' '.join(aggregate[header])))
    else:
        out.write('samples.txt: not found\n'.format(batch_name))
    # could also show data, target regions

def main():
    '''
        parse command line
    '''
    parser = argparse.ArgumentParser(description='Manage batches')
    parser.add_argument('command', help='command to execute', choices=['add_batch', 'show_batches', 'show_batch', 'add_sample'])
    parser.add_argument('--batch', required=False, help='batch name')
    parser.add_argument('--profile', required=False, help='analysis profile')
    parser.add_argument('--exome', required=False, help='target regions')
    parser.add_argument('--data', required=False, nargs='*', help='fastq files (relative to batch directory)')
    parser.add_argument('--force', action="store_true", required=False, help='force action')
    args = parser.parse_args()

    if args.command == 'show_batches': # list all batches
        show_batches(out=sys.stdout)
    else:
        if not args.batch:
            parser.print_help()
            sys.exit(1)
        if args.command == 'show_batch': # show details of a batch
            show_batch(args.batch, out=sys.stdout)
        elif args.command == 'add_batch': # add a new batch
            if not args.profile:
                write_log(sys.stderr, "ERROR: please provide a profile")
                parser.print_help()
                sys.exit(1)
            add_batch(args.batch, args.profile, args.exome, args.data, args.force, log=sys.stderr)
        elif args.command == 'add_sample': # add additional samples to batch
            add_sample(args.batch, args.profile, args.data, log=sys.stderr)

if __name__ == '__main__':
    main()
