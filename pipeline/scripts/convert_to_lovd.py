#!/usr/bin/env python
'''
===========================================================================
  convert vcf table output by GATK to lovd format

  Usage:
    python convert_to_lovd.py --vcf source_vcf < table.tsv > table.lovd
===========================================================================
'''

import argparse
import collections
import datetime
import sys

# what to consider an unknown or empty value
EMPTY_BEFORE = ('', 'NA', 'unknown')

# what to replace empty fields with 
EMPTY_AFTER = ''

# max number of warnings to show
MAX_WARNINGS = 100

def write_log_factory(params):
    '''
      write a datestamped log message
    '''
    def write_log_fn(fh, msg):
        if msg.startswith('WARNING'):
            params[0] -= 1
        if params[0] >= 0 or not msg.startswith('WARNING'):
            now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
            fh.write('{0}: {1}\n'.format(now, msg))
        else:
            if params[0] % 100 == 0:
                fh.write('.')
    return write_log_fn

write_log = write_log_factory([MAX_WARNINGS])

def get_ann_fields(vcf, log):
    '''
        extract ANN fields from vcf header
    '''
    write_log(log, 'parsing vcf...')
    result = None
    for line in vcf:
        if line.startswith('##INFO=<ID=ANN'):
           info = line.split('Format: ')
           if len(info) == 2:
               description = info[1].split('"')[0]
               result = description.split('|')
               write_log(log, 'parsing vcf: success: {0} fields: {1}'.format(len(result), result))
               return result
    write_log(log, 'WARNING: parsing vcf: ANN INFO field not found. Variants have not been annotated.')

def process_table(source, target, ann_headers, log):
    '''
        write tab separated columns, expanding ANN fields and writing one transcript per line
    '''
    write_log(log, 'processing table...')
    count = 0
    warnings = 0
    transcript_counts = collections.defaultdict(int)
    for count, line in enumerate(source):
        fields = line.strip('\n').split('\t')
        if count == 0: # header
            header = fields
            ann = fields.index('ANN')
            output = fields[:ann] + ann_headers + fields[ann + 1:]
            target.write('{0}\n'.format('\t'.join(output)))
        else:
            transcripts = fields[ann].split(',') 
            transcript_counts[len(transcripts)] += 1
            for transcript in transcripts:
                ann_fields = transcript.split('|') 
                if len(ann_fields) != len(ann_headers): # ann column count mismatch
                    error = []
                    common = min(len(ann_fields), len(ann_headers))
                    for idx in range(0, common):
                        if ann_fields[idx] in EMPTY_BEFORE:
                            ann_fields[idx] = EMPTY_AFTER
                        error.append("{0}: '{1}'".format(ann_headers[idx], ann_fields[idx]))
                    write_log(log, 'WARNING: line {0}: {1} ANN fields, {2} expected: {3}'.format(count, len(ann_fields), len(ann_headers), ', '.join(error)))
                    warnings += 1
                # if ann_fields is shorter than expected, extend with empty
                while len(ann_fields) < len(ann_headers):
                    ann_fields.append(EMPTY_AFTER)
                output = []
                # map empty to unknown
                for col in fields[:ann] + ann_fields + fields[ann + 1:]:
                    if col in EMPTY_BEFORE:
                        col = EMPTY_AFTER
                    output.append(col)
                target.write('{0}\n'.format('\t'.join(output)))

    write_log(log, 'processing table: finished {0} records. {1} warnings. transcript counts: {2}'.format(count, warnings, transcript_counts))

def main():
    '''
        run via command line
    '''
    parser = argparse.ArgumentParser(description='Convert table to lovd')
    parser.add_argument('--vcf', required=True, help='original vcf')
    args = parser.parse_args() 
    # get ANN fields
    ann_fields = get_ann_fields(open(args.vcf, 'r'), sys.stderr)
    # process table
    process_table(sys.stdin, sys.stdout, ann_fields, sys.stderr)

if __name__ == '__main__':
    main()
