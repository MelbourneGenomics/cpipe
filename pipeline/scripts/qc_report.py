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
#
# Purpose:
#   Generate coverage stats
# Usage:
#
# arguments:
# --report_cov: intersected coverage file with genes from bedtools
# --gene_cov: coverage of each gene
# --exome_cov: exome coverage file with genes
# --ontarget: target reads count file
# --metrics: metrics output from Picard
# --study: ID of study
# --meta: meta data file
# --threshold: threshold for satisfactory coverage
# --classes: how to categorise results
# --exome: target regions for entire exome
# --gc: gene categories of genes
# --anonymous: do not show study ID
# --write_karyotype: write karyotype details to this file
# --fragments: file containing fragment details
#
##############################################################################
'''

import collections
import datetime
import re
import sys

MEAN_RANGE = 0.8 # calculate proportion of coverage within this fraction of the mean

def write_log(log, msg):
    '''
        write a date stamped message to log
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))

def calculate_karyotype(exome_cov, log=None):
    '''
        * calculate the mean coverage on chr1 and chr22, chrX, and chrY
              (calculated by number of reads overlapping each position in exome)
        * if chrY coverage < 5 and chrX coverage > 30: sample is female
        * else if chrX coverage / autosome coverage < 0.7: sample is male
        * else sample is other
    '''
    write_log(log, 'calculating karyotype...')
    stats = collections.defaultdict(int)
    for idx, line in enumerate(exome_cov):
        fields = line.strip().split('\t') # chr, start, end, offset, coverage
        if len(fields) > 4:
            target = None
            if fields[0].lower() == 'chrx':
                target = 'x'
            if fields[0].lower() == 'chry':
                target = 'y'
            if fields[0].lower() == 'chr1' or fields[0].lower() == 'chr22':
                target = 'a'
            if target is not None:
                stats['{0}c'.format(target)] += int(fields[4])
                stats['{0}n'.format(target)] += 1
        else:
            write_log(log, 'skipped line {0}: {1}'.format(idx, line.strip()))
        if idx % 100000 == 0:
            write_log(log, 'processed {0} lines...'.format(idx))
    write_log(log, 'processed {0} lines'.format(idx))

    result = {'sex': 'OTHER', 'x_mean_coverage': 0.0, 'y_mean_coverage': 0.0, \
        'autosome_mean_coverage': 0.0}

    if stats['xn'] > 0:
        result['x_mean_coverage'] = 1. * stats['xc'] / stats['xn']
    if stats['yn'] > 0:
        result['y_mean_coverage'] = 1. * stats['yc'] / stats['yn']
    if stats['an'] > 0:
        result['autosome_mean_coverage'] = 1. * stats['ac'] / stats['an']

    # determine gender
    if result['y_mean_coverage'] < 5 and result['x_mean_coverage'] > 30:
        result['sex'] = 'FEMALE'
    elif result['autosome_mean_coverage'] == 0:
        write_log(log, 'WARNING: no autosome coverage')
        result['sex'] = 'OTHER'
    elif result['x_mean_coverage'] / result['autosome_mean_coverage'] < 0.7:
        result['sex'] = 'MALE'

    write_log(log, 'calculating karyotype: done: {0}'.format(result))

    return result

def write_karyotype(target, karyotype, meta):
    '''
        write karyotype details to target
    '''
    target.write('Sex\t{0}\n'.format(meta['sex']))
    target.write('Inferred Sex\t{0}\n'.format(karyotype['sex']))
    target.write('xCoverage\t{0}\n'.format(karyotype['x_mean_coverage']))
    target.write('yCoverage\t{0}\n'.format(karyotype['y_mean_coverage']))
    target.write('autosomeCoverage\t{0}\n'.format(karyotype['autosome_mean_coverage']))

def mean(items):
    '''
        compute the mean of a list of items
    '''
    if len(items) == 0:
        return 0
    else:
        return sum(items) / float(len(items))

def median(items):
    '''
      return the median of a list of items
    '''
    sorted_list = sorted(items)
    if len(items) % 2 == 0:
        high = len(items) / 2
        return (sorted_list[high] + sorted_list[high-1]) / 2.
    else:
        mid = (len(items) - 1) / 2
        return sorted_list[mid]

def parse_metadata(meta, study):
    '''
       reads metadata file and returns info about a study
    '''
    header = None
    for line in meta:
        if header is None:
            header = line.strip('\n').split('\t')
        else:
            fields = line.strip('\n').split('\t')
            # convert to dict
            candidate = {}
            for idx, field in enumerate(fields):
                candidate[header[idx].strip().lower()] = field.strip()
            if 'sample_id' in candidate and candidate['sample_id'] == study:
                return candidate
    # problem
    return None

def calculate_mean_stats(overall_stats, overall_mean, log):
    '''
        calculate coverage stats of the form: [in mean range, cov >1, >10, >20, >50]
    '''
    write_log(log, 'calculating coverage stats...')
    mean_stats = [0, 0, 0, 0, 0]

    for coverage in overall_stats:
        if coverage > overall_mean * (1.0 - MEAN_RANGE) and coverage < overall_mean * (1.0 + MEAN_RANGE):
            mean_stats[0] += 1
        if coverage >= 1:
            mean_stats[1] += 1
        if coverage >= 10:
            mean_stats[2] += 1
        if coverage >= 20:
            mean_stats[3] += 1
        if coverage >= 50:
            mean_stats[4] += 1

    mean_stats = [100. * x / len(overall_stats) for x in mean_stats]
    write_log(log, 'calculating: done')

    return mean_stats
 
def calculate_summary(report_cov, threshold, log):
    '''
      calculate a summary of coverage across genes in report_cov
    '''
    write_log(log, 'calculating gene summaries...')
    stats = collections.defaultdict(list)
    total_ok = collections.defaultdict(int)
    total = collections.defaultdict(int)
    overall_stats = []
    for idx, line in enumerate(report_cov):
        fields = line.strip().split('\t') # chr, start, end, gene, offset, cov
        if len(fields) > 5:
            gene = fields[3]
            cov = int(fields[5])
            stats[gene].append(cov)
            overall_stats.append(cov)
            total[gene] += 1
            if cov > threshold:
                total_ok[gene] += 1
        if idx % 100000 == 0:
            #write_log(log, 'processed {0} lines... {1} > {2}'.format(idx, cov, threshold))
            write_log(log, 'processed {0} lines...'.format(idx))

    overall_mean = mean(overall_stats)
    mean_stats = calculate_mean_stats(overall_stats, overall_mean, log)

    # now record medians of each gene
    gene_results = {}
    for gene in stats:
        gene_results[gene] = {'ok': 100. * total_ok[gene] / total[gene], 'median': int(median(stats[gene]))}

    return {'mean': overall_mean, 'median': median(overall_stats), 'genes': gene_results, 'mean_stats': mean_stats}

def is_ok(percent, conversion):
    '''return the status for a gene'''
    for cand in conversion.split(','):
        fields = cand.split(':')
        if percent >= int(fields[1]):
            return '<span style="color:{0};">{1}</span>'.format(fields[2], fields[0])

def category(gene, categories):
    '''return the category for a gene'''
    if gene.lower() in categories:
        return int(categories[gene.lower()])
    else:
        return 0

def build_metrics(picard, ontarget, log):
    '''
        parse metric details from picard file
    '''
    stage = 0
    result = {}
    for line in picard:
        if stage == 0:
            if line.startswith('## METRICS CLASS'):
                write_log(log, "found metrics line")
                stage += 1
        elif stage == 1:
            headers = line.strip('\n').split('\t')
            write_log(log, "found headers line: {0}".format(headers))
            stage += 1
        else:
            values = line.strip('\n').split('\t')
            write_log(log, "found values line: {0}".format(values))
            for idx, value in enumerate(values):
                result[headers[idx].lower()] = value
                # ontarget
            for line in ontarget:
                result['on_target'] = int(line.strip())
            write_log(log, result)
            return result
    write_log(log, "ERROR: failed to parse metrics file")

def parse_tsv(tsv):
    '''
        turn lines of the form key\tvalue into a dictionary
    '''
    result = {}
    for line in tsv:
        fields = line.strip('\n').split()
        if len(fields) > 1:
            result[fields[0]] = fields[1]
    return result

def parse_date(date):
    '''
      convert yyyymmdd into something nicer
    '''
    dates = date.strip('"')
    result = []
    for date in dates.split(' '):
        if len(date) >= 8:
            result.append('{0}-{1}-{2}'.format(date[:4], date[4:6], date[6:]))
    if len(result) > 0:
        return ' '.join(result)
    else:
        return 'N/A'

def parse_genes(genes):
    '''
      convert prioritized genes into a list
    '''
    genes = genes.strip('"')
    genes = re.sub('G?[0-9]:', '', genes)
    return re.sub(', *', ' ', genes)

def find_first_of(meta, keys, default='N/A'):
    '''
        return the value of the first key from keys found in meta, otherwise return default
    '''
    for key in keys:
        if key in meta and meta[key] != '':
            return meta[key]
    return default

def generate_report(summary, karyotype, meta, threshold, categories, conversion, metrics, capture, anonymous, fragments, padding, out):
    '''
        generate a report from the provided summary
    '''
    if anonymous:
        meta['sample_id'] = 'ANONYMOUS'

    # headline
    out.write('# Sequencing Summary Report for Study {0}\n\n'.format(meta['sample_id']))

    # summary data
    out.write('\n## Summary Data\n')
    out.write('|\n')
    out.write('-------------|---\n')
    out.write('**Batch**    | {0}\n'.format(find_first_of(meta, ['batch'])))
    out.write('**Study ID** | {0}\n'.format(find_first_of(meta, ['sample_id'])))
    out.write('**Sex**      | {0}\n'.format(find_first_of(meta, ['sex']).upper()))
    if meta['sex'].lower() == karyotype['sex'].lower():
        out.write('**Inferred Sex** | {0}\n'.format(karyotype['sex'].upper()))
    else:
        out.write('**Inferred Sex** | **{0}**\n'.format(karyotype['sex'].upper()))
    out.write('**Disease Cohort** | {0}\n'.format(find_first_of(meta, ['cohort'])))
    out.write('**Hospital/Institution** | {0}\n'.format(find_first_of(meta, ['sequencing_lab', 'hospital_centre'])))
    out.write('**Ethnicity** | {0}\n'.format(find_first_of(meta, ['ethnicity'])))
    out.write('**Prioritized Genes** | {0}\n'.format(parse_genes(meta['prioritised_genes'])))
    out.write('**Consanguinity** | {0}\n'.format(find_first_of(meta, ['consanguinity'])))
    out.write('**Sample Type** | {0}\n'.format(find_first_of(meta, ['sample_type'])))
    out.write('**Sequencing Dates** | {0}\n'.format(parse_date(find_first_of(meta, ['sequencing_date']))))
    out.write('**DNA Collection Dates** | {0}\n'.format(parse_date(find_first_of(meta, ['capture_date']))))
    out.write('**Sequencing Machines** | {0}\n'.format(find_first_of(meta, ['machine_id'])))

    # coverage summary
    out.write('\n## Coverage Summary\n')
    out.write('|\n')
    out.write('-----------------------------|---\n')
    out.write('**Mean Coverage Reported by Lab**   | {0}\n'.format(find_first_of(meta, ['mean_coverage'])))
    out.write('**Observed Mean Coverage**   | {0:.1f}\n'.format(summary['mean']))
    out.write('**Observed Median Coverage** | {0}\n'.format(summary['median']))

    mapped_reads = int(metrics['read_pairs_examined']) * 2
    on_target_reads = metrics['on_target']
    unmapped_reads = int(metrics['unmapped_reads'])
    total_reads = unmapped_reads + mapped_reads + int(metrics['unpaired_reads_examined'])

    out.write('**Total Reads** | {0:,d}\n'.format(total_reads))
    out.write('**Unmapped Reads (% of total)** | {0:,d} ({1:.1f}%)\n'.format(unmapped_reads, 100. * unmapped_reads / total_reads))
    out.write('**Mapped Paired Reads (% of total)** | {0:,d} ({1:.1f}%)\n'.format(mapped_reads, 100. * mapped_reads / total_reads))
    out.write('**% Mapped On Target (% off target)** | {0:.1f}% ({1:.1f}%)\n'.format(100. * on_target_reads / mapped_reads, 100. * (1. - 1. * on_target_reads / mapped_reads)))

    # coverage uniformity
    out.write('**% Coverage within 80% of Mean** | {0:.1f}%\n'.format(summary['mean_stats'][0]))
    out.write('**% Coverage at 1x, 10x, 20x, 50x** | {0:.1f}%, {1:.1f}%, {2:.1f}%, {3:.1f}%\n'.format(summary['mean_stats'][1], summary['mean_stats'][2], summary['mean_stats'][3], summary['mean_stats'][4]))

    # fragments
    if fragments is not None:
        out.write('**Mean fragment size (Std. dev)** | {0:.1f} ({1:.1f})\n'.format(float(fragments['fragment_mean']), float(fragments['fragment_sd'])))
        out.write('**Mean read length (Std. dev)** | {0:.1f} ({1:.1f})\n'.format(float(fragments['read_mean']), float(fragments['read_sd'])))
        out.write('**% Bases >= Q30** | {0:.1f}%\n'.format(100. * float(fragments['base_pass']) / float(fragments['base_count'])))

    # padding
    if padding is not None:
        padding_items = padding.split(',')
        out.write('**Exon Padding: Overall, Indels, SNVs** | {0}, {1}, {2}\n'.format(padding_items[0], padding_items[1], padding_items[2]))

    # gene summary
    out.write('\n## Gene Summary\n')
    out.write('Gene | Category | % > {0}x  | Median | OK? | % in capture\n'.format(threshold))
    out.write('-----|----------|-----------|--------|-----|-------------\n'.format(threshold))
    for gene in sorted(summary['genes'].keys(), key=lambda gene: "{0}{1}".format(9 - category(gene, categories), gene)):
        record = summary['genes'][gene]
        in_capture = capture[gene.lower()] if gene.lower() in capture else 0.0
        out.write('{0} | {1} | {2:.1f} | {3} | {4} | {5:.1f}\n'.format(gene, category(gene, categories), record['ok'], record['median'], is_ok(record['ok'], conversion), in_capture))

    # definitions
    out.write('\n## Definitions\n')

    out.write('\n\n**Mean Coverage Reported by Lab**: the mean coverage reported by the sequencing lab')
    out.write('\n\n**Observed Mean Coverage**: the mean coverage across the capture region')
    out.write('\n\n**Observed Median Coverage**: the median coverage across the capture region')
    out.write('\n\n**Total Reads**: the total number of reads generated by the sequencer')
    out.write('\n\n**Unmapped Reads**: reads that were not mapped to the genome')
    out.write('\n\n**Mapped Paired Reads**: paired reads that were mapped to the genome')
    out.write('\n\n**% Mapped on Target**: the proportion of mapped reads that have any part align to any part of the capture region')
    out.write('\n\n**% Coverage within 20% of Mean**: bases in the capture region with coverage within 20% of the observed mean coverage')
    out.write('\n\n**Mean Fragment Size**: the average distance between correctly mapped and paired reads')
    out.write('\n\n**Mean Read Length**: the average length of all sequenced reads')
    out.write('\n\n**% Bases >= Q30**: what percentage of bases were given a recalibrated quality of at least 30')
    out.write('\n\n**Exon Padding**: when calling variants, how much padding is given to the exon boundary, overall, for indels, and for SNVs')
    out.write('\n\n**Perc**: the percentage of the gene overlapping the capture region with acceptable coverage')
    out.write('\n\n**Median**: the median coverage across the gene overlapping the capture region')
    out.write('\n\n**% in capture**: the proportion of the gene that overlaps the capture region')


def build_categories(categories, prioritized, log):
    '''
        build a dictionary that maps genes to categories
    '''
    result = {}
    for line in categories:
        fields = line.strip().split('\t')
        if len(fields) > 1:
            result[fields[0].lower()] = int(fields[1])

    # override with prioritized
    prioritized = prioritized.strip('"')
    override = 0
    for item in prioritized.split(' '):
        priority, genes = item.split(':')
        for gene in genes.split(','):
            result[gene.strip().lower()] = priority
            override += 1

    write_log(log, 'Added {0} prioritized genes to category list'.format(override))

    return result

def build_capture(coverage, log):
    '''
        data for percent in capture
    '''
    result = {}
    for idx, line in enumerate(coverage):
        field = line.strip('\n').split()
        result[field[0].lower()] = float(field[1])
        if idx % 100000 == 0:
            write_log(log, 'exome_cov: {0} processed {1} -> {2}...'.format(idx, field[0], field[1]))
    return result

def main():
    '''
        parse command line and execute
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate coverage report')
    parser.add_argument('--report_cov', required=True, help='intersected coverage file with genes from bedtools')
    parser.add_argument('--gene_cov', required=True, help='coverage of each gene')
    parser.add_argument('--exome_cov', required=True, help='exome coverage file with genes')
    parser.add_argument('--ontarget', required=False, help='target reads count file')
    parser.add_argument('--metrics', required=False, help='metrics output from Picard')
    #parser.add_argument('--bam', required=False, help='alignment of reads')
    parser.add_argument('--study', required=True, help='ID of study')
    parser.add_argument('--meta', required=True, help='meta data file')
    parser.add_argument('--threshold', type=int, required=True, help='threshold for satisfactory coverage')
    parser.add_argument('--classes', required=True, help='how to categorise results')
    #parser.add_argument('--exome', required=True, #   help='target regions for entire exome')
    parser.add_argument('--gc', required=True, help='gene categories of genes')
    parser.add_argument('--anonymous', action='store_true', required=False, help='do not show study ID')
    parser.add_argument('--write_karyotype', required=False, help='write karyotype details to specified file')
    parser.add_argument('--fragments', required=False, help='file containing fragment statistics')
    parser.add_argument('--padding', required=False, help='comma separated padding stats for all,indel,snv')
    args = parser.parse_args()
    write_log(sys.stderr, 'opening {0} for karyotype'.format(args.exome_cov))
    karyotype = calculate_karyotype(open(args.exome_cov, 'r'), log=sys.stderr)
    summary = calculate_summary(open(args.report_cov, 'r'), args.threshold, log=sys.stderr)
    sample = parse_metadata(open(args.meta, 'r'), args.study)
    if args.write_karyotype:
        write_karyotype(open(args.write_karyotype, 'w'), karyotype, sample)
    categories = build_categories(open(args.gc, 'r'), sample['prioritised_genes'], log=sys.stderr)
    metrics = build_metrics(open(args.metrics, 'r'), open(args.ontarget, 'r'), log=sys.stderr)
    capture = build_capture(open(args.gene_cov, 'r'), log=sys.stderr)
    if args.fragments:
        fragments = parse_tsv(open(args.fragments, 'r'))
    else:
        fragments = None
    generate_report(summary, karyotype, sample, args.threshold, categories, args.classes, metrics, capture, args.anonymous, fragments, args.padding, out=sys.stdout)

if __name__ == '__main__':
    main()
