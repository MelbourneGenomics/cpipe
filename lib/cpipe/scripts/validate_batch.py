#!/usr/bin/env python2.7
"""
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.    If not, see <http:#www.gnu.org/licenses/>.
#
#############################################################################
# Purpose:
# * Validate batch results and generate a markdown flavoured report
#
##############################################################################
"""

import argparse
import collections
import glob
import operator
import os
#import subprocess

#def pdf_to_text(file):
#    return subprocess.Popen(["pdftotext",file,"-"], stdout=subprocess.PIPE).communicate()[0]

def md_to_text(md_file):
    """
        given md file, return all lines
    """
    return open(md_file, 'r').readlines()

def extract_sample(sample_file):
    """
        determine sample ID from filename
    """
    sample_file = os.path.basename(sample_file) # remove leading directories
    if '_' in sample_file:
        sample_file = sample_file.rsplit("_", 1)[-1] # remove pipeline run id
    return sample_file.split(".")[0] # remove extensions

def check_sex(dir_name):
    """
        compare inferred sex from karyotype file to sample sex
    """
    print("\n# Gender Validation")
    print("Sample     | Outcome  | Sex    | Inferred")
    print("-----------|----------|--------|---------")
    for karyotype_file in glob.glob(os.path.join(dir_name, '*.karyotype.tsv')):
        result = {}
        for line in open(karyotype_file, 'r'):
            key, value = line.strip().split('\t')
            result[key] = value

        sample = extract_sample(karyotype_file)
        if "Sex" not in result:
            print(("%s | **Sex not found** | | " % sample))
        if "Inferred Sex" not in result:
            print(("%s | **Inferred Sex not found** | | " % sample))
        if "Sex" in result and "Inferred Sex" in result:
            outcome = "OK" if result["Sex"].upper() == result["Inferred Sex"].upper() else "**FAIL*"
            print(("%s | %s | %s | %s" % (sample.ljust(10), outcome.ljust(8), result["Sex"].ljust(6), result["Inferred Sex"].ljust(8))))

def check_gene_coverage(dir_name, bad_threshold=15):
    """
        calculate gene coverage from qc reports
    """
    print(("\n# Gene coverage by sample (flagged if >%i%% fail)" % bad_threshold))
    print("Sample     | Outcome  | % Fail | Good | Pass | Fail | Total")
    print("-----------|----------|--------|------|------|------|------")
    for summary_file in glob.glob(os.path.join(dir_name, '*.summary.htm')):
        out = md_to_text(summary_file)
        result = collections.defaultdict(int)
        for line in out:
            if 'GOOD' in line:
                result['GOOD'] += 1
            if 'FAIL' in line:
                result['FAIL'] += 1
            if 'PASS' in line:
                result['PASS'] += 1
        total = result['GOOD'] + result['FAIL'] + result['PASS']
        bad_percent = 100. * result['FAIL'] / total if total > 0 else 100
        sample = extract_sample(summary_file)
        outcome = "OK" if bad_percent < bad_threshold else "**FAIL**"
        print(("%s | %s | %s | %4i | %4i | %4i | %5i" % (sample.ljust(10), outcome.ljust(8), str('%.1f' % bad_percent).rjust(6), result['GOOD'], result['PASS'], result['FAIL'], total)))

def check_observed_mean_coverage(dir_name, bad_threshold=90):
    """
        extract observed mean coverage from each sample qc report
    """
    print(("\n# Observed mean coverage by sample (flagged if coverage <%i)" % bad_threshold))
    print("Sample     | Outcome  | OMC")
    print("-----------|----------|------")
    for summary_file in glob.glob(os.path.join(dir_name, '*.summary.md')):
        out = md_to_text(summary_file)
        sample = extract_sample(summary_file)
        for line in out:
            if 'Observed Mean Coverage' in line:
                try:
                    omc = float(line.split('|')[1])
                    outcome = "OK\t" if omc > bad_threshold else "**FAIL**"
                    print(("%s | %s | %s" % (sample.ljust(10), outcome.ljust(8), str('%.1f' % omc).rjust(4))))
                    break
                except ValueError:
                    print(("%s | %s | Unexpected string: %s" % (sample, "**FAIL**", line)))

def check_individual_genes(dir_name, bad_threshold=75):
    """
        find genes that fail across multiple samples
    """
    print(("\n# Individual genes with >%i%% fail across samples" % bad_threshold))
    print("Gene     | Outcome  | % Fail | Good | Pass | Fail | Total")
    print("---------|----------|--------|------|------|------|------")
    genes = {}

    for summary_file in glob.glob(os.path.join(dir_name, '*.summary.md')):
        out = md_to_text(summary_file)
        for line in out:
            if 'GOOD' in line or 'FAIL' in line or 'PASS' in line:
                gene = line.split('|')[0].strip()
                if gene not in genes:
                    genes[gene] = {'GOOD': 0, 'FAIL': 0, 'PASS': 0}
                    if 'GOOD' in line:
                        genes[gene]['GOOD'] += 1
                    if 'FAIL' in line:
                        genes[gene]['FAIL'] += 1
                    if 'PASS' in line:
                        genes[gene]['PASS'] += 1

    bad_percent = {}
    for gene in genes:
        bad_percent[gene] = 100. * genes[gene]['FAIL'] / sum([genes[gene][status] for status in genes[gene]])

    for key, value in sorted(list(bad_percent.items()), key=operator.itemgetter(1)):
        if value > bad_threshold:
            outcome = 'OK' if value <= bad_threshold else '**FAIL**'
            print(("%s | %s | %s | %4i | %4i | %4i | %4i" % (key.ljust(8), outcome.ljust(8), str('%.1f' % value).rjust(6), genes[key]['GOOD'], genes[key]['PASS'], genes[key]['FAIL'], sum([genes[key][status] for status in genes[key]]))))

def show_not_found(handle, title):
    """
        list genes not found
    """
    print(("# Requested Genes not found in %s" % title))
    found = False
    for line in handle:
        print(("* %s" % line.strip()))
        found = True
    if not found:
        print("None")

def main():
    """
        validate batch from command line
    """
    parser = argparse.ArgumentParser(description='Validate cpipe output')
    parser.add_argument('--dir', default='./results', help='results directory')
    parser.add_argument('--gene_coverage', default=15, help='report genes with coverage below this')
    parser.add_argument('--mean_coverage', default=90, help='report batches with mean coverage below this')
    parser.add_argument('--gene_sample_fail', default=80, help='report genes that fail in more than this proportion of samples')
    parser.add_argument('--missing_exons', required=False, help='file containing genes not in exons')
    parser.add_argument('--missing_annovar', required=False, help='file containing genes not in annovar')
    parser.add_argument('--excluded_genes', required=False, help='file containing excluded genes')
    args = parser.parse_args()
    check_sex(args.dir)
    check_gene_coverage(args.dir, bad_threshold=args.gene_coverage)
    check_observed_mean_coverage(args.dir, bad_threshold=args.mean_coverage)
    check_individual_genes(args.dir, bad_threshold=args.gene_sample_fail)
    print("")
    if args.missing_exons and os.path.isfile(args.missing_exons):
        show_not_found(open(args.missing_exons, 'r'), 'Reference')
    else:
        print("* No missing gene information at exon level")
    print("")
    if args.missing_annovar and os.path.isfile(args.missing_annovar):
        show_not_found(open(args.missing_annovar, 'r'), 'Annovar')
    else:
        print("* No missing gene information at annovar level")
    print("")
    if args.excluded_genes and os.path.isfile(args.excluded_genes):
        print("# Excluded genes found in gene lists")
        for line in open(args.excluded_genes, 'r'):
            print((line.strip()))
if __name__ == '__main__':
    main()
