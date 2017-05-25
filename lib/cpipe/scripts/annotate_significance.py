#!/usr/bin/env python
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################

####################################################################################
#
# Purpose:
#
# This script adds annotations to an Annovar output to rate the
# clinical significance of each variant as a "priority index". These rules are
# defined by the Bioinformatics working group and are currently
# specified as follows:
#
#     (1)     Index 1 - missense
#     (1)     Index 2 - rare missense
#             Variant has MAF <0.01 in dbSNP or ESP or 1000G databases
#     (2)     Index 3 - novel missense
#             Variant is not in dbSNP or ESP or 1000G databases
#     (3)     Index 4 - highly conserved missense
#             Variant has condel score > 0.7 (other possible filters - for discussion)
#     (4)     Index 5 - truncating
#             Out of frame indel (non-recurrent*), nonsense, splice site +/-2bp
#
#     NOTE: these rules were updated 27/5/2014:
#
#             - include 'very rare' variants <0.0005 MAF in cat 2,3
#             - include truncating variants into priority 1 if rare (but not novel)
#
#     NOTE: these rules were updated 2015-6-10:
#
#             - all priorities shifted up by 1
#             - new priority 1 introduced to capture all misense variants
#                (even non-rare)
#
# Author:   Simon Sadedin, simon.sadedin@mcri.edu.au
# Date:     23/1/2014
#
####################################################################################
"""

import argparse
import collections
import csv
import logging as log
import sys

log.basicConfig(level=log.INFO)

class Annovar(object):
    """
    Helper class to map Annovar column names to their fields parsed from CSV,
    and to implement logic surrounding categorization.

    The function init_columns() must be called, passing the header row as
    returned by csv.reader() to initialise the class before use.
    """

    # Column names of Annovar file
    columns = []

    # Default MAF threshold for considering a variant 'rare'
    MAF_THRESHOLD = 0.01

    # Default MAF threshold for considering a variant 'very rare'
    MAF_THRESHOLD_VERY_RARE = 0.0005

    # Default Condel Threshold
    CONDEL_THRESHOLD = 0.7

    # Categories of variants as specified by Annovar, mapped to functional categories
    # defined for Melbourne Genomics
    ANNOVAR_EXONIC_FUNCS = {
        "truncating" : ["frameshift insertion", "frameshift deletion", "frameshift substitution", "stopgain SNV", "stoploss SNV", "stoploss", "stopgain"],
        "missense" : ["nonframeshift insertion", "nonframeshift deletion", "nonframeshift substitution", "nonsynonymous SNV"],
        "synonymous" : ["synonymous SNV"],
        "noncoding" : ["intronic", "intergenic", "ncRNA_intronic", "ncRNA_exonic", "upstream", "downstream", "UTR5", "UTR3", "ncRNA_splicing", "ncRNA_exonic;splicing", "upstream;downstream", "UTR5;UTR3"]

    }

    # These are the Annovar fields that contain population frequency estimates
    # Note we do some fooling around in the maf_value() method to maintain
    # compatibility with different versions of Annovar
    POPULATION_FREQ_FIELDS = ["esp6500siv2_all", "1000g2014oct_all", "exac03"]

    def __init__(self, line, synonymous=None):
        self.line = line
        self.synonymous = synonymous

    def priority(self):
        """
            Main logic describing how to map any given variant to a clinical significance
            priority index. See the main header for the definition of these categories.

            Note: unknown categories are returned as 9 - that is, extremely high.
        """

        if self.is_missense(): # nonframeshift...
            if self.is_rare():
                if self.is_novel() or self.is_very_rare():
                    if self.is_conserved():
                        return 4 # Missense, novel and conserved => category 4
                    else:
                        return 3 # Missesnse, novel but not highly conserved => category 3
                else:
                    return 2 # Missense & rare but not novel => category 2
            else:
                return 1 # Missense but not even rare => category 1

        elif self.is_truncating(): # frameshift, stopgain, stoploss
            # From Natalie, 27/5/2014:
            # With regard to priority 5 truncating variants:
            #  novel should stay in priority 5
            #  rare should be priority 2
            if self.is_novel():
                return 5
            elif self.is_rare():
                log.debug("%s:%s is rare", self.Chr, self.Start)
                return 2
            else:
                return 1

        elif self.is_noncoding():
            return 0

        elif self.ExonicFunc == "synonymous SNV":
            # From Natalie, 18/11/15
            if '{0},{1}'.format(self.Chr, self.Start) in self.synonymous:
                if self.is_novel():
                    log.info("variant %s:%s %s/%s func=%s not filtered due to exon boundary proximity", self.Chr, self.Start, self.Ref, self.Obs, self.ExonicFunc)
                    return 5
                elif self.is_rare():
                    log.info("variant %s:%s %s/%s func=%s not filtered due to exon boundary proximity", self.Chr, self.Start, self.Ref, self.Obs, self.ExonicFunc)
                    return 2
                else:
                    return 0
            else:
                return 0
        elif self.ExonicFunc == "unknown":
            return 0
        else:
            log.warning("variant %s:%s %s/%s func=%s failed to be categorized", self.Chr, self.Start, self.Ref, self.Alt, self.ExonicFunc)
            return 9

    def is_noncoding(self):
        """is the variant non coding"""
        return self.Func in self.ANNOVAR_EXONIC_FUNCS["noncoding"]

    def is_missense(self):
        """is the variant missense"""
        return self.ExonicFunc in self.ANNOVAR_EXONIC_FUNCS["missense"]

    def is_truncating(self):
        """is the variant truncating or splicing"""
        return self.ExonicFunc in self.ANNOVAR_EXONIC_FUNCS["truncating"] or self.Func in ["splicing", "exonic;splicing"]

    def is_rare(self):
        """Return true iff at least one database has the variant at > the MAF_THRESHOLD"""
        log.debug("MAF values for %s:%s are %s", self.Chr, self.Start, [self.maf_value(f) for f in self.POPULATION_FREQ_FIELDS])
        return not any([self.maf_value(f) > self.MAF_THRESHOLD for f in self.POPULATION_FREQ_FIELDS])

    def is_very_rare(self):
        """Return true iff at least one database has the variant at > the MAF_THRESHOLD_VERY_RARE"""
        return not any([self.maf_value(f) > self.MAF_THRESHOLD_VERY_RARE for f in self.POPULATION_FREQ_FIELDS])

    def is_novel(self):
        """return true iff the variant has no MAF in any database AND no DBSNP ID"""
        return not any([self.maf_value(f) > 0.0 for f in self.POPULATION_FREQ_FIELDS]) and (self.snp138 in ["", "."])

    def is_conserved(self):
        """
            Clarification 27/5/2014:
            ONLY if condel score is missing, then it can categorised as a 3 if CONSERVED by Annovar
        """
        condel_str = self.Condel
        if condel_str != "":
            return float(condel_str) >= 0.7
        else:
            return self.phastConsElements46way != ""

    @staticmethod
    def init_columns(cols):
        """
            prepare column names in annovar
        """
        Annovar.columns = cols #+ ["MapQ", "QD"]

    def maf_value(self, name):
        """
            Trying to be compatible with multiple versions of Annovar, each having different
            names for this column
        """
        if name == "exac03" and "ExAC_Freq" in self.columns:
            name = "ExAC_Freq"
        if name == "exac03" and "ExAC_ALL" in self.columns:
            name = "ExAC_ALL"
        value = self.line[self.columns.index(name)]
        if value == "" or value == ".":
            return 0
        else:
            return float(value)

    def __getattr__(self, name):
        return self.line[self.columns.index(name)]

    def set_value(self, name, value):
        """
            update value of a given field at the current line
        """
        self.line[self.columns.index(name)] = value

def process_annovar(annovar, output, synonymous=None):
    """
        annotate priority
    """
    log.info("started processing...")
    # prepare synonymous set
    synonymous_set = set()
    if synonymous is not None:
        for line in synonymous:
            fields = line.strip().split('\t')
            if len(fields) > 2:
                for field in range(int(fields[1]), int(fields[2])):
                    key = '{0},{1}'.format(fields[0], field)
                    synonymous_set.add(key)
    log.info("finished reading synonymous set: {0} positions.".format(len(synonymous_set)))

    # Read the file
    reader = csv.reader(annovar, delimiter=',', quotechar='"', doublequote=True)

    # Open CSV writer to standard output, first for header (for body comes in the loop below)
    header_out = csv.writer(output, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONE)
    is_header = True
    priorities = collections.defaultdict(int)
    log.info("calculating priorities with thresholds: rare {0}, very rare {1}, condel {2}".format(Annovar.MAF_THRESHOLD, Annovar.MAF_THRESHOLD_VERY_RARE, Annovar.CONDEL_THRESHOLD))
    for line in reader:

        if is_header:
            is_header = False
            # Note: Annovar does not seem to provide Qual and Depth headings itself
            if "Qual" not in line:
                line += ["Qual"]

            if "Depth" not in line:
                line += ["Depth"]

            Annovar.init_columns(line)

            header_out.writerow(Annovar.columns + ["Priority_Index"])
            output.flush()
            csv_output = csv.writer(output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            continue

        annovar = Annovar(line, synonymous_set)

        while len(line) < len(Annovar.columns):
            line.append("")

        priority = annovar.priority()
        priorities[priority] += 1
        csv_output.writerow(line + [priority])
    log.info("priority distribution: {0}".format(priorities))

####################################################################################
#
# Main body
#
####################################################################################

def main():
    """
        Parse command line options
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--annovar', required=True, help='annovar file')
    parser.add_argument('--rare', required=False, help='threshold for rare')
    parser.add_argument('--very_rare', required=False, help='threshold for very rare')
    parser.add_argument('--condel', required=False, help='threshold for condel')
    parser.add_argument('--synonymous', required=False, help='bed file allowing synonymous variants')
    args = parser.parse_args()

    if args.rare:
        Annovar.MAF_THRESHOLD = float(args.rare)

    if args.condel:
        Annovar.CONDEL_THRESHOLD = float(args.condel)

    if args.very_rare:
        Annovar.MAF_THRESHOLD_VERY_RARE = float(args.very_rare)

    if args.synonymous:
        process_annovar(open(args.annovar, 'r'), sys.stdout, synonymous=open(args.synonymous, 'r'))
    else:
        process_annovar(open(args.annovar, 'r'), sys.stdout)

if __name__ == "__main__":
    main()
