#!/usr/bin/env python
# vim: expandtab:ts=4:sw=4:
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
#     (1)     Index 1 - rare missense
#             Variant has MAF <0.01 in dbSNP or ESP or 1000G databases
#     (2)     Index 2 - novel missense
#             Variant is not in dbSNP or ESP or 1000G databases 
#     (3)     Index 3 - highly conserved missense
#             Variant has condel score > 0.7 (other possible filters - for discussion)
#     (4)     Index 4 - truncating
#             Out of frame indel (non-recurrent*), nonsense, splice site +/-2bp
#
#     NOTE: these rules were updated 27/5/2014:
#
#             - include 'very rare' variants <0.0005 MAF in cat 2,3
#             - include truncating variants into priority 1 if rare (but not novel)
#
# Author:   Simon Sadedin, simon.sadedin@mcri.edu.au
# Date:     23/1/2014
#
####################################################################################

import csv, getopt, sys, logging as log

log.basicConfig(level=log.DEBUG)

class Annovar:
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
        "truncating" : ["frameshift insertion","frameshift deletion","frameshift substitution","stopgain SNV","stoploss SNV"],
        "missense" : ["nonframeshift insertion","nonframeshift deletion","nonframeshift substitution","nonsynonymous SNV"],
        "synonymous" : ["synonymous SNV"],
        "noncoding" : ["intronic","intergenic","ncRNA_intronic"]
    }

    # These are the Annovar fields that contain population frequency estimates
    POPULATION_FREQ_FIELDS = ["esp5400_all", "1000g2010nov_all","exac03"]

    def __init__(self, line):
        self.line = line

    def priority(self):
        """
            Main logic describing how to map any given variant to a clinical significance
            priority index. See the main header for the definition of these categories.

            Note: unknown categories are returned as 9 - that is, extremely high.
        """

        if self.is_missense():
           if self.is_rare():
               if self.is_novel() or self.is_very_rare():
                   if self.is_conserved():
                       return 3 # Missense, novel and conserved => category 3
                   else:
                       return 2 # Missesnse, novel but not highly conserved => category 2
               else:
                    return 1 # Missense & rare but not novel => category 1
           else:
               return 0 # Missense but not even rare => no category

        elif self.is_truncating():
            # From Natalie, 27/5/2014:
            # With regard to priority 4 truncating variants:
            #  novel should stay in priority 4
            #  rare should be priority 1
            if self.is_novel():
                return 4
            elif self.is_rare():
                log.debug("%s:%s is rare" % (self.Chr,self.Start))
                return 1
            else:
                return 0

        elif self.is_noncoding():
            return 0

        elif self.ExonicFunc in ["synonymous SNV", "unknown"]:
            return 0
        else:
            print >>sys.stderr, "WARNING: variant %s:%s %s/%s func=%s failed to be categorized" % \
                    (self.Chr, self.Start, self.Ref, self.Alt, self.ExonicFunc)
            return 9
        
    def is_noncoding(self):
        return self.Func in self.ANNOVAR_EXONIC_FUNCS["noncoding"]

    def is_missense(self):
        return self.ExonicFunc in self.ANNOVAR_EXONIC_FUNCS["missense"]

    def is_truncating(self):
        return self.ExonicFunc in self.ANNOVAR_EXONIC_FUNCS["truncating"] or self.Func in ["splicing","exonic;splicing"]

    def is_rare(self):
        # Return true iff at least one database has the variant at > the MAF_THRESHOLD
        log.debug("MAF values for %s:%s are %s", self.Chr, self.Start, map(lambda f: self.maf_value(f),self.POPULATION_FREQ_FIELDS))
        return not any(map(lambda f: self.maf_value(f)>self.MAF_THRESHOLD, self.POPULATION_FREQ_FIELDS))
        return not any(map(lambda f: self.maf_value(f)>self.MAF_THRESHOLD, self.POPULATION_FREQ_FIELDS))

    def is_very_rare(self):
        # Return true iff at least one database has the variant at > the MAF_THRESHOLD_VERY_RARE
        return not any(map(lambda f: self.maf_value(f)>self.MAF_THRESHOLD_VERY_RARE, self.POPULATION_FREQ_FIELDS))

    def is_novel(self):
        # return true iff the variant has no MAF in any database AND no DBSNP ID
        return not any(map(lambda f: self.maf_value(f) > 0.0, self.POPULATION_FREQ_FIELDS)) and (self.snp138 in ["","."])

    def is_conserved(self):
        # Clarification 27/5/2014:
        # ONLY if condel score is missing, then it can categorised as a 3 if CONSERVED by Annovar
        condel_str = self.Condel 
        if condel_str != "":
            return float(condel_str) >= 0.7
        else:
            return self.phastConsElements46way != ""

    @staticmethod
    def init_columns(cols):
        Annovar.columns = cols #+ ["MapQ","QD"]

    def maf_value(self, name):
        value = self.line[self.columns.index(name)]
        if value == "" or value ==".":
            return 0
        else:
            return float(value)

    def __getattr__(self,name):
        return self.line[self.columns.index(name)]
    
    def set_value(self,name,value):
        self.line[self.columns.index(name)]=value
    
####################################################################################
#
# Main body
#
####################################################################################

def main():
    # Parse command line options
    optstring = "a:"
    opts,args = getopt.getopt(sys.argv[1:],optstring)
    
    options = {}
    for opt in opts:
       options[opt[0]] = opt[1]
    
    def usage(msg):
        print >>sys.stderr, "\nERROR: %s\n\nUsage: annotate_significance.py -a <annovar file>\n" % msg
        sys.exit(1)
            
    if not '-a' in options:
        usage("Please provide -a option.")
    
    # Read the file
    reader = csv.reader(open(options["-a"]), delimiter=',', quotechar='"', doublequote=True)
    
    # Open CSV writer to standard output, first for header (for body comes in the loop below)
    header_out = csv.writer(sys.stdout, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONE)
    is_header = True
    for line in reader:
    
        if is_header:
            is_header = False
            # Note: Annovar does not seem to provide Qual and Depth headings itself
            if "Qual" not in line:
                line = line + ["Qual"]
            
            if "Depth" not in line:
                line = line + ["Depth"]
    
            Annovar.init_columns(line)
    
            header_out.writerow(Annovar.columns + ["Priority_Index"])
            sys.stdout.flush()
            output = csv.writer(sys.stdout, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            continue
    
        av = Annovar(line)
    
        while len(line)<len(Annovar.columns):
                line.append("")
          
        output.writerow(line + [av.priority()])
    
if __name__ == "__main__":    
    main()
