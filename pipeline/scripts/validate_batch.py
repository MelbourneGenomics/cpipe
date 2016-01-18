#!/usr/bin/env python2.7
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
####################################################################################
# Purpose:
# * Validate batch results and generate a markdown flavoured report
#
# Dependencies:
# * pdftotext
####################################################################################

import argparse
import collections
import glob
import operator
import os
import subprocess

def pdf_to_text(file):
  return subprocess.Popen(["pdftotext",file,"-"], stdout=subprocess.PIPE).communicate()[0]

def extract_sample( file ):
  file = os.path.basename(file) # remove leading directories
  if '_' in file:
    file = file.rsplit("_", 1)[-1] # remove pipeline run id
  return file.split(".")[0] # remove extensions

def check_sex( dir ):
  print "\n# Gender Validation"
  print "Sample     | Outcome  | Sex    | Inferred"
  print "-----------|----------|--------|---------"
  for file in glob.glob( os.path.join( dir, '*.karyotype.tsv' ) ):
    result = {}
    for line in open( file, 'r' ):
      key, value = line.strip().split('\t')
      result[key] = value
      
    sample = extract_sample( file )
    if "Sex" not in result:
      print "%s | **Sex not found** | | " % sample
    if "Inferred Sex" not in result:
      print "%s | **Inferred Sex not found** | | " % sample
    if "Sex" in result and "Inferred Sex" in result:
      outcome = "OK" if result["Sex"] == result["Inferred Sex"] else "**FAIL*"
      print "%s | %s | %s | %s" % ( sample.ljust(10), outcome.ljust(8), result["Sex"].ljust(6), result["Inferred Sex"].ljust(8) )

def check_gene_coverage( dir, bad_threshold=15 ):
  print "\n# Gene coverage by sample (flagged if >%i%% fail)" % bad_threshold
  print "Sample     | Outcome  | % Fail | Good | Pass | Fail | Total"
  print "-----------|----------|--------|------|------|------|------"
  for file in glob.glob( os.path.join( dir, '*.summary.pdf' ) ):
    out = pdf_to_text(file)
    result = collections.defaultdict(int)
    for line in out.split('\n'):
      result[line.strip()] += 1
    total = result['GOOD'] + result['FAIL'] + result['PASS']
    bad_percent = 100. * result['FAIL'] / total if total > 0 else 100
    sample = extract_sample( file )
    outcome = "OK" if bad_percent < bad_threshold else "**FAIL**"
    print "%s | %s | %s | %4i | %4i | %4i | %5i" % ( sample.ljust(10), outcome.ljust(8), str( '%.1f' % bad_percent ).rjust(6), result['GOOD'], result['PASS'], result['FAIL'], total )

def check_observed_mean_coverage( dir, bad_threshold=90 ):
  print "\n# Observed mean coverage by sample (flagged if coverage <%i)" % bad_threshold
  print "Sample     | Outcome  | OMC"
  print "-----------|----------|------"
  for file in glob.glob( os.path.join( dir, '*.summary.pdf' ) ):
    out = pdf_to_text(file)
    result = collections.defaultdict(int)
    in_omc = False
    sample = extract_sample( file )
    for line in out.split('\n'):
      line = line.strip()
      if in_omc:
        if len(line) > 0:
          try:
            omc = float(line)
            outcome = "OK\t" if omc > bad_threshold else "**FAIL**"
            print "%s | %s | %s" % ( sample.ljust(10), outcome.ljust(8), str( '%.1f' % omc ).rjust(4) )
            break
          except:
            print "%s | %s | Unexpected string: %s" % ( sample, "**FAIL**", line )
      else:
        if line == 'Observed Mean Coverage':
          in_omc = True

def check_individual_genes( dir, bad_threshold=75 ):
  print "\n# Individual genes with >%i%% fail across samples" % bad_threshold
  print "Gene     | Outcome  | % Fail | Good | Pass | Fail | Total"
  print "---------|----------|--------|------|------|------|------"
  genes = {}

  for file in glob.glob( os.path.join( dir, '*.summary.pdf' ) ):
    out = pdf_to_text(file)
    result = collections.defaultdict(int)
    last = None
    for line in out.split('\n'):
      line = line.strip()
      if line in ('GOOD', 'FAIL', 'PASS') and last is not None:
        if last not in genes:
          genes[last] = { 'GOOD': 0, 'FAIL': 0, 'PASS': 0 }
        genes[last][line] += 1
      elif len(line) > 0 and '%' not in line:
        try:
          _ = float(line)
        except:
          last = line

  bad_percent = {}
  for gene in genes:
    bad_percent[gene] = 100. * genes[gene]['FAIL'] / sum( [ genes[gene][status] for status in genes[gene] ] )

  for key, value in sorted( bad_percent.items(), key=operator.itemgetter(1)):
    if value > bad_threshold:
      outcome = 'OK' if value <= bad_threshold else '**FAIL**'
      print "%s | %s | %s | %4i | %4i | %4i | %4i" % ( key.ljust(8), outcome.ljust(8), str( '%.1f' % value ).rjust(6), genes[key]['GOOD'], genes[key]['PASS'], genes[key]['FAIL'], sum( [ genes[key][status] for status in genes[key] ] ) )

def show_not_found( fh, title ):
  print "# Requested Genes not found in %s" % title
  found = False
  for line in fh:
    print "* %s" % line.strip()
    found = True
  if not found:
    print "None"

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Validate cpipe output') 
  parser.add_argument('--dir', default='./results', help='results directory')
  parser.add_argument('--gene_coverage', default=15, help='report genes with coverage below this')
  parser.add_argument('--mean_coverage', default=90, help='report batches with mean coverage below this')
  parser.add_argument('--gene_sample_fail', default=80, help='report genes that fail in more than this proportion of samples')
  parser.add_argument('--missing_exons', required=False, help='file containing genes not in exons' )
  parser.add_argument('--missing_annovar', required=False, help='file containing genes not in annovar' )
  parser.add_argument('--excluded_genes', required=False, help='file containing excluded genes' )
  args = parser.parse_args()
  check_sex( args.dir )
  check_gene_coverage( args.dir, bad_threshold=args.gene_coverage )
  check_observed_mean_coverage( args.dir, bad_threshold=args.mean_coverage )
  check_individual_genes( args.dir, bad_threshold=args.gene_sample_fail )
  print ""
  if args.missing_exons and os.path.isfile(args.missing_exons):
    show_not_found( open( args.missing_exons, 'r' ), 'Reference' )
  else:
    print "* No missing gene information at exon level"
  print ""
  if args.missing_annovar and os.path.isfile(args.missing_annovar):
    show_not_found( open( args.missing_annovar, 'r' ), 'Annovar' )
  else:
    print "* No missing gene information at annovar level"
  print ""
  if args.excluded_genes and os.path.isfile(args.excluded_genes):
    print "# Excluded genes found in gene lists"
    for line in open( args.excluded_genes, 'r' ):
      print line.strip()
