
# convert a bed file into a list of genes

import sys

genes = set()
for line in sys.stdin:
  if line.startswith( '#' ):
    continue
  fields = line.strip().split()
  if len(fields) > 3:
    gene = fields[3].strip().upper()
    genes.add( gene )

for gene in sorted( list( genes ) ):
  sys.stdout.write( '{0}\t{1}\n'.format( gene, 1 ) )
