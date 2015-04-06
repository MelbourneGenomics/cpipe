# vim: sw=4:expandtab:cindent:ts=4
###################################################################
#
# Reads the output summary from Annovar and adds a column to indicate 
# whether the mutation in the row affects one of the transcripts that
# are identified in the prioritised transcripts file.
#
###################################################################

import csv,sys,re

# Transcripts of interest
transcripts = sys.argv[1]

# Summary exome_summary.csv file from Annovar
summary = sys.argv[2]

# Full exonic_variant_function file from Annovar
full = sys.argv[3]

# Read all the transcripts of interest
txs = {}
for i in open(sys.argv[1]):
        txs[re.sub('\\.[0-9].*$', '', i.strip())] = True

w = csv.writer(sys.stdout)

# Read the CSV summary file and examine each variant
reader = csv.reader(open(summary))

# First read the header and since by default some columns don't have headers,
# as a side benefit we fix those
header = reader.next()

# Fix missing column headings
if len(header)<28:
        header.append('Qual')
        header.append('Depth')
        header.append('PTY_TX')

w.writerow(header)

debug = False

for l in reader:
        gene = l[1]
        chr = l[21]
        start = l[22]
        aachange = l[3]

        if gene == 'unknown':
                continue

        # Search for the aachange in the full file
        # to get the full list of isoforms / transcripts and see if any are 
        # flagged as of interest
        f = open(full)
        row = l
        found_vtx = ''
        for v in csv.reader(f,delimiter='\t'):
                if v[1] == 'unknown':
                        continue
                
                # Same location
                if v[3] != chr:
                        continue

                if v[4] != start:
                        continue

                # Has to be the same gene
                if gene not in map(lambda x: x.split(':')[0], v[2].split(',')):
                        continue

                # Column 2 is in the following format:
                # TTN:NM_003319:exon73:c.G18133A:p.D6045N,TTN:NM_133432:exon74:c.G18508A:p.D6170N
                # We want to report only the transcript and the AA change
        
                try:
                    vtxs = zip(map(lambda x: x.split(':')[1], v[2].split(',')), # the NM_.. transcript id
                           map(lambda x: ':'.join(x.split(':')[3:5]), v[2].split(','))) # the AA change
                    
                    if debug:
                        print "Full transcripts for %s are %s" % (aachange, vtxs)

                    # vtxs is a list of transcripts, each element is another list
                    # of 2 elements, (tx name, aa change)
                    vtxs_flag = filter(lambda x: x[0] in txs, vtxs)
                    if vtxs_flag:
                            if debug:
                                print "Full transcripts for %s are %s" % (aachange, vtxs)
                                print "FLAGGED: %s" % vtxs_flag
                            found_vtx = ";".join(map(lambda f: ":".join(f), vtxs_flag))
                except Exception,e:
                    found_vtx = 'Error: Please check manually'

        row.append(found_vtx)
        w.writerow(row)
        f.close()
