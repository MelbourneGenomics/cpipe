# vim: ts=4:expandtab:sw=4:cindent
########################################################################
#
# Reads a bed file containing genes and expands it to exon definitions
# by relating the genes to exons defined in the hg19 refseq gene file
# 
# The output is written to the third / last command line argument
#
########################################################################
import sys, csv, getopt, re

# Whether to include UTR regions 
include_utr = True

opts,args = getopt.getopt(sys.argv[1:],"c",)
for opt,value in opts:
        if opt == '-c':
               include_utr = False 

if not args or len(args)<4:
        print >>sys.stderr, "\nUsage: python %s [-c] <gene bed file> <hg19 UCSC RefSeq genes file> <transcript file> <output file>\n" % sys.argv[0]
        print >>sys.stderr, "\t-c   do not include UTR\n"
        sys.exit(1)


# Read all of the transcripts
tx_file = args[2]
vcgs_txes = [ re.sub('\.[0-9]*$', '',line).strip() for line in open(tx_file) ]

print >>sys.stderr, "The vcgs txes are %s" % vcgs_txes


def check_overlap(existing_exons,newexon):            
    """
    Return an overlapping existing exon if the new exon overlaps an existing exon
    """
    start,end = newexon[0],newexon[1]
    i=0
    for existing in existing_exons:
        # exact match
        if newexon == existing:
           return i
        # overlaps start point
        elif start<existing[0] and end > existing[0]:
            return i
        # contained within existing
        elif start>=existing[0] and end <= existing[1]:
            return i
        # overlaps end point
        elif start<=existing[1] and end > existing[1]:
            return i
        elif start<=existing[0] and end >= existing[1]:
            return i
        i+=1

    return None


# Start by reading the genes we are interested in (don't care about
# coordinates)
gene_file = csv.reader(open(args[0]), delimiter='\t')
genes = {}
gene_ranges = {}
gene_chr = {}
for g in gene_file:
    genes[g[3]] = []
    gene_chr[g[3]] = g[0]
    gene_ranges[g[3]] = [int(g[1])-1,int(g[2])+1]

print >>sys.stderr,  "Read %d genes from gene file" % len(genes)

cds_starts = {}
cds_ends = {}

vcgs_prioritized_txes = { }

# Now create a hash in memory of exons indexed by gene
# We are reading the whole hg19 RefSeq gene annotation database into memory here
refseq_genes = csv.reader(open(args[1]), delimiter='\t')
for g in refseq_genes:
    gene = g[12]
    # debug: limit to specific gene
    # if gene != 'DSP':
    #  continue

    if gene in genes:

        # Exon start positions in column 9
        starts = map(lambda x: int(x), filter(lambda x: x != '', g[9].split(",")))

        # Exon end positions in column 10
        ends = map(lambda x: int(x), filter(lambda x: x != '', g[10].split(",")))

        exons = map(lambda x: list(x), zip(starts,ends))
        print >>sys.stderr,  "Found gene %s with transcript %s (%d exons)" % (g[12],g[1],len(exons))

        if g[2] != gene_chr[gene]:
            raise Exception("Gene %s is annotated to multiple chromosomes: %s vs %s" % (gene, gene_chr[gene], g[2]))

        cds_start = int(g[6])
        cds_end = int(g[7])

        # Update the index of coding sequence starts / ends
        # This is to be used later when we adjust exons based on 
        # UTR inclusion
        cds_starts[gene] = min(cds_starts.get(gene,sys.maxint),cds_start)
        cds_ends[gene] = max(cds_ends.get(gene,0),cds_end)

        existing_exons = genes[gene]

        #print "Existing exons = %s" % existing_exons
        gene_range = gene_ranges[gene]

        # Merge the exons with existing ones
        for e in exons:
            if e[0] > gene_range[1] or e[1] < gene_range[0]:
                print >>sys.stderr,  "Exon %s in gene %s outside range %s" % (e,gene,gene_range)
                continue

            overlapping = check_overlap(existing_exons,e)
            if overlapping is not None:

                old_start = existing_exons[overlapping][0]
                new_start = min(e[0],old_start)

                old_end = existing_exons[overlapping][1]
                new_end = max(e[1],old_end)

                if new_start == old_start and new_end == old_end:
                    continue

                if g[1] in vcgs_txes:
                    print >>sys.stderr,  "Exon %s in gene %s has multiple potential splice regions! Selecting VCGS prioritized tx=%s : will use %d-%d as coding sequence" % (e,gene,g[1],e[0],e[1])
                    existing_exons[overlapping][0] = e[0]
                    existing_exons[overlapping][1] = e[1]
                    vcgs_prioritized_txes[gene] = g[1]
                elif gene in vcgs_prioritized_txes:
                    print >>sys.stderr,\
                        "Exon %s in gene %s has multiple potential splice regions! Not using longest sequence because a VCGS prioritized tx (%s) exists vs tx=%s : will use %d-%d as coding sequence" \
                        % (e,gene,vcgs_prioritized_txes[gene], g[1],new_start,new_end)
                else:
                    print >>sys.stderr,  "Exon %s in gene %s has multiple potential splice regions! Current tx=%s Will use %d-%d as longest coding sequence" % (e,gene,g[1],new_start,new_end)
                    existing_exons[overlapping][0] = new_start
                    existing_exons[overlapping][1] = new_end
            else:
                existing_exons.append(e)

# Exons include the UTR by default, so if it should not be included,
# trim the first and last exons 
if not include_utr:
    for g in genes:
        if g in cds_starts:
            exons = genes[g]

            # NOTE: in some cases the entire first / last exon are
            # non-coding. In that case they seem to be NOT considered
            # to be UTR, so they are not trimmed

            # Move start of first exon to start of CDS for gene
            if len(exons) == 0:
                print "WARNING: no exons overlapped by coordinates for gene %s: " % g
            else:
                first_exon = exons[min(range(len(exons)), key=lambda i: exons[i][0])]

                if first_exon[1] > cds_starts[g]:
                    first_exon[0] = cds_starts[g]

                # Move the end of the last exon to the end of the CDS
                # Note that exons are not sorted, so we have to find 
                # the last exon ...
                last_exon = exons[max(range(len(exons)), key=lambda i: exons[i][1])]
                if last_exon[0] < cds_ends[g]:
                    last_exon[1] = cds_ends[g]
            

# Write out result bed file
output = csv.writer(open(args[3],'wb'), delimiter='\t', lineterminator='\n')
for g in genes:
    exon_count = 0
    for e in genes[g]:
        exon_count += 1
        # print ','.join(map(lambda x: str(x), [ gene_chr[g], e[0], e[1], "%s|%d" % (g, exon_count)])) 
        output.writerow( [ gene_chr[g], e[0], e[1], "%s|%d" % (g, exon_count)] )

