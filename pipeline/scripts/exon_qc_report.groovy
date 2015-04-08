// vim: ts=4:sw=4:expandtab:cindent
java.util.logging.LogManager.getLogManager().reset();
Cli cli = new Cli()
cli.with {
    refgene 'UCSC RefGene annotation file (download from UCSC website)', args:1, required:true
    targets 'Target regions: coverage will be annotated for regions in the given targets', args:1, required:true
    cov 'BEDTools coverage output', args:1, required:true
    o 'TSV file to write output to', args:1
    x 'XLSX file to write output to (optional)', args:1
}

opts = cli.parse(args)
if(!opts)
    System.exit(1)

if(!opts.o && !opts.x) {
    System.err.println "Please supply at least one of -o or -x to set output file\n"
    cli.usage()
    System.exit(1)
}
    
// Open the coverage file
println "Reading coverage ..."
BEDToolsCoverage covs = new BEDToolsCoverage(opts.cov)
        
// Read the UCSC file
println "Reading RefSeq genes ..."
RefGenes refgene = new RefGenes(opts.refgene)

// Read in the target regions to annotate for
println "Reading Target regions ..."
BED targetRegions = new BED(opts.targets, withExtra:true).load()

// Finding unique set of genes
List<String> allGenes = targetRegions*.extra.unique()

COVERAGE_LEVELS=[1,15,30,50,100]

println "Annotating coverage for $targetRegions.numberOfRanges regions (${allGenes.size()} genes) ..."

ProgressCounter c = new ProgressCounter(withRate:true)


class ExonCov {

    Region exon // the exon region referred to

    List<Map> indices // tx, index pairs for the transcripts that the exon is in

    List<Float> fracs // List of coverage 
}

coverageResults = allGenes.collectEntries { String gene ->
    
    c.count()
    println "Processing $gene"
    
    Regions exons = refgene.getExons(gene)
    
    // Find the overlapping transcripts
    def txes = refgene.getTranscripts(gene)
    
    def coverage = exons.collect { exon ->
        // The exon may have different numbers in different transcripts
        // For each transcript, figure out what index the exon is in that transcript
        def indices = txes.collect {  tx ->
            Integer i = refgene.getTranscriptExons(tx.tx).findIndexOf { 
                it.overlaps(exon) 
            }

            if(i >= 0) {
                [tx: tx.tx, index:i+1]
            }
        }.grep { it != null }

        if(indices.empty) {
            System.err.println "WARNING: No exons in $exons overlap any exons from any transcript in : " + txes.collect { it.toString() }.join(",")
            indices = [ [ tx:txes.collect { it.tx }.join(','), index:-1] ]
        }
      
        // This gets us the per-base coverage over the exon intersected
        // with the target region (because the target region is only part that coverage
        // scores are computed for
        def cov = covs.coverage(exon)

        // println "Coverage for exon $exon = $cov"
        def fracs = COVERAGE_LEVELS.collect { threshold ->
            cov.count { it > threshold }
        }*.div((float)(Math.max(cov.size(),1)))

        new ExonCov(exon: exon, indices:indices, fracs:fracs)
    }

    [ gene, coverage ]
}

format_indices = { ExonCov e ->
    // Only one exon? just exon number
    if(e.indices*.index.unique().size()<=1) {
        e.indices[0].index
    }
    else { // join them all together
        e.indices.collect { it.tx + ':' + it.index }.join(',')
    }
}

if(opts.o) {
    new File(opts.o).withWriter { w ->
        allGenes.each { gene ->
            def exons = coverageResults[gene]
            for(e in exons) {
                w.println(([ gene, e.exon.chr, e.exon.from, e.exon.to+1, format_indices(e) ] + e.fracs).join('\t'))
            }
        }
    }
}

if(opts.x) {
    new ExcelBuilder().build {
        // Summary for all samples in the batch
        def s = sheet("PerExonQC ") {
            row { bold { cells("Gene","Region","Exon"); COVERAGE_LEVELS.each { cell(">"+it + "x") } } }
            allGenes.each { gene ->
                def exons = coverageResults[gene]
                boolean first = true
                for(e in exons) {
                    row { 
                      if(first) {
                          bold { cell(gene) }
                      }
                      else {
                        cell("")
                      }
                      cell("$e.exon.chr:$e.exon.from-${e.exon.to+1}").link("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$e.exon.chr%3A$e.exon.from-$e.exon.to&refGene=pack")
                      applyStyle("left") { cell(format_indices(e)) }
                      e.fracs.each { cell(String.format('%2.2f',it*100)) }
                    }
                    first=false
                }
            }
        }
        s.autoSize()
        s.setColumnWidth(2,20*256)
    }.save(opts.x)
}


