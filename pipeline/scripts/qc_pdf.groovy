// vim: shiftwidth=4:ts=4:expandtab:
/////////////////////////////////////////////////////////////////////////////////
//
// This file is part of Cpipe.
// 
// Cpipe is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, under version 3 of the License, subject
// to additional terms compatible with the GNU General Public License version 3,
// specified in the LICENSE file that is part of the Cpipe distribution.
//
// Cpipe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Cpipe.  If not, see <http://www.gnu.org/licenses/>.
//
/////////////////////////////////////////////////////////////////////////////////
//
// Sequencing Summary Report
//
// This script reads the coverage and sample information and
// creates a summary report suitable for high level evaluation
// of the sequencing quality and success in PDF form.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////
import java.awt.Color

/////////////////////////////////////////////////////////////////////////
//
// Parse command line options
//
/////////////////////////////////////////////////////////////////////////
CliBuilder cli = new CliBuilder(usage: "qc_pdf.groovy <options>")

cli.with { 
    cov "Coverage file from Bedtools", args:1, required:true
    ontarget "file containing count of on target reads", args: 1
    metrics "metrics output from Picard", args: 1
    bam "Final aligment of sample reads", args:1, required:true
    study "ID of study for which this QC report is being generated", args:1, required: true
    meta "Meta data for the sample to produce the QC report for", args:1, required:true
    threshold "Coverage threshold for signalling a region as having satisfactory coverage", args:1, required: true
    classes "Percentages of bases and corresponding classes of regions in the form: GOOD:95:GREEN,PASS:80:ORANGE,FAIL:0:RED", args:1, required:true
    exome "BED file containing target regions for the whole exome", args:1, required: true
    gc "Gene categories indicating the categories of genes in the cohort / target / flagship", args:1, required:true
    capture "Statistics for percent gene covered by capture", args:1, required:true
    o    "Output file name (PDF format)", args: 1, required:true
    anonymous    "Don't included study ID in PDF", args: 0, required:false
}

opts = cli.parse(args)
if(!opts)
  System.exit(1)

err = { System.err.println "\nERROR: $it\n"; System.exit(1) }
if(!new File(opts.cov).exists()) 
  err "The input coverage file ($opts.cov) does not exist or could not be accessed"

if(!opts.threshold.isInteger())
  err "Please provide the coverage threshold as an integer between 1 and 1000"

int coverageThreshold = opts.threshold.toInteger()

if(!new File(opts.meta).exists())
  err "The sample meta data file ($opts.meta) does not exist or could not be accessed"

// try both metadata formats
Map samples;
try {
  samples = SampleInfo.parse_mg_sample_info(opts.meta)
}
catch (RuntimeException e) {
  samples = SampleInfo.parse_sample_info(opts.meta)
}

if(!samples.containsKey(opts.study))
  err "The provided meta data file ($opts.meta) did not contain meta information for study $opts.study"

SampleInfo meta = samples[opts.study]

// Read the gene categories for the target / cohort / flagship
// geneCategories = new File(opts.gc).readLines()*.split('\t').collect { [it[0],it[1]] }.collectEntries()
geneCategories = new File(opts.gc).readLines().findAll( { it =~ /^[^#]/ } )*.split('\t').collect { [it[0],it[1]] }.collectEntries()
geneCapture = new File(opts.capture).readLines().findAll( { it =~ /^[^#]/ } )*.split('\t').collect { [it[0],it[1]] }.collectEntries()

// Update gene categories with sample specific data
if(meta.geneCategories) {
    meta.geneCategories.each { gene, category ->
        geneCategories[gene] = category
    }
}

String onTargetCapture = "Unknown"
String offTargetCapture = "Unknown"
int totalReads = -1
int mappedPairedReads = -1
int unmappedReads = -1
if(opts.ontarget && opts.ontarget != "") {
  String onTargetText = new File(opts.ontarget).text
  if (onTargetText != "") {
    int onTargetCount = onTargetText.toInteger()
    if(opts.metrics) {
        Map metrics = PicardMetrics.parse(opts.metrics)
        int totalCount = metrics.READ_PAIRS_EXAMINED.toInteger() * 2
        float onTargetPerc = ((float)onTargetCount / ((float)totalCount))
        onTargetCapture = String.format("%2.2f%%", onTargetPerc * 100)
        offTargetCapture = String.format("%2.2f%%", 100 - onTargetPerc * 100)
        unpairedReads = metrics.UNPAIRED_READS_EXAMINED.toInteger()
        mappedPairedReads = (metrics.READ_PAIRS_EXAMINED.toInteger()) * 2
        unmappedReads = metrics.UNMAPPED_READS.toInteger()
        totalReads = unmappedReads + mappedPairedReads + unpairedReads
    }
    else {
        onTarget = String.valueOf(onTargetCount)
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Parse user specified classes into colors and levels
//
/////////////////////////////////////////////////////////////////////////
def classes
try {
  classes = opts.classes.split(",")*.trim().collect { it.split(":") }.collect { 
    [it[0], it[1].toInteger(),it[2]] 
  }.sort{ -it[1] }
}
catch(Exception e) {
  err "The class string $opts.classes couldn't be parsed. Please use the format: <class1>:<percentage>:<color>,<class2>:<percentage>..."
}

println "Percentage thresholds are: $classes"

println "Meta info = $meta"


/////////////////////////////////////////////////////////////////////////
//
// Compute Karyotype
//
/////////////////////////////////////////////////////////////////////////    
BED exomeBed = new BED(opts.exome).load()
SAM sam = new SAM(opts.bam)
SexKaryotyper karyotyper = new SexKaryotyper(sam, exomeBed) 
Utils.time("Running Karyotyping") { 
    karyotyper.run()
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate Coverage Statistics
//
/////////////////////////////////////////////////////////////////////////    
CoverageStats stats = null
int totalBP = 0
int totalOK = 0
String currentGene = null

List<Map> geneReport = []

CoverageStats allGeneStats = new CoverageStats(1000)

def stream = new FileInputStream(opts.cov)
if(opts.cov.endsWith(".gz")) 
    stream = new java.util.zip.GZIPInputStream(stream) 

ProgressCounter.withProgress {
  stream.eachLine { line ->
    def (chr,start,end,gene,offset,cov) = line.split("\t")
    (start,end,offset,cov) = [start,end,offset,cov]*.toInteger()
    if(gene != currentGene) {
      if(stats != null) {
        def geneSummary = [
              gene: currentGene, 
              fracOK: totalOK / (float)totalBP,
              totalOK: totalOK,
              totalBP: totalBP,
              median: stats.median,
              stats: stats
            ]

        geneReport.add(geneSummary)
      }
      
      def prev = geneReport.find { it.gene == gene }
      if(prev) { 
        totalOK = prev.totalOK; totalBP = prev.totalBP; stats = prev.stats; 
        geneReport.removeAll { it.gene == gene }
      } else {
          totalOK = 0; totalBP = 0; stats = new CoverageStats(1000) 
      }
      currentGene = gene
    }
    if(cov > coverageThreshold) 
      ++totalOK

    stats.addValue(cov)
    allGeneStats.addValue(cov)

    ++totalBP
    count()
  }
  if(stats != null) {
    def geneSummary = [
          gene: currentGene, 
          fracOK: totalOK / (float)totalBP,
          totalOK: totalOK,
          totalBP: totalBP,
          median: stats.median,
          stats: stats
        ]

    geneReport.add(geneSummary)
  }
}

// Sort the gene report by category
if(geneCategories) {
    geneReport.sort { String.format( "%d%s", geneCategories[it.gene] ? 9 - geneCategories[it.gene].toInteger() : 9, it.gene ) } // no category is lowest; highest cat first
}

/////////////////////////////////////////////////////////////////////////
//
// Generate PDF
//
/////////////////////////////////////////////////////////////////////////
// Read the coverage file and write out a line in the PDF for each gene containing how many 
// bp are below the threshold
new PDF().document(opts.o) {

  if (opts.anonymous) {
    title("Sequencing Summary Report for Study XXX")
  }
  else {
    title("Sequencing Summary Report for Study $opts.study")
  }

  p("For definitions of calculations, refer to the section at the end of this document.")

  br()

  bold { p("Summary Data") }

  def hd = { text -> bg("#eeeeee") { bold { cell(text) } } }
  table(cols:2,padding:4) {
    hd("Batch"); cell(meta.batch);
    if (opts.anonymous) {
      hd("Study ID"); cell("XXX");
    }
    else {
      hd("Study ID"); cell(meta.sample);
    }
    hd("Sex"); cell(meta.sex);

    hd("Inferred Sex"); 
    if(meta.sex != karyotyper.sex) {
        color("red") { cell(karyotyper.sex.name()) }
    }
    else
        cell(karyotyper.sex.name()) 

    hd("Disease Cohort"); cell(meta.target);
    hd("Hospital / Institution"); cell(meta.institution);
    hd("Ethnicity"); cell(meta.ethnicity);
    hd("Prioritized Genes"); cell(meta.geneCategories.collect { it.key }.join(","));
    hd("Consanguinity Status"); cell(meta.consanguinity.name().replaceAll("_"," "));
    hd("Sample Type (tumor/normal)"); cell(meta.sampleType);
    hd("Sequencing Dates"); cell(meta.sequencingDates*.format("yyyy-MM-dd")?.unique()?.join(", "));
    hd("DNA Collection Dates"); cell(meta.dnaDates*.format("yyyy-MM-dd")?.unique()?.join(", "));
    hd("Sequencing Machines"); cell(meta.machineIds?.unique()?.join(","));
  }
  br()

  bold { p("Coverage Summary") }
  table(cols:2,padding:4) {
    hd("Reported Mean Coverage"); cell(meta.meanCoverage);
    hd("Observed Mean Coverage"); cell(String.format("%2.2f",allGeneStats.mean));
    hd("Observed Median Coverage"); cell(allGeneStats.median);
    hd("Total Reads"); cell(totalReads);
    hd("Unmapped Reads (% of total)"); cell( String.format( "%d (%2.2f%%)", unmappedReads, 100.0 * unmappedReads / totalReads ) );
    hd("Mapped Paired Reads (% of total)"); cell( String.format( "%d (%2.2f%%)", mappedPairedReads, 100.0 * mappedPairedReads / totalReads ) );
    hd("% Mapped On Target (Off Target)"); cell( String.format( "%s (%s)", onTargetCapture, offTargetCapture ) );
  }

  br()

  bold { p("Gene Summary") }
  table {
    bg("#eeeeee") { head {
      cells("Gene", "Category", "Perc > ${coverageThreshold}X","Median", "OK?", "% in Capture")
    } }

    for(geneSummary in geneReport) {

        if(geneSummary.gene.startsWith("Intergenic"))
            continue

        // gene
        cell(geneSummary.gene)
        // category
        cell(geneCategories[geneSummary.gene]?:"")
        // percentage over threshold, median
        align("center") {
          cells(String.format("%2.1f%%", 100 * geneSummary.fracOK), geneSummary.median)
        }

        // ok?
        // Color depends on class
        def clazz = classes.find { geneSummary.fracOK * 100 >= it[1] }
        if(clazz != null) {
            color(clazz[2]) {
              cell(clazz[0])
            }
        }
        else {
            color("RED") { cell("FAIL (*)") } // Less than any provided category, assume fail
        }

        // % in capture
        cell( String.format( "%2.1f%%", (double) Double.parseDouble(geneCapture[geneSummary.gene.toLowerCase()]) ?: -100.0 ) )
    }
  }

  br()
  bold { p("Coverage Summary Definitions") }
  
  p("Reported Mean Coverage: the mean coverage reported by the sequencing lab")
  p("Observed Mean Coverage: the mean coverage across the disease cohort region")
  p("Observed Median Coverage: the median coverage across the disease cohort region")
  p("Total Reads: the total number of reads generated by the sequencer")
  p("Unmapped Reads: reads that were not mapped to the genome")
  p("Mapped Paired Reads: paired reads that were mapped to the genome")
  p("% Mapped on Target: the proportion of mapped reads that have any part align to any part of the capture region")

  br()
  bold { p("Gene Summary Definitions") }
  p("Perc: the percentage of the gene overlapping the capture region with acceptable coverage")
  p("Median: the median coverage across the gene overlapping the capture region")
  p("% in capture: the proportion of the gene that overlaps the capture region")
  
}

// Write out the karyotyping statistics for later reference
new File(opts.o.replaceAll('.pdf','.karyotype.tsv')).text = [
        'Sex' : meta.sex,
        'Inferred Sex' : karyotyper.sex.name(),
        'xCoverage' : karyotyper.xCoverage.mean,
        'yCoverage' : karyotyper.yCoverage.mean,
        'autosomeCoverage' : karyotyper.autosomeCoverage.mean
    ].collect { [it.key, String.valueOf(it.value)].join('\t') }.join('\n') 

