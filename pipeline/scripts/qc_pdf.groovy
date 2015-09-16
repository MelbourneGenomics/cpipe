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
    o    "Output file name (PDF format)", args: 1, required:true
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
geneCategories = new File(opts.gc).readLines()*.split('\t').collect { [it[0],it[1]] }.collectEntries()

// Update gene categories with sample specific data
if(meta.geneCategories) {
    meta.geneCategories.each { gene, category ->
        geneCategories[gene] = category
    }
}

String onTarget = "Unknown"
String totalReads = "Unknown"
if(opts.ontarget) {
    int onTargetCount = new File(opts.ontarget).text.toInteger()
    if(opts.metrics) {
        Map metrics = PicardMetrics.parse(opts.metrics)
        int totalCount = metrics.READ_PAIRS_EXAMINED.toInteger() * 2
        float onTargetPerc = ((float)onTargetCount / ((float)totalCount))
        onTarget = String.format("%2.1f",onTargetPerc)
        totalReads = (metrics.READ_PAIRS_EXAMINED.toInteger())*2
    }
    else {
        onTarget = String.valueOf(onTargetCount)
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
    geneReport.sort { geneCategories[it.gene] ? geneCategories[it.gene].toInteger() : -1 } // no category is lowest
    geneReport = geneReport.reverse() // highest category first
}

/////////////////////////////////////////////////////////////////////////
//
// Generate PDF
//
/////////////////////////////////////////////////////////////////////////
// Read the coverage file and write out a line in the PDF for each gene containing how many 
// bp are below the threshold
new PDF().document(opts.o) {

  title("Sequencing Summary Report for Study $opts.study")

  bold { p("Summary Data") }

  def hd = { text -> bg("#eeeeee") { bold { cell(text) } } }
  table(cols:2,padding:4) {
    hd("Batch"); cell(meta.batch);
    hd("Study ID"); cell(meta.sample);
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
    hd("Prioritzed Genes"); cell(meta.geneCategories.collect { it.key }.join(","));
    hd("Consanguinity Status"); cell(meta.consanguinity.name().replaceAll("_"," "));
    hd("Sample Type (tumor/normal)"); cell(meta.sampleType);
    hd("Sequencing Dates"); fontSize(10) { cell(meta.sequencingDates*.format("yyyy-MM-dd")?.unique()?.join(", ")); }
    hd("DNA Collection Dates"); fontSize(10) { cell(meta.dnaDates*.format("yyyy-MM-dd")?.unique()?.join(", ")); }
    hd("Sequencing Machines"); cell(meta.machineIds?.unique()?.join(","));
  }
  br()

  bold { p("Coverage Summary") }
  table(cols:2,padding:4) {
    hd("Reported Mean Coverage"); cell(meta.meanCoverage);
    hd("Observed Mean Coverage"); cell(String.format("%2.2f",allGeneStats.mean));
    hd("Observed Median Coverage"); cell(allGeneStats.median);
    hd("Total Reads"); cell(totalReads);
    hd("Fraction on Target"); cell(onTarget);
  }

  br()

  bold { p("Gene Summary") }
  table {
    bg("#eeeeee") { head {
      cells("Gene", "Category", "Perc > ${coverageThreshold}X","Median", "OK?")
    } }

    for(geneSummary in geneReport) {

        if(geneSummary.gene.startsWith("Intergenic"))
            continue

        cell(geneSummary.gene)
        cell(geneCategories[geneSummary.gene]?:"")
        align("center") {
          cells(String.format("%2.1f%%",100*geneSummary.fracOK), geneSummary.median)
        }

        // Color depends on class
        def clazz = classes.find { geneSummary.fracOK*100 >= it[1] }
        if(clazz != null) {
            color(clazz[2]) {
              cell(clazz[0])
            }
        }
        else {
            color("RED") { cell("FAIL (*)") } // Less than any provided category, assume fail
        }
    }
  }
}

// Write out the karyotyping statistics for later reference
new File(opts.o.replaceAll('.pdf','.karyotype.tsv')).text = [
        'Sex' : meta.sex,
        'Inferred Sex' : karyotyper.sex.name(),
        'xCoverage' : karyotyper.xCoverage.mean,
        'yCoverage' : karyotyper.yCoverage.mean,
        'autosomeCoverage' : karyotyper.autosomeCoverage.mean
    ].collect { [it.key, String.valueOf(it.value)].join('\t') }.join('\n') 

