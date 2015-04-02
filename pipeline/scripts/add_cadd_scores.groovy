// vim: shiftwidth=4:ts=4:expandtab:cindent
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// Script to augment standard Annovar exome summary output with CADD scores
//
// This script reads both the standard exome summary AND also the 
// CADD score output and then adds a column at the end of the Annovar
// output containing the CADD score for each variant, where one is
// available.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
//
/////////////////////////////////////////////////////////////////////////

import com.xlson.groovycsv.*
import au.com.bytecode.opencsv.*

System.err.println("="*80)
System.err.println("\nCADD Annotation Merge Tool\n")
System.err.println("="*80)
System.err.println()

CliBuilder cli = new CliBuilder(usage: "add_cadd_scores.groovy <options>",writer:new PrintWriter(System.err))
cli.with {
    c "file containing cadd scores (*.hg19_cadd_dropped)", longOpt: "cadd", args:1, required:true
    a "Annovar exome_summary output", longOpt: "annovar", required:true, args:1
    o "Output file (CSV format)", args:1, required:true
}

opts = cli.parse(args)

CADD_COLUMNS = [ "annotation",  "score","chr","start","end","ref","alt","dosage", "qual", "depth" ]
def cadd_csv  = ExcelCategory.parseCSV("",opts.c, CADD_COLUMNS)

// Let's index them by chr:pos (they are all SNVs so this only needs 1 possition)
Map cadd_index = [:]
for(cadd in cadd_csv) {
  cadd_index[cadd.chr + ":" + cadd.start + ":" + cadd.alt] = cadd
}

println "Read ${cadd_index.size()} CADD annotations from Annovar output"

// We need the annovar column headings for later, so we read
// them separately, even though the CSV reader does it too
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

def annovar_csv = ExcelCategory.parseCSV("",opts.a, ',')

ProgressCounter c = new ProgressCounter()
int countFound = 0
int countMissing = 0

Stats scores = new Stats()
Stats quals = new Stats()
Stats depths = new Stats()

new File(opts.o).withWriter { out ->
    // Write out the column headers first, because CSVWriter doesn't do that
    out.println((ANNOVAR_FIELDS + "CADD").join(","))
    CSVWriter csvWriter = new CSVWriter(out);
    for(av in annovar_csv) {

        c.count()

        // Find the CADD score
        // Important Annovar columns are: Chr,Start,End,Ref,Obs
        String key = "$av.Chr:$av.Start:$av.Obs"
        def score = ""
        if(cadd_index.containsKey(key)) {
            ++countFound
            cadd = cadd_index[key]
            score = cadd.score
            scores << cadd.score.toDouble()
            quals << cadd.qual.toDouble()
            depths << cadd.depth.toDouble()
        }
        else {
            ++countMissing
        }
        csvWriter.writeNext((av.values + [score]) as String[])
    }
//c.end()
}

System.err.println "Processed $c.count Annovar variants"
System.err.println "Found $countFound CADD scores (${((float)countFound/(c.count+1))*100}%)"
System.err.println "$countMissing CADD scores were missing (${((float)countMissing/c.count)*100}%)"
System.err.println "Mean CADD score ${scores.mean} (Median=${scores.getPercentile(50)})"
System.err.println "Mean Quality Score for variants with CADD scores ${quals.mean} (Median=${quals.getPercentile(50)})"
System.err.println "Mean Coverage Depth for variants with CADD scores ${depths.mean} (Median=${depths.getPercentile(50)})"
