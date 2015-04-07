// vim: shiftwidth=4:ts=4:expandtab
/////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Demonstration Project
//
// Script to augment annovar output with a Condel score, until
// Annovar gets around to including it in their annotations.
//
// This script reads the VCF file and Annovar exome summary file
// and then writes out the Annovar file with one extra column, that
// being the Condel score for the Annovar variant.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////

import com.xlson.groovycsv.*
import au.com.bytecode.opencsv.*

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "merge_condel.groovy [options]\n")
cli.with {
  i 'VCF file containing Condel scores', args:1
  a 'Annovar file containing annotations', args:1
  t 'Trim <integer> columns from end (different versions of Annovar put unwanted columns on the end)', args:1
  o 'Output file (*.csv)', args:1
}
opts = cli.parse(args)

// Quick and simple way to exit with a message
err = { msg ->
  System.err.println("\nERROR: " + msg + "\n")
  cli.usage()
  System.err.println()
  System.exit(1)
}

if(!opts) {
   cli.usage()
   err "Failed to parse command line options"
}

args = opts.arguments()
if(!opts.i) 
    err "Please provide -i option to specify VCF file to process"
if(!opts.o) 
    err "Please provide -o option to specify output file name"
if(!opts.a)
    err "Please provide -a option to specify Annovar annotation file"

trim_columns = opts.t ? opts.t.toInteger() : 0

// First sniff the header line from the Annovar file to get the columns
ANNOVAR_FIELDS = null
new File(opts.a).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }

// Not all the fields have headers (why?)
if(!("Qual" in ANNOVAR_FIELDS))
    ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCFIndex vcf = new VCFIndex(opts.i)

//
// Function to find an Annovar variant in the original VCF file
//
def find_vcf_variant(vcf, av, lineIndex) {
  try {
      int pos = av.Start.toInteger()
      int end = av.End.toInteger()
      return vcf.find(av.Chr,pos-10,end) { variant -> variant.equalsAnnovar(av.Chr,pos,av.Alt) }
  }
  catch(Exception e) {
      try { println "WARNING: unable to locate annovar variant at $lineIndex in VCF ($e)" } catch(Exception e2) { e.printStackTrace() }
  }
}

// Output file
def writer = new FileWriter(opts.o)
writer.println(ANNOVAR_FIELDS.join(",")+",Condel")

CSVWriter csvWriter = new CSVWriter(writer);
def annovar_csv = new CsvParser().parse(new File(opts.a).text, separator:',')
int lineIndex = 0
for(av in annovar_csv) {
    ++lineIndex
    if(lineIndex%5000==0)
        println new Date().toString() + "\tProcessed $lineIndex lines"

    def values = trim_columns ? av.values[0..(trim_columns-1)] : av.values

    def variant = find_vcf_variant(vcf,av,lineIndex)
    if(!variant) {
        println "WARNING: Variant $av.Chr:$av.Start at line $lineIndex could not be found in the original VCF file"
        csvWriter.writeNext((values + [""]) as String[])
        continue
    }

    // Parse out the VEP annotations - the only one we want is Condel
    // def veps = variant.info.CSQ.split(",").collect { csq -> [VEP_FIELDS,csq.split("\\|")].transpose().collectEntries() }
    def veps = variant.vepInfo

    def gene = av.Gene.replaceAll(/\(.*$/,"")
    vep = veps.grep { it.SYMBOL == gene }.max { it.Condel ? it.Condel.toFloat() : 0 }
    if(!vep) {
        println "WARNING: No vep consequence was for the same gene: ${veps*.SYMBOL} vs $av.Gene"
        vep = veps[0]
    }
    csvWriter.writeNext((values + [vep.Condel]) as String[])
}
writer.close()
println "Processed $lineIndex rows"
