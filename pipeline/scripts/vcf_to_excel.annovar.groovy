// vim: shiftwidth=4:ts=4:expandtab:cindent:number
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
// VCF to Excel Format Conversion
//
// This script reads the VCF file and Annovar exome summary file
// and produces an Excel file that can be easily viewed and filtered
// by less technical users.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////

import groovy.sql.Sql
import com.xlson.groovycsv.*
import au.com.bytecode.opencsv.*
import org.apache.commons.cli.Option

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  s 'comma separated list of samples to include', args:1, required:true
  vcf 'VCF file to convert to Excel format', args: Option.UNLIMITED_VALUES, required:true
  a 'Annovar file containing annotations', args: Option.UNLIMITED_VALUES, required:true
  o 'Name of output file', args:1, required:true
  x 'Comma separated list of functional types to exclude', args:1
  si 'sample meta data file for the pipeline', args:1, required:true
  db 'Sqlite database containing known variants. If known, a column will be populated with the count of times observed.', args:1
  gc 'File listing genes and categories', args:1, required:true
  oocf 'Count at which out-of-cohort variants will cause variants to be filtered from report', args:1, required: true
  pgx 'VCF file containing variants to treat as pharmacogenomic variants (always report)', args:1
  bam 'BAM file for annotating coverage depth where not available from VCF files', args: Option.UNLIMITED_VALUES
  pgxcov 'Coverage threshold below which a pharmocogenomic site is considered untested (15)', args: 1
  annox 'Directory to send Annovar style per-sample summaries to', args: 1, required: true
  xprof 'Analysis profiles to exclude from contributing variant counts in variant filtering by internal database', args:1
  log 'Log file for writing information about variants filtered out', args: 1
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

Writer log = new File(opts.log?:'/dev/null').newWriter()

int out_of_cohort_variant_count_threshold = opts.oocf.toInteger() 

def pg_variants = []
if(opts.pgx) 
    pg_variants = VCF.parse(opts.pgx)

int pgx_coverage_threshold = opts.pgxcov ? opts.pgxcov.toInteger() : 15

sample_info = SampleInfo.parse_sample_info(opts.si)

// println "sample_info = $sample_info"




exclude_types = opts.x ? opts.x.split(",") : []
excluded_profiles_from_counts = opts.xprof ? opts.xprof.split(",") as List : ["AML"]

samples = opts.s.split(",")

Map<String,SAM> bams = null
if(opts.bams) {
    bams = opts.bams.collectEntries { def bam = new SAM(it); [ bam.samples[0], bam ] }
    println "="*80
    println "Read ${bams.size()} BAM files for querying read depth:"
    bams.each {  println "Sample: $it.key => $it.value.samFile " }
    println "="*80
}

// Read the header from the first annovar file to find out the column names
ANNOVAR_FIELDS = null
new File(opts.as[0]).withReader { r -> ANNOVAR_FIELDS = r.readLine().split(",") as List }
// Not all the fields have headers (why?)
if(!("Qual" in ANNOVAR_FIELDS))
    ANNOVAR_FIELDS += ["Qual","Depth"]

println "Annovar fields are " + ANNOVAR_FIELDS

// Then read all the annovar files
println "Processing ${opts.as.size()} Annovar files"

// connect to database, if specified
db = null
if(opts.db) 
    db = new VariantDB(opts.db)

// Read the gene categories
geneCategories = new File(opts.gc).readLines()*.split('\t').collect { [it[0],it[1]] }.collectEntries()

AACHANGE_FIELDS = ANNOVAR_FIELDS.grep { it.startsWith("AAChange") }

EXAC_FIELDS=["exac03","ExAC_ALL","ExAC_Freq"]

EXAC_FIELD = EXAC_FIELDS.find { it in ANNOVAR_FIELDS }
if(EXAC_FIELD == null) 
    EXAC_FIELD = "exac03"

ONEKG_FIELD="1000g2014oct_all"

ESP_FIELD = ANNOVAR_FIELDS.find { it =~ /^esp/ }
if(ESP_FIELD == null)
    ESP_FIELD = "esp5400_all"

LJB_FIELDS = [ "SIFT_score", "SIFT_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "GERP++_RS", "phyloP100way_vertebrate"]

// The LJB fields changed column headings. To preserve backwards compatibility with
// downstream scripts, we replace them with the old headings in output
OLD_LJB_FIELDS = [ "LJB_SIFT", "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred", "LJB_GERP++", "LJB_PhyloP"]

// Order preferred if clinicians need to review output directly
OUTPUT_FIELDS = ["Func", "Gene", "ExonicFunc"] + 
                AACHANGE_FIELDS + 
                ["Gene Category", "Priority_Index", "CADD_raw", "CADD_phred", "Condel", "phastConsElements46way", ESP_FIELD, ONEKG_FIELD, "snp138", EXAC_FIELD] +
                LJB_FIELDS + 
                [ "genomicSuperDups", "Chr", "Start", "End", "Ref", "Alt", "Otherinfo", "Qual", "Depth", "#Obs", "RefCount", "AltCount", "PRIORITY_TX"]

OUTPUT_CSV_FIELDS = ["Func","Gene","ExonicFunc"] +
                    AACHANGE_FIELDS + 
                    ["phastConsElements46way","genomicSuperDups",ESP_FIELD,ONEKG_FIELD,EXAC_FIELD,"snp138"] +
                    LJB_FIELDS +
                    ["Chr","Start","End","Ref","Alt","Otherinfo","Qual","Depth","Condel","Priority_Index","CADD_raw","CADD_phred", "Gene Category","Priority_Index","#Obs","RefCount","AltCount","PRIORITY_TX"]

CENTERED_COLUMNS = ["Gene Category", "Priority_Index", ONEKG_FIELD, ESP_FIELD, "LJB_PhyloP_Pred","LJB_SIFT_Pred","LJB_PolyPhen2","LJB_PolyPhen2_Pred"]

// The output headings are sometimes different to the input headings
// this is done to preserve compatibility as annovar headings change
// occasionally with the software
HEADING_MAP = OUTPUT_FIELDS.collectEntries{[it,it]} + [
   "phastConsElements46way" : "Conserved",
   "esp5400_All" : "ESP5400_ALL",
   "1000g2010nov_all"  : "1000g2010nov_ALL",
   ONEKG_FIELD  : "1000g",
   "snp138" : "dbSNP138",
   "genomicSuperDups" : "SegDup",
   "ExAC_ALL" : "exac03",
   "ExAC_Freq" : "exac03",
   "CADD_raw": "CADD"
] + [ LJB_FIELDS, OLD_LJB_FIELDS ].transpose().collectEntries()

println HEADING_MAP

extractAAChange = { gene, aaChange ->
    if(gene.indexOf("(")>=0) {
        def geneParts = (gene =~ /(.*)\((.*)\)/)[0]
        gene = geneParts[1].toString()
        return geneParts[2].toString()
    }
    else
	return aaChange
}

//
// Utility function to collect information about a variant into the columns
// required for export.
//
collectOutputValues = { lineIndex, funcGene, variant, sample, variant_counts, av ->

    // Build up values for the row in a map with the column name as the key
    def outputValues = [:]

    (func,gene) = funcGene

    for(aaChange in AACHANGE_FIELDS) {
        outputValues[aaChange] = extractAAChange(gene, av[aaChange])
    }
    outputValues.ExonicFunc = func=="splicing"?"":av.ExonicFunc

    def geneCategory = geneCategories[gene]
    if(sample_info[sample].geneCategories[gene])
        geneCategory = sample_info[sample].geneCategories[gene]

    outputValues["Gene Category"] = geneCategory?:1
    
    outputValues["Gene"] = gene
    outputValues["Func"] = func

    for(af in ANNOVAR_FIELDS) {
        if(av.columns.containsKey(af)) {
            outputValues[af] = av[af]
        }
    }

    // New version of Annovar puts het/hom, Qual and Depth all in one tab separated field called Otherinfo
    def otherInfo = av.Otherinfo.split("\t")
    outputValues["Otherinfo"] = otherInfo[0] // the original value that was called Otherinfo
    outputValues["Qual"] = otherInfo[1]
    outputValues["Depth"] = otherInfo[2]

    outputValues.CADD = av.columns.CADD != null ? av.CADD: ""

    if(db) {
        outputValues["#Obs"] = variant_counts.in_target
    }

    outputValues.RefCount=outputValues.AltCount="";
    if(!variant) 
        return outputValues // Cannot annotate allele depths for this variant

    // Try to annotate allele frequencies
    def gt = variant.sampleGenoType(sample)
    if(gt) {
        // Reference depth
        if(gt.containsKey('AD')) {
            outputValues.RefCount=gt.AD[0]

            // Alternate depth depends on which allele
            int altAllele = (variant.alts.size()==1)?1:variant.equalsAnnovar(av.Chr, av.Start.toInteger(), av.Alt)
            outputValues.AltCount = gt.AD[altAllele]
        }
        else {
          System.err.println("WARNING: variant $variant.chr:$variant.pos ($variant.ref/$variant.alt) had no AD info for sample $sample at line $lineIndex")
        }
    }
    else {
      System.err.println("WARNING: variant $variant.chr:$variant.pos ($variant.ref/$variant.alt) had no genotype for sample $sample at line $lineIndex")
    }
    return outputValues
}

// Because excel can only handle up to 30 chars in the worksheet name,
// we may have to shorten them
int sampleNumber = 1
MAX_SAMPLE_NAME_LENGTH=30
sheet_samples = [samples, samples.collect { 
    it.size() > MAX_SAMPLE_NAME_LENGTH ? "S_" + (sampleNumber++) + "_" + it.substring(0,MAX_SAMPLE_NAME_LENGTH-10) : it 
}].transpose().collectEntries()

//
// Now build our spreadsheet, and export CSV in the same loop
//
try {
    new ExcelBuilder().build {

        for(sample in samples) { // one sample per spreadsheet tab
            def s = sheet(sheet_samples[sample]) { 
                lineIndex = 0
                sampleCount = 0
                includeCount=0

                // Read the CSV file entirely
                // Sort the annovar output by Priority Index
                String samplePrefix = sample+"."
                String annovarName = opts.as.find{new File(it).name.startsWith(samplePrefix)}
                if(annovarName == null)
                    err "The following samples did not have an Annovar file provided: $sample in Annovar files:\n${opts.as.join('\n')}"

                println "Processing $annovarName ..."
                def annovar_csv = parseCSV(annovarName,',').grep { it.Priority_Index.toInteger()>0 }.sort { -it.Priority_Index.toInteger() }

                // Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
                String vcfName = opts.vcfs.find { new File(it).name.startsWith(samplePrefix) }
                if(vcfName == null)
                    err "The following samples were not found in the VCF file provided: ${sample}"

                VCFIndex vcf = new VCFIndex(vcfName)

                // Write out header row
                bold { row {
                        cells(OUTPUT_FIELDS.collect {HEADING_MAP[it]} )
                } }

                println "Priority genes for $sample are ${sample_info[sample].geneCategories.keySet()}"

                // We are going to write out a CSV that is identical to the original annovar output
                // but which includes our custom fields on the end
                // Start by writing the headers
                def writer = new FileWriter("${opts.annox}/${sample}.annovarx.csv")
                writer.println(OUTPUT_CSV_FIELDS.collect{HEADING_MAP[it]}.join(","))
                CSVWriter csvWriter = new CSVWriter(writer);
                for(av in annovar_csv) {
                    ++lineIndex
                    if(lineIndex%5000==0)
                        println new Date().toString() + "\tProcessed $lineIndex lines"

                    // note: check for exonic, because splice events show up as synonymous but with 
                    // Func="exonic;splicing", and should not be filtered out this way
                    if(av.ExonicFunc in exclude_types && av.Func=="exonic") { 
                        log.println "Variant $av.Chr:$av.Start-$av.End excluded by being an excluded type: $av.ExonicFunc"
                        continue
                    }

                    def variantInfo = vcf.findAnnovarVariant(av.Chr, av.Start, av.End, av.Alt)
                    if(!variantInfo) {
                        println "WARNING: Variant $av.Chr:$av.Start at line $lineIndex could not be found in the original VCF file"
                        log.println "Variant $av.Chr:$av.Start excluded because it could not be identified in the source VCF file"
                        continue
                    }

                    Variant variant = variantInfo.variant
                    if(variant.sampleDosage(sample, variantInfo.allele)==0)
                        continue

                    Map variant_counts = [total: 0, other_target:0]
                    if(db) {
                        variant_counts = db.queryVariantCounts(variant, 
                                                               variant.alleles[variantInfo.allele], 
                                                               sample, 
                                                               sample_info[sample].target, 
                                                               excludeCohorts: excluded_profiles_from_counts,
                                                               batch: sample_info[sample].batch)
                        if(variant_counts.other_target>out_of_cohort_variant_count_threshold) {
                            log.println "Variant $variant excluded by presence ${variant_counts.other_target} times in other targets"
                            continue
                        }
                    }
                    ++includeCount

                    ++sampleCount
                    def funcs = av.Func.split(";")
                    def genes = av.Gene.split(";")

                    [funcs,genes].transpose().each { funcGene ->
                        
                        def outputValues = collectOutputValues(lineIndex, funcGene, variant, sample, variant_counts, av)

                        // Write the row into the spreadsheet
                        row {
                            OUTPUT_FIELDS.each { fieldName ->
                                if(fieldName in CENTERED_COLUMNS) { 
                                    center {cell(outputValues[fieldName])}
                                }
                                else {
                                    cell(outputValues[fieldName]) 
                                }
                            }
                      }

                      // Write Annovar CSV format
                      csvWriter.writeNext( OUTPUT_CSV_FIELDS.collect { fieldName ->
                        outputValues[fieldName] == null ? "" : outputValues[fieldName]
                      } as String[])
                  }
                } // End Annovar variants

                // Now add pharmacogenomic variants
                for(pvx in pg_variants) {

                    values = OUTPUT_FIELDS.collect { "" }

                    // Check if the unfiltered VCF has the variant
                    // Note: there's an issue here about canonicalizing the variant
                    // representation. For now, it's being ignored.
                    def vx = vcf.contains(pvx)

                    List<Map> vepInfos = pvx.vepInfo
                    def genes = vepInfos*.SYMBOL.grep { it != null }.join(",")

                    def state = "Untested"
                    int depth = bams[sample].coverage(pvx.chr, pvx.pos)
                    // println "Queried depth $depth at $pvx.chr:$pvx.pos"
                    // TODO: why is vx sometimes null? should always be genotyped
                    if(vx && depth >= pgx_coverage_threshold) {
                        int allele = vx.findAlleleIndex(pvx.alleles[0])
                        state = vx.sampleDosage(sample, allele) > 0 ? "Present" : "Absent"
                    }
                    else { // Variant not called - but is there coverage?
                        if(depth >= pgx_coverage_threshold)
                            state = "Absent"
                        System.err.println "WARNING : PGX variant $pvx was not genotyped for sample $sample"
                    }

                    // Convert to annovar form since we are using Annovar annotations in the 
                    // rest of the report
                    def annovarVx = vx ? vx.toAnnovar() : pvx.toAnnovar()

                    // Now set the fields that we can
                    def output = [ 
                        'Gene Category': 1, 
                        'Priority Index': 1, 
                        Func: "pharma", 
                        ExonicFunc: state,
                        Gene: genes,
                        snp138: pvx.id,
                        Chr: pvx.chr,
                        Start: pvx.pos,
                        End: pvx.pos + pvx.size(),
                        Ref: annovarVx.ref,
                        Obs: annovarVx.obs,
                        Otherinfo: vx ? (vx.dosages[0] == 1 ? "het" : "hom") : ""
                    ]
                    output.each { k, v ->
                        values[OUTPUT_FIELDS.indexOf(k)] = v 
                    }

                    csv_out = []
                    nvlcell = { cell(it == null ? "" : it ) }
                    row { 
                        OUTPUT_FIELDS.each { fieldName ->
                            //println "Export $fieldName = ${outputValues[fieldName]}"
                            if(fieldName in CENTERED_COLUMNS) { 
                                center { nvlcell(output[fieldName]) } 
                            }
                            else
                            if(fieldName == "ExonicFunc" && state == 'Present') {
                                red {nvlcell(output[fieldName])}  
                            }
                            else {
                                nvlcell(output[fieldName]) 
                            }
                        }
                    }
                    csvWriter.writeNext( ANNOVAR_FIELDS.collect { f -> if(output[f] != null) { String.valueOf(output[f]) } else "" } as String[] )

                }
                csvWriter.close()
            }
            println "Sample $sample has ${sampleCount} / ${includeCount} of included variants"
            try { s/*.autoFilter("A:"+(char)(65+6+samples.size()))*/.autoSize() } catch(e) { println "WARNING: Unable to autosize columns: " + String.valueOf(e) }
            s.setColumnWidth(OUTPUT_FIELDS.indexOf("Gene"),60*256) // 30 chars wide for Gene column
            AACHANGE_FIELDS.each { aaChange -> 
                s.setColumnWidth(OUTPUT_FIELDS.indexOf(aaChange),30*256) // 60 chars wide for AAChange column
            }
            // s.setColumnWidth(OUTPUT_FIELDS.indexOf("AAChange_RefSeq"),30*256) // 60 chars wide for AAChange column
            // s.setColumnWidth(OUTPUT_FIELDS.indexOf("AAChange_UCSC"),30*256) // 60 chars wide for AAChange column
            s.setColumnWidth(OUTPUT_FIELDS.indexOf("Gene Category"),14*256) // 14 chars wide for Gene category column
        }

        sheet("README") {
            row { }
            row { cell("This sheet contains explanations of columns in the previous sheet(s)") }
            row {}
            row { cell("Gene").bold(); cell("The gene affected. A mutation may occur on multiple rows if more than one gene or transcript is affected") }
            row { cell("ESP5400").bold(); cell("Frequency of allele in ESP project (5400 exomes)") }
            row { cell("1000g2010nov_all").bold(); cell("Frequency of allele in 1000 Genomes project 2010 Nov release") }
            row { cell("LJB_XXX").bold(); cell("DBNSFP annotations indicating predictions of pathogenicity") }
            row { cell(""); cell("Numeric: 0 = low impact, 1.0 = high impact") }
            row { cell(""); cell("D=Damaging") }
            row { cell(""); cell("P=Probably Damaging") }
            row { cell(""); cell("T=Tolerated") }
            row { cell(""); cell("N=Neutral") }
            row { cell(""); cell("B=Benign") }
            row { cell(""); cell("See link below for more information") }
            row { cell(""); cell("http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b4.readme.txt").link("http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b4.readme.txt") }
            row {}
            if(opts.x) {
                    row{ cell("NOTE:").bold(); cell("The following categories of variant are excluded from this spreadsheet:")}
                    row{ cell(""); cell( opts.x ) }
            }
        }.autoSize()
    }.save(opts.o)
}
finally {
    if(db)
       db.close()

    log.close()
}
