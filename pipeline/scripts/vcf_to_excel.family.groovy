// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:cindent:nowrap
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

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

help = """
--------------------------------------------------------------------
Produces an excel report containing a single Excel sheet per VCF file 
provided. Each VCF file may contain multiple samples, in which case 
the samples are each represented as a column within the same excel
sheet / tab. This allows for comparison between samples for family
based diagnosis.
--------------------------------------------------------------------
"""

// Parse command line args
Cli cli = new Cli(usage: "vcf_to_excel.groovy [-s <sample order>] [...csv / vcf / exoncoverage.txt files...]", header: help)

cli.with {
  s 'define order in which samples should appear in spreadsheet', args:1
  p 'define a list of phenotypically relevant genes to highlight (first column used in tab delimited file)', args:1
  l 'linkify elements (only works for small data sets)'
  unique 'Only output the most significant effect (1 row) for each variant'
  g 'annotate coverage over genes'
  meta 'meta data file for samples. If provided it will be used to resolve pedigree information', args:1
  ped 'Pedigree file in PED format. If provided, will be used to annotate affected individuals', args:1
  mindp 'filter out variants when all samples have coverage depth < this value', args:1
  mingq 'filter out variants when all samples have genotype quality < this value', args:1
  db    'database containing historical variant observations. If provided, frequency information within the database will be annotated', args:1
  o     'output file name (results.xlsx)', args:1
}

opts = cli.parse(args)
if(!opts) {
   cli.usage()
   err "Failed to parse command line options"
}
args = opts.arguments()

// Quick and simple way to exit with a message
def err(msg) {
  System.err.println("ERROR: " + msg)
  System.exit(1)
}

LJB_FIELDS = [ "SIFT_score", "SIFT_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "GERP++_RS"]


// Note: 'AAChange' comes from Annovar but is handled separately
def EXPORTED_ANNOVAR_COLUMNS =
        [ 'phastConsElements46way','genomicSuperDups','esp5400_all','1000g2014oct_all','ExAC_Freq','SIFT_score'] +
        LJB_FIELDS

// Genes that will be highlighted in spreadsheet as linked to 
// phenotype of interest
def phenotype_genes = [] as Set
if(opts.p) {
  // PHENOTYPE_GENES = "/home/simons/work/trio/dsdgenes.txt"
  if(!new File(opts.p).exists()) 
        err "Failed to find phenotype gene file $opts.p"
  phenotype_genes = new File(opts.p).readLines().collect { it.split("\t")[0] } as Set
}

List affected = []
if(opts.ped)  {
    List pedigrees = Pedigree.parse(opts.ped)*.value
    affected = pedigrees.collect { it.affected }.flatten().unique()
}

int minDP = 3
if(opts.mindp) 
    minDP = opts.mindp.toInteger()

float minGQ = 5.0f
if(opts.mingq) 
    minGQ = opts.mingq.toFloat()

VariantDB variantDB = null
if(opts.db)
    variantDB = new VariantDB(opts.db)

if(args.size()<2) {
	cli.usage()
    System.exit(1)
}

csvs = args.grep { it.endsWith('.csv') }
vcfs = args.grep { it.endsWith('.vcf') }

fileNameParser = new IlluminaFileNameParser()

outputFileName = "results.xlsx"
if(opts.o)
	outputFileName = opts.o

allSamples = []

annovar_summary = csvs[0]

new ExcelBuilder().build {

  vcfs.each { vcf_file ->

        println "Reading Annovar summary $annovar_summary (vcf file = $vcf_file) ..."

        // Read the exome annotations from annovar
        def annovar = new CSV(annovar_summary)

        // println "Read ${annovar_variants.size()} annotations"
        def sampleName = "unknown_sample"
        VCF vcfHeader = new VCF(vcf_file)
        sampleName = vcfHeader.samples[0]

        println "Sample name is $sampleName"

        def exportedColumns 

        def sht = sheet(sampleName) {

            // All lines, including ones we ignored
            int lineCount = 0

            // Count of lines we actually wrote
            int lineIndex = 0

            // Current line being processed in current file
            def currentLine  = "Unknown"

            try {
                // Find the line with the column names. We're doing this to figure out the sample names
                println "Processing $vcf_file"

                def sampleNames = vcfHeader.samples
                def exportSamples = sampleNames
                if(opts.s) {
                    exportSamples = opts.s.split(",").toList()
                    println "User defined export order: $exportSamples"
                    if(!exportSamples.every { it in sampleNames })
                        throw new Exception("Could not find one or more samples $opts.s in sample names from VCF file: $sampleNames")
                }
                else {
                    exportSamples = vcfHeader.samples
                }


                println "Exporting samples $exportSamples from $vcf_file"

                // If the file contains CLR (samtools constraint likelihood ratios), add a column for those

                // Write the header line and make it bold
                exportedColumns = ['gene','chr','start','end','id','rank','effect','qual','depth'] + 
                    (exportSamples.size()>1?["scount"]:[]) +
                    (affected.size()>1?["acount"]:[]) +
                    (vcfHeader.hasInfo("FC")?["mut fm cnt"]:[]) +
                    (vcfHeader.hasInfo("GC")?["gene fm cnt"]:[]) +
                    (opts.db?["db cnt","db fam cnt"]:[]) +
                    exportSamples*.replaceAll('_$','') + 
                    (exportSamples.size()!=sampleNames.size()?["allcount"]:[]) +
                    (vcfHeader.hasInfo("CLR")?["CLR"]:[]) +
                    (vcfHeader.hasInfo("SSDNP")?["SSDNP"]:[]) +
                    (phenotype_genes ? ["pheno match"] : []) +
                    [ 'Length' ] +
                    [ 'DNAChange' ] +
                    ['Annotation'] +
                    EXPORTED_ANNOVAR_COLUMNS 


                allSamples.addAll(exportSamples*.replaceAll('_$','')) 

                row { bold {
                    if(opts.sex) { 
                      bottomBorder { cells(exportedColumns) }
                    }
                    else {
                      cells(exportedColumns)
                    }
                  }
                }

                // Reasons why variants filtered out
                Map reasons = [minGQ : 0, minDP : 0, dosage:0]

                VCFIndex vcfIndex = new VCFIndex(vcf_file)
                annovar.each { variant ->

                  Map vInfo = vcfIndex.findAnnovarVariant(variant.Chr, variant.Start, variant.End, variant.Alt)
                  if(!vInfo) {
                      println "WARNING: Could not find Annovar variant $variant.Chr:$variant.Start-$variant.End $variant.Ref/$variant.Alt from $annovar_summary in VCF file $vcf_file"
                      return
                  }

                  Variant v = vInfo.variant
                  Allele allele = v.alleles[vInfo.allele]

                  currentLine = v.line // for debugging, save this so the exception can see it

                  ++ lineCount

                  def snpEffInfo 
                  if(opts.unique) {
                      snpEffInfo = [v.maxEffect]
                  }
                  else
                      snpEffInfo = v.snpEffInfo

                  genes = v.snpEffInfo*.gene
                  effs = v.snpEffInfo*.type
                  ranks = v.snpEffInfo*.impact

                  if(genes == null)
                      err "No SnpEFF annotations are available on the variants in VCF file $vcf_file"

                  if(lineCount % 1000 == 0)
                      println lineCount

                  String urlPos = v.pos
                  String urlChr = URLEncoder.encode(v.chr.replaceAll('chr',''))

                  // For trio calling, Samtools outputs CLR in the INFO field.
                  // We would like to extract it to a separate column, IF it exists 
                  def clr = vcfHeader.hasInfo("CLR") ? v.info.CLR : "0"

                  def dosages =exportSamples.collect { v.sampleDosage(it, vInfo.allele) } 
                  if(dosages.every { it == 0 }) {
                      reasons.dosage++
                      return false
				  }

                  def depths = exportSamples.collect { v.genoTypes[sampleNames.indexOf(it)] }*.AD*.sum()*.toInteger().collect { it ?: 0 }
                  if(depths.every { it < minDP }) {
                    reasons.minDP++
                    return false
                  }

                  def gqs = exportSamples.collect { v.genoTypes[sampleNames.indexOf(it)] }*.GQ.collect { it ?: 0.0f }
                  if(gqs.every { it < minGQ }) {
                    reasons.minGQ++
                    return false
                  }

                  // There are places where there are overlapping genes
                  // To enable easy filtering in the spreadsheet, we actually 
                  // write out the same row for each gene
                  for(geneAndEffect in [genes,effs,ranks].transpose().unique()) {

                      def (gene,effect,rank) = geneAndEffect
                      // if(!(rank in ['HIGH','MODERATE']))
                      //    return

                      def row = row {}
                      def genesCell = row.addCell(gene)
                      if(opts.l && v.id == '.') {
                            genesCell.link("http://www.genecards.org/cgi-bin/carddisp.pl?gene=${gene}&search=${gene}")
                            if(phenotype_genes.contains(gene))
                                genesCell.red()
                      }

                      // if(v.id == '.')
                      //    genesCell.link("http://asia.ensembl.org/Homo_sapiens/Gene/Phenotype?g=${URLEncoder.encode(genes[0])};r=$urlChr:$urlPos")

                      // search pubmed
                      // row.addCell(genes.join(',')).link("http://www.ncbi.nlm.nih.gov/pubmed?term=${genes[0]}")

                      row.addCell(v.chr)

                      // UCSC
                      // row.addCell(v.pos).link("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$v.chr%3A$v.pos&dgv=pack&knownGene=pack&omimGene=pack")

                      boolean isDeletion = v.ref.size() > v.alt.size()
                      boolean isInsertion = v.alt.size() > v.ref.size()
                      int len = (allele.alt.size() == v.ref.size()) ? 1 : Math.abs(allele.alt.size() - v.ref.size())

                      int start = v.pos
                      int end = v.pos
                      if(isDeletion) {
                          start = v.pos+1
                          end = start+len
                      }

                      def posCell = row.addCell(start)

                      // Length
                      row.add(v.pos + len)

                      // ENSEMBL
                      // if(v.id == '.')
                      //    posCell.link("http://asia.ensembl.org/Homo_sapiens/Location/View?g=${URLEncoder.encode(genes[0])};r=$urlChr:$urlPos")

                      row.add(v.id, rank=="HIGH"?2:1, effect, gqs.grep { it > 0 }.min(), depths.grep { it > 0 }.min())
                      // row.add(v.id, rank=="HIGH"?2:1, effect, v.qual, v.info.DP)

                      // If there is more than one sample then it's useful to have the count of 
                      // how many total samples the variant was observed in
                      if(exportSamples.size()>1) {
                          row.add(dosages.count { it > 0 })
                      }

                      if(affected) {
                          row.add(affected.count { v.sampleDosage(it) > 0 })
                      }

                      if(vcfHeader.hasInfo("FC")) {
                          row.addCell(v.info.FC)
                      }

                      if(vcfHeader.hasInfo("GC")) {
                          row.addCell(v.info.GC?:"0")
                      }

                      if(opts.db) {
                          def counts = variantDB.countObservations(v, allele)
                          row.add(counts.samples, counts.families)
                      }

                      dosages.eachWithIndex { dosage,index ->
                          boolean isAffected = affected.contains(exportSamples[index])
                          withStyle { 
                              if((dosage > 0) && ((gqs[index]<minGQ) || (depths[index]<minDP))) {
                                  applyStyle("gray",null) 
                              }

                              if((dosage > 0) && isAffected) {
                                  applyStyle('red',null) 
                              }

                              if(isAffected) {
                                  applyStyle('bgOrange',null)
                              }

                              cell(dosage)
                          }
                      }

                      if(exportSamples.size()!=sampleNames.size()) {
                          row.add(exportSamples.count { v.sampleDosage(it)>0} )
                      }

                      if(vcfHeader.hasInfo("SSDNP")) 
                          row.addCell(v.info.SSDNP)

                      if(vcfHeader.hasInfo("CLR")) {
                          row.addCell(clr?Integer.parseInt(clr):0)
                      }

                      if(phenotype_genes) {
                          row.add(phenotype_genes.contains(gene) ? 1 : 0)
                      }

                      row.add(len) 

                      row.add(v.ref + ' / ' + allele.alt)

                      // The ESP5400 frequencies can sometimes have ridiculously long precisions: round them
                      // to a few digits
                      def esp5400 = variant.esp5400_all
                      if(esp5400 && (esp5400 != ".")) {
						  if(esp5400 instanceof String)
							esp5400 = Float.parseFloat(esp5400)
                          esp5400 = Float.parseFloat(String.format("%2.3f",esp5400))
                      }
                      else {
                          esp5400 = 0f;
                      }

                      def g1000 = 0f;
                      if(variant['1000g2014oct_all'] && variant['1000g2014oct_all'] != ".")
                          g1000 = variant['1000g2014oct_all'] 

                      // def key = v.chr+'_'+v.pos+'_'+v.ref+'_'+v.alt
                      // if(annovar_variants.containsKey(key)) {
                          // def variant = annovar_variants[key]
                          row.add(variant.AAChange) //.replaceAll('^.*?:',''))
                          List exported = EXPORTED_ANNOVAR_COLUMNS.collect { variant[it] }
                          exported[EXPORTED_ANNOVAR_COLUMNS.indexOf('esp5400_all')] = esp5400;
                          exported[EXPORTED_ANNOVAR_COLUMNS.indexOf('1000g2014oct_all')] = g1000;
                          row.add(exported)
                      // }

                      // row.addCell(info.find { it.startsWith('EFF=') }?.substring(4))
                      ++lineIndex
                  }
                  return false
                }
                println "Variants filtered from $sampleName : " + reasons
                // sheet.setAutoFilter(CellRangeAddress.valueOf("A1:AK"+lineIndex))
            }
            catch(Exception e) {
                //  System.err.println "Failed processing at line $lineCount:\n\n$currentLine\n\n"
                System.err.println "Failed processing at line $lineCount"
                throw e
            }
            println "$lineIndex / $lineCount rows were exported for sample $sampleName"

    }

    for(int i=0; i<exportedColumns.size(); ++i) {
        sht.autoSizeColumn((short)i)
    }
  }

}.save(outputFileName)
println "Output is $outputFileName"

