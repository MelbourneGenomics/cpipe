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
// This script reads the VCF file containing annotations by VEP 
// and produces an Excel file that can be easily viewed and filtered
// by less technical users.
//
// Requires: Groovy NGS Utils (https://github.com/ssadedin/groovy-ngs-utils)
//           ExcelCategory    (https://github.com/ssadedin/excelcatgory)
//
// Author: Simon Sadedin, simon.sadedin@mcri.edu.au
//
/////////////////////////////////////////////////////////////////////////

// Quick and simple way to exit with a message
err = { msg ->
  System.err.println("\nERROR: " + msg + "\n")
  cli.usage()
  System.err.println()
  System.exit(1)
}

// Parse command line args
CliBuilder cli = new CliBuilder(usage: "vcf_to_excel.groovy [options]\n")
cli.with {
  s 'comma separated list of samples to include', args:1
  i 'VCF file to convert to Excel format', args:1
  t 'Name for main spreadsheet tab (should reflect batch, sample group, etc)', args:1
  o 'Name of output file', args:1
}

opts = cli.parse(args)
if(!opts) {
   cli.usage()
   err "Failed to parse command line options"
}

args = opts.arguments()
if(!opts.s) 
    err "Please provide -s option to specify samples to export"
if(!opts.i) 
    err "Please provide -i option to specify VCF file to process"
if(!opts.o) 
    err "Please provide -o option to specify output file name"
if(!opts.t)
    err "Please provide -t option to specify title for spreadsheet"

// These are copied directly from the ##INFO section of an example VCF
// that was processed by VEP. If the flags to VEP are changed, then they
// may need to be updated
VEP_FIELDS = "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|CLIN_SIG|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|DISTANCE|ENSP|HGVSc|HGVSp|PUBMED|SYMBOL|SYMBOL_SOURCE|CANONICAL|SIFT|PolyPhen|AA_MAF|EA_MAF".split("\\|")

samples = opts.s.split(",")

// Parse the VCF. It is assumed that all the samples to be exported are included in the VCF
VCF vcf = VCF.parse(opts.i)

missing_samples = samples.grep { !(it in vcf.samples) }
if(missing_samples)
    err "The following samples were not found in the VCF file provided: ${missing_samples.join(',')}"

new ExcelBuilder().build {

    sheet(opts.t) { 

        // Write out the header columns
        bold {
            row {
                cells("Chr","Pos","Gene","Type","Cons","Id","Qual","DP")
                cells(samples)
                cells("DNA Chg","Protein Chg","1000g","ESP", "SIFT", "PolyPhen")
            }
        }

        for(Variant v in vcf) {

            // Annotations are coming from VEP, which puts them in the CSQ info field
            // A bit of clever groovy code can parse them out
            vep = [VEP_FIELDS,v.info.CSQ.split(",")[0].split("\\|")].transpose().collectEntries()

            row {
                cells(v.chr, v.pos, vep.SYMBOL, v.type, vep.Consequence, v.id, v.qual, v.info.DP)
                
                // Output the genotype for each sample, for this variant
                for(sample in samples) {
                    cell(v.sampleDosage(sample))
                }
                cells(vep.HGVSc, vep.HGVSp, vep.EUR_MAF?:"", vep.EA_MAF?:"", vep.SIFT?:"", vep.PolyPhen?:"")
            }
        }
    }/*.autoFilter("A:"+(char)(65+6+samples.size()))*/.autoSize()
}.save(opts.o)
