//vim: shiftwidth=4:ts=4:expandtab:
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
// The provenance report, produced for each sample in a Cpipe run
//
/////////////////////////////////////////////////////////////////////////////////

// Start by getting the basic information we will need
si = sample_info[sample]

// Pull out the output files that are for this sample specifically
println "Searching for sample $sample in output files"
sampleFiles = outputGraph.findAllOutputsBy { 
    def branchParts = it.branchPath.split("/")
    return (branchParts.contains(sample)|| branchParts[-1]==si.target) 
}

// Pull out the output files that are for this target region (flagship/cohort)
targetFiles = outputGraph.findAllOutputsBy { it.branchPath.split("/").contains(si.target) }

files = [
    rawbam: sampleFiles.grep { it.stageName == "align_bwa" },
    finalbam: sampleFiles.grep { it.stageName.endsWith("recal") },
    vcf: sampleFiles.grep { it.stageName.startsWith("call_variants")  && it.outputFile.name.endsWith("vcf") },
    annovarx: sampleFiles.grep { it.stageName == "vcf_to_excel" && it.outputFile.name.endsWith("annovarx.csv") },
    summary: sampleFiles.grep { it.stageName == "summary_pdf" && it.outputFile.name.endsWith(".pdf") }
].collectEntries { key, fs -> [ key, fs.unique { it.outputFile.absolutePath }[0] ] }

tools = sampleFiles.grep { it.tools }                   // Only files with tools
                   .collect { it.tools.split("\n") }    // A single file has multiple tools -> split/ flatten
                   .flatten()*.trim()*.replaceAll(',$','')*.replaceAll('^,','') // remove leading or trailing commas
                   .unique()                            // Don't list tools multiple times
                   .collect { it.split(":") }           // Split the tool:version for each tool 
                   .collectEntries { it.size()>1 ? [it[0], it[1]] : [it[0], "UNKNOWN"] }  // Create a map of tool => version from above split

// Read the version of the pipeline from the version file in the root directory
versionFile = new File(BASE,"version.txt")
if(versionFile.exists())
    pipelineVersion = versionFile.text.trim()
else
    pipelineVersion = "Unknown"

new PDF().document(outputFile.absolutePath) {

    title "Provenance Report for Study $sample"
    br()

    bold { p "Pipeline" }

    table(cols:2, padding:4) {
        head { cells("Property","Value") }
        cells("Version", pipelineVersion + " / " + new File("revision.txt").text )
        cells("Revision Date", "git log -1 --format=%cd".execute().text)
        cells("Run By", System.properties["user.name"])
        cells("Date", (new Date()).toString())
    }

    br()
    bold { p "Tools" }

    table(cols:2, padding:4, widths: [0.2f,0.4f,0.2f,0.2f]) {
        head {
            cells("Tool","Version")
        }

        tools.each {  tool, version ->
            // BAM files from alignment
            cell tool
            align('center') { cell version }
        }
    }

    br()
    bold { p "Files" }

    table(cols:4,padding:4) {
        head {
            cells("Description","File(s)","Timestamp","File Size")
        }

        fontSize(9) {

            // BAM files from alignment
            cell("Raw Alignment")
            cell(files.rawbam.outputFile.name)
            cell(files.rawbam.timestamp)
            cell(files.rawbam.outputFile.length())

            cell("Final Alignment")
            cell(files.finalbam.outputFile.name)
            cell(files.finalbam.timestamp)
            cell(files.finalbam.outputFile.length())

            cell("Variant Calls")
            cell(files.vcf.outputFile.name)
            cell(files.vcf.timestamp)
            cell(files.vcf.outputFile.length())

            cell("Annotated Variants")
            cell(files.annovarx.outputFile.name)
            cell(files.annovarx.timestamp)
            cell(files.annovarx.outputFile.length())

            cell("Summary PDF")
            cell(files.summary.outputFile.name)
            cell(files.summary.timestamp)
            cell(files.summary.outputFile.length())
        }
    }
}
