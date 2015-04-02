/**
 * Filter to find variants that differ between a number of VCFs
 * Writes an output VCF with a "DS" tag specifying which input files
 * the variant was found in.
 *
 * NOTE: Ignores sample information. Assumes each VCF contains 1 sample.
 * Requires Groovy-NGS-Utils (https://github.com/ssadedin/groovy-ngs-utils)
 */
import org.apache.commons.cli.Option

Cli cli = new Cli()
cli.with { 
        vcf 'VCFs to compare (supply multiple times)', args: Option.UNLIMITED_VALUES, required:true
}

opts = cli.parse(args)
if(!opts) 
        System.exit(1)

// Read all the VCFs
vcfs = opts.vcfs.collect { VCF.parse(it) }

// Create a super-VCF containing all the variants from all the VCFs
combined = new VCF(headerLines: vcfs[0].headerLines)
vcfs.each { vcf ->
  for(Variant v in vcf) { 
        combined.add(v)
  }
}

// Now write out the VCF containing only variants that differ between the three
// Note: ignore heterozygosity here.
VCF output = new VCF(headerLines:vcfs[0].headerLines)
for(Variant v in combined) {
    def foundIn = vcfs.collect { v in it }
    if(!foundIn.every()) {
        v.update {
          v.info.DS=foundIn.collect { it ? 1 : 0 }.join(",")     
        }
        output.add(v)
    } 
}

output.print()

