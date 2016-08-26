// vim: ts=4:expandtab:sw=4:cindent

load 'config.groovy'

init_size = {
    branch.variant_size = branch.name.toInteger()
    branch.size_min = size_ranges[variant_size].from
    branch.size_max = size_ranges[variant_size].to
}

@filter("filt")
filter_vcf_to_target = {

    output.dir="eval"
    msg "Evaluating $variant_types of size $size_min - $size_max"

    produce(file(input.vcf).name.replaceAll('\\.vcf$',".${size_min}_${size_max}.vcf")) {
        exec """
            $BEDTOOLS/bin/bedtools intersect -header -u -a $input.vcf -b $EXOME_TARGET | awk '{ if(\$1 != "chrY" && \$1 != "chrM" ) print \$0 }' |
                JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY 
                    -cp $GROOVY_NGS/groovy-ngs-utils.jar
                    -e 'VCF.filter { Variant v -> (v.filter in ["PASS","."]) && !v.isSV() && (v.size() >= $size_min) && (v.size() <= $size_max) }'
                    > $output.vcf
        """
    }
}

@filter("sample") 
extract_sample = {
    exec """
         cat $input.vcf | groovy -e 'VCF.filter(samples:["$sample"])' > $output.vcf
         """
}

@filter("rename")
rename_vcf_sample = {
    output.dir="eval"
    exec """
         cat $input.vcf | awk 'BEGIN { OFS="\\t" } { if(\$1 == "#CHROM") \$10="$sample"; print \$0 }' > $output.vcf
         """
}

rtg_bgzip = {

    var set_baseline : false

    output.dir="eval"
    transform(".vcf") to(".vcf.gz") {
        exec """
         $RTG bgzip -c $input.vcf > $output.gz

         $RTG index -f vcf $output.gz
        """

        if(set_baseline)
            branch.parent.baseline_zip = output.gz
    }
}

rtg_eval = {

    requires sample : "The sample to compare"

    // Name the output directory that RTG writes to based on the
    // original VCF
    output.dir = "eval/roc/"+file(input.gz).name.replaceAll('\\.merge.*\\.vcf\\.gz','').replaceAll('\\.','_') + "_${size_min}_${size_max}"

    produce("heterozygous_roc.tsv") {
        // Note: Bpipe creates the "roc" directory but then RTG refuses to use it
        // because it insists on creating the directory itself. So remove it first.
        exec """
            rm -rf $output.dir;
            $RTG vcfeval --no-gzip --sample $sample -b $baseline_zip -c $input1.gz -t ${REF}.sdf -o $output.dir
        """
    }
}

rtg_roc_plot = {
    exec """
         $RTG rocplot --png $output.png -t "ROC Curve for Sensitivity for Different Median Coverages" $inputs.tsv
         """
}

requires baseline : "VCF of the baseline sample to compare with",
         sample : "ID of baseline sample to compare against. Must be a valid sample id in the baseline VCF.",
         EXOME_TARGET : "The BED file describing target regions for analysis"

testsamples = args.collectEntries { [it,it] }

size_ranges = [0..0, 1..1, 2..3, 4..6]  

sizes = (0..<size_ranges.size())*.toString()

run {
   sizes * [
     init_size +

     // Step 1: filter the gold standard calls to our target region
     [ baseline : baseline ] * [ filter_vcf_to_target + rtg_bgzip.using(set_baseline:true) ]  +

     // Step 2: filter all the test samples also to the target region, 
     //         also rename the sample within the VCF to have the same
     //         sample id as the gold standard, and then compare
     //         each sample to baseline using RTG eval
     testsamples * [ 
         filter_vcf_to_target + rename_vcf_sample + rtg_bgzip + rtg_eval 
     ]
 ] 
}


