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


variant_annotation = segment {
   annotate_vep + 
   index_vcf +
   annovar_table +
   [ 
       add_to_database, 
       augment_condel + 
       annotate_significance
   ]
}

family_vcf = segment { 
    merge_target_vcfs + 
    annotate_snpeff + 
    index_vcf.using(sort_vcf:false) + 
    vcf_to_family_excel 
}


@filter("vep")
annotate_vep = {
    doc "Annotate variants using VEP to add Ensemble annotations"
    output.dir="variants"

    // Note: if the input VCF file is empty, VEP will not create an output.
    // To avoid this causing the pipeline to fail, we first copy only the
    // headers from the input file to the output file so that if VEP does not
    // overwrite the file, we end up with an empty file
    exec """
        grep '^#' $input.vcf > $output.vcf 

        PERL5LIB="$CONDEL:\$PERL5LIB"
        perl $VEP/variant_effect_predictor.pl --cache --dir $VEP/../vep_cache 
            -i $input.vcf 
            --vcf -o $output.vcf 
            -species homo_sapiens 
            --canonical --per_gene --protein 
            --sift=b --polyphen=b
            --symbol hgnc --force_overwrite --hgvs  --maf_1kg --maf_esp --pubmed
            --plugin Condel,$CONDEL/config,s
            --offline
            --verbose
    """, "vep"
}

index_vcf = {
    output.dir="variants"
    var sort_vcf : true
    if(sort_vcf) {
        transform("vcf") to("sort.vcf") {
            exec """

                $IGVTOOLS/igvtools sort $input.vcf $output.vcf 

                if [ ! -e $output.vcf ];
                then
                    grep '^#' $input.vcf > $output.vcf ;
                fi

                $IGVTOOLS/igvtools index $output.vcf
            """
        }
    }
    else {
        transform("vcf") to("vcf.idx") {
            exec """
                $IGVTOOLS/igvtools index $input.vcf
            """
        }
    }
}


annovar_table = {

    output.dir="variants"

    transform("vcf","vcf","vcf") to("av", "hg19_multianno.csv","refGene.exonic_variant_function") {
        exec """
            $ANNOVAR/convert2annovar.pl $input.vcf -format vcf4 > $output.av

            $ANNOVAR/table_annovar.pl $output.av $ANNOVAR_DB/  -buildver hg19 
            -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,exac03,snp138,ljb26_all
            -operation g,r,r,f,f,f,f,f 
            -nastring . 
            --otherinfo   
            --csvout
            --outfile $output.csv.prefix.prefix
            --argument '-exonicsplicing -splicing $INTERVAL_PADDING_CALL',,,,,,,

            sed -i '/^Chr,/ s/\\.refGene//g' $output.csv
        """
    }
}

add_to_database = {

    doc "Add discovered variants to a database to enable annotation of number of observations of the variant"

    requires UPDATE_VARIANT_DB : "SQLite database to store variants in"

    output.dir="variants"

    uses(variantdb:1) {
        exec """

            echo "====> Adding variants for flaship $target_name to database"

            JAVA_OPTS="-Xmx24g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/vcf_to_db.groovy 
                   -v $input.recal.vcf 
                   -a $input.csv 
                   -db $UPDATE_VARIANT_DB 
                   -cohort $target_name
                   -idmask '$SAMPLE_ID_MASK'
                   -b "$batch"

            echo "<==== Finished adding variants for flaship $target_name to database"

            echo "Variants from $input.recal.vcf were added to database $VARIANT_DB on ${new Date()}" > $output.txt
        """, "add_to_database"
    }
}

augment_condel = {

    doc "Extract Condel scores from VEP annotated VCF files and add them to Annovar formatted CSV output"

    output.dir="variants"
    from("*.hg19_multianno*.csv") filter("con") {
        exec """
            JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY 
                -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar 
                $SCRIPTS/merge_condel.groovy
                    -i $input.vcf
                    -a $input.csv
                    -o $output.csv
        """
    }
}

@filter("sig")
annotate_significance = {
    doc "Add clinical significance category annotations as defined by Melbourne Genomics"
        var MAF_THRESHOLD_RARE : 0.01,
            MAF_THRESHOLD_VERY_RARE : 0.0005,
            CONDEL_THRESHOLD : 0.7

        output.dir="variants"
        from("con.csv") {
            exec """
                python $SCRIPTS/annotate_significance.py 
                --annovar $input.csv
                --rare $MAF_THRESHOLD_RARE
                --very_rare $MAF_THRESHOLD_VERY_RARE
                --condel $CONDEL_THRESHOLD
                --synonymous $COMBINED_SYNONYMOUS
                > $output.csv
                """
        }
}

merge_target_vcfs = {

    doc "Merges only the VCFs belonging to the current analysis profile"

    output.dir="variants"

    var enable_family_excel : false
    if(!enable_family_excel)
        succeed "Family VCF output not configured for $target_name"

    // Find the sample names that are in the current analysis profile
    def target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample

    // Work with raw filtered VCFs, not ones that were already annotated
    def variant_files = target_samples*.plus(".*.filter.vcf")

    from(variant_files) produce(target_name + ".merge.vcf") {
        msg "Merging vcf files: " + inputs.vcf
        exec """
                $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
                -T CombineVariants
                -R $REF
                ${inputs.vcf.withFlag("--variant")}
                --out $output.vcf
             """
    }
}

@filter("snpeff")
annotate_snpeff = {
    output.dir="variants"

    var enable_snpeff:false,
        SNPEFF : false

    if(!enable_snpeff)
        succeed "Snpeff support not enabled"

    if(!SNPEFF)
        fail "Please define the SNPEFF variable to point to the location of your SNPEFF installation"

    exec """
            $JAVA -Xmx2g -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -treatAllAsProteinCoding false -a 2 hg19 $input.vcf  > $output.vcf
    ""","snpeff"
}

vcf_to_family_excel = {

    var with_sex_karyotype : false,
        unique_variants : false,
        min_priority : 0

    doc """
         Convert variants annotated with SnpEff and Annovar to an excel based format
         designed for diagnostics from family based sequencing.
        """

    requires ANNOTATION_VARIANT_DB : "File name of SQLite variant database for storing variants"

    output.dir = "results"

    def UNIQUE = unique_variants.toBoolean() ? " -unique " : ""

    def target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample
    
    def annovar_files = target_samples*.plus(".*.hg19_multianno.con.sig.csv") 

    println "Annovar files for target $target_name are " + annovar_files
    from(annovar_files) {
        produce(target_name + ".family.xlsx") {
            exec """
                JAVA_OPTS="-Xmx12g -Djava.awt.headless=true" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/vcf_to_excel.family.groovy 
                    -targets $target_bed_file ${with_sex_karyotype ? "-sex" : ""}
                    -p "" 
                    -db $ANNOTATION_VARIANT_DB
                    -minpri $min_priority
                    -gc $target_gene_file
                    -ped $input.ped ${inputs.bam.withFlag("-bam")}
                    -o $output.xlsx
                    -p $run_id
                    $UNIQUE $input.vcf $inputs.csv 
            """, "vcf_to_family_excel"
        }
    }
}

