import groovy.util.logging.Log

/**
 * 
 * Creates a directory structure with all the necessary files ready to import into Seqr
 * 
 * @author Simon Sadedin
 */
@Log
class CreateSeqrProject {
    
    private Set allParents
    
    String batchId
    
    Map<String,SampleInfo> sampleInfo
    
    File seqrDir = new File("seqr")
        
    CreateSeqrProject(OptionAccessor opts) {
        this.batchId = opts.batch ?: new File(".").absoluteFile.parentFile.parentFile.name
        this.sampleInfo = SampleInfo.parse_mg_sample_info(opts.sample_info ?: 'samples.txt')
        
        log.info "Creating Seqr project for batch $batchId based on $opts.sample_info"
    }
    
    /**
     * Create all parts of the Seqr project (main method)
     */
    void create() {
        seqrDir.mkdir()
        
        createSampleListFile()
        
        createPEDFile()
        
        createYAMLFile()
    }
    
    /**
     * Discover all the sample ids and write them to a file for inclusion in project.yaml
     */
    void createSampleListFile() {
        new File(seqrDir,"all_samples.txt").text = sampleInfo*.key.join("\n") +"\n"
    }
    
    /**
     * Create the final YAML file, ready for import
     */
    void createYAMLFile() {
        // Create the YAML file
        new File(seqrDir,"project.yaml").text = 
        """
            project_id: '${batchId}'

            project_name: 'Cpipe Batch ${batchId}'

            sample_id_list: 'all_samples.txt'
            
            ped_files:
              - 'samples.ped'
            
            vcf_files:
            ${getVCFFiles().collect {"  - " + it}*.plus('\n').join('')}

        """.stripIndent()
    }
    
    /**
     * Discover all the appropriate VCF files for inclusion in upload.
     * <p>
     * This is complicated by the fact htat trios and singletons are handled differently -
     * when a trio VCF exists, this is used as it should contain the most accurate 
     * calls for each of the samples concerned. When no trio VCF exists for a sample,
     * the singleton VCF is identified.
     * 
     * @return  List of absolute paths to VCF files
     */
    List<String> getVCFFiles() {
        
        List<File> vcfFiles = new File("variants").listFiles()
        
        // We expect to find 1 VEP annotated VCF file in the variants directory for each proband
        // When the proband has parents, it ends with "trio.genotype.soi.vep.vcf"
        // When the proband has no parents, it just ends with "individual.genotype.soi.vep.vcf"
        sampleInfo*.key.grep { !(it in allParents) }.collect { sampleId ->
            SampleInfo si = sampleInfo[sampleId]
            String prefix = null
            if(si.pedigree) {
                prefix = "trio"
            }
            else {
                prefix = "individual"
            }
            
            File vcfFile = vcfFiles.find {  File f ->
                f.name.startsWith(sampleId) && f.name.endsWith("${prefix}.genotype.soi.vep.ens.vcf")
            }
            
            if(vcfFile == null) {
                throw new RuntimeException("No VCF could be found for sample $sampleId")
            }
            
            vcfFile
        }*.absolutePath
    }
    
    /**
     * Create a PED file containing all the samples included in the batch
     * <p>
     * Note that although a PED file is created by Cpipe, that PED file is not
     * quite usable. First, it only exists for trio samples (I think), and second
     * it doesn't include all samples combined in one file. This method creates
     * a single PED file for all samples that is consumable by Seqr.
     */
    void createPEDFile() {
        
        // Parents are stored in the pedigree field
        Map<String,List<String>> parents = sampleInfo.collectEntries { e ->
            [e, e.value.pedigree.replaceAll('familyID=','').tokenize(',')]
        }
        
        allParents = parents*.value.flatten() as Set
        
        List<Pedigree> peds = []
        
        // Generate a pedigree for each proband (for now, single proband families are all we support)
        // At this point, I am assuming that a sample is a proband if 
        // - there is info in the Pedigree field
        // or
        // - no other sample mentions this one has being a parent
        //
        for(String sample in sampleInfo.keySet()) {
            SampleInfo info = sampleInfo[sample]
            
            if(sample in allParents) { // parent not a proband, skip
                log.info "Sample $sample is a parent: assuming unaffected"
                continue
            }
            
            log.info "Inferring pedigree for $sample"
            peds << getSamplePedigree(sample, parents[sample])
        }
        
        new File(seqrDir,"samples.ped").withWriter { w ->
            peds.each { it.toPed(w); }
        }
    }
    
    /**
     * Construct a pedigree object for the given sample and parent set
     * 
     * @param sampleId
     * @param parents
     * @return
     */
    Pedigree getSamplePedigree(String sampleId, List<String> parents) {
        Pedigree ped = new Pedigree()
        
        SampleInfo si = sampleInfo[sampleId]
        ped.id = sampleId
        
        Subject proband = new Subject()
        proband.id = sampleId
        proband.sex = si.sex
        ped.individuals << proband
        ped.phenoTypes << 2
        
        RelationshipType probandRel = proband.sex == Sex.MALE ? RelationshipType.SON : RelationshipType.DAUGHTER 
        
        // Find mum, if there is one
        List<Subject> mum = parents.grep { sampleInfo[it].sex == Sex.FEMALE }.collect { SampleInfo pi ->
            Subject s = new Subject()
            s.id = pi.sample
            s.sex = Sex.FEMALE
            s.relationships << new Relationship(type:RelationshipType.MOTHER, to:proband.id)
            proband.relationships << new Relationship(type:probandRel, to:proband.id)
        }
        
        if(mum) {
            ped.individuals << mum[0]
            ped.phenoTypes << 1
        }
        
        // Find dad, if there is one
        List<Subject> dad = parents.grep { sampleInfo[it].sex == Sex.MALE }.collect { SampleInfo pi ->
            Subject s = new Subject()
            s.id = pi.sample
            s.sex = Sex.MALE
            s.relationships << new Relationship(type:RelationshipType.FATHER, to:proband.id)
            proband.relationships << new Relationship(type:probandRel, to:proband.id)
        }
        
        if(dad) {
            ped.individuals << dad[0]
            ped.phenoTypes << 1
        }
        
        return ped
    }

    static void main(String [] args) {
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli()
        cli.with { 
            batch "Batch id", args:1, required: false
            sample_info "Sample Info File", args:1, required:true
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) 
            System.exit(1)
        new CreateSeqrProject(opts).create()
        
        log.info "Done."
    }

}
