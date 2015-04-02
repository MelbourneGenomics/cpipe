====================== Melbourne Genomics Variant Calling Pipeline =======================

Directory structure under here:

 - pipeline : the source code for running the pipeline
        - scripts : various support scripts

 - batches  
     <batch id> : directories numbered by batch id containing source data
           - data     : source data (FastQ)
           - design   : BED files specifying target regions for the batch
           - analysis : the outputs and working files from execution of the pipeline

 - designs
     
    <design id / name> : BED files and other information associated with a design
                         Each sequencing run will have used a specific design file

 - hg19

   All the HG19 reference files for alignment and downstream processing.
   These should be downloaded from GATK resource bundle for compatibility with GATK.

For information on running and setting up the pipeline, please refer to the
documentation on the Melbourne Genomics Development site. Briefly:

  -  you have to create the file pipeline/config.groovy by customizing 
     pipeline/config.groovy.template

  -  you have to download reference data (and put correct paths to it
     in the above config.groovy), and index it (.fai, .bwt, .dict)

  -  you have to download annovar and its data (see tools/annovar/README)

  -  you have to download VEP data (see tools/vep/README)

  -  Run ./pipeline/scripts/check_install.sh to check everything is in place.
