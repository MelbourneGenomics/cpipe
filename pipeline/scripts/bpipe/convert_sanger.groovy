convert_sanger = {
        doc """
                Convert old-format illumina reads to Sanger (Phred+33) format.
                Requires manually patched version of MAQ.
        """
        var MAQ : "maq"

        // filter("phred") {
        transform("fastq.gz") to("sanger.fastq.gz") {
            exec """
                    gunzip -c  $input.gz | $MAQ ill2sanger  - - | gzip -c > $output.gz
            """
        }
}

run {
   "%.gz" * [ convert_sanger ]
}
