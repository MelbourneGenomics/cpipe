#!/usr/bin/env groovy
CliBuilder cli = new CliBuilder()
cli.with {
    batch "batch to which samples belong", args:1, required: true
    disease "disease cohort to which samples belong", args:1, required: true
}

opts = cli.parse(args)
if(!opts)
        System.exit(0)

samples = SampleInfo.fromFiles(opts.arguments() as List)
samples.each { it.value.batch = opts.batch; it.value.target = opts.disease }

println(
    samples*.value*.toTsv().join("\n")
)
