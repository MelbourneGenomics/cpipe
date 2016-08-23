#!/usr/bin/env bash
ARGS=$(getopt -o e:c:g: --long "ensembl-release:,cache-type:l,usage,help,genome-build" -n $(basename $0) -- "$@")
eval set -- "$ARGS"

case "$1" in
    batch)
        #  e.g. docker run cpipe batch add_batch --batch batch_identifier --profile profile_name
        shift 1
        python pipeline/scripts/manage_batch.py "$@" < /dev/stdin
    ;;
    genelist)
        # e.g. docker run cpipe genelist show_bed --profile profile_name
        shift 1
        python pipeline/scripts/manage_genelists.py "$@" < /dev/stdin
    ;;
    metadata)
        shift 1
        case "$1" in
            check)
                #e.g. docker run cpipe metadata check < ./batches/batch_identifier/samples.txt
                python pipeline/scripts/check_metadata.py "$@" < /dev/stdin
            ;;
            update)
                #e.g. docker run cpipe metadata update --sample sample_name --name prioritised_genes --value “4:ABC1,ABC2” --target ./batches/batch_identifier/samples.txt
                python pipeline/scripts/update_metadata.py "$@" < /dev/stdin
            ;;
        esac
    # e.g. docker run cpipe genelist show_bed --profile profile_name
    ;;
    pipeline)
        cd batches/batch/analysis
        ../../../bpipe run ../../../pipeline/pipeline.groovy ../samples.txt "$@" < /dev/stdin
    ;;
esac