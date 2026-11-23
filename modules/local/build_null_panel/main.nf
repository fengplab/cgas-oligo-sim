// STUB ONLY. Generates the composition- and dinucleotide-preserving shuffle
// panel that every occupancy threshold is set against. Not implemented: the
// shuffle method (mononucleotide vs dinucleotide-preserving, and which
// algorithm) is a scientific choice, and the sequence-preference test is only
// interpretable if it is made deliberately.
process BUILD_NULL_PANEL {
    label 'analysis'
    tag "${sequence_id}"

    container params.container_python

    input:
    tuple val(sequence_id), path(fasta)

    output:
    path "null_${sequence_id}.tsv", emit: null_occupancy
    path "versions.yml"           , emit: versions

    script:
    error "BUILD_NULL_PANEL is not implemented. Supply --null_occupancy directly. " +
          "See README section 'Not implemented'."

    stub:
    """
    set -euo pipefail
    printf 'sequence_id\\tshuffle_id\\treplicate\\tposition\\tdwell_frac\\n' > null_${sequence_id}.tsv
    printf '${sequence_id}\\tshuf01\\t1\\t1\\t0.01\\n' >> null_${sequence_id}.tsv
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
