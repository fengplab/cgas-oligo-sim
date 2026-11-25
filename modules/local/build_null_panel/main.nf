// STUB ONLY: not yet wired into the workflow, and not yet implemented.
//
// DECIDED 2026-11-25: MONONUCLEOTIDE shuffle. It preserves base composition and
// randomises arrangement. The rationale is that composition here is fixed by the
// prior ChIP data, so this null supports the specificity of that input and of
// the motifs discovered in that workflow - it asks whether the arrangement
// matters given the composition, not whether composition matters.
//
// Known limitation to state in the methods: 3SPN.2C encodes sequence-dependent
// DNA geometry and flexibility at the base-step level, so a mononucleotide
// shuffle also changes dinucleotide step statistics. An occupancy difference
// between observed and null therefore cannot, on its own, separate motif
// arrangement from altered bendability. A dinucleotide-preserving shuffle would
// isolate arrangement but is conservative, absorbing real signal if preference
// is itself a dinucleotide-level shape effect. Reporting both nulls would
// bound the answer from each side if the question becomes contested.
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
