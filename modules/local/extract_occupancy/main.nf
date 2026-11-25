// Turn one CG trajectory into an occupancy profile along the oligo.
//
// Sits between GENESIS_CG_RUN and CALL_SITES. The bead map comes from
// GENESIS_SETUP, joined on sequence_id rather than inferred here: bead identity
// must never be guessed from position in the file, and a bead-count mismatch
// between map and trajectory is a hard failure in extract_occupancy.py.
//
// The arm label is carried through and written into every row. Arms are
// different force fields or different structures, so pooling across them would
// average over distinct physics - extract_occupancy.py says so itself.
process EXTRACT_OCCUPANCY {
    label 'analysis'
    tag  "${sequence_id}/${arm}/rep${replicate}"

    // Python, not GENESIS - so this follows container_python like the other
    // analysis modules. Keying it on genesis_bin would run python natively just
    // because atdyn happens to be a local build.
    container params.container_python

    input:
    tuple val(sequence_id), val(lineage_id), val(arm), val(replicate),
          path(traj), path(bead_map)

    output:
    tuple val(sequence_id), val(arm), path("${prefix}_occupancy.tsv")   , emit: occupancy
    tuple val(sequence_id), val(arm), path("${prefix}_contacts.tsv")    , emit: contacts
    tuple val(sequence_id), val(arm), path("${prefix}_orientation.tsv") , emit: orientation
    path "${prefix}_occupancy_summary.json"                             , emit: summary
    path "versions.yml"                                                 , emit: versions

    script:
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}_rep${replicate}"
    """
    set -euo pipefail

    extract_occupancy.py \\
        --trajectory ${traj} \\
        --bead-map ${bead_map} \\
        --sequence-id ${sequence_id} \\
        --replicate ${replicate} \\
        --arm '${arm}' \\
        --contact-cutoff ${params.contact_cutoff} \\
        --footprint-min ${params.footprint_min} \\
        --equil-frac ${params.equil_frac} \\
        --min-frames ${params.min_frames_occ} \\
        --out-occupancy ${prefix}_occupancy.tsv \\
        --out-contacts ${prefix}_contacts.tsv \\
        --out-orientation ${prefix}_orientation.tsv \\
        --out-summary ${prefix}_occupancy_summary.json

    # A replicate in which cGAS never bound is a VALID result and must survive
    # to the null comparison, so emptiness is not failed here. What is not valid
    # is a missing or header-only file, which would mean the extractor wrote
    # nothing at all.
    for f in ${prefix}_occupancy.tsv ${prefix}_contacts.tsv ${prefix}_orientation.tsv; do
        if [ ! -s "\$f" ]; then
            echo "ERROR: \$f was not written" >&2
            exit 1
        fi
    done

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 --version | awk '{print \$2}')
	    contact_cutoff: ${params.contact_cutoff}
	    footprint_min: ${params.footprint_min}
	    equil_frac: ${params.equil_frac}
	VERSIONS
    """

    stub:
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}_rep${replicate}"
    """
    set -euo pipefail
    printf 'sequence_id\\tarm\\treplicate\\tposition\\tdwell_frac\\n' > ${prefix}_occupancy.tsv
    printf '${sequence_id}\\t${arm}\\t${replicate}\\t1\\t0.0\\n' >> ${prefix}_occupancy.tsv
    printf 'stub\\n' > ${prefix}_contacts.tsv
    printf 'stub\\n' > ${prefix}_orientation.tsv
    echo '{"stub": true}' > ${prefix}_occupancy_summary.json
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
