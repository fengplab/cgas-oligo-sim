// Call occupancy sites for ONE arm.
//
// Arm-scoped on purpose. call_sites.py takes
// OCC_COLS = [sequence_id, replicate, position, dwell_frac] and ignores an arm
// column, so handing it a table that pools arms averages dwell_frac across
// different force fields and structures. Worse, n_rep is
// df["replicate"].nunique() over replicate LABELS, which repeat across arms, so
// a position occupied in any one arm counts as occupied for that replicate and
// replicate_frac inflates toward 1. Running once per arm keeps the statistics
// inside one arm without touching the script's maths.
process CALL_SITES {
    label 'analysis'
    tag  "${arm}/null_p${params.null_percentile}"

    container params.container_python

    // The occupancy parts arrive as a list for THIS arm and are concatenated
    // here. Collapsing them earlier with collectFile would have meant encoding
    // the arm in a filename and parsing it back out, and the sanitised name is
    // not invertible.
    input:
    tuple val(arm), path(occupancy, stageAs: 'occ/*')
    path null_occupancy

    // Output names carry the arm, or several arms publishing to the same
    // directory would overwrite one another.
    output:
    tuple val(arm), path("${prefix}_sites.tsv")         , emit: sites
    tuple val(arm), path("${prefix}_profile.tsv")       , emit: profile
    path "${prefix}_tiers.tsv"                          , emit: tiers
    path "${prefix}_summary.json"                       , emit: summary
    path "versions.yml"                                 , emit: versions

    script:
    prefix = arm.replaceAll(/[^A-Za-z0-9._-]/, '_')
    def rep_min = params.replicate_min == null ? '' : "${params.replicate_min}"
    """
    set -euo pipefail

    # one header, then every part; fail if a part has a different header
    first=\$(ls occ/* | head -1)
    head -1 "\$first" > ${prefix}_occupancy.tsv
    for f in occ/*; do
        if [ "\$(head -1 "\$f")" != "\$(head -1 "\$first")" ]; then
            echo "ERROR: \$f has a different header from \$first" >&2
            exit 1
        fi
        tail -n +2 "\$f" >> ${prefix}_occupancy.tsv
    done

    call_sites.py \\
        --occupancy       ${prefix}_occupancy.tsv \\
        --null-occupancy  ${null_occupancy} \\
        --null-percentile ${params.null_percentile} \\
        --dwell-min       ${params.dwell_min} \\
        --bootstrap-n     ${params.bootstrap_n} \\
        --seed            ${params.seed} \\
        --replicate-min   '${rep_min}' \\
        --out-sites       ${prefix}_sites.tsv \\
        --out-profile     ${prefix}_profile.tsv \\
        --out-tiers       ${prefix}_tiers.tsv \\
        --out-summary     ${prefix}_summary.json

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
	    numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
	    pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
	VERSIONS
    """

    stub:
    prefix = arm.replaceAll(/[^A-Za-z0-9._-]/, '_')
    """
    set -euo pipefail
    printf 'sequence_id\\tsite_id\\tpeak_position\\tbed_start\\tbed_end\\trun_start\\trun_end\\twidth\\toccupancy\\toccupancy_ci_lo\\toccupancy_ci_hi\\treplicate_frac\\tn_occupying_replicates\\tnull_threshold\\tnull_percentile\\tn_shuffles\\treplicate_gate\\trendered\\n' > ${prefix}_sites.tsv
    printf 'SEQ_STUB\\tSEQ_STUB_s1\\t1\\t0\\t1\\t1\\t1\\t1\\t0.5\\t0.4\\t0.6\\t0.9\\t9\\t0.02\\t${params.null_percentile}\\t10\\tunset\\tTrue\\n' >> ${prefix}_sites.tsv
    printf 'sequence_id\\tposition\\toccupancy\\tn_replicates_seen\\tn_occupying\\treplicate_frac\\tnull_threshold\\tabove_threshold\\n' > ${prefix}_profile.tsv
    printf 'sequence_id\\treplicate_frac_tier\\tobserved_positions\\tnull_positions_per_shuffle\\n' > ${prefix}_tiers.tsv
    echo '{"stub": true}' > ${prefix}_summary.json
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
