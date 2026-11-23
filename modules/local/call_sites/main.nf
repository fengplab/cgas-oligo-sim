process CALL_SITES {
    label 'analysis'
    tag  "null_p${params.null_percentile}"

    container params.container_python

    input:
    path occupancy
    path null_occupancy

    output:
    path "sites.tsv"              , emit: sites
    path "occupancy_profile.tsv"  , emit: profile
    path "tiers.tsv"              , emit: tiers
    path "call_sites_summary.json", emit: summary
    path "versions.yml"           , emit: versions

    script:
    def rep_min = params.replicate_min == null ? '' : "${params.replicate_min}"
    """
    set -euo pipefail

    call_sites.py \\
        --occupancy       ${occupancy} \\
        --null-occupancy  ${null_occupancy} \\
        --null-percentile ${params.null_percentile} \\
        --dwell-min       ${params.dwell_min} \\
        --bootstrap-n     ${params.bootstrap_n} \\
        --seed            ${params.seed} \\
        --replicate-min   '${rep_min}' \\
        --out-sites       sites.tsv \\
        --out-profile     occupancy_profile.tsv \\
        --out-tiers       tiers.tsv \\
        --out-summary     call_sites_summary.json

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
	    numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
	    pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    printf 'sequence_id\\tsite_id\\tpeak_position\\tbed_start\\tbed_end\\trun_start\\trun_end\\twidth\\toccupancy\\toccupancy_ci_lo\\toccupancy_ci_hi\\treplicate_frac\\tn_occupying_replicates\\tnull_threshold\\tnull_percentile\\tn_shuffles\\treplicate_gate\\trendered\\n' > sites.tsv
    printf 'SEQ_STUB\\tSEQ_STUB_s1\\t1\\t0\\t1\\t1\\t1\\t1\\t0.5\\t0.4\\t0.6\\t0.9\\t9\\t0.02\\t${params.null_percentile}\\t10\\tunset\\tTrue\\n' >> sites.tsv
    printf 'sequence_id\\tposition\\toccupancy\\tn_replicates_seen\\tn_occupying\\treplicate_frac\\tnull_threshold\\tabove_threshold\\n' > occupancy_profile.tsv
    printf 'sequence_id\\treplicate_frac_tier\\tobserved_positions\\tnull_positions_per_shuffle\\n' > tiers.tsv
    echo '{"stub": true}' > call_sites_summary.json
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
