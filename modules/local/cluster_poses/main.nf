process CLUSTER_POSES {
    label 'cluster'
    tag  "${site_id}"

    container params.container_python

    input:
    tuple val(sequence_id), val(site_id), path(npz)
    path contacts

    output:
    tuple val(site_id), path("${site_id}.pdb")          , emit: pdb
    tuple val(site_id), path("${site_id}_cluster.json") , emit: summary
    path "${site_id}_clusters.tsv"                      , emit: clusters
    path "versions.yml"                                 , emit: versions

    script:
    if (params.rmsd_cutoff == null) {
        error "params.rmsd_cutoff is not set. A cutoff below the data's noise floor " +
              "silently yields one cluster per frame. Run cluster_poses.py --report-only " +
              "on a representative site first and choose from the RMSD sweep."
    }
    """
    set -euo pipefail

    cluster_poses.py \\
        --npz                    ${npz} \\
        --site-id                ${site_id} \\
        --sequence-id            ${sequence_id} \\
        --contacts               ${contacts} \\
        --rmsd-cutoff            ${params.rmsd_cutoff} \\
        --min-frames             ${params.min_frames} \\
        --min-top-population-frac ${params.min_top_pop_frac} \\
        --out-pdb                ${site_id}.pdb \\
        --out-clusters           ${site_id}_clusters.tsv \\
        --out-summary            ${site_id}_cluster.json

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
	    numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    printf 'TITLE     ${site_id} stub\\nEND\\n' > ${site_id}.pdb
    echo '{"site_id": "${site_id}", "stub": true, "top_cluster_population_frac": 1.0}' > ${site_id}_cluster.json
    printf 'sequence_id\\tsite_id\\tcluster_rank\\tcentre_frame\\tpopulation\\tpopulation_frac\\tn_replicates_in_cluster\\tmean_rmsd_to_centre\\n' > ${site_id}_clusters.tsv
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
