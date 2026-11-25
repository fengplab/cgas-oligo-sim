// Pool bound frames per called site into the pose ensembles CLUSTER_POSES reads.
//
// Runs once per (arm, sequence_id), after that arm's sites are called, and sees
// every replicate of that pair at once - the ensemble is the pooled frames, and
// params.min_frames is checked against the pool.
//
// Replicate identity is written into a manifest from the channel's own values,
// never taken from the order in which the trajectories were staged.
process EXTRACT_POSES {
    label 'analysis'
    tag  "${sequence_id}/${arm}"

    container params.container_python

    input:
    tuple val(arm), val(sequence_id), val(replicates), path(trajectories),
          path(bead_map), path(sites)

    output:
    tuple val(arm), val(sequence_id), path("poses/*.npz"), emit: poses, optional: true
    tuple val(arm), path("${prefix}_poses.tsv")          , emit: manifest
    path "${prefix}_poses_summary.json"                  , emit: summary
    path "versions.yml"                                  , emit: versions

    script:
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}"
    // Zip the replicate labels with the staged filenames positionally HERE,
    // where both lists come from the same channel item and are guaranteed to
    // correspond. This is the only place the pairing exists, which is why it is
    // written to a manifest the script reads by column.
    def reps  = replicates instanceof List ? replicates : [replicates]
    def trajs = trajectories instanceof List ? trajectories : [trajectories]
    if (reps.size() != trajs.size()) {
        error "EXTRACT_POSES ${sequence_id}/${arm}: ${reps.size()} replicate " +
              "labels for ${trajs.size()} trajectories; the pairing is unsafe."
    }
    def manifest_rows = [reps, trajs].transpose()
        .collect { r, t -> "${r}\\t${t}" }.join('\\n')
    """
    set -euo pipefail
    mkdir -p poses

    printf 'replicate\\tpath\\n${manifest_rows}\\n' > ${prefix}_manifest.tsv

    extract_poses.py \\
        --manifest ${prefix}_manifest.tsv \\
        --bead-map ${bead_map} \\
        --sites ${sites} \\
        --sequence-id ${sequence_id} \\
        --arm '${arm}' \\
        --contact-cutoff ${params.contact_cutoff} \\
        --footprint-min ${params.footprint_min} \\
        --equil-frac ${params.equil_frac} \\
        --site-overlap ${params.pose_site_overlap} \\
        --pose-beads ${params.pose_beads} \\
        --min-frames ${params.min_frames} \\
        --out-dir poses \\
        --out-manifest ${prefix}_poses.tsv \\
        --out-summary ${prefix}_poses_summary.json

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 --version | awk '{print \$2}')
	    pose_beads: ${params.pose_beads}
	    site_overlap: ${params.pose_site_overlap}
	    contact_cutoff: ${params.contact_cutoff}
	VERSIONS
    """

    stub:
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}"
    """
    set -euo pipefail
    mkdir -p poses
    # Emit one real pose so a stub run exercises CLUSTER_POSES and BUILD_INDEX.
    # A stub that produces an empty manifest leaves both unstarted and the run
    # still exits 0, which hides the rest of the chain.
    python3 - <<'PYSTUB'
import numpy as np
rng = np.random.default_rng(0)
np.savez_compressed("poses/${sequence_id}_stub_poses.npz",
                    coords=rng.normal(size=(60, 8, 3)).astype("float32"),
                    bead_labels=np.array([f"pro:{i}" for i in range(8)], dtype="U32"),
                    replicate=np.ones(60, dtype=int))
PYSTUB
    printf 'sequence_id\\tsite_id\\tnpz_path\\n' > ${prefix}_poses.tsv
    printf '${sequence_id}\\t${sequence_id}_stub\\t%s/poses/${sequence_id}_stub_poses.npz\\n' "\$PWD" \\
        >> ${prefix}_poses.tsv
    echo '{"stub": true}' > ${prefix}_poses_summary.json
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
