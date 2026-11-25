// Build the flat job matrix for the whole campaign.
//
// One row per trajectory, so the workflow can fan out over a single channel
// instead of nesting loops: arm x oligo x panel x replicate. Each row carries
// its own seed, derived from (global seed, arm, oligo, panel, replicate), which
// is what keeps replicates independent across arms and oligos.
//
// PWMcos arms appear only when motifs are supplied, so a no-motif run produces
// the baseline and structure arms and nothing else.
process ENUMERATE_ARMS {
    label 'analysis'

    container params.container_python

    input:
    path panel
    path structures
    path motifs

    output:
    path "jobs.tsv"              , emit: jobs
    path "arms.tsv"              , emit: arms
    path "arms_summary.json"     , emit: summary
    path "versions.yml"          , emit: versions

    script:
    // An empty motif table is the no-motif case, not an error: enumerate_arms.py
    // takes an absent --motifs to mean "no PWMcos arms".
    // Must match assets/NO_FILE exactly. An earlier version tested a different
    // name, so the placeholder was passed in as a real motif table and
    // enumerate_arms.py tried to read its comment line as a header.
    def motif_opt = motifs.name != 'NO_FILE' ? "--motifs ${motifs}" : ''
    def mb_opt    = params.multibasin ? '--multibasin' : '--no-multibasin'
    """
    set -euo pipefail

    enumerate_arms.py \\
        --panel ${panel} \\
        --structures ${structures} \\
        ${motif_opt} \\
        --partners '${params.partners ?: ''}' \\
        --replicates ${params.cg_replicates} \\
        --shuffles ${params.null_shuffles} \\
        --null-replicates ${params.null_replicates} \\
        ${mb_opt} \\
        --seed ${params.seed} \\
        --out-jobs jobs.tsv \\
        --out-arms arms.tsv \\
        --out-summary arms_summary.json

    # An empty matrix means the panel or structure sheet was empty, which would
    # otherwise produce a pipeline that runs nothing and reports success.
    if [ "\$(wc -l < jobs.tsv)" -lt 2 ]; then
        echo "ERROR: jobs.tsv has no job rows" >&2
        exit 1
    fi

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 --version | awk '{print \$2}')
	    replicates: ${params.cg_replicates}
	    seed: ${params.seed}
	VERSIONS
    """

    stub:
    """
    set -euo pipefail

    # Derive ids from the real input sheets rather than hardcoding them, so a
    # stub run exercises the downstream joins against the SAME identifiers a
    # real run would use. A stub that invents ids proves only that the
    # processes execute, not that the channels line up.
    printf 'job_id\\tarm\\tmode\\tpartner\\tmotif\\tstructures\\tbasins\\tsequence_id\\tpanel\\tshuffle\\treplicate\\tseed\\n' > jobs.tsv

    struct0=\$(awk -F'\\t' 'NR==2 {print \$1}' ${structures})
    if [ -z "\$struct0" ]; then
        echo "ERROR: ${structures} has no data row" >&2; exit 1
    fi

    awk -F'\\t' -v s="\$struct0" -v reps=${params.cg_replicates} '
        NR > 1 && \$1 != "" {
            for (r = 1; r <= reps; r++)
                printf "baseline|%s|obs|r%d\\tbaseline\\tfree\\t\\t\\t%s\\t1\\t%s\\tobserved\\t0\\t%d\\t%d\\n",
                       \$1, r, s, \$1, r, 1000 + NR * 100 + r
        }' ${panel} >> jobs.tsv

    if [ "\$(wc -l < jobs.tsv)" -lt 2 ]; then
        echo "ERROR: stub produced no job rows from ${panel}" >&2; exit 1
    fi

    printf 'arm\\tmode\\tstructures\\nbaseline\\tfree\\t%s\\n' "\$struct0" > arms.tsv
    echo '{"stub": true}' > arms_summary.json
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
