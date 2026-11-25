// Combine the per-arm site calls into one table and compare the arms.
//
// Arms are called separately so call_sites.py never pools across them; this
// puts them back together with `arm` as a column, which is what makes the
// single combined report possible.
//
// Arm labels are passed alongside their files as arm=path pairs, so the script
// never has to infer which arm a file belongs to from its name.
process CROSS_ARM_REPORT {
    label 'render'

    container params.container_python

    input:
    val site_specs
    val profile_specs
    path sites_files
    path profile_files

    output:
    path "sites_all_arms.tsv"      , emit: sites
    path "cross_arm_comparison.tsv", emit: comparison
    path "figures/*"               , emit: figures
    path "cross_arm_summary.json"  , emit: summary
    path "versions.yml"            , emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p figures

    cross_arm_report.py \\
        --sites ${site_specs.collect { "'" + it + "'" }.join(' ')} \\
        --profiles ${profile_specs.collect { "'" + it + "'" }.join(' ')} \\
        --out-sites sites_all_arms.tsv \\
        --out-comparison cross_arm_comparison.tsv \\
        --out-figures figures/cross_arm.png \\
        --out-summary cross_arm_summary.json

    # A combined table with no rows means every arm called nothing, which is a
    # result worth seeing rather than an empty report published as success.
    if [ "\$(wc -l < sites_all_arms.tsv)" -lt 2 ]; then
        echo "ERROR: no site rows across any arm" >&2
        exit 1
    fi

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 --version | awk '{print \$2}')
	    matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    mkdir -p figures && touch figures/cross_arm.png
    printf 'arm\\tsequence_id\\tsite_id\\n' > sites_all_arms.tsv
    printf 'baseline\\tSTUB\\tSTUB_s1\\n' >> sites_all_arms.tsv
    printf 'sequence_id\\tsite_id\\tn_arms_calling\\n' > cross_arm_comparison.tsv
    echo '{"stub": true}' > cross_arm_summary.json
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
