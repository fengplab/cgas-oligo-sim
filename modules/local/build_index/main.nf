process BUILD_INDEX {
    label 'render'

    container params.container_python

    input:
    path sites
    path cluster_summaries
    path pdbs
    path figures

    output:
    path "report/index.html", emit: index
    path "report/figures/*" , emit: figures, optional: true
    path "versions.yml"     , emit: versions

    script:
    def vendor = params.vendor_js ? 'true' : 'false'
    def bundle = params.js_bundle ? "--js-bundle ${params.js_bundle}" : ''
    def runp = groovy.json.JsonOutput.toJson([
        null_percentile: params.null_percentile,
        dwell_min      : params.dwell_min,
        replicate_min  : params.replicate_min == null ? 'unset' : params.replicate_min,
        rmsd_cutoff    : params.rmsd_cutoff,
        seed           : params.seed
    ])
    """
    set -euo pipefail
    mkdir -p report/figures
    for f in ${figures}; do cp "\$f" report/figures/ ; done

    build_index.py \\
        --sites             ${sites} \\
        --cluster-summaries ${cluster_summaries} \\
        --pdbs              ${pdbs} \\
        --figdir            report/figures \\
        --out               report/index.html \\
        --vendor-js         ${vendor} ${bundle} \\
        --run-params        '${runp}'

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    mkdir -p report/figures
    echo '<!DOCTYPE html><html><head><title>stub</title></head><body>stub</body></html>' > report/index.html
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
