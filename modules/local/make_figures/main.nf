process MAKE_FIGURES {
    label 'render'

    container params.container_python

    input:
    path profile
    path sites
    path occupancy
    path contacts

    output:
    path "figures/*"   , emit: figures
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p figures

    make_figures.py \\
        --profile   ${profile} \\
        --sites     ${sites} \\
        --occupancy ${occupancy} \\
        --contacts  ${contacts} \\
        --outdir    figures

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
	    matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    mkdir -p figures
    touch figures/fig1_occupancy_STUB.png figures/fig1_occupancy_STUB.svg figures/fig1_occupancy_STUB.csv
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
