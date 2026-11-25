// Build CG topology and coordinates for one oligo with the GENESIS-CG-tool.
//
// Flags are taken from the published GENESIS CG tutorials, not from memory:
//   build_dna.jl -s <fasta> -C -o <prefix>        (DNA, 3SPN.2C)
//   aa_2_cg.jl <pdb> --force-field-protein AICG2+ --force-field-DNA 3SPN.2C
//   aa_2_cg.jl <pdb> --pwmcos -p <pfm> --patch <itp> --pwmcos-scale --pwmcos-shift
//
// UNVERIFIED: nothing here was executed during authoring - GENESIS was not
// available. Run `aa_2_cg.jl -h` and `build_dna.jl -h` on the cluster and
// reconcile every flag before a production run.
process GENESIS_SETUP {
    label 'analysis'
    tag  "${sequence_id}"

    // Local install when julia_exe is supplied, container otherwise. Verified
    // in Nextflow 25.10.2: a null container runs natively even with
    // singularity.enabled = true.
    container { params.julia_exe ? null : params.container_genesis }

    input:
    tuple val(sequence_id), val(lineage_id), path(seq_fasta)
    path cgas_pdb
    path pfm

    output:
    tuple val(sequence_id), val(lineage_id),
          path("${sequence_id}_system.top"), path("${sequence_id}_init.gro"),
          path("top/*")                                   , emit: system
    path "${sequence_id}_setup.log"                        , emit: log
    path "versions.yml"                                    , emit: versions

    script:
    if (params.genesis_cg_tool == null) {
        error "params.genesis_cg_tool is not set: path to the genesis_cg_tool " +
              "checkout (ships with GENESIS >1.7 under src/analysis/cg_tools)."
    }
    if (params.julia_exe == null && params.container_genesis == null) {
        error "GENESIS_SETUP needs Julia. Set params.julia_exe to a local julia " +
              "binary, or params.container_genesis to an image that provides " +
              "julia and genesis_cg_tool. Neither is set."
    }
    def tool = params.genesis_cg_tool
    // The scripts carry a `#!/usr/bin/env julia` shebang, which needs julia on
    // PATH and ignores --project. Calling julia explicitly removes both
    // problems; inside a container the shebang is used as before.
    def jl_project = params.julia_project ?: tool
    def jl = params.julia_exe ? "${params.julia_exe} --project=${jl_project}" : ''
    def jl_env = params.julia_depot ?
        "export JULIA_DEPOT_PATH='${params.julia_depot}'\n    export JULIA_PROJECT='${jl_project}'" :
        "export JULIA_PROJECT='${jl_project}'"
    def jl_version = params.julia_exe ? params.julia_exe : 'julia'
    def pwmcos_block = params.pwmcos ? """
    # PWMcos ON. The PFM becomes a force-field term, so this arm imposes the
    # motif and cannot test whether cGAS prefers it.
    ${jl} ${tool}/src/aa_2_cg.jl ${cgas_pdb} --pwmcos -p ${pfm} \\
        --patch cgas_cg.itp \\
        --pwmcos-scale ${params.pwmcos_scale} \\
        --pwmcos-shift ${params.pwmcos_shift} 2>&1 | tee -a ${sequence_id}_setup.log
""" : """
    echo "PWMcos OFF: any sequence preference must emerge from 3SPN.2C" \\
        | tee -a ${sequence_id}_setup.log
"""
    """
    set -euo pipefail

    ${jl_env}

    ${jl} ${tool}/tools/modeling/DNA_general/build_dna.jl \\
        -s ${seq_fasta} -C -o ${sequence_id}_dna 2>&1 | tee ${sequence_id}_setup.log

    ${jl} ${tool}/src/aa_2_cg.jl ${cgas_pdb} \\
        --force-field-protein AICG2+ \\
        --force-field-DNA 3SPN.2C \\
        --output-name cgas 2>&1 | tee -a ${sequence_id}_setup.log
${pwmcos_block}
    mkdir -p top
    cp ${sequence_id}_dna_cg.itp cgas_cg.itp top/

    # UNVERIFIED: the combined .top is assembled here, not by the tool.
    # [molecules] names must match [moleculetype] in each .itp, and their order
    # must match the .gro. Chain indices follow that order: 1,2 = DNA strands,
    # 3 = cGAS.
    {
      echo '#include "atom_types.itp"'
      echo '#include "flexible_local_angle.itp"'
      echo '#include "flexible_local_dihedral.itp"'
      echo '#include "pair_energy_MJ_96.itp"'
      echo '#include "./top/${sequence_id}_dna_cg.itp"'
      echo '#include "./top/cgas_cg.itp"'
      echo
      echo '[ system ]'
      echo '${sequence_id} cGAS-DNA'
      echo
      echo '[ molecules ]'
      echo '${sequence_id}_dna_cg 1'
      echo 'cgas_cg 1'
      echo
      echo '[ cg_ele_chain_pairs ]'
      echo 'ON 1 - 3 : 1 - 3'
    } > ${sequence_id}_system.top

    if [ "${params.pwmcos}" = "true" ]; then
      { echo; echo '[ pwmcos_chain_pairs ]'; echo 'ON 1 - 2 : 3 - 3'; } \\
        >> ${sequence_id}_system.top
    fi

    cat ${sequence_id}_dna_cg.gro cgas_cg.gro > ${sequence_id}_init.gro

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    genesis_cg_tool: "UNVERIFIED - record the checkout commit"
	    julia: \$(${jl_version} --version | awk '{print \$3}')
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    mkdir -p top && touch top/stub.itp
    echo stub > ${sequence_id}_system.top
    echo stub > ${sequence_id}_init.gro
    echo stub > ${sequence_id}_setup.log
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
