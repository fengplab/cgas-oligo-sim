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

    // Everything that varies per arm travels IN the tuple. Declaring cgas_pdb
    // and pfm as separate inputs made them broadcast value channels, so every
    // arm would have received the same structure - which silently defeats the
    // whole point of a structure arm.
    input:
    tuple val(sequence_id), val(lineage_id), val(arm), val(mode), val(motif),
          path(seq_fasta), path(cgas_pdb), path(pfm)

    output:
    tuple val(sequence_id), val(lineage_id), val(arm),
          path("${prefix}_system.top"), path("${prefix}_init.gro"),
          path("top/*")                                   , emit: system
    tuple val(sequence_id), val(arm), path("${prefix}_bead_map.tsv"), emit: bead_map
    path "${prefix}_setup.log"                             , emit: log
    path "versions.yml"                                    , emit: versions

    script:
    // The same oligo appears in several arms, so the arm goes in the filename.
    // Sanitised because arm labels contain ':' (e.g. "structure:6CT9").
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}"
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
        --pwmcos-shift ${params.pwmcos_shift} 2>&1 | tee -a ${prefix}_setup.log
""" : """
    echo "PWMcos OFF: any sequence preference must emerge from 3SPN.2C" \\
        | tee -a ${prefix}_setup.log
"""
    """
    set -euo pipefail

    ${jl_env}

    ${jl} ${tool}/tools/modeling/DNA_general/build_dna.jl \\
        -s ${seq_fasta} -C -o ${prefix}_dna 2>&1 | tee ${prefix}_setup.log

    ${jl} ${tool}/src/aa_2_cg.jl ${cgas_pdb} \\
        --force-field-protein AICG2+ \\
        --force-field-DNA 3SPN.2C \\
        --output-name cgas 2>&1 | tee -a ${prefix}_setup.log
${pwmcos_block}
    mkdir -p top
    cp ${prefix}_dna_cg.itp cgas_cg.itp top/

    # UNVERIFIED: the combined .top is assembled here, not by the tool.
    # [molecules] names must match [moleculetype] in each .itp, and their order
    # must match the .gro. Chain indices follow that order: 1,2 = DNA strands,
    # 3 = cGAS.
    {
      echo '#include "atom_types.itp"'
      echo '#include "flexible_local_angle.itp"'
      echo '#include "flexible_local_dihedral.itp"'
      echo '#include "pair_energy_MJ_96.itp"'
      echo '#include "./top/${prefix}_dna_cg.itp"'
      echo '#include "./top/cgas_cg.itp"'
      echo
      echo '[ system ]'
      echo '${sequence_id} cGAS-DNA'
      echo
      echo '[ molecules ]'
      echo '${prefix}_dna_cg 1'
      echo 'cgas_cg 1'
      echo
      echo '[ cg_ele_chain_pairs ]'
      echo 'ON 1 - 3 : 1 - 3'
    } > ${prefix}_system.top

    if [ "${params.pwmcos}" = "true" ]; then
      { echo; echo '[ pwmcos_chain_pairs ]'; echo 'ON 1 - 2 : 3 - 3'; } \\
        >> ${prefix}_system.top
    fi

    # Concatenating .gro files does not work: each has a title line, an atom
    # count and a trailing box line, so `cat` buries those in the middle and
    # leaves atom indices restarting. merge_cg_system.py writes one valid .gro
    # with continuous indices, and emits the matching bead map from the same
    # parse so the two cannot disagree about bead order.
    merge_cg_system.py \\
        --dna ${prefix}_dna_cg.gro \\
        --protein cgas_cg.gro \\
        --box-angstrom ${params.cg_box_size} \\
        --title "${sequence_id} cGAS-DNA CG system" \\
        --out-gro ${prefix}_init.gro \\
        --out-bead-map ${prefix}_bead_map.tsv 2>&1 | tee -a ${prefix}_setup.log

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    genesis_cg_tool: "UNVERIFIED - record the checkout commit"
	    julia: \$(${jl_version} --version | awk '{print \$3}')
	VERSIONS
    """

    stub:
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}"
    """
    set -euo pipefail
    mkdir -p top && touch top/stub.itp
    echo stub > ${prefix}_system.top
    echo stub > ${prefix}_init.gro
    printf 'bead_index\\trole\\tnt_index\\tstrand\\tcopy\\n0\\tdna\\t1\\t1\\t1\\n' \\
        > ${prefix}_bead_map.tsv
    echo stub > ${prefix}_setup.log
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
