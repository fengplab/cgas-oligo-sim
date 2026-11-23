// One independent CG free-association replicate for one oligo.
//
// Control-file keys are taken verbatim from the published TBP-DNA tutorial:
//   [ENERGY]  forcefield = RESIDCG, electrostatic = CUTOFF, cg_* cutoffs,
//             cg_sol_ionic_strength, cg_pro_DNA_ele_scale_Q
//   [DYNAMICS] integrator = VVER_CG
//   [ENSEMBLE] NVT / LANGEVIN
//   [BOUNDARY] type = NOBC  (see the boundary warning below)
// The control file never names AICG2+ or 3SPN.2C; the model comes from the
// topology files.
//
// The executable is atdyn, run under MPI with OpenMP threads. GENESIS needs
// more than one task, so this process REPLACES the profile clusterOptions and
// therefore repeats --ntasks and --cpus-per-task itself.
//
// UNVERIFIED: not executed during authoring. Confirm against the GENESIS
// manual for the installed version before a production run.
process GENESIS_CG_RUN {
    label 'genesis'
    tag  "${sequence_id}/rep${replicate}"

    container params.container_genesis

    input:
    tuple val(sequence_id), val(lineage_id), path(top), path(gro), path(topdir)
    each replicate

    output:
    tuple val(sequence_id), val(lineage_id), val(replicate),
          path("${sequence_id}_rep${replicate}.dcd")    , emit: traj
    path "${sequence_id}_rep${replicate}.log"           , emit: log
    path "versions.yml"                                 , emit: versions

    script:
    // Seed is deterministic in the sequence and replicate, and is recorded.
    def seed = "${params.seed}".toInteger() + replicate
    if (params.cg_boundary == 'NOBC' && params.cg_free_association) {
        log.warn "GENESIS_CG_RUN: boundary is NOBC with free association enabled. " +
                 "cGAS can diffuse away and never return, and no concentration is " +
                 "defined. See params.cg_boundary."
    }
    """
    set -euo pipefail

    cat > md.inp <<'INPEOF'
[INPUT]
grotopfile = ${top}
grocrdfile = ${gro}

[OUTPUT]
dcdfile = ${sequence_id}_rep${replicate}.dcd
rstfile = ${sequence_id}_rep${replicate}.rst

[ENERGY]
forcefield             = RESIDCG
electrostatic          = CUTOFF
cg_cutoffdist_ele      = ${params.cg_cutoffdist_ele}
cg_cutoffdist_DNAbp    = ${params.cg_cutoffdist_DNAbp}
cg_pairlistdist_ele    = ${params.cg_pairlistdist_ele}
cg_pairlistdist_PWMcos = ${params.cg_pairlistdist_PWMcos}
cg_pairlistdist_DNAbp  = ${params.cg_pairlistdist_DNAbp}
cg_pairlistdist_exv    = ${params.cg_pairlistdist_exv}
cg_sol_ionic_strength  = ${params.cg_ionic_strength}
cg_pro_DNA_ele_scale_Q = ${params.cg_pro_DNA_ele_scale_Q}

[DYNAMICS]
integrator     = VVER_CG
nsteps         = ${params.cg_nsteps}
timestep       = ${params.cg_timestep}
rstout_period  = ${params.cg_nsteps}
crdout_period  = ${params.cg_crdout_period}
eneout_period  = ${params.cg_crdout_period}
nbupdate_period = 100
stoptr_period  = 100
iseed          = ${seed}

[CONSTRAINTS]
rigid_bond = NO

[ENSEMBLE]
ensemble    = NVT
tpcontrol   = LANGEVIN
temperature = ${params.cg_temperature}
gamma_t     = ${params.cg_gamma_t}

[BOUNDARY]
type = ${params.cg_boundary}
INPEOF

    export OMP_NUM_THREADS=${params.cg_omp_threads}
    mpirun -np ${params.cg_mpi_ranks} atdyn md.inp \\
        > ${sequence_id}_rep${replicate}.log 2>&1

    # A trajectory that is absent or empty is a failure, not an empty result.
    if [ ! -s "${sequence_id}_rep${replicate}.dcd" ]; then
        echo "ERROR: atdyn produced no trajectory for ${sequence_id} rep ${replicate}" >&2
        tail -40 ${sequence_id}_rep${replicate}.log >&2
        exit 1
    fi

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    atdyn: \$(atdyn --version 2>&1 | head -1 || echo "UNVERIFIED")
	    seed: ${seed}
	VERSIONS
    """

    stub:
    """
    set -euo pipefail
    echo stub > ${sequence_id}_rep${replicate}.dcd
    echo stub > ${sequence_id}_rep${replicate}.log
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
