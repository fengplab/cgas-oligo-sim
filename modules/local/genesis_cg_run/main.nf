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
// The executable is atdyn, run SERIALLY with OpenMP threads. GENESIS is built
// with --disable-mpi, so there are no ranks: parallelism is Slurm fan-out over
// independent trajectories. This process still REPLACES the profile
// clusterOptions, so it repeats --ntasks and --cpus-per-task itself.
//
// Do not reintroduce mpirun here against a --disable-mpi build. `mpirun -np N`
// on a serial atdyn starts N independent processes that all write the same
// dcd, rst and log. That does not fail; it produces a plausible, corrupt
// trajectory.
//
// UNVERIFIED: not executed during authoring. Confirm against the GENESIS
// manual for the installed version before a production run.
process GENESIS_CG_RUN {
    label 'genesis'
    tag  "${sequence_id}/rep${replicate}"

    // Local install when genesis_bin is supplied, container otherwise.
    container { params.genesis_bin ? null : params.container_genesis }

    // Replicate and seed come from the job matrix row, not from `each replicate`
    // and not from arithmetic here. enumerate_arms.py derives each seed from
    // (global seed, arm, oligo, panel, replicate); the previous
    // `params.seed + replicate` gave every oligo and every arm at the same
    // replicate index an IDENTICAL random stream, correlating thermal noise
    // across exactly the comparison this campaign is making.
    input:
    tuple val(sequence_id), val(lineage_id), val(arm), val(replicate), val(seed),
          path(top), path(gro), path(topdir)

    output:
    tuple val(sequence_id), val(lineage_id), val(arm), val(replicate),
          path("${prefix}.dcd")                         , emit: traj
    path "${prefix}.log"                                , emit: log
    path "versions.yml"                                 , emit: versions

    script:
    if (params.genesis_bin == null && params.container_genesis == null) {
        error "GENESIS_CG_RUN needs atdyn. Set params.genesis_bin to the bin/ " +
              "of a local GENESIS install, or params.container_genesis to an " +
              "image that provides atdyn. Neither is set."
    }
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}_rep${replicate}"
    // atdyn reads box_size_x/y/z only under type = PBC, and does not fail when
    // they are missing - it silently omits the box_size line from its setup log,
    // because a box may legitimately come from a restart file. So refuse here
    // rather than start a run with an undefined volume.
    // Validate numerically, not just against Groovy null: a value supplied on the
    // command line arrives as a String, so `--cg_box_size null` is the four-letter
    // string and would sail past a null check straight into the control file.
    if (params.cg_boundary == 'PBC') {
        def box_ok = params.cg_box_size != null &&
                     params.cg_box_size.toString().isNumber() &&
                     params.cg_box_size.toString().toDouble() > 0
        if (!box_ok) {
            error "params.cg_boundary is 'PBC' but params.cg_box_size is not a " +
                  "positive number (got: ${params.cg_box_size}). The box volume sets " +
                  "the effective concentration, so it cannot be unset or zero. " +
                  "A 300 A cube is about 61 uM for one copy of each species."
        }
    }
    def boundary_box = params.cg_boundary == 'PBC' ? """box_size_x = ${params.cg_box_size}
box_size_y = ${params.cg_box_size}
box_size_z = ${params.cg_box_size}
""" : ''
    def atdyn_resolve = params.genesis_bin ?
        "ATDYN='${params.genesis_bin}/atdyn'" :
        'ATDYN="\$(command -v atdyn)"'
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
dcdfile = ${prefix}.dcd
rstfile = ${prefix}.rst

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
${boundary_box}INPEOF

    ${atdyn_resolve}
    if [ ! -x "\$ATDYN" ]; then
        echo "ERROR: atdyn not executable at '\$ATDYN'" >&2
        exit 1
    fi

    # Serial: one process, OpenMP threads only. No mpirun - see the header.
    export OMP_NUM_THREADS=${params.cg_omp_threads}
    "\$ATDYN" md.inp > ${prefix}.log 2>&1

    # A trajectory that is absent or empty is a failure, not an empty result.
    if [ ! -s "${prefix}.dcd" ]; then
        echo "ERROR: atdyn produced no trajectory for ${sequence_id} rep ${replicate}" >&2
        tail -40 ${prefix}.log >&2
        exit 1
    fi

    cat <<-VERSIONS > versions.yml
	"${task.process}":
	    atdyn_path: "\$ATDYN"
	    atdyn_sha256: \$(sha256sum "\$ATDYN" | awk '{print \$1}')
	    omp_threads: ${params.cg_omp_threads}
	    seed: ${seed}
	VERSIONS
    """

    stub:
    prefix = "${sequence_id}_${arm.replaceAll(/[^A-Za-z0-9._-]/, '_')}_rep${replicate}"
    """
    set -euo pipefail
    echo stub > ${prefix}.dcd
    echo stub > ${prefix}.log
    echo '"${task.process}": {stub: true}' > versions.yml
    """
}
