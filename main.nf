#!/usr/bin/env nextflow

/*
 * cgas-oligo-sim - CG simulation, occupancy site calling, pose clustering and
 * visualisation, in one execution.
 *
 * This is the CLUSTER pipeline. A laptop branch is maintained separately.
 *
 * Two entry points into the same downstream analysis:
 *
 *   params.oligos set  -> the simulation tier runs. The arm matrix is generated
 *                         here, every (arm, oligo, replicate) is submitted as an
 *                         independent task, and occupancy and contacts are
 *                         produced by the pipeline.
 *   params.oligos null -> the simulation tier is skipped and params.occupancy
 *                         and params.contacts are required instead, which is the
 *                         visualisation-only path.
 *
 * Still external by design: the null panel (params.null_occupancy) and the pose
 * ensembles (params.poses). See README 'Not implemented'.
 */

include { ENUMERATE_ARMS    } from './modules/local/enumerate_arms/main.nf'
include { GENESIS_SETUP     } from './modules/local/genesis_setup/main.nf'
include { GENESIS_CG_RUN    } from './modules/local/genesis_cg_run/main.nf'
include { EXTRACT_OCCUPANCY } from './modules/local/extract_occupancy/main.nf'
include { EXTRACT_POSES     } from './modules/local/extract_poses/main.nf'
include { CALL_SITES        } from './modules/local/call_sites/main.nf'
include { CROSS_ARM_REPORT  } from './modules/local/cross_arm_report/main.nf'
include { CLUSTER_POSES     } from './modules/local/cluster_poses/main.nf'
include { MAKE_FIGURES      } from './modules/local/make_figures/main.nf'
include { BUILD_INDEX       } from './modules/local/build_index/main.nf'

def requireParam(name, value, hint) {
    if (value == null) {
        error "params.${name} is required. ${hint}"
    }
    return value
}

// A samplesheet row is only usable if every column it is read for is present.
// Failing here names the row; failing later names a null path.
def requireCols(row, cols, what) {
    def missing = cols.findAll { c -> !row[c] }
    if (missing) {
        error "${what} row is missing ${missing.join(', ')}: ${row}"
    }
    return row
}

workflow {

    main:

    requireParam('outdir', params.outdir,
                 'Publish directory on NFS, so results open from the laptop.')
    requireParam('null_occupancy', params.null_occupancy,
                 'Long-format TSV of the matched shuffle null. The occupancy ' +
                 'threshold is defined against it, so it is not optional. ' +
                 'BUILD_NULL_PANEL is not wired in yet, so supply it directly.')
    if (!params.oligos) {
        requireParam('poses', params.poses,
                     'Samplesheet TSV: sequence_id, site_id, npz_path. ' +
                     'Required only when params.oligos is unset, because then ' +
                     'nothing in the run extracts poses.')
    }
    requireParam('rmsd_cutoff', params.rmsd_cutoff,
                 'Pose-clustering cutoff in Angstrom. No default: run ' +
                 'cluster_poses.py --report-only and choose from the sweep.')

    ch_null = channel.fromPath(params.null_occupancy, checkIfExists: true)

    if (params.oligos) {

        requireParam('panel', params.panel,
                     'TSV with a sequence_id column: the oligo panel that ' +
                     'defines the arm matrix.')
        requireParam('structures', params.structures,
                     'Samplesheet TSV with a structure_id column; see ' +
                     'conf/structures.template.tsv')
        requireParam('prepared_structures', params.prepared_structures,
                     'TSV: structure_id, pdb. Output of prepare_structures.py, ' +
                     'which is run by hand - it is not in the workflow.')

        // ---- arm matrix ----
        // An absent motif table is the no-motif case and must still stage a
        // file, because the process declares a path input. ENUMERATE_ARMS tests
        // for this placeholder by name.
        ch_no_file = channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)
        ch_motifs  = params.motif_pfms
            ? channel.fromPath(params.motif_pfms, checkIfExists: true)
            : ch_no_file

        ENUMERATE_ARMS(
            channel.fromPath(params.panel,      checkIfExists: true),
            channel.fromPath(params.structures, checkIfExists: true),
            ch_motifs
        )

        // One row per trajectory. The null half of the matrix is filtered out
        // here rather than never generated, so the matrix stays the single
        // description of the campaign and the null can be switched on later by
        // relaxing this filter alone.
        ch_jobs = ENUMERATE_ARMS.out.jobs
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.panel == 'observed' }

        // ---- inputs addressed by identity, never by file order ----
        ch_oligos = channel.fromPath(params.oligos, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                requireCols(row, ['sequence_id', 'fasta'], 'oligos samplesheet')
                tuple(row.sequence_id,
                      row.lineage_id ?: row.sequence_id,
                      file(row.fasta, checkIfExists: true))
            }

        ch_structures = channel.fromPath(params.prepared_structures, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                requireCols(row, ['structure_id', 'pdb'], 'prepared_structures samplesheet')
                tuple(row.structure_id, file(row.pdb, checkIfExists: true))
            }

        ch_pfms = params.motif_pfms
            ? channel.fromPath(params.motif_pfms, checkIfExists: true)
                .splitCsv(header: true, sep: '\t')
                .map { row ->
                    requireCols(row, ['family', 'pfm'], 'motif_pfms samplesheet')
                    tuple(row.family, file(row.pfm, checkIfExists: true))
                }
            : channel.empty()

        // The multibasin arm names SEVERAL structures in one field
        // ("a,b"), and GENESIS_SETUP builds from a single reference structure.
        // That id matches no structure_id, so combine(by:) below would delete
        // the whole arm without a word. Drop it loudly instead.
        ch_jobs_split = ch_jobs.branch { r ->
            multibasin: r.structures.contains(',')
            single:     true
        }

        ch_jobs_split.multibasin
            .map { r -> r.arm }
            .unique()
            .subscribe { arm ->
                log.warn "Skipping arm '${arm}': it names several reference " +
                         "structures (multi-basin), which GENESIS_SETUP does not " +
                         "build. Pass --multibasin false to stop generating it."
            }

        ch_jobs_obs = ch_jobs_split.single

        // ---- setup runs once per (oligo, arm), not per replicate ----
        ch_setup_keys = ch_jobs_obs
            .map { r -> tuple(r.sequence_id, r.arm, r.mode, r.motif, r.structures) }
            .unique()

        // combine(by: 0), not join(by: 0): an oligo appears in several arms, so
        // the keys are not unique on this side and join would silently keep one
        // pairing per key and drop the rest.
        ch_setup = ch_setup_keys
            .map { seq, arm, mode, motif, struct -> tuple(seq, arm, mode, motif, struct) }
            .combine(ch_oligos, by: 0)
            .map { seq, arm, mode, motif, struct, lineage, fasta ->
                tuple(struct, seq, arm, mode, motif, lineage, fasta)
            }
            .combine(ch_structures, by: 0)
            .map { _struct, seq, arm, mode, motif, lineage, fasta, pdb ->
                tuple(motif, seq, arm, mode, lineage, fasta, pdb)
            }
            // A motif arm needs its PFM; every other arm stages the placeholder.
            .branch { motif, _seq, _arm, _mode, _lineage, _fasta, _pdb ->
                with_motif: motif as boolean
                without:    true
            }

        ch_setup_motif = ch_setup.with_motif
            .combine(ch_pfms, by: 0)
            .map { _motif, seq, arm, mode, lineage, fasta, pdb, pfm ->
                tuple(seq, lineage, arm, mode, _motif, fasta, pdb, pfm)
            }

        ch_setup_plain = ch_setup.without
            .combine(ch_no_file)
            .map { _motif, seq, arm, mode, lineage, fasta, pdb, no_file ->
                tuple(seq, lineage, arm, mode, '', fasta, pdb, no_file)
            }

        ch_setup_all = ch_setup_motif.mix(ch_setup_plain)

        // Every join above is a combine(by:), which drops a key that the other
        // side lacks instead of failing. One guard covers all of them: if the
        // number of (oligo, arm) setups does not survive the joins, a
        // sequence_id, structure_id or motif family referenced by the arm
        // matrix is missing from a samplesheet, and those arms would simply
        // never run while the pipeline reported success.
        ch_setup_keys.count()
            .combine(ch_setup_all.count())
            .subscribe { expected, got ->
                if (got != expected) {
                    error "setup lost ${expected - got} of ${expected} " +
                          "(oligo, arm) combinations while joining the " +
                          "samplesheets. A sequence_id, structure_id or motif " +
                          "family named in the arm matrix is absent from " +
                          "params.oligos, params.prepared_structures or " +
                          "params.motif_pfms."
                }
            }

        GENESIS_SETUP(ch_setup_all)

        // ---- one task per job row, all submittable at once ----
        // Keyed on [sequence_id, arm] so a replicate can never be attached to
        // another arm's topology.
        ch_runs = GENESIS_SETUP.out.system
            .map { seq, lineage, arm, top, gro, topdir ->
                tuple([seq, arm], lineage, top, gro, topdir)
            }
            .combine(
                ch_jobs_obs.map { r -> tuple([r.sequence_id, r.arm], r.replicate, r.seed) },
                by: 0
            )
            .map { key, lineage, top, gro, topdir, replicate, seed ->
                tuple(key[0], lineage, key[1], replicate as Integer, seed as Integer,
                      top, gro, topdir)
            }

        GENESIS_CG_RUN(ch_runs)

        ch_occ_in = GENESIS_CG_RUN.out.traj
            .map { seq, lineage, arm, replicate, dcd ->
                tuple([seq, arm], lineage, replicate, dcd)
            }
            .combine(
                GENESIS_SETUP.out.bead_map.map { seq, arm, bm -> tuple([seq, arm], bm) },
                by: 0
            )
            .map { key, lineage, replicate, dcd, bead_map ->
                tuple(key[0], lineage, key[1], replicate, dcd, bead_map)
            }

        EXTRACT_OCCUPANCY(ch_occ_in)

        // Occupancy is grouped BY ARM and called per arm. call_sites.py takes
        // no arm column and computes n_rep from replicate labels, which repeat
        // across arms, so one pooled table would average different force
        // fields and inflate replicate_frac toward 1. The single combined
        // report is assembled afterwards by CROSS_ARM_REPORT.
        ch_occ_by_arm = EXTRACT_OCCUPANCY.out.occupancy
            .map { _seq, arm, f -> tuple(arm, f) }
            .groupTuple()

        // Contacts are not arm-partitioned downstream (cluster_poses reads one
        // table keyed by sequence_id and site_id), so they are pooled, with the
        // arm retained as a column in the file itself.
        ch_contacts = EXTRACT_OCCUPANCY.out.contacts
            .map { _seq, _arm, f -> f }
            .collectFile(name: 'contacts_all.tsv', keepHeader: true, skip: 1,
                         sort: true, storeDir: "${params.outdir}/occupancy")

        // poses, cross-arm report and the figure steps are wired after the
        // branch, because CALL_SITES must be instantiated exactly once.
        ch_sim_versions = ENUMERATE_ARMS.out.versions
            .mix(GENESIS_SETUP.out.versions)
            .mix(GENESIS_CG_RUN.out.versions)
            .mix(EXTRACT_OCCUPANCY.out.versions)

    } else {

        requireParam('occupancy', params.occupancy,
                     'Long-format TSV: sequence_id, replicate, position, ' +
                     'dwell_frac. Required when params.oligos is not set, ' +
                     'because then nothing in the run produces it.')
        requireParam('contacts', params.contacts,
                     'Long-format TSV: sequence_id, site_id, residue, ' +
                     'position, contact_frac. Required when params.oligos is ' +
                     'not set.')

        // One pseudo-arm, so the arm-scoped CALL_SITES signature serves both
        // paths and the process is still instantiated exactly once.
        ch_occ_by_arm = channel.fromPath(params.occupancy, checkIfExists: true)
            .map { f -> tuple('external', [f]) }
        ch_contacts   = channel.fromPath(params.contacts, checkIfExists: true)
        ch_sim_versions = channel.empty()
    }

    // Called once for both paths: a process may be invoked only once per
    // workflow, so this cannot sit inside either branch.
    CALL_SITES(ch_occ_by_arm, ch_null.first())

    if (params.oligos) {

        // ---- poses, per (arm, sequence_id), after that arm's sites exist ----
        // Trajectories are grouped so one task sees every replicate of a pair,
        // and the replicate labels travel with them so the pairing is explicit
        // rather than positional.
        ch_traj_grouped = GENESIS_CG_RUN.out.traj
            .map { seq, _lineage, arm, replicate, dcd ->
                tuple(tuple(arm, seq), replicate, dcd)
            }
            .groupTuple()

        ch_pose_in = ch_traj_grouped
            .combine(
                GENESIS_SETUP.out.bead_map
                    .map { seq, arm, bm -> tuple(tuple(arm, seq), bm) },
                by: 0
            )
            .map { key, reps, dcds, bm -> tuple(key[0], key[1], reps, dcds, bm) }
            // Sites are per arm, so join on the arm alone: every sequence in
            // that arm reads the same sites table.
            .combine(CALL_SITES.out.sites, by: 0)

        EXTRACT_POSES(ch_pose_in)

        ch_poses = EXTRACT_POSES.out.manifest
            .map { _arm, f -> f }
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                requireCols(row, ['sequence_id', 'site_id', 'npz_path'],
                            'extracted poses manifest')
                tuple(row.sequence_id, row.site_id,
                      file(row.npz_path, checkIfExists: true))
            }

        // ---- one report across all arms ----
        // arm=path pairs, so the script never infers an arm from a filename.
        CROSS_ARM_REPORT(
            CALL_SITES.out.sites.map   { arm, f -> "${arm}=${f.name}" }.collect(),
            CALL_SITES.out.profile.map { arm, f -> "${arm}=${f.name}" }.collect(),
            CALL_SITES.out.sites.map   { _arm, f -> f }.collect(),
            CALL_SITES.out.profile.map { _arm, f -> f }.collect()
        )

        ch_sites       = CROSS_ARM_REPORT.out.sites
        ch_extra_vers  = EXTRACT_POSES.out.versions
            .mix(CROSS_ARM_REPORT.out.versions)

    } else {

        ch_poses = channel.fromPath(params.poses, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                requireCols(row, ['sequence_id', 'site_id', 'npz_path'],
                            'poses samplesheet')
                tuple(row.sequence_id, row.site_id,
                      file(row.npz_path, checkIfExists: true))
            }
        // .first() here, not at the consumers: this makes ch_sites a value
        // channel in BOTH branches, so it can feed two processes.
        ch_sites      = CALL_SITES.out.sites.map { _arm, f -> f }.first()
        ch_extra_vers = channel.empty()
    }

    // One occupancy and one profile file for the figure steps, which take a
    // single table each. The per-arm detail lives in the cross-arm comparison.
    ch_occupancy_one = ch_occ_by_arm.map { _arm, fs -> fs }.flatten().first()
    ch_profile_one   = CALL_SITES.out.profile.map { _arm, f -> f }.first()

    CLUSTER_POSES(ch_poses, ch_contacts.first())

    MAKE_FIGURES(
        ch_profile_one,
        ch_sites,
        ch_occupancy_one,
        ch_contacts.first()
    )

    BUILD_INDEX(
        ch_sites,
        CLUSTER_POSES.out.summary.map { _id, f -> f }.collect(),
        CLUSTER_POSES.out.pdb.map     { _id, f -> f }.collect(),
        MAKE_FIGURES.out.figures.collect()
    )

    // No emit section: this is the entry workflow, and only a NAMED workflow may
    // emit. Outputs reach the user via the process-scope publishDir default in
    // nextflow.config (mode 'copy', to params.outdir on NFS). That default is
    // global, so grepping this module for publishDir finds nothing.
    ch_sim_versions
        .mix(CALL_SITES.out.versions)
        .mix(ch_extra_vers)
        .mix(CLUSTER_POSES.out.versions)
        .mix(MAKE_FIGURES.out.versions)
        .mix(BUILD_INDEX.out.versions)
        .collectFile(name: 'versions.yml', storeDir: "${params.outdir}/pipeline_info")
}
