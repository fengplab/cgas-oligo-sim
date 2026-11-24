#!/usr/bin/env nextflow

/*
 * cgas-oligo-sim - occupancy site calling, pose clustering and visualisation.
 *
 * This is the CLUSTER pipeline. A laptop branch is maintained separately.
 *
 * Implemented here: the visualisation branch (CALL_SITES -> CLUSTER_POSES ->
 * MAKE_FIGURES -> BUILD_INDEX). The upstream simulation tiers are stubs; see
 * README 'Not implemented'.
 */

include { CALL_SITES    } from './modules/local/call_sites/main.nf'
include { CLUSTER_POSES } from './modules/local/cluster_poses/main.nf'
include { MAKE_FIGURES  } from './modules/local/make_figures/main.nf'
include { BUILD_INDEX   } from './modules/local/build_index/main.nf'

def requireParam(name, value, hint) {
    if (value == null) {
        error "params.${name} is required. ${hint}"
    }
    return value
}

workflow {

    main:

    requireParam('occupancy',      params.occupancy,
                 'Long-format TSV: sequence_id, replicate, position, dwell_frac')
    requireParam('null_occupancy', params.null_occupancy,
                 'Long-format TSV of the matched shuffle panel. The occupancy ' +
                 'threshold is defined against it, so it is not optional.')
    requireParam('poses',          params.poses,
                 'Samplesheet TSV: sequence_id, site_id, npz_path')
    requireParam('contacts',       params.contacts,
                 'Long-format TSV: sequence_id, site_id, residue, position, contact_frac')
    requireParam('outdir',         params.outdir,
                 'Publish directory on NFS, so results open from the laptop.')

    ch_occupancy = channel.fromPath(params.occupancy,      checkIfExists: true)
    ch_null      = channel.fromPath(params.null_occupancy, checkIfExists: true)
    ch_contacts  = channel.fromPath(params.contacts,       checkIfExists: true)

    CALL_SITES(ch_occupancy, ch_null)

    // Site identity comes from the samplesheet, never from filename order.
    ch_poses = channel
        .fromPath(params.poses, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.sequence_id || !row.site_id || !row.npz_path) {
                error "poses samplesheet row is missing a required column " +
                      "(sequence_id, site_id, npz_path): ${row}"
            }
            tuple(row.sequence_id, row.site_id, file(row.npz_path, checkIfExists: true))
        }

    CLUSTER_POSES(ch_poses, ch_contacts.first())

    MAKE_FIGURES(
        CALL_SITES.out.profile,
        CALL_SITES.out.sites,
        ch_occupancy.first(),
        ch_contacts.first()
    )

    BUILD_INDEX(
        CALL_SITES.out.sites,
        CLUSTER_POSES.out.summary.map { _id, f -> f }.collect(),
        CLUSTER_POSES.out.pdb.map     { _id, f -> f }.collect(),
        MAKE_FIGURES.out.figures.collect()
    )

    // No emit section: this is the entry workflow, and only a NAMED workflow may
    // emit. The outputs this used to declare - CALL_SITES.out.sites and
    // BUILD_INDEX.out.index - still reach the user, via the process-scope
    // publishDir default in nextflow.config (mode 'copy', to params.outdir on
    // NFS). That default is global, so grepping this module for publishDir
    // finds nothing; it is set once in nextflow.config for every process.
    // If the pipeline is ever composed into another workflow or given a
    // workflow-level test, move this body into a named workflow and put the
    // emit section there.
    CALL_SITES.out.versions
        .mix(CLUSTER_POSES.out.versions)
        .mix(MAKE_FIGURES.out.versions)
        .mix(BUILD_INDEX.out.versions)
        .collectFile(name: 'versions.yml', storeDir: "${params.outdir}/pipeline_info")
}
