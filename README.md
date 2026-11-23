# cgas-oligo-sim

Occupancy site calling, representative-pose selection and visualisation for
cGAS binding on candidate oligos.

**This is the cluster pipeline.** A laptop branch is maintained separately and
is not part of this tree.

---

## Status

| Component | State |
| - | - |
| `CALL_SITES` | implemented, tested on synthetic fixtures |
| `CLUSTER_POSES` | implemented, tested on synthetic fixtures |
| `MAKE_FIGURES` | implemented (figures 1, 2, 5), tested |
| `BUILD_INDEX` | implemented, tested |
| `CG_FREE_ASSOC` | **stub only** |
| `BUILD_NULL_PANEL` | **stub only** |
| Figure 4 (fold geometry) | not implemented - needs bend-angle output from the CG tier |
| Figure 6 (ChIP overlay) | not implemented - needs the ChIP dataset |

### Not implemented, and why

`CG_FREE_ASSOC` and `BUILD_NULL_PANEL` raise an error if run. Both depend on
open scientific decisions, and neither command will be written from memory:

- **CG engine.** GENESIS (AICG2+ / 3SPN.2C), CafeMol, COFFEE, or an
  IDR-capable model if the full-length construct is used. The truncated and
  full-length constructs cannot share a potential, so they are separate
  experiments with separate nulls.
- **Shuffle method.** Mononucleotide vs dinucleotide-preserving, and which
  algorithm. The sequence-preference test is only interpretable if this is
  chosen deliberately: a preference that survives a dinucleotide-preserving
  shuffle is evidence of base readout; one that does not is mechanics.

Until then, supply `--occupancy` and `--null_occupancy` directly.

---

## Before the first run

1. **Set `params.container_python`.** It has no default. Container tags could
   not be verified when this was written, so none is asserted. Pull an image
   with `python3`, `numpy`, `pandas` and `matplotlib`, pin it to an exact tag
   or digest, and set it. Never `:latest`.
2. **Confirm the GPU GRES string** in the `gpu` label in `nextflow.config`.
   It is marked `# UNVERIFIED:` and is a guess at `--gres=gpu:1`. The `gpu`
   label is unused by this branch but will be wrong for the simulation tiers
   if left unchecked.
3. **Choose `params.rmsd_cutoff` from your own data.** It has no default, on
   purpose. A cutoff below the data's noise floor silently yields one cluster
   per frame; above it, one cluster holding everything. Both produce a
   confident-looking PDB. Run the diagnostic first:

   ```
   cluster_poses.py --npz <a representative site>.npz --site-id X \
       --sequence-id Y --contacts contacts.tsv --min-frames 50 \
       --rmsd-cutoff 1 --report-only \
       --out-pdb /dev/null --out-clusters /dev/null --out-summary rmsd_report.json
   ```

   Pick a cutoff on the plateau where cluster count is low *and* the top
   cluster is a real majority.

---

## Parameters

| Param | Default | Meaning |
| - | - | - |
| `null_percentile` | 90 | percentile of the matched shuffle null used as the occupancy threshold |
| `dwell_min` | 0.05 | within one replicate, the fraction of post-equilibration frames with cGAS resident at a position for that replicate to count as occupying it |
| `replicate_min` | **unset** | optional gate: fraction of replicates that must occupy a site before it is rendered. Unset means no gate - every site clearing the null is rendered, and `replicate_frac` is reported per site so a value can be chosen from real data |
| `rmsd_cutoff` | **unset** | gromos clustering cutoff, Angstrom. Must be set; see above |
| `min_top_pop_frac` | 0.05 | fail a site whose largest cluster is smaller than this |
| `min_frames` | 50 | fail a site with fewer pooled frames |
| `bootstrap_n` | 2000 | bootstrap resamples for the occupancy CI |
| `seed` | 20261120 | recorded in the run report |
| `vendor_js` | false | false = 3Dmol.js from CDN; true = inline it (needs `js_bundle`) |

### Choosing `replicate_min`

`tiers.tsv` reports, for each sequence, how many positions clear the null at
replicate-fraction tiers 0.2 / 0.4 / 0.6 / 0.8, **and how many the shuffle null
itself produces at each tier**. That second column is the specificity curve:
pick the tier where the null's count collapses. On the synthetic fixture the
null yields ~6 positions per shuffle at 0.2, 0.2 at 0.4, and 0.0 at 0.6+.
Do this against real data, not the fixture.

---

## Input contracts

`occupancy` (TSV): `sequence_id, replicate, position, dwell_frac`
`null_occupancy` (TSV): `sequence_id, shuffle_id, replicate, position, dwell_frac`
`poses` (TSV samplesheet): `sequence_id, site_id, npz_path`
`contacts` (TSV): `sequence_id, site_id, residue, position, contact_frac`

Pose NPZ must contain `coords` (frames, beads, 3), `bead_labels` (DNA beads
prefixed `DNA`), and `replicate` (one per frame).

### Coordinates

`position` is a **1-based index into the oligo sequence** naming the anchor of
the cGAS register. It is not a genomic coordinate. `sites.tsv` also carries
`bed_start` / `bed_end` as **0-based half-open** so the table can be
intersected with ChIP output without an implicit conversion. Anything crossing
into genomic space must be lifted onto the project's genome build first; that
is not done here.

### Site identity

`site_id` comes from the sites table, never from filename order. Pose
extraction must be keyed on it. `BUILD_INDEX` fails if any rendered site has no
structure, or any structure matches no site.

---

## Failure semantics

Hard failures (the upstream step did not produce what it claims):
empty input; `dwell_frac` outside [0,1]; a sequence with no matched shuffle
null; all-zero occupancy; zero contact rows for a site; fewer pooled frames
than `min_frames`; a clustering cutoff below the data's noise floor; a
site/pose join that drops a site.

**Not** a failure: zero called sites. That is the null result. It exits 0,
writes a header-only `sites.tsv`, says so on stdout, and the HTML index
explains it.

Flagged, not failed: a single cluster holding every frame. That can be genuine
convergence or a cutoff that absorbed all variation, so it is surfaced in the
summary, on stdout, and on the card.

---

## Running

```bash
# From the project directory on BeeGFS. Never from NFS: no file locking.
sbatch run_pipeline.sh
```

Head-job script: see the project's standard header. Publish `--outdir` to NFS
so the report opens from the laptop.

First run, before trusting anything:

```bash
nextflow lint .
nextflow run main.nf -profile slurm,stub_local -stub-run \
    --occupancy x --null_occupancy x --poses x --contacts x --outdir results
head -20 work/<hash>/.command.run   # expect -p, --ntasks, --cpus-per-task, --mem; no -t
```

## Tests

```bash
bash tests/run_tests.sh
```

Generates synthetic fixtures, then runs the full downstream chain on them:
site calling against a shuffle null, gromos pose clustering, figures and the
interactive report. Two sites are planted in `SEQ_A` and `SEQ_B` is left as
background, so both the positive and the null path are exercised.

This exercises the `bin/` scripts directly. It is **not** a test of Nextflow,
Slurm or Singularity. Synthetic data calibrates the method, never the biology.

`tests/examples/` holds hand-written exemplars of the expected input shapes -
an oligo FASTA, the lineage metadata TSV, a minimal-MEME motif file and an
IUPAC consensus list. Those are tracked; generated fixtures are not.

---

## Building the panel

The usual entry point. Your FASTA carries its own coordinates, so no lookup
against another project is needed.

```bash
python3 bin/build_panel.py \
    --fasta   my_oligos.fasta \
    --mutants mutants.tsv \
    --out-fasta oligos.fasta --out-metadata metadata.tsv --out-summary panel.json
```

`--mutants` is optional; without it no mutant is generated and the panel is
the FASTA alone.

### FASTA header

The oligo name is the first token; the rest is annotation.

```
>rank1 chr1:130352560-130352639 summit=130352599 track=HSV-1 Infected THP1 Cells (Cytosolic Fraction - 1 HPI)
```

| field | meaning |
| - | - |
| `rank1` | the oligo name used everywhere downstream |
| `chr1:130352560-130352639` | locus, **0-based half-open**, + strand |
| `summit=` | optional, genomic coordinate |
| `track=` | optional free text, runs to end of line |

`end - start` must equal the sequence length, and that is asserted per record.
All sequences are taken to be + strand, so no reverse-complement handling
exists anywhere in the pipeline. Contigs must be UCSC-style: a RefSeq
(`NC_060925.1`) or GenBank (`CP068277.2`) name is rejected, because it would
fail to join the T2T BAMs later.

`lineage=`, `role=` and `mutation=` are also honoured as header keys, so a
control or a pre-grouped lineage can be declared in the FASTA without a
separate table. Any other `key=value` is kept verbatim in `notes`.

`build_panel.py` holds the only FASTA-header parser in the repository. Do not
add a second one: two parsers disagree eventually, and the disagreement
surfaces as a coordinate that is wrong rather than as an error.

### Mutant table

Three columns. One row per mutant.

| column | meaning |
| - | - |
| `name` | label for the mutant in results |
| `oligo` | name of the parent sequence it derives from |
| `mutation` | `30T>G`, or several separated by `;`: `30T>G;50A>C` |

**Mutation positions are 1-based into the oligo**, because "the 30th base
pair" is an ordinal, while locus coordinates are 0-based half-open. Two
conventions in one panel is where off-by-one errors breed, so every
mutation's reference base is checked against the parent. If `30T>G` names a
position that is not `T`, the run fails and reports whether the 0-based
reading would have matched:

```
ERROR: bad1: 30T>A says the reference base at position 30 of rank1 is 'T',
but it is 'A'. Position 30 is 'A', but 31 is 'T': this looks like a 0-based
position in a 1-based field.
```

Also fatal: an unknown parent, a position past the end of the parent, the
same position mutated twice (the result would depend on order), a mutation
changing a base to itself, and a mutant name colliding with a parent.

Mutants carry no coordinates of their own. They take `lineage_id` = the
parent and inherit its locus, which is the only correct placement for a
sequence that no longer matches the genome. The genomic position of every
change is still recorded in the summary JSON for audit.
