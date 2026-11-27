# rosetta_gaps

Rebuild the unmodelled loops of the structures named in `config/structures.tsv`,
and write the repaired coordinates back as mmCIF for `structures.tsv` to point
at. Standalone: nothing here is a Nextflow module and `main.nf` does not call
it. It reads the same samplesheet the pipeline reads and writes a file the
pipeline can consume, and that is the whole coupling.

## Why two Rosetta stages

Next-generation KIC cannot build a residue that is not already in the pose. The
`loopmodel` documentation is explicit that the input is a *"Template pdb file
with real coordinates for all residues plus the first and last residue of each
loop region"*, and KIC moves are local perturbations in which the C-alpha of a
loop's first and last residue stays fixed. A residue absent from the crystal is
absent from the pose, so there is nothing to perturb.

So:

| stage | tool | job |
| - | - | - |
| 1 | `remodel` | insert the missing residues from a blueprint, close the chain |
| 2 | `loopmodel` | resample each inserted span with NGK, `extend` set |

Stage 2's loop file sets the extend flag, which randomises the span with
idealised geometry before sampling. Stage 1's conformation is therefore
discarded: Remodel supplies connectivity, NGK does the modelling.

All of a chain's gaps go in **one** blueprint and **one** loop file. Pose
numbering shifts by the number of residues inserted below a position, so
inserting gaps in separate runs would leave every later loop file pointing at
the wrong residues — which Rosetta cannot detect. One gap-free pose also means
one output structure carries every rebuilt loop, instead of several structures
that each rebuilt a different one.

## Three numbering systems

Never mixed implicitly, because this is where this kind of code goes wrong:

| name | what it is | who uses it |
| - | - | - |
| author | the numbering in the deposited file; may be non-contiguous | `config/structures.tsv`, `gaps.tsv`, the repaired mmCIF |
| input pose | 1..N over **observed** residues only | the blueprint, Remodel |
| full pose | 1..M over observed **+** inserted | the loop file, loopmodel |

`numbering.tsv` records all three for every residue, and every loop file is
checked against it before it is written. Residue indices in a Rosetta loop file
are in Rosetta numbering, which runs continuously across chains — which is why
the protein chain is written first.

## Run order

Everything runs on the **cluster**, from the project directory on BeeGFS. Only
`promote.py` writes somewhere the laptop can read.

```bash
# 0. What is actually in the image. Do this first; the flag files are written
#    against its answer, not against a guess.
module load singularity
bash rosetta_gaps/bin/probe_container.sh containers/rosetta.sif | tee logs/probe.txt

# 1. Enumerate the gaps of the UNCOMMENTED rows only (commented rows are skipped).
python3 rosetta_gaps/bin/find_gaps.py \
    --samplesheet config/structures.tsv \
    --out rosetta_gaps/work/gaps.tsv

# 2. Generate the job tree. --rosetta-db and the binary names come from step 0.
python3 rosetta_gaps/bin/make_inputs.py \
    --gaps rosetta_gaps/work/gaps.tsv \
    --out-dir rosetta_gaps/work \
    --rosetta-db <database path from step 0> \
    --loopmodel-bin <loopmodel binary from step 0> \
    --remodel-bin   <remodel binary from step 0> \
    --nstruct 500
#    Read rosetta_gaps/work/plan.tsv before submitting anything.

# 3. Stage 1, one task per structure.
N=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="work_dir") c=i; next} {print $c}' \
      rosetta_gaps/work/jobs.tsv | sort -u | wc -l)
sbatch --array=1-$N rosetta_gaps/submit_insert.sbatch \
       rosetta_gaps/work containers/rosetta.sif

# 4. Stage 1 check (the submit script runs it; this is the re-read).
python3 rosetta_gaps/bin/check_models.py \
    --work-dir rosetta_gaps/work/<structure_id> --stage insert

# 5. Stage 2, one task per jobs.tsv row.
sbatch --array=1-$(( $(wc -l < rosetta_gaps/work/jobs.tsv) - 1 )) \
       rosetta_gaps/submit_ngk.sbatch rosetta_gaps/work containers/rosetta.sif

# 6. Validate, rank, and read the spread.
python3 rosetta_gaps/bin/check_models.py \
    --work-dir rosetta_gaps/work/<structure_id> --stage ngk

# 7. Promote ONE model. Lands where structures.tsv can point at it.
python3 rosetta_gaps/bin/promote.py \
    --work-dir rosetta_gaps/work/<structure_id> \
    --model rosetta_gaps/work/<structure_id>/runs/<task>/<model>.pdb.gz \
    --out data/structures/<structure_id>.repaired.cif
```

Then put the printed line into `config/structures.tsv` in place of the row it
replaces, and run the pipeline as usual. **The WT reversion is not applied
here**: `prepare_structures.py` still does it, and that order is deliberate —
reverting before the repair would leave any substitution the rebuild introduced
undetected.

## What to check, and in what order

1. `plan.tsv` — the gaps, their merged loop definitions and their lengths.
   A gap you did not expect means `find_gaps.py` and your reading of the entry
   disagree; settle that before spending cluster time.
2. `insert_report.tsv` — every Remodel candidate, and why any was rejected.
3. `ngk_report.tsv` — every model, valid or not, with its total score.
4. `spread_report.tsv` — **the one to read before promoting.** Pairwise
   C-alpha RMSD of each rebuilt span across the top models, after superposing
   on the residues that were never rebuilt.

On the spread: in AICG2+ the reference structure defines the native contacts
and the equilibrium bond, angle and dihedral values, so whichever loop you
promote *becomes* the definition of the native state and nothing downstream can
penalise it. A gap whose top models agree is forced by the surrounding
structure. A gap whose top models disagree by several Angstrom is a modelling
assumption, and `promote.py` refuses it unless you pass `--force`. If you force
it, it belongs in the limitations of any claim that rests on that loop.

Three checks run on every model: the sequence must match
`expected_sequence.txt`, every consecutive C-N distance must be inside
`--bond-range`, and residues outside every loop must not have moved
(`--fixed-tol`). The third exists because KIC is supposed to hold them fixed,
so a nonzero value means the loop file or the numbering was wrong rather than
that the chain relaxed.

## Decisions baked in, and where to change them

| decision | default | flag |
| - | - | - |
| DNA and ions stay in the pose | keep | `--no-keep-nucleic`, `--no-keep-het` |
| models per structure | 500 | `--nstruct` |
| flank each side of a rebuilt span | 2 | `--flank` |
| longest run to attempt | 12 residues | `--max-len` (on `find_gaps.py`) |
| method | NGK | `--compare-gap <id>` also runs KIC-with-fragments on one loop |

DNA is kept by default because a loop at a DNA interface modelled with the DNA
removed is built in the wrong state — the same reason an apo donor was ruled
out for the KRKR loop in the gap-closing plan. For a structure with no DNA the
flag does nothing.

Two things this does **not** attempt: terminal runs, because KIC holds a loop's
end C-alphas fixed and a terminal run has no anchor on one side; and runs longer
than `--max-len`. `find_gaps.py` lists both with the reason in the `reason`
column rather than dropping them silently.

The dimer is not modelled directly. Repair the monomer, then let
`prepare_structures.py` expand the assembly from the repaired file — which is
why `promote.py` preserves the deposited `_pdbx_struct_assembly*` categories.
A two-copy pose would double every gap, and Remodel's `-run:chain` takes one
chain at a time.

## Tests

```bash
bash rosetta_gaps/tests/run_tests.sh
```

42 assertions over gap enumeration and its refusal paths, job generation and
its numbering, the three model-validity checks, the spread arithmetic,
promotion with metadata preservation, and acceptance of the promoted file by
`bin/prepare_structures.py`. The fixture is a synthetic 20-residue peptide with
planted gaps and an author-numbering offset, so a script that confuses the
numbering systems fails visibly.

Rosetta itself is not exercised: no `remodel`, no `loopmodel`, no container, no
Slurm. The flag names come from the current Rosetta documentation, not from
`-help` on your image — which is what step 0 is for.
