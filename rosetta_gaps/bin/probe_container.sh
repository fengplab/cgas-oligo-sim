#!/bin/bash
# Report what is actually inside the Rosetta image, before any flag file is
# written against it. Read-only: this script never writes into the container
# and never launches a modelling run.
#
# Run on the CLUSTER, from the project directory on BeeGFS:
#     module load singularity
#     bash rosetta_gaps/bin/probe_container.sh containers/rosetta.sif
#
# Everything the generator needs is unverifiable from outside the image:
# the binary suffix (linuxgccrelease / default.linuxgccrelease / static),
# the release, the database path, and whether the blueprint helper ships.
# Paste the output back before trusting any generated flag file.
set -euo pipefail

SIF="${1:-containers/rosetta.sif}"
[ -f "$SIF" ] || { echo "no such image: $SIF" >&2; exit 1; }

run() { singularity exec "$SIF" "$@"; }

echo "=== image"
ls -l "$SIF"
singularity inspect "$SIF" 2>&1 | sed 's/^/    /' || true

echo
echo "=== loopmodel / remodel binaries on PATH and in common prefixes"
# `command -v` only finds what is on PATH; Rosetta images often leave the
# binaries outside it, so look in the usual prefixes too. -maxdepth keeps this
# from walking the whole filesystem.
for app in loopmodel remodel relax score_jd2; do
    printf '%-12s PATH: ' "$app"
    run bash -lc "command -v $app || command -v $app.linuxgccrelease || command -v $app.default.linuxgccrelease || echo '(not on PATH)'"
done
echo "--- find, depth-limited"
run bash -lc 'for d in /opt /usr/local /rosetta /software /app; do
    [ -d "$d" ] || continue
    find "$d" -maxdepth 6 -type f \( -name "loopmodel*" -o -name "remodel*" \) 2>/dev/null
done | head -40'

echo
echo "=== version / release"
# -version is the documented way to get the release string. If the binary name
# guessed above is wrong this prints nothing, which is itself the answer.
run bash -lc 'for b in loopmodel loopmodel.linuxgccrelease loopmodel.default.linuxgccrelease \
                       loopmodel.static.linuxgccrelease; do
    command -v "$b" >/dev/null 2>&1 || continue
    echo "binary: $b"
    "$b" -version 2>&1 | head -20
    break
done'

echo
echo "=== database"
# ROSETTA3 / ROSETTA_DATABASE are conventions, not guarantees. A database
# directory is identifiable by its chemical/ subtree.
run bash -lc 'echo "ROSETTA3=${ROSETTA3:-unset}"
              echo "ROSETTA_DATABASE=${ROSETTA_DATABASE:-unset}"
              for d in /opt /usr/local /rosetta /software /app; do
                  [ -d "$d" ] || continue
                  find "$d" -maxdepth 6 -type d -name database 2>/dev/null
              done | head -20'
echo "--- does a candidate database hold chemical/ ?"
run bash -lc 'for d in $(for p in /opt /usr/local /rosetta /software /app; do
                  [ -d "$p" ] && find "$p" -maxdepth 6 -type d -name database 2>/dev/null
              done | head -20); do
                  [ -d "$d/chemical" ] && echo "YES $d"
              done'

echo
echo "=== blueprint helper (needed to seed blueprints from coordinates)"
# The docs spell this both getBlueprintFromCoords.pl and
# getBluePrintFromCoords.pl. Match case-insensitively rather than guessing.
run bash -lc 'for p in /opt /usr/local /rosetta /software /app; do
                  [ -d "$p" ] || continue
                  find "$p" -maxdepth 8 -iname "getblueprintfromcoords.pl" 2>/dev/null
              done | head
              command -v perl >/dev/null 2>&1 && echo "perl: $(perl --version | sed -n 3p)" || echo "perl: ABSENT (the helper cannot run)"'

echo
echo "=== do the NGK options exist in this build?"
# -help lists the options this binary was compiled with. If a flag the
# generator writes is absent here, the run will die on an unrecognised option
# rather than silently ignoring it - but find out now, not in an array job.
run bash -lc 'b=$(command -v loopmodel || command -v loopmodel.linuxgccrelease \
                  || command -v loopmodel.default.linuxgccrelease || true)
    [ -n "$b" ] || { echo "(no loopmodel binary found; skipping)"; exit 0; }
    "$b" -help 2>&1 | grep -E -- "-loops:(remodel|refine|loop_file|kic_rama2b|kic_omega_sampling|ramp_fa_rep|ramp_rama|taboo_sampling|frag_files|frag_sizes)|-allow_omega_move" \
        || echo "(no matching options in -help output - check the version above)"'

echo
echo "=== threading / MPI"
# Rosetta is single-threaded per process unless built with extras. If this is
# an MPI build the array job should still run one rank per task.
run bash -lc 'b=$(command -v loopmodel || command -v loopmodel.linuxgccrelease \
                  || command -v loopmodel.default.linuxgccrelease || true)
    [ -n "$b" ] && ldd "$b" 2>/dev/null | grep -iE "mpi|omp" || echo "(no mpi/omp libraries linked)"'

echo
echo "=== done. Paste everything above back before any flag file is trusted."
