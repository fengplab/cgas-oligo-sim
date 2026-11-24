#!/usr/bin/env python3
"""Fail if environment.yml and the code's imports have drifted apart.

An environment file is only useful while it matches what the code imports. The
two drift silently: a new import works for whoever added it, because their
shell already has the package, and breaks for everyone who builds the
environment from the file. This check closes that gap.

environment.yml is parsed with a regex rather than PyYAML on purpose - adding
a YAML parser to run this check would itself be a new dependency that the
environment file would then have to carry.
"""
from __future__ import annotations

import glob
import pathlib
import re
import sys

# Modules in the standard library, which need no conda package.
STDLIB = {
    "__future__", "argparse", "collections", "csv", "dataclasses", "glob",
    "hashlib", "html", "io", "itertools", "json", "math", "os", "pathlib",
    "random", "re", "shutil", "subprocess", "sys", "textwrap", "typing",
    "warnings",
}
# Import name -> conda package name, where they differ.
IMPORT_TO_CONDA = {
    "Bio": "biopython",
    "MDAnalysis": "mdanalysis",
    "matplotlib": "matplotlib-base",
    "yaml": "pyyaml",
    "sklearn": "scikit-learn",
    "PIL": "pillow",
}
# Packages allowed to be present without being imported anywhere.
ALLOWED_UNUSED = {"python"}


def conda_deps(path: str) -> set[str]:
    text = pathlib.Path(path).read_text()
    block = re.search(r"^dependencies:\s*$(.*?)(?=^\S|\Z)", text, re.M | re.S)
    if not block:
        sys.exit(f"FAIL: {path} has no dependencies: block")
    out = set()
    for line in block.group(1).splitlines():
        line = line.strip()
        if not line.startswith("- "):
            continue
        spec = line[2:].strip()
        if spec.endswith(":"):          # a nested key such as 'pip:'
            continue
        name = re.split(r"[=<>!~ ]", spec, maxsplit=1)[0].strip()
        if name:
            out.add(name)
    return out


def imported(paths: list[str]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for f in paths:
        for m in re.findall(r"^\s*(?:import|from)\s+([A-Za-z_][A-Za-z0-9_]*)",
                            pathlib.Path(f).read_text(), re.M):
            if m in STDLIB:
                continue
            out.setdefault(IMPORT_TO_CONDA.get(m, m), set()).add(f)
    return out


def main() -> None:
    files = sorted(glob.glob("bin/*.py") + glob.glob("tests/*.py"))
    if not files:
        sys.exit("FAIL: no Python files found; run this from the repo root")
    deps = conda_deps("environment.yml")
    used = imported(files)

    missing = sorted(set(used) - deps)
    unused = sorted(deps - set(used) - ALLOWED_UNUSED)

    if missing:
        for m in missing:
            where = ", ".join(sorted(used[m])[:3])
            print(f"  imported but absent from environment.yml: {m}  ({where})")
    if unused:
        print(f"  in environment.yml but imported nowhere: {', '.join(unused)}")
    if missing or unused:
        sys.exit("FAIL: environment.yml and the code's imports have drifted. "
                 "Anyone building the environment from this file would get a "
                 "different result from whoever last ran the tests.")
    print(f"  OK environment.yml matches {len(used)} third-party import(s) "
          f"across {len(files)} file(s)")


if __name__ == "__main__":
    main()
