#!/usr/bin/env python3
"""Fixture-only helper: rekey synthetic pose/contact files onto called site ids.

Matched by position, never by sort order. The real pipeline does not need this:
pose extraction runs downstream of CALL_SITES and is keyed on site_id directly.
"""
import glob, os, re, shutil, sys
import pandas as pd

sites = pd.read_csv('tests/out/sites.tsv', sep='\t')
if len(sites) == 0:
    sys.exit("no called sites in tests/out/sites.tsv; nothing to rekey")
called = {int(r.peak_position): r.site_id for r in sites.itertuples()}
ct = pd.read_csv('tests/fixtures/contacts.tsv', sep='\t')

rows, mapping = [], {}
for f in sorted(glob.glob('tests/fixtures/*_poses.npz')):
    old = os.path.basename(f).replace('_poses.npz', '')
    pos = int(re.search(r'_s(\d+)$', old).group(1))
    nearest = min(called, key=lambda p: abs(p - pos))
    if abs(nearest - pos) > 5:
        sys.exit(f"no called site within 5 bp of fixture {old}")
    new = called[nearest]
    mapping[old] = new
    dst = os.path.abspath(f'tests/fixtures/{new}_poses.npz')
    if os.path.abspath(f) != dst:
        shutil.move(f, dst)
    seq = sites.loc[sites.site_id == new, 'sequence_id'].iloc[0]
    rows.append({'sequence_id': seq, 'site_id': new, 'npz_path': dst})

ct['site_id'] = ct['site_id'].replace(mapping)
ct.to_csv('tests/fixtures/contacts.tsv', sep='\t', index=False)
pd.DataFrame(rows).to_csv('tests/fixtures/poses.tsv', sep='\t', index=False)
print(f"[rekey] {mapping}")
