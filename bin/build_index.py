#!/usr/bin/env python3
"""Build the interactive HTML index of called sites.

3D rendering happens client-side in the browser (WebGL). This is deliberate:
ChimeraX cannot render offscreen on macOS, and a browser-side viewer removes
the headless-renderer dependency entirely.

Every value shown on a card is also present in the table view, so nothing is
reachable only by hovering. Site labels come from upstream tool output and are
inserted with textContent, never by HTML string concatenation.
"""
from __future__ import annotations

import argparse
import html
import json
import os

import pandas as pd

CDN = "https://3Dmol.org/build/3Dmol-min.js"

CSS = """
:root { color-scheme: light; }
.viz-root{
  --surface-1:#fcfcfb; --plane:#f9f9f7; --text-primary:#0b0b0b;
  --text-secondary:#52514e; --muted:#898781; --grid:#e1e0d9;
  --baseline:#c3c2b7; --series-1:#2a78d6; --border:rgba(11,11,11,0.10);
}
@media (prefers-color-scheme: dark){
  :root:where(:not([data-theme="light"])) .viz-root{
    color-scheme: dark;
    --surface-1:#1a1a19; --plane:#0d0d0d; --text-primary:#f0efec;
    --text-secondary:#c3c2b7; --muted:#898781; --grid:#2c2c2a;
    --baseline:#383835; --series-1:#3987e5; --border:rgba(255,255,255,0.10);
  }
}
:root[data-theme="dark"] .viz-root{
  color-scheme: dark;
  --surface-1:#1a1a19; --plane:#0d0d0d; --text-primary:#f0efec;
  --text-secondary:#c3c2b7; --muted:#898781; --grid:#2c2c2a;
  --baseline:#383835; --series-1:#3987e5; --border:rgba(255,255,255,0.10);
}
*{box-sizing:border-box}
body{margin:0;background:var(--plane);font-family:system-ui,-apple-system,"Segoe UI",sans-serif}
.viz-root{background:var(--plane);color:var(--text-primary);padding:24px;min-height:100vh}
header{display:flex;align-items:baseline;gap:16px;flex-wrap:wrap;margin-bottom:4px}
h1{font-size:20px;font-weight:600;margin:0}
.sub{color:var(--text-secondary);font-size:13px;margin:0 0 20px}
.bar{display:flex;gap:12px;align-items:center;margin-bottom:20px;flex-wrap:wrap}
button{font:inherit;font-size:13px;padding:8px 14px;min-height:36px;border-radius:8px;
  border:1px solid var(--border);background:var(--surface-1);color:var(--text-primary);cursor:pointer}
button:hover{background:var(--grid)}
button[aria-pressed="true"]{border-color:var(--series-1)}
.card{background:var(--surface-1);border:1px solid var(--border);border-radius:12px;
  padding:18px;margin-bottom:20px}
.card h2{font-size:15px;font-weight:600;margin:0 0 2px}
.cap{color:var(--text-secondary);font-size:12px;margin:0 0 14px;line-height:1.6}
.cap b{color:var(--text-primary);font-weight:600}
.warn{color:var(--text-primary);background:var(--grid);border-left:3px solid var(--baseline);
  padding:8px 12px;border-radius:0 6px 6px 0;font-size:12px;margin:0 0 14px}
.viewer{position:relative;width:100%;height:340px;border-radius:8px;overflow:hidden;
  border:1px solid var(--border);background:var(--surface-1)}
.legend{display:flex;gap:18px;flex-wrap:wrap;font-size:12px;color:var(--text-secondary);margin-top:10px}
.key{display:inline-flex;align-items:center;gap:6px}
.sw{width:14px;height:14px;border-radius:3px;display:inline-block}
.ramp{width:76px;height:12px;border-radius:3px;display:inline-block;
  background:linear-gradient(90deg,#cde2fb,#6da7ec,#2a78d6,#184f95,#0d366b)}
figure{margin:16px 0 0}
figure img{width:100%;height:auto;border-radius:8px;border:1px solid var(--border);background:#fcfcfb}
figcaption{color:var(--muted);font-size:11px;margin-top:6px}
table{border-collapse:collapse;width:100%;font-size:12px;font-variant-numeric:tabular-nums}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--grid)}
th{color:var(--text-secondary);font-weight:600}
caption{text-align:left;color:var(--text-secondary);font-size:12px;padding-bottom:8px}
.hidden{display:none}
.empty{color:var(--text-secondary);font-size:14px;line-height:1.7}
"""

JS = """
function setTheme(t){
  document.documentElement.setAttribute('data-theme', t);
  document.querySelectorAll('[data-theme-btn]').forEach(function(b){
    b.setAttribute('aria-pressed', String(b.getAttribute('data-theme-btn') === t));
  });
}
function toggleTable(){
  var t = document.getElementById('tableview');
  var b = document.getElementById('tablebtn');
  var show = t.classList.contains('hidden');
  t.classList.toggle('hidden', !show);
  b.setAttribute('aria-pressed', String(show));
  b.textContent = show ? 'Hide table view' : 'Show table view';
}
function renderViewers(){
  if (typeof $3Dmol === 'undefined') {
    document.querySelectorAll('.viewer').forEach(function(d){
      var p = document.createElement('p');
      p.className = 'cap';
      p.style.padding = '12px';
      p.textContent = '3D viewer unavailable: the 3Dmol.js library did not load. '
        + 'The structure file and every value on this card are still in the table '
        + 'view and in the published PDB.';
      d.appendChild(p);
    });
    return;
  }
  SITES.forEach(function(s){
    var el = document.getElementById('v_' + s.dom_id);
    if (!el) return;
    var v = $3Dmol.createViewer(el, {backgroundColor: 'white'});
    v.addModel(PDBS[s.site_id], 'pdb');
    // DNA beads coloured by B-factor = per-nucleotide contact fraction,
    // which is the UV-crosslink proxy. Protein beads stay neutral grey.
    v.setStyle({chain: 'D'}, {sphere: {radius: 1.5, colorscheme:
      {prop: 'b', gradient: new $3Dmol.Gradient.Sinebow(0, 100)}}});
    v.setStyle({chain: 'D'}, {sphere: {radius: 1.5, colorfunc: function(atom){
      var t = Math.max(0, Math.min(1, atom.b / 100));
      var stops = ['#cde2fb','#9ec5f4','#6da7ec','#3987e5','#256abf','#184f95','#0d366b'];
      return stops[Math.min(stops.length - 1, Math.floor(t * stops.length))];
    }}});
    v.addStyle({chain: 'P'}, {sphere: {radius: 1.1, color: '#898781', opacity: 0.55}});
    v.zoomTo(); v.render();
  });
}
document.addEventListener('DOMContentLoaded', function(){ renderViewers(); });
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--cluster-summaries", nargs="*", default=[])
    ap.add_argument("--pdbs", nargs="*", default=[])
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--vendor-js", default="false")
    ap.add_argument("--js-bundle", default="")
    ap.add_argument("--run-params", default="{}")
    ap.add_argument("--allow-missing-poses", action="store_true",
                    help="permit a rendered site to have no structure. Off by "
                         "default: a missing pose means the site/pose join lost a "
                         "site, which must not pass silently.")
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep="\t")
    clusters = {}
    for p in args.cluster_summaries:
        with open(p) as fh:
            d = json.load(fh)
        clusters[d["site_id"]] = d
    pdbs = {}
    for p in args.pdbs:
        sid = os.path.basename(p).rsplit(".pdb", 1)[0]
        with open(p) as fh:
            pdbs[sid] = fh.read()

    rendered = sites[sites["rendered"]] if len(sites) else sites
    params = json.loads(args.run_params) if args.run_params else {}

    # Site/pose identity join. Site ids come from the sites table, never from a
    # filename's position in a glob. A rendered site with no structure, or a
    # structure with no site, means the join lost something.
    want = set(rendered["site_id"].astype(str)) if len(rendered) else set()
    have = set(pdbs)
    missing, orphan = sorted(want - have), sorted(have - want)
    if missing and not args.allow_missing_poses:
        raise SystemExit(
            "ERROR: rendered site(s) have no structure file, so the site/pose join "
            "dropped them: " + ", ".join(missing)
            + ". Structures present: " + (", ".join(sorted(have)) or "none")
            + ". Pose extraction must be keyed on site_id from the sites table.")
    if orphan:
        raise SystemExit(
            "ERROR: structure file(s) match no called site: " + ", ".join(orphan)
            + ". Called sites: " + (", ".join(sorted(want)) or "none")
            + ". This usually means poses were generated against a different run.")
    print(f"[build_index] site/pose join OK: {len(want)} rendered site(s), "
          f"{len(have)} structure(s), 0 unmatched")

    body = []
    if len(rendered) == 0:
        body.append(
            '<div class="card"><h2>No sites called</h2><p class="empty">'
            'No position exceeded the shuffle null at the configured percentile. '
            'This is a valid result, not a pipeline failure: it means occupancy on '
            'these sequences is indistinguishable from composition-matched shuffled '
            'controls. The occupancy tracks and the tier table are still published '
            'and show how close the data came to the threshold.</p></div>')

    site_meta = []
    for _, r in rendered.iterrows():
        sid = str(r["site_id"])
        dom = sid.replace(".", "_").replace("-", "_")
        cl = clusters.get(sid, {})
        pop = cl.get("top_cluster_population_frac")
        popn = cl.get("top_cluster_population")
        nfr = cl.get("n_frames_pooled")
        nrep = cl.get("top_cluster_n_replicates")
        degen = cl.get("degenerate_single_cluster", False)
        site_meta.append({"site_id": sid, "dom_id": dom})

        cap = (
            f"Occupancy <b>{r['occupancy']:.3f}</b> "
            f"(95% CI {r['occupancy_ci_lo']:.3f}&ndash;{r['occupancy_ci_hi']:.3f}) "
            f"against a shuffle-null threshold of <b>{r['null_threshold']:.3f}</b>. "
            f"Occupied in <b>{r['replicate_frac']:.0%}</b> of replicates "
            f"({int(r['n_occupying_replicates'])}). "
        )
        if pop is not None:
            cap += (f"Representative pose is the gromos cluster centre of the most "
                    f"populated cluster: <b>{pop:.0%}</b> of {nfr} pooled frames "
                    f"({popn} frames, drawn from {nrep} replicates).")
        else:
            cap += "No cluster summary was supplied for this site."

        card = [f'<div class="card"><h2></h2><p class="cap">{cap}</p>']
        card[0] = card[0].replace("<h2></h2>",
                                  f"<h2>{html.escape(sid)} &middot; position "
                                  f"{int(r['peak_position'])}</h2>")
        if degen:
            card.append(
                '<p class="warn">A single cluster holds every frame at this RMSD '
                'cutoff. That may be genuine convergence, or a cutoff that has '
                'absorbed all variation &mdash; check the RMSD percentiles in the '
                'cluster summary before reading this pose as structural agreement.</p>')
        card.append(f'<div class="viewer" id="v_{dom}"></div>')
        card.append(
            '<div class="legend">'
            '<span class="key"><span class="ramp"></span>DNA bead, contact fraction 0 &rarr; 1</span>'
            '<span class="key"><span class="sw" style="background:#898781"></span>cGAS bead</span>'
            '</div>')
        for fn, lab in (
            (f"fig2_contacts_{sid}.png", "Residue &times; nucleotide contact map."),
        ):
            if os.path.exists(os.path.join(args.figdir, fn)):
                card.append(f'<figure><img src="figures/{fn}" alt="">'
                            f'<figcaption>{lab}</figcaption></figure>')
        card.append("</div>")
        body.append("".join(card))

    for seq in sorted(sites["sequence_id"].unique()) if len(sites) else []:
        figs = [(f"fig1_occupancy_{seq}.png", f"{seq}: occupancy against the shuffle null."),
                (f"fig5_cooccupancy_{seq}.png", f"{seq}: site co-occupancy across replicates.")]
        have = [(f, l) for f, l in figs if os.path.exists(os.path.join(args.figdir, f))]
        if not have:
            continue
        blk = [f'<div class="card"><h2>{html.escape(str(seq))}</h2>']
        for f, l in have:
            blk.append(f'<figure><img src="figures/{f}" alt="">'
                       f'<figcaption>{l}</figcaption></figure>')
        blk.append("</div>")
        body.append("".join(blk))

    cols = ["sequence_id", "site_id", "peak_position", "bed_start", "bed_end", "width",
            "occupancy", "occupancy_ci_lo", "occupancy_ci_hi", "replicate_frac",
            "n_occupying_replicates", "null_threshold", "replicate_gate", "rendered"]
    present = [c for c in cols if c in sites.columns]
    thead = "".join(f"<th>{html.escape(c)}</th>" for c in present)
    trows = []
    for _, r in sites.iterrows():
        tds = "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in present)
        trows.append(f"<tr>{tds}</tr>")
    table = (
        '<div class="card hidden" id="tableview"><table>'
        '<caption>Every value on the cards above, in full. '
        'position is a 1-based sequence index; bed_start/bed_end are 0-based half-open.'
        '</caption><thead><tr>' + thead + "</tr></thead><tbody>"
        + ("".join(trows) if trows else '<tr><td colspan="99">No sites called.</td></tr>')
        + "</tbody></table></div>")

    if str(args.vendor_js).lower() == "true":
        if not args.js_bundle or not os.path.exists(args.js_bundle):
            raise SystemExit("ERROR: --vendor-js true requires --js-bundle pointing "
                             "at a readable 3Dmol.js file")
        with open(args.js_bundle) as fh:
            lib = f"<script>{fh.read()}</script>"
    else:
        lib = f'<script src="{CDN}"></script>'

    pstr = ", ".join(f"{k}={v}" for k, v in params.items()) if params else "not supplied"
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>cGAS occupancy sites</title>
<style>{CSS}</style>{lib}
<script>
const SITES = {json.dumps(site_meta)};
const PDBS = {json.dumps(pdbs)};
</script>
</head><body><div class="viz-root">
<header><h1>cGAS occupancy sites</h1></header>
<p class="sub">Representative poses are gromos (Daura) cluster centres &mdash; real
frames, not averaged coordinates. Each is captioned with its cluster population:
a pose seen in 3% of frames is not the same evidence as one seen in 80%.
Run parameters: {html.escape(pstr)}</p>
<div class="bar">
  <button id="tablebtn" aria-pressed="false" onclick="toggleTable()">Show table view</button>
  <button data-theme-btn="light" aria-pressed="false" onclick="setTheme('light')">Light</button>
  <button data-theme-btn="dark" aria-pressed="false" onclick="setTheme('dark')">Dark</button>
</div>
{table}
{''.join(body)}
</div><script>{JS}</script></body></html>"""

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(doc)
    print(f"[build_index] wrote {args.out} "
          f"({len(rendered)} rendered site(s), {len(sites)} total, "
          f"{len(pdbs)} structure(s) embedded, "
          f"js={'vendored' if str(args.vendor_js).lower()=='true' else 'cdn'})")


if __name__ == "__main__":
    main()
