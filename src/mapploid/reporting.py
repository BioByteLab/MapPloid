"""Offline, self-contained reports; no changes to ploidy decisions."""
from __future__ import annotations
import csv
import html
import io
import math
from pathlib import Path
import uuid

OUTPUT_SCHEMA = 4
RESULT_FIELDS = ("sample", "predicted_ploidy", "status", "selected_sites", "delta_bic",
                 "analyzed_region_mb")
PRIMARY_WARNINGS = {
    "LOW_BOOTSTRAP_SUPPORT", "INSUFFICIENT_CHROMOSOME_EVIDENCE",
    "CHROMOSOME_HETEROGENEITY",
}
MESSAGES = {
    "LOW_SITE_COUNT": "Too few fitting SNPs. Check coverage, retained-site counts and usable reference regions.",
    "LOW_DELTA_BIC": "The best two models are insufficiently separated. Inspect the frequency distribution and regional consistency.",
    "MODEL_NOT_CONVERGED": "At least one main model did not converge. Inspect model diagnostics before interpreting the candidate.",
    "NO_COVERAGE": "No valid A/C/G/T observations in the analysis regions. Check the BAM, reference and coverage.",
    "NO_VALID_SITES": "Coverage was present, but no sites passed the filters. Inspect the actual depth and allele filters.",
    "NO_ELIGIBLE_REGIONS": "No eligible analysis regions. Check reference preparation and excluded regions.",
    "LOW_BOOTSTRAP_SUPPORT": "The global candidate was not consistently selected across block-bootstrap replicates.",
    "INSUFFICIENT_CHROMOSOME_EVIDENCE": "Fewer than three reference sequences provided informative chromosome-level fits.",
    "CHROMOSOME_HETEROGENEITY": "Fewer than 70% of informative reference sequences supported the global candidate.",
    "TOO_FEW_BLOCKS": "Fewer than five occupied genomic blocks are available for resampling.",
    "RHO_AT_BOUNDARY": "The fitted dispersion parameter is close to an allowed boundary.",
    "A_AT_BOUNDARY": "The tetraploid outer-peak weight is close to its lower bound.",
    "A_AT_UPPER_BOUNDARY": "The tetraploid central-peak weight is close to zero.",
    "MODEL_ASSUMPTION_SENSITIVE": "Relaxing the tetraploid outer-peak constraint changes the best candidate.",
}
STYLE = """
:root{color-scheme:light;--ink:#18313d;--muted:#526873;--line:#dae3e7;--accent:#007c83}
*{box-sizing:border-box}body{margin:0;background:#f4f7f8;color:var(--ink);font:15px/1.6 system-ui,sans-serif}
main{max-width:1160px;margin:32px auto;padding:0 24px 40px}header{margin:24px 0}
h1{font-size:30px;margin:4px 0}h2{font-size:20px;margin:0 0 14px}h3{font-size:16px}
.eyebrow,.muted{color:var(--muted)}.eyebrow{letter-spacing:.08em;font-size:12px}
section,.hero{background:white;border:1px solid var(--line);border-radius:12px;padding:22px;margin:18px 0}
.hero{border-left:5px solid var(--accent)}.result{font-size:34px;font-weight:700}
.metrics{display:flex;flex-wrap:wrap;gap:14px 32px}.metric strong{display:block;font-size:21px}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
th,td{padding:11px 12px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}
th{background:#f2f6f7;font-size:13px}td{overflow-wrap:anywhere}
.scroll{overflow-x:auto}code{font-size:12px;overflow-wrap:anywhere}
.PASS,.CALLED{color:#00664b}.FAIL,.FAILED{color:#b13e12}.WARNING,.UNCERTAIN{color:#916400}
.badge{font-weight:700}details{margin:12px 0}summary{cursor:pointer;font-weight:600}
svg{width:100%;height:auto}a{color:#006e8c}input,select{font:inherit;padding:8px;border:1px solid #aabcc5;border-radius:5px}
.filters{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}.notice{background:#fff7e5;padding:12px 16px;border-radius:6px}
footer{font-size:12px;color:var(--muted)}@media(max-width:650px){main{padding:0 12px}.hero,section{padding:15px}h1{font-size:24px}}
"""
def esc(value):
    return html.escape(str(value), quote=True)

def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False

def num(value, digits=2):
    if not finite(value):
        return "Not calculated"
    if isinstance(value, int):
        return f"{value:,}"
    return f"{float(value):,.{digits}f}"

def flags(value):
    return [v for v in str(value or "").split(";") if v not in ("", "NONE", "NA")]

def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value

def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)

def page(title, body):
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{esc(title)}</title><style>{STYLE}</style></head><body><main>'
            + body + '</main></body></html>')

def table(headers, rows):
    return '<div class="scroll"><table><thead><tr>' + ''.join(f'<th>{esc(x)}</th>' for x in headers) +            '</tr></thead><tbody>' + ''.join('<tr>' + ''.join(f'<td>{esc(v)}</td>' for v in r) + '</tr>' for r in rows) + '</tbody></table></div>'

def diagnostics(stage, settings):
    e = stage.get("evidence", {})
    fits = stage.get("model", {}).get("fits", {})
    computed = all(str(p) in fits for p in (2, 3, 4))
    sites = e.get("selected_sites", 0)
    delta = e.get("delta_bic")
    minimum = settings.get("min_sites", 25000)
    separation = settings.get("delta_bic", 10)
    convergence = [f"{p}x: {'yes' if fits[str(p)].get('converged') else 'no'}" for p in (2, 3, 4)] if computed else []
    gates = [
        ["Fitting SNPs", num(sites), f"At least {minimum:,}", "PASS" if sites >= minimum else "FAIL"],
        ["BIC separation", num(delta), f"At least {separation:g}",
         ("PASS" if float(delta) >= separation else "FAIL") if finite(delta) else "NOT_CALCULATED"],
        ["Main-model convergence", "; ".join(convergence) or "Not calculated", "All three models",
         ("PASS" if all(fits[str(p)].get("converged") for p in (2, 3, 4)) else "FAIL") if computed else "NOT_CALCULATED"],
    ]
    planned = e.get("bootstrap_planned")
    support = e.get("bootstrap_support")
    boot = f"{float(support):.1%}" if planned and finite(support) else "Not calculated"
    boot_context = f"Agreement at least {settings.get('bootstrap_support', 0.85):.0%}; all planned replicates"
    chromosomes = e.get("informative_chromosomes", 0)
    agreement = e.get("chromosome_agreement")
    chrom = f"{float(agreement):.1%}" if chromosomes and finite(agreement) else "Not assessable"
    w = set(flags(e.get("warnings")))
    rows = [
        ["Bootstrap stability", boot, boot_context, "WARNING" if "LOW_BOOTSTRAP_SUPPORT" in w else ("CONSISTENT" if planned and finite(support) else "NOT_CALCULATED")],
        ["Chromosome consistency", f"{chromosomes} informative; {chrom}", "At least 3 informative; agreement at least 70%",
         "WARNING" if w & {"INSUFFICIENT_CHROMOSOME_EVIDENCE", "CHROMOSOME_HETEROGENEITY"} else ("CONSISTENT" if computed and chromosomes and finite(agreement) else "NOT_ASSESSABLE")],
    ]
    return gates, rows

def evidence_figure(stage, export_dir=None):
    """One figure: observed and fitted distributions for 2x, 3x and 4x."""
    model = stage.get("model", {})
    fits = model.get("fits", {})
    profile = model.get("frequency_profile", {})
    curves = profile.get("fitted_models", {})
    edges = profile.get("edges", [])
    if not edges or not curves:
        return "", []
    import matplotlib as mpl
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    import numpy as np

    observed = np.asarray(profile["observed"]) * 100
    valid_bics = {p: float(f["bic"]) for p, f in fits.items() if finite(f.get("bic"))}
    best = min(valid_bics, key=valid_bics.get) if valid_bics else None
    colors = {2: "#597181", 3: "#d79246", 4: "#82475f"}
    files = []
    with mpl.rc_context({"font.size": 10, "svg.fonttype": "none", "svg.hashsalt": "MapPloid",
                         "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"}):
        fig = Figure(figsize=(6.2, 7.2), layout="constrained")
        axes = fig.subplots(3, 1, sharex=True, sharey=True)
        top = fig.suptitle("Observed and model-fitted distributions", fontsize=14, fontweight="bold")
        peak = max([float(observed.max()), *[max(values) * 100 for values in curves.values()]])
        for ax, ploidy in zip(axes, (2, 3, 4)):
            key = str(ploidy)
            color = colors[ploidy]
            ax.stairs(observed, edges, fill=True, color="#c7d1d8", linewidth=0, label="Observed")
            predicted = curves.get(key)
            if predicted is not None:
                ax.stairs(np.asarray(predicted) * 100, edges, color=color, linewidth=1.5, label="Fitted model")
            delta = valid_bics[key] - valid_bics[best] if key in valid_bics and best is not None else None
            title = f"{ploidy}x model   ΔBIC = {delta:,.1f}" if delta is not None else f"{ploidy}x model   ΔBIC = NA"
            if key == best:
                title += "   BEST"
            if not fits.get(key, {}).get("converged"):
                title += "   (not converged)"
            ax.set_title(title, loc="left", color=color, fontsize=11, pad=9)
            ax.set_ylim(0, max(peak * 1.15, 1))
            ax.set_xlim(edges[0], edges[-1])
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].legend(handles=[Patch(facecolor="#c7d1d8", label="Observed"),
                                Line2D([0], [0], color="#597181", linewidth=1.5, label="Fitted model")],
                       loc="lower left", bbox_to_anchor=(0, 1.18), ncol=2, frameon=False,
                       fontsize=9, borderaxespad=0)
        ticks = [(edges[0], f"{edges[0]:.2f}"), (.25, "1/4"), (1/3, "1/3"), (.40, "0.40"), (.5, "1/2")]
        ticks = [(x, label) for x, label in ticks if edges[0] <= x <= edges[-1]]
        # A custom lower bound may coincide with one of the standard tick marks.
        ticks = list(dict(ticks).items())
        axes[-1].set_xticks([x for x, _ in ticks], [label for _, label in ticks])
        axes[-1].set_xlabel("Minor allele fraction")
        fig.supylabel("Sites (%)", fontsize=11)
        svg = io.StringIO()
        fig.savefig(svg, format="svg", metadata={"Date": None})
        image = svg.getvalue()
        image = image[image.index("<svg"):]
        if export_dir:
            export_dir = Path(export_dir)
            export_dir.mkdir(parents=True, exist_ok=True)
            for fmt in ("pdf", "png"):
                p = export_dir / ("evidence." + fmt)
                tmp = p.with_name(p.name + ".tmp")
                fig.savefig(tmp, format=fmt, dpi=300,
                            metadata={"CreationDate": None} if fmt == "pdf" else None)
                tmp.replace(p)
                files.append(p)
    return image, files

def write_sample_report(root, record, export_figures=False):
    root = Path(root)
    row = record["result"]
    stage = record.get("stages", [{}])[-1] if record.get("stages") else {}
    evidence = stage.get("evidence", {})
    state = row["status"]
    reasons, warnings = flags(row.get("reason")), flags(row.get("warnings"))
    settings = record.get("model_settings", {})
    gates, diagnostic_rows = diagnostics(stage, settings)
    body = f'<header><div class="eyebrow">MAPPLOID {esc(record["version"])} · SAMPLE REPORT</div><h1>{esc(row["sample"])}</h1></header>'
    body += f'<div class="hero"><div class="result">{esc(row["predicted_ploidy"])}</div><span class="badge {esc(state)}">{esc(state)}</span>'
    if state == "UNCERTAIN" and evidence.get("best_hypothesis"):
        body += f'<p>Best candidate: {esc(evidence["best_hypothesis"])}. No definitive ploidy call.</p>'
    if record.get("state") == "INTERRUPTED":
        body += '<p>Run interrupted. Rerun the same command to resume.</p>'
    if state in ("UNCERTAIN", "FAILED"):
        body += '<p>See <a href="details/run.json">details/run.json</a> for the reason.</p>'
    body += '<div class="metrics">'
    for label, val in [("Analyzed region (Mb)", row.get("analyzed_region_mb")),
                       ("Available qualifying SNPs", evidence.get("available_sites")),
                       ("Fitting SNPs", row.get("selected_sites"))]:
        body += f'<div class="metric">{esc(label)}<strong>{esc(num(val))}</strong></div>'
    body += '</div></div>'
    if state != "FAILED":
        body += '<section><h2>Decision requirements and supporting diagnostics</h2>'
        body += '<p class="muted">Only the three decision requirements determine the ploidy call. The two supporting diagnostics do not veto it. Bootstrap support is resampling agreement, not the probability that the call is correct.</p>'
        combined = [["Decision requirement", *row] for row in gates]
        combined += [["Supporting diagnostic", *row] for row in diagnostic_rows]
        body += table(["Role", "Item", "Observed", "Criterion or context", "Result"], combined) + '</section>'
        svg, files = evidence_figure(stage, root / "figures" if export_figures else None)
        if svg:
            body += '<section><h2>Observed and model-fitted distributions</h2>' + svg
            body += '<p class="muted">Frequency = minor count / sum of the two leading allele counts. All three rows use the same fitting SNPs and bins. Each model is shown with its own fitted parameters, the observed depth distribution and the same discrete truncation as the likelihood. No smoothing. ΔBIC is relative to the best model; BEST identifies the best candidate, not necessarily a definitive call.</p></section>'
    else:
        files = []
    technical = reasons + warnings
    body += '<section><details><summary>Detailed diagnostics</summary>'
    body += ''.join(f'<p><code>{esc(c)}</code>: {esc(MESSAGES.get(c, c))}</p>' for c in technical) or '<p>No additional warnings.</p>'
    chromosomes = stage.get("chromosomes", [])
    if chromosomes:
        body += '<h3>Reference-sequence fits</h3>' + table(["Sequence", "SNPs", "Candidate", "Delta BIC", "Status"],
            [[c.get("chrom"), c.get("sites"), c.get("best_hypothesis"), c.get("delta_bic"), c.get("status")] for c in chromosomes])
    body += '</details><details><summary>Inputs and analysis settings</summary>'
    request = record.get("request", {})
    reference = record.get("reference", {})
    body += table(["Setting", "Value"], [
        ["BAM", request.get("bam", "NA")], ["Reference assets", request.get("assets", "NA")],
        ["K-mer length", reference.get("kmer_length", "NA")],
        ["Minimum region length", reference.get("min_region_length", "NA")],
        ["Eligible reference region (Mb)", num(reference.get("eligible_bases", 0) / 1e6)],
        ["Seed", settings.get("seed", "NA")],
    ] + [[k, v] for k, v in request.get("filters", {}).items()])
    actual = stage.get("counts", {}).get("thresholds", {})
    if actual:
        body += '<h3>Actual site filters</h3>' + table(["Setting", "Value"], list(actual.items()))
    if record.get("input_warnings"):
        body += '<h3>Input notices</h3>' + ''.join(f'<p>{esc(x)}</p>' for x in record["input_warnings"])
    body += '</details></section>'
    body += f'<footer>Output schema {OUTPUT_SCHEMA} · <a href="result.tsv">Result table</a> · <a href="details/run.json">Run record</a> · <a href="details/run.log">Run log</a>'
    if (root / "regions.bed").exists():
        body += ' · <a href="regions.bed">Analysis regions</a>'
    body += '</footer>'
    p = root / "report.html"
    atomic_text(p, page(f'MapPloid — {row["sample"]}', body))
    return [p, *files]

def write_batch_report(root, rows, total, state, version):
    from collections import Counter
    root = Path(root)
    counts = Counter(r["status"] for r in rows)
    ploidies = Counter(r["predicted_ploidy"] for r in rows if r["status"] == "CALLED")
    body = f'<header><div class="eyebrow">MAPPLOID {esc(version)} · BATCH REPORT</div><h1>Sample overview</h1></header>'
    body += f'<div class="hero"><strong>{esc(state)}</strong><p>{len(rows)} / {total} samples processed · {counts["CALLED"]} called · {counts["UNCERTAIN"]} uncertain · {counts["FAILED"]} failed</p>'
    body += '<p>Called ploidies: ' + ' · '.join(f'{p}: {ploidies[p]}' for p in ("2x", "3x", "4x")) + '</p></div>'
    body += '<section><h2>Results</h2><div class="filters"><input id="search" aria-label="Search sample" placeholder="Search sample"><select id="status" aria-label="Filter status"><option value="">All statuses</option><option>CALLED</option><option>UNCERTAIN</option><option>FAILED</option></select></div>'
    body += '<div class="scroll"><table id="samples"><thead><tr>' + ''.join(f'<th>{esc(x)}</th>' for x in ("Sample", "Ploidy", "Status", "Fitting SNPs", "Delta BIC")) + '</tr></thead><tbody>'
    for r in rows:
        name = r["sample"]
        path = root / "samples" / name / "report.html"
        label = f'<a href="samples/{esc(name)}/report.html">{esc(name)}</a>' if path.is_file() else esc(name)
        body += f'<tr data-sample="{esc(name.lower())}" data-status="{esc(r["status"])}"><td>{label}</td>'
        for x in [r["predicted_ploidy"], r["status"], num(r["selected_sites"], 0), num(r["delta_bic"])]:
            body += f'<td>{esc(x)}</td>'
        body += '</tr>'
    body += '</tbody></table></div></section><footer><a href="result.tsv">Download result table</a> · <a href="details/run.json">Run record</a>. Ploidy distributions are predictions, not accuracy estimates.</footer>'
    body += """<script>
const q=document.getElementById('search'),s=document.getElementById('status');
function filter(){for(const r of document.querySelectorAll('#samples tbody tr')){
r.hidden=!(r.dataset.sample.includes(q.value.toLowerCase())&&(!s.value||r.dataset.status===s.value));}}
q.addEventListener('input',filter);s.addEventListener('change',filter);
</script>"""
    atomic_text(root / "report.html", page("MapPloid — batch results", body))
