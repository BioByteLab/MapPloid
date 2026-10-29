"""Generate a small synthetic 3x/empty-BAM batch to exercise the public outputs."""
from pathlib import Path
import sys
import csv
import json
import hashlib
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mapploid.pipeline import PrepareConfig, prepare
from mapploid.cli import main

def run():
    root = ROOT / "validation" / "demo_schema3"
    root.mkdir(exist_ok=True)
    sys.path.insert(0, str(ROOT / "examples/quickstart"))
    from run import make_inputs
    ref, scores, covered = make_inputs(root)
    assets = root / "assets"
    prepare(PrepareConfig(ref, assets, 50, 1, genmap_bedgraph=scores))
    with pysam.AlignmentFile(str(covered), "rb") as handle:
        header = handle.header.to_dict()
    empty = root / "empty.bam"
    if not empty.exists():
        with pysam.AlignmentFile(str(empty), "wb", header=header):
            pass
        pysam.index(str(empty))
    listing = root / "bams.txt"
    if not listing.exists():
        listing.write_text("synthetic_3x.bam\nempty.bam\nmissing.bam\n")
    output = root / "batch"
    args = ["call", "--bam-list", str(listing), "--assets", str(assets), "--out", str(output), "--export-figures"]
    assert main(args) == 1  # One intentionally missing BAM; other samples continue.
    rows = list(csv.DictReader((output / "result.tsv").open(), delimiter="\t"))
    assert [r["status"] for r in rows] == ["CALLED", "UNCERTAIN", "FAILED"], rows
    assert rows[0]["predicted_ploidy"] == "3x"
    sample = output / "samples/synthetic_3x"
    assert sorted(p.name for p in (sample / "details").iterdir()) == ["run.json", "run.log"]
    record = json.loads((sample / "details/run.json").read_text())
    assert list(record["result"]) == ["sample", "predicted_ploidy", "status", "selected_sites",
                                      "delta_bic", "analyzed_region_mb", "reason", "warnings"]
    assert record["output_schema"] == 3
    profile = record["stages"][-1]["model"]["frequency_profile"]
    assert abs(sum(profile["observed"]) - 1) < 1e-10
    assert set(profile["fitted_models"]) == {"2", "3", "4"}
    for fitted in profile["fitted_models"].values():
        assert len(fitted) == len(profile["observed"]) == len(profile["edges"]) - 1
        assert abs(sum(fitted) - 1) < 1e-10
    fields = ["sample", "predicted_ploidy", "status", "selected_sites", "delta_bic", "analyzed_region_mb"]
    assert list(rows[0]) == fields + ["result_file"]
    for row in rows:
        with (output / row["result_file"]).open() as handle:
            assert next(csv.reader(handle, delimiter="\t")) == fields
    report = (sample / "report.html").read_text()
    svg = report[report.index("<svg"):report.index("</svg>") + len("</svg>")]
    assert all(f"{p}x model" in svg for p in (2, 3, 4))
    assert svg.count("BEST") == 1
    assert "Relative BIC" not in svg
    for fmt in ("pdf", "png"):
        assert (sample / "figures" / f"evidence.{fmt}").is_file()
    assert not (sample / "figures/evidence.svg").exists()
    from PIL import Image
    with Image.open(sample / "figures/evidence.png") as image:
        assert image.format == "PNG" and image.size == (1860, 2160)
        assert abs(image.info["dpi"][0] - 300) < 1
        image.verify()
    tracked = [sample / "result.tsv", sample / "report.html", sample / "details/run.json"]
    stamps = {p: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest()) for p in tracked}
    assert main(args) == 1
    assert stamps == {p: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest()) for p in tracked}
    # Missing-data reports must not contain empty charts.
    assert "<svg" not in (output / "samples/empty/report.html").read_text()
    assert "FAILED" in (output / "samples/missing/report.html").read_text()
    (ROOT / "validation/output_smoke.json").write_text(json.dumps({
        "passed": True, "case": "Synthetic 3x, empty BAM and intentionally missing BAM",
        "batch_report": str(output / "report.html"), "tests": [
            "3x call", "UNCERTAIN no coverage", "FAILED missing input", "batch continuation",
            "six sample columns and seven batch columns", "compact cleanup", "normalized fitting profile",
            "PDF and 300-dpi PNG export without SVG file", "completed reuse without rewriting", "no empty charts"
        ]}, indent=2) + "\n")
    print("Output smoke check passed:", output / "report.html")

if __name__ == "__main__":
    run()
