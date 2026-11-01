#!/usr/bin/env python3
"""Export chromosome support from an existing sample run.json without refitting."""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mapploid.reporting import chromosome_support_figures, evidence_figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-json", type=Path, required=True, help="sample details/run.json")
    parser.add_argument("--out", type=Path, required=True, help="new or empty figure directory")
    parser.add_argument("--combined", action="store_true", help="export distributions with the support matrix underneath")
    args = parser.parse_args()
    try:
        record = json.loads(args.run_json.read_text())
        stages = record.get("stages", [])
        if record.get("result", {}).get("status") == "FAILED" or not stages or not stages[-1].get("chromosomes"):
            raise ValueError("No completed sample chromosome results to plot; use a sample run.json, not a batch record")
        if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
            raise ValueError("Use a new or empty output directory")
        if args.combined:
            _, files = evidence_figure(stages[-1], args.out)
            if not files:
                raise ValueError("No fitted distribution data available for the combined figure")
        else:
            _, files = chromosome_support_figures(stages[-1], args.out)
        for path in files:
            print(path.resolve())
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
