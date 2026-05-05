#!/usr/bin/env python3

"""
Standalone scanner for additional central-pair pT-based selections.

It scans, separately:
 - scaled leading-jet pT:    central_jet1_pt / central_pair_mass > threshold
 - scaled subleading-jet pT: central_jet2_pt / central_pair_mass > threshold
 - pT imbalance:             (central_jet1_pt - central_jet2_pt) / (central_jet1_pt + central_jet2_pt) < threshold

The script reads the reco study tree produced by studyVBFReco_scoutNano.py and
reports how each threshold changes the selected, matched, and unmatched yields.
"""

import os
import re
import csv
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Scan pT-based central-pair cuts on reco VBF study trees")
parser.add_argument(
    "--input",
    type=str,
    required=True,
    help="Input ROOT file containing the Events tree from studyVBFReco_scoutNano.py",
)
parser.add_argument(
    "--tree",
    type=str,
    default="Events",
    help="Tree name inside the input ROOT file",
)
parser.add_argument(
    "--outdir",
    type=str,
    default=None,
    help="Directory for optional CSV outputs. Defaults to scan_ptCuts_recoVBF_M<MASS> when the mass is parsable from the input filename.",
)
parser.add_argument(
    "--baseSelection",
    type=str,
    default="pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1",
    help="Base selection used before applying the scanned cut",
)
parser.add_argument(
    "--leadScaledMin",
    type=float,
    default=0.10,
    help="Minimum threshold in the leading scaled-pT scan",
)
parser.add_argument(
    "--leadScaledMax",
    type=float,
    default=0.60,
    help="Maximum threshold in the leading scaled-pT scan",
)
parser.add_argument(
    "--leadScaledStep",
    type=float,
    default=0.02,
    help="Step size in the leading scaled-pT scan",
)
parser.add_argument(
    "--subleadScaledMin",
    type=float,
    default=0.05,
    help="Minimum threshold in the subleading scaled-pT scan",
)
parser.add_argument(
    "--subleadScaledMax",
    type=float,
    default=0.50,
    help="Maximum threshold in the subleading scaled-pT scan",
)
parser.add_argument(
    "--subleadScaledStep",
    type=float,
    default=0.01,
    help="Step size in the subleading scaled-pT scan",
)
parser.add_argument(
    "--imbalanceMin",
    type=float,
    default=0.10,
    help="Minimum threshold in the pT-imbalance scan",
)
parser.add_argument(
    "--imbalanceMax",
    type=float,
    default=0.90,
    help="Maximum threshold in the pT-imbalance scan",
)
parser.add_argument(
    "--imbalanceStep",
    type=float,
    default=0.02,
    help="Step size in the pT-imbalance scan",
)
parser.add_argument(
    "--writeCSV",
    action="store_true",
    help="Write one CSV file per scan in --outdir",
)
args = parser.parse_args()


def default_outdir():
    if args.outdir:
        return args.outdir

    input_base = os.path.basename(args.input)
    mass_match = re.search(r"_M(\d+)", input_base)
    if mass_match:
        return f"scan_ptCuts_recoVBF_M{mass_match.group(1)}"
    return "scan_ptCuts_recoVBF"


def frange(start, stop, step):
    values = []
    current = start
    while current <= stop + 1e-9:
        values.append(round(current, 6))
        current += step
    return values


def parse_tree(tree, base_selection):
    selector = ROOT.TTreeFormula("base_selection_formula", base_selection, tree)
    events = []

    for entry in tree:
        if selector.EvalInstance() == 0:
            continue

        mass = float(entry.central_pair_mass)
        pt1 = float(entry.central_jet1_pt)
        pt2 = float(entry.central_jet2_pt)
        matched = int(entry.central_pair_genmatched) == 1

        if mass <= 0.0 or pt1 <= 0.0 or pt2 <= 0.0:
            continue

        scaled_pt1 = pt1 / mass
        scaled_pt2 = pt2 / mass
        pt_imbalance = (pt1 - pt2) / (pt1 + pt2) if (pt1 + pt2) > 0.0 else 999.0

        events.append(
            {
                "mass": mass,
                "pt1": pt1,
                "pt2": pt2,
                "scaled_pt1": scaled_pt1,
                "scaled_pt2": scaled_pt2,
                "pt_imbalance": pt_imbalance,
                "matched": matched,
            }
        )

    return events


def evaluate_scan(events, thresholds, variable_name, direction):
    base_total = len(events)
    base_matched = sum(1 for event in events if event["matched"])
    base_unmatched = base_total - base_matched

    rows = []
    for threshold in thresholds:
        if direction == "min":
            selected = [event for event in events if event[variable_name] > threshold]
        elif direction == "max":
            selected = [event for event in events if event[variable_name] < threshold]
        else:
            raise ValueError(f"Unknown scan direction: {direction}")

        n_selected = len(selected)
        n_matched = sum(1 for event in selected if event["matched"])
        n_unmatched = n_selected - n_matched

        matched_eff = 100.0 * n_matched / base_matched if base_matched > 0 else 0.0
        unmatched_eff = 100.0 * n_unmatched / base_unmatched if base_unmatched > 0 else 0.0
        purity = 100.0 * n_matched / n_selected if n_selected > 0 else 0.0
        rej = 100.0 * (base_unmatched - n_unmatched) / base_unmatched if base_unmatched > 0 else 0.0

        rows.append(
            {
                "threshold": threshold,
                "selected": n_selected,
                "matched": n_matched,
                "unmatched": n_unmatched,
                "matched_eff_pct": matched_eff,
                "unmatched_eff_pct": unmatched_eff,
                "matched_purity_pct": purity,
                "unmatched_rejection_pct": rej,
            }
        )

    return {
        "base_total": base_total,
        "base_matched": base_matched,
        "base_unmatched": base_unmatched,
        "rows": rows,
    }


def print_scan(title, summary, comparator_label):
    print(f"\n[SCAN] {title}")
    print(f"  Base selected events : {summary['base_total']}")
    print(f"  Base matched events  : {summary['base_matched']}")
    print(f"  Base unmatched events: {summary['base_unmatched']}")
    print(f"  {'Cut':>11s} {'Selected':>10s} {'Matched':>10s} {'Unmatched':>11s} {'MatchEff%':>10s} {'UnmatchEff%':>12s} {'Purity%':>9s} {'UnmatchRej%':>12s}")
    for row in summary["rows"]:
        cut_label = f"{comparator_label}{row['threshold']:.3f}"
        print(
            f"  {cut_label:>11s}"
            f" {row['selected']:10d}"
            f" {row['matched']:10d}"
            f" {row['unmatched']:11d}"
            f" {row['matched_eff_pct']:10.2f}"
            f" {row['unmatched_eff_pct']:12.2f}"
            f" {row['matched_purity_pct']:9.2f}"
            f" {row['unmatched_rejection_pct']:12.2f}"
        )


def write_csv(path, summary):
    fieldnames = [
        "threshold",
        "selected",
        "matched",
        "unmatched",
        "matched_eff_pct",
        "unmatched_eff_pct",
        "matched_purity_pct",
        "unmatched_rejection_pct",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary["rows"]:
            writer.writerow(row)


def main():
    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open input file: {args.input}")

    tree = root_file.Get(args.tree)
    if tree is None:
        raise RuntimeError(f"Could not find tree '{args.tree}' in {args.input}")

    events = parse_tree(tree, args.baseSelection)
    if not events:
        raise RuntimeError("No events passed the base selection. Check --input or --baseSelection.")

    print(f"[INFO] Input file:           {args.input}")
    print(f"[INFO] Tree:                 {args.tree}")
    print(f"[INFO] Base selection:       {args.baseSelection}")
    print(f"[INFO] Events after base:    {len(events)}")

    lead_summary = evaluate_scan(
        events,
        frange(args.leadScaledMin, args.leadScaledMax, args.leadScaledStep),
        "scaled_pt1",
        "min",
    )
    sublead_summary = evaluate_scan(
        events,
        frange(args.subleadScaledMin, args.subleadScaledMax, args.subleadScaledStep),
        "scaled_pt2",
        "min",
    )
    imbalance_summary = evaluate_scan(
        events,
        frange(args.imbalanceMin, args.imbalanceMax, args.imbalanceStep),
        "pt_imbalance",
        "max",
    )

    print_scan("Leading scaled pT: central_jet1_pt / central_pair_mass", lead_summary, ">")
    print_scan("Subleading scaled pT: central_jet2_pt / central_pair_mass", sublead_summary, ">")
    print_scan("pT imbalance: (pt1 - pt2) / (pt1 + pt2)", imbalance_summary, "<")

    if args.writeCSV:
        outdir = default_outdir()
        os.makedirs(outdir, exist_ok=True)
        write_csv(os.path.join(outdir, "scan_scaledPt_leading.csv"), lead_summary)
        write_csv(os.path.join(outdir, "scan_scaledPt_subleading.csv"), sublead_summary)
        write_csv(os.path.join(outdir, "scan_ptImbalance.csv"), imbalance_summary)
        print(f"\n[INFO] Wrote CSV summaries to: {outdir}")

    root_file.Close()


if __name__ == "__main__":
    main()
