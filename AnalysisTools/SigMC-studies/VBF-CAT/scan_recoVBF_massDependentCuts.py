#!/usr/bin/env python3

"""
Standalone comparator for mass-dependent central-pair selections.

It evaluates a mass-dependent working point based on:
 - central_pair_dR > dR_min(mass)
 - pT imbalance = (pt1 - pt2) / (pt1 + pt2) < A_pt_max(mass)

and compares offline scenarios on the stored central pair.

Recommended usage:
 - run on trees produced with loose central-pair cuts, e.g.
   --centralPairDEtaMax 999 --centralPairDPhiMin -999
 - then use this script to compare offline angular and dR/A_pt scenarios
   on the chosen pair.
"""

import os
import re
import csv
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)


DEFAULT_WORKING_POINTS = {
    300: {"dr_min": 1.5, "apt_max": 0.60},
    500: {"dr_min": 2.0, "apt_max": 0.45},
    1000: {"dr_min": 2.5, "apt_max": 0.35},
}


parser = argparse.ArgumentParser(description="Compare offline mass-dependent dR + pT-imbalance selections on reco VBF trees")
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
    "--mass",
    type=int,
    default=None,
    help="Mass hypothesis used to choose default working points. If omitted, it is parsed from the input filename.",
)
parser.add_argument(
    "--drMin",
    type=float,
    default=None,
    help="Override the mass-dependent minimum central-pair dR.",
)
parser.add_argument(
    "--aptMax",
    type=float,
    default=None,
    help="Override the mass-dependent maximum pT imbalance.",
)
parser.add_argument(
    "--currentDEtaMax",
    type=float,
    default=1.3,
    help="Reference maximum |DeltaEta| used for the current-style angular requirement.",
)
parser.add_argument(
    "--currentDPhiMin",
    type=float,
    default=1.5,
    help="Reference minimum |DeltaPhi| used for the current-style angular requirement.",
)
parser.add_argument(
    "--baseSelection",
    type=str,
    default="pass_trigger_baseline == 1 && has_vbf_tag == 1 && central_pair_mass > 0 && central_jet1_pt > 0 && central_jet2_pt > 0",
    help="Base selection before applying the tested offline central-pair requirements.",
)
parser.add_argument(
    "--outdir",
    type=str,
    default=None,
    help="Directory for optional CSV output. Defaults to scan_massDependentCuts_M<MASS>.",
)
parser.add_argument(
    "--writeCSV",
    action="store_true",
    help="Write the comparison table to CSV.",
)
args = parser.parse_args()


def infer_mass():
    if args.mass is not None:
        return args.mass
    match = re.search(r"_M(\d+)", os.path.basename(args.input))
    if match:
        return int(match.group(1))
    raise RuntimeError("Could not infer mass from input filename. Please pass --mass explicitly.")


def get_working_point(mass):
    if args.drMin is not None and args.aptMax is not None:
        return args.drMin, args.aptMax

    if mass not in DEFAULT_WORKING_POINTS:
        raise RuntimeError(
            f"No default working point configured for mass {mass}. "
            "Pass both --drMin and --aptMax explicitly."
        )

    wp = DEFAULT_WORKING_POINTS[mass]
    dr_min = args.drMin if args.drMin is not None else wp["dr_min"]
    apt_max = args.aptMax if args.aptMax is not None else wp["apt_max"]
    return dr_min, apt_max


def default_outdir(mass):
    if args.outdir:
        return args.outdir
    return f"scan_massDependentCuts_M{mass}"


def parse_events(tree, base_selection):
    selector = ROOT.TTreeFormula("base_selection_formula_massdep", base_selection, tree)
    events = []

    for entry in tree:
        if selector.EvalInstance() == 0:
            continue

        pt1 = float(entry.central_jet1_pt)
        pt2 = float(entry.central_jet2_pt)
        if pt1 <= 0.0 or pt2 <= 0.0:
            continue

        apt = (pt1 - pt2) / (pt1 + pt2) if (pt1 + pt2) > 0.0 else 999.0
        events.append(
            {
                "matched": int(entry.central_pair_genmatched) == 1,
                "deta": abs(float(entry.central_pair_dEta)),
                "dphi": abs(float(entry.central_pair_dPhi)),
                "dr": float(entry.central_pair_dR),
                "apt": apt,
            }
        )

    return events


def count_summary(events, predicate):
    selected = [event for event in events if predicate(event)]
    n_selected = len(selected)
    n_matched = sum(1 for event in selected if event["matched"])
    n_unmatched = n_selected - n_matched
    return n_selected, n_matched, n_unmatched


def make_row(name, events, predicate, base_matched, base_unmatched):
    n_selected, n_matched, n_unmatched = count_summary(events, predicate)
    return {
        "selection": name,
        "selected": n_selected,
        "matched": n_matched,
        "unmatched": n_unmatched,
        "matched_eff_pct": 100.0 * n_matched / base_matched if base_matched > 0 else 0.0,
        "unmatched_eff_pct": 100.0 * n_unmatched / base_unmatched if base_unmatched > 0 else 0.0,
        "matched_purity_pct": 100.0 * n_matched / n_selected if n_selected > 0 else 0.0,
        "unmatched_rejection_pct": 100.0 * (base_unmatched - n_unmatched) / base_unmatched if base_unmatched > 0 else 0.0,
    }


def print_table(rows, mass, dr_min, apt_max):
    print(f"[INFO] Mass hypothesis:        {mass}")
    print(f"[INFO] Working point dR min:  {dr_min}")
    print(f"[INFO] Working point A_pt max:{apt_max}")
    print(f"[INFO] Reference |dEta| max:  {args.currentDEtaMax}")
    print(f"[INFO] Reference |dPhi| min:  {args.currentDPhiMin}")
    print("[INFO] Interpretation: this is most meaningful on trees produced with loose pair-finding cuts")
    print("[INFO]                 such as --centralPairDEtaMax 999 --centralPairDPhiMin -999.")
    print()
    print(
        f"{'Selection':34s} {'Selected':>10s} {'Matched':>10s} {'Unmatched':>11s} "
        f"{'MatchEff%':>10s} {'UnmatchEff%':>12s} {'Purity%':>9s} {'UnmatchRej%':>12s}"
    )
    for row in rows:
        print(
            f"{row['selection']:34s} "
            f"{row['selected']:10d} "
            f"{row['matched']:10d} "
            f"{row['unmatched']:11d} "
            f"{row['matched_eff_pct']:10.2f} "
            f"{row['unmatched_eff_pct']:12.2f} "
            f"{row['matched_purity_pct']:9.2f} "
            f"{row['unmatched_rejection_pct']:12.2f}"
        )


def write_csv(path, rows):
    fieldnames = [
        "selection",
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
        for row in rows:
            writer.writerow(row)


def main():
    mass = infer_mass()
    dr_min, apt_max = get_working_point(mass)

    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open input file: {args.input}")

    tree = root_file.Get(args.tree)
    if tree is None:
        raise RuntimeError(f"Could not find tree '{args.tree}' in {args.input}")

    events = parse_events(tree, args.baseSelection)
    if not events:
        raise RuntimeError("No events passed the base selection. Check --input or --baseSelection.")

    base_matched = sum(1 for event in events if event["matched"])
    base_unmatched = len(events) - base_matched

    current_angular = lambda event: event["deta"] < args.currentDEtaMax and event["dphi"] > args.currentDPhiMin
    dr_apt = lambda event: event["dr"] > dr_min and event["apt"] < apt_max

    rows = [
        make_row("Stored loose pair", events, lambda event: True, base_matched, base_unmatched),
        make_row("Offline current angular", events, current_angular, base_matched, base_unmatched),
        make_row("Offline dR(m) + A_pt(m)", events, dr_apt, base_matched, base_unmatched),
        make_row(
            "Offline current angular + dR + A_pt",
            events,
            lambda event: current_angular(event) and dr_apt(event),
            base_matched,
            base_unmatched,
        ),
        make_row(
            "Offline drop dEta, keep dPhi + dR + A_pt",
            events,
            lambda event: event["dphi"] > args.currentDPhiMin and dr_apt(event),
            base_matched,
            base_unmatched,
        ),
        make_row(
            "Offline keep dEta, drop dPhi + dR + A_pt",
            events,
            lambda event: event["deta"] < args.currentDEtaMax and dr_apt(event),
            base_matched,
            base_unmatched,
        ),
        make_row(
            "Offline drop dEta and dPhi, keep dR + A_pt",
            events,
            dr_apt,
            base_matched,
            base_unmatched,
        ),
    ]

    print_table(rows, mass, dr_min, apt_max)

    if args.writeCSV:
        outdir = default_outdir(mass)
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, f"compare_massDependentCuts_M{mass}.csv")
        write_csv(outpath, rows)
        print(f"\n[INFO] Wrote CSV summary to: {outpath}")

    root_file.Close()


if __name__ == "__main__":
    main()
