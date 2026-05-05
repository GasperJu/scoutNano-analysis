#!/usr/bin/env python3

"""
Summary plotter for scan_recoVBF_massDependentCuts.py CSV outputs.

It makes one canvas per mass point with:
 - left pad: unmatched rejection vs matched efficiency
 - right pad: matched purity vs matched efficiency

Each point corresponds to one offline scenario from the comparison CSV.
"""

import os
import re
import csv
import glob
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Plot summary canvases for mass-dependent reco VBF cut comparisons")
parser.add_argument(
    "--inputs",
    nargs="*",
    default=None,
    help="Comparison CSV files or directories. Defaults to scan_massDependentCuts_M*/compare_massDependentCuts_M*.csv in the current working directory.",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="plots_scan_massDependentCuts_recoVBF",
    help="Output directory for summary canvases",
)
parser.add_argument(
    "--label",
    type=str,
    default="",
    help="Optional extra label shown in the canvas header",
)
parser.add_argument(
    "--preferredScenario",
    type=str,
    default="Offline dR(m) + A_pt(m)",
    help="Scenario to highlight with a circle and bold legend entry.",
)
args = parser.parse_args()


SCENARIO_STYLES = {
    "Stored loose pair": (ROOT.TColor.GetColor("#B8B8C8"), 20),
    "Offline current angular": (ROOT.TColor.GetColor("#7FA8D1"), 21),
    "Offline dR(m) + A_pt(m)": (ROOT.TColor.GetColor("#E8A07C"), 22),
    "Offline current angular + dR + A_pt": (ROOT.TColor.GetColor("#8DBB8B"), 23),
    "Offline drop dEta, keep dPhi + dR + A_pt": (ROOT.TColor.GetColor("#C8A2C8"), 33),
    "Offline keep dEta, drop dPhi + dR + A_pt": (ROOT.TColor.GetColor("#D6B46E"), 34),
    "Offline drop dEta and dPhi, keep dR + A_pt": (ROOT.TColor.GetColor("#5FA8A6"), 29),
}

SCENARIO_LABELS = {
    "Stored loose pair": "Loose pair",
    "Offline current angular": "Current angular",
    "Offline dR(m) + A_pt(m)": "dR(m) + A_{pT}(m)",
    "Offline current angular + dR + A_pt": "Current + dR + A_{pT}",
    "Offline drop dEta, keep dPhi + dR + A_pt": "Drop dEta",
    "Offline keep dEta, drop dPhi + dR + A_pt": "Drop dPhi",
    "Offline drop dEta and dPhi, keep dR + A_pt": "Drop dEta,dPhi",
}

DEFAULT_SCENARIOS = [
    "Stored loose pair",
    "Offline current angular",
    "Offline dR(m) + A_pt(m)",
    "Offline drop dEta, keep dPhi + dR + A_pt",
]


def cms_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.1, "X")
    ROOT.gStyle.SetTitleOffset(1.15, "Y")


def draw_label(extra_text, cms_x=0.075, cms_y=0.9, extra_x=0.96, extra_y=0.9):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.042)
    latex.DrawLatex(cms_x, cms_y, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.034)
    latex2.DrawLatex(cms_x + 0.06, cms_y, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.036)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(extra_x, extra_y, extra_text)
    return latex, latex2, latex3


def resolve_csv_inputs():
    if not args.inputs:
        return sorted(glob.glob(os.path.join(os.getcwd(), "scan_massDependentCuts_M*", "compare_massDependentCuts_M*.csv")))

    resolved = []
    for item in args.inputs:
        if os.path.isdir(item):
            resolved.extend(sorted(glob.glob(os.path.join(item, "compare_massDependentCuts_M*.csv"))))
        else:
            resolved.append(item)
    return sorted(resolved)


def load_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "selection": row["selection"],
                    "matched_eff_pct": float(row["matched_eff_pct"]),
                    "matched_purity_pct": float(row["matched_purity_pct"]),
                    "unmatched_rejection_pct": float(row["unmatched_rejection_pct"]),
                }
            )
    return rows


def filter_rows(rows):
    ordered = []
    for scenario in DEFAULT_SCENARIOS:
        for row in rows:
            if row["selection"] == scenario:
                ordered.append(row)
                break
    return ordered if ordered else rows


def extract_mass_label(path):
    match = re.search(r"_M(\d+)", os.path.basename(path))
    if match:
        return f"VBFH{match.group(1)}"
    return os.path.basename(path)


def make_graph_point(row):
    color, marker_style = SCENARIO_STYLES.get(row["selection"], (ROOT.kBlack, 20))
    graph = ROOT.TGraph(1)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(1.9)
    graph.SetLineWidth(3)
    return graph


def draw_pad(pad, rows, y_key, y_title, legend_header, preferred_scenario):
    pad.cd()
    frame = ROOT.TH2F(
        f"frame_{pad.GetName()}_{y_key}",
        f";Matched efficiency [%];{y_title}",
        100,
        0.0,
        101.0,
        100,
        0.0,
        100.0,
    )
    frame.Draw()

    legend = ROOT.TLegend(0.14, 0.14, 0.56, 0.54)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.040)
    legend.SetHeader(legend_header, "C")

    graphs = []
    highlight_objects = []
    for row in rows:
        graph = make_graph_point(row)
        graph.SetPoint(0, row["matched_eff_pct"], row[y_key])
        graph.Draw("P SAME")
        graphs.append(graph)
        legend_label = SCENARIO_LABELS.get(row["selection"], row["selection"])
        legend_value = f"{row[y_key]:.1f}%"
        if row["selection"] == preferred_scenario:
            legend_label = f"#bf{{{legend_label}}}"
            legend_value = f"#bf{{{legend_value}}}"
        legend.AddEntry(graph, f"{legend_label} ({legend_value})", "p")

        if row["selection"] == preferred_scenario:
            circle = ROOT.TEllipse(row["matched_eff_pct"], row[y_key], 2.0, 3.3)
            circle.SetFillStyle(0)
            circle.SetLineColor(ROOT.kBlack)
            circle.SetLineWidth(2)
            circle.Draw()
            highlight_objects.append(circle)

    legend.Draw()
    return frame, legend, graphs, highlight_objects


def main():
    cms_style()
    os.makedirs(args.outdir, exist_ok=True)

    csv_inputs = resolve_csv_inputs()
    if not csv_inputs:
        raise RuntimeError("No comparison CSVs found. Pass them with --inputs or run from the scan directory area.")

    for path in csv_inputs:
        rows = load_rows(path)
        if not rows:
            print(f"[WARN] Empty CSV: {path}")
            continue
        rows = filter_rows(rows)

        mass_label = extract_mass_label(path)
        canvas_objects = []

        canvas = ROOT.TCanvas(f"c_massdep_{mass_label}", "", 1700, 850)
        canvas.Divide(2, 1)

        left_pad = canvas.cd(1)
        left_pad.SetLeftMargin(0.14)
        left_pad.SetRightMargin(0.05)
        left_pad.SetBottomMargin(0.12)
        left_pad.SetTopMargin(0.12)

        right_pad = canvas.cd(2)
        right_pad.SetLeftMargin(0.14)
        right_pad.SetRightMargin(0.05)
        right_pad.SetBottomMargin(0.12)
        right_pad.SetTopMargin(0.12)

        canvas_objects.extend(draw_pad(
            left_pad,
            rows,
            "unmatched_rejection_pct",
            "Unmatched rejection [%]",
            "Background Rejection",
            args.preferredScenario,
        ))
        canvas_objects.extend(draw_pad(
            right_pad,
            rows,
            "matched_purity_pct",
            "Matched purity [%]",
            "Matched Purity",
            args.preferredScenario,
        ))

        canvas.cd()
        header_text = mass_label if not args.label else f"{args.label}  {mass_label}"
        canvas_objects.extend(draw_label(header_text))
        canvas.Modified()
        canvas.Update()

        outbase = os.path.join(args.outdir, f"summary_massDependentCuts_{mass_label}")
        canvas.SaveAs(f"{outbase}.png")
        canvas.SaveAs(f"{outbase}.pdf")

    print(f"[INFO] Saved mass-dependent scan summary plots in: {args.outdir}")


if __name__ == "__main__":
    main()
