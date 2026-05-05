#!/usr/bin/env python3

"""
Summary plotter for scan_recoVBF_ptCuts.py CSV outputs.

It overlays the different scan families on the same canvas, one canvas per
mass point:
 - leading scaled pT
 - subleading scaled pT
 - pT imbalance

Each canvas contains two pads:
 - unmatched rejection vs matched efficiency
 - matched purity vs matched efficiency
"""

import os
import re
import csv
import glob
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Plot summary canvases for reco VBF pT-cut scans")
parser.add_argument(
    "--inputs",
    nargs="*",
    default=None,
    help="Scan directories to plot. Defaults to scan_ptCuts_recoVBF_M* in the current working directory.",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="plots_scan_ptCuts_recoVBF",
    help="Output directory for summary canvases",
)
parser.add_argument(
    "--label",
    type=str,
    default="",
    help="Global label shown on plots",
)
args = parser.parse_args()


SCAN_SPECS = [
    ("scan_scaledPt_leading.csv", "Leading p_{T}/m", ROOT.TColor.GetColor("#7FA8D1"), 20),
    ("scan_scaledPt_subleading.csv", "Subleading p_{T}/m", ROOT.TColor.GetColor("#E8A07C"), 21),
    ("scan_ptImbalance.csv", "p_{T} imbalance", ROOT.TColor.GetColor("#8DBB8B"), 22),
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


def draw_label(extra_text, cms_x=0.078, cms_y=0.9, extra_x=0.96, extra_y=0.9):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.042)
    latex.DrawLatex(cms_x, cms_y, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.034)
    latex2.DrawLatex(cms_x + 0.05, cms_y, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.036)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(extra_x, extra_y, extra_text)
    return latex, latex2, latex3


def resolve_input_dirs():
    if args.inputs:
        return sorted(args.inputs)
    return sorted(glob.glob(os.path.join(os.getcwd(), "scan_ptCuts_recoVBF_M*")))


def load_csv_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "threshold": float(row["threshold"]),
                    "matched_eff_pct": float(row["matched_eff_pct"]),
                    "matched_purity_pct": float(row["matched_purity_pct"]),
                    "unmatched_rejection_pct": float(row["unmatched_rejection_pct"]),
                }
            )
    return rows


def make_graph(rows, x_key, y_key, color, marker_style):
    graph = ROOT.TGraph(len(rows))
    for idx, row in enumerate(rows):
        graph.SetPoint(idx, row[x_key], row[y_key])
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(1.1)
    graph.SetLineWidth(3)
    return graph


def extract_mass_label(scan_dir):
    base = os.path.basename(os.path.normpath(scan_dir))
    match = re.search(r"_M(\d+)", base)
    return f"VBFH{match.group(1)}" if match else base


def draw_pad(pad, graphs, title, xtitle, ytitle, legend_coords):
    pad.cd()
    mg = ROOT.TMultiGraph()
    legend = ROOT.TLegend(*legend_coords)
    legend.SetTextSize(0.040)
    legend.SetFillStyle(0)
    legend.SetHeader(title, "C")

    for label, graph in graphs:
        mg.Add(graph, "LP")
        legend.AddEntry(graph, label, "lp")

    mg.Draw("A")
    mg.SetTitle(f";{xtitle};{ytitle}")
    mg.GetXaxis().SetLimits(0.0, 101.0)
    mg.SetMinimum(0.0)
    mg.SetMaximum(100.0)
    mg.GetXaxis().SetTitleSize(0.050)
    mg.GetYaxis().SetTitleSize(0.050)
    mg.GetXaxis().SetLabelSize(0.042)
    mg.GetYaxis().SetLabelSize(0.042)
    mg.GetYaxis().SetTitleOffset(1.15)

    legend.Draw()
    return mg, legend


def main():
    cms_style()
    os.makedirs(args.outdir, exist_ok=True)

    scan_dirs = resolve_input_dirs()
    if not scan_dirs:
        raise RuntimeError("No scan directories found. Pass them with --inputs or run from the scan directory area.")

    for scan_dir in scan_dirs:
        mass_label = extract_mass_label(scan_dir)

        roc_graphs = []
        purity_graphs = []
        missing = []
        canvas_objects = []

        for filename, label, color, marker_style in SCAN_SPECS:
            path = os.path.join(scan_dir, filename)
            if not os.path.exists(path):
                missing.append(filename)
                continue

            rows = load_csv_rows(path)
            roc_graphs.append(
                (label, make_graph(rows, "matched_eff_pct", "unmatched_rejection_pct", color, marker_style))
            )
            purity_graphs.append(
                (label, make_graph(rows, "matched_eff_pct", "matched_purity_pct", color, marker_style))
            )

        if not roc_graphs:
            print(f"[WARN] No scan CSVs found in {scan_dir}")
            continue

        if missing:
            print(f"[WARN] Missing CSVs in {scan_dir}: {', '.join(missing)}")

        canvas = ROOT.TCanvas(f"c_summary_{mass_label}", "", 1600, 800)
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
            roc_graphs,
            "Background Rejection",
            "Matched efficiency [%]",
            "Unmatched rejection [%]",
            (0.17, 0.16, 0.53, 0.40),
        ))
        canvas_objects.extend(draw_pad(
            right_pad,
            purity_graphs,
            "Matched Purity",
            "Matched efficiency [%]",
            "Matched purity [%]",
            (0.17, 0.16, 0.53, 0.40),
        ))

        canvas.cd()
        header_text = mass_label if not args.label else f"{args.label}  {mass_label}"
        canvas_objects.extend(draw_label(header_text, cms_x=0.078, cms_y=0.9, extra_x=0.96, extra_y=0.9))
        canvas.Modified()
        canvas.Update()

        outbase = os.path.join(args.outdir, f"summary_ptCutScans_{mass_label}")
        canvas.SaveAs(f"{outbase}.png")
        canvas.SaveAs(f"{outbase}.pdf")

    print(f"[INFO] Saved scan-summary plots in: {args.outdir}")


if __name__ == "__main__":
    main()
