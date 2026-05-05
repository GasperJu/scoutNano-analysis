#!/usr/bin/env python3

"""
Reco-level summary plotter for the ScoutNano VBF study tree.

It produces:
 - inclusive reco kinematics
 - matched vs not-matched comparisons
 - two reco cutflows ending in matched / not-matched categories
"""

import os
import argparse
import re

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Plot reco-level ScoutNano VBF summary")
parser.add_argument(
    "--input",
    type=str,
    default="VBFHTo2B_M300_recoVBF_scoutNano.root",
    help="Input ROOT file containing the Events tree",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="plots_recoVBF_scoutNano",
    help="Output directory for plots",
)
parser.add_argument(
    "--label",
    type=str,
    default="VBFHTo2B M300 ScoutNano",
    help="Label shown on plots",
)
parser.add_argument(
    "--fitOutput",
    type=str,
    default=None,
    help="Output ROOT file with post-selection mass histograms for fitting. Defaults to VBF-dijetMass-Histos_ForFIT/vbf-m<MASS>-<ALGO>.root when parsable from the input filename.",
)
args = parser.parse_args()


def cms_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.1, "X")
    ROOT.gStyle.SetTitleOffset(1.3, "Y")


def draw_label(extra_text, cms_x=0.14, cms_y=0.93, extra_x=0.95, extra_y=0.93):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    latex.DrawLatex(cms_x, cms_y, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.04)
    latex2.DrawLatex(cms_x + 0.10, cms_y, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(extra_x, extra_y, extra_text)


def make_hist(tree, name, expr, selection, nbins, xmin, xmax, xtitle, ytitle="Events"):
    hist = ROOT.TH1F(name, f";{xtitle};{ytitle}", nbins, xmin, xmax)
    hist.Sumw2()
    tree.Draw(f"{expr}>>{name}", selection, "goff")
    return hist


def style_hist(hist, color):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineWidth(3)


def save_single(hist, outdir, outname, label, logy=False, note_lines=None, plain_text=None, plain_text_coords=None, max_scale=1.35):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()
        hist.SetMinimum(0.5)
        hist.SetMaximum(hist.GetMaximum() * 20.0)
    else:
        hist.SetMaximum(hist.GetMaximum() * max_scale)
    hist.Draw("hist")
    note_box = None
    if note_lines:
        note_box = draw_note(note_lines)
    plain = None
    if plain_text:
        if plain_text_coords is None:
            plain = draw_text_label(plain_text)
        else:
            plain = draw_text_label(plain_text, *plain_text_coords)
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def draw_note(lines, x1=0.58, y1=0.62, x2=0.90, y2=0.84):
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.18)
    box.SetFillStyle(1001)
    box.SetLineColor(ROOT.kGray + 1)
    box.SetLineWidth(2)
    box.SetShadowColor(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.032)
    for line in lines:
        box.AddText(line)
    box.Draw()
    return box


def draw_text_label(text, x=0.95, y=0.84, align=31, font=62, size=0.032):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(align)
    latex.SetTextFont(font)
    latex.SetTextSize(size)
    latex.DrawLatex(x, y, text)
    return latex


def save_overlay(hists, labels, outdir, outname, label, logy=False, normalize=False, note_lines=None, legend_coords=None, note_coords=None):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()

    if legend_coords is None:
        legend_coords = (0.56, 0.72, 0.90, 0.89)
    legend = ROOT.TLegend(*legend_coords)
    legend.SetTextSize(0.035)
    legend.SetFillStyle(0)

    max_y = 0.0
    for hist in hists:
        if normalize and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())
        max_y = max(max_y, hist.GetMaximum())

    for idx, hist in enumerate(hists):
        hist.SetMaximum(max_y * (20.0 if logy else 1.35))
        if logy:
            hist.SetMinimum(1e-4 if normalize else 0.5)
        hist.Draw("hist" if idx == 0 else "hist same")
        legend.AddEntry(hist, labels[idx], "l")

    legend.Draw()
    note_box = None
    if note_lines:
        if note_coords is None:
            note_box = draw_note(note_lines)
        else:
            note_box = draw_note(note_lines, *note_coords)
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def save_overlay_with_vline(hists, labels, outdir, outname, label, vline_x, vline_label=None, logy=False, normalize=False, note_lines=None, legend_coords=None, note_coords=None, x_range=None):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()

    if legend_coords is None:
        legend_coords = (0.52, 0.72, 0.88, 0.89)
    legend = ROOT.TLegend(*legend_coords)
    legend.SetTextSize(0.035)
    legend.SetFillStyle(0)

    max_y = 0.0
    for hist in hists:
        if normalize and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())
        max_y = max(max_y, hist.GetMaximum())
        if x_range is not None:
            hist.GetXaxis().SetRangeUser(*x_range)

    for idx, hist in enumerate(hists):
        hist.SetMaximum(max_y * (22.0 if logy else 1.40))
        if logy:
            hist.SetMinimum(1e-4 if normalize else 0.5)
        hist.Draw("hist" if idx == 0 else "hist same")
        legend.AddEntry(hist, labels[idx], "l")

    line = ROOT.TLine(vline_x, hists[0].GetMinimum() if logy else 0.0, vline_x, max_y * (22.0 if logy else 1.40))
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(7)
    line.SetLineWidth(3)
    line.Draw()

    legend.AddEntry(line, vline_label if vline_label else f"x = {vline_x}", "l")
    legend.Draw()

    note_box = None
    if note_lines:
        if note_coords is None:
            note_box = draw_note(note_lines)
        else:
            note_box = draw_note(note_lines, *note_coords)

    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def save_2d_with_thresholds(tree, outdir, outname, label, x_expr, y_expr, selection, nbins_x, xmin, xmax, nbins_y, ymin, ymax, xtitle, ytitle, threshold=0.4):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 850)
    hist = ROOT.TH2F(outname, f";{xtitle};{ytitle}", nbins_x, xmin, xmax, nbins_y, ymin, ymax)
    tree.Draw(f"{y_expr}:{x_expr}>>{outname}", selection, "colz")
    hist = ROOT.gDirectory.Get(outname)
    hist.SetTitle(f";{xtitle};{ytitle}")
    hist.GetZaxis().SetTitle("Events")
    hist.GetZaxis().SetTitleOffset(1.2)

    line_v = ROOT.TLine(threshold, ymin, threshold, ymax)
    line_h = ROOT.TLine(xmin, threshold, xmax, threshold)
    for line in [line_v, line_h]:
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(7)
        line.SetLineWidth(3)
        line.Draw()

    note_box = draw_note(
        [
            "Matched-event region",
            "#DeltaR_{1} < 0.4 and #DeltaR_{2} < 0.4",
        ],
        0.56, 0.75, 0.89, 0.87
    )
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def build_cutflow_hist(name, counts):
    hist = ROOT.TH1F(name, ";Selection;Events", len(counts), 0.5, len(counts) + 0.5)
    for idx, (label, value) in enumerate(counts, start=1):
        hist.SetBinContent(idx, value)
        hist.GetXaxis().SetBinLabel(idx, label)
    hist.SetFillColor(ROOT.kAzure - 9)
    hist.SetLineColor(ROOT.kAzure + 2)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetXaxis().SetLabelOffset(0.02)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleOffset(2.1)
    return hist


def save_cutflow(hist, counts, outdir, outname, label, title):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1200, 1000)
    canvas.SetBottomMargin(0.24)
    canvas.SetTopMargin(0.10)
    canvas.SetRightMargin(0.06)
    hist.SetBarWidth(0.8)
    hist.SetBarOffset(0.1)
    hist.SetMaximum(hist.GetMaximum() * 1.35)
    hist.Draw("hist text0")
    hist.LabelsOption("u", "X")

    total = float(counts[0][1]) if counts else 0.0
    box = ROOT.TPaveText(0.56, 0.54, 0.91, 0.82, "NDC")
    box.SetFillColor(ROOT.kGray)
    box.SetFillColorAlpha(ROOT.kGray, 0.15)
    box.SetFillStyle(1001)
    box.SetLineColor(ROOT.kGray + 2)
    box.SetLineWidth(2)
    box.SetShadowColor(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(0.026)
    title_text = box.AddText(title)
    title_text.SetTextFont(62)
    for step_label, value in counts:
        frac = 100.0 * value / total if total > 0 else 0.0
        box.AddText(f"{step_label}: {value} ({frac:.2f}%)")
    box.Draw()

    draw_label(label, cms_x=0.14, cms_y=0.845, extra_x=0.955, extra_y=0.905)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def default_fit_output(outdir):
    if args.fitOutput:
        return args.fitOutput
    fit_dir = os.path.join(os.getcwd(), "VBF-dijetMass-Histos_ForFIT")
    input_base = os.path.basename(args.input)
    input_stem = os.path.splitext(input_base)[0]

    mass_match = re.search(r"_M(\d+)", input_stem)
    mass_token = f"m{mass_match.group(1)}" if mass_match else "mUnknown"

    algo_match = re.search(r"(leading|minDR|minDPhi|maxPtSum|massClosest)", input_stem)
    algo_token = algo_match.group(1) if algo_match else "unknownAlgo"

    return os.path.join(fit_dir, f"vbf-{mass_token}-{algo_token}.root")


cms_style()
os.makedirs(args.outdir, exist_ok=True)

root_file = ROOT.TFile.Open(args.input)
if not root_file or root_file.IsZombie():
    raise RuntimeError(f"Could not open input file: {args.input}")

tree = root_file.Get("Events")
if tree is None:
    raise RuntimeError(f"Could not find tree 'Events' in {args.input}")

base_reco = "pass_trigger_baseline == 1"
has_forward_pair = "pass_trigger_baseline == 1 && has_forward_pair == 1"
has_vbf_tag = "pass_trigger_baseline == 1 && has_vbf_tag == 1"
has_central_pair = "pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1"
matched = "pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1 && central_pair_genmatched == 1"
unmatched = "pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1 && central_pair_genmatched == 0"

# Inclusive reco plots
h_n_reco_total_raw = make_hist(tree, "h_n_reco_total_raw", "n_recojets_total_raw", base_reco, 16, -0.5, 15.5, "Reco-jet multiplicity", "Normalized events")
h_n_reco_central_raw = make_hist(tree, "h_n_reco_central_raw", "n_recojets_central_raw", base_reco, 16, -0.5, 15.5, "Reco-jet multiplicity", "Normalized events")
h_n_reco_forward_raw = make_hist(tree, "h_n_reco_forward_raw", "n_recojets_forward_raw", base_reco, 16, -0.5, 15.5, "Reco-jet multiplicity", "Normalized events")
h_n_reco_total = make_hist(tree, "h_n_reco_total", "n_recojets_total", base_reco, 16, -0.5, 15.5, "Reco-jet multiplicity", "Normalized events")
h_n_reco_central = make_hist(tree, "h_n_reco_central", "n_recojets_central", base_reco, 16, -0.5, 15.5, "Reco-jet multiplicity", "Normalized events")
h_n_reco_forward = make_hist(tree, "h_n_reco_forward", "n_recojets_forward", base_reco, 16, -0.5, 15.5, "Reco-jet multiplicity", "Normalized events")
style_hist(h_n_reco_total, ROOT.kSpring + 3)
style_hist(h_n_reco_central, ROOT.kBlue + 1)
style_hist(h_n_reco_forward, ROOT.kRed + 1)
style_hist(h_n_reco_total_raw, ROOT.kGreen + 1)
style_hist(h_n_reco_central_raw, ROOT.kCyan + 2)
style_hist(h_n_reco_forward_raw, ROOT.kOrange + 7)
h_n_reco_total.SetLineWidth(4)
h_n_reco_central.SetLineWidth(4)
h_n_reco_forward.SetLineWidth(4)
h_n_reco_total_raw.SetLineWidth(4)
h_n_reco_central_raw.SetLineWidth(4)
h_n_reco_forward_raw.SetLineWidth(4)
h_n_reco_total_raw.SetLineStyle(2)
h_n_reco_central_raw.SetLineStyle(2)
h_n_reco_forward_raw.SetLineStyle(2)
save_overlay(
    [h_n_reco_total, h_n_reco_central, h_n_reco_forward],
    ["All reco jets", "Central reco jets", "Forward reco jets"],
    args.outdir,
    "recojet_multiplicities",
    args.label,
    normalize=True,
    legend_coords=(0.20, 0.70, 0.52, 0.88),
)
save_overlay(
    [h_n_reco_total_raw, h_n_reco_total],
    ["Raw all reco jets", "Corrected all reco jets"],
    args.outdir,
    "compare_recojet_mult_total_raw_vs_corr",
    args.label,
    normalize=True,
    legend_coords=(0.20, 0.72, 0.52, 0.88),
)
save_overlay(
    [h_n_reco_central_raw, h_n_reco_central],
    ["Raw central reco jets", "Corrected central reco jets"],
    args.outdir,
    "compare_recojet_mult_central_raw_vs_corr",
    args.label,
    normalize=True,
    legend_coords=(0.20, 0.72, 0.52, 0.88),
)
raw_fwd_ge2 = tree.GetEntries("n_recojets_forward_raw >= 2")
corr_fwd_ge2 = tree.GetEntries("n_recojets_forward >= 2")
fwd_reduction_pct = 100.0 * (raw_fwd_ge2 - corr_fwd_ge2) / raw_fwd_ge2 if raw_fwd_ge2 > 0 else 0.0
save_overlay(
    [h_n_reco_forward_raw, h_n_reco_forward],
    ["Raw forward reco jets", "Corrected forward reco jets"],
    args.outdir,
    "compare_recojet_mult_forward_raw_vs_corr",
    args.label,
    normalize=True,
    legend_coords=(0.25, 0.72, 0.6, 0.88),
    note_coords=(0.5, 0.4, 0.88, 0.65),
    note_lines=[
        "VBF-like proxy: >= 2 forward jets",
        f"Raw: {raw_fwd_ge2} events",
        f"Corrected: {corr_fwd_ge2} events",
        f"Reduction after HLT JECs: {fwd_reduction_pct:.1f}%",
    ],
)

h_n_reco_total_diff = make_hist(tree, "h_n_reco_total_diff", "n_recojets_total-n_recojets_total_raw", base_reco, 21, -10.5, 10.5, "#DeltaN reco jets (corr - raw)")
h_n_reco_central_diff = make_hist(tree, "h_n_reco_central_diff", "n_recojets_central-n_recojets_central_raw", base_reco, 21, -10.5, 10.5, "#DeltaN central reco jets (corr - raw)")
h_n_reco_forward_diff = make_hist(tree, "h_n_reco_forward_diff", "n_recojets_forward-n_recojets_forward_raw", base_reco, 21, -10.5, 10.5, "#DeltaN forward reco jets (corr - raw)")
style_hist(h_n_reco_total_diff, ROOT.kBlack)
style_hist(h_n_reco_central_diff, ROOT.kBlue + 1)
style_hist(h_n_reco_forward_diff, ROOT.kRed + 1)
save_overlay(
    [h_n_reco_total_diff, h_n_reco_central_diff, h_n_reco_forward_diff],
    ["All reco jets", "Central reco jets", "Forward reco jets"],
    args.outdir,
    "compare_recojet_mult_diff_raw_vs_corr",
    args.label,
    normalize=True,
    legend_coords=(0.20, 0.70, 0.52, 0.88),
)

inclusive_specs = [
    ("forward_pair_dEta", has_forward_pair, 60, 0, 10, "Forward reco pair |#Delta#eta|"),
    ("forward_pair_mass", has_forward_pair, 100, 0, 5000, "Forward reco pair mass [GeV]"),
    ("central_jet1_pt", has_central_pair, 100, 0, 1000, "Leading central reco-jet p_{T} [GeV]"),
    ("central_jet2_pt", has_central_pair, 100, 0, 1000, "Subleading central reco-jet p_{T} [GeV]"),
    ("central_jet1_eta", has_central_pair, 60, -3.2, 3.2, "Leading central reco-jet #eta"),
    ("central_jet2_eta", has_central_pair, 60, -3.2, 3.2, "Subleading central reco-jet #eta"),
    ("central_pair_dR", has_central_pair, 60, 0, 6, "Central reco pair #DeltaR"),
    ("central_pair_mass", has_central_pair, 100, 0, 3000, "Central reco pair mass [GeV]"),
]

for expr, selection, nbins, xmin, xmax, xtitle in inclusive_specs:
    hist = make_hist(tree, f"h_{expr}", expr, selection, nbins, xmin, xmax, xtitle)
    style_hist(hist, ROOT.kBlue + 1)
    save_single(hist, args.outdir, expr, args.label, logy=("mass" in expr or "_pt" in expr))

# Inclusive plots for the final reco-selected sample used for analysis-like studies
selected_specs = [
    ("forward_pair_dEta", 60, 0, 10, "Forward reco pair |#Delta#eta|"),
    ("forward_pair_mass", 100, 0, 5000, "Forward reco pair mass [GeV]"),
    ("central_jet1_pt", 100, 0, 1000, "Leading central reco-jet p_{T} [GeV]"),
    ("central_jet2_pt", 100, 0, 1000, "Subleading central reco-jet p_{T} [GeV]"),
    ("central_jet1_eta", 60, -3.2, 3.2, "Leading central reco-jet #eta"),
    ("central_jet2_eta", 60, -3.2, 3.2, "Subleading central reco-jet #eta"),
    ("central_pair_dR", 60, 0, 6, "Central reco pair #DeltaR"),
    ("central_pair_dEta", 60, 0, 6, "Central reco pair |#Delta#eta|"),
    ("central_pair_dPhi", 64, 0, 3.2, "Central reco pair |#Delta#phi|"),
    ("central_pair_mass", 100, 0, 3000, "Central reco pair mass [GeV]"),
]

for expr, nbins, xmin, xmax, xtitle in selected_specs:
    hist = make_hist(tree, f"h_selected_{expr}", expr, has_central_pair, nbins, xmin, xmax, xtitle)
    style_hist(hist, ROOT.kMagenta + 1)
    save_single(
        hist,
        args.outdir,
        f"selected_{expr}",
        args.label,
        logy=False,
        plain_text="pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1",
        plain_text_coords=(0.95, 0.84, 31, 62, 0.030),
        max_scale=1.55,
    )

# Matched vs unmatched comparisons
compare_specs = [
    ("forward_pair_dEta", 60, 0, 10, "Forward reco pair |#Delta#eta|"),
    ("forward_pair_mass", 100, 0, 5000, "Forward reco pair mass [GeV]"),
    ("central_jet1_pt", 100, 0, 1000, "Leading central reco-jet p_{T} [GeV]"),
    ("central_jet2_pt", 100, 0, 1000, "Subleading central reco-jet p_{T} [GeV]"),
    ("central_jet1_eta", 60, -3.2, 3.2, "Leading central reco-jet #eta"),
    ("central_jet2_eta", 60, -3.2, 3.2, "Subleading central reco-jet #eta"),
    ("central_pair_dR", 60, 0, 6, "Central reco pair #DeltaR"),
    ("central_pair_dEta", 60, 0, 6, "Central reco pair |#Delta#eta|"),
    ("central_pair_dPhi", 64, 0, 3.2, "Central reco pair |#Delta#phi|"),
    ("central_pair_mass", 100, 0, 3000, "Central reco pair mass [GeV]"),
]

for idx, (expr, nbins, xmin, xmax, xtitle) in enumerate(compare_specs):
    h_match = make_hist(tree, f"h_match_{idx}", expr, matched, nbins, xmin, xmax, xtitle, "Normalized events")
    h_unmatch = make_hist(tree, f"h_unmatch_{idx}", expr, unmatched, nbins, xmin, xmax, xtitle, "Normalized events")
    style_hist(h_match, ROOT.kRed + 1)
    style_hist(h_unmatch, ROOT.kBlack)
    save_overlay(
        [h_match, h_unmatch],
        ["Gen matched", "Not matched"],
        args.outdir,
        f"compare_{expr}",
        args.label,
        normalize=True,
        logy=False,
    )

# Matching distance plots
h_match_dr1 = make_hist(tree, "h_match_dr1", "central_pair_match_dR1", "pass_trigger_baseline == 1 && has_central_pair == 1 && has_two_decay_products == 1", 60, 0, 3, "Best-match #DeltaR(j, q)")
h_match_dr2 = make_hist(tree, "h_match_dr2", "central_pair_match_dR2", "pass_trigger_baseline == 1 && has_central_pair == 1 && has_two_decay_products == 1", 60, 0, 3, "Best-match #DeltaR(j, q)")
style_hist(h_match_dr1, ROOT.kAzure + 2)
style_hist(h_match_dr2, ROOT.kOrange + 7)
h_match_dr1.SetLineWidth(4)
h_match_dr2.SetLineWidth(4)
save_overlay_with_vline(
    [h_match_dr1, h_match_dr2],
    ["Leading central reco jet", "Subleading central reco jet"],
    args.outdir,
    "compare_matching_dR",
    args.label,
    vline_x=0.4,
    vline_label="Matching threshold #DeltaR = 0.4",
    normalize=True,
    logy=True,
    note_lines=[
        "Matching criterion",
        "central_pair_genmatched == 1 if",
        "#DeltaR_{1} < 0.4 and #DeltaR_{2} < 0.4",
    ],
    note_coords=(0.52, 0.46, 0.90, 0.70),
    legend_coords=(0.52, 0.72, 0.90, 0.89),
    x_range=(0.0, 1.5),
)
save_2d_with_thresholds(
    tree,
    args.outdir,
    "compare_matching_dR2D",
    args.label,
    "central_pair_match_dR1",
    "central_pair_match_dR2",
    "pass_trigger_baseline == 1 && has_central_pair == 1 && has_two_decay_products == 1",
    50, 0.0, 1.5,
    50, 0.0, 1.5,
    "Best-match #DeltaR(j, q) for leading central reco jet",
    "Best-match #DeltaR(j, q) for subleading central reco jet",
    threshold=0.4,
)

# Cutflows
n_all = tree.GetEntries()
n_fwd = tree.GetEntries(has_forward_pair)
n_vbf = tree.GetEntries(has_vbf_tag)
n_cent = tree.GetEntries(has_central_pair)
n_match = tree.GetEntries(matched)
n_unmatch = tree.GetEntries(unmatched)

matched_cutflow = [
    ("Reco events", n_all),
    (">=2 forward jets", n_fwd),
    ("+ VBF tag", n_vbf),
    ("+ central pair", n_cent),
    ("+ gen matched", n_match),
]
unmatched_cutflow = [
    ("Reco events", n_all),
    (">=2 forward jets", n_fwd),
    ("+ VBF tag", n_vbf),
    ("+ central pair", n_cent),
    ("+ not matched", n_unmatch),
]

h_cutflow_match = build_cutflow_hist("h_cutflow_match", matched_cutflow)
h_cutflow_unmatch = build_cutflow_hist("h_cutflow_unmatch", unmatched_cutflow)
save_cutflow(h_cutflow_match, matched_cutflow, args.outdir, "cutflow_matched", args.label, "Reco-level study: matched")
save_cutflow(h_cutflow_unmatch, unmatched_cutflow, args.outdir, "cutflow_unmatched", args.label, "Reco-level study: not matched")

print("[INFO] Matched cutflow")
for step, value in matched_cutflow:
    frac = 100.0 * value / n_all if n_all > 0 else 0.0
    print(f"  {step:18s}: {value:8d} ({frac:6.2f}%)")

print("[INFO] Not-matched cutflow")
for step, value in unmatched_cutflow:
    frac = 100.0 * value / n_all if n_all > 0 else 0.0
    print(f"  {step:18s}: {value:8d} ({frac:6.2f}%)")

fit_output = default_fit_output(args.outdir)
fit_output_dir = os.path.dirname(fit_output)
if fit_output_dir:
    os.makedirs(fit_output_dir, exist_ok=True)
fit_file = ROOT.TFile(fit_output, "RECREATE")

h_signal_peak_selected = make_hist(
    tree,
    "h_signal_peak_selected_5GeV",
    "central_pair_mass",
    has_central_pair,
    600,
    0,
    3000,
    "Central reco pair mass [GeV]",
)
h_signal_peak_selected.SetTitle(";Central reco pair mass [GeV];Events / 5 GeV")
h_signal_peak_selected.Write()

h_signal_peak_matched = make_hist(
    tree,
    "h_signal_peak_matched_5GeV",
    "central_pair_mass",
    matched,
    600,
    0,
    3000,
    "Central reco pair mass [GeV]",
)
h_signal_peak_matched.SetTitle(";Central reco pair mass [GeV];Events / 5 GeV")
h_signal_peak_matched.Write()

fit_file.Close()

root_file.Close()
print(f"[INFO] Saved plots in: {args.outdir}")
print(f"[INFO] Saved fit histograms in: {fit_output}")
