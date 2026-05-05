#!/usr/bin/env python3

"""
Gen-level summary plotter for the ScoutNano Higgs study tree.

It produces:
 - generated Higgs kinematics
 - decay-product comparisons before/after minimal acceptance cuts
 - gen-jet multiplicity comparisons
 - forward gen-jet pair properties
 - a simple VBF acceptance cutflow based on forward gen-jets
"""

import os
import argparse

import ROOT

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Plot gen-level ScoutNano Higgs summary")
parser.add_argument(
    "--input",
    type=str,
    default="VBFHTo2B_M300_genHiggs_scoutNano.root",
    help="Input ROOT file containing the Events tree",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="plots_genVBFHiggs_scoutNano",
    help="Output directory for plots",
)
parser.add_argument(
    "--label",
    type=str,
    default="VBFHTo2B M300 (2024)",
    help="Label shown on plots",
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
    draw_expr = f"{expr}>>{name}"
    tree.Draw(draw_expr, selection, "goff")
    return hist


def style_hist(hist, color):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetLineWidth(3)


def save_single(hist, outdir, outname, label, logy=False):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()
        hist.SetMinimum(0.5)
    hist.Draw("hist")
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def save_overlay(hists, labels, outdir, outname, label, logy=False, normalize=False):
    canvas = ROOT.TCanvas(f"c_{outname}", "", 1000, 800)
    if logy:
        canvas.SetLogy()

    legend = ROOT.TLegend(0.5, 0.72, 0.88, 0.89)
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
        option = "hist" if idx == 0 else "hist same"
        hist.Draw(option)
        legend.AddEntry(hist, labels[idx], "l")

    legend.Draw()
    draw_label(label)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


def build_cutflow_hist(counts):
    hist = ROOT.TH1F("h_vbf_cutflow", ";Selection;Events", len(counts), 0.5, len(counts) + 0.5)
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


def save_cutflow(hist, counts, outdir, outname, label):
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
    title = box.AddText("Gen-level study: summary")
    title.SetTextFont(62)
    for idx, (step_label, value) in enumerate(counts):
        frac = 100.0 * value / total if total > 0 else 0.0
        line = f"{step_label}: {value} ({frac:.2f}%)"
        box.AddText(line)

    box.Draw()
    draw_label(label, cms_x=0.14, cms_y=0.845, extra_x=0.945, extra_y=0.912)
    #draw_label(label, cms_x=0.14, cms_y=0.845, extra_x=0.955, extra_y=0.905)
    canvas.SaveAs(os.path.join(outdir, f"{outname}.png"))
    canvas.SaveAs(os.path.join(outdir, f"{outname}.pdf"))


cms_style()
os.makedirs(args.outdir, exist_ok=True)

root_file = ROOT.TFile.Open(args.input)
if not root_file or root_file.IsZombie():
    raise RuntimeError(f"Could not open input file: {args.input}")

tree = root_file.Get("Events")
if tree is None:
    raise RuntimeError(f"Could not find tree 'Events' in {args.input}")

base_higgs = "higgs_mass > -900"
base_decay = "decay_mass > -900"
accepted_decay = (
    "decay_mass > -900"
    " && lead_decay_pt > 10 && sublead_decay_pt > 10"
    " && abs(lead_decay_eta) < 3 && abs(sublead_decay_eta) < 3"
)
has_forward_pair = "has_forward_genjet_pair == 1"
vbf_deta = "has_forward_genjet_pair == 1 && forward_genjet_pair_dEta > 5"
vbf_deta_pt10 = (
    "has_forward_genjet_pair == 1 && forward_genjet_pair_dEta > 5"
    " && forward_genjet1_pt > 10 && forward_genjet2_pt > 10"
)
vbf_deta_pt20 = (
    "has_forward_genjet_pair == 1 && forward_genjet_pair_dEta > 5"
    " && forward_genjet1_pt > 20 && forward_genjet2_pt > 20"
)
vbf_deta_pt20_decayacc = (
    "has_forward_genjet_pair == 1 && forward_genjet_pair_dEta > 5"
    " && forward_genjet1_pt > 20 && forward_genjet2_pt > 20"
    " && lead_decay_pt > 10 && sublead_decay_pt > 10"
    " && abs(lead_decay_eta) < 3 && abs(sublead_decay_eta) < 3"
)

# Higgs properties
h_higgs_pt = make_hist(tree, "h_higgs_pt", "higgs_pt", base_higgs, 80, 0, 800, "Generated Higgs p_{T} [GeV]")
h_higgs_eta = make_hist(tree, "h_higgs_eta", "higgs_eta", base_higgs, 60, -6, 6, "Generated Higgs #eta")
h_higgs_phi = make_hist(tree, "h_higgs_phi", "higgs_phi", base_higgs, 64, -3.2, 3.2, "Generated Higgs #phi")
h_higgs_mass = make_hist(tree, "h_higgs_mass", "higgs_mass", base_higgs, 80, 0, 1200, "Generated Higgs mass [GeV]")

for hist, outname in [
    (h_higgs_pt, "higgs_pt"),
    (h_higgs_eta, "higgs_eta"),
    (h_higgs_phi, "higgs_phi"),
    (h_higgs_mass, "higgs_mass"),
]:
    style_hist(hist, ROOT.kBlue + 1)
    save_single(hist, args.outdir, outname, args.label)

# Decay comparisons
decay_specs = [
    ("lead_decay_pt", 80, 0, 800, "Leading decay-product p_{T} [GeV]"),
    ("sublead_decay_pt", 80, 0, 800, "Subleading decay-product p_{T} [GeV]"),
    ("lead_decay_eta", 60, -6, 6, "Leading decay-product #eta"),
    ("sublead_decay_eta", 60, -6, 6, "Subleading decay-product #eta"),
    ("decay_dR", 60, 0, 6, "Decay-product #DeltaR"),
    ("decay_dEta", 60, 0, 6, "Decay-product |#Delta#eta|"),
    ("decay_dPhi", 64, 0, 3.2, "Decay-product |#Delta#phi|"),
    ("decay_mass", 80, 0, 1200, "Decay-product invariant mass [GeV]"),
]

for idx, (expr, nbins, xmin, xmax, xtitle) in enumerate(decay_specs):
    h_all = make_hist(tree, f"h_decay_all_{idx}", expr, base_decay, nbins, xmin, xmax, xtitle, "Normalized events")
    h_acc = make_hist(tree, f"h_decay_acc_{idx}", expr, accepted_decay, nbins, xmin, xmax, xtitle, "Normalized events")
    style_hist(h_all, ROOT.kBlack)
    style_hist(h_acc, ROOT.kRed + 1)
    save_overlay(
        [h_all, h_acc],
        ["No acceptance cuts", "gen p_{T} > 10 GeV, |#eta| < 3"],
        args.outdir,
        f"compare_{expr}",
        args.label,
        normalize=True,
    )

# Gen-jet multiplicities
h_n_genjets_total = make_hist(tree, "h_n_genjets_total", "n_genjets_total", "", 16, -0.5, 15.5, "GenJet multiplicity", "Normalized events")
h_n_genjets_central = make_hist(tree, "h_n_genjets_central", "n_genjets_central", "", 16, -0.5, 15.5, "GenJet multiplicity", "Normalized events")
h_n_genjets_forward = make_hist(tree, "h_n_genjets_forward", "n_genjets_forward", "", 16, -0.5, 15.5, "GenJet multiplicity", "Normalized events")
style_hist(h_n_genjets_total, ROOT.kBlack)
style_hist(h_n_genjets_central, ROOT.kBlue + 1)
style_hist(h_n_genjets_forward, ROOT.kRed + 1)
save_overlay(
    [h_n_genjets_total, h_n_genjets_central, h_n_genjets_forward],
    ["All GenJets", "|#eta| < 3", "3 #leq |#eta| #leq 5"],
    args.outdir,
    "genjet_multiplicities",
    args.label,
    normalize=True,
)

# Forward-pair properties
forward_specs = [
    ("forward_genjet1_pt", 80, 0, 300, "Forward GenJet 1 p_{T} [GeV]"),
    ("forward_genjet2_pt", 80, 0, 300, "Forward GenJet 2 p_{T} [GeV]"),
    ("forward_genjet1_eta", 50, -5.2, 5.2, "Forward GenJet 1 #eta"),
    ("forward_genjet2_eta", 50, -5.2, 5.2, "Forward GenJet 2 #eta"),
    ("forward_genjet_pair_dEta", 60, 0, 10, "Forward GenJet pair |#Delta#eta|"),
    ("forward_genjet_pair_dPhi", 64, 0, 3.2, "Forward GenJet pair |#Delta#phi|"),
    ("forward_genjet_pair_dR", 80, 0, 12, "Forward GenJet pair #DeltaR"),
    ("forward_genjet_pair_mass", 80, 0, 3000, "Forward GenJet pair mass [GeV]"),
]

for expr, nbins, xmin, xmax, xtitle in forward_specs:
    hist = make_hist(tree, f"h_{expr}", expr, has_forward_pair, nbins, xmin, xmax, xtitle)
    style_hist(hist, ROOT.kGreen + 2)
    save_single(hist, args.outdir, expr, args.label, logy=("mass" in expr or "_pt" in expr))

# VBF acceptance summary
n_all = tree.GetEntries()
n_pair = tree.GetEntries(has_forward_pair)
n_deta = tree.GetEntries(vbf_deta)
n_deta_pt10 = tree.GetEntries(vbf_deta_pt10)
n_deta_pt20 = tree.GetEntries(vbf_deta_pt20)
n_deta_pt20_decayacc = tree.GetEntries(vbf_deta_pt20_decayacc)

cutflow_counts = [
    ("Generated events", n_all),
    (">=2 fwd GenJets", n_pair),
    ("+ |#Delta#eta| > 5", n_deta),
    ("+ p_{T}^{1,2} > 10", n_deta_pt10),
    ("+ p_{T}^{1,2} > 20", n_deta_pt20),
    ("+ central pair acceptance", n_deta_pt20_decayacc),
]
h_cutflow = build_cutflow_hist(cutflow_counts)
save_cutflow(h_cutflow, cutflow_counts, args.outdir, "vbf_genjet_cutflow", args.label)

print("[INFO] Summary counts")
for label, value in cutflow_counts:
    frac = 100.0 * value / n_all if n_all > 0 else 0.0
    print(f"  {label:20s}: {value:8d} ({frac:6.2f}%)")

root_file.Close()
print(f"[INFO] Saved plots in: {args.outdir}")
