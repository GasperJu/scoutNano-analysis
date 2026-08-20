#!/usr/bin/env python3

import argparse
import csv
import math
import os
import re

import ROOT

from fit_signalPeak_VBF import (
    COLORS,
    compute_fwhm_from_roofit_pdf,
    discover_files,
    draw_cms_label,
    effective_sigma,
    get_draw_range,
    get_fit_window,
    sanitize_token,
    set_style,
)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

FLAVOR_ORDER = ["bb", "cc", "qq", "gluglu"]
FLAVOR_COLORS = {
    "bb": ROOT.kBlue + 1,
    "cc": ROOT.kRed + 1,
    "qq": ROOT.kGreen + 2,
    "gluglu": ROOT.kMagenta + 2,
    "unknown": ROOT.kGray + 2,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare VBF signal masses across flavors and run per-flavor fits."
    )
    parser.add_argument("--indir", type=str, default="VBF-dijetMass-Histos_ForFIT")
    parser.add_argument("--outdir", type=str, default="plots_VBFFlavorMassFits")
    parser.add_argument("--hist", type=str, default="h_signal_peak_selected_5GeV")
    parser.add_argument("--algo", type=str, default=None)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument(
        "--flavors",
        type=str,
        default="bb,cc,qq,gluglu",
        help="Comma-separated flavors to process (bb,cc,qq,gluglu).",
    )
    parser.add_argument("--norm", action="store_true", help="Normalize comparison overlays.")
    parser.add_argument(
        "--reference-flavor",
        type=str,
        default="bb",
        help="Reference flavor used for ratio panels.",
    )
    parser.add_argument(
        "--no-ratio",
        action="store_true",
        help="Disable ratio panel in flavor comparison plots.",
    )
    parser.add_argument("--ratio-ymin", type=float, default=0.5)
    parser.add_argument("--ratio-ymax", type=float, default=1.5)
    parser.add_argument(
        "--summary-name",
        type=str,
        default="flavor_mass_fit_summary.csv",
        help="Output summary table filename.",
    )
    return parser.parse_args()


def infer_flavor(token):
    lowered = token.lower()
    checks = [
        ("gluglu", [r"(?:^|[^a-z0-9])(?:2glu|gluglu|gg)(?:[^a-z0-9]|$)"]),
        ("bb", [r"(?:^|[^a-z0-9])(?:2b|bb)(?:[^a-z0-9]|$)"]),
        ("cc", [r"(?:^|[^a-z0-9])(?:2c|cc)(?:[^a-z0-9]|$)"]),
        ("qq", [r"(?:^|[^a-z0-9])(?:2q|qq)(?:[^a-z0-9]|$)"]),
    ]
    for flavor, patterns in checks:
        for pattern in patterns:
            if re.search(pattern, lowered):
                return flavor
    return "unknown"


def load_histogram(path, hist_name, suffix):
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        print(f"[WARNING] Could not open file: {path}")
        return None
    hist = root_file.Get(hist_name)
    if not hist or not hist.InheritsFrom("TH1"):
        print(f"[WARNING] Missing TH1 histogram {hist_name} in {path}")
        root_file.Close()
        return None
    cloned = hist.Clone(f"{hist.GetName()}_{suffix}")
    cloned.SetDirectory(0)
    root_file.Close()
    if cloned.Integral() <= 0:
        print(f"[WARNING] Empty histogram in {path}")
        return None
    return cloned


def count_events_in_window(hist, xlow, xhigh):
    if xlow is None or xhigh is None or xhigh <= xlow:
        return 0.0
    axis = hist.GetXaxis()
    lo_bin = axis.FindBin(max(axis.GetXmin(), xlow + 1e-6))
    hi_bin = axis.FindBin(min(axis.GetXmax(), xhigh - 1e-6))
    if hi_bin < lo_bin:
        return 0.0
    return hist.Integral(lo_bin, hi_bin)


def fit_flavor_mass_histogram(info, hist, outdir):
    sigma_eff, x_eff_low, x_eff_high = effective_sigma(hist, frac=0.683)
    hist_xmin = hist.GetXaxis().GetXmin()
    hist_xmax = hist.GetXaxis().GetXmax()
    fit_xmin, fit_xmax = get_fit_window(info["mass"], hist_xmin, hist_xmax)
    draw_xmin, draw_xmax = get_draw_range(info["mass"], hist_xmin, hist_xmax)

    token = f"{sanitize_token(info['flavor'])}_{info['mass']}_{sanitize_token(info['algo'])}_{sanitize_token(info['tag'])}"
    xvar = ROOT.RooRealVar(f"x_{token}", "m_{jj} [GeV]", fit_xmin, fit_xmax)
    datahist = ROOT.RooDataHist(
        f"datahist_{token}",
        f"datahist_{token}",
        ROOT.RooArgList(xvar),
        ROOT.RooFit.Import(hist),
    )

    mean = ROOT.RooRealVar(f"mean_{token}", "mean", info["mass"], 0.75 * info["mass"], 1.25 * info["mass"])
    sigma_left = ROOT.RooRealVar(f"sigmaL_{token}", "sigmaL", max(10.0, 0.08 * info["mass"]), 10.0, 600.0)
    sigma_right = ROOT.RooRealVar(f"sigmaR_{token}", "sigmaR", max(12.0, 0.12 * info["mass"]), 1.0, 500.0)
    alpha_left = ROOT.RooRealVar(f"alphaL_{token}", "alphaL", 5.0, 0.05, 30.0)
    n_left = ROOT.RooRealVar(f"nL_{token}", "nL", 1.0, 0.5, 50.0)
    alpha_right = ROOT.RooRealVar(f"alphaR_{token}", "alphaR", 1.0, 0.01, 5.0)
    n_right = ROOT.RooRealVar(f"nR_{token}", "nR", 2.0, 0.1, 50.0)

    model = ROOT.RooCrystalBall(
        f"dscb_{token}",
        f"dscb_{token}",
        xvar,
        mean,
        sigma_left,
        sigma_right,
        alpha_left,
        n_left,
        alpha_right,
        n_right,
    )

    fitres = model.fitTo(
        datahist,
        ROOT.RooFit.Save(True),
        ROOT.RooFit.PrintLevel(-1),
        ROOT.RooFit.Strategy(1),
    )

    fwhm_info = compute_fwhm_from_roofit_pdf(model, xvar, fit_xmin, fit_xmax)
    sigma_eq = fwhm_info["sigma_eq"] if fwhm_info else None
    rel_width_pct = None
    if sigma_eq is not None and mean.getVal() != 0:
        rel_width_pct = 100.0 * sigma_eq / mean.getVal()

    total_events = hist.Integral(1, hist.GetNbinsX())
    n_eff = count_events_in_window(hist, x_eff_low, x_eff_high)
    sigma_for_windows = sigma_eq if sigma_eq is not None else 0.5 * (sigma_left.getVal() + sigma_right.getVal())
    x1sig_low = mean.getVal() - sigma_for_windows
    x1sig_high = mean.getVal() + sigma_for_windows
    x3sig_low = mean.getVal() - 3.0 * sigma_for_windows
    x3sig_high = mean.getVal() + 3.0 * sigma_for_windows
    n_1sig = count_events_in_window(hist, x1sig_low, x1sig_high)
    n_3sig = count_events_in_window(hist, x3sig_low, x3sig_high)
    eff_68 = (100.0 * n_eff / total_events) if total_events > 0 else None
    eff_1sig = (100.0 * n_1sig / total_events) if total_events > 0 else None
    eff_3sig = (100.0 * n_3sig / total_events) if total_events > 0 else None

    frame = xvar.frame()
    datahist.plotOn(
        frame,
        ROOT.RooFit.Name("data"),
        ROOT.RooFit.MarkerStyle(20),
        ROOT.RooFit.MarkerSize(0.9),
        ROOT.RooFit.LineColor(ROOT.kBlack),
    )
    model.plotOn(
        frame,
        ROOT.RooFit.Name("model"),
        ROOT.RooFit.LineColor(FLAVOR_COLORS.get(info["flavor"], COLORS.get(info["mass"], ROOT.kBlue + 1))),
        ROOT.RooFit.LineWidth(3),
    )

    canvas = ROOT.TCanvas(f"c_fit_{token}", "", 1200, 1000)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)

    frame.SetTitle("")
    frame.GetXaxis().SetTitle("m_{jj} [GeV]")
    frame.GetYaxis().SetTitle("Events / 5 GeV")
    frame.GetXaxis().SetRangeUser(draw_xmin, draw_xmax)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.Draw()

    legend = ROOT.TLegend(0.62, 0.68, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.SetHeader(f"{info['flavor']} M={info['mass']} GeV", "C")
    legend.AddEntry(frame.findObject("data"), "Data", "lep")
    legend.AddEntry(frame.findObject("model"), "Double-CB fit", "l")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.030)
    x0 = 0.58
    y0 = 0.60
    dy = 0.043
    text.DrawLatex(x0, y0, f"#mu = {mean.getVal():.2f} #pm {mean.getError():.2f} GeV")
    text.DrawLatex(x0, y0 - dy, f"#sigma_{{L}} = {sigma_left.getVal():.2f} #pm {sigma_left.getError():.2f}")
    text.DrawLatex(x0, y0 - 2 * dy, f"#sigma_{{R}} = {sigma_right.getVal():.2f} #pm {sigma_right.getError():.2f}")
    text.DrawLatex(x0, y0 - 3 * dy, f"#sigma_{{eff}} = {sigma_eff:.1f} GeV")
    text.DrawLatex(x0, y0 - 4 * dy, f"fit status/cov = {fitres.status()}/{fitres.covQual()}")

    if x_eff_low is not None and x_eff_high is not None:
        text.DrawLatex(0.58, 0.38, f"68.3% = [{x_eff_low:.0f}, {x_eff_high:.0f}] GeV")
    if fwhm_info is not None:
        text.DrawLatex(0.58, 0.34, f"FWHM = {fwhm_info['fwhm']:.1f} GeV")
    if rel_width_pct is not None:
        text.DrawLatex(0.58, 0.30, f"#sigma/m = {rel_width_pct:.2f}%")

    draw_cms_label()

    stem = f"flavor_{sanitize_token(info['flavor'])}_m{info['mass']}_{sanitize_token(info['algo'])}_{sanitize_token(info['tag'])}"
    outbase = os.path.join(outdir, stem + "_fit")
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")

    result = {
        "flavor": info["flavor"],
        "mass": info["mass"],
        "algo": info["algo"],
        "tag": info["tag"],
        "fit_status": fitres.status(),
        "cov_qual": fitres.covQual(),
        "mean": mean.getVal(),
        "mean_err": mean.getError(),
        "sigmaL": sigma_left.getVal(),
        "sigmaL_err": sigma_left.getError(),
        "sigmaR": sigma_right.getVal(),
        "sigmaR_err": sigma_right.getError(),
        "alphaL": alpha_left.getVal(),
        "alphaL_err": alpha_left.getError(),
        "nL": n_left.getVal(),
        "nL_err": n_left.getError(),
        "alphaR": alpha_right.getVal(),
        "alphaR_err": alpha_right.getError(),
        "nR": n_right.getVal(),
        "nR_err": n_right.getError(),
        "sigma_eff": sigma_eff,
        "eff_low": x_eff_low,
        "eff_high": x_eff_high,
        "n_events_total": total_events,
        "n_events_68": n_eff,
        "eff_68_pct": eff_68,
        "fwhm": fwhm_info["fwhm"] if fwhm_info else None,
        "sigma_eq": sigma_eq,
        "rel_resolution_pct": rel_width_pct,
        "window_sigma_used": sigma_for_windows,
        "n_events_1sigma": n_1sig,
        "eff_1sigma_pct": eff_1sig,
        "n_events_3sigma": n_3sig,
        "eff_3sigma_pct": eff_3sig,
        "three_sigma_low": x3sig_low,
        "three_sigma_high": x3sig_high,
    }

    print(f"\n[INFO] {os.path.basename(info['path'])} ({info['flavor']})")
    print(f"       mean       = {result['mean']:.3f} +/- {result['mean_err']:.3f}")
    print(f"       sigmaL     = {result['sigmaL']:.3f} +/- {result['sigmaL_err']:.3f}")
    print(f"       sigmaR     = {result['sigmaR']:.3f} +/- {result['sigmaR_err']:.3f}")
    print(f"       fit status = {result['fit_status']}, covQual = {result['cov_qual']}")
    print(f"       sigma_eff  = {result['sigma_eff']:.3f} GeV")
    if result["fwhm"] is not None:
        print(f"       FWHM       = {result['fwhm']:.3f} GeV")
    print(f"       N(3sigma)  = {result['n_events_3sigma']:.1f} ({result['eff_3sigma_pct']:.2f}%)")

    return result


def split_pads_for_ratio(canvas):
    pad_top = ROOT.TPad("pad_top", "", 0.0, 0.30, 1.0, 1.0)
    pad_bot = ROOT.TPad("pad_bot", "", 0.0, 0.0, 1.0, 0.30)
    pad_top.SetBottomMargin(0.03)
    pad_top.SetLeftMargin(0.14)
    pad_top.SetRightMargin(0.05)
    pad_bot.SetTopMargin(0.02)
    pad_bot.SetBottomMargin(0.34)
    pad_bot.SetLeftMargin(0.14)
    pad_bot.SetRightMargin(0.05)
    pad_top.Draw()
    pad_bot.Draw()
    return pad_top, pad_bot


def draw_flavor_mass_overlay(group_key, hist_entries, outdir, args):
    if len(hist_entries) < 2:
        return

    tag, algo, mass = group_key
    draw_ratio = (not args.no_ratio)

    prepared = []
    for info, hist in hist_entries:
        clone = hist.Clone(f"cmp_{hist.GetName()}")
        clone.SetDirectory(0)
        if args.norm and clone.Integral() > 0:
            clone.Scale(1.0 / clone.Integral())
        color = FLAVOR_COLORS.get(info["flavor"], ROOT.kGray + 2)
        clone.SetLineColor(color)
        clone.SetMarkerColor(color)
        clone.SetLineWidth(3)
        prepared.append((info, clone))

    prepared.sort(key=lambda item: FLAVOR_ORDER.index(item[0]["flavor"]) if item[0]["flavor"] in FLAVOR_ORDER else 999)
    ref_entry = next((item for item in prepared if item[0]["flavor"] == args.reference_flavor), None)
    if ref_entry is None:
        draw_ratio = False

    canvas = ROOT.TCanvas(
        f"c_cmp_{sanitize_token(tag)}_{sanitize_token(algo)}_{mass}",
        "",
        1200,
        1000,
    )

    pad_top = canvas
    pad_bot = None
    drawable_refs = []
    if draw_ratio:
        pad_top, pad_bot = split_pads_for_ratio(canvas)

    pad_top.cd()
    draw_xmin = max(prepared[0][1].GetXaxis().GetXmin(), 0.25 * mass)
    draw_xmax = min(prepared[0][1].GetXaxis().GetXmax(), 2.0 * mass)

    ymax = 0.0
    for _, hist in prepared:
        bmin = hist.GetXaxis().FindBin(draw_xmin + 1e-6)
        bmax = hist.GetXaxis().FindBin(draw_xmax - 1e-6)
        local_max = max(hist.GetBinContent(i) for i in range(bmin, bmax + 1))
        ymax = max(ymax, local_max)

    ymin = 5e-8 if args.norm else 0.0
    ymax *= 1.35

    legend = ROOT.TLegend(0.62, 0.66, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.SetHeader(f"M={mass} GeV ({algo})", "C")

    first = True
    for info, hist in prepared:
        hist.SetTitle("")
        hist.GetXaxis().SetRangeUser(draw_xmin, draw_xmax)
        hist.SetMaximum(ymax)
        hist.SetMinimum(ymin)
        hist.GetXaxis().SetTitle("m_{jj} [GeV]")
        hist.GetYaxis().SetTitle("A.U." if args.norm else "Events / 5 GeV")
        hist.Draw("hist" if first else "hist same")
        legend.AddEntry(hist, info["flavor"], "l")
        first = False

    legend.Draw()
    draw_cms_label()

    if draw_ratio and pad_bot:
        pad_bot.cd()
        ref_hist = ref_entry[1]
        first_ratio = True
        for info, hist in prepared:
            if info["flavor"] == args.reference_flavor:
                continue
            ratio = hist.Clone(f"ratio_{hist.GetName()}")
            ratio.SetDirectory(0)
            ratio.Divide(ref_hist)
            ratio.SetTitle("")
            ratio.GetXaxis().SetRangeUser(draw_xmin, draw_xmax)
            ratio.GetYaxis().SetRangeUser(args.ratio_ymin, args.ratio_ymax)
            ratio.GetYaxis().SetTitle(f"/{args.reference_flavor}")
            ratio.GetXaxis().SetTitle("m_{jj} [GeV]")
            ratio.GetYaxis().SetNdivisions(505)
            ratio.GetYaxis().SetTitleSize(0.11)
            ratio.GetYaxis().SetTitleOffset(0.45)
            ratio.GetYaxis().SetLabelSize(0.09)
            ratio.GetXaxis().SetTitleSize(0.12)
            ratio.GetXaxis().SetLabelSize(0.10)
            ratio.Draw("hist" if first_ratio else "hist same")
            drawable_refs.append(ratio)
            first_ratio = False

        line = ROOT.TLine(draw_xmin, 1.0, draw_xmax, 1.0)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray + 2)
        line.Draw()
        drawable_refs.append(line)

    suffix = "_norm" if args.norm else ""
    suffix += "_ratio" if draw_ratio else ""
    outbase = os.path.join(
        outdir,
        f"flavorCompare_m{mass}_{sanitize_token(algo)}_{sanitize_token(tag)}{suffix}",
    )
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")


def save_summary_table(rows, outpath):
    columns = [
        "flavor",
        "mass",
        "algo",
        "tag",
        "fit_status",
        "cov_qual",
        "mean",
        "mean_err",
        "sigmaL",
        "sigmaL_err",
        "sigmaR",
        "sigmaR_err",
        "alphaL",
        "alphaL_err",
        "nL",
        "nL_err",
        "alphaR",
        "alphaR_err",
        "nR",
        "nR_err",
        "sigma_eff",
        "eff_low",
        "eff_high",
        "n_events_total",
        "n_events_68",
        "eff_68_pct",
        "fwhm",
        "sigma_eq",
        "rel_resolution_pct",
        "window_sigma_used",
        "n_events_1sigma",
        "eff_1sigma_pct",
        "three_sigma_low",
        "three_sigma_high",
        "n_events_3sigma",
        "eff_3sigma_pct",
    ]

    with open(outpath, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            serialized = {}
            for key in columns:
                val = row.get(key)
                if isinstance(val, float) and (math.isinf(val) or math.isnan(val)):
                    val = None
                serialized[key] = val
            writer.writerow(serialized)


def main():
    args = parse_args()
    set_style()
    os.makedirs(args.outdir, exist_ok=True)

    requested_flavors = [token.strip().lower() for token in args.flavors.split(",") if token.strip()]
    file_infos = discover_files(args.indir, args.recursive, args.algo)
    if not file_infos:
        raise RuntimeError(f"No matching ROOT files found in {args.indir}")

    enriched_infos = []
    for info in file_infos:
        flavor = infer_flavor("/".join([info["path"], info["tag"], info["algo"]]))
        if requested_flavors and flavor not in requested_flavors:
            continue
        record = dict(info)
        record["flavor"] = flavor
        enriched_infos.append(record)

    if not enriched_infos:
        raise RuntimeError("No files matched the requested flavor filter.")

    summary_rows = []
    overlay_groups = {}

    for info in enriched_infos:
        suffix = f"{info['flavor']}_{info['mass']}_{sanitize_token(info['algo'])}_{sanitize_token(info['tag'])}"
        hist = load_histogram(info["path"], args.hist, suffix)
        if not hist:
            continue

        fit_row = fit_flavor_mass_histogram(info, hist, args.outdir)
        summary_rows.append(fit_row)

        key = (info["tag"], info["algo"], info["mass"])
        overlay_groups.setdefault(key, []).append((info, hist))

    if not summary_rows:
        raise RuntimeError(f"No valid histograms named {args.hist} were processed.")

    for key, entries in overlay_groups.items():
        draw_flavor_mass_overlay(key, entries, args.outdir, args)

    summary_path = os.path.join(args.outdir, args.summary_name)
    save_summary_table(summary_rows, summary_path)
    print(f"\n[INFO] Wrote summary table: {summary_path}")
    print(f"[INFO] Done. Outputs saved in: {args.outdir}")


if __name__ == "__main__":
    main()
