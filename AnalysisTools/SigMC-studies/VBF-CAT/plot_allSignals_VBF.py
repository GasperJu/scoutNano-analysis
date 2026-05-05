#!/usr/bin/env python3

import argparse
import os
import re

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay VBF signal peak histograms for all available masses."
    )
    parser.add_argument(
        "--indir",
        type=str,
        default="VBF-dijetMass-Histos_ForFIT",
        help="Directory containing vbf-m<MASS>-<ALGO>.root files.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots_VBFSignalFits",
        help="Directory where the overlay plots will be saved.",
    )
    parser.add_argument(
        "--hist",
        type=str,
        default="h_signal_peak_selected_5GeV",
        help="Histogram name to overlay.",
    )
    parser.add_argument(
        "--algo",
        type=str,
        default=None,
        help="Only use files matching this algorithm token, e.g. leading or minDR.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Scan subdirectories of --indir as well.",
    )
    parser.add_argument(
        "--norm",
        action="store_true",
        help="Normalize each histogram to unit area before plotting.",
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="Use a logarithmic y axis.",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=None,
        help="Optional x-axis minimum.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Optional x-axis maximum.",
    )
    return parser.parse_args()


def hex_to_rootcolor(hexcode):
    hexcode = hexcode.lstrip("#")
    r = int(hexcode[0:2], 16) / 255.0
    g = int(hexcode[2:4], 16) / 255.0
    b = int(hexcode[4:6], 16) / 255.0
    return ROOT.TColor.GetColor(r, g, b)


COLORS = [
    hex_to_rootcolor("#3f90da"),
    hex_to_rootcolor("#ffa90e"),
    hex_to_rootcolor("#bd1f01"),
    hex_to_rootcolor("#832db6"),
    hex_to_rootcolor("#94a4a2"),
    hex_to_rootcolor("#e76300"),
    hex_to_rootcolor("#b9ac70"),
]


def set_style():
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
    ROOT.gStyle.SetTitleOffset(1.05, "X")
    ROOT.gStyle.SetTitleOffset(1.35, "Y")
    ROOT.gStyle.SetEndErrorSize(0)
    ROOT.gStyle.SetLegendFont(42)


def draw_cms_label(extra="2024 (13.6 TeV)"):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.132, 0.932, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.04)
    latex2.DrawLatex(0.21, 0.932, "Simulation Preliminary")

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(0.95, 0.932, extra)


def parse_signal_file(path, relpath):
    match = re.search(r"vbf-m(?P<mass>\d+)-(?P<algo>[^/]+)\.root$", relpath)
    if not match:
        return None
    parent_dir = os.path.dirname(relpath)
    return {
        "path": path,
        "mass": int(match.group("mass")),
        "algo": match.group("algo"),
        "tag": os.path.basename(parent_dir) if parent_dir not in ("", ".") else "base",
    }


def discover_files(indir, recursive, algo_filter):
    candidates = []
    if recursive:
        for root_dir, _, filenames in os.walk(indir):
            for filename in filenames:
                if not filename.endswith(".root"):
                    continue
                full_path = os.path.join(root_dir, filename)
                relpath = os.path.relpath(full_path, indir)
                info = parse_signal_file(full_path, relpath)
                if info is None:
                    continue
                if algo_filter and info["algo"] != algo_filter:
                    continue
                candidates.append(info)
    else:
        for filename in os.listdir(indir):
            if not filename.endswith(".root"):
                continue
            full_path = os.path.join(indir, filename)
            info = parse_signal_file(full_path, filename)
            if info is None:
                continue
            if algo_filter and info["algo"] != algo_filter:
                continue
            candidates.append(info)
    return sorted(candidates, key=lambda item: (item["tag"], item["algo"], item["mass"]))


def load_histograms(file_infos, hist_name, normalize):
    groups = {}
    for idx, info in enumerate(file_infos):
        root_file = ROOT.TFile.Open(info["path"])
        if not root_file or root_file.IsZombie():
            print(f"[WARNING] Could not open file: {info['path']}")
            continue

        hist = root_file.Get(hist_name)
        if not hist:
            print(f"[WARNING] Missing histogram {hist_name} in {info['path']}")
            root_file.Close()
            continue

        hist = hist.Clone(f"{hist.GetName()}_{info['tag']}_{info['algo']}_{info['mass']}")
        hist.SetDirectory(0)
        root_file.Close()

        if normalize and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

        color = COLORS[idx % len(COLORS)]
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetLineWidth(3)

        key = (info["tag"], info["algo"])
        groups.setdefault(key, []).append((info, hist))

    for key in groups:
        groups[key].sort(key=lambda pair: pair[0]["mass"])
    return groups


def compute_xrange(entries, args):
    xmin = args.xmin
    xmax = args.xmax
    if xmin is not None and xmax is not None:
        return xmin, xmax

    masses = [entry[0]["mass"] for entry in entries]
    if xmin is None:
        xmin = max(0.0, 0.35 * min(masses))
    if xmax is None:
        xmax = min(entries[0][1].GetXaxis().GetXmax(), 1.9 * max(masses))
    return xmin, xmax


def sanitize_token(token):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", token)


def draw_group(group_key, entries, outdir, hist_name, args):
    tag, algo = group_key
    canvas_name = f"c_{sanitize_token(tag)}_{sanitize_token(algo)}"
    canvas = ROOT.TCanvas(canvas_name, "", 1400, 1000)
    if args.logy:
        canvas.SetLogy()

    xmin, xmax = compute_xrange(entries, args)
    ymax = 0.0
    for _, hist in entries:
        bmin = hist.GetXaxis().FindBin(xmin + 1e-6)
        bmax = hist.GetXaxis().FindBin(xmax - 1e-6)
        local_max = max(hist.GetBinContent(bin_idx) for bin_idx in range(bmin, bmax + 1))
        ymax = max(ymax, local_max)

    ymax *= 1.35
    ymin = 5e-6 if args.logy and args.norm else (0.5 if args.logy else 0.0)

    legend = ROOT.TLegend(0.60, 0.62, 0.88, 0.88)
    legend.SetTextSize(0.033)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetHeader(f"VBFHTo2B")
    #legend.SetHeader(f"VBF signal ({algo})")

    first = True
    for info, hist in entries:
        hist.SetTitle("")
        hist.GetXaxis().SetRangeUser(xmin, xmax)
        hist.SetMaximum(ymax)
        hist.SetMinimum(ymin)
        hist.GetXaxis().SetTitle("m_{jj} [GeV]")
        hist.GetYaxis().SetTitle("A.U." if args.norm else "Events / 5 GeV")
        draw_opt = "hist" if first else "hist same"
        hist.Draw(draw_opt)
        legend.AddEntry(hist, f"M = {info['mass']} GeV", "l")
        first = False

    legend.Draw()

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.036)
    subdir_text = "mass-dependent selection" if tag == "base" else tag
    #label.DrawLatex(0.16, 0.84, f"{subdir_text}")
    #subdir_text = "latest mass-dependent selection" if tag == "base" else tag
    #label.DrawLatex(0.16, 0.84, f"{subdir_text}, hist: {hist_name}")

    draw_cms_label()

    suffix = "_norm" if args.norm else ""
    suffix += "_logy" if args.logy else ""
    base_name = f"allSignals_{sanitize_token(tag)}_{sanitize_token(algo)}{suffix}"
    outbase = os.path.join(outdir, base_name)
    canvas.SaveAs(outbase + ".png")
    canvas.SaveAs(outbase + ".pdf")
    print(f"[INFO] Saved overlay: {outbase}.(png|pdf)")


def main():
    args = parse_args()
    set_style()
    os.makedirs(args.outdir, exist_ok=True)

    file_infos = discover_files(args.indir, args.recursive, args.algo)
    if not file_infos:
        raise RuntimeError(f"No matching ROOT files found in {args.indir}")

    groups = load_histograms(file_infos, args.hist, args.norm)
    if not groups:
        raise RuntimeError(f"No histograms named {args.hist} could be loaded.")

    for group_key, entries in groups.items():
        draw_group(group_key, entries, args.outdir, args.hist, args)

    print(f"[INFO] Done. Outputs saved in: {args.outdir}")


if __name__ == "__main__":
    main()
