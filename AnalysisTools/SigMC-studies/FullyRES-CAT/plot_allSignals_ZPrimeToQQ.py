#!/usr/bin/env python3
import os
import ROOT
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------
parser = argparse.ArgumentParser(description="Overlay resolved dijet mass signal templates")

parser.add_argument(
    "--indir",
    type=str,
    default=".",
    help="Directory containing resolved_dijetMass_M*.root"
)
parser.add_argument(
    "--outdir",
    type=str,
    default="/eos/user/e/elfontan/www/dijetAnaRun3/SIGModelling/",
    help="Output directory"
)
parser.add_argument(
    "--norm",
    action="store_true",
    help="Normalize histograms to unit area"
)
parser.add_argument(
    "--logy",
    action="store_true",
    help="Use logarithmic y axis"
)
parser.add_argument(
    "--xmin",
    type=float,
    default=0.0,
    help="x-axis minimum"
)
parser.add_argument(
    "--xmax",
    type=float,
    default=1500.0,   # updated from 2000 -> 1500
    help="x-axis maximum"
)

args = parser.parse_args()

MASSES = [100, 200, 500, 1000]
#MASSES = [50, 100, 200, 500, 1000]

# ------------------------------------------------------------
# CMS-style, color-vision-friendly palette
# Inspired by CMS/CAT cmsstyle direction.
# ------------------------------------------------------------
def hex_to_rootcolor(hexcode):
    hexcode = hexcode.lstrip("#")
    r = int(hexcode[0:2], 16) / 255.0
    g = int(hexcode[2:4], 16) / 255.0
    b = int(hexcode[4:6], 16) / 255.0
    return ROOT.TColor.GetColor(r, g, b)

COLORS = { #94a4a2: grey-green; #3f90da: blue; #ffa90e: orange; #bd1f01: red; #832db6:purple        
    #50:   hex_to_rootcolor("#94a4a2"),  # grey-green                     
    100:  hex_to_rootcolor("#3f90da"),  # blue                          
    200:  hex_to_rootcolor("#bd1f01"),  # red                                       
    500:  hex_to_rootcolor("#ffa90e"),  # orange 
    1000: hex_to_rootcolor("#832db6"),  # purple      
}

LINESTYLES = {
    #50: 1,
    100: 1,
    200: 1,
    500: 1,
    1000: 1,
}

# ------------------------------------------------------------
# CMS style helpers
# ------------------------------------------------------------
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

def draw_cms_label():
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
    latex3.DrawLatex(0.95, 0.932, "2024 (13.6 TeV)")

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    set_style()
    os.makedirs(args.outdir, exist_ok=True)

    histos = {}

    # Read histograms
    for mass in MASSES:
        infile = os.path.join(args.indir, f"resolved_dijetMass_M{mass}.root")
        hname = f"dijet_m_{mass}"

        if not os.path.exists(infile):
            print(f"[WARNING] Missing file: {infile}")
            continue

        f = ROOT.TFile.Open(infile)
        if not f or f.IsZombie():
            print(f"[WARNING] Could not open file: {infile}")
            continue

        h = f.Get(hname)
        if not h:
            print(f"[WARNING] Missing histogram {hname} in {infile}")
            f.Close()
            continue

        h = h.Clone(f"{hname}_clone")
        h.SetDirectory(0)
        f.Close()

        if args.norm and h.Integral() > 0:
            h.Scale(1.0 / h.Integral())

        h.SetLineColor(COLORS[mass])
        h.SetMarkerColor(COLORS[mass])
        h.SetLineStyle(LINESTYLES[mass])
        h.SetLineWidth(3)

        histos[mass] = h
        print(f"[INFO] Loaded {hname} from {infile} | integral = {h.Integral():.3f}")

    if len(histos) == 0:
        raise RuntimeError("No histograms loaded.")

    # Determine y-range after restricting x-range
    ymax = 0.0
    for h in histos.values():
        bmin = h.GetXaxis().FindBin(args.xmin + 1e-6)
        bmax = h.GetXaxis().FindBin(args.xmax - 1e-6)
        local_max = max(h.GetBinContent(b) for b in range(bmin, bmax + 1))
        ymax = max(ymax, local_max)

    ymax *= 1.35
    ymin = 5e-6 if args.logy and args.norm else (0.5 if args.logy else 0.0)

    # Canvas
    c = ROOT.TCanvas("c_signals", "", 1500, 1000)
    if args.logy:
        c.SetLogy()

    leg = ROOT.TLegend(0.60, 0.62, 0.88, 0.87)
    leg.SetTextSize(0.033)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetHeader("ZPrimeToQQ signal")

    first = True
    for mass in MASSES:
        if mass not in histos:
            continue

        h = histos[mass]
        h.GetXaxis().SetRangeUser(args.xmin, args.xmax)
        h.SetMaximum(ymax)
        h.SetMinimum(ymin)
        h.GetXaxis().SetTitle("m_{jj} [GeV]")
        h.GetYaxis().SetTitle("A.U." if args.norm else "Events")
        h.GetXaxis().SetMoreLogLabels(True)
        h.GetXaxis().SetNoExponent(False)

        if first:
            h.SetTitle("")
            h.Draw("hist")
            first = False
        else:
            h.Draw("hist same")

        leg.AddEntry(h, f"M = {mass} GeV", "l")

    leg.Draw()
    draw_cms_label()

    suffix = "_norm" if args.norm else ""
    suffix += "_log" if args.logy else ""
    #suffix += f"_xmax{int(args.xmax)}"

    outbase = os.path.join(args.outdir, f"resolved_dijetMass_signals{suffix}")
    c.SaveAs(outbase + ".png")
    c.SaveAs(outbase + ".pdf")

    print(f"[INFO] Saved plots:")
    print(f"       {outbase}.png")
    print(f"       {outbase}.pdf")


if __name__ == "__main__":
    main()
