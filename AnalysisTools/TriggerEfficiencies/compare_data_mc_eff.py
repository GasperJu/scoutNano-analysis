# ******************************************** #
# Run 3 dijet scouting analysis:               #
# trigger efficiency studies (data-MC)         #
# ******************************************** #
# HowToRun:
# python3 compare_data_mc_eff.py --data /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_inclTrgEffData_2024H_wJECs.root --mc histos_TT4Q_JECs.root:762.1 histos_Wto2Q_HT-100to400_JECs.root:16100 histos_Wto2Q_HT-400to800_JECs.root:357.5 histos_Zto2Q-HT100to400_JECs.root:6290.0  histos_QCD-HT200to400_JECs.root:1951000.0 histos_QCD-HT400to600_JECs.root:96660.0 histos_QCD-HT600to800_JECs.root:13684.0 histos_QCD-HT800to1000_JECs.root:3047.0 --mc-indir /eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff --lumi-pb 1000 --rebin 2 --variables ht_inclusive --vline 280

#SAMPLES = {
    #"histos_QCD-HT100to200_JECs.root": 25360000.0,
    #"histos_QCD-HT200to400_JECs.root": 1951000.0,
    #"histos_QCD-HT400to600_JECs.root": 96660.0,
    #"histos_QCD-HT600to800_JECs.root": 13684.0,
    #"histos_QCD-HT800to1000_JECs.root": 3047.0,
    #"histos_TT4Q_JECs.root": 762.1,
    #"histos_Wto2Q_HT-100to400_JECs.root": 16100.0,
    #"histos_Wto2Q_HT-400to800_JECs.root": 357.5,
    #"histos_Zto2Q-HT100to400_JECs.root": 6290.0,
#}


#!/usr/bin/env python3
import os
from argparse import ArgumentParser

import ROOT
import numpy as np
import matplotlib.pyplot as plt

try:
    import mplhep as hep
except ImportError as exc:
    hep = None
    print(f"[WARN] Could not import mplhep ({exc}). Falling back to plain matplotlib style.")

ROOT.gROOT.SetBatch(True)
if hep is not None:
    hep.style.use("CMS")

def add_cms_label(ax, year):
    if hep is not None:
        hep.cms.label("Preliminary", data=True, year=year, com="13.6", ax=ax)
        return

    ax.text(0.0, 1.01, "CMS Preliminary", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=18, fontweight="bold")
    ax.text(1.0, 1.01, f"{year} (13.6 TeV)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=16)

# -----------------------------
# Helpers (ROOT -> numpy)
# -----------------------------
def get_sumw(fdir, hname="h_sumw"):
    h = fdir.Get(hname)
    if not h:
        raise RuntimeError(f"Missing {hname} in directory {fdir.GetName()}")
    sw = h.Integral()
    if sw <= 0:
        raise RuntimeError(f"Invalid sumw={sw} from {hname} in {fdir.GetName()}")
    return sw

def clone_detach(h, name):
    out = h.Clone(name)
    out.SetDirectory(0)
    out.Sumw2()
    return out

def apply_rebin(h, rebin=None, new_edges=None):
    """
    Either:
      - rebin = integer factor (TH1::Rebin)
      - new_edges = list/np.array of bin edges (variable rebin)
    """
    if new_edges is not None:
        # ROOT expects an array('d') or python list is usually ok via PyROOT
        nb = len(new_edges) - 1
        h2 = h.Rebin(nb, f"{h.GetName()}_rebin", new_edges)
        h2.SetDirectory(0)
        h2.Sumw2()
        return h2

    if rebin is not None and int(rebin) > 1:
        h.Rebin(int(rebin))
    return h

def load_hist(fdir, hname, label, rebin=None, new_edges=None):
    hin = fdir.Get(hname)
    if not hin:
        raise RuntimeError(f"Missing {hname} in {label}")
    h = clone_detach(hin, f"{hname}_{label}")
    h = apply_rebin(h, rebin=rebin, new_edges=new_edges)
    return h

def check_same_binning(href, hnew, label=""):
    if href.GetNbinsX() != hnew.GetNbinsX():
        raise RuntimeError(f"[{label}] nbins mismatch: {href.GetNbinsX()} vs {hnew.GetNbinsX()}")
    ax1, ax2 = href.GetXaxis(), hnew.GetXaxis()
    if ax1.GetXmin() != ax2.GetXmin() or ax1.GetXmax() != ax2.GetXmax():
        raise RuntimeError(f"[{label}] x-range mismatch: [{ax1.GetXmin()},{ax1.GetXmax()}] vs [{ax2.GetXmin()},{ax2.GetXmax()}]")
    # catches variable binning differences
    for i in range(1, href.GetNbinsX() + 2):
        if ax1.GetBinLowEdge(i) != ax2.GetBinLowEdge(i):
            raise RuntimeError(f"[{label}] bin edge mismatch at edge {i}")

def teff_to_arrays(num, den, stat_opt=ROOT.TEfficiency.kFCP):
    """
    Build TEfficiency and return arrays:
    centers, eff, err_low, err_up
    """
    eff_obj = ROOT.TEfficiency(num, den)
    eff_obj.SetStatisticOption(stat_opt)

    centers, eff, elow, eup = [], [], [], []
    nb = den.GetNbinsX()

    for b in range(1, nb + 1):
        total = den.GetBinContent(b)
        if total <= 0:
            continue
        centers.append(den.GetBinCenter(b))
        eff.append(eff_obj.GetEfficiency(b))
        elow.append(eff_obj.GetEfficiencyErrorLow(b))
        eup.append(eff_obj.GetEfficiencyErrorUp(b))

    return np.array(centers), np.array(eff), np.array(elow), np.array(eup)

def ratio_asymm(data_eff, data_lo, data_up, mc_eff, mc_lo, mc_up, eps=1e-12):
    """
    Asymmetric ratio uncertainty (simple propagation).
    For r = d/m:
      r_low uses d_low and m_up
      r_up  uses d_up  and m_low
    """
    d = np.clip(data_eff, eps, None)
    m = np.clip(mc_eff,   eps, None)

    r = d / m

    # protect division by zero in relative errors
    rel_d_lo = data_lo / np.clip(d, eps, None)
    rel_d_up = data_up / np.clip(d, eps, None)
    rel_m_lo = mc_lo   / np.clip(m, eps, None)
    rel_m_up = mc_up   / np.clip(m, eps, None)

    r_lo = r * np.sqrt(rel_d_lo**2 + rel_m_up**2)
    r_up = r * np.sqrt(rel_d_up**2 + rel_m_lo**2)

    return r, r_lo, r_up

# -----------------------------
# MC combination logic
# -----------------------------
def parse_mc_spec(items):
    """
    Parse --mc like:
      /path/histos_TT.root:762.1 /path/histos_W.root:16100
    Returns list of tuples: (path, xsec_pb, label)
    """
    out = []
    for it in items:
        if ":" not in it:
            raise ValueError(f"MC spec must be 'file.root:xsec_pb' but got '{it}'")
        path, xs = it.rsplit(":", 1)
        xs = float(xs)
        label = os.path.basename(path).replace(".root", "")
        out.append((path, xs, label))
    return out

def build_mc_totals(mc_samples, indir, directory_name, var, trg, lumi_pb, rebin=None, new_edges=None, verbose=False):
    den_tot = None
    num_tot = None

    for (path, xs, label) in mc_samples:
        fpath = path
        if indir and not os.path.isabs(path):
            fpath = os.path.join(indir, path)

        f = ROOT.TFile.Open(fpath)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open MC file: {fpath}")

        fdir = f.Get(directory_name)
        if not fdir:
            raise RuntimeError(f"Missing directory '{directory_name}' in {fpath}")

        sumw = get_sumw(fdir)
        scale = xs * lumi_pb / sumw

        # Load raw (genWeight-filled) hists, then rebin, then scale
        den = load_hist(fdir, f"h_{var}_all", label, rebin=rebin, new_edges=new_edges)
        num = load_hist(fdir, f"h_{var}_{trg}", label, rebin=rebin, new_edges=new_edges)

        den.Scale(scale)
        num.Scale(scale)

        if verbose:
            print(f"[MC] {label}: xs={xs} pb, sumw={sumw:.6g}, scale={scale:.6g}, "
                  f"denInt={den.Integral():.6g}, numInt={num.Integral():.6g}")

        if den_tot is None:
            den_tot = den.Clone("den_tot")
            den_tot.SetDirectory(0)
        else:
            check_same_binning(den_tot, den, label=f"MC den {var}")
            den_tot.Add(den)

        if num_tot is None:
            num_tot = num.Clone("num_tot")
            num_tot.SetDirectory(0)
        else:
            check_same_binning(num_tot, num, label=f"MC num {var} {trg}")
            num_tot.Add(num)

        f.Close()

    if not ROOT.TEfficiency.CheckConsistency(num_tot, den_tot):
        # You can optionally print offending bins here if needed.
        raise RuntimeError(f"TEfficiency consistency failed for MC totals (num>den somewhere?) for var={var}, trg={trg}")

    return num_tot, den_tot

# -----------------------------
# Data loading
# -----------------------------
def build_data_totals(data_files, directory_name, var, trg, rebin=None, new_edges=None, verbose=False):
    """
    If multiple data files are provided, just add them (no scaling).
    """
    den_tot = None
    num_tot = None

    for fpath in data_files:
        f = ROOT.TFile.Open(fpath)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open Data file: {fpath}")

        fdir = f.Get(directory_name)
        if not fdir:
            raise RuntimeError(f"Missing directory '{directory_name}' in {fpath}")

        label = os.path.basename(fpath).replace(".root", "")

        den = load_hist(fdir, f"h_{var}_all", label, rebin=rebin, new_edges=new_edges)
        num = load_hist(fdir, f"h_{var}_{trg}", label, rebin=rebin, new_edges=new_edges)

        if verbose:
            print(f"[Data] {label}: denInt={den.Integral():.6g}, numInt={num.Integral():.6g}")

        if den_tot is None:
            den_tot = den.Clone("den_data_tot")
            den_tot.SetDirectory(0)
        else:
            check_same_binning(den_tot, den, label=f"Data den {var}")
            den_tot.Add(den)

        if num_tot is None:
            num_tot = num.Clone("num_data_tot")
            num_tot.SetDirectory(0)
        else:
            check_same_binning(num_tot, num, label=f"Data num {var} {trg}")
            num_tot.Add(num)

        f.Close()

    if not ROOT.TEfficiency.CheckConsistency(num_tot, den_tot):
        raise RuntimeError(f"TEfficiency consistency failed for Data totals (num>den somewhere?) for var={var}, trg={trg}")

    return num_tot, den_tot

# -----------------------------
# Plotting
# -----------------------------
def plot_eff_and_ratio(ax_eff, ax_ratio, x_d, y_d, dlo, dup, x_m, y_m, mlo, mup, label_data, label_mc):
    # efficiencies
    ax_eff.errorbar(x_d, y_d, yerr=[dlo, dup], fmt='o', ms=4, lw=2.3, capsize=2, label=label_data)
    ax_eff.errorbar(x_m, y_m, yerr=[mlo, mup], fmt='s', ms=4, lw=2.3, capsize=2, label=label_mc, color='orchid')

    # Ratio only where both samples have a populated denominator bin.
    x_common, idx_d, idx_m = np.intersect1d(x_d, x_m, return_indices=True)
    if len(x_common) == 0:
        raise RuntimeError("Data/MC have no common populated x bins for the ratio panel.")
    if len(x_common) != len(x_d) or len(x_common) != len(x_m):
        print(f"[WARN] Ratio uses {len(x_common)} common populated bins "
              f"(data={len(x_d)}, MC={len(x_m)}).")

    r, rlo, rup = ratio_asymm(
        y_d[idx_d], dlo[idx_d], dup[idx_d],
        y_m[idx_m], mlo[idx_m], mup[idx_m],
    )
    ax_ratio.errorbar(x_common, r, yerr=[rlo, rup], fmt='o', ms=4, lw=2.0, capsize=2)
    ax_ratio.axhline(1.0, linestyle='--', linewidth=1.2, color='black')

def main(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    directory_name = args.dir

    # Variables / triggers
    variables = args.variables
    triggers = args.triggers

    # Parse MC samples
    mc_samples = parse_mc_spec(args.mc)
    print(mc_samples)

    # Rebin configuration:
    # - if --edges is provided, it overrides --rebin
    new_edges = None
    if args.edges:
        # edges provided as comma-separated floats
        new_edges = [float(x) for x in args.edges.split(",")]
        if len(new_edges) < 2:
            raise ValueError("--edges must contain at least two numbers")
    rebin = args.rebin if new_edges is None else None

    for var in variables:
        for trg in triggers:
            # Build totals
            num_data, den_data = build_data_totals(
                data_files=args.data,
                directory_name=directory_name,
                var=var,
                trg=trg,
                rebin=rebin,
                new_edges=new_edges,
                verbose=args.verbose
            )

            num_mc, den_mc = build_mc_totals(
                mc_samples=mc_samples,
                indir=args.mc_indir,
                directory_name=directory_name,
                var=var,
                trg=trg,
                lumi_pb=args.lumi_pb,
                rebin=rebin,
                new_edges=new_edges,
                verbose=args.verbose
            )

            # Convert to arrays (Clopper-Pearson by default)
            x_d, y_d, dlo, dup = teff_to_arrays(num_data, den_data, stat_opt=ROOT.TEfficiency.kFCP)
            x_m, y_m, mlo, mup = teff_to_arrays(num_mc,   den_mc,   stat_opt=ROOT.TEfficiency.kFCP)

            # Figure with ratio panel
            fig = plt.figure(figsize=(10, 10))
            gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.0], hspace=0.05)
            ax_eff = fig.add_subplot(gs[0])
            ax_ratio = fig.add_subplot(gs[1], sharex=ax_eff)

            plot_eff_and_ratio(
                ax_eff, ax_ratio,
                x_d, y_d, dlo, dup,
                x_m, y_m, mlo, mup,
                label_data=args.data_label,
                label_mc=args.mc_label
            )

            # Labels / style
            add_cms_label(ax_eff, args.year)

            xlabel_map = {
                "ht_inclusive": r"$H_{T}$ [GeV]",
                "pt_leading": r"AK4 leading jet $p_{T}$ [GeV]",
                "dijetMass": r"Dijet mass [GeV]",
            }
            ax_ratio.set_xlabel(xlabel_map.get(var, var))
            ax_eff.set_ylabel("Trigger efficiency")
            ax_ratio.set_ylabel("Data/MC")
            ax_eff.set_ylim(0, 1.1)
            ax_ratio.set_ylim(args.ratio_ymin, args.ratio_ymax)

            ax_eff.legend(fontsize=24, loc="lower right")

            # Optional vertical line (your HT threshold)
            if args.vline is not None:
                ax_eff.axvline(args.vline, color='gray', linestyle='--', linewidth=1.2)
                ax_ratio.axvline(args.vline, color='gray', linestyle='--', linewidth=1.2)

                # text label near the line (top of eff panel)
                ax_eff.text(
                    args.vline * 1.015, 0.3,              # x in data coords, y in axis fraction
                    "L1 threshold",
                    rotation=90,
                    va="top",
                    ha="left",
                    fontsize=18,
                    color='gray',
                    transform=ax_eff.get_xaxis_transform()
                )

            # Remove x tick labels on top panel
            plt.setp(ax_eff.get_xticklabels(), visible=False)

            # Save
            tag = f"{var}_{trg}"
            if rebin is not None:
                tag += f"_rebin{rebin}"
            if new_edges is not None:
                tag += f"_varbins{len(new_edges)-1}"

            for fmt in args.formats:
                outname = os.path.join(outdir, f"TrigEff_DataVsMC_{tag}{fmt}")
                fig.savefig(outname, dpi=300)
                print(f"Saved: {outname}")

            plt.close(fig)

if __name__ == "__main__":
    parser = ArgumentParser(description="Data vs combined-MC trigger efficiency with ratio (mplhep)")

    # Inputs
    parser.add_argument("--data", nargs="+", required=True,
                        help="Data ROOT file(s). If multiple, they will be added.")
    parser.add_argument("--mc", nargs="+", required=True,
                        help="MC samples as 'file.root:xsec_pb' (multiple allowed).")
    parser.add_argument("--mc-indir", default="",
                        help="Optional base directory for MC files (used if MC paths are relative).")

    # ROOT directory inside files
    parser.add_argument("--dir", default="InclusiveTrigNanoAOD",
                        help="TDirectory name inside ROOT files (default: InclusiveTrigNanoAOD)")

    # Physics config
    parser.add_argument("--lumi-pb", type=float, default=1000.0,
                        help="Luminosity in pb^-1 for MC normalization (default: 1000)")
    parser.add_argument("--year", default="2024", help="Year label for CMS text")

    parser.add_argument("--variables", nargs="+", default=["ht_inclusive", "pt_leading"],
                        help="Variables to plot (default: ht_inclusive pt_leading)")
    parser.add_argument("--triggers", nargs="+", default=["passed"],
                        help="Trigger numerator suffixes (default: passed)")

    # Rebinning
    parser.add_argument("--rebin", type=int, default=1,
                        help="Integer rebin factor (applies to both data and MC). Ignored if --edges is set.")
    parser.add_argument("--edges", default="",
                        help="Comma-separated variable bin edges, e.g. '0,50,100,150,200'. Overrides --rebin.")

    # Plot config
    parser.add_argument("--outdir", default="/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/INCLUSIVE_TrgEff/",
                        help="Output directory")
    parser.add_argument("--formats", nargs="+", default=[".png", ".pdf"], help="Output formats")

    parser.add_argument("--data-label", default="Data", help="Legend label for data")
    parser.add_argument("--mc-label", default="MC (weighted sum)", help="Legend label for combined MC")

    parser.add_argument("--vline", type=float, default=280.0,
                        help="Optional vertical line x-position (set to -1 to disable)")
    parser.add_argument("--ratio-ymin", type=float, default=0.7)
    parser.add_argument("--ratio-ymax", type=float, default=1.3)

    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()

    if args.vline is not None and args.vline < 0:
        args.vline = None
    if args.edges.strip() == "":
        args.edges = ""

    main(args)
