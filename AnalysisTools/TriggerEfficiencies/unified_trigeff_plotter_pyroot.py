#!/usr/bin/env python3
"""
unified_trigeff_plotter_pyroot.py

Framework-compatible sibling of unified_trigeff_plotter.py -- v2.

What changed from v1
---------------------
v1 copy-pasted compare_data_mc_eff.py's functions into this file as text.
That's still "the original logic," but it's a fork: if compare_data_mc_eff.py
is fixed or extended in the real repo, this script would silently drift out
of sync. v2 instead *imports the actual source files* (via importlib, since
their filenames -- e.g. "plotDijetHTEff_AN_Ortogonal-Mjj-junction.py" --
aren't valid Python module names for a normal `import`) and calls their real
functions directly. The original files are the single source of truth; this
script is a thin orchestration/auxiliary layer on top of them.

Two source files are imported:
  - compare_data_mc_eff.py            -> I/O + statistics primitives
                                          (get_sumw, load_hist, teff_to_arrays,
                                          ratio_asymm, parse_mc_spec, ...)
  - Mods Tested/Plot Mods/
      plotDijetHTEff_AN_Ortogonal-Mjj-junction.py
                                       -> create_ratio() -- a complete,
                                          correctly-implemented ratio
                                          function in the ORIGINAL file...
                                          that its own main() never actually
                                          calls. Checked with `grep -n
                                          "create_ratio("` on the source:
                                          one match, the definition itself.
                                          The panel main() actually renders
                                          just redraws both curves' raw
                                          values on a zoomed axis
                                          (ratio_print), not a real ratio --
                                          likely unfinished/orphaned code.
                                          This script resurrects and calls
                                          the working create_ratio() (proper
                                          asymmetric-error propagation)
                                          rather than replicating the
                                          apparent bug in what main() ships.

Closing the -junctionTRG.py gap
--------------------------------
v1's junction_ratio style hardcoded numerator suffix "passed" -- it couldn't
reproduce -junctionTRG.py's combination of a weighted background stack WITH
a multi-hypothesis L1-seed overlay. junction_ratio now accepts --hypotheses
too: with it set, one full stack-vs-primary-plus-ratio figure is produced
per hypothesis (looping the numerator suffix), instead of only "passed".

What's still new glue (not in any original script), clearly marked below:
  - import_original()                 (loads a .py file as a module by path)
  - graph_to_arrays() / point_to_xy()  (ROOT TGraph -> numpy, for mplhep)
  - find_target_eff_crossing()         (Step 4a item 5's new feature)
  - discover_variables()               (needed for --var all)
  - build_weighted_stack()             (glue: loops samples calling the
                                          IMPORTED get_sumw/load_hist/etc.,
                                          since neither original script
                                          factored its background-stacking
                                          loop into a standalone function)
  - the --style dispatch / CLI / file_cache

IMPORTANT: requires real PyROOT (CMSSW). Syntax- and import-graph-checked
(`python3 -m py_compile`, `pyflakes`) in a sandbox without ROOT -- NOT
execution-tested. Validate in the actual framework, ideally through Step 5
staging on red_list_data2024H.txt, before trusting its output.
"""
from __future__ import annotations

import os
import sys
import importlib.util
from argparse import ArgumentParser
from ctypes import c_double

import numpy as np
import matplotlib.pyplot as plt

import ROOT

try:
    import mplhep as hep
    HAVE_MPLHEP = True
except ImportError as exc:
    hep = None
    HAVE_MPLHEP = False
    print(f"[WARN] Could not import mplhep ({exc}). Falling back to plain matplotlib style.")

ROOT.gROOT.SetBatch(True)
if HAVE_MPLHEP:
    hep.style.use("CMS")

# Shared with unified_trigeff_plotter.py (the uproot backend): XLABEL_MAP,
# get_xlabel, and find_target_eff_crossing are backend-agnostic (pure dict
# lookups / numpy), so both scripts import the SAME definitions from here
# instead of each carrying its own copy that could drift out of sync.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common.efficiency import get_xlabel, find_target_eff_crossing, match_variable_pairs, parse_weight_spec, resolve_scale, PETROFF10


# =====================================================================
# NEW: load the original source files as modules and call their real code
# =====================================================================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_COMPARE_SCRIPT = os.path.join(_THIS_DIR, "compare_data_mc_eff.py")
DEFAULT_JUNCTION_SCRIPT = os.path.join(
    _THIS_DIR, "Mods Tested", "Plot Mods", "plotDijetHTEff_AN_Ortogonal-Mjj-junction.py"
)


def import_original(path, module_name):
    """
    Load a .py file as a live module by path (importlib, not `import`),
    because these scripts have filenames with hyphens/spaces that aren't
    valid Python identifiers. Their module-level code (ROOT.gROOT.SetBatch,
    gStyle tweaks) runs harmlessly on import; their CLI blocks are guarded
    by `if __name__ == "__main__":` in both source files, so importing them
    does not trigger their argparse/main().
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Cannot find original script '{path}'. Pass --compare-script / "
            f"--junction-script explicitly if this unity script isn't sitting "
            f"in AnalysisTools/TriggerEfficiencies/ alongside the originals."
        )
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_compare = None   # set in main(), after --compare-script is known
_junction = None  # set in main(), after --junction-script is known


# =====================================================================
# NEW (small): bridge ROOT TGraph objects to numpy, for mplhep plotting.
# Needed because create_ratio() (imported from the original junction
# script, below) returns a TGraphAsymmErrors, not arrays.
# =====================================================================
def graph_to_arrays(g):
    n = g.GetN()
    x = np.zeros(n); y = np.zeros(n)
    ylo = np.zeros(n); yhi = np.zeros(n)
    xv, yv = c_double(), c_double()
    for i in range(n):
        g.GetPoint(i, xv, yv)
        x[i] = xv.value
        y[i] = yv.value
        ylo[i] = g.GetErrorYlow(i)
        yhi[i] = g.GetErrorYhigh(i)
    return x, y, ylo, yhi


# =====================================================================
# NEW: auxiliary layer (Step 4/4a additions -- not in any original script)
# =====================================================================
# XLABEL_MAP, get_xlabel, find_target_eff_crossing: imported from
# common/efficiency.py above -- see that module for their definitions.


def get_root_file(path, file_cache):
    """Open each ROOT file at most once per invocation -- fixes the N*M
    redundant-reopen cost in the originals' var x trigger loops, where
    ROOT.TFile.Open is called fresh on every combination."""
    if path not in file_cache:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open ROOT file: {path}")
        file_cache[path] = f
    return file_cache[path]


def discover_variables(fdir):
    """PyROOT-backend wrapper around the shared match_variable_pairs():
    extract this TDirectory's key names, then apply the one shared matching
    rule (same rule the uproot backend uses via common/efficiency.py)."""
    key_names = {k.GetName() for k in fdir.GetListOfKeys()}
    return match_variable_pairs(key_names)


def add_cms_label(ax, year):
    """mplhep CMS styling applied uniformly to every style this script
    produces (per the last round of changes), regardless of whether the
    ORIGINAL script that inspired that style used mplhep or a native
    ROOT TCanvas."""
    if HAVE_MPLHEP:
        hep.cms.label("Preliminary", data=True, year=year, com="13.6", ax=ax)
        return
    ax.text(0.0, 1.01, "CMS Preliminary", transform=ax.transAxes,
             ha="left", va="bottom", fontsize=18, fontweight="bold")
    ax.text(1.0, 1.01, f"{year} (13.6 TeV)", transform=ax.transAxes,
             ha="right", va="bottom", fontsize=16)


def draw_target_eff_marker(ax_list, centers, eff, target_eff, plateau_n, color="black"):
    x_cross, plateau = find_target_eff_crossing(centers, eff, target_eff, plateau_n)
    if x_cross is None:
        print(f"[WARN] target-eff={target_eff} never reached (plateau={plateau:.3f}); no marker drawn.")
        return None
    for ax in ax_list:
        ax.axvline(x_cross, color=color, linestyle=":", linewidth=1.6)
    ax_list[0].text(x_cross, 0.05, f"{target_eff*100:.0f}% eff. @ {x_cross:.1f}",
                     rotation=90, va="bottom", ha="right", fontsize=13, color=color,
                     transform=ax_list[0].get_xaxis_transform())
    print(f"[INFO] target-eff={target_eff}: crossing at x={x_cross:.3f} (plateau={plateau:.4f})")
    return x_cross, plateau


PALETTE = PETROFF10  # real CMS petroff10, verified against -CMS.py's own source (see common/efficiency.py)


def resolve_edges(args):
    if args.edges:
        new_edges = np.array([float(x) for x in args.edges.split(",")], dtype=float)
        if len(new_edges) < 2:
            raise ValueError("--edges must contain at least two numbers")
        return None, new_edges
    return args.rebin, None


def resolve_bkg_items_and_mode(args):
    """
    Two ways to specify background weights, both still supported:
    the legacy --bkg-weights positional list (always 'raw' mode, for
    backward compatibility with existing commands), or inline
    'file.root:value' specs in --bkg itself, interpreted per
    --bkg-weight-mode.
    """
    if args.bkg_weights:
        if len(args.bkg_weights) != len(args.bkg):
            raise ValueError("--bkg-weights must have exactly one entry per --bkg file.")
        items = [f"{path}:{w}" for path, w in zip(args.bkg, args.bkg_weights)]
        return items, "raw"
    return list(args.bkg), args.bkg_weight_mode


def resolve_variables(args, file_cache):
    if args.var != ["all"]:
        return args.var
    ref_paths = []
    if args.style == "cms_ratio":
        raw = list(args.data) + list(args.mc)
        ref_paths = [parse_weight_spec(item)[0] for item in raw]
    elif args.style == "multicurve":
        ref_paths = list(args.inputs)
    elif args.style == "junction_ratio":
        bkg_items, _ = resolve_bkg_items_and_mode(args)
        raw = [args.primary] + bkg_items
        ref_paths = [parse_weight_spec(item)[0] for item in raw]

    found = set()
    for p in ref_paths:
        fpath = p if (os.path.isabs(p) or not args.mc_indir) else os.path.join(args.mc_indir, p)
        f = get_root_file(fpath, file_cache)
        fdir = f.Get(args.dir)
        if fdir:
            found.update(discover_variables(fdir))
    if not found:
        raise RuntimeError(f"--var all found no h_{{var}}_all/_passed pairs in directory '{args.dir}'.")
    return sorted(found)


# =====================================================================
# GLUE (new): loops samples calling the IMPORTED primitives underneath.
# Neither original script factored background-stacking into a standalone,
# reusable function -- -junction.py's main() does it inline -- so this is
# new orchestration code, but every histogram read/scale/add call inside
# it is the imported original (_compare.load_hist, .get_sumw, etc.), not
# a reimplementation.
# =====================================================================
def build_weighted_stack(items, directory_name, var, trg, file_cache, rebin=None, new_edges=None,
                          mode="raw", lumi_pb=1000.0):
    """
    items: list of 'file.root' or 'file.root:value' weight-spec strings
    (parsed via the SHARED parse_weight_spec/resolve_scale from
    common/efficiency.py -- same parser the uproot script uses).
    mode="raw" (default) matches -junction.py's actual convention: value
    used directly as Scale() multiplier, no sumw/lumi division. mode=
    "xsec_sumw" instead computes scale = xsec*lumi/sumw via _compare.get_sumw,
    for when the more physically correct normalization is wanted instead.
    """
    den_tot, num_tot = None, None
    for item in items:
        path, value = parse_weight_spec(item, mode=mode)
        f = get_root_file(path, file_cache)
        fdir = f.Get(directory_name)
        if not fdir:
            raise RuntimeError(f"Missing directory '{directory_name}' in {path}")
        label = os.path.basename(path).replace(".root", "")
        den = _compare.load_hist(fdir, f"h_{var}_all", label, rebin=rebin, new_edges=new_edges)
        num = _compare.load_hist(fdir, f"h_{var}_{trg}", label, rebin=rebin, new_edges=new_edges)

        sumw = _compare.get_sumw(fdir) if mode == "xsec_sumw" else None
        scale = resolve_scale(value, mode, sumw=sumw, lumi_pb=lumi_pb)
        den.Scale(scale)
        num.Scale(scale)

        if den_tot is None:
            den_tot = den.Clone("den_bkg_tot"); den_tot.SetDirectory(0)
        else:
            _compare.check_same_binning(den_tot, den, label=f"bkg den {var}"); den_tot.Add(den)
        if num_tot is None:
            num_tot = num.Clone("num_bkg_tot"); num_tot.SetDirectory(0)
        else:
            _compare.check_same_binning(num_tot, num, label=f"bkg num {var} {trg}"); num_tot.Add(num)

    if not ROOT.TEfficiency.CheckConsistency(num_tot, den_tot):
        raise RuntimeError(f"TEfficiency consistency failed for background stack, var={var}, trg={trg}")
    return num_tot, den_tot


def combine_cached(items, directory_name, var, trg, file_cache, rebin=None, new_edges=None, mode="none"):
    """
    New (small) glue for cms_ratio's --mc-weight-mode "none"/"raw" paths:
    unlike _compare.build_mc_totals (the real original, which always
    requires an xsec and always opens files uncached), this respects
    file_cache and allows weight=1 samples for a quick test without an
    xsec. mode="xsec_sumw" still routes through the real
    _compare.build_mc_totals/build_data_totals instead, to keep that path
    exactly the original code.
    """
    return build_weighted_stack(items, directory_name, var, trg, file_cache,
                                 rebin=rebin, new_edges=new_edges, mode=mode)


# =====================================================================
# Styles
# =====================================================================
def style_cms_ratio(args, var, file_cache):
    """
    Every histogram-handling call below is the imported compare_data_mc_eff.py
    code when --mc-weight-mode is left at its default (xsec_sumw) -- this
    function's own job is only to hand it the right arguments and route the
    resulting arrays into an mplhep figure. --mc-weight-mode none/raw routes
    through combine_cached() instead (new glue), trading strict fidelity to
    the original's mandatory-xsec, uncached TFile.Open for the ability to
    test with a plain path and for keeping the file_cache benefit -- the
    original build_mc_totals opens every file fresh regardless of
    file_cache, so only the default xsec_sumw path carries that known cost.
    """
    rebin, new_edges = resolve_edges(args)

    if args.mc_weight_mode == "xsec_sumw":
        mc_samples = _compare.parse_mc_spec(args.mc)
        num_data, den_data = _compare.build_data_totals(
            args.data, args.dir, var, args.trigger,
            rebin=rebin, new_edges=new_edges, verbose=args.verbose,
        )
        num_mc, den_mc = _compare.build_mc_totals(
            mc_samples, args.mc_indir, args.dir, var, args.trigger, args.lumi_pb,
            rebin=rebin, new_edges=new_edges, verbose=args.verbose,
        )
    else:
        num_data, den_data = combine_cached(args.data, args.dir, var, args.trigger, file_cache,
                                             rebin=rebin, new_edges=new_edges, mode="none")
        num_mc, den_mc = combine_cached(args.mc, args.dir, var, args.trigger, file_cache,
                                         rebin=rebin, new_edges=new_edges, mode=args.mc_weight_mode)

    x_d, y_d, dlo, dup = _compare.teff_to_arrays(num_data, den_data)
    x_m, y_m, mlo, mup = _compare.teff_to_arrays(num_mc, den_mc)

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.0], hspace=0.05)
    ax_eff = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax_eff)

    ax_eff.errorbar(x_d, y_d, yerr=[dlo, dup], fmt='o', ms=4, lw=2.3, capsize=2, label=args.data_label)
    ax_eff.errorbar(x_m, y_m, yerr=[mlo, mup], fmt='s', ms=4, lw=2.3, capsize=2, label=args.mc_label, color='orchid')

    x_common, idx_d, idx_m = np.intersect1d(x_d, x_m, return_indices=True)
    if len(x_common) == 0:
        raise RuntimeError("Data/MC have no common populated x bins for the ratio panel.")
    r, rlo, rup = _compare.ratio_asymm(y_d[idx_d], dlo[idx_d], dup[idx_d], y_m[idx_m], mlo[idx_m], mup[idx_m])
    ax_ratio.errorbar(x_common, r, yerr=[rlo, rup], fmt='o', ms=4, lw=2.0, capsize=2)
    ax_ratio.axhline(1.0, linestyle='--', linewidth=1.2, color='black')

    add_cms_label(ax_eff, args.year)
    ax_ratio.set_xlabel(get_xlabel(var, args.xlabel))
    ax_eff.set_ylabel("Trigger efficiency")
    ax_ratio.set_ylabel("Data/MC")
    ax_eff.set_ylim(0, 1.1)
    ax_ratio.set_ylim(args.ratio_ymin, args.ratio_ymax)
    ax_eff.legend(fontsize=16, loc="lower right")
    plt.setp(ax_eff.get_xticklabels(), visible=False)

    if args.vline is not None:
        for ax in (ax_eff, ax_ratio):
            ax.axvline(args.vline, color='gray', linestyle='--', linewidth=1.2)
    if args.target_eff is not None:
        draw_target_eff_marker([ax_eff, ax_ratio], x_d, y_d, args.target_eff, args.plateau_n)

    return fig


def style_multicurve(args, var, file_cache):
    """Covers plotDijetHTEff.py/-trigEffIncl.py's category overlay, -CMS.py's
    petroff+mplhep multi-sample overlay (same palette, mplhep now default
    everywhere), and -TRG.py's multi-hypothesis overlay (--hypotheses)."""
    rebin, new_edges = resolve_edges(args)
    fig, ax = plt.subplots(figsize=(9, 8))

    if args.hypotheses:
        combos = [(path, hyp, hyp if len(args.inputs) == 1 else f"{lbl}: {hyp}")
                  for path, lbl in zip(args.inputs, args.input_labels) for hyp in args.hypotheses]
    else:
        combos = [(path, "passed", lbl) for path, lbl in zip(args.inputs, args.input_labels)]

    curves = []
    for i, (path, trg, label) in enumerate(combos):
        f = get_root_file(path, file_cache)
        fdir = f.Get(args.dir)
        if not fdir:
            raise RuntimeError(f"Missing directory '{args.dir}' in {path}")
        den = _compare.load_hist(fdir, f"h_{var}_all", label, rebin=rebin, new_edges=new_edges)
        num = _compare.load_hist(fdir, f"h_{var}_{trg}", label, rebin=rebin, new_edges=new_edges)
        x, y, elo, eup = _compare.teff_to_arrays(num, den)
        ax.errorbar(x, y, yerr=[elo, eup], fmt='o', ms=5, lw=2.0, capsize=2,
                    color=PALETTE[i % len(PALETTE)], label=label)
        curves.append((x, y))

    add_cms_label(ax, args.year)
    ax.set_xlabel(get_xlabel(var, args.xlabel))
    ax.set_ylabel("Trigger efficiency")
    ax.set_ylim(-0.02, 1.1)
    ax.legend(fontsize=13, loc="lower right", frameon=False, title=args.legend_title or None)

    if args.vline is not None:
        ax.axvline(args.vline, color='navy', linestyle='--', linewidth=1.4)
    if args.target_eff is not None and curves:
        draw_target_eff_marker([ax], curves[0][0], curves[0][1], args.target_eff, args.plateau_n)

    return fig


def style_junction_ratio(args, var, file_cache, trg="passed"):
    """
    Covers -junction.py's weighted-background-stack + ratio, now reusing the
    IMPORTED create_ratio() from that original file for the ratio panel
    itself, instead of a separate reimplementation.

    create_ratio(graph1, graph2) assumes matching point-for-point ordering
    (it indexes both graphs by position i, not by x-value) -- that's an
    assumption in the ORIGINAL code, carried over here deliberately rather
    than silently patched. If the primary and background-stack curves don't
    have the same number of populated bins (e.g. very different statistics
    knock out different bins), this falls back to the intersect1d-based
    ratio_asymm approach instead, with a printed warning, rather than
    crashing or silently misaligning points.

    trg: numerator suffix, defaults to "passed"; looped externally when
    --hypotheses is set (closing the -junctionTRG.py gap -- v1 could not
    combine background-stacking with hypothesis overlay at all).
    """
    rebin, new_edges = resolve_edges(args)
    bkg_items, bkg_mode = resolve_bkg_items_and_mode(args)

    f = get_root_file(args.primary, file_cache)
    fdir = f.Get(args.dir)
    if not fdir:
        raise RuntimeError(f"Missing directory '{args.dir}' in {args.primary}")
    den_p = _compare.load_hist(fdir, f"h_{var}_all", "primary", rebin=rebin, new_edges=new_edges)
    num_p = _compare.load_hist(fdir, f"h_{var}_{trg}", "primary", rebin=rebin, new_edges=new_edges)
    eff_p = ROOT.TEfficiency(num_p, den_p)
    eff_p.SetStatisticOption(ROOT.TEfficiency.kFCP)
    xp, yp, plo, pup = _compare.teff_to_arrays(num_p, den_p)

    num_b, den_b = build_weighted_stack(bkg_items, args.dir, var, trg, file_cache,
                                         rebin=rebin, new_edges=new_edges, mode=bkg_mode, lumi_pb=args.lumi_pb)
    eff_b = ROOT.TEfficiency(num_b, den_b)
    eff_b.SetStatisticOption(ROOT.TEfficiency.kFCP)
    xb, yb, blo, bup = _compare.teff_to_arrays(num_b, den_b)

    fig = plt.figure(figsize=(9, 9))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.05)
    ax_eff = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax_eff)

    ax_eff.errorbar(xp, yp, yerr=[plo, pup], fmt='o', ms=5, lw=2.0, capsize=2,
                     color=PALETTE[0], label=args.primary_label)
    ax_eff.errorbar(xb, yb, yerr=[blo, bup], fmt='s', ms=5, lw=2.0, capsize=2,
                     color=PALETTE[3], label=args.bkg_label)

    graph_p = eff_p.CreateGraph()
    graph_b = eff_b.CreateGraph()
    if graph_p.GetN() == graph_b.GetN():
        # Original code path: literally the imported create_ratio().
        ratio_graph = _junction.create_ratio(graph_p, graph_b)
        x_common, r, rlo_signed_lo, rlo_signed_hi = graph_to_arrays(ratio_graph)
        # create_ratio's ratio is Bkg/Primary already (graph1=primary is y1,
        # graph2=bkg is y2 -> ratio=y1/y2 in the original source... note this
        # is Primary/Bkg, the OPPOSITE convention of v1's Bkg/Primary label;
        # kept faithful to the original rather than silently flipped)
        rlo, rup = rlo_signed_lo, rlo_signed_hi
        ax_ratio.set_ylabel("Eff Primary/Eff Bkg")  # matches create_ratio()'s own axis title
    else:
        print(f"[WARN] var={var}, trg={trg}: primary has {graph_p.GetN()} populated bins, "
              f"background stack has {graph_b.GetN()} -- create_ratio()'s point-by-point "
              f"assumption doesn't hold here. Falling back to intersect1d + ratio_asymm.")
        x_common, idx_p, idx_b = np.intersect1d(xp, xb, return_indices=True)
        if len(x_common) == 0:
            raise RuntimeError("No common populated bins between primary and background stack.")
        r, rlo, rup = _compare.ratio_asymm(yb[idx_b], blo[idx_b], bup[idx_b], yp[idx_p], plo[idx_p], pup[idx_p])
        ax_ratio.set_ylabel("Bkg/Primary")

    ax_ratio.errorbar(x_common, r, yerr=[rlo, rup], fmt='o', ms=4, lw=2.0, capsize=2, color="black")
    ax_ratio.axhline(1.0, linestyle='--', linewidth=1.2, color='gray')

    add_cms_label(ax_eff, args.year)
    ax_ratio.set_xlabel(get_xlabel(var, args.xlabel))
    ax_eff.set_ylabel("Trigger efficiency")
    ax_eff.set_ylim(0, 1.1)
    ax_ratio.set_ylim(args.ratio_ymin, args.ratio_ymax)
    ax_eff.legend(fontsize=14, loc="lower right", frameon=False)
    plt.setp(ax_eff.get_xticklabels(), visible=False)

    if args.target_eff is not None:
        draw_target_eff_marker([ax_eff, ax_ratio], xp, yp, args.target_eff, args.plateau_n)

    return fig


def run_junction_ratio(args, var, file_cache):
    """
    Dispatch wrapper: if --hypotheses is set, loop the numerator suffix and
    produce one figure per hypothesis (closing the -junctionTRG.py gap --
    background stack x hypothesis overlay together). Otherwise, one figure
    with the default "passed" suffix, same as v1 and as -junction.py itself.
    Returns a list of (tag_suffix, fig) so main() can save each separately.
    """
    if not args.hypotheses:
        return [("", style_junction_ratio(args, var, file_cache, trg="passed"))]
    out = []
    for hyp in args.hypotheses:
        fig = style_junction_ratio(args, var, file_cache, trg=hyp)
        out.append((f"_{hyp}", fig))
    return out


STYLES = {"cms_ratio": style_cms_ratio, "multicurve": style_multicurve}  # junction_ratio handled via run_junction_ratio


def build_parser():
    p = ArgumentParser(description="Framework-compatible unified trigger-efficiency plotter (real PyROOT backend)")
    p.add_argument("--style", choices=["cms_ratio", "multicurve", "junction_ratio"], required=True)
    p.add_argument("--dir", default="InclusiveTrigNanoAOD")
    p.add_argument("--var", nargs="+", default=["ht_inclusive"])
    p.add_argument("--rebin", type=int, default=1)
    p.add_argument("--edges", default="")
    p.add_argument("--year", default="2024")
    p.add_argument("--xlabel", default="")
    p.add_argument("--vline", type=float, default=None)
    p.add_argument("--target-eff", type=float, default=None)
    p.add_argument("--plateau-n", type=int, default=5)
    p.add_argument("--ratio-ymin", type=float, default=0.5)
    p.add_argument("--ratio-ymax", type=float, default=1.5)
    p.add_argument("--outdir", default="./plots")
    p.add_argument("--formats", nargs="+", default=[".png"])
    p.add_argument("-v", "--verbose", action="store_true")

    p.add_argument("--compare-script", default=DEFAULT_COMPARE_SCRIPT,
                    help="Path to the original compare_data_mc_eff.py to import.")
    p.add_argument("--junction-script", default=DEFAULT_JUNCTION_SCRIPT,
                    help="Path to the original -junction.py to import.")

    # cms_ratio
    p.add_argument("--data", nargs="+",
                    help="One or more data .root files (plain paths, no xsec). Combined by summing.")
    p.add_argument("--mc", nargs="+",
                    help="MC samples. Format depends on --mc-weight-mode: 'xsec_sumw' (default) "
                         "requires 'file.root:xsec_pb' for every sample; 'none'/'raw' accept plain "
                         "paths too.")
    p.add_argument("--mc-weight-mode", choices=["xsec_sumw", "none", "raw"], default="xsec_sumw",
                    help="'xsec_sumw' (default): routes through the REAL imported "
                         "compare_data_mc_eff.py build_mc_totals (mandatory ':xsec_pb', proper "
                         "xsec*lumi/sumw normalization) -- exact original behavior, but that "
                         "original function opens files uncached regardless of file_cache. "
                         "'none'/'raw': lighter new glue accepting plain paths (weight=1) or "
                         "'file.root:value' as a direct multiplier -- trades strict fidelity to the "
                         "original function for testability without an xsec, while keeping the "
                         "file-cache benefit.")
    p.add_argument("--mc-indir", default="", help="Optional directory prefix prepended to relative --mc paths.")
    p.add_argument("--lumi-pb", type=float, default=1000.0,
                    help="Integrated luminosity in pb^-1, used only in xsec_sumw weight mode.")
    p.add_argument("--trigger", default="passed")
    p.add_argument("--data-label", default="Data")
    p.add_argument("--mc-label", default="MC (weighted sum)")

    # multicurve
    p.add_argument("--inputs", nargs="+", default=[])
    p.add_argument("--input-labels", nargs="+", default=[])
    p.add_argument("--hypotheses", nargs="+", default=[])
    p.add_argument("--legend-title", default="")

    # junction_ratio
    p.add_argument("--primary")
    p.add_argument("--primary-label", default="Data")
    p.add_argument("--bkg", nargs="+", default=[],
                    help="Background samples: 'file.root' or 'file.root:value' (see --bkg-weight-mode). "
                         "If --bkg-weights is also given, it takes precedence (legacy positional form).")
    p.add_argument("--bkg-weights", nargs="+", type=float, default=[],
                    help="Legacy: one weight per --bkg file, positionally paired, always a direct "
                         "multiplier (equivalent to --bkg-weight-mode raw via inline ':value').")
    p.add_argument("--bkg-weight-mode", choices=["raw", "xsec_sumw", "none"], default="raw",
                    help="How to interpret ':value' in --bkg when --bkg-weights isn't used. Default "
                         "'raw' matches -junction.py's actual convention (xsec used directly as the "
                         "scale, verified via the original source -- no sumw/lumi division). Use "
                         "'xsec_sumw' for the more physically correct xsec*lumi/sumw normalization "
                         "instead, via the imported _compare.get_sumw().")
    p.add_argument("--bkg-label", default="Bkg (weighted sum)")

    return p


def main():
    global _compare, _junction
    args = build_parser().parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    _compare = import_original(args.compare_script, "compare_data_mc_eff")
    _junction = import_original(args.junction_script, "junction_original")

    if args.style == "multicurve" and not args.input_labels:
        args.input_labels = [os.path.basename(p) for p in args.inputs]

    file_cache = {}
    variables = resolve_variables(args, file_cache)
    print(f"[INFO] Plotting {len(variables)} variable(s): {variables}")

    n_ok, n_failed = 0, 0
    for var in variables:
        try:
            if args.style == "junction_ratio":
                results = run_junction_ratio(args, var, file_cache)
            else:
                results = [("", STYLES[args.style](args, var, file_cache))]
        except Exception as exc:
            print(f"[ERROR] style='{args.style}' var='{var}' failed: {exc}")
            n_failed += 1
            continue

        for suffix, fig in results:
            tag = f"{args.style}_{var}{suffix}"
            for fmt in args.formats:
                outname = os.path.join(args.outdir, f"TrigEff_{tag}{fmt}")
                fig.savefig(outname, dpi=200, bbox_inches="tight")
                print(f"Saved: {outname}")
            plt.close(fig)
        n_ok += 1

    for f in file_cache.values():
        f.Close()

    print(f"[INFO] Done: {n_ok} variable(s) succeeded, {n_failed} failed.")
    if n_failed and not n_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
