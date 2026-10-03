#!/usr/bin/env python3
"""
unified_trigeff_plotter.py

One entry point that can reproduce the distinct plotting *aesthetics* found
across the mainstream repo and the branch's "Mods Tested/Plot Mods", so the
same input histograms can be rendered either way for a side-by-side visual
recheck, instead of maintaining N separate scripts.

Catalogued source styles -> --style value
------------------------------------------
  "cms_ratio"     <- compare_data_mc_eff.py
                     mplhep CMS style, data vs MC overlay + Data/MC ratio panel,
                     Clopper-Pearson errors, optional fixed --vline.

  "multicurve"    <- plotDijetHTEff.py / plotDijet_trigEffIncl.py / plotDijet_mcTrigEffIncl.py
                     Native-ROOT-canvas look (no ratio panel): several turn-on
                     curves overlaid on one axis, CMS/Preliminary text corners,
                     legend box -- reproduced here in matplotlib. This is the
                     "categorized" style: one curve per category
                     (VBF / Boosted / Resolved / Rest) or per background sample.

  "junction_ratio" <- Mods Tested/Plot Mods/plotDijetHTEff_AN_Ortogonal-Mjj-junction.py
                     Background-stack style: a weighted sum of several
                     background samples drawn as one combined curve alongside
                     a primary (data or signal-MC) curve, with a TPad-style
                     ratio panel underneath (background/primary).

New, style-independent addition (Step 4a, item 5)
--------------------------------------------------
  --target-eff (default 0.90): draws a *computed* vertical marker at the
  x-position where the curve first crosses target_eff * plateau (plateau =
  mean efficiency of the last N bins), annotated with the crossing value.
  This replaces guessing a fixed --vline with a number read directly off
  the curve, for every style.

Notes on scope
--------------
This intentionally does not chase pixel-parity with the original ROOT-canvas
scripts (exact font sizes, margins, etc.) -- the goal is a single, testable
code path that reproduces each style's *structure* (what's plotted, how
panels relate, what's overlaid) so it can be validated against synthetic
data first, then pointed at real histogram files.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    HAVE_MPLHEP = True
except ImportError:
    HAVE_MPLHEP = False

try:
    import uproot
except ImportError:
    print("[FATAL] uproot is required (pip install uproot).", file=sys.stderr)
    raise

from common.efficiency import (
    read_th1, rebin_counts, bin_centers,
    teff_to_arrays, check_consistency, ratio_asymm,
    find_target_eff_crossing, get_xlabel, discover_variables,
    parse_weight_spec, resolve_scale, get_sumw, PETROFF10,
)

PALETTE = PETROFF10  # real CMS petroff10, verified against -CMS.py's own source (see common/efficiency.py)


# -----------------------------------------------------------------------
# Shared data loading
# -----------------------------------------------------------------------
def get_file(file_path, file_cache):
    """
    Open a ROOT file with uproot exactly once per invocation, no matter how
    many variables/styles use it. This is the fix for the N*M redundant
    reopen cost confirmed in compare_data_mc_eff.py's original variable x
    trigger loop (each combination reopened every input file from scratch).
    """
    if file_path not in file_cache:
        file_cache[file_path] = uproot.open(file_path)
    return file_cache[file_path]


def load_efficiency_curve(file_path, dirname, var, num_suffix, rebin, file_cache):
    """Build (centers, eff, elo, eup) from an already-open (or now-cached) file."""
    f = get_file(file_path, file_cache)
    edges, den = read_th1(f, dirname, f"h_{var}_all")
    _, num = read_th1(f, dirname, f"h_{var}_{num_suffix}")
    edges, den = rebin_counts(edges, den, rebin)
    _, num = rebin_counts(edges, num, rebin)
    check_consistency(num, den, label=f"{file_path}:{var}_{num_suffix}")
    centers = bin_centers(edges)
    return teff_to_arrays(num, den, centers)


def resolve_variables(args, file_cache):
    """
    Expand --var into a concrete variable list:
      - explicit list of names -> used as-is
      - ["all"] -> auto-discover every h_{var}_all/_passed pair present in
        the style's reference file(s), so "plot everything in this file"
        doesn't require a hardcoded variable list that may not match what
        this particular file actually contains.
    """
    if args.var != ["all"]:
        return args.var

    ref_paths = []
    if args.style == "cms_ratio":
        # strip any 'file.root:value' weight-spec suffix before opening
        raw = list(args.data) + list(args.mc)
        ref_paths = [parse_weight_spec(item)[0] for item in raw]
    elif args.style == "multicurve":
        ref_paths = list(args.inputs)
    elif args.style == "junction_ratio":
        raw = [args.primary] + list(args.bkg)
        ref_paths = [parse_weight_spec(item)[0] for item in raw]

    found = set()
    for p in ref_paths:
        if not p:
            continue
        f = get_file(p, file_cache)
        found.update(discover_variables(f, args.dir))
    if not found:
        raise RuntimeError(
            "--var all requested, but no h_{var}_all/_passed pairs were found "
            f"in directory '{args.dir}' of the input file(s)."
        )
    return sorted(found)


def add_cms_label(ax, year, use_mplhep=True):
    """
    mplhep CMS styling is now the default for every style produced by this
    plotter (matching -CMS.py's convention), not just cms_ratio. The plain
    fallback only fires if mplhep genuinely isn't importable in this
    environment, or if a caller explicitly opts out.
    """
    if use_mplhep and HAVE_MPLHEP:
        hep.style.use("CMS")
        hep.cms.label("Preliminary", data=True, year=year, com="13.6", ax=ax)
        return
    ax.text(0.02, 0.96, "CMS", transform=ax.transAxes, ha="left", va="top",
             fontsize=18, fontweight="bold")
    ax.text(0.02, 0.90, "Preliminary", transform=ax.transAxes, ha="left", va="top",
             fontsize=14, style="italic")
    ax.text(0.98, 0.96, f"{year} (13.6 TeV)", transform=ax.transAxes,
             ha="right", va="top", fontsize=14)


def draw_target_eff_marker(ax_list, centers, eff, target_eff, plateau_n, color="black"):
    """Style-independent Step 4a addition: computed crossing marker + annotation."""
    x_cross, plateau = find_target_eff_crossing(centers, eff, target_eff=target_eff, plateau_n=plateau_n)
    if x_cross is None:
        print(f"[WARN] target-eff={target_eff} never reached (plateau={plateau:.3f}); no marker drawn.")
        return None
    for ax in ax_list:
        ax.axvline(x_cross, color=color, linestyle=":", linewidth=1.6)
    ax_list[0].text(
        x_cross, 0.05, f"{target_eff*100:.0f}% eff. @ {x_cross:.1f}",
        rotation=90, va="bottom", ha="right", fontsize=13, color=color,
        transform=ax_list[0].get_xaxis_transform(),
    )
    print(f"[INFO] target-eff={target_eff}: crossing at x={x_cross:.3f} (plateau={plateau:.4f})")
    return x_cross, plateau


def combine_curve(items, dirname, var, num_suffix, rebin, file_cache, mode="none", lumi_pb=1000.0):
    """
    Combine one or more weight-spec'd samples ('file.root' or
    'file.root:value', see common/efficiency.parse_weight_spec) into one
    (num, den, edges) triple, honoring whichever weighting convention this
    style/flag combination calls for (none / raw / xsec_sumw). Used for
    --data (always mode="none" -- none of the originals weighted data),
    --mc (mode configurable, default "xsec_sumw" matching
    compare_data_mc_eff.py), and --bkg (mode configurable, default "raw"
    matching -junction.py's actual convention).
    """
    den_stack, num_stack, edges_ref = None, None, None
    for item in items:
        path, value = parse_weight_spec(item, mode=mode)
        f = get_file(path, file_cache)
        edges, den = read_th1(f, dirname, f"h_{var}_all")
        _, num = read_th1(f, dirname, f"h_{var}_{num_suffix}")
        edges, den = rebin_counts(edges, den, rebin)
        _, num = rebin_counts(edges, num, rebin)

        sumw = get_sumw(f, dirname) if mode == "xsec_sumw" else None
        scale = resolve_scale(value, mode, sumw=sumw, lumi_pb=lumi_pb)
        den, num = den * scale, num * scale

        if den_stack is None:
            den_stack, num_stack, edges_ref = den, num, edges
        else:
            den_stack += den
            num_stack += num
    check_consistency(num_stack, den_stack, label=f"{dirname}/{var}_{num_suffix}")
    return num_stack, den_stack, edges_ref


# -----------------------------------------------------------------------
# Style: cms_ratio  (compare_data_mc_eff.py aesthetic)
# -----------------------------------------------------------------------
def style_cms_ratio(args, var, file_cache):
    num_d, den_d, edges_d = combine_curve(args.data, args.dir, var, "passed", args.rebin, file_cache, mode="none")
    xd, yd, dlo, dup = teff_to_arrays(num_d, den_d, bin_centers(edges_d))

    num_m, den_m, edges_m = combine_curve(args.mc, args.dir, var, "passed", args.rebin, file_cache,
                                           mode=args.mc_weight_mode, lumi_pb=args.lumi_pb)
    xm, ym, mlo, mup = teff_to_arrays(num_m, den_m, bin_centers(edges_m))

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.0], hspace=0.05)
    ax_eff = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax_eff)

    ax_eff.errorbar(xd, yd, yerr=[dlo, dup], fmt='o', ms=4, lw=2.3, capsize=2, label=args.data_label)
    ax_eff.errorbar(xm, ym, yerr=[mlo, mup], fmt='s', ms=4, lw=2.3, capsize=2, label=args.mc_label, color='orchid')

    x_common, idx_d, idx_m = np.intersect1d(xd, xm, return_indices=True)
    if len(x_common) == 0:
        raise RuntimeError("No common populated bins between data and MC for ratio panel.")
    r, rlo, rup = ratio_asymm(yd[idx_d], dlo[idx_d], dup[idx_d], ym[idx_m], mlo[idx_m], mup[idx_m])
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
        draw_target_eff_marker([ax_eff, ax_ratio], xd, yd, args.target_eff, args.plateau_n)

    return fig


# -----------------------------------------------------------------------
# Style: multicurve  (plotDijetHTEff.py / plotDijet_trigEffIncl.py aesthetic,
#                      now also covering -CMS.py's mplhep+petroff multi-sample
#                      overlay and -TRG.py's multi-hypothesis L1-seed overlay)
# -----------------------------------------------------------------------
def style_multicurve(args, var, file_cache):
    fig, ax = plt.subplots(figsize=(9, 8))

    if args.hypotheses:
        # -TRG.py / -junctionTRG.py pattern: vary the NUMERATOR SUFFIX on one
        # (or each) input file instead of varying the file. One curve per
        # (input, hypothesis) combination; if there's only one input, curves
        # are labeled by hypothesis alone (matches the original's usage).
        combos = []
        for path, in_label in zip(args.inputs, args.input_labels):
            for hyp in args.hypotheses:
                label = hyp if len(args.inputs) == 1 else f"{in_label}: {hyp}"
                combos.append((path, hyp, label))
    else:
        combos = [(path, "passed", label) for path, label in zip(args.inputs, args.input_labels)]

    curves = []
    for i, (path, num_suffix, label) in enumerate(combos):
        x, y, elo, eup = load_efficiency_curve(path, args.dir, var, num_suffix, args.rebin, file_cache)
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
        # marker computed off the first curve (conventionally the primary/data one)
        draw_target_eff_marker([ax], curves[0][0], curves[0][1], args.target_eff, args.plateau_n)

    return fig


# -----------------------------------------------------------------------
# Style: junction_ratio  (Mods Tested "-junction.py" aesthetic)
# -----------------------------------------------------------------------
def style_junction_ratio(args, var, file_cache):
    num_p, den_p, edges_p = combine_curve([args.primary], args.dir, var, "passed", args.rebin, file_cache, mode="none")
    xp, yp, plo, pup = teff_to_arrays(num_p, den_p, bin_centers(edges_p))

    # Weighted background stack -- mirrors -junction.py's Clone/Scale/Add
    # pattern. Two ways to specify weights, kept both for backward
    # compatibility: the original --bkg-weights positional list (if given),
    # or inline 'file.root:value' specs in --bkg itself (consistent with
    # --mc's syntax), interpreted per --bkg-weight-mode.
    if args.bkg_weights:
        if len(args.bkg_weights) != len(args.bkg):
            raise ValueError("--bkg-weights must have exactly one entry per --bkg file.")
        bkg_items = [f"{path}:{w}" for path, w in zip(args.bkg, args.bkg_weights)]
        bkg_mode = "raw"  # --bkg-weights has always meant a direct multiplier
    else:
        bkg_items = args.bkg
        bkg_mode = args.bkg_weight_mode

    num_b, den_b, edges_b = combine_curve(bkg_items, args.dir, var, "passed", args.rebin, file_cache,
                                           mode=bkg_mode, lumi_pb=args.lumi_pb)
    xb, yb, blo, bup = teff_to_arrays(num_b, den_b, bin_centers(edges_b))

    fig = plt.figure(figsize=(9, 9))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.05)
    ax_eff = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax_eff)

    ax_eff.errorbar(xp, yp, yerr=[plo, pup], fmt='o', ms=5, lw=2.0, capsize=2,
                     color=PALETTE[0], label=args.primary_label)
    ax_eff.errorbar(xb, yb, yerr=[blo, bup], fmt='s', ms=5, lw=2.0, capsize=2,
                     color=PALETTE[3], label=args.bkg_label)

    if args.ratio_style == "overlay":
        # Replicates what -junction.py's main() ACTUALLY renders in its lower
        # panel: both curves' raw efficiency values redrawn on a zoomed axis,
        # not a computed ratio (its create_ratio() function exists but is
        # never called -- see the PyROOT script's docstring for the verified
        # grep check). Offered as an option for exact visual parity with what
        # the original script outputs today, even though it isn't a ratio.
        ax_ratio.errorbar(xp, yp, yerr=[plo, pup], fmt='o', ms=4, lw=2.0, capsize=2, color=PALETTE[0])
        ax_ratio.errorbar(xb, yb, yerr=[blo, bup], fmt='s', ms=4, lw=2.0, capsize=2, color=PALETTE[3])
        ax_ratio.set_ylabel("Efficiency (zoom)")
    else:
        x_common, idx_p, idx_b = np.intersect1d(xp, xb, return_indices=True)
        if len(x_common) == 0:
            raise RuntimeError("No common populated bins between primary and background stack for ratio panel.")
        r, rlo, rup = ratio_asymm(yb[idx_b], blo[idx_b], bup[idx_b], yp[idx_p], plo[idx_p], pup[idx_p])
        # Neutral color -- NOT PALETTE[3] (the background's own color), which
        # previously made the ratio panel look like "only the background is
        # shown" since it visually matched the top panel's bkg marker exactly.
        ax_ratio.errorbar(x_common, r, yerr=[rlo, rup], fmt='o', ms=4, lw=2.0, capsize=2, color="black")
        ax_ratio.axhline(1.0, linestyle='--', linewidth=1.2, color='gray')
        ax_ratio.set_ylabel("Bkg/Primary")

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


STYLES = {
    "cms_ratio": style_cms_ratio,
    "multicurve": style_multicurve,
    "junction_ratio": style_junction_ratio,
}


def build_parser():
    p = argparse.ArgumentParser(description="Unified trigger-efficiency plotter (multi-style)")
    p.add_argument("--style", choices=STYLES.keys(), required=True)
    p.add_argument("--dir", default="InclusiveTrigNanoAOD")
    p.add_argument("--var", nargs="+", default=["ht_inclusive"],
                    help="One or more variable names, or the single value 'all' to "
                         "auto-discover and plot every h_{var}_all/_passed pair "
                         "present in the input file(s) for the chosen style.")
    p.add_argument("--rebin", type=int, default=1)
    p.add_argument("--year", default="2024")
    p.add_argument("--xlabel", default="")
    p.add_argument("--vline", type=float, default=None)
    p.add_argument("--target-eff", type=float, default=None,
                    help="Draw computed plateau-crossing marker at this efficiency (e.g. 0.90). Omit to disable.")
    p.add_argument("--plateau-n", type=int, default=5,
                    help="Number of trailing bins averaged to estimate the plateau.")
    p.add_argument("--ratio-ymin", type=float, default=0.5)
    p.add_argument("--ratio-style", choices=["ratio", "overlay"], default="ratio",
                    help="junction_ratio only. 'ratio' (default): the lower panel shows the computed "
                         "Bkg/Primary ratio, one series, neutral color. 'overlay': the lower panel "
                         "instead re-draws both Primary and Bkg raw efficiency curves on a zoomed axis "
                         "-- matches what -junction.py's main() actually renders today (its own "
                         "create_ratio() function is defined but never called there).")
    p.add_argument("--ratio-ymax", type=float, default=1.5)
    p.add_argument("--outdir", default="./plots")
    p.add_argument("--formats", nargs="+", default=[".png"])

    # cms_ratio
    p.add_argument("--data", nargs="+",
                    help="One or more data files, plain paths -- always summed unweighted.")
    p.add_argument("--mc", nargs="+",
                    help="One or more MC samples: 'file.root' (weight=1) or 'file.root:value' "
                         "(value's meaning depends on --mc-weight-mode).")
    p.add_argument("--mc-weight-mode", choices=["none", "raw", "xsec_sumw"], default="xsec_sumw",
                    help="How to interpret ':value' in --mc. Default 'xsec_sumw' (value=xsec_pb, "
                         "scale=xsec*lumi/sumw) matches compare_data_mc_eff.py. 'raw' treats value "
                         "as a direct Scale() multiplier (matches -junction.py's convention). "
                         "'none' ignores any ':value' and uses weight=1 for every sample.")
    p.add_argument("--lumi-pb", type=float, default=1000.0,
                    help="Integrated luminosity in pb^-1, used only in xsec_sumw weight mode.")
    p.add_argument("--data-label", default="Data")
    p.add_argument("--mc-label", default="MC")

    # multicurve
    p.add_argument("--inputs", nargs="+", default=[])
    p.add_argument("--input-labels", nargs="+", default=[])
    p.add_argument("--hypotheses", nargs="+", default=[],
                    help="Numerator histogram suffixes to overlay from the SAME input file(s) "
                         "(e.g. passed passedEXC passedL1HTT280 passedL1Jet passedL1HTDouble "
                         "passedL1HTSingle), reproducing -TRG.py's multi-hypothesis L1-seed "
                         "overlay. If omitted, one curve per --inputs file is drawn using 'passed'.")
    p.add_argument("--legend-title", default="", help="Optional legend title (e.g. 'Scenarios').")

    # junction_ratio
    p.add_argument("--primary")
    p.add_argument("--primary-label", default="Data")
    p.add_argument("--bkg", nargs="+", default=[],
                    help="Background samples: 'file.root' or 'file.root:value' (see --bkg-weight-mode). "
                         "If --bkg-weights is also given, it takes precedence (legacy positional form).")
    p.add_argument("--bkg-weights", nargs="+", type=float, default=[],
                    help="Legacy: one weight per --bkg file, positionally paired, always used as a "
                         "direct multiplier (equivalent to --bkg-weight-mode raw via inline ':value').")
    p.add_argument("--bkg-weight-mode", choices=["none", "raw", "xsec_sumw"], default="raw",
                    help="How to interpret ':value' in --bkg when --bkg-weights isn't used. Default "
                         "'raw' matches -junction.py's actual convention (xsec used directly as the "
                         "scale, no sumw/lumi division). Use 'xsec_sumw' for the more physically "
                         "correct xsec*lumi/sumw normalization instead.")
    p.add_argument("--bkg-label", default="Bkg (weighted sum)")

    return p


def main():
    args = build_parser().parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.style == "multicurve" and not args.input_labels:
        args.input_labels = [os.path.basename(p) for p in args.inputs]

    # Shared across every variable plotted in this invocation: each distinct
    # input file is opened with uproot exactly once here, no matter how many
    # variables get looped over below.
    file_cache = {}

    variables = resolve_variables(args, file_cache)
    print(f"[INFO] Plotting {len(variables)} variable(s): {variables}")

    n_ok, n_failed = 0, 0
    for var in variables:
        try:
            fig = STYLES[args.style](args, var, file_cache)
        except Exception as exc:
            # Loud, per-variable failure -- one missing/malformed histogram
            # shouldn't silently drop out of a multi-variable batch, nor
            # should it be swallowed the way the original scripts'
            # `except Exception: continue` blocks did.
            print(f"[ERROR] style='{args.style}' var='{var}' failed: {exc}")
            n_failed += 1
            continue

        tag = f"{args.style}_{var}"
        for fmt in args.formats:
            outname = os.path.join(args.outdir, f"TrigEff_{tag}{fmt}")
            fig.savefig(outname, dpi=200, bbox_inches="tight")
            print(f"Saved: {outname}")
        plt.close(fig)
        n_ok += 1

    print(f"[INFO] Done: {n_ok} plot(s) written, {n_failed} failed.")
    if n_failed and not n_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
