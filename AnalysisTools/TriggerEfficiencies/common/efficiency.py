#!/usr/bin/env python3
"""
common/efficiency.py

ROOT-independent efficiency core for the unified trigger-efficiency plotter.

Design rationale
-----------------
The existing scripts (compare_data_mc_eff.py, plotDijetHTEff.py, the
"Mods Tested/Plot Mods" variants) all lean on PyROOT's ROOT.TEfficiency for
the Clopper-Pearson interval. PyROOT is a heavy, environment-specific
dependency (CMSSW/LCG stack) that isn't available in every place someone
might want to smoke-test this code (e.g. a plain CI runner or this sandbox).

This module re-implements exactly what's needed -- Clopper-Pearson efficiency
+ asymmetric errors -- using scipy.stats.beta, which is the same underlying
formula ROOT uses for kFCP. Histograms are read with uproot (pure Python,
no ROOT install required) instead of ROOT.TFile.

If a real CMSSW environment *does* have PyROOT available, nothing here
needs to change: this module never imports ROOT, and the plotting script
only depends on this module's functions, not on TEfficiency directly. That
means the exact same plotting code path is used whether you're running
against real CERN-hosted ROOT files with PyROOT elsewhere, or smoke-testing
here with synthetic files read via uproot.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

# ROOT's default confidence level for TEfficiency (1 sigma, ~68.27%)
DEFAULT_CL = 0.682689492137

# -----------------------------------------------------------------------
# Clopper-Pearson interval (matches ROOT::TEfficiency::ClopperPearson)
# -----------------------------------------------------------------------
def clopper_pearson(k, n, cl=DEFAULT_CL):
    """
    Vectorized Clopper-Pearson confidence interval for a binomial efficiency.

    k, n : array-like (passed counts, total counts)
    cl   : central confidence level (ROOT default: 1-sigma equivalent)

    Returns (eff, err_low, err_up) as numpy arrays, matching
    TEfficiency::GetEfficiency / GetEfficiencyErrorLow / GetEfficiencyErrorUp
    for StatisticOption kFCP.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    alpha = 1.0 - cl
    eff = np.divide(k, n, out=np.zeros_like(k), where=n > 0)

    lo = np.where(
        k == 0, 0.0,
        stats.beta.ppf(alpha / 2.0, np.clip(k, 1e-12, None), n - k + 1)
    )
    hi = np.where(
        k == n, 1.0,
        stats.beta.ppf(1.0 - alpha / 2.0, k + 1, np.clip(n - k, 1e-12, None))
    )
    lo = np.nan_to_num(lo, nan=0.0)
    hi = np.nan_to_num(hi, nan=1.0)

    err_lo = eff - lo
    err_up = hi - eff
    return eff, err_lo, err_up


# -----------------------------------------------------------------------
# Histogram I/O (uproot instead of ROOT.TFile)
# -----------------------------------------------------------------------
def read_th1(uproot_file, dirname, hist_name):
    """
    Read a 1D histogram from a (possibly nested) TDirectory inside a ROOT
    file opened with uproot. Returns (edges, values) as numpy arrays.
    Raises KeyError with a clear message if missing -- deliberately loud,
    same philosophy as the original scripts' `raise RuntimeError` on a
    missing histogram, since a silently-skipped histogram is exactly the
    kind of structural bug this harness exists to catch.
    """
    key = f"{dirname}/{hist_name}" if dirname else hist_name
    try:
        h = uproot_file[key]
    except KeyError as exc:
        raise KeyError(
            f"Missing histogram '{hist_name}' in directory '{dirname}' "
            f"of {uproot_file.file_path}"
        ) from exc
    values, edges = h.to_numpy()
    return np.asarray(edges, dtype=float), np.asarray(values, dtype=float)


def get_sumw(uproot_file, dirname, hist_name="h_sumw"):
    """uproot-backend sumw reader (PyROOT backend reuses the imported
    original's own get_sumw() instead of this one)."""
    key = f"{dirname}/{hist_name}" if dirname else hist_name
    try:
        h = uproot_file[key]
    except KeyError as exc:
        raise KeyError(f"Missing {hist_name} in directory '{dirname}' of {uproot_file.file_path}") from exc
    values, _ = h.to_numpy()
    sumw = float(np.sum(values))
    if sumw <= 0:
        raise RuntimeError(f"Invalid sumw={sumw} from {hist_name} in {dirname}")
    return sumw


def rebin_counts(edges, values, factor=1):
    """Integer rebin by summing adjacent bins (mirrors TH1::Rebin)."""
    factor = int(factor)
    if factor <= 1:
        return edges, values
    nbins = len(values)
    if nbins % factor != 0:
        usable = (nbins // factor) * factor
        values = values[:usable]
        edges = edges[: usable + 1]
    new_values = values.reshape(-1, factor).sum(axis=1)
    new_edges = edges[::factor]
    return new_edges, new_values


def bin_centers(edges):
    edges = np.asarray(edges, dtype=float)
    return 0.5 * (edges[:-1] + edges[1:])


# -----------------------------------------------------------------------
# x-axis label map (mirrors compare_data_mc_eff.py / comp_plotDijetHTEff.py's
# xlabel_map dicts, unified into one place instead of duplicated per script)
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# CMS-approved qualitative color palettes (Petroff scheme), copied exactly
# from the values literally defined in Last Dev/plotDijetHTEff_AN_Ortogonal
# -Mjj-CMS.py -- NOT approximated. PETROFF6 also matches mplhep's own
# hep.style.use("CMS") default axes.prop_cycle exactly (verified against
# mplhep 1.3.3), confirming these are the real values, not a guess.
# -----------------------------------------------------------------------
PETROFF6 = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
PETROFF10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
             "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]

XLABEL_MAP = {
    "ht_inclusive": r"$H_{T}$ Inclusive [GeV]",
    "pt_leading":   r"AK4 leading jet $p_{T}$ [GeV]",
    "dijetMass":    r"Dijet mass [GeV]",
    "minv_1":       r"$m_{jj}$ VBF [GeV]",
    "ht_1":         r"$H_{T}$ VBF [GeV]",
    "minv_2":       r"$m_{jj}$ Boosted ISR [GeV]",
    "ht_2":         r"$H_{T}$ Boosted ISR [GeV]",
    "ptak8_2":      r"AK8 $p_{T}$ Boosted ISR [GeV]",
    "minv_3":       r"$m_{jj}$ Resolved ISR [GeV]",
    "ht_3":         r"$H_{T}$ Resolved ISR [GeV]",
    "minv_4":       r"$m_{jj}$ Rest [GeV]",
    "ht_4":         r"$H_{T}$ Rest [GeV]",
}


def get_xlabel(var, override=""):
    """Explicit --xlabel always wins; otherwise use the shared map; otherwise fall back to the raw name."""
    if override:
        return override
    return XLABEL_MAP.get(var, var)


# -----------------------------------------------------------------------
# Shared weight-spec parsing (used by cms_ratio's --mc and junction_ratio's
# --bkg in both backend scripts). Three modes found across the original
# scripts, verified by reading their actual source rather than assumed:
#   "none"      -- plotDijetHTEff.py / -trigEffIncl.py / -CMS.py / -TRG.py:
#                  each curve is a pre-made histogram, no scaling at all.
#   "raw"       -- -junction.py: `weights = [1.0, 25360000.0, ...]` (a
#                  comment literally says "xsec took from samples codes"),
#                  applied directly as histog.Scale(weights[k]) -- NO
#                  division by sumw or multiplication by lumi. This is an
#                  approximation in the original, kept faithfully as an
#                  available mode rather than silently "corrected".
#   "xsec_sumw" -- compare_data_mc_eff.py / plotDijet_mcTrigEffIncl.py:
#                  scale = xsec_pb * lumi_pb / sumw (get_sumw() reads
#                  h_sumw.Integral()). The physically correct normalization
#                  for combining multiple MC samples with different
#                  generator cross sections.
# -----------------------------------------------------------------------
def parse_weight_spec(item, mode="none"):
    """
    item: 'file.root' or 'file.root:value'. The colon suffix is always
    OPTIONAL -- a bare path never errors, regardless of mode (fixes the
    friction of --mc requiring ':xsec_pb' even for a quick single-file
    test). Returns (path, value_or_None).

    mode="none":      value is ignored if present (printed as a no-op note).
    mode="raw":       value is used directly as a Scale() multiplier;
                      defaults to 1.0 if omitted.
    mode="xsec_sumw": value is REQUIRED (a cross-section in pb) -- there's
                      no safe default for an unknown xsec; raises if a bare
                      path is given in this mode.
    """
    if ":" in item:
        path, val_str = item.rsplit(":", 1)
        value = float(val_str)
    else:
        path, value = item, None

    if mode == "xsec_sumw" and value is None:
        raise ValueError(
            f"'{item}': --*-weight-mode xsec_sumw requires 'file.root:xsec_pb' "
            f"for every sample (no safe default for an unknown cross section)."
        )
    if mode == "none" and value is not None:
        print(f"[INFO] '{item}': weight-mode is 'none', ignoring the ':{val_str}' suffix.")
    return path, value


def resolve_scale(value, mode, sumw=None, lumi_pb=1000.0):
    """Turn a parsed weight-spec value into the actual histogram Scale() factor."""
    if mode == "none":
        return 1.0
    if mode == "raw":
        return 1.0 if value is None else value
    if mode == "xsec_sumw":
        if sumw is None or sumw <= 0:
            raise RuntimeError(f"Invalid or missing sumw ({sumw}) for xsec_sumw weighting.")
        return value * lumi_pb / sumw
    raise ValueError(f"Unknown weight mode: {mode}")


def match_variable_pairs(key_names):
    """
    Backend-agnostic core of variable discovery: given any flat set/list of
    histogram name strings (from uproot's file.keys(), or from
    {k.GetName() for k in TDirectory.GetListOfKeys()} on the PyROOT side),
    return every var with both h_{var}_all and h_{var}_passed present.
    Both unified_trigeff_plotter.py and unified_trigeff_plotter_pyroot.py
    call this after extracting their own backend's key names, so the
    matching rule itself has exactly one definition.
    """
    keys = set(key_names)
    found = []
    for k in keys:
        if not k.startswith("h_") or not k.endswith("_all"):
            continue
        var = k[2:-4]
        if f"h_{var}_passed" in keys:
            found.append(var)
    return sorted(found)


def discover_variables(uproot_file, dirname):
    """
    uproot-backend convenience wrapper around match_variable_pairs(): list
    every variable that has both h_{var}_all and h_{var}_passed in the given
    directory, instead of relying on a hardcoded variable list that may not
    match what a particular file actually contains (inclusive files have
    ht_inclusive/pt_leading/dijetMass; categorized files have
    minv_1..4/ht_1..4/ptak8_2; a file may have only some of these).
    """
    prefix = f"{dirname}/" if dirname else ""
    raw_keys = {k.split(";")[0] for k in uproot_file.keys()}
    # strip the directory prefix so match_variable_pairs sees bare h_*_all names
    stripped = {k[len(prefix):] for k in raw_keys if k.startswith(prefix + "h_")}
    return match_variable_pairs(stripped)


# -----------------------------------------------------------------------
# num/den -> efficiency arrays (replaces compare_data_mc_eff.teff_to_arrays)
# -----------------------------------------------------------------------
def teff_to_arrays(num_values, den_values, centers, cl=DEFAULT_CL):
    """
    Mirrors the original teff_to_arrays(): only keeps bins where the
    denominator is populated (den > 0), same as the ROOT-based version's
    `if total <= 0: continue` guard.
    """
    num_values = np.asarray(num_values, dtype=float)
    den_values = np.asarray(den_values, dtype=float)
    centers = np.asarray(centers, dtype=float)

    mask = den_values > 0
    if not np.any(mask):
        raise RuntimeError("No populated denominator bins -- cannot build efficiency curve.")

    eff, elo, eup = clopper_pearson(num_values[mask], den_values[mask], cl=cl)
    return centers[mask], eff, elo, eup


def check_consistency(num_values, den_values, label=""):
    """Mirrors ROOT.TEfficiency.CheckConsistency: numerator must never exceed denominator."""
    num_values = np.asarray(num_values, dtype=float)
    den_values = np.asarray(den_values, dtype=float)
    bad = num_values > den_values + 1e-9
    if np.any(bad):
        bad_bins = np.where(bad)[0]
        raise RuntimeError(
            f"TEfficiency consistency check failed{f' for {label}' if label else ''}: "
            f"numerator exceeds denominator in bin(s) {bad_bins.tolist()}"
        )


def ratio_asymm(data_eff, data_lo, data_up, mc_eff, mc_lo, mc_up, eps=1e-12):
    """Unchanged from compare_data_mc_eff.py -- same asymmetric-ratio propagation."""
    d = np.clip(data_eff, eps, None)
    m = np.clip(mc_eff, eps, None)
    r = d / m

    rel_d_lo = data_lo / np.clip(d, eps, None)
    rel_d_up = data_up / np.clip(d, eps, None)
    rel_m_lo = mc_lo / np.clip(m, eps, None)
    rel_m_up = mc_up / np.clip(m, eps, None)

    r_lo = r * np.sqrt(rel_d_lo ** 2 + rel_m_up ** 2)
    r_up = r * np.sqrt(rel_d_up ** 2 + rel_m_lo ** 2)
    return r, r_lo, r_up


# -----------------------------------------------------------------------
# NEW (Step 4a): computed 90%-of-plateau crossing, replacing the fixed --vline
# -----------------------------------------------------------------------
def find_target_eff_crossing(centers, eff, target_eff=0.90, plateau_n=5):
    """
    Estimate the plateau value as the mean efficiency over the last
    `plateau_n` populated bins, then find the x-position where the curve
    first crosses (target_eff * plateau) via linear interpolation between
    the two bracketing bins.

    Returns (x_cross, plateau_value) or (None, plateau_value) if the curve
    never reaches the target within the available range.
    """
    centers = np.asarray(centers, dtype=float)
    eff = np.asarray(eff, dtype=float)
    if len(centers) == 0:
        raise ValueError("Empty efficiency curve -- cannot estimate plateau.")

    order = np.argsort(centers)
    centers, eff = centers[order], eff[order]

    n_tail = min(plateau_n, len(eff))
    plateau = float(np.mean(eff[-n_tail:]))
    threshold = target_eff * plateau

    x_cross = None
    for i in range(1, len(eff)):
        y0, y1 = eff[i - 1], eff[i]
        if (y0 - threshold) * (y1 - threshold) <= 0 and y1 != y0:
            x0, x1 = centers[i - 1], centers[i]
            frac = (threshold - y0) / (y1 - y0)
            x_cross = float(x0 + frac * (x1 - x0))
            break

    return x_cross, plateau
