"""
plot_mcAna_miniAOD.py

Plotter for PyROOT plus FWLite analyzer to:
 - compare dijet mass distribution for signal sample with different acceptane cuts

Usage:
  cmsenv
  python3 plot_mcAna_massComp.py --mHyp 200
"""

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import argparse
from ROOT import TFile

# ------------------------
# --- Argument parsing ---
# ------------------------
parser = argparse.ArgumentParser(description="MiniAOD gluon+jet analysis")
parser.add_argument("--mHyp", type=int, default=200,
                    help="Mass hypothesis (default: 200 GeV)")
args = parser.parse_args()

# --------------------------
# --- User configuration ---
# --------------------------
mHyp = args.mHyp
print(f"[INFO] Running with mHyp = {mHyp} GeV")

outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/"

# Input ROOT files and labels
files = {
    "pt > 30 GeV": "HISTOS_MCSIG/recoDijetMass_RSGraviton_M300_wPUPPIJets_FSR0p8_jetPt30Eta5.root",
    "pt > 72 GeV": "HISTOS_MCSIG/recoDijetMass_RSGraviton_M300_wPUPPIJets_FSR0p8_jetPt72Eta5.root",
    "pt > 100 GeV": "HISTOS_MCSIG/recoDijetMass_RSGraviton_M300_wPUPPIJets_FSR0p8_jetPt100Eta5.root",
}

colors = {
    "pt > 30 GeV": "mediumseagreen",
    "pt > 72 GeV": "darkorange",
    "pt > 100 GeV": "royalblue",
}

# Function to extract numpy arrays from ROOT histogram
def get_hist_arrays(hist):
    nbins = hist.GetNbinsX()
    edges = np.array([hist.GetBinLowEdge(i+1) for i in range(nbins+1)])
    values = np.array([hist.GetBinContent(i+1) for i in range(nbins)])
    return values, edges


# Apply CMS style
plt.style.use(hep.style.CMS)

# ----------------------------------- #
# --- Compare all acceptance cuts --- #
# ----------------------------------- #
fig, ax = plt.subplots(figsize=(10, 9))

for label, path in files.items():
    f = TFile.Open(path)
    h = f.Get("h_dijet_m")
    if not h:
        print(f"Histogram 'h_dijet_m' not found in {path}")
        continue
    values, edges = get_hist_arrays(h)
    hep.histplot(values, edges, ax=ax, histtype="step", label=label, color=colors[label])

ax.axvline(mHyp, linestyle="--", color="black", label=f"m = {mHyp} GeV (truth)")
ax.set_xlabel(r"$m_{jj}$ [GeV]")
ax.set_ylabel("Events")
ax.legend(title="Varying pt threshold:", fontsize=16)
hep.cms.label("Preliminary", data=False, ax=ax, year="2024", com="13.6")
ax.set_ylim(0, ax.get_ylim()[1]*1.3)

plt.tight_layout()
plt.savefig(f"{outdir}/compare_dijetMass_RSGraviton_M300_acceptances.png", dpi=300)
plt.savefig(f"{outdir}/compare_dijetMass_RSGraviton_M300_acceptances.pdf", dpi=300)


# --------------------------------- #
# --- Mass plot and resolution  --- #
# --------------------------------- #
file_pt100 = files["pt > 100 GeV"]
f = TFile.Open(file_pt100)
h = f.Get("h_dijet_m")
values, edges = get_hist_arrays(h)

fig2, ax2 = plt.subplots(figsize=(10, 9))
hep.histplot(values, edges, ax=ax2, histtype="step", color="dodgerblue", label="Jet pt > 100 GeV")

# --- Estimate resolution w/ RMS ---
#bin_centers = 0.5 * (edges[1:] + edges[:-1])
#peak_idx = np.argmax(values)
#peak = bin_centers[peak_idx]
#mean = h.GetMean()
#rms = h.GetRMS()
#resolution = rms / mean * 100  # in %

# --- Estimate peak and FWHM (Full Width at Half Maximum) ---
# -----------------------------------------------
bin_centers = 0.5 * (edges[1:] + edges[:-1])
peak_idx = np.argmax(values)
peak = bin_centers[peak_idx]
halfmax = values[peak_idx] / 2

left_idx = np.where(values[:peak_idx] < halfmax)[0]
right_idx = np.where(values[peak_idx:] < halfmax)[0]

if len(left_idx) > 0:
    left_edge = bin_centers[left_idx[-1]]
else:
    left_edge = bin_centers[0]

if len(right_idx) > 0:
    right_edge = bin_centers[peak_idx + right_idx[0]]
else:
    right_edge = bin_centers[-1]

fwhm = right_edge - left_edge
resolution = fwhm / (2.35 * peak) * 100  # approximate sigma from FWHM assuming Gaussian

# Draw double arrow at half max width around peak
# -----------------------------------------------
ax2.annotate("", xy=(left_edge, halfmax), xytext=(right_edge, halfmax),
             arrowprops=dict(arrowstyle='<->', color='black', lw=2))
ax2.text(peak*1.6, halfmax * 1.02,
         f"$\sigma$ = {fwhm/2.35:.1f} GeV ({resolution:.1f}%)",
         color="black", fontsize=18, ha="center", va="bottom", fontstyle="italic")

# CMS label
# -----------------------------------------------
hep.cms.label("Preliminary", data=False, ax=ax2, year="2024", com="13.6")
#ax.axvline(mHyp, linestyle="--", color="red", label=f"m = {mHyp} GeV (truth)")
ax2.axvline(peak, color="darkturquoise", linestyle="--", label=f"Mean = {peak:.1f} GeV")
ax2.set_xlabel(r"$m_{jj}$ [GeV]")
ax2.set_ylabel("Events")
ax2.legend(loc="upper right")
ax2.set_ylim(0, max(values) * 1.3)
ax2.set_xlim(0, 600)

plt.tight_layout()
plt.savefig(f"{outdir}/dijetMass_RSGraviton_M300_jetPt100_resolution.png", dpi=300)
plt.savefig(f"{outdir}/dijetMass_RSGraviton_M300_jetPt100_resolution.pdf", dpi=300)


# ----------------------------------------------- #
# --- Mass plot and resolution: jpt > 72 GeV  --- #
# ----------------------------------------------- #
file_pt72 = files["pt > 72 GeV"]
f = TFile.Open(file_pt72)
h = f.Get("h_dijet_m")
values, edges = get_hist_arrays(h)

fig2, ax2 = plt.subplots(figsize=(10, 9))
hep.histplot(values, edges, ax=ax2, histtype="step", color="dodgerblue", label="Jet pt > 72 GeV")

# Estimate peak and FWHM (Full Width at Half Maximum) 
# ---------------------------------------------------
bin_centers = 0.5 * (edges[1:] + edges[:-1])
peak_idx = np.argmax(values)
peak = bin_centers[peak_idx]
halfmax = values[peak_idx] / 2

left_idx = np.where(values[:peak_idx] < halfmax)[0]
right_idx = np.where(values[peak_idx:] < halfmax)[0]

if len(left_idx) > 0:
    left_edge = bin_centers[left_idx[-1]]
else:
    left_edge = bin_centers[0]

if len(right_idx) > 0:
    right_edge = bin_centers[peak_idx + right_idx[0]]
else:
    right_edge = bin_centers[-1]

fwhm = right_edge - left_edge
resolution = fwhm / (2.35 * peak) * 100  # approximate sigma from FWHM assuming Gaussian

# Draw double arrow at half max width around peak
# -----------------------------------------------
ax2.annotate("", xy=(left_edge, halfmax), xytext=(right_edge, halfmax),
             arrowprops=dict(arrowstyle='<->', color='black', lw=2))
ax2.text(peak*1.7, halfmax * 1.02,
         f"$\sigma$ = {fwhm/2.35:.1f} GeV ({resolution:.1f}%)",
         color="black", fontsize=18, ha="center", va="bottom", fontstyle="italic")

# CMS label
# -----------------------------------------------
hep.cms.label("Preliminary", data=False, ax=ax2, year="2024", com="13.6")
#ax.axvline(mHyp, linestyle="--", color="red", label=f"m = {mHyp} GeV (truth)")
ax2.axvline(peak, color="darkturquoise", linestyle="--", label=f"Mean = {peak:.1f} GeV")
ax2.set_xlabel(r"$m_{jj}$ [GeV]")
ax2.set_ylabel("Events")
ax2.legend(loc="upper right")
ax2.set_ylim(0, max(values) * 1.3)
ax2.set_xlim(0, 600)

plt.tight_layout()
plt.savefig(f"{outdir}/dijetMass_RSGraviton_M300_jetPt72_resolution.png", dpi=300)
plt.savefig(f"{outdir}/dijetMass_RSGraviton_M300_jetPt72_resolution.pdf", dpi=300)

print("---------------------------------------------------------------------------------------")
print(f"All plots saved! Location: {outdir}")
