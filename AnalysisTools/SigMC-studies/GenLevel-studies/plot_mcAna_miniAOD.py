"""
plot_mcAna_miniAOD.py

Plotter for PyROOT plus FWLite analyzer to:
 - read MINIAODSIM files
 - find generator gluons (pdgId=21) coming from mother with pdgId=5100039
 - plot generator level quantities

Usage:
  cmsenv
  python3 plot_mcAna_miniAOD.py --mHyp 200
"""

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import argparse

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

infile = f"HISTOS_MCSIG/recoDijetMass_RSGraviton_M{mHyp}_wPUPPIJets_FSR0p8_jetPt72Eta5.root"
outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/"

# Apply CMS style
plt.style.use(hep.style.CMS)

# -----------------------------
# --- Helper for plotting -----
# -----------------------------
def plot_hist(values, edges, xlabel, ylabel, outname, color="dodgerblue",
              label=None, vline=None, yrange=None, logy=False):
    fig, ax = plt.subplots(figsize=(10, 9))
    hep.histplot(values, edges, ax=ax, histtype="step", color=color, label=label)

    if vline:
        ax.axvline(vline, linestyle="--", color="red")

    hep.cms.label("Preliminary", data=False, ax=ax, year="2024", com="13.6")
    ax.set_ylim(0, max(values)*1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if yrange:
        ax.set_ylim(*yrange)
    if label:
        ax.legend()
    if logy:
        ax.set_yscale("log")

    fig.savefig(f"{outdir}/{outname}_M{mHyp}_wPUPPIJets_FSR0p8_jetPt100Eta5.png", dpi=300)
    fig.savefig(f"{outdir}/{outname}_M{mHyp}_wPUPPIJets_FSR0p8_jetPt100Eta5.pdf", dpi=300)
    plt.close(fig)


# -----------------------------
# --- Load ROOT histograms ----
# -----------------------------
with uproot.open(infile) as f:
    h_mjj            = f["h_dijet_m"]
    h_mjj_preFSR     = f["h_dijet_m_preFSR"]
    h_njets          = f["h_njetsTotal"]
    h_njets_ISRCands = f["h_njetsAdditional"]
    h_dRgg           = f["h_dRgg"]

    h_ISR_pt  = [f[f"h_ISR{i}_pt"]  for i in range(1, 4)]
    h_ISR_eta = [f[f"h_ISR{i}_eta"] for i in range(1, 4)]
    h_ISR_phi = [f[f"h_ISR{i}_phi"] for i in range(1, 4)]

# Extract values and bin edges
def get_vals_edges(h):
    return h.values(flow=False), h.axes[0].edges()

values_mjj, edges_mjj = get_vals_edges(h_mjj)
values_mjj_preFSR, edges_mjj_preFSR = get_vals_edges(h_mjj_preFSR)
values_njets, edges_njets = get_vals_edges(h_njets)
values_dRgg, edges_dRgg = get_vals_edges(h_dRgg)
values_njets_ISRCands, edges_njets_ISRCands = get_vals_edges(h_njets_ISRCands)

# -----------------------------
# --- Core plots --------------
# -----------------------------
plot_hist(h_njets.values(flow=False), h_njets.axes[0].edges(),
          xlabel="Number of slimmedJets per event", ylabel="Events",
          outname="njets", color="dodgerblue")

plot_hist(h_dRgg.values(flow=False), h_dRgg.axes[0].edges(),
          xlabel=r"$\Delta R(g,g)$ ", ylabel="Events",
          outname="dR_gg", color="dodgerblue")

#plot_hist(h_mjj.values(flow=False), h_mjj.axes[0].edges(),
#          xlabel=r"$m_{jj}$ [GeV]", ylabel="Events",
#          outname="dijet_mass_wFSR", label="Reco dijet mass (matched)", vline=mHyp)

#plot_hist(h_mjj_preFSR.values(flow=False), h_mjj_preFSR.axes[0].edges(),
#          xlabel=r"$m_{jj}$ [GeV]", ylabel="Events",
#          outname="dijet_mass_preFSR", label="Reco dijet mass (pre-FSR)", vline=mHyp)

#plot_hist(h_njets_ISRCands.values(flow=False), h_njets_ISRCands.axes[0].edges(),
#          xlabel="Number of ISR candidates", ylabel="Events",
#          outname="njets_ISRCands_wFSR", color="darkgreen")


# ----------------------------------------------
# --- Dijet mass before/after FSR comparison ---
# ----------------------------------------------
fig, ax = plt.subplots(figsize=(10,9))
hep.histplot(values_mjj_preFSR, edges_mjj_preFSR, ax=ax, histtype="step", color="darkorange", label="Pre-FSR")
hep.histplot(values_mjj, edges_mjj, ax=ax, histtype="step", color="royalblue", label="Post-FSR")
ax.axvline(mHyp, linestyle="--", color="red", label=f"m = {mHyp} GeV (truth)")
hep.cms.label("Preliminary", data=False, ax=ax, year="2024", com="13.6")
ax.set_ylim(0, max(values_mjj_preFSR)*1.2)
ax.set_xlabel(r"$m_{gg}$ [GeV]")
ax.set_ylabel("Events")
ax.legend()
fig.savefig(f"{outdir}/dijet_mass_prePostFSR_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.png", dpi=300)
fig.savefig(f"{outdir}/dijet_mass_prePostFSR_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.pdf", dpi=300)

# -----------------------------
# --- ISR comparisons ---------
# -----------------------------
colors = ["dodgerblue", "orange", "hotpink"]

# ISR pT comparison
fig, ax = plt.subplots(figsize=(10, 9))
for i, h in enumerate(h_ISR_pt):
    values, edges = h.values(flow=False), h.axes[0].edges()
    hep.histplot(values, edges, ax=ax, histtype="step", color=colors[i], label=f"Jet {i+1}")
hep.cms.label("Preliminary", data=False, ax=ax, year="2024", com="13.6")
#ax.set_ylim(0, max(values)*1.5)
ax.set_xlabel(r"ISR jet $p_T$ [GeV]")
ax.set_ylabel("Events")
ax.set_yscale("log")  
ax.legend(title="ISR candidates")
fig.savefig(f"{outdir}/ISR_pt_comparison_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.png", dpi=300)
fig.savefig(f"{outdir}/ISR_pt_comparison_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.pdf", dpi=300)
plt.close(fig)

# ISR eta comparison
fig, ax = plt.subplots(figsize=(10, 9))
for i, h in enumerate(h_ISR_eta):
    values, edges = h.values(flow=False), h.axes[0].edges()
    hep.histplot(values, edges, ax=ax, histtype="step", color=colors[i], label=f"Jet {i+1}")
hep.cms.label("Preliminary", data=False, ax=ax, year="2024", com="13.6")
#ax.set_ylim(0, max(values)*1.5)
ax.set_xlabel(r"ISR jet $\eta$")
ax.set_ylabel("Events")
ax.set_yscale("log")  
ax.legend(title="ISR candidates")
fig.savefig(f"{outdir}/ISR_eta_comparison_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.png", dpi=300)
fig.savefig(f"{outdir}/ISR_eta_comparison_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.pdf", dpi=300)
plt.close(fig)

# --- Number of ISR candidates plot (with percentages) ---
fig, ax = plt.subplots(figsize=(10,9))
hep.histplot(values_njets_ISRCands, edges_njets_ISRCands, ax=ax, histtype="step", color="dodgerblue")
hep.cms.label("Preliminary", data=False, ax=ax, year="2024", com="13.6")
ax.set_ylim(0, max(values_njets_ISRCands)*1.2)
ax.set_xlabel("Number of slimmedJets (ISR candidates)")
ax.set_ylabel("Events")
# Compute percentages
total_events = np.sum(values_njets_ISRCands)
percentages = values_njets_ISRCands / total_events * 100
# Annotate on canvas
for i, (val, perc) in enumerate(zip(values_njets_ISRCands, percentages)):
    if val > 20:
        ax.text(edges_njets_ISRCands[i] + 0.22, val + max(values_njets_ISRCands)*0.025,
                f"{perc:.1f}%", rotation=0, fontsize=14, color="black")

fig.savefig(f"{outdir}/njets_ISRCands_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.png", dpi=300)
fig.savefig(f"{outdir}/njets_ISRCands_M{mHyp}_wPUPPIJets_newFSR0p8_jetPt100Eta5.pdf", dpi=300)

print("---------------------------------------------------------------------------------------")
print(f"All plots saved! Location: {outdir}")
