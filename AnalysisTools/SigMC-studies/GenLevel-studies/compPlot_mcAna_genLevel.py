import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

# -----------------------------
# --- Config ------------------
# -----------------------------
#mHyp = 50
#mHyp = 200
mHyp = 300

file_incl = uproot.open(f"HISTOS_MCSIG/histos_genLevel_GluonPairs_M{mHyp}.root")       
file_acc  = uproot.open(f"HISTOS_MCSIG/histos_genLevel_GluonPairs_M{mHyp}_wAccCuts.root") 
outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/GenLevel/wAccCuts/"

# --- Histograms dictionary ---
hist_1D = {
    "Graviton pT [GeV]": ("h_gravPt", file_incl["h_gravPt"], file_acc["h_gravPt"]),
    "Dijet mass [GeV]": ("h_dijetM", file_incl["h_dijetM"], file_acc["h_dijetM"]),
    r"$\Delta R(g,g)$": ("h_dR",     file_incl["h_dR"],     file_acc["h_dR"])
}

# -----------------------------
# --- CMS style 1D plots -----
# -----------------------------
plt.style.use(hep.style.CMS)

for xlabel, (hname, h_incl, h_acc) in hist_1D.items():
    # Extract values and bin edges
    values_incl = h_incl.values(flow=False)
    values_acc  = h_acc.values(flow=False)
    edges = h_incl.axes[0].edges()

    # Integrals (sum of bin contents)
    total_incl = np.sum(values_incl)
    total_acc  = np.sum(values_acc)
    acc_pct = (total_acc / total_incl * 100) if total_incl > 0 else 0.0

    fig, ax = plt.subplots(figsize=(12, 9))

    # Plot both histos
    hep.histplot(values_incl, edges, ax=ax, histtype="step", color="royalblue",
                 label="inclusive")
    hep.histplot(values_acc, edges, ax=ax, histtype="step", color="tomato",
                 label="w/ acceptance cuts")

    # CMS label
    hep.cms.label("Preliminary", data=False, ax=ax, lumi=None, year="2024", com="13.6")

    # Axis styling
    ymax = max(max(values_incl), max(values_acc)) * 1.3
    ax.set_ylim(0, ymax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events")

    # Legend with title
    ax.legend(title=f"GravitonToGG (M={mHyp} GeV)", loc="upper right")

    # Add acceptance text on canvas (top left corner)
    ax.text(0.65, 0.6, f"Sig acc = {acc_pct:.1f}%", transform=ax.transAxes,
            fontsize=24, fontweight="bold", color="black")

    # Save
    fig.savefig(f"{outdir}/genLevelAcceptance_{hname}_M{mHyp}.png", dpi=300)
    fig.savefig(f"{outdir}/genLevelAcceptance_{hname}_M{mHyp}.pdf")
    plt.close(fig)
    print(f"Saved 1D comparison plot: genLevel_{hname}_M{mHyp}, Sig acc = {acc_pct:.1f}%")
