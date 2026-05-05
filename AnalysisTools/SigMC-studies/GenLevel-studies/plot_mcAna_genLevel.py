import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

# -----------------------------
# --- Load ROOT histograms ----
# -----------------------------
mHyp = 100
file = uproot.open(f"HISTOS_MCSIG/histos_genLevel_GluonPairs_M{mHyp}.root")
outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/GenLevel/"
#file = uproot.open(f"HISTOS_MCSIG/histos_genLevel_GluonPairs_M{mHyp}_wAccCuts.root")
#outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/GenLevel/wAccCuts"

hist_1D = {
    "Graviton pT [GeV]": file["h_gravPt"],
    "Dijet mass [GeV]": file["h_dijetM"],
    r"$\Delta R(g,g)$": file["h_dR"]
}

h2 = file["h_dRvsPt"]  # 2D histogram: DeltaR vs Graviton Pt

# -----------------------------
# --- CMS style 1D plots -----
# -----------------------------
plt.style.use(hep.style.CMS)

for xlabel, h in hist_1D.items():
    values = h.values(flow=False)
    edges = h.axes[0].edges()

    fig, ax = plt.subplots(figsize=(12,9))
    hep.histplot(values, edges, ax=ax, histtype="step", color="royalblue",
                 label=f"GravitonToGG (M={mHyp} GeV)")

    hep.cms.label("Preliminary", data=False, ax=ax, lumi=None, year="2024", com="13.6")

    ax.set_ylim(0, max(values)*1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events")
    ax.legend(loc = "upper right")
    
    # Save
    fig.savefig(f"{outdir}/genLevel_{h.name}_M{mHyp}.png", dpi=300)
    fig.savefig(f"{outdir}/genLevel_{h.name}_M{mHyp}.pdf")
    plt.close(fig)
    print(f"Saved 1D plot: genLevel_{h.name}")


# -----------------------------
# --- CMS style 2D plot -------
# -----------------------------
# Extract bin edges and contents
x_edges = h2.axis(0).edges()  # Graviton pT
y_edges = h2.axis(1).edges()  # DeltaR
contents = h2.values()   # 2D array

X, Y = np.meshgrid(x_edges, y_edges)

fig, ax = plt.subplots(figsize=(12,8))
pcm = ax.pcolormesh(X, Y, contents.T, cmap="viridis")  
cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label("Events")

ax.set_xlabel(r"Graviton $p_T$ [GeV]")
ax.set_ylabel(r"$\Delta R(g,g)$")
#ax.set_title("Gluon pair DR vs Graviton pT")
hep.cms.label("Preliminary", data=False, ax=ax, lumi=None, year="2024", com="13.6")

plt.tight_layout()
fig.savefig(f"{outdir}/genLevel_dR_vs_GravPt_M{mHyp}.png", dpi=300)
fig.savefig(f"{outdir}/genLevel_dR_vs_GravPt_M{mHyp}.pdf")
plt.close(fig)
print("Saved 2D DR vs Graviton Pt plot")
