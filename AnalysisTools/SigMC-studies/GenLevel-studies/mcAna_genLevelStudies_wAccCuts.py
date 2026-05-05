#!/usr/bin/env python3

"""
mcAna_genLevelStudies.py

PyROOT plus FWLite analyzer to:
 - read MINIAODSIM files
 - find generator gluons (pdgId=21) coming from mother with pdgId=5100039
 - plot generator level quantities

Usage:
  cmsenv
  python3 mcAna_genLevelStudies.py
"""


import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import glob, os, math


# --------------
# --- Config ---
# --------------
#mHyp = 50
mHyp = 100
#mHyp = 200
#mHyp = 300

file_pattern = f"/eos/user/e/elfontan/SampleFactory/RSGravitonToGluonGluon_M{mHyp}_TuneCP5_13p6TeV_pythia8/chain_RunIII2024Summer24GS-RunIII2024Summer24MiniAOD/*/RunIII2024Summer24MiniAOD_*root"
maxEvents = -1
print_every = 500
outdir = "/eos/user/e/elfontan/www/dijetAnaRun3/MCSignal/GenLevel/"
os.makedirs(outdir, exist_ok=True)

# ------------------
# --- Histograms ---
# ------------------
ROOT.gROOT.SetBatch(True)
h_gravPt  = ROOT.TH1F("h_gravPt",  "Graviton pT; pT [GeV]; Events", 50, 0, 500)
h_dijetM  = ROOT.TH1F("h_dijetM",  "Dijet mass from gluons; m_{jj} [GeV]; Events", 50, 0, 500)
h_dR      = ROOT.TH1F("h_dR",      "#DeltaR between gluons; #DeltaR; Events", 50, 0, 6)
h_dRvsPt  = ROOT.TH2F("h_dRvsPt",  "#DeltaR vs Graviton pT; pT [GeV]; #DeltaR", 50, 0, 500, 50, 0, 6)


# ---------------
# --- Helpers ---
# ---------------
def deltaR(eta1, phi1, eta2, phi2):
    dphi = abs(phi1 - phi2)
    while dphi > math.pi:
        dphi = (2*math.pi - dphi)
    deta = eta1 - eta2
    return math.sqrt(deta*deta + dphi*dphi)

def find_mother_with_pdg(genparticle, target_pdg):
    """
    Traverse mothers (recursively through immediate mothers) and return True if any ancestor has |pdgId| == target_pdg.
    prunedGenParticles -> has mother(index) accessor.
    """
    stack = [genparticle]
    visited = set()
    while stack:
        gp = stack.pop()
        gid = id(gp)
        if gid in visited:
            continue
        visited.add(gid)
        nmat = gp.numberOfMothers()
        for im in range(nmat):
            try:
                mom = gp.mother(im)
            except Exception:
                continue
            if mom is None:
                continue
            try:
                mpid = abs(mom.pdgId())
            except Exception:
                continue
            if mpid == target_pdg:
                return True
            stack.append(mom)
    return False


# -----------------------------
# --- Prepare FWLite events ---
# -----------------------------
from DataFormats.FWLite import Events, Handle

files = sorted(glob.glob(file_pattern))
if not files:
    raise RuntimeError(f"No files found for pattern: {file_pattern}")
print(f"-----------------------------")
print(f"------ Found {len(files)} files")
print(f"-----------------------------")

events = Events(files)

# Handles
# ---------------
gen, genLabel = Handle("GenEventInfoProduct"), ("generator")
handlePruned, prunedLabel  = Handle("std::vector<reco::GenParticle>"), ("prunedGenParticles")


# ------------
# --- Loop ---
# ------------
evt_count = 0
selected_pairs = 0

for event in events:
    evt_count += 1
    if 0 < maxEvents <= evt_count:
        break
    if evt_count % print_every == 0:
        print(f"Processed {evt_count} events")

    # Get pruned gen particles
    # ------------------------
    event.getByLabel(prunedLabel, handlePruned) #event.getByLabel(genLabel, gen)
    pruned = handlePruned.product()

    # Collect gluons and graviton
    # ---------------------------    
    gluons = []
    graviton_p4 = None

    for gp in pruned:
        if gp.pdgId() == 5100039:  # Graviton
            graviton_p4 = ROOT.TLorentzVector()
            graviton_p4.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy())
        if (abs(gp.pdgId()) == 21):

            if find_mother_with_pdg(gp, 5100039):
                p4 = ROOT.TLorentzVector()
                p4.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy())
                if (p4.Pt() > 30 and abs(p4.Eta()) < 5.0):
                    gluons.append(p4)
                
    # Find pair closest to Graviton mass
    # ----------------------------------
    best_pair = None
    best_mass_diff = 1e6
    if (len(gluons) < 2): continue
    
    for i in range(len(gluons)):
        for j in range(i + 1, len(gluons)):
            dijet = gluons[i] + gluons[j]
            mass_diff = abs(dijet.M() - mHyp)
            if mass_diff < best_mass_diff:
                best_mass_diff = mass_diff
                best_pair = (gluons[i], gluons[j])

    if best_pair is None:
        continue

    g1, g2 = best_pair
    dijet = g1 + g2
    dr = deltaR(g1.Eta(), g1.Phi(), g2.Eta(), g2.Phi())

    if (dr < 0.4): continue
        
    # Fill histograms
    # ---------------
    h_gravPt.Fill(graviton_p4.Pt())
    h_dijetM.Fill(dijet.M())
    h_dR.Fill(dr)
    h_dRvsPt.Fill(graviton_p4.Pt(), dr)

    selected_pairs += 1

print(f"Finished loop: processed {evt_count} events, selected {selected_pairs} gluon pairs")

# -----------------------------
# --- Save ROOT histograms ---
# -----------------------------
outfile = f"HISTOS_MCSIG/histos_genLevel_GluonPairs_M{mHyp}_wAccCuts.root"
fout = ROOT.TFile(outfile, "RECREATE")
h_gravPt.Write()
h_dijetM.Write()
h_dR.Write()
h_dRvsPt.Write()
fout.Close()
print(f"Saved histograms to {outfile}")
