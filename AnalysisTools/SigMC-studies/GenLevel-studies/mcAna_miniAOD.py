#!/usr/bin/env python3

"""
mcAna_miniAOD.py

PyROOT plus FWLite analyzer to:
 - read MINIAODSIM files
 - find generator gluons (pdgId=21) coming from mother with pdgId=5100039
 - match reconstructed jets (slimmedJets) to those gen gluons (by DeltaR)
 - compute dijet invariant mass of the two matched jets
 - plot histogram with CMS-like style and save canvas/png + root file

Usage:
  cmsenv
  python3 mcAna_miniAOD.py --mHyp 200
"""

import sys
import glob
import math
from array import array
import ROOT
from DataFormats.FWLite import Events, Handle
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

file_pattern = f"/eos/user/e/elfontan/SampleFactory/RSGravitonToGluonGluon_M{mHyp}_TuneCP5_13p6TeV_pythia8/chain_RunIII2024Summer24GS-RunIII2024Summer24MiniAOD/*/RunIII2024Summer24MiniAOD_*root"
maxEvents = -1          # -1 means process all
match_dr = 0.2          # DeltaR for jets <---> matching
hist_name = "h_dijet_m"
hist_bins = 100
hist_min = 0.0
hist_max = 500.0
print_every = 500

# ------------------------
# --- Helper functions --- 
# ------------------------
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


# --- Prepare file list ---
files = sorted(glob.glob(file_pattern))
if not files:
    print("No files found with pattern:", file_pattern)
    sys.exit(1)

print(f"Found {len(files)} files; first file: {files[0]}")

# --- Prepare FWLite handles ---
events = Events(files)

handle_jets = Handle("std::vector<pat::Jet>")
#label_jets = ("slimmedJets", "")
label_jets = ("slimmedJetsPuppi", "")

handle_pruned = Handle("std::vector<reco::GenParticle>")
label_pruned = ("prunedGenParticles", "")

# --------------------------
# --- Histogram
# --------------------------
ROOT.gROOT.SetBatch(False)  # set True to run batch (no GUI)
h_mjj = ROOT.TH1F("h_dijet_m", "Reco dijet mass (matched gluons from 5100039);m_{jj} (w/ FSR) [GeV];Events", 100, 0, 500)
h_mjj_preFSR = ROOT.TH1F("h_dijet_m_preFSR", "Reco dijet mass preFSR (matched gluons from 5100039);m_{jj} (w/o FSR) [GeV];Events", 100, 0, 500)
h_dRgg = ROOT.TH1F("h_dRgg", "#DeltaR between two gluons from 5100039;#DeltaR(g,g);Events", 60, 0, 6)
h_njetsTotal = ROOT.TH1F("h_njetsTotal", "Number of slimmedJets per event;N_{jets};Events", 20, -0.5, 19.5)
h_njetsAdditional = ROOT.TH1F("h_njetsAdditional", "Number of ISR candidates per event;N_{jets} ISR candidates;Events", 20, -0.5, 19.5)

# ISR candidate histograms (first 3 only)
h_ISR_pt  = [ROOT.TH1F(f"h_ISR{i}_pt", f"ISR jet {i} pT; pT [GeV]; Events", 50, 0, 500) for i in range(1, 4)]
h_ISR_eta = [ROOT.TH1F(f"h_ISR{i}_eta", f"ISR jet {i} #eta; #eta; Events", 50, -5, 5) for i in range(1, 4)]
h_ISR_phi = [ROOT.TH1F(f"h_ISR{i}_phi", f"ISR jet {i} #phi; #phi; Events", 64, -3.2, 3.2) for i in range(1, 4)]

evt_count = 0
filled = 0

print("----------------------")
print("Starting event loop...")
print("----------------------")
for event in events:
    evt_count += 1
    if 0 < maxEvents <= evt_count:
        break
    if evt_count % print_every == 0:
        print(f"Processed {evt_count} events, filled {filled} entries so far...")

    # --- jets ---
    event.getByLabel(label_jets, handle_jets)
    jets = handle_jets.product()
    h_njetsTotal.Fill(len(jets))  # fill N_jets histogram before selections

    # --- gen particles ---
    event.getByLabel(label_pruned, handle_pruned)
    pruned = handle_pruned.product()

    # collect gen gluons that come from 5100039
    target_gen_gluons = []
    for gp in pruned:
        try:
            if abs(gp.pdgId()) == 21:
                if find_mother_with_pdg(gp, 5100039):                    
                    p4 = ROOT.TLorentzVector()
                    p4.SetPtEtaPhiM(gp.pt(), gp.eta(), gp.phi(), gp.mass() if hasattr(gp, "mass") else 0.0)
                    target_gen_gluons.append( (gp, p4) )
        except Exception:
            continue

    if len(target_gen_gluons) >= 2:
        print("----------------------------------")
        print("Found two gen gluons from GRAVITON!")
        g0, g1 = target_gen_gluons[0][1], target_gen_gluons[1][1]
        dRgg = deltaR(g0.Eta(), g0.Phi(), g1.Eta(), g1.Phi())
        h_dRgg.Fill(dRgg)

    # --- reco jet matching as before ---
    reco_jets = []
    for j in jets:
        if (j.pt() < 72.0 or abs(j.eta()) > 5):
            #print("--- REJECTED! Acceptance cuts not passed.")
            continue
        p4j = ROOT.TLorentzVector()
        p4j.SetPtEtaPhiM(j.pt(), j.eta(), j.phi(), j.mass())
        reco_jets.append((j, p4j))

    matches = []
    for ig, (genobj, genp4) in enumerate(target_gen_gluons):
        best_dr, best_j = 999.0, -1
        for ij, (jobj, jp4) in enumerate(reco_jets):
            dr = deltaR(genp4.Eta(), genp4.Phi(), jp4.Eta(), jp4.Phi())
            if dr < best_dr:
                best_dr, best_j = dr, ij
        if best_j >= 0 and best_dr <= match_dr:
            matches.append((ig, best_j, best_dr))

    if len(matches) < 2:
        #print("--- REJECTED! Less than two matches.")
        continue

    # Build reco --> gen gluon mapping (for FSR checks)
    # -------------------------------------------------
    reco_to_gen = {}
    for ij, (jobj, jp4) in enumerate(reco_jets):
        best_dr, best_gp = 999.0, None
        for gp, genp4 in target_gen_gluons:
            dr = deltaR(jp4.Eta(), jp4.Phi(), genp4.Eta(), genp4.Phi())
            if dr < best_dr:
                best_dr, best_gp = dr, gp
        if best_gp and best_dr < 0.3:  # tighter matching threshold
            reco_to_gen[ij] = best_gp
    
    jet_to_best = {}
    for (ig, ij, dr) in matches:
        if ij not in jet_to_best or dr < jet_to_best[ij][1]:
            jet_to_best[ij] = (ig, dr)

    if len(jet_to_best) < 2:
        #print("--- REJECTED! jet_to_best < 2.")
        continue

    sorted_jets = sorted(jet_to_best.items(), key=lambda kv: kv[1][1])
    jet_indices = [sorted_jets[0][0], sorted_jets[1][0]]

    # --- Dijet mass before FSR recovery ---
    p4_0 = reco_jets[jet_indices[0]][1]
    p4_1 = reco_jets[jet_indices[1]][1]
    mjj_before = (p4_0 + p4_1).M()
    h_mjj_preFSR.Fill(mjj_before)
    
    # ----------------------
    # ---- FSR recovery ----
    # ----------------------
    matched_jets = [reco_jets[i][1].Clone() for i in jet_indices]  
    for i, mj in enumerate(matched_jets):
        for ij, (jobj, jp4) in enumerate(reco_jets):
            if ij in jet_indices:
                continue  
            dr = deltaR(mj.Eta(), mj.Phi(), jp4.Eta(), jp4.Phi())
            if dr < 0.8:
                is_from_graviton = False
                if ij in reco_to_gen:
                    gp = reco_to_gen[ij]
                    is_from_graviton = find_mother_with_pdg(gp, 5100039)
                if not is_from_graviton:
                    print(f"[DEBUG] Rejected FSR candidate: dR={dr:.3f}, pT={jp4.Pt():.1f}, not from Graviton gluon")
                    continue
                
                print(f"[DEBUG] Adding FSR jet to matched jet {i}: dR={dr:.3f}, FSR jet pT={jp4.Pt():.1f}")
                mj += jp4
                print(f"[DEBUG] Updated matched jet {i} pT after FSR = {mj.Pt():.1f}")
                
    # --- Dijet mass after FSR ---
    mjj_after = (matched_jets[0] + matched_jets[1]).M()
    print(f"[DEBUG] Dijet mass before FSR = {mjj_before:.1f} GeV, after FSR = {mjj_after:.1f} GeV")
    h_mjj.Fill(mjj_after)
    filled += 1

    # --- Remaining jets (not matched, not used in FSR) ---
    remaining = []
    for ij, (jobj, jp4) in enumerate(reco_jets):
        if ij not in jet_indices:
            # check it wasn't used in FSR
            used_in_fsr = False
            for mj in matched_jets:
                if deltaR(mj.Eta(), mj.Phi(), jp4.Eta(), jp4.Phi()) < 0.8:
                    used_in_fsr = True
                    break
            if not used_in_fsr:
                remaining.append(jp4)

    # sort by pt
    remaining_sorted = sorted(remaining, key=lambda v: v.Pt(), reverse=True)
    h_njetsAdditional.Fill(len(remaining_sorted))

    if remaining_sorted:
        print(f"[DEBUG] ISR candidates (sorted by pT): {[f'{v.Pt():.1f}' for v in remaining_sorted]}")

    # fill ISR histos for first 3 jets
    for i, v in enumerate(remaining_sorted[:3]):
        h_ISR_pt[i].Fill(v.Pt())
        h_ISR_eta[i].Fill(v.Eta())
        h_ISR_phi[i].Fill(v.Phi())

# --- done loop ---
print(f"Event loop finished: processed {evt_count} events, filled {filled} dijet entries.")


# Save outputs
outfile_root = f"recoDijetMass_RSGraviton_M{mHyp}_wPUPPIJets_FSR0p8_jetPt72Eta5.root"
fout = ROOT.TFile(outfile_root, "RECREATE")
h_mjj.Write()
h_mjj_preFSR.Write()
h_dRgg.Write()
h_njetsTotal.Write()
h_njetsAdditional.Write()
for h in h_ISR_pt + h_ISR_eta + h_ISR_phi:
    h.Write()
fout.Close()

print(f"Saved histogram to {outfile_root}")
print("Done.")
