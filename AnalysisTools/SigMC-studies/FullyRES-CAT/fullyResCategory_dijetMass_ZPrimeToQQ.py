#!/usr/bin/env python3
import os
import math
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------
parser = argparse.ArgumentParser(description="Build resolved dijet mass histogram from ScoutNano ZPrimeToQQ signal")

parser.add_argument(
    "--mass",
    type=int,
    required=True,
    choices=[50, 100, 200, 500, 1000],
    help="Signal mass point"
)
parser.add_argument(
    "--baseDir",
    type=str,
    default="/eos/cms/store/cmst3/user/elfontan/scoutAna/samples/mc/2024/Signal/ZPrimeToQQ",
    help="Base directory containing ZprimeToQQ ROOT files"
)
parser.add_argument(
    "--tree",
    type=str,
    default="Events",
    help="Tree name"
)
parser.add_argument(
    "--xmin",
    type=float,
    default=0.0,
    help="Histogram minimum"
)
parser.add_argument(
    "--xmax",
    type=float,
    default=2000.0,
    help="Histogram maximum"
)
parser.add_argument(
    "--printEvery",
    type=int,
    default=5000,
    help="Print progress every N events"
)

args = parser.parse_args()

TREE_NAME = args.tree
XMIN = args.xmin
XMAX = args.xmax
NBINS = int(XMAX - XMIN)   # 1 GeV binning
PRINT_EVERY = args.printEvery

INPUT_FILE = os.path.join(args.baseDir, f"ZprimeToQQ_Par-M-{args.mass}.root")
OUTFILE = f"resolved_dijetMass_M{args.mass}.root"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def JetID_from_arrays(i, event):
    eta = event.ScoutingPFJet_eta[i]
    aeta = abs(eta)

    totalE = (
        event.ScoutingPFJet_neutralHadronEnergy[i]
        + event.ScoutingPFJet_HFHadronEnergy[i]
        + event.ScoutingPFJet_photonEnergy[i]
        + event.ScoutingPFJet_HFEMEnergy[i]
        + event.ScoutingPFJet_muonEnergy[i]
        + event.ScoutingPFJet_electronEnergy[i]
        + event.ScoutingPFJet_chargedHadronEnergy[i]
    )

    if totalE <= 0:
        return False

    NHF = (
        event.ScoutingPFJet_neutralHadronEnergy[i]
        + event.ScoutingPFJet_HFHadronEnergy[i]
    ) / float(totalE)

    NEMF = (
        event.ScoutingPFJet_photonEnergy[i]
        + event.ScoutingPFJet_HFEMEnergy[i]
    ) / float(totalE)

    muFrac = event.ScoutingPFJet_muonEnergy[i] / float(totalE)

    chargedMult = (
        event.ScoutingPFJet_chargedHadronMultiplicity[i]
        + event.ScoutingPFJet_HFHadronMultiplicity[i]
    )

    neutralMult = (
        event.ScoutingPFJet_neutralHadronMultiplicity[i]
        + event.ScoutingPFJet_HFEMMultiplicity[i]
    )

    nconst = (
        event.ScoutingPFJet_chargedHadronMultiplicity[i]
        + event.ScoutingPFJet_neutralHadronMultiplicity[i]
        + event.ScoutingPFJet_muonMultiplicity[i]
        + event.ScoutingPFJet_electronMultiplicity[i]
        + event.ScoutingPFJet_photonMultiplicity[i]
    )

    # -------- |eta| < 2.6 --------
    if aeta < 2.6:
        if NHF >= 0.99: return False
        if NEMF >= 0.90: return False
        if nconst <= 1: return False
        if chargedMult <= 0: return False
        if muFrac >= 0.80: return False
        return True

    # --- |eta| = [2.6,2.7] --------
    if 2.6 <= aeta < 2.7:
        if NEMF >= 0.99: return False
        if muFrac >= 0.80: return False
        return True

    # --- |eta| = [2.7,3.0] --------
    if 2.7 <= aeta < 3.0:
        if NEMF >= 0.99: return False
        if neutralMult <= 1: return False
        return True

    # --- |eta| = [3.0,5.0] --------
    if 3.0 <= aeta < 5.0:
        if NEMF >= 0.10: return False
        return True

    return False


def dijet_mass(j1, j2):
    pt1, eta1, phi1, m1 = j1
    pt2, eta2, phi2, m2 = j2

    px = pt1 * math.cos(phi1) + pt2 * math.cos(phi2)
    py = pt1 * math.sin(phi1) + pt2 * math.sin(phi2)
    pz = pt1 * math.sinh(eta1) + pt2 * math.sinh(eta2)

    e1 = math.sqrt((pt1 * math.cosh(eta1))**2 + m1*m1)
    e2 = math.sqrt((pt2 * math.cosh(eta2))**2 + m2*m2)
    E = e1 + e2

    m2jj = E*E - px*px - py*py - pz*pz
    return math.sqrt(max(m2jj, 0.0))

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    if not os.path.exists(INPUT_FILE):
        raise RuntimeError(f"Input file not found:\n  {INPUT_FILE}")

    f = ROOT.TFile.Open(INPUT_FILE)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open input file:\n  {INPUT_FILE}")

    tree = f.Get(TREE_NAME)
    if not tree:
        raise RuntimeError(f"Tree '{TREE_NAME}' not found in:\n  {INPUT_FILE}")

    nentries = tree.GetEntries()

    print(f"[INFO] Input file   : {INPUT_FILE}")
    print(f"[INFO] Output file  : {OUTFILE}")
    print(f"[INFO] Mass point   : {args.mass} GeV")
    print(f"[INFO] Entries      : {nentries}")
    print(f"[INFO] Histogram    : dijet_m_{args.mass}")
    print(f"[INFO] Binning      : {NBINS} bins from {XMIN} to {XMAX} GeV (1 GeV/bin)")

    hname = f"dijet_m_{args.mass}"
    htitle = f"dijet_m_{args.mass};m_{{jj}} [GeV];Events"
    h_mjj = ROOT.TH1F(hname, htitle, NBINS, XMIN, XMAX)
    h_mjj.Sumw2()

    selected_events = 0

    for iev, event in enumerate(tree):
        if iev % PRINT_EVERY == 0:
            print(f"[INFO] Processing event {iev}/{nentries}")

        jets = []
        njet = int(getattr(event, "nScoutingPFJet", 0))

        for i in range(njet):
            pt = float(event.ScoutingPFJet_pt[i])
            eta = float(event.ScoutingPFJet_eta[i])

            if pt <= 30.0:
                continue
            if abs(eta) >= 2.5:
                continue
            if not JetID_from_arrays(i, event):
                continue

            phi = float(event.ScoutingPFJet_phi[i])
            massj = float(event.ScoutingPFJet_m[i])

            jets.append((pt, eta, phi, massj))

        if len(jets) < 2:
            continue

        jets.sort(key=lambda x: x[0], reverse=True)
        j1 = jets[0]
        j2 = jets[1]

        # Require the two leading jets to be close in eta
        deta = abs(j1[1] - j2[1])
        if deta >= 1.3:
            continue

        mjj = dijet_mass(j1, j2)
        h_mjj.Fill(mjj)
        selected_events += 1

    print(f"[INFO] Selected events    : {selected_events}")
    print(f"[INFO] Histogram integral: {h_mjj.Integral()}")

    fout = ROOT.TFile(OUTFILE, "RECREATE")
    h_mjj.Write()
    fout.Close()
    f.Close()

    print(f"[INFO] Saved histogram to {OUTFILE}")


if __name__ == "__main__":
    main()
