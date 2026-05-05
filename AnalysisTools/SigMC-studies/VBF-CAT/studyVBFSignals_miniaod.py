#!/usr/bin/env python3

"""
Minimal FWLite + PyROOT analyzer for generator-level Higgs studies in MiniAODSIM.

For each event, the script:
 - finds a generated Higgs boson (pdgId = 25)
 - stores the Higgs pt/eta/phi/mass
 - identifies its final non-self decay products
 - stores a two-body decay summary using the two daughters sorted by pt

The default input is the VBFHTo2B M300 sample.
"""

import sys
import glob
import math
import argparse
from array import array

import ROOT
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description="Generator-level Higgs study in MiniAOD")
parser.add_argument(
    "--input",
    type=str,
    default="/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/PrivateProduction/SampleFactory//VBFHTo2B_Par-M-200-W-0p014__chain_RunIII2024Summer24wmLHEGS-RunIII2024Summer24MiniAOD/SampleFactory/260416_155518/0000/*root",
    #default="/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/PrivateProduction/SampleFactory/VBFHTo2B_Par-M-300-W-0p014__chain_RunIII2024Summer24wmLHEGS-RunIII2024Summer24MiniAOD/SampleFactory/260416_160237/0000/*root",
    help="Input file pattern",
)
parser.add_argument(
    "--maxEvents",
    type=int,
    default=-1,
    help="Maximum number of events to process (-1 = all)",
)
parser.add_argument(
    "--out",
    type=str,
    default="VBFHTo2B_M300_genHiggs.root",
    help="Output ROOT file",
)
parser.add_argument(
    "--printEvery",
    type=int,
    default=1000,
    help="Print progress every N events",
)
args = parser.parse_args()


def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2.0 * math.pi
    while dphi <= -math.pi:
        dphi += 2.0 * math.pi
    return dphi


def delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = delta_phi(phi1, phi2)
    return math.sqrt(deta * deta + dphi * dphi)


def p4(obj):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(obj.pt(), obj.eta(), obj.phi(), obj.mass())
    return vec


def get_nonself_daughters(particle):
    daughters = []
    for idx in range(particle.numberOfDaughters()):
        dau = particle.daughter(idx)
        if dau is None:
            continue
        if dau.pdgId() == particle.pdgId():
            daughters.extend(get_nonself_daughters(dau))
        else:
            daughters.append(dau)
    return daughters


def find_last_higgs(pruned):
    higgs_candidates = []
    for gp in pruned:
        try:
            if abs(gp.pdgId()) != 25:
                continue
            daughters = get_nonself_daughters(gp)
            if len(daughters) < 2:
                continue
            higgs_candidates.append((gp, daughters))
        except Exception:
            continue

    if not higgs_candidates:
        return None, []

    higgs_candidates.sort(key=lambda item: item[0].pt(), reverse=True)
    return higgs_candidates[0]


files = sorted(glob.glob(args.input))
if not files:
    print(f"[ERROR] No files found for pattern:\n  {args.input}")
    sys.exit(1)

print(f"[INFO] Found {len(files)} files")
print(f"[INFO] First file: {files[0]}")

events = Events(files)
handle_pruned = Handle("std::vector<reco::GenParticle>")
label_pruned = ("prunedGenParticles", "")

fout = ROOT.TFile(args.out, "RECREATE")
tree = ROOT.TTree("Events", "Generator-level Higgs tree")

run = array("I", [0])
lumi = array("I", [0])
event_num = array("L", [0])
n_decay = array("i", [0])

higgs_pt = array("f", [0.0])
higgs_eta = array("f", [0.0])
higgs_phi = array("f", [0.0])
higgs_mass = array("f", [0.0])

dau1_pdgId = array("i", [0])
dau2_pdgId = array("i", [0])
dau_dR = array("f", [0.0])
dau_dEta = array("f", [0.0])
dau_dPhi = array("f", [0.0])
dau_mass = array("f", [0.0])

lead_pt = array("f", [0.0])
lead_eta = array("f", [0.0])
lead_phi = array("f", [0.0])
sublead_pt = array("f", [0.0])
sublead_eta = array("f", [0.0])
sublead_phi = array("f", [0.0])

tree.Branch("run", run, "run/i")
tree.Branch("lumi", lumi, "lumi/i")
tree.Branch("event", event_num, "event/l")
tree.Branch("n_decay_products", n_decay, "n_decay_products/I")

tree.Branch("higgs_pt", higgs_pt, "higgs_pt/F")
tree.Branch("higgs_eta", higgs_eta, "higgs_eta/F")
tree.Branch("higgs_phi", higgs_phi, "higgs_phi/F")
tree.Branch("higgs_mass", higgs_mass, "higgs_mass/F")

tree.Branch("decay1_pdgId", dau1_pdgId, "decay1_pdgId/I")
tree.Branch("decay2_pdgId", dau2_pdgId, "decay2_pdgId/I")
tree.Branch("decay_dR", dau_dR, "decay_dR/F")
tree.Branch("decay_dEta", dau_dEta, "decay_dEta/F")
tree.Branch("decay_dPhi", dau_dPhi, "decay_dPhi/F")
tree.Branch("decay_mass", dau_mass, "decay_mass/F")

tree.Branch("lead_decay_pt", lead_pt, "lead_decay_pt/F")
tree.Branch("lead_decay_eta", lead_eta, "lead_decay_eta/F")
tree.Branch("lead_decay_phi", lead_phi, "lead_decay_phi/F")
tree.Branch("sublead_decay_pt", sublead_pt, "sublead_decay_pt/F")
tree.Branch("sublead_decay_eta", sublead_eta, "sublead_decay_eta/F")
tree.Branch("sublead_decay_phi", sublead_phi, "sublead_decay_phi/F")


def reset_branches():
    n_decay[0] = 0

    higgs_pt[0] = -999.0
    higgs_eta[0] = -999.0
    higgs_phi[0] = -999.0
    higgs_mass[0] = -999.0

    dau1_pdgId[0] = 0
    dau2_pdgId[0] = 0
    dau_dR[0] = -999.0
    dau_dEta[0] = -999.0
    dau_dPhi[0] = -999.0
    dau_mass[0] = -999.0

    lead_pt[0] = -999.0
    lead_eta[0] = -999.0
    lead_phi[0] = -999.0
    sublead_pt[0] = -999.0
    sublead_eta[0] = -999.0
    sublead_phi[0] = -999.0


evt_count = 0
n_higgs_found = 0
n_two_body = 0

print("[INFO] Starting event loop")
for event in events:
    evt_count += 1
    if 0 < args.maxEvents <= evt_count:
        break
    if evt_count % args.printEvery == 0:
        print(f"[INFO] Processed {evt_count} events")

    aux = event.eventAuxiliary()
    run[0] = aux.run()
    lumi[0] = aux.luminosityBlock()
    event_num[0] = aux.event()
    reset_branches()

    event.getByLabel(label_pruned, handle_pruned)
    pruned = handle_pruned.product()

    higgs, daughters = find_last_higgs(pruned)
    if higgs is None:
        tree.Fill()
        continue

    n_higgs_found += 1
    n_decay[0] = len(daughters)
    higgs_pt[0] = higgs.pt()
    higgs_eta[0] = higgs.eta()
    higgs_phi[0] = higgs.phi()
    higgs_mass[0] = higgs.mass()

    if len(daughters) >= 2:
        daughters = sorted(daughters, key=lambda x: x.pt(), reverse=True)
        lead = daughters[0]
        sublead = daughters[1]

        n_two_body += 1
        dau1_pdgId[0] = lead.pdgId()
        dau2_pdgId[0] = sublead.pdgId()
        dau_dEta[0] = abs(lead.eta() - sublead.eta())
        dau_dPhi[0] = abs(delta_phi(lead.phi(), sublead.phi()))
        dau_dR[0] = delta_r(lead.eta(), lead.phi(), sublead.eta(), sublead.phi())
        dau_mass[0] = (p4(lead) + p4(sublead)).M()

        lead_pt[0] = lead.pt()
        lead_eta[0] = lead.eta()
        lead_phi[0] = lead.phi()
        sublead_pt[0] = sublead.pt()
        sublead_eta[0] = sublead.eta()
        sublead_phi[0] = sublead.phi()

    tree.Fill()

print("\n================ SUMMARY ================")
print(f"Processed events:         {evt_count}")
print(f"Events with a gen Higgs:  {n_higgs_found}")
print(f"Events with >=2 daughters:{n_two_body}")
print("========================================\n")

fout.cd()
tree.Write()
fout.Close()
print(f"[INFO] Saved ROOT file: {args.out}")
