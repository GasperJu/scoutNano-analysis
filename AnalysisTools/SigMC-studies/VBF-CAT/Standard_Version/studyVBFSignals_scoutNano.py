#!/usr/bin/env python3

"""
Minimal generator-level Higgs study for ScoutNano signal samples.

The script:
 - reads ScoutNano files with NanoAODTools
 - finds a generated Higgs boson (pdgId = 25)
 - stores the Higgs pt/eta/phi/mass
 - identifies its final non-self decay products
 - stores a two-body decay summary using the two daughters sorted by pt

Input can be provided either explicitly with --input or via --decay and --mass.
"""

import os
import math
import glob
import argparse
from array import array

import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

ROOT.gROOT.SetBatch(True)

BASE_DIR = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFHToXX"


parser = argparse.ArgumentParser(description="Generator-level ScoutNano Higgs study")
parser.add_argument(
    "--input",
    type=str,
    default=None,
    help="Input ScoutNano ROOT file or glob pattern. Overrides --decay/--mass.",
)
parser.add_argument(
    "--decay",
    type=str,
    default="2B",
    help="Signal decay label used in the directory/file naming, e.g. 2B, 2C, 2Q, 2Glu.",
)
parser.add_argument(
    "--mass",
    type=int,
    default=300,
    help="Signal mass hypothesis used in the file naming.",
)
parser.add_argument(
    "--output",
    type=str,
    default=None,
    help="Output ROOT file. If omitted, a name is built from decay and mass.",
)
parser.add_argument(
    "--maxEntries",
    type=int,
    default=None,
    help="Maximum number of entries to process.",
)
parser.add_argument(
    "--firstEntry",
    type=int,
    default=0,
    help="First entry to process.",
)
parser.add_argument(
    "--prefetch",
    action="store_true",
    help="Enable NanoAODTools prefetching.",
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


def p4_from_genpart(gp):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(gp.pt, gp.eta, gp.phi, gp.mass)
    return vec


def p4_from_genjet(gj):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(gj.pt, gj.eta, gj.phi, gj.mass)
    return vec


def get_input_files():
    if args.input:
        files = sorted(glob.glob(args.input))
        if files:
            return files
        if os.path.exists(args.input):
            return [args.input]
        raise RuntimeError(f"No files found for explicit input pattern: {args.input}")

    decay_dir = f"VBFHTo{args.decay}"
    file_name = f"hadd_VBFTo{args.decay}_M{args.mass}.root"
    candidate = os.path.join(BASE_DIR, decay_dir, file_name)
    if not os.path.exists(candidate):
        raise RuntimeError(f"Default sample not found: {candidate}")
    return [candidate]


def default_output_name():
    if args.output:
        return args.output
    return f"VBFHTo{args.decay}_M{args.mass}_genHiggs_scoutNano.root"


def final_nonself_daughters(genparts, higgs_idx):
    daughters = []
    queue = [higgs_idx]
    visited = set()

    while queue:
        idx = queue.pop()
        if idx in visited:
            continue
        visited.add(idx)

        parent = genparts[idx]
        child_indices = [
            i for i, gp in enumerate(genparts)
            if getattr(gp, "genPartIdxMother", -1) == idx
        ]

        for child_idx in child_indices:
            child = genparts[child_idx]
            if child.pdgId == parent.pdgId:
                queue.append(child_idx)
            else:
                daughters.append(child)

    return daughters


def find_last_higgs(genparts):
    higgs_indices = [i for i, gp in enumerate(genparts) if abs(gp.pdgId) == 25]
    if not higgs_indices:
        return None, []

    higgs_indices.sort(key=lambda idx: genparts[idx].pt, reverse=True)
    for idx in higgs_indices:
        daughters = final_nonself_daughters(genparts, idx)
        if len(daughters) >= 2:
            return genparts[idx], daughters

    best_idx = higgs_indices[0]
    return genparts[best_idx], final_nonself_daughters(genparts, best_idx)


def best_forward_genjet_pair(genjets):
    forward = [gj for gj in genjets if 3.0 <= abs(gj.eta) <= 5.0]
    if len(forward) < 2:
        return None, forward

    best_pair = None
    best_deta = -1.0
    for i in range(len(forward)):
        for j in range(i + 1, len(forward)):
            deta = abs(forward[i].eta - forward[j].eta)
            if deta > best_deta:
                best_deta = deta
                best_pair = (forward[i], forward[j])

    return best_pair, forward


class GenHiggsScoutNanoStudy(Module):
    def __init__(self, output_name):
        self.output_name = output_name
        self.writeHistFile = False
        self.total = 0
        self.n_higgs_found = 0
        self.n_two_body = 0

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.outfile = ROOT.TFile(self.output_name, "RECREATE")
        self.tree = ROOT.TTree("Events", "Generator-level Higgs tree from ScoutNano")

        self.run = array("I", [0])
        self.lumi = array("I", [0])
        self.event = array("L", [0])
        self.n_decay_products = array("i", [0])
        self.n_genjets_total = array("i", [0])
        self.n_genjets_central = array("i", [0])
        self.n_genjets_forward = array("i", [0])
        self.has_forward_genjet_pair = array("i", [0])

        self.higgs_pt = array("f", [0.0])
        self.higgs_eta = array("f", [0.0])
        self.higgs_phi = array("f", [0.0])
        self.higgs_mass = array("f", [0.0])

        self.decay1_pdgId = array("i", [0])
        self.decay2_pdgId = array("i", [0])
        self.decay_dR = array("f", [0.0])
        self.decay_dEta = array("f", [0.0])
        self.decay_dPhi = array("f", [0.0])
        self.decay_mass = array("f", [0.0])

        self.lead_decay_pt = array("f", [0.0])
        self.lead_decay_eta = array("f", [0.0])
        self.lead_decay_phi = array("f", [0.0])
        self.sublead_decay_pt = array("f", [0.0])
        self.sublead_decay_eta = array("f", [0.0])
        self.sublead_decay_phi = array("f", [0.0])

        self.forward_genjet1_pt = array("f", [0.0])
        self.forward_genjet1_eta = array("f", [0.0])
        self.forward_genjet1_phi = array("f", [0.0])
        self.forward_genjet1_mass = array("f", [0.0])
        self.forward_genjet2_pt = array("f", [0.0])
        self.forward_genjet2_eta = array("f", [0.0])
        self.forward_genjet2_phi = array("f", [0.0])
        self.forward_genjet2_mass = array("f", [0.0])
        self.forward_genjet_pair_dEta = array("f", [0.0])
        self.forward_genjet_pair_dPhi = array("f", [0.0])
        self.forward_genjet_pair_dR = array("f", [0.0])
        self.forward_genjet_pair_mass = array("f", [0.0])

        self.tree.Branch("run", self.run, "run/i")
        self.tree.Branch("lumi", self.lumi, "lumi/i")
        self.tree.Branch("event", self.event, "event/l")
        self.tree.Branch("n_decay_products", self.n_decay_products, "n_decay_products/I")
        self.tree.Branch("n_genjets_total", self.n_genjets_total, "n_genjets_total/I")
        self.tree.Branch("n_genjets_central", self.n_genjets_central, "n_genjets_central/I")
        self.tree.Branch("n_genjets_forward", self.n_genjets_forward, "n_genjets_forward/I")
        self.tree.Branch("has_forward_genjet_pair", self.has_forward_genjet_pair, "has_forward_genjet_pair/I")

        self.tree.Branch("higgs_pt", self.higgs_pt, "higgs_pt/F")
        self.tree.Branch("higgs_eta", self.higgs_eta, "higgs_eta/F")
        self.tree.Branch("higgs_phi", self.higgs_phi, "higgs_phi/F")
        self.tree.Branch("higgs_mass", self.higgs_mass, "higgs_mass/F")

        self.tree.Branch("decay1_pdgId", self.decay1_pdgId, "decay1_pdgId/I")
        self.tree.Branch("decay2_pdgId", self.decay2_pdgId, "decay2_pdgId/I")
        self.tree.Branch("decay_dR", self.decay_dR, "decay_dR/F")
        self.tree.Branch("decay_dEta", self.decay_dEta, "decay_dEta/F")
        self.tree.Branch("decay_dPhi", self.decay_dPhi, "decay_dPhi/F")
        self.tree.Branch("decay_mass", self.decay_mass, "decay_mass/F")

        self.tree.Branch("lead_decay_pt", self.lead_decay_pt, "lead_decay_pt/F")
        self.tree.Branch("lead_decay_eta", self.lead_decay_eta, "lead_decay_eta/F")
        self.tree.Branch("lead_decay_phi", self.lead_decay_phi, "lead_decay_phi/F")
        self.tree.Branch("sublead_decay_pt", self.sublead_decay_pt, "sublead_decay_pt/F")
        self.tree.Branch("sublead_decay_eta", self.sublead_decay_eta, "sublead_decay_eta/F")
        self.tree.Branch("sublead_decay_phi", self.sublead_decay_phi, "sublead_decay_phi/F")

        self.tree.Branch("forward_genjet1_pt", self.forward_genjet1_pt, "forward_genjet1_pt/F")
        self.tree.Branch("forward_genjet1_eta", self.forward_genjet1_eta, "forward_genjet1_eta/F")
        self.tree.Branch("forward_genjet1_phi", self.forward_genjet1_phi, "forward_genjet1_phi/F")
        self.tree.Branch("forward_genjet1_mass", self.forward_genjet1_mass, "forward_genjet1_mass/F")
        self.tree.Branch("forward_genjet2_pt", self.forward_genjet2_pt, "forward_genjet2_pt/F")
        self.tree.Branch("forward_genjet2_eta", self.forward_genjet2_eta, "forward_genjet2_eta/F")
        self.tree.Branch("forward_genjet2_phi", self.forward_genjet2_phi, "forward_genjet2_phi/F")
        self.tree.Branch("forward_genjet2_mass", self.forward_genjet2_mass, "forward_genjet2_mass/F")
        self.tree.Branch("forward_genjet_pair_dEta", self.forward_genjet_pair_dEta, "forward_genjet_pair_dEta/F")
        self.tree.Branch("forward_genjet_pair_dPhi", self.forward_genjet_pair_dPhi, "forward_genjet_pair_dPhi/F")
        self.tree.Branch("forward_genjet_pair_dR", self.forward_genjet_pair_dR, "forward_genjet_pair_dR/F")
        self.tree.Branch("forward_genjet_pair_mass", self.forward_genjet_pair_mass, "forward_genjet_pair_mass/F")

    def reset(self):
        self.n_decay_products[0] = 0
        self.n_genjets_total[0] = 0
        self.n_genjets_central[0] = 0
        self.n_genjets_forward[0] = 0
        self.has_forward_genjet_pair[0] = 0

        self.higgs_pt[0] = -999.0
        self.higgs_eta[0] = -999.0
        self.higgs_phi[0] = -999.0
        self.higgs_mass[0] = -999.0

        self.decay1_pdgId[0] = 0
        self.decay2_pdgId[0] = 0
        self.decay_dR[0] = -999.0
        self.decay_dEta[0] = -999.0
        self.decay_dPhi[0] = -999.0
        self.decay_mass[0] = -999.0

        self.lead_decay_pt[0] = -999.0
        self.lead_decay_eta[0] = -999.0
        self.lead_decay_phi[0] = -999.0
        self.sublead_decay_pt[0] = -999.0
        self.sublead_decay_eta[0] = -999.0
        self.sublead_decay_phi[0] = -999.0

        self.forward_genjet1_pt[0] = -999.0
        self.forward_genjet1_eta[0] = -999.0
        self.forward_genjet1_phi[0] = -999.0
        self.forward_genjet1_mass[0] = -999.0
        self.forward_genjet2_pt[0] = -999.0
        self.forward_genjet2_eta[0] = -999.0
        self.forward_genjet2_phi[0] = -999.0
        self.forward_genjet2_mass[0] = -999.0
        self.forward_genjet_pair_dEta[0] = -999.0
        self.forward_genjet_pair_dPhi[0] = -999.0
        self.forward_genjet_pair_dR[0] = -999.0
        self.forward_genjet_pair_mass[0] = -999.0

    def analyze(self, event):
        self.total += 1

        self.run[0] = getattr(event, "run", 0)
        self.lumi[0] = getattr(event, "luminosityBlock", 0)
        self.event[0] = getattr(event, "event", 0)
        self.reset()

        genparts = Collection(event, "GenPart")
        genjets = Collection(event, "GenJet")
        higgs, daughters = find_last_higgs(genparts)

        self.n_genjets_total[0] = len(genjets)
        self.n_genjets_central[0] = len([gj for gj in genjets if abs(gj.eta) < 3.0])
        pair, forward_genjets = best_forward_genjet_pair(genjets)
        self.n_genjets_forward[0] = len(forward_genjets)

        if higgs is not None:
            self.n_higgs_found += 1
            self.n_decay_products[0] = len(daughters)
            self.higgs_pt[0] = higgs.pt
            self.higgs_eta[0] = higgs.eta
            self.higgs_phi[0] = higgs.phi
            self.higgs_mass[0] = higgs.mass

        if higgs is not None and len(daughters) >= 2:
            daughters = sorted(daughters, key=lambda x: x.pt, reverse=True)
            lead = daughters[0]
            sublead = daughters[1]

            self.n_two_body += 1
            self.decay1_pdgId[0] = lead.pdgId
            self.decay2_pdgId[0] = sublead.pdgId
            self.decay_dEta[0] = abs(lead.eta - sublead.eta)
            self.decay_dPhi[0] = abs(delta_phi(lead.phi, sublead.phi))
            self.decay_dR[0] = delta_r(lead.eta, lead.phi, sublead.eta, sublead.phi)
            self.decay_mass[0] = (p4_from_genpart(lead) + p4_from_genpart(sublead)).M()

            self.lead_decay_pt[0] = lead.pt
            self.lead_decay_eta[0] = lead.eta
            self.lead_decay_phi[0] = lead.phi
            self.sublead_decay_pt[0] = sublead.pt
            self.sublead_decay_eta[0] = sublead.eta
            self.sublead_decay_phi[0] = sublead.phi

        if pair is not None:
            gj1, gj2 = pair
            self.has_forward_genjet_pair[0] = 1
            self.forward_genjet1_pt[0] = gj1.pt
            self.forward_genjet1_eta[0] = gj1.eta
            self.forward_genjet1_phi[0] = gj1.phi
            self.forward_genjet1_mass[0] = gj1.mass
            self.forward_genjet2_pt[0] = gj2.pt
            self.forward_genjet2_eta[0] = gj2.eta
            self.forward_genjet2_phi[0] = gj2.phi
            self.forward_genjet2_mass[0] = gj2.mass
            self.forward_genjet_pair_dEta[0] = abs(gj1.eta - gj2.eta)
            self.forward_genjet_pair_dPhi[0] = abs(delta_phi(gj1.phi, gj2.phi))
            self.forward_genjet_pair_dR[0] = delta_r(gj1.eta, gj1.phi, gj2.eta, gj2.phi)
            self.forward_genjet_pair_mass[0] = (p4_from_genjet(gj1) + p4_from_genjet(gj2)).M()

        self.tree.Fill()
        return True

    def endJob(self):
        print("\n================ ScoutNano Gen Higgs Study ===============")
        print(f"Total events:          {self.total}")
        print(f"Events with gen Higgs: {self.n_higgs_found}")
        print(f"Events with >=2 daughters: {self.n_two_body}")
        print("==========================================================\n")

        self.outfile.cd()
        self.tree.Write()
        self.outfile.Close()
        Module.endJob(self)


if __name__ == "__main__":
    input_files = get_input_files()
    output_name = default_output_name()

    print(f"[INFO] Using {len(input_files)} input file(s)")
    print(f"[INFO] First file: {input_files[0]}")
    print(f"[INFO] Output file: {output_name}")

    processor = PostProcessor(
        ".",
        input_files,
        modules=[GenHiggsScoutNanoStudy(output_name)],
        noOut=True,
        maxEntries=args.maxEntries,
        firstEntry=args.firstEntry,
        prefetch=args.prefetch,
    )
    processor.run()
