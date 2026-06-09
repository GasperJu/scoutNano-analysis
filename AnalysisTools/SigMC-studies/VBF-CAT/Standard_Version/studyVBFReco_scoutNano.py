# ---------------------------------------------------------
# ----- Run 3 scouting dijet analysis: signal studies -----
# ---------------------------------------------------------
# Optimising VBF selection
# ---------------------------------------------------------

# --------
# HowToRun
# --------
#python3 studyVBFReco_scoutNano.py --decay 2B --mass 300 --requireTriggerBaseline


#!/usr/bin/env python3

"""
Minimal reco-level VBF signal study for ScoutNano samples.

The script:
 - reads ScoutNano files with NanoAODTools
 - stores trigger bits for JetHT and hadronic L1 seeds
 - identifies a forward reco-jet pair to tag VBF
 - studies the leading/subleading central reco-jet pair
 - stores reco kinematics with and without gen matching to Higgs decay products

Input can be provided either explicitly with --input or via --decay and --mass.
"""

import os
import math
import glob
import argparse
from array import array
from types import SimpleNamespace

import ROOT
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

ROOT.gROOT.SetBatch(True)

BASE_DIR = "/eos/cms/store/cmst3/group/run3Scouting/LowMassDijetSearch/samples/2024/mc/Signals-ScoutNano/VBFHToXX"
MASS_DEPENDENT_PAIR_WP = [
    (300.0, 1.5, 0.60),
    (500.0, 2.0, 0.45),
    (1000.0, 2.5, 0.35),
]


parser = argparse.ArgumentParser(description="Reco-level ScoutNano VBF signal study")
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
parser.add_argument(
    "--jetPtMin",
    type=float,
    default=30.0,
    help="Minimum reco-jet pT.",
)
parser.add_argument(
    "--centralEtaMax",
    type=float,
    default=3.0,
    help="Maximum |eta| for central reco jets.",
)
parser.add_argument(
    "--forwardEtaMin",
    type=float,
    default=3.0,
    help="Minimum |eta| for forward reco jets.",
)
parser.add_argument(
    "--forwardEtaMax",
    type=float,
    default=5.0,
    help="Maximum |eta| for forward reco jets.",
)
parser.add_argument(
    "--forwardPairDEtaMin",
    type=float,
    default=5.0,
    help="Minimum |DeltaEta| for the forward reco-jet pair.",
)
parser.add_argument(
    "--centralPairDEtaMax",
    type=float,
    default=1.3,
    help="Maximum |DeltaEta| required for the central reco-jet pair.",
)
parser.add_argument(
    "--centralPairDPhiMin",
    type=float,
    default=1.5,
    help="Minimum |DeltaPhi| required for the central reco-jet pair.",
)
parser.add_argument(
    "--centralPairSelection",
    type=str,
    default="massDependent",
    choices=["massDependent", "angular"],
    help="Central-pair acceptance mode: mass-dependent dR + pT imbalance or legacy angular cuts.",
)
parser.add_argument(
    "--centralPairDRMin",
    type=float,
    default=None,
    help="Override minimum DeltaR required for the central reco-jet pair in massDependent mode.",
)
parser.add_argument(
    "--centralPairPtImbalanceMax",
    type=float,
    default=None,
    help="Override maximum pT imbalance required for the central reco-jet pair in massDependent mode.",
)
parser.add_argument(
    "--centralPairAlgo",
    type=str,
    default="leading",
    choices=["leading", "minDR", "minDPhi", "maxPtSum", "massClosest"],
    help="Algorithm used to choose the central jet pair.",
)
parser.add_argument(
    "--matchDR",
    type=float,
    default=0.4,
    help="DeltaR matching threshold between reco central jets and Higgs decay products.",
)
parser.add_argument(
    "--requireTriggerBaseline",
    action="store_true",
    help="Require DST_PFScouting_JetHT and (L1_HTT280er or L1_SingleJet180).",
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


def interpolate_wp(signal_mass, x0, y0, x1, y1):
    if x1 == x0:
        return y0
    frac = (signal_mass - x0) / float(x1 - x0)
    return y0 + frac * (y1 - y0)


def mass_dependent_pair_wp(signal_mass):
    if args.centralPairDRMin is not None and args.centralPairPtImbalanceMax is not None:
        return args.centralPairDRMin, args.centralPairPtImbalanceMax

    anchors = MASS_DEPENDENT_PAIR_WP
    if signal_mass <= anchors[0][0]:
        dr_min = anchors[0][1]
        apt_max = anchors[0][2]
    elif signal_mass >= anchors[-1][0]:
        dr_min = anchors[-1][1]
        apt_max = anchors[-1][2]
    else:
        dr_min = None
        apt_max = None
        for idx in range(len(anchors) - 1):
            m0, dr0, apt0 = anchors[idx]
            m1, dr1, apt1 = anchors[idx + 1]
            if m0 <= signal_mass <= m1:
                dr_min = interpolate_wp(signal_mass, m0, dr0, m1, dr1)
                apt_max = interpolate_wp(signal_mass, m0, apt0, m1, apt1)
                break

    if args.centralPairDRMin is not None:
        dr_min = args.centralPairDRMin
    if args.centralPairPtImbalanceMax is not None:
        apt_max = args.centralPairPtImbalanceMax
    return dr_min, apt_max


def pair_pt_imbalance(j1, j2):
    ptsum = j1.pt + j2.pt
    if ptsum <= 0.0:
        return 999.0
    return abs(j1.pt - j2.pt) / float(ptsum)


def p4_from_recojet(jet):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, recojet_mass(jet))
    return vec


def recojet_mass(jet):
    for attr in ["mass", "m"]:
        try:
            return float(getattr(jet, attr))
        except Exception:
            continue
    return 0.0


def corrected_recojet(jet, pt_corr):
    return SimpleNamespace(
        original=jet,
        pt=pt_corr,
        eta=jet.eta,
        phi=jet.phi,
        mass=recojet_mass(jet),
        neutralHadronEnergy=jet.neutralHadronEnergy,
        HFHadronEnergy=jet.HFHadronEnergy,
        photonEnergy=jet.photonEnergy,
        HFEMEnergy=jet.HFEMEnergy,
        muonEnergy=jet.muonEnergy,
        electronEnergy=jet.electronEnergy,
        chargedHadronEnergy=jet.chargedHadronEnergy,
        chargedHadronMultiplicity=jet.chargedHadronMultiplicity,
        HFHadronMultiplicity=jet.HFHadronMultiplicity,
        neutralHadronMultiplicity=jet.neutralHadronMultiplicity,
        HFEMMultiplicity=jet.HFEMMultiplicity,
        muonMultiplicity=jet.muonMultiplicity,
        electronMultiplicity=jet.electronMultiplicity,
        photonMultiplicity=jet.photonMultiplicity,
    )


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
    return f"VBFHTo{args.decay}_M{args.mass}_recoVBF_scoutNano.root"


def get_trigger_bits(event):
    pass_dst = bool(getattr(event, "DST_PFScouting_JetHT", False))
    pass_htt280 = bool(getattr(event, "L1_HTT280er", False))
    pass_singlejet180 = bool(getattr(event, "L1_SingleJet180", False))
    pass_baseline = pass_dst and (pass_htt280 or pass_singlejet180)
    return {
        "pass_dst": pass_dst,
        "pass_htt280": pass_htt280,
        "pass_singlejet180": pass_singlejet180,
        "pass_baseline": pass_baseline,
    }


def jet_id(jet):
    aeta = abs(jet.eta)

    total_e = (
        jet.neutralHadronEnergy + jet.HFHadronEnergy +
        jet.photonEnergy + jet.HFEMEnergy +
        jet.muonEnergy + jet.electronEnergy +
        jet.chargedHadronEnergy
    )
    if total_e <= 0:
        return False

    nhf = (jet.neutralHadronEnergy + jet.HFHadronEnergy) / float(total_e)
    nemf = (jet.photonEnergy + jet.HFEMEnergy) / float(total_e)
    mu_frac = jet.muonEnergy / float(total_e)

    charged_mult = jet.chargedHadronMultiplicity + jet.HFHadronMultiplicity
    neutral_mult = jet.neutralHadronMultiplicity + jet.HFEMMultiplicity
    nconst = (
        jet.chargedHadronMultiplicity + jet.neutralHadronMultiplicity +
        jet.muonMultiplicity + jet.electronMultiplicity +
        jet.photonMultiplicity
    )

    if aeta < 2.6:
        if nhf >= 0.99 or nemf >= 0.90 or nconst <= 1:
            return False
        if charged_mult <= 0 or mu_frac >= 0.80:
            return False
        return True

    if 2.6 <= aeta < 2.7:
        if nemf >= 0.99 or mu_frac >= 0.80:
            return False
        return True

    if 2.7 <= aeta < 3.0:
        if nemf >= 0.99 or neutral_mult <= 1:
            return False
        return True

    if 3.0 <= aeta < 5.0:
        if nemf >= 0.10:
            return False
        return True

    return False


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


def best_forward_pair(jets):
    if len(jets) < 2:
        return None

    best_pair = None
    best_deta = -1.0
    for i in range(len(jets)):
        for j in range(i + 1, len(jets)):
            deta = abs(jets[i].eta - jets[j].eta)
            if deta > best_deta:
                best_deta = deta
                best_pair = (jets[i], jets[j])
    return best_pair


def choose_central_pair(jets, algo, max_deta=None, min_dphi=None, target_mass=None):
    if len(jets) < 2:
        return None

    candidates = []
    for i in range(len(jets)):
        for j in range(i + 1, len(jets)):
            j1, j2 = jets[i], jets[j]
            deta = abs(j1.eta - j2.eta)
            if max_deta is not None and deta >= max_deta:
                continue
            dr = delta_r(j1.eta, j1.phi, j2.eta, j2.phi)
            dphi = abs(delta_phi(j1.phi, j2.phi))
            if min_dphi is not None and dphi <= min_dphi:
                continue
            ptsum = j1.pt + j2.pt
            mass = (p4_from_recojet(j1) + p4_from_recojet(j2)).M()
            candidates.append((j1, j2, deta, dr, dphi, ptsum, mass))

    if not candidates:
        return None

    if algo == "leading":
        best = candidates[0]
    elif algo == "minDR":
        best = min(candidates, key=lambda item: item[3])
    elif algo == "minDPhi":
        best = min(candidates, key=lambda item: item[4])
    elif algo == "maxPtSum":
        best = max(candidates, key=lambda item: item[5])
    elif algo == "massClosest":
        if target_mass is None:
            target_mass = 0.0
        best = min(candidates, key=lambda item: abs(item[6] - target_mass))
    else:
        raise ValueError(f"Unknown centralPairAlgo: {algo}")

    return best[0], best[1]


def central_pair_passes_selection(j1, j2):
    dR = delta_r(j1.eta, j1.phi, j2.eta, j2.phi)
    dEta = abs(j1.eta - j2.eta)
    dPhi = abs(delta_phi(j1.phi, j2.phi))
    apt = pair_pt_imbalance(j1, j2)

    if args.centralPairSelection == "angular":
        passes = dEta < args.centralPairDEtaMax and dPhi > args.centralPairDPhiMin
    elif args.centralPairSelection == "massDependent":
        dr_min, apt_max = mass_dependent_pair_wp(float(args.mass))
        passes = dR > dr_min and apt < apt_max
    else:
        raise ValueError(f"Unknown centralPairSelection: {args.centralPairSelection}")

    return passes, dR, dEta, dPhi, apt


def best_match_distances(j1, j2, daughters):
    if len(daughters) < 2:
        return False, -999.0, -999.0

    daughters = sorted(daughters, key=lambda x: x.pt, reverse=True)
    d1, d2 = daughters[0], daughters[1]

    dr11 = delta_r(j1.eta, j1.phi, d1.eta, d1.phi)
    dr22 = delta_r(j2.eta, j2.phi, d2.eta, d2.phi)
    dr12 = delta_r(j1.eta, j1.phi, d2.eta, d2.phi)
    dr21 = delta_r(j2.eta, j2.phi, d1.eta, d1.phi)

    if max(dr11, dr22) <= max(dr12, dr21):
        return True, dr11, dr22
    return True, dr12, dr21


class RecoVBFScoutNanoStudy(Module):
    def __init__(self, output_name):
        self.output_name = output_name
        self.writeHistFile = False
        self.total = 0
        self.after_trigger = 0
        self.n_forward_pair = 0
        self.n_vbf_tag = 0
        self.n_central_pair = 0
        self.n_genmatched = 0

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        jec_json = os.path.join(
            os.getenv("CMSSW_BASE", ""),
            "src",
            "2024_UtilsDataQuality",
            "jetHLT_jerc.json"
        )
        if not os.path.exists(jec_json):
            raise RuntimeError(f"HLT JEC JSON not found: {jec_json}")

        cset = correctionlib.CorrectionSet.from_file(jec_json)
        self.jec = cset.compound["HLT_Winter24_V1_MC_L1L2L3Res_AK4PFHLT"]

        self.outfile = ROOT.TFile(self.output_name, "RECREATE")
        self.tree = ROOT.TTree("Events", "Reco-level VBF study tree from ScoutNano")

        self.run = array("I", [0])
        self.lumi = array("I", [0])
        self.event = array("L", [0])

        self.pass_dst_jetht = array("i", [0])
        self.pass_l1_htt280 = array("i", [0])
        self.pass_l1_singlejet180 = array("i", [0])
        self.pass_trigger_baseline = array("i", [0])

        self.n_recojets_total_raw = array("i", [0])
        self.n_recojets_central_raw = array("i", [0])
        self.n_recojets_forward_raw = array("i", [0])
        self.n_recojets_total = array("i", [0])
        self.n_recojets_central = array("i", [0])
        self.n_recojets_forward = array("i", [0])

        self.has_forward_pair = array("i", [0])
        self.has_vbf_tag = array("i", [0])
        self.has_central_pair = array("i", [0])
        self.has_gen_higgs = array("i", [0])
        self.has_two_decay_products = array("i", [0])
        self.central_pair_genmatched = array("i", [0])

        self.forward_jet1_pt = array("f", [0.0])
        self.forward_jet1_eta = array("f", [0.0])
        self.forward_jet1_phi = array("f", [0.0])
        self.forward_jet1_mass = array("f", [0.0])
        self.forward_jet2_pt = array("f", [0.0])
        self.forward_jet2_eta = array("f", [0.0])
        self.forward_jet2_phi = array("f", [0.0])
        self.forward_jet2_mass = array("f", [0.0])
        self.forward_pair_dEta = array("f", [0.0])
        self.forward_pair_dPhi = array("f", [0.0])
        self.forward_pair_dR = array("f", [0.0])
        self.forward_pair_mass = array("f", [0.0])

        self.central_jet1_pt = array("f", [0.0])
        self.central_jet1_eta = array("f", [0.0])
        self.central_jet1_phi = array("f", [0.0])
        self.central_jet1_mass = array("f", [0.0])
        self.central_jet2_pt = array("f", [0.0])
        self.central_jet2_eta = array("f", [0.0])
        self.central_jet2_phi = array("f", [0.0])
        self.central_jet2_mass = array("f", [0.0])
        self.central_pair_dEta = array("f", [0.0])
        self.central_pair_dPhi = array("f", [0.0])
        self.central_pair_dR = array("f", [0.0])
        self.central_pair_mass = array("f", [0.0])
        self.central_pair_match_dR1 = array("f", [0.0])
        self.central_pair_match_dR2 = array("f", [0.0])

        self.forward_jet1_pt_raw = array("f", [0.0])
        self.forward_jet2_pt_raw = array("f", [0.0])
        self.central_jet1_pt_raw = array("f", [0.0])
        self.central_jet2_pt_raw = array("f", [0.0])

        self.tree.Branch("run", self.run, "run/i")
        self.tree.Branch("lumi", self.lumi, "lumi/i")
        self.tree.Branch("event", self.event, "event/l")
        self.tree.Branch("pass_dst_jetht", self.pass_dst_jetht, "pass_dst_jetht/I")
        self.tree.Branch("pass_l1_htt280", self.pass_l1_htt280, "pass_l1_htt280/I")
        self.tree.Branch("pass_l1_singlejet180", self.pass_l1_singlejet180, "pass_l1_singlejet180/I")
        self.tree.Branch("pass_trigger_baseline", self.pass_trigger_baseline, "pass_trigger_baseline/I")

        self.tree.Branch("n_recojets_total_raw", self.n_recojets_total_raw, "n_recojets_total_raw/I")
        self.tree.Branch("n_recojets_central_raw", self.n_recojets_central_raw, "n_recojets_central_raw/I")
        self.tree.Branch("n_recojets_forward_raw", self.n_recojets_forward_raw, "n_recojets_forward_raw/I")
        self.tree.Branch("n_recojets_total", self.n_recojets_total, "n_recojets_total/I")
        self.tree.Branch("n_recojets_central", self.n_recojets_central, "n_recojets_central/I")
        self.tree.Branch("n_recojets_forward", self.n_recojets_forward, "n_recojets_forward/I")

        self.tree.Branch("has_forward_pair", self.has_forward_pair, "has_forward_pair/I")
        self.tree.Branch("has_vbf_tag", self.has_vbf_tag, "has_vbf_tag/I")
        self.tree.Branch("has_central_pair", self.has_central_pair, "has_central_pair/I")
        self.tree.Branch("has_gen_higgs", self.has_gen_higgs, "has_gen_higgs/I")
        self.tree.Branch("has_two_decay_products", self.has_two_decay_products, "has_two_decay_products/I")
        self.tree.Branch("central_pair_genmatched", self.central_pair_genmatched, "central_pair_genmatched/I")

        self.tree.Branch("forward_jet1_pt", self.forward_jet1_pt, "forward_jet1_pt/F")
        self.tree.Branch("forward_jet1_pt_raw", self.forward_jet1_pt_raw, "forward_jet1_pt_raw/F")
        self.tree.Branch("forward_jet1_eta", self.forward_jet1_eta, "forward_jet1_eta/F")
        self.tree.Branch("forward_jet1_phi", self.forward_jet1_phi, "forward_jet1_phi/F")
        self.tree.Branch("forward_jet1_mass", self.forward_jet1_mass, "forward_jet1_mass/F")
        self.tree.Branch("forward_jet2_pt", self.forward_jet2_pt, "forward_jet2_pt/F")
        self.tree.Branch("forward_jet2_pt_raw", self.forward_jet2_pt_raw, "forward_jet2_pt_raw/F")
        self.tree.Branch("forward_jet2_eta", self.forward_jet2_eta, "forward_jet2_eta/F")
        self.tree.Branch("forward_jet2_phi", self.forward_jet2_phi, "forward_jet2_phi/F")
        self.tree.Branch("forward_jet2_mass", self.forward_jet2_mass, "forward_jet2_mass/F")
        self.tree.Branch("forward_pair_dEta", self.forward_pair_dEta, "forward_pair_dEta/F")
        self.tree.Branch("forward_pair_dPhi", self.forward_pair_dPhi, "forward_pair_dPhi/F")
        self.tree.Branch("forward_pair_dR", self.forward_pair_dR, "forward_pair_dR/F")
        self.tree.Branch("forward_pair_mass", self.forward_pair_mass, "forward_pair_mass/F")

        self.tree.Branch("central_jet1_pt", self.central_jet1_pt, "central_jet1_pt/F")
        self.tree.Branch("central_jet1_pt_raw", self.central_jet1_pt_raw, "central_jet1_pt_raw/F")
        self.tree.Branch("central_jet1_eta", self.central_jet1_eta, "central_jet1_eta/F")
        self.tree.Branch("central_jet1_phi", self.central_jet1_phi, "central_jet1_phi/F")
        self.tree.Branch("central_jet1_mass", self.central_jet1_mass, "central_jet1_mass/F")
        self.tree.Branch("central_jet2_pt", self.central_jet2_pt, "central_jet2_pt/F")
        self.tree.Branch("central_jet2_pt_raw", self.central_jet2_pt_raw, "central_jet2_pt_raw/F")
        self.tree.Branch("central_jet2_eta", self.central_jet2_eta, "central_jet2_eta/F")
        self.tree.Branch("central_jet2_phi", self.central_jet2_phi, "central_jet2_phi/F")
        self.tree.Branch("central_jet2_mass", self.central_jet2_mass, "central_jet2_mass/F")
        self.tree.Branch("central_pair_dEta", self.central_pair_dEta, "central_pair_dEta/F")
        self.tree.Branch("central_pair_dPhi", self.central_pair_dPhi, "central_pair_dPhi/F")
        self.tree.Branch("central_pair_dR", self.central_pair_dR, "central_pair_dR/F")
        self.tree.Branch("central_pair_mass", self.central_pair_mass, "central_pair_mass/F")
        self.tree.Branch("central_pair_match_dR1", self.central_pair_match_dR1, "central_pair_match_dR1/F")
        self.tree.Branch("central_pair_match_dR2", self.central_pair_match_dR2, "central_pair_match_dR2/F")

    def reset(self):
        self.pass_dst_jetht[0] = 0
        self.pass_l1_htt280[0] = 0
        self.pass_l1_singlejet180[0] = 0
        self.pass_trigger_baseline[0] = 0

        self.n_recojets_total_raw[0] = 0
        self.n_recojets_central_raw[0] = 0
        self.n_recojets_forward_raw[0] = 0
        self.n_recojets_total[0] = 0
        self.n_recojets_central[0] = 0
        self.n_recojets_forward[0] = 0

        self.has_forward_pair[0] = 0
        self.has_vbf_tag[0] = 0
        self.has_central_pair[0] = 0
        self.has_gen_higgs[0] = 0
        self.has_two_decay_products[0] = 0
        self.central_pair_genmatched[0] = 0

        for branch in [
            self.forward_jet1_pt, self.forward_jet1_eta, self.forward_jet1_phi, self.forward_jet1_mass,
            self.forward_jet2_pt, self.forward_jet2_eta, self.forward_jet2_phi, self.forward_jet2_mass,
            self.forward_pair_dEta, self.forward_pair_dPhi, self.forward_pair_dR, self.forward_pair_mass,
            self.central_jet1_pt, self.central_jet1_eta, self.central_jet1_phi, self.central_jet1_mass,
            self.central_jet2_pt, self.central_jet2_eta, self.central_jet2_phi, self.central_jet2_mass,
            self.central_pair_dEta, self.central_pair_dPhi, self.central_pair_dR, self.central_pair_mass,
            self.central_pair_match_dR1, self.central_pair_match_dR2,
            self.forward_jet1_pt_raw, self.forward_jet2_pt_raw,
            self.central_jet1_pt_raw, self.central_jet2_pt_raw,
        ]:
            branch[0] = -999.0

    def corrected_jets(self, event, jets_raw):
        rho = float(getattr(event, "ScoutingRho_fixedGridRhoFastjetAll", 0.0))
        run = int(getattr(event, "run", 1))

        corrected = []
        for jet in jets_raw:
            variables = {
                "JetPt": jet.pt,
                "JetEta": jet.eta,
                "JetPhi": jet.phi,
                "JetA": getattr(jet, "jetArea", 0.0),
                "Rho": rho,
                "run": run,
            }
            inputs = [variables[inp.name] for inp in self.jec.inputs]
            corr_factor = self.jec.evaluate(*inputs)
            corrected.append(corrected_recojet(jet, jet.pt * corr_factor))
        return corrected

    def analyze(self, event):
        self.total += 1

        self.run[0] = getattr(event, "run", 0)
        self.lumi[0] = getattr(event, "luminosityBlock", 0)
        self.event[0] = getattr(event, "event", 0)
        self.reset()

        trigger_bits = get_trigger_bits(event)
        self.pass_dst_jetht[0] = int(trigger_bits["pass_dst"])
        self.pass_l1_htt280[0] = int(trigger_bits["pass_htt280"])
        self.pass_l1_singlejet180[0] = int(trigger_bits["pass_singlejet180"])
        self.pass_trigger_baseline[0] = int(trigger_bits["pass_baseline"])

        if args.requireTriggerBaseline and not trigger_bits["pass_baseline"]:
            self.tree.Fill()
            return True

        if trigger_bits["pass_baseline"]:
            self.after_trigger += 1

        genparts = Collection(event, "GenPart")
        jets_raw = Collection(event, "ScoutingPFJet")
        jets = self.corrected_jets(event, jets_raw)

        higgs, daughters = find_last_higgs(genparts)
        if higgs is not None:
            self.has_gen_higgs[0] = 1
        if len(daughters) >= 2:
            self.has_two_decay_products[0] = 1

        good_jets_raw = [
            jet for jet in jets_raw
            if jet_id(jet) and jet.pt > args.jetPtMin and abs(jet.eta) < args.forwardEtaMax
        ]
        forward_jets_raw = [
            jet for jet in good_jets_raw
            if args.forwardEtaMin <= abs(jet.eta) <= args.forwardEtaMax
        ]
        central_jets_raw = [
            jet for jet in good_jets_raw
            if abs(jet.eta) < args.centralEtaMax
        ]
        central_jets_raw = sorted(central_jets_raw, key=lambda x: x.pt, reverse=True)

        good_jets = [
            jet for jet in jets
            if jet_id(jet) and jet.pt > args.jetPtMin and abs(jet.eta) < args.forwardEtaMax
        ]
        forward_jets = [
            jet for jet in good_jets
            if args.forwardEtaMin <= abs(jet.eta) <= args.forwardEtaMax
        ]
        central_jets = [
            jet for jet in good_jets
            if abs(jet.eta) < args.centralEtaMax
        ]

        central_jets = sorted(central_jets, key=lambda x: x.pt, reverse=True)

        self.n_recojets_total_raw[0] = len(good_jets_raw)
        self.n_recojets_forward_raw[0] = len(forward_jets_raw)
        self.n_recojets_central_raw[0] = len(central_jets_raw)
        self.n_recojets_total[0] = len(good_jets)
        self.n_recojets_forward[0] = len(forward_jets)
        self.n_recojets_central[0] = len(central_jets)

        forward_pair = best_forward_pair(forward_jets)
        if forward_pair is not None:
            self.n_forward_pair += 1
            self.has_forward_pair[0] = 1
            jf1, jf2 = forward_pair
            self.forward_jet1_pt[0] = jf1.pt
            self.forward_jet1_pt_raw[0] = getattr(jf1.original, "pt", -999.0)
            self.forward_jet1_eta[0] = jf1.eta
            self.forward_jet1_phi[0] = jf1.phi
            self.forward_jet1_mass[0] = recojet_mass(jf1)
            self.forward_jet2_pt[0] = jf2.pt
            self.forward_jet2_pt_raw[0] = getattr(jf2.original, "pt", -999.0)
            self.forward_jet2_eta[0] = jf2.eta
            self.forward_jet2_phi[0] = jf2.phi
            self.forward_jet2_mass[0] = recojet_mass(jf2)
            self.forward_pair_dEta[0] = abs(jf1.eta - jf2.eta)
            self.forward_pair_dPhi[0] = abs(delta_phi(jf1.phi, jf2.phi))
            self.forward_pair_dR[0] = delta_r(jf1.eta, jf1.phi, jf2.eta, jf2.phi)
            self.forward_pair_mass[0] = (p4_from_recojet(jf1) + p4_from_recojet(jf2)).M()

            if self.forward_pair_dEta[0] > args.forwardPairDEtaMin:
                self.has_vbf_tag[0] = 1
                self.n_vbf_tag += 1

        central_pair = choose_central_pair(
            central_jets,
            args.centralPairAlgo,
            max_deta=(args.centralPairDEtaMax if args.centralPairSelection == "angular" else None),
            min_dphi=(args.centralPairDPhiMin if args.centralPairSelection == "angular" else None),
            target_mass=args.mass,
        )
        if central_pair is not None:
            jc1, jc2 = central_pair
            passes_central_pair, dR, dEta, dPhi, apt = central_pair_passes_selection(jc1, jc2)

            self.central_jet1_pt[0] = jc1.pt
            self.central_jet1_pt_raw[0] = getattr(jc1.original, "pt", -999.0)
            self.central_jet1_eta[0] = jc1.eta
            self.central_jet1_phi[0] = jc1.phi
            self.central_jet1_mass[0] = recojet_mass(jc1)
            self.central_jet2_pt[0] = jc2.pt
            self.central_jet2_pt_raw[0] = getattr(jc2.original, "pt", -999.0)
            self.central_jet2_eta[0] = jc2.eta
            self.central_jet2_phi[0] = jc2.phi
            self.central_jet2_mass[0] = recojet_mass(jc2)
            self.central_pair_dEta[0] = dEta
            self.central_pair_dPhi[0] = dPhi
            self.central_pair_dR[0] = dR
            self.central_pair_mass[0] = (p4_from_recojet(jc1) + p4_from_recojet(jc2)).M()

            if passes_central_pair:
                self.has_central_pair[0] = 1
                self.n_central_pair += 1

                has_match_info, dr1, dr2 = best_match_distances(jc1, jc2, daughters)
                self.central_pair_match_dR1[0] = dr1
                self.central_pair_match_dR2[0] = dr2
                if has_match_info and dr1 < args.matchDR and dr2 < args.matchDR:
                    self.central_pair_genmatched[0] = 1
                    self.n_genmatched += 1

        self.tree.Fill()
        return True

    def endJob(self):
        dr_min, apt_max = mass_dependent_pair_wp(float(args.mass))
        n_after_trigger = self.tree.GetEntries("pass_trigger_baseline == 1")
        n_forward_pair_after_trigger = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && has_forward_pair == 1"
        )
        n_vbf_tag_after_trigger = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && has_vbf_tag == 1"
        )
        n_central_pair_after_vbf = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1"
        )
        n_genmatched_after_cascade = self.tree.GetEntries(
            "pass_trigger_baseline == 1 && has_vbf_tag == 1 && has_central_pair == 1 && central_pair_genmatched == 1"
        )

        print("\n================ ScoutNano Reco VBF Study ================")
        print(f"Total events:                {self.total}")
        print(f"After trigger baseline:      {n_after_trigger}")
        print(f"Central pair algorithm:      {args.centralPairAlgo}")
        print(f"Central pair selection:      {args.centralPairSelection}")
        print(f"Central jet |eta| max:       {args.centralEtaMax}")
        if args.centralPairSelection == "angular":
            print(f"Central pair |DeltaEta| max: {args.centralPairDEtaMax}")
            print(f"Central pair |DeltaPhi| min: {args.centralPairDPhiMin}")
        else:
            print(f"Central pair DeltaR min:     {dr_min:.3f}")
            print(f"Central pair A_pt max:       {apt_max:.3f}")
        print(f"Events with forward pair:    {n_forward_pair_after_trigger}")
        print(f"Events with VBF tag pair:    {n_vbf_tag_after_trigger}")
        print(f"Events with central pair:    {n_central_pair_after_vbf}")
        print(f"Events with gen-matched pair:{n_genmatched_after_cascade}")
        print(
            f"Forward pair definition:     >=2 forward jets, choose pair with largest |DeltaEta|"
        )
        print(
            f"VBF tag pair definition:     forward pair with |DeltaEta| > {args.forwardPairDEtaMin}"
        )
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
        modules=[RecoVBFScoutNanoStudy(output_name)],
        noOut=True,
        maxEntries=args.maxEntries,
        firstEntry=args.firstEntry,
        prefetch=args.prefetch,
    )
    processor.run()
