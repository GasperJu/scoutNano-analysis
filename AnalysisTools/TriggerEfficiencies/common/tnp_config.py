"""Shared tag-and-probe configuration for trigger-efficiency producers."""

from __future__ import annotations

import math
from typing import Iterable, List

import ROOT
from ROOT import TLorentzVector

REFERENCE_PATHS = ["PFScouting_SingleMuon"]
SIGNAL_PATHS = ["PFScouting_JetHT"]

L1_SEED_HYPOTHESES = [
    {
        "name": "baseline",
        "require_any": ["L1_HTT280er", "L1_SingleJet180"],
        "exclude_all": [],
    },
    {
        "name": "exclusive",
        "require_any": [
            "L1_HTT280er",
            "L1_SingleJet180",
            "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5",
            "L1_ETT2000",
        ],
        "exclude_all": ["L1_HTT200er", "L1_HTT255er"],
    },
    {
        "name": "HTT280_only",
        "require_any": ["L1_HTT280er"],
        "exclude_all": [],
    },
    {
        "name": "anyJet",
        "require_any": ["L1_SingleJet180"],
        "exclude_all": [],
    },
    {
        "name": "HTDouble",
        "require_any": ["L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"],
        "exclude_all": [],
    },
    {
        "name": "HTSingle",
        "require_any": ["L1_HTT280er", "L1_SingleJet180"],
        "exclude_all": [],
    },
]


def pass_jet_id(jet) -> bool:
    aeta = abs(jet.eta)

    total_e = (
        jet.neutralHadronEnergy
        + jet.HFHadronEnergy
        + jet.photonEnergy
        + jet.HFEMEnergy
        + jet.muonEnergy
        + jet.electronEnergy
        + jet.chargedHadronEnergy
    )
    if total_e <= 0:
        return False

    nhf = (jet.neutralHadronEnergy + jet.HFHadronEnergy) / float(total_e)
    nemf = (jet.photonEnergy + jet.HFEMEnergy) / float(total_e)
    mu_frac = jet.muonEnergy / float(total_e)

    charged_mult = jet.chargedHadronMultiplicity + jet.HFHadronMultiplicity
    neutral_mult = jet.neutralHadronMultiplicity + jet.HFEMMultiplicity
    nconst = (
        jet.chargedHadronMultiplicity
        + jet.neutralHadronMultiplicity
        + jet.muonMultiplicity
        + jet.electronMultiplicity
        + jet.photonMultiplicity
    )

    if aeta < 2.6:
        if nhf >= 0.99:
            return False
        if nemf >= 0.90:
            return False
        if nconst <= 1:
            return False
        if charged_mult <= 0:
            return False
        if mu_frac >= 0.80:
            return False
        return True

    if 2.6 <= aeta < 2.7:
        if nemf >= 0.99:
            return False
        if mu_frac >= 0.80:
            return False
        return True

    if 2.7 <= aeta < 3.0:
        if nemf >= 0.99:
            return False
        if neutral_mult <= 1:
            return False
        return True

    if 3.0 <= aeta < 5.0:
        if nemf >= 0.10:
            return False
        return True

    return False


def pass_muon_id(mu) -> bool:
    if mu.pt <= 30:
        return False
    if abs(mu.eta) >= 0.8:
        return False
    if abs(mu.trk_dxy) >= 0.2:
        return False
    if abs(mu.trk_dz) >= 0.5:
        return False
    if mu.normchi2 >= 3:
        return False
    if mu.nValidRecoMuonHits <= 0:
        return False
    if mu.nRecoMuonMatchedStations <= 1:
        return False
    if mu.nValidPixelHits <= 0:
        return False
    if mu.nTrackerLayersWithMeasurement <= 5:
        return False
    return True


def delta_r(jet, mu) -> float:
    deta = jet.eta - mu.eta
    dphi = math.fabs(math.atan2(math.sin(jet.phi - mu.phi), math.cos(jet.phi - mu.phi)))
    return math.sqrt(deta * deta + dphi * dphi)


def clean_jets(jets: Iterable, muons: Iterable, dr_cut: float = 0.4) -> List[TLorentzVector]:
    cleaned = []
    for jet in jets:
        jvec = TLorentzVector()
        jvec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.m)

        for mu in muons:
            if delta_r(jet, mu) >= dr_cut:
                continue
            muvec = TLorentzVector()
            muvec.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, 0.105)
            jvec -= muvec

        if jvec.Pt() < 1:
            continue
        cleaned.append(jvec)

    return cleaned


def evaluate_l1_hypothesis(event, hypothesis) -> bool:
    required = hypothesis.get("require_any", [])
    excluded = hypothesis.get("exclude_all", [])
    if required and not any(bool(getattr(event, bit, False)) for bit in required):
        return False
    if excluded and any(bool(getattr(event, bit, False)) for bit in excluded):
        return False
    return True
