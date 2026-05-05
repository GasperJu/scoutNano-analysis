# **********************************************
# Run 3 dijet scouting analysis:               #
# trigger efficiency studies                   #
# **********************************************

#!/usr/bin/env python3

import os, sys, math, copy
import ROOT
from ROOT import TLorentzVector, TMath

ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
import json
import array
import numpy as np
import fastjet as fj
#print(fj.__version__) 

jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.4)

# Importing tools from nanoAOD processing set up to store the ratio histograms in a root file
# -------------------------------------------------------------------------------------------
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def PseudoJ(list_jets):
    pseudojets = []
    for j in list_jets:

        # Convert Collection Objects (pt, eta, phi, mass) to (px, py, pz, E)
        px = j.Pt() * TMath.Cos(j.Phi())
        py = j.Pt() * TMath.Sin(j.Phi())
        pz = j.Pt() * TMath.SinH(j.Eta())
        E = TMath.Sqrt(px**2 + py**2 + pz**2 + j.M()**2) #j.mass**2)
        
        # Input jets (px, py, pz, E)
        pj = fj.PseudoJet(px, py, pz, E)
        pseudojets.append(pj)

    return pseudojets

def Lorentz(pseudojets):    
    # Convert PseudoJets into TLorentz vectors (able to calculate Eta)
    return [TLorentzVector(pj.px(), pj.py(), pj.pz(), pj.E()) for pj in pseudojets]

def deltaR(jet, mu):
    deta = jet.eta - mu.eta
    dphi = math.fabs(math.atan2(math.sin(jet.phi - mu.phi), math.cos(jet.phi - mu.phi)))
    return math.sqrt(deta*deta + dphi*dphi)

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2.0 * math.pi
    while dphi <= -math.pi:
        dphi += 2.0 * math.pi
    return dphi

# ---------------------------------------------------------------------------------------------------
def JetID(jet):
    eta = jet.eta
    aeta = abs(eta)
    
    totalE = (
        jet.neutralHadronEnergy + jet.HFHadronEnergy + 
        jet.photonEnergy + jet.HFEMEnergy +
        jet.muonEnergy + jet.electronEnergy +
        jet.chargedHadronEnergy
    )
    if totalE <= 0:
        return False

    NHF = (jet.neutralHadronEnergy + jet.HFHadronEnergy) / float(totalE)
    NEMF = (jet.photonEnergy + jet.HFEMEnergy) / float(totalE)
    muFrac = jet.muonEnergy / float(totalE)

    chargedMult = jet.chargedHadronMultiplicity + jet.HFHadronMultiplicity
    neutralMult = jet.neutralHadronMultiplicity + jet.HFEMMultiplicity
    nconst = jet.chargedHadronMultiplicity + jet.neutralHadronMultiplicity + jet.muonMultiplicity + jet.electronMultiplicity + jet.photonMultiplicity 

    # -------- |Eta| < 2.6 --------
    if aeta < 2.6:
        if NHF >= 0.99: return False
        if NEMF >= 0.90: return False
        if nconst <= 1: return False
        if chargedMult <= 0: return False
        if muFrac >= 0.80: return False
        return True

    # --- |Eta| = [2.6,2.7] --------
    if (aeta >= 2.6 and aeta < 2.7):
        if NEMF >= 0.99: return False
        if muFrac >= 0.80: return False
        return True
    
    # --- |Eta| = [2.7,3.0] --------
    if (aeta >= 2.7 and aeta < 3.0):

        if NEMF >= 0.99: return False
        if neutralMult <= 1: return False
        return True

    # --- |Eta| = [3.0,5.0] --------
    if (aeta >= 3.0 and aeta < 5.0):

        if NEMF >= 0.10: return False
        return True

    return False
# ---------------------------------------------------------------------------------------------------
def MuonID(mu):
    if mu.pt <= 30: return False
    if abs(mu.eta) >= 0.8: return False

    if abs(mu.trk_dxy) >= 0.2: return False
    if abs(mu.trk_dz) >= 0.5: return False
    #if mu.trackIso >= 0.15: return False
    if mu.normchi2 >= 3: return False

    if mu.nValidRecoMuonHits <= 0: return False
    if mu.nRecoMuonMatchedStations <= 1: return False
    if mu.nValidPixelHits <= 0: return False
    if mu.nTrackerLayersWithMeasurement <= 5: return False

    return True

def get_trigger_bits(event):
    dst = Object(event, "DST")
    pass_dst = bool(getattr(dst, "PFScouting_JetHT", False))

    pass_htt200 = bool(getattr(event, "L1_HTT200er", False))
    pass_htt255 = bool(getattr(event, "L1_HTT255er", False))
    pass_htt280 = bool(getattr(event, "L1_HTT280er", False))
    pass_singlejet180 = bool(getattr(event, "L1_SingleJet180", False))
    pass_doublejet = bool(getattr(event, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", False))
    pass_ett2000 = bool(getattr(event, "L1_ETT2000", False))

    pass_any_unprescaled = pass_htt280 or pass_singlejet180 or pass_doublejet or pass_ett2000
    pass_nonhtt_unprescaled = pass_singlejet180 or pass_doublejet or pass_ett2000
    prescaled_triggers = pass_htt200 or pass_htt255
    pass_prescaled_only = (pass_htt200 or pass_htt255) and not pass_any_unprescaled

    pass_baseline = (
        pass_dst and
        (pass_htt280 or pass_singlejet180)
    )
    pass_baseline_L1 = (
        pass_htt280 or pass_singlejet180
    )

    return {
        "pass_dst": pass_dst,
        "pass_htt200": pass_htt200,
        "pass_htt255": pass_htt255,
        "pass_htt280": pass_htt280,
        "pass_singlejet180": pass_singlejet180,
        "pass_doublejet": pass_doublejet,
        "pass_ett2000": pass_ett2000,
        "pass_any_unprescaled": pass_any_unprescaled,
        "pass_nonhtt_unprescaled": pass_nonhtt_unprescaled,
        "prescaled_triggers": prescaled_triggers,
        "pass_prescaled_only": pass_prescaled_only,
        "pass_baseline": pass_baseline,
        "pass_baseline_L1": pass_baseline_L1,
    }
# ---------------------------------------------------------------------------------------------------

class TrigDijetHTAnalysis(Module):
    def __init__(self):
        self.writeHistFile=True
        self.reference_paths=reference_paths
        print("[INFO] Reference path: ", reference_paths)
        self.signal_paths=signal_paths

    def is_good_lumi(self, run, lumi):
        """Return True if (run,lumi) is contained in golden JSON self.good_ls."""
        if not self.good_ls:
            print(">> LS discarded!")
            return True  # If no JSON file, accept everything
        run = int(run)
        lumi = int(lumi)
        if run not in self.good_ls:
            return False
        for (a, b) in self.good_ls[run]:
            if a <= lumi <= b:
                #print(">> Good run and LS!")
                return True
        return False
        
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

        # Counters
        # --------
        self.n_total_events = 0 
        self.n_VBF = 0
        self.n_Boosted = 0
        self.n_Resolved = 0
        self.n_Rest = 0
        self.n_Unclassified = 0

        # Histos
        # --------
        self.h_passreftrig = ROOT.TH1F("h_passreftrig" , "; passed ref trigger", 2, 0. , 2.)
        self.h_ht_inclusive_all = ROOT.TH1F("h_ht_inclusive_all", "; H_{T} Inclusive [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_ht_inclusive_passed = ROOT.TH1F("h_ht_inclusive_passed", "; H_{T} Inclusive [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_minv_1_all = ROOT.TH1F("h_minv_1_all", "; M_{jj} VBF [GeV]; Efficiency;", 100, 30., 630.)
        self.h_minv_1_passed = ROOT.TH1F("h_minv_1_passed", "; M_{jj} VBF [GeV]; Efficiency;", 100, 30., 630.)
        self.h_ht_1_all = ROOT.TH1F("h_ht_1_all", "; H_{T} VBF [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_ht_1_passed = ROOT.TH1F("h_ht_1_passed", "; H_{T} VBF [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_minv_2_all = ROOT.TH1F("h_minv_2_all", "; M_{AK8} Boosted ISR [GeV]; Efficiency;", 100, 30., 630.)
        self.h_minv_2_passed = ROOT.TH1F("h_minv_2_passed", "; M_{AK8} Boosted ISR [GeV]; Efficiency;", 100, 30., 630.)
        self.h_ht_2_all = ROOT.TH1F("h_ht_2_all", "; H_{T} Boosted ISR [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_ht_2_passed = ROOT.TH1F("h_ht_2_passed", "; H_{T} Boosted ISR [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_ptak8_2_all = ROOT.TH1F("h_ptak8_2_all", "; p_{T} AK8 Boosted ISR [GeV]; Efficiency;", 100, 30., 630.)
        self.h_ptak8_2_passed = ROOT.TH1F("h_ptak8_2_passed", "; p_{T} AK8 Boosted ISR [GeV]; Efficiency;", 100, 30., 630.)
        self.h_minv_3_all = ROOT.TH1F("h_minv_3_all", "; M_{jj} Resolved ISR [GeV]; Efficiency;", 100, 30., 630.)
        self.h_minv_3_passed = ROOT.TH1F("h_minv_3_passed", "; M_{jj} Resolved ISR [GeV]; Efficiency;", 100, 30., 630.)
        self.h_ht_3_all = ROOT.TH1F("h_ht_3_all", "; H_{T} Resolved ISR [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_ht_3_passed = ROOT.TH1F("h_ht_3_passed", "; H_{T} Resolved ISR [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_minv_4_all = ROOT.TH1F("h_minv_4_all", "; M_{jj} pp [GeV]; Efficiency;", 100, 30., 630.)
        self.h_minv_4_passed = ROOT.TH1F("h_minv_4_passed", "; M_{jj} pp [GeV]; Efficiency;", 100, 30., 630.)
        self.h_ht_4_all = ROOT.TH1F("h_ht_4_all", "; H_{T} pp [GeV]; Efficiency;", 100, 0., 1500.)
        self.h_ht_4_passed = ROOT.TH1F("h_ht_4_passed", "; H_{T} pp [GeV]; Efficiency;", 100, 0., 1500.)

        for h in [
                self.h_passreftrig, self.h_ht_inclusive_all, self.h_ht_inclusive_passed,
                self.h_minv_1_all, self.h_minv_1_passed, self.h_ht_1_all, self.h_ht_1_passed,
                self.h_minv_2_all, self.h_minv_2_passed, self.h_ht_2_all, self.h_ht_2_passed, self.h_ptak8_2_all, self.h_ptak8_2_passed,
                self.h_minv_3_all, self.h_minv_3_passed,self.h_ht_3_all, self.h_ht_3_passed,
                self.h_minv_4_all, self.h_minv_4_passed, self.h_ht_4_all, self.h_ht_4_passed,
        ]:
            self.addObject(h)

        self.golden_json_path = "/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/TRIGGER_EFF/2024_UtilsDataQuality/Cert_Collisions2024_378981_386951_Golden.json"  
        if os.path.exists(self.golden_json_path):
            with open(self.golden_json_path, "r") as f:
                gj = json.load(f)
            # convert to dict[int, list[(start,end),...]]
            self.good_ls = {int(run): [(int(a), int(b)) for (a, b) in ranges] for run, ranges in gj.items()}
            print(f"[INFO] Loaded golden JSON ---{self.golden_json_path}--- with {len(self.good_ls)} runs")
        else:
            self.good_ls = {}
            print(f"[WARN] Golden JSON not found at {self.golden_json_path} - no run/lumi filtering will be applied.")

            
    def analyze(self, event):

        isBoosted = False
        isResolved = False
        HT = 0.0
        Mjj_1 = 0.0
        HT_1 = 0.0
        Mjj_2 = 0.0
        HT_2 = 0.0
        PtAk8_2 = 0.0
        Mjj_3 = 0.0
        HT_3 = 0.0
        Mjj_4 = 0.0
        HT_4 = 0.0

        # -------------------------------------------------------
        # --- Event-level run/lumi selection based on Golden Json
        # -------------------------------------------------------
        run = getattr(event, "run", None)
        lumi = getattr(event, "luminosityBlock", None)  # in NanoAOD this is usually 'luminosityBlock'
        if run is None or lumi is None:
            pass
        else:
            if not self.is_good_lumi(run, lumi):
                print("RUN = ", run)
                print("LS = ", lumi)
                return False           

        # -------------------------
        # --- Reference trigger ---
        # -------------------------
        dst = Object(event, "DST")

        refAccept = any(getattr(dst, path) for path in self.reference_paths)
        self.h_passreftrig.Fill(refAccept)
        if not refAccept:
            return False

        self.n_total_events += 1


        # ---------------------------------
        # ------ Including HLT JECs -------  
        # ---------------------------------
        import correctionlib
        import os
        
        # Path to JSON
        # --------------------------------
        jec_json = os.path.join(
            os.getenv("CMSSW_BASE"),
            "src",
            "2024_UtilsDataQuality",
            "jetHLT_jerc.json"
        )
        
        # Load correction set
        # --------------------------------
        cset = correctionlib.CorrectionSet.from_file(jec_json)
        #print("--- Regular corrections (printing keys): ")
        #for k in cset.keys():
        #    if "AK4" in k:
        #        print(k)
        #print("--- Compound corrections (printing keys): ")
        #print(cset.compound.keys())

        # Get the specific JEC
        # --------------------------------
        jec = cset.compound["HLT_Winter24_V1_MC_L1L2L3Res_AK4PFHLT"]

        # Get Rho in the event
        # --------------------------------
        run = getattr(event, "run", None)                                                                                                            
        rho = getattr(event, "ScoutingRho_fixedGridRhoFastjetAll", 0.0)
        rho = event.ScoutingRho_fixedGridRhoFastjetAll
        #rho = getattr(event, "ScoutingRho", 0.0)
        #rho = event.ScoutingRho
        
        # ---------------------------------
        # --- Preselection requirements ---
        # ---------------------------------
        jets = Collection(event, "ScoutingPFJet")
        njets = getattr(event, "n" + jets._prefix)

        corrected_pts = []
        for j in jets:
            variables = {
                "JetPt": j.pt,
                "JetEta": j.eta,
                "JetPhi": j.phi,
                "JetA": j.jetArea,
                "Rho": rho,
                "run": run,
            }
            inputs = [variables[inp.name] for inp in jec.inputs]
            corr_factor = jec.evaluate(*inputs)
            corrected_pt = j.pt * corr_factor
            corrected_pts.append(corrected_pt)    
            #pt_corr = jec.evaluate(*inputs)
            #print("## PT = ", j.pt, " &&& PT_corr = ", j.pt*pt_corr, " ETA = ", j.eta)

        # --- Jet preselection ---
        njetAcc = [
            j for j, pt_corr in zip(jets, corrected_pts) #j for j in jets 
            if pt_corr > 30 #if j.pt * jec.evaluate(*inputs) > 30 #if j.pt > 30
            and abs(j.eta) < 5
            and JetID(j)          
        ]

        # --- Muon veto or cleaning (ScoutingMuonVtx with dR(mu,jet) < 0.4 ---
        # --------------------------------------------------------------------
        muonsVtx = Collection(event, "ScoutingMuonVtx")
        muons = Collection(event, "ScoutingMuonVtx")
        ngood_muonsVtx = [mu for mu in muonsVtx if MuonID(mu)]

        if (len(ngood_muonsVtx) < 1):
            return False

        clean_jets = []

        # ---  Muon CLEANING:
        # -----------------------------
        for j in njetAcc:
            # Build jet 4-vector
            jvec = TLorentzVector()
            jvec.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.m)
            
            # Check if too close to any muon
            close_muons = [mu for mu in muons if deltaR(j, mu) < 0.4]
           
            if close_muons:
               # Subtract overlapping muons
                for mu in close_muons:
                    muvec = TLorentzVector()
                    muvec.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, 0.105)
                    jvec -= muvec
                   
               # Discard if jet is fully removed
                if jvec.Pt() < 1:
                    continue
                
            clean_jets.append(jvec)
            
        njetAcc = clean_jets

        # --- Global HT definition ---
        cutHT = 0.0
        njetHt = [j for j in clean_jets if j.Pt() > 30 and abs(j.Eta()) < 2.5]
        HT = sum(j.Pt() for j in njetHt)

        if len(njetAcc) < 2:
            return False

        # -------------------------------
        # ------- CATEGORISATION  -------
        # -------------------------------
        sort_jets = sorted(njetAcc, key=lambda x: x.Pt(), reverse=True)
        
        jetdef_ak4 = fj.JetDefinition(fj.antikt_algorithm, 0.4)
        jetdef_ak8 = fj.JetDefinition(fj.antikt_algorithm, 0.8)
        
        # -------------------------------
        # [1] VBF CATEGORY
        # -------------------------------
        VBF_jets = []

        forward_pairs = []
        for i in range(len(sort_jets)):
            for j in range(i + 1, len(sort_jets)):
                eta_i, eta_j = sort_jets[i].Eta(), sort_jets[j].Eta()
                dEta = abs(eta_i - eta_j)
                if (
                        dEta > 5.0
                        and sort_jets[i].Pt() > 30
                        and sort_jets[j].Pt() > 30
                        and abs(eta_i) > 3.0
                        and abs(eta_j) > 3.0
                        and eta_i * eta_j < 0
                ):
                    forward_pairs.append((i, j, dEta))
                    
        if forward_pairs:
            i_fwd, j_fwd, _ = max(forward_pairs, key=lambda x: x[2])
            eta_min, eta_max = sorted([sort_jets[i_fwd].Eta(), sort_jets[j_fwd].Eta()])
            central_jets = [
                j for k, j in enumerate(sort_jets)
                if k not in (i_fwd, j_fwd)
                and abs(j.Eta()) < 2.5
                and (eta_min < j.Eta() < eta_max)
            ]
            if len(central_jets) >= 2:
                central_sorted = sorted(central_jets, key=lambda x: x.Pt(), reverse=True)
                if abs(central_sorted[0].Eta() - central_sorted[1].Eta()) < 1.3:
                    HT_1 = HT

                    v0, v1 = Lorentz(PseudoJ(central_sorted[:2]))
                    if (HT_1 > cutHT and (v0 + v1).M() > 100):
                        Mjj_1 = (v0 + v1).M()
                        self.n_VBF += 1
                    
        else:
            # -----------------
            # Tuning thresholds
            # -----------------
            boosted_pt_cut = 170
            boosted_eta_cut = 3.0
            isr_pt_cut_boosted = 200
            isr_eta_cut_boosted = 3.0
            isr_min_dr_boosted = 2.0
            isr_min_dphi_boosted = 2.5
            resJ_pt_cut = 70.0  
            ISR_pt_cut_resolved = 70.0
   
            # ---------------------------------
            # [2] BOOSTED w/ ISR CATEGORY (AK8)
            # ---------------------------------
            cluster_boosted = fj.ClusterSequence(PseudoJ(sort_jets), jetdef_ak8)
            ak8jets = fj.sorted_by_pt(cluster_boosted.inclusive_jets())

            if len(ak8jets) >= 1:
                boosted_lv = Lorentz(ak8jets)[0]
                if boosted_lv.Pt() > boosted_pt_cut and abs(boosted_lv.Eta()) < boosted_eta_cut:

                    # Now build the ISR jet candidates from AK4 jets that are NOT in the AK8 jet
                    isr_candidates = []
                    for j in sort_jets:
                        j_lv = TLorentzVector()
                        j_lv.SetPtEtaPhiM(j.Pt(), j.Eta(), j.Phi(), j.M())

                        dR = boosted_lv.DeltaR(j_lv)
                        dPhi = abs(delta_phi(boosted_lv.Phi(), j_lv.Phi()))

                        # Select ISR jet if it's well separated and above threshold
                        if (
                            j_lv.Pt() > isr_pt_cut_boosted and
                            abs(j_lv.Eta()) < isr_eta_cut_boosted and
                            dR > isr_min_dr_boosted and
                            dPhi > isr_min_dphi_boosted
                        ):
                            isr_candidates.append(j_lv)

                    # Require at least 1 ISR candidate: if any ISR candidate exists, select the highest-pT one
                    if len(isr_candidates) > 0:
                        ISR_jet = max(isr_candidates, key=lambda x: x.Pt())
                        HT_2 = HT

                        PtAk8_2 = boosted_lv.Pt()
                        if (HT_2 > cutHT and boosted_lv.M() > 50):
                            Mjj_2 = boosted_lv.M()
                            isBoosted = True
                            self.n_Boosted += 1

            if not (isBoosted):
                # -----------------------------------
                # [3] RESOLVED w/ ISR CATEGORY (AK4)
                # -----------------------------------
                cluster_resolved = fj.ClusterSequence(PseudoJ(sort_jets), jetdef_ak4)
                j_ak4 = Lorentz(fj.sorted_by_pt(cluster_resolved.inclusive_jets()))

                # Chosing as dijet pair the one with 0-1 jets and highest ISR pT
                # --------------------------------------------------------------
                if len(j_ak4) >= 3:
                    # Take the two leading AK4 jets (highest pT)
                    i, j = 0, 1
                    min_dEta = abs(j_ak4[i].Eta() - j_ak4[j].Eta())

                    if (min_dEta < 1.3 and
                        j_ak4[i].Pt() >= resJ_pt_cut and
                        j_ak4[j].Pt() >= resJ_pt_cut):
            
                        leadj = TLorentzVector(j_ak4[i])
                        subleadj = TLorentzVector(j_ak4[j])
                        res_pair = leadj + subleadj
                        
                        # Pick ISR jet as the highest-pT one not in the pair
                        ISR_jet = max(
                            (jet for k, jet in enumerate(j_ak4) if k not in (i, j)),
                            key=lambda jet: jet.Pt(),
                            default=None
                        )

                        if ISR_jet and ISR_jet.Pt() >= ISR_pt_cut_resolved:
                            HT_3 = HT
                            if (HT_3 > cutHT and res_pair.M() > 170):
                                Mjj_3 = res_pair.M()
                                isResolved = True
                                self.n_Resolved += 1
                
                # Chosing as dijet pair the one closest in dPhi
                # ---------------------------------------------
                #if len(j_ak4) >= 3:
                #    best_pair = None
                #    min_dPhi = 999
                #    for i in range(len(j_ak4)):
                #        for j in range(i + 1, len(j_ak4)):
                #            dPhi_pair = abs(j_ak4[i].Phi() - j_ak4[j].Phi())
                #            if dPhi_pair < min_dPhi:
                #                min_dPhi = dPhi_pair
                #                best_pair = (i, j)

                #    if best_pair is not None:
                #        i, j = best_pair
                #        min_dEta = abs(j_ak4[i].Eta() - j_ak4[j].Eta())

                #         if min_dEta < 1.3 and j_ak4[i].Pt() >= resJ_pt_cut and j_ak4[j].Pt() >= resJ_pt_cut:
                #            leadj = TLorentzVector(j_ak4[i])
                #            subleadj = TLorentzVector(j_ak4[j])
                #            res_pair = leadj + subleadj
                #            # Pick ISR jet as the highest-pT one not in the pair
                #            ISR_jet = max(
                #                (jet for k, jet in enumerate(j_ak4) if k not in best_pair),
                #                key=lambda jet: jet.Pt(),
                #                default=None
                #            )
                            # Pick ISR jet as the farest in phi
                            #ISR_jet = max(
                            #    (jet for k, jet in enumerate(j_ak4) if k not in best_pair),
                            #    key=lambda jet: abs(res_pair.DeltaPhi(jet)),
                            #    default=None
                            #)
                #            if ISR_jet and ISR_jet.Pt() >= ISR_pt_cut_resolved:
                #            #if ISR_jet and ISR_jet.Pt() >= ISR_pt_cut_resolved and res_pair.M() > 150:
                #                HT_3 = HT
                #                if (HT_3 > cutHT and res_pair.M() > 170):
                #                    Mjj_3 = res_pair.M()
                #                    isResolved = True
                #                    self.n_Resolved += 1

                if not(isResolved):
                    # ===============================
                    # [4] REST CATEGORY
                    # ===============================
                    if len(j_ak4) >= 2 and abs(j_ak4[0].Eta() - j_ak4[1].Eta()) < 1.3:
                        HT_4 = HT

                        if (HT_4 > cutHT and (j_ak4[0] + j_ak4[1]).M() > 100):
                            Mjj_4 = (j_ak4[0] + j_ak4[1]).M()
                            self.n_Rest += 1


        # *************
        # Set NUMERATOR
        # *************
        self.h_ht_inclusive_all.Fill(HT) ## Inclusive
        self.h_minv_1_all.Fill(Mjj_1) ## VBF
        self.h_ht_1_all.Fill(HT_1) ## VBF
        self.h_minv_2_all.Fill(Mjj_2) ## BOOSTED
        self.h_ht_2_all.Fill(HT_2) ## BOOSTED
        self.h_ptak8_2_all.Fill(PtAk8_2) ## BOOSTED
        self.h_minv_3_all.Fill(Mjj_3) ## RESOLVED
        self.h_ht_3_all.Fill(HT_3) ## RESOLVED
        self.h_minv_4_all.Fill(Mjj_4) ## REST
        self.h_ht_4_all.Fill(HT_4) ## REST

        # -------------------------
        # --- L1 bits selection ---
        # -------------------------

        # --- Unprescaled L1 bits in JetHT scouting trigger
        #unprescaled_l1Triggers = [
        #    "L1_HTT280er",
        #    "L1_SingleJet180",
        #    "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5",
        #    "L1_ETT2000",
        #]
        # --- Prescaled L1 bits in JetHT scouting trigger
        #prescaled_l1Triggers = [
        #    "L1_HTT200er",
        #    "L1_HTT255er",
        #]

        trig = get_trigger_bits(event)
        if trig["pass_baseline"]:
            self.h_ht_inclusive_passed.Fill(HT) ## Inclusive
            self.h_minv_1_passed.Fill(Mjj_1) ## VBF
            self.h_ht_1_passed.Fill(HT_1) ## VBF
            self.h_minv_2_passed.Fill(Mjj_2) ## BOOSTED
            self.h_ht_2_passed.Fill(HT_2) ## BOOSTED
            self.h_ptak8_2_passed.Fill(PtAk8_2) ## BOOSTED
            self.h_minv_3_passed.Fill(Mjj_3) ## RESOLVED
            self.h_ht_3_passed.Fill(HT_3) ## RESOLVED
            self.h_minv_4_passed.Fill(Mjj_4) ## REST
            self.h_ht_4_passed.Fill(HT_4) ## REST

        return True
    

    def endJob(self):
        print(f"------------------------------------------------------------")
        print("Summary after processing:")
        print(f"------------------------------------------------------------")
        print(f"Total events after reference trigger: {self.n_total_events}")
        print(f"VBF:      {self.n_VBF}")
        print(f"Boosted:  {self.n_Boosted}")
        print(f"Resolved: {self.n_Resolved}")
        print(f"Rest:     {self.n_Rest}")
        print(f"Sum: {self.n_VBF + self.n_Boosted + self.n_Resolved + self.n_Rest}")
        print(f"------------------------------------------------------------")
        Module.endJob(self)

#reference_paths = ["PFScouting_SingleMuon", "PFScouting_SinglePhotonEB", "PFScouting_ZeroBias"]
reference_paths = ["PFScouting_SingleMuon"]
signal_paths    = ["PFScouting_JetHT"]

preselection= "" #DST_PFScouting_SingleMuon == 1 && DST_PFScouting_JetHT == 1"

### --------------- ###
### Parse arguments ###
### --------------- ###
args = dict(arg.split('=') for arg in sys.argv[1:] if '=' in arg)
#inputFile = args.get('inputFile', 'root://eoscms.cern.ch//eos/cms/store/data/Run2024I/ScoutingPFRun3/NANOAOD/PromptReco-v2/000/386/864/00000/c0920fe5-182a-47c4-bf4c-aa8b4b4c7518.root')
#inputFile = args.get('inputFile', 'root://cms-xrd-global.cern.ch/')
#inputFile = args.get('inputFile', 'root://eoscms.cern.ch//eos/cms/store/data/Run2024I/ScoutingPFRun3/NANOAOD/PromptReco-v2/000/386/951/00000/8547f361-c4a3-468e-97ab-3b58535911ca.root')
inputFile = args.get('inputFile', 'root://cms-xrd-global.cern.ch//store/data/Run2024H/ScoutingPFRun3/NANOAOD/ScoutNano-v1/2810000/6c46f11a-aba9-49c9-ab9d-df054244b544.root')
outputFile = args.get('outputFile', 'histos_DijetHTTrigNanoAOD.root')

### ------- ###
### Running ###
### ------- ###
p = PostProcessor(
    ".",
    [inputFile],
    cut=preselection,
    branchsel=None,
    modules=[TrigDijetHTAnalysis()],
    noOut=True,
    histFileName=outputFile,
    histDirName="DijethtTrigAnalyzerNanoAOD"
)

p.run()
