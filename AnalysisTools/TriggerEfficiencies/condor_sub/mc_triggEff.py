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
        px = j.pt * TMath.Cos(j.phi)
        py = j.pt * TMath.Sin(j.phi)
        pz = j.pt * TMath.SinH(j.eta)
        E = TMath.Sqrt(px**2 + py**2 + pz**2 + j.m**2) #j.mass**2)
        
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
    #if mu.nValidRecoMuonHits <= 0: return False
    #if mu.nRecoMuonMatchedStations <= 3: return False
    #if mu.nValidPixelHits <= 1: return False
    #if mu.nTrackerLayersWithMeasurement <= 7: return False

    return True


# -------------------- TRIGGER SELECTION ------------------------#
def get_trigger_bits(event):                                                                                                                                                                           
    """
    Check for:
      - Hadronic L1 seed selection
    """

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

    pass_baseline_L1 = (
        pass_htt280 or pass_singlejet180
    )
    pass_baseline = pass_dst and pass_baseline_L1

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
        "pass_any_unprescaled": pass_any_unprescaled,
        "pass_baseline": pass_baseline,
        "pass_baseline_L1": pass_baseline_L1
    }

def build_pp_dijet_mass(clean_jets, max_deta=1.3):
    if len(clean_jets) < 2:
        return None
    sort_jets = sorted(clean_jets, key=lambda x: x.Pt(), reverse=True)
    lead, sublead = sort_jets[0], sort_jets[1]
    if abs(lead.Eta() - sublead.Eta()) >= max_deta:
        return None
    return (lead + sublead).M()

class TrigDijetHTAnalysis(Module):
    def __init__(self):
        self.isData=False
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

    def get_event_weight(self, event):
        if self.isData:
            return 1.0
        return event.genWeight if hasattr(event, "genWeight") else 1.0

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # inputFile is a TFile
        runs = inputFile.Get("Runs")
        if not runs:
            return

        sumw  = 0.0
        sumw2 = 0.0
        cnt   = 0.0
        
        for r in runs:
            if hasattr(r, "genEventSumw"):
                sumw += r.genEventSumw
            if hasattr(r, "genEventSumw2"):
                sumw2 += r.genEventSumw2
            if hasattr(r, "genEventCount"):
                cnt += r.genEventCount

        self.h_sumw.Fill(0.5, sumw)
        self.h_sumw2.Fill(0.5, sumw2)
        self.h_cnt.Fill(0.5, cnt)

    
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

        # Metadata
        # --------
        self.h_sumw  = ROOT.TH1D("h_sumw",  "sumw",  1, 0, 1)
        self.h_sumw2 = ROOT.TH1D("h_sumw2", "sumw2", 1, 0, 1)
        self.h_cnt   = ROOT.TH1D("h_cnt",   "count", 1, 0, 1)
        
        self.addObject(self.h_sumw)
        self.addObject(self.h_sumw2)
        self.addObject(self.h_cnt)
        
        # Counters
        # --------
        self.n_totEvents_refTrig = 0 
        self.n_totEvents_refTrigJetId = 0 

        # Histos
        # --------
        self.h_passreftrig = ROOT.TH1F("h_passreftrig" , "; passed ref trigger", 2, 0. , 2.)
        self.h_ht_inclusive_all = ROOT.TH1F("h_ht_inclusive_all", "; H_{T} Inclusive [GeV]; Efficiency;", 150, 0., 1500.)
        self.h_ht_inclusive_passed = ROOT.TH1F("h_ht_inclusive_passed", "; H_{T} Inclusive [GeV]; Efficiency;", 150, 0., 1500.)
        self.h_pt_leading_all = ROOT.TH1F("h_pt_leading_all", "; AK4 Leading jet p_{T} [GeV]; Efficiency;", 50, 0., 500.)
        self.h_pt_leading_passed = ROOT.TH1F("h_pt_leading_passed", "; AK4 Leading jet p_{T} [GeV]; Efficiency;", 50, 0., 500.)
        self.h_dijetMass_all = ROOT.TH1F("h_dijetMass_all", "; Dijet mass [GeV]; Efficiency;", 200, 0., 2000.)
        self.h_dijetMass_passed = ROOT.TH1F("h_dijetMass_passed", "; Dijet mass [GeV]; Efficiency;", 200, 0., 2000.)

        for h in [
                self.h_passreftrig,
                self.h_ht_inclusive_all, self.h_ht_inclusive_passed,
                self.h_pt_leading_all, self.h_pt_leading_passed,
                self.h_dijetMass_all, self.h_dijetMass_passed,
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

        eventWeight = self.get_event_weight(event)

        # -------------------------------------------------------
        # --- Event-level run/lumi selection based on Golden Json
        # -------------------------------------------------------
        run = getattr(event, "run", None)
        lumi = getattr(event, "luminosityBlock", None)  # in NanoAOD this is usually 'luminosityBlock'
        #if run is None or lumi is None:
        #    pass
        #else:
        #    if not self.is_good_lumi(run, lumi):
        #        return False           
        
        # -------------------------
        # --- Reference trigger ---
        # -------------------------
        HT = 0.0        
        dst = Object(event, "DST")

        refAccept = any(getattr(dst, path) for path in self.reference_paths)
        self.h_passreftrig.Fill(refAccept, eventWeight)

        if not refAccept:
            return False

        self.n_totEvents_refTrig += 1

        #print("*** ----------- ***")
        #print("*** Event = ", event)
        #print("*** ----------- ***")
        #print("*** n_total_events_refTrigger: ", self.n_totEvents_refTrig)
        #print("*** ----------- ***")

        jets = Collection(event, "ScoutingPFJet")
        njets = getattr(event, "n" + jets._prefix)

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

        # --- Muon veto (ScoutingMuonNoVtx with dR(mu,jet) < 0.4 ---
        muonsVtx = Collection(event, "ScoutingMuonVtx")
        muons = Collection(event, "ScoutingMuonVtx")
        #print("Number of muons = ", len(muons))
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

         
        # ----------------------------
        # --- Global HT definition ---
        # ----------------------------
        njetHt = [j for j in clean_jets if j.Pt() > 30 and abs(j.Eta()) < 2.5]
        HT = sum(j.Pt() for j in njetHt)


        # **********************************
        # Check trigger as a function of HT:
        # **********************************
        self.h_ht_inclusive_all.Fill(HT, eventWeight) ## Inclusive

        self.n_totEvents_refTrigJetId += 1 

        # **************************************
        # Check trigger as a function of lead pT
        # **************************************
        if len(njetAcc) >= 1:
            sort_jets = sorted(njetAcc, key=lambda x: x.Pt(), reverse=True)
            self.h_pt_leading_all.Fill(sort_jets[0].Pt(), eventWeight) 

        dijet_mass = build_pp_dijet_mass(njetAcc, max_deta=1.3)
        if dijet_mass is not None:
            self.h_dijetMass_all.Fill(dijet_mass, eventWeight)

        #if len(njetAcc) >= 2:
        #    sort_jets = sorted(njetAcc, key=lambda x: x.Pt(), reverse=True)
        #    self.h_pt_leading_all.Fill(sort_jets[0].Pt(), eventWeight) 

        # --- Trigger requirement ---                         
        trig = get_trigger_bits(event)                                                                                    

        if trig["pass_baseline"]:
            self.h_ht_inclusive_passed.Fill(HT, eventWeight) ## Inclusive
            if len(njetAcc) >= 1:
                self.h_pt_leading_passed.Fill(sort_jets[0].Pt(), eventWeight)
            if dijet_mass is not None:
                self.h_dijetMass_passed.Fill(dijet_mass, eventWeight)

        return True
    

    def endJob(self):
        print(f"----------------------------------------------------------------")
        print("Summary after processing:")
        print(f"----------------------------------------------------------------")
        print(f"Total events after reference trigger: {self.n_totEvents_refTrig}")
        print(f"----------------------------------------------------------------")
        Module.endJob(self)

reference_paths = ["PFScouting_SingleMuon"]
signal_paths    = ["PFScouting_JetHT"]

preselection= "" #DST_PFScouting_SingleMuon == 1 && DST_PFScouting_JetHT == 1"

### --------------- ###
### Parse arguments ###
### --------------- ###
args = dict(arg.split('=') for arg in sys.argv[1:] if '=' in arg)
inputFile = args.get('inputFile', 'file:/eos/cms/store/cmst3/group/vhcc/ScoutingNanoAOD/2024/mc/QCD-4Jets_Bin-HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/ScoutNanoTuples-16May2025_Run3_RunIII2024Summer24MiniAOD-140X_v26-v2/250517_120052/0000/nano_918.root')
#inputFile = args.get('inputFile', 'file:/eos/cms/store/cmst3/group/vhcc/ScoutingNanoAOD/2024/mc/QCD-4Jets_Bin-HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/ScoutNanoTuples-16May2025_Run3_RunIII2024Summer24MiniAOD-140X_v26-v2/250517_120110/0000/nano_109.root')
outputFile = args.get('outputFile', 'histosMC_InclusiveTrigNanoAOD.root')

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
    histDirName="InclusiveTrigNanoAOD"
)

p.run()
