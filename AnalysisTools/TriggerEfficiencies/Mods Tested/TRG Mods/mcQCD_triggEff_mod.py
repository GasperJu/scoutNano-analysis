# **********************************************
# Run 3 dijet scouting analysis:               #
# trigger efficiency studies                   #
# **********************************************

import concurrent.futures
import time

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

def process_data(input_value, output_value):
    """A sample function that processes an input value."""
    print(f"Processing started for input: {input_value}")
    p.run()
    print(f"Processing finished for input: {input_value}, result: {output_value}")
    #return result

jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.4)

# Importing tools from nanoAOD processing set up to store the ratio histograms in a root file
# -------------------------------------------------------------------------------------------
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def read_list_from_file(file_path):
    print(file_path)
    try:
        with open(file_path, 'r') as file:
            # Choose one of these methods:
            # lines = file.readlines()          # Keeps newline characters
            #lines = [line.strip() for line in file]  # Removes newline characters
            lines = file.read().splitlines()

        #print(f"File contents as list: {lines}"
        print(f"Number of lines: {len(lines)}")
        return lines

    except FileNotFoundError:
        print("Error: File not found!")
        return []
    except IOError:
        print("Error: Could not read the file!")
        return []

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
    if mu.nRecoMuonMatchedStations <= 3: return False
    if mu.nValidPixelHits <= 1: return False
    if mu.nTrackerLayersWithMeasurement <= 7: return False

    return True


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
        
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)

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

        for h in [
                self.h_passreftrig,
                self.h_ht_inclusive_all, self.h_ht_inclusive_passed,
                self.h_pt_leading_all, self.h_pt_leading_passed,
        ]:
            self.addObject(h)

        self.golden_json_path = "/eos/user/j/jleite/SecFAILING/CMSSW_14_0_12/src/Boosted-Elisa/TriggerEfficiencies/GoldenJSON/Cert_Collisions2024_378981_386951_Golden.json"
        #"/afs/cern.ch/work/e/elfontan/private/dijetAnalysis_ScoutingRun3/TRIGGER_EFF/2024_UtilsDataQuality/Cert_Collisions2024_378981_386951_Golden.json"  
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
        # --- Preselection requirements ---
        # ---------------------------------

        # --- Jet preselection ---
        njetAcc = [
            j for j in jets
            if j.pt > 30
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

        input_l1Triggers = [
            "L1_HTT200er",
            "L1_HTT255er",
            "L1_HTT280er",
            "L1_SingleJet180",
            "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5",
            "L1_ETT2000",
        ]

        fired = []
        for trig in input_l1Triggers:
            if getattr(event, trig, False):
                fired.append(trig)
                
        # --- DEBUG PRINTOUT to check the L1 bit requirement ---
        #print(f"[DEBUG] Event {event.event}: L1 fired -> {fired}")

        # --- Keep only events in which one of the unprescaled triggers in the JetHT scouting path fired ---          
        if not (                                                                                              
                getattr(event, "L1_HTT280er", False) or                                                        
                getattr(event, "L1_SingleJet180", False) or                                                        
                getattr(event, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", False) or                                                        
                getattr(event, "L1_ETT2000", False)                                                       
        ):                                                                                    
            return False                                                                                                          

        # --- Discard events in which one of the prescaled/disabled triggers in the JetHT scouting path fired ---          
        noHTT_triggers = (
            getattr(event, "L1_SingleJet180", False) or
            getattr(event, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", False) or
            getattr(event, "L1_ETT2000", False)
        )
        
        disabled_triggers = (
            getattr(event, "L1_HTT200er", False) or
            getattr(event, "L1_HTT255er", False)
        )

        # Discard events triggered ONLY by bad triggers
        if disabled_triggers and not noHTT_triggers:
            return False

        if dst.PFScouting_JetHT == 1:
            self.h_ht_inclusive_passed.Fill(HT, eventWeight) ## Inclusive
            if len(njetAcc) >= 1:
                self.h_pt_leading_passed.Fill(sort_jets[0].Pt(), eventWeight) 
                #self.h_pt_leading_passed.Fill(sort_jets[0].pt) 

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

## Sample NanoAOD Scouting Format
Input = "/eos/user/j/jleite/SecFAILING/CMSSW_14_0_12/src/ScoutingDijetAN/"

##### BKG
## QCD-4Jets (HT=100-200 GeV)
QCD_4Jets_Bin_HT_100to200 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=200-400 GeV)
QCD_4Jets_Bin_HT_200to400 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=400-600 GeV)
QCD_4Jets_Bin_HT_400to600 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=600-800 GeV)
QCD_4Jets_Bin_HT_600to800 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=800-1000 GeV)
QCD_4Jets_Bin_HT_800to1000 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=1000-1200 GeV)
QCD_4Jets_Bin_HT_1000to1200 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=1200-1500 GeV)
QCD_4Jets_Bin_HT_1200to1500 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT=1500-2000 GeV)
QCD_4Jets_Bin_HT_1500to2000 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## QCD-4Jets (HT>2000 GeV)
QCD_4Jets_Bin_HT_2000 = read_list_from_file(str(Input)+"bkg/QCD-4Jets_Bin-HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt")
## TTto4Q (2024)
## Full 7740
TTto4Q = read_list_from_file(str(Input)+"bkg/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8.txt")
## Reduced (2000)
TTto4Q_red = read_list_from_file(str(Input)+"bkg/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8-reduced.txt")
## Small (1000)
TTto4Q_s =  read_list_from_file(str(Input)+"bkg/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8-small.txt")
## WWto4Q (2024)
WWto4Q = read_list_from_file(str(Input)+"bkg/WWto4Q_TuneCP5_13p6TeV_powheg-pythia8.txt")

### --------------- ###
### Parse arguments ###
### --------------- ###
#args = dict(arg.split('=') for arg in sys.argv[1:] if '=' in arg)
inputFile = QCD_4Jets_Bin_HT_800to1000 #TTto4Q
#args.get('inputFile', WWto4Q)
#'file:/eos/cms/store/cmst3/group/vhcc/ScoutingNanoAOD/2024/mc/QCD-4Jets_Bin-HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/ScoutNanoTuples-16May2025_Run3_RunIII2024Summer24MiniAOD-140X_v26-v2/250517_120110/0000/nano_109.root')
outputFile = 'histosQCD_HT_800to1000_InclusiveTrigNanoAOD.root'
#args.get('outputFile', 'histosWWto4Q_InclusiveTrigNanoAOD.root')

### ------- ###
### Running ###
### ------- ###
p = PostProcessor(
    ".",
    inputFile,
    cut=preselection,
    branchsel=None,
    modules=[TrigDijetHTAnalysis()],
    noOut=True,
    histFileName=outputFile,
    histDirName="InclusiveTrigNanoAOD"
)

p.run()

'''
# List of different inputs to process
inputs_list = [1, 2, 3, 4, 5]

# use Parse
# Use ProcessPoolExecutor to run the function with different inputs in parallel
with concurrent.futures.ProcessPoolExecutor() as executor:
    # Map the function to the list of inputs
    results = executor.map(process_data, inputs_list, output_list)

print("\nAll tasks completed.")
'''
