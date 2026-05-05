#############################################################
# Run 3 scouting analysis: Boosted (+ ISR) category - Utils #
#############################################################

#!/usr/bin/env python3

import math
import os
import glob
import json
import ROOT
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object

# ------------------ Event-related functions ------------------- #
def load_golden_json(path):
    if not path or not os.path.exists(path):
        print(f"[WARN] Golden JSON not found: {path}")
        return {}

    with open(path, "r") as f:
        gj = json.load(f)

    good_ls = {
        int(run): [(int(a), int(b)) for (a, b) in ranges]
        for run, ranges in gj.items()
    }
    print(f"[INFO] Loaded golden JSON {path} with {len(good_ls)} runs")
    return good_ls


def is_good_lumi(good_ls, run, lumi):
    if not good_ls:
        return True
    run = int(run)
    lumi = int(lumi)
    if run not in good_ls:
        return False
    for a, b in good_ls[run]:
        if a <= lumi <= b:
            return True
    return False

def get_event_weight(event, is_data):
    if is_data:
        return 1.0
    return float(getattr(event, "genWeight", 1.0))

# -------------------- TRIGGER SELECTION ------------------------#

def get_trigger_bits(event):                                                                                                                                                                           
    """
    Check for:
      - DST.PFScouting_JetHT
      - Hadronic L1 seed selection
      - Veto specific muon-seeded L1 bits
    """

    dst = Object(event, "DST")
    pass_dst = bool(getattr(dst, "PFScouting_JetHT", False))

    pass_mu_veto = not (
        getattr(event, "L1_SingleMu11_SQ14_BMTF", False) or
        getattr(event, "L1_SingleMu10_SQ14_BMTF", False) or
        getattr(event, "L1_SingleMu9_SQ14_BMTF", False) or
        getattr(event, "L1_SingleMu8_SQ14_BMTF", False) or
        getattr(event, "L1_SingleMu7_SQ14_BMTF", False) or
        getattr(event, "L1_SingleMu6_SQ14_BMTF", False) or
        getattr(event, "L1_SingleMu5_SQ14_BMTF", False)
    )

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
        pass_mu_veto and
        (pass_htt280 or pass_singlejet180)
        #pass_any_unprescaled and 
        #not pass_prescaled_only
    )

    return {
        "pass_dst": pass_dst,
        "pass_mu_veto": pass_mu_veto,
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
        "pass_baseline": pass_baseline
    }

# --------------------- Jet ID UTILITIES ------------------------ #
def jetid_ak4(i, eta_arr,
              nh_arr, hfh_arr, pho_arr, hfem_arr,
              muE_arr, eleE_arr, chE_arr,
              chMult_arr, hfhMult_arr, neuMult_arr, hfemMult_arr,
              muMult_arr, eleMult_arr, phoMult_arr):
    eta = float(eta_arr[i])
    aeta = abs(eta)

    nh  = float(nh_arr[i])
    hfh = float(hfh_arr[i])
    pho = float(pho_arr[i])
    hfem = float(hfem_arr[i])
    muE = float(muE_arr[i])
    eleE = float(eleE_arr[i])
    chE = float(chE_arr[i])

    totalE = nh + hfh + pho + hfem + muE + eleE + chE
    if totalE <= 0:
        return False

    NHF = (nh + hfh) / totalE
    NEMF = (pho + hfem) / totalE
    muFrac = muE / totalE

    chargedMult = int(chMult_arr[i]) + int(hfhMult_arr[i])
    neutralMult = int(neuMult_arr[i]) + int(hfemMult_arr[i])

    nconst = (
        int(chMult_arr[i]) +
        int(neuMult_arr[i]) +
        int(muMult_arr[i]) +
        int(eleMult_arr[i]) +
        int(phoMult_arr[i])
    )

    if aeta < 2.6:
        if NHF >= 0.99: return False
        if NEMF >= 0.90: return False
        if nconst <= 1: return False
        if chargedMult <= 0: return False
        if muFrac >= 0.80: return False
        return True

    if 2.6 <= aeta < 2.7:
        if NEMF >= 0.99: return False
        if muFrac >= 0.80: return False
        return True

    if 2.7 <= aeta < 3.0:
        if NEMF >= 0.99: return False
        if neutralMult <= 1: return False
        return True

    if 3.0 <= aeta < 5.0:
        if NEMF >= 0.10: return False
        return True

    return False

def jetid_ak4_recluster(i, eta_arr, nhf_arr, nemf_arr, muf_arr,
                    nch_arr, nnh_arr, nconst_arr):
    eta = float(eta_arr[i])
    aeta = abs(eta)

    NHF = float(nhf_arr[i])
    NEMF = float(nemf_arr[i])
    muFrac = float(muf_arr[i])

    chargedMult = int(nch_arr[i])
    neutralMult = int(nnh_arr[i])
    nconst = int(nconst_arr[i])

    # -------- |eta| < 2.6 --------
    if aeta < 2.6:
        if NHF >= 0.99: return False
        if NEMF >= 0.90: return False
        if nconst <= 1: return False
        if chargedMult <= 0: return False
        if muFrac >= 0.80: return False
        return True

    # -------- 2.6 <= |eta| < 2.7 --------
    if 2.6 <= aeta < 2.7:
        if NEMF >= 0.99: return False
        if muFrac >= 0.80: return False
        return True

    # -------- 2.7 <= |eta| < 3.0 --------
    if 2.7 <= aeta < 3.0:
        if NEMF >= 0.99: return False
        if neutralMult <= 0: return False
        return True

    return False

# --------------------- HLT JECs UTILITIES ------------------------ #
def load_ak4_jec(jec_json_path, jec_name):
    if not jec_json_path or not os.path.exists(jec_json_path):
        print(f"[WARN] JEC JSON not found: {jec_json_path}")
        return None, None

    cset = correctionlib.CorrectionSet.from_file(jec_json_path)
    jec = cset.compound[jec_name]
    input_names = [inp.name for inp in jec.inputs]
    print(f"[INFO] Loaded JEC: {jec_name} from {jec_json_path}")
    return jec, input_names

def evaluate_ak4_jec(jec, input_names, pt, eta, phi, area, rho, run):
    if jec is None:
        return 1.0

    inputs = []
    for name in input_names:
        if name == "JetPt":
            inputs.append(pt)
        elif name == "JetEta":
            inputs.append(eta)
        elif name == "JetPhi":
            inputs.append(phi)
        elif name == "JetA":
            inputs.append(area)
        elif name == "Rho":
            inputs.append(rho)
        elif name == "run":
            inputs.append(run)
        else:
            raise RuntimeError(f"Unsupported JEC input: {name}")

    return float(jec.evaluate(*inputs))
# ------------------------------------------------------------------ #

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2.0 * math.pi
    while dphi <= -math.pi:
        dphi += 2.0 * math.pi
    return dphi

def delta_r(eta1, phi1, eta2, phi2):
    dphi = delta_phi(phi1, phi2)
    deta = eta1 - eta2
    return math.sqrt(deta * deta + dphi * dphi)

def score_sum(pxqq, pxbb, pxcc, pxgg, pqcd):
#def score_sum(pxqq, pxbb, pxcc,  pqcd):
    xsum = pxqq + pxbb + pxcc + pxgg
    denom = xsum + pqcd
    if denom <= 0:
        return -1.0
    return xsum / denom

def get_fatjet_prefix(rjet):
    if rjet == "15":
        return "ScoutingFatPFJet15Recluster"
    if rjet == "8":
        return "ScoutingFatPFJetRecluster"
    raise ValueError(f"Unknown rjet={rjet}")

def choose_fatjet(jets, wp=None, mode="bestScore"):
    if len(jets) == 0:
        return None

    if mode == "bestScore":
        passing = [j for j in jets if j["score_xsum"] > wp]
        if len(passing) == 0:
            return None
        return max(passing, key=lambda x: x["score_xsum"])

    if mode == "leading":
        lead = max(jets, key=lambda x: x["pt"])
        if wp is not None and lead["score_xsum"] <= wp:
            return None
        return lead

    if mode == "leadingTagged":
        jets_pt_sorted = sorted(jets, key=lambda x: x["pt"], reverse=True)
        for j in jets_pt_sorted:
            if j["score_xsum"] > wp:
                return j
        return None

    raise ValueError(f"Unknown jet choice mode: {mode}")

def build_good_ak4_jets(event, ak4_max_eta, jec=None, jec_input_names=None, apply_jec_eta_max=3.0):
    prefix = "ScoutingPFJet"
    nAK4 = int(getattr(event, f"n{prefix}", 0))
    if nAK4 == 0:
        return []

    ak4_pt   = getattr(event, f"{prefix}_pt", [])
    ak4_eta  = getattr(event, f"{prefix}_eta", [])
    ak4_phi  = getattr(event, f"{prefix}_phi", [])
    ak4_area = getattr(event, f"{prefix}_jetArea", [])

    nh_arr    = getattr(event, f"{prefix}_neutralHadronEnergy", [])
    hfh_arr   = getattr(event, f"{prefix}_HFHadronEnergy", [])
    pho_arr   = getattr(event, f"{prefix}_photonEnergy", [])
    hfem_arr  = getattr(event, f"{prefix}_HFEMEnergy", [])
    muE_arr   = getattr(event, f"{prefix}_muonEnergy", [])
    eleE_arr  = getattr(event, f"{prefix}_electronEnergy", [])
    chE_arr   = getattr(event, f"{prefix}_chargedHadronEnergy", [])

    chMult_arr   = getattr(event, f"{prefix}_chargedHadronMultiplicity", [])
    hfhMult_arr  = getattr(event, f"{prefix}_HFHadronMultiplicity", [])
    neuMult_arr  = getattr(event, f"{prefix}_neutralHadronMultiplicity", [])
    hfemMult_arr = getattr(event, f"{prefix}_HFEMMultiplicity", [])
    muMult_arr   = getattr(event, f"{prefix}_muonMultiplicity", [])
    eleMult_arr  = getattr(event, f"{prefix}_electronMultiplicity", [])
    phoMult_arr  = getattr(event, f"{prefix}_photonMultiplicity", [])

    run = int(getattr(event, "run", 1))
    rho = float(getattr(event, "ScoutingRho_fixedGridRhoFastjetAll", 0.0))

    jets = []
    for j in range(nAK4):
        eta = float(ak4_eta[j])
        if abs(eta) > ak4_max_eta:
            continue

        if not jetid_ak4(
            j, ak4_eta,
            nh_arr, hfh_arr, pho_arr, hfem_arr,
            muE_arr, eleE_arr, chE_arr,
            chMult_arr, hfhMult_arr, neuMult_arr, hfemMult_arr,
            muMult_arr, eleMult_arr, phoMult_arr
        ):
            continue

        pt_raw = float(ak4_pt[j])
        phi    = float(ak4_phi[j])
        area   = float(ak4_area[j])

        # Apply JEC only in central region, according to ISR selection
        # ------------------------------------------------------------
        if abs(eta) < apply_jec_eta_max and jec is not None:
            corr = evaluate_ak4_jec(jec, jec_input_names, pt_raw, eta, phi, area, rho, run)
            pt = pt_raw * corr
        else:
            corr = 1.0
            pt = pt_raw

        #corr = evaluate_ak4_jec(jec, jec_input_names, pt_raw, eta, phi, area, rho, run)
        #pt = pt_raw * corr

        jets.append({
            "idx": j,
            "pt": pt,
            "pt_raw": pt_raw,
            "eta": eta,
            "phi": phi,
            "corr": corr,
        })

    return jets

def build_good_isr_recluster_jets(event, isr_max_eta):
    prefix = "ScoutingPFJetRecluster"
    njet = int(getattr(event, f"n{prefix}", 0))
    if njet == 0:
        return []

    pt_arr   = getattr(event, f"{prefix}_pt", [])
    eta_arr  = getattr(event, f"{prefix}_eta", [])
    phi_arr  = getattr(event, f"{prefix}_phi", [])

    nhf_arr    = getattr(event, f"{prefix}_neHEF", [])
    nemf_arr   = getattr(event, f"{prefix}_neEmEF", [])
    muf_arr    = getattr(event, f"{prefix}_muEF", [])
    nch_arr    = getattr(event, f"{prefix}_nCh", [])
    nnh_arr    = getattr(event, f"{prefix}_nNh", [])
    nconst_arr = getattr(event, f"{prefix}_nConstituents", [])

    jets = []
    for j in range(njet):
        eta = float(eta_arr[j])
        if abs(eta) > isr_max_eta:
            continue

        if not jetid_ak4_recluster(
            j, eta_arr, nhf_arr, nemf_arr, muf_arr,
            nch_arr, nnh_arr, nconst_arr
        ):
            continue

        jets.append({
            "idx": j,
            "pt": float(pt_arr[j]),
            "eta": eta,
            "phi": float(phi_arr[j]),
        })

    return jets

def get_best_isr_from_list(ak4jets, fat_eta, fat_phi, min_dr, min_dphi, ak4_min_pt):
    candidates = []
    for j in ak4jets:
        if j["pt"] < ak4_min_pt:
            continue
        if delta_r(fat_eta, fat_phi, j["eta"], j["phi"]) < min_dr:
            continue
        if abs(delta_phi(fat_phi, j["phi"])) < min_dphi:
            continue
        candidates.append(j)

    if not candidates:
        return None

    return max(candidates, key=lambda x: x["pt"])

def get_isr_category(pt):
    if 200.0 <= pt < 300.0:
        return "200to300"
    if 300.0 <= pt < 400.0:
        return "300to400"
    if 400.0 <= pt < 500.0:
        return "400to500"
    if pt >= 500.0:
        return "ge500"
    return None

def build_chain(tree_name, pattern):
    files = sorted(glob.glob(pattern))
    if len(files) == 0 and os.path.isfile(pattern):
        files = [pattern]
    if len(files) == 0:
        raise RuntimeError(f"No files matched input: {pattern}")

    chain = ROOT.TChain(tree_name)
    for f in files:
        chain.Add(f)
    return chain, files




# ------ Plots -------------
def cms_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.1, "X")
    ROOT.gStyle.SetTitleOffset(1.3, "Y")

def draw_cms_label(extra_text="Simulation Preliminary", lumi_text="2024 (13.6 TeV)"):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.14, 0.93, "CMS")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextFont(52)
    latex2.SetTextSize(0.04)
    latex2.DrawLatex(0.24, 0.93, extra_text)

    latex3 = ROOT.TLatex()
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.04)
    latex3.SetTextAlign(31)
    latex3.DrawLatex(0.95, 0.93, lumi_text)

def ratio(x, qcd):
    denom = x + qcd
    if denom <= 0:
        return -1.0
    return x / denom

def make_output_root_name(outdir, rjet, plusISR, jetChoice, massVar):
    suffix = "_plusISR" if plusISR else ""
    return os.path.join(outdir, f"fatjet_ZPrime_M50_AK{rjet}_{jetChoice}_{massVar}{suffix}.root")

def draw_overlaid_hists(hmap, outname, xtitle, ytitle="Events", vertical_line=None, header=None, logy=False):
    c = ROOT.TCanvas(f"c_{outname}", "", 900, 800)
    c.cd()

    if logy:
        c.SetLogy()

    colors = [
        ROOT.kBlack,
        ROOT.kBlue + 1,
        ROOT.kRed + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 2,
    ]

    ymax = max(h.GetMaximum() for h in hmap.values()) if len(hmap) else 0.0
    ymax = 50 * ymax if logy else 1.35 * ymax
    ymin = 0.5 if logy else 0.0

    leg = ROOT.TLegend(0.15, 0.480, 0.55, 0.88)
    leg.SetTextSize(0.032)
    if header:
        leg.SetHeader(header)

    first = True
    for i, (label, h) in enumerate(hmap.items()):
        h.SetLineColor(colors[i % len(colors)])
        h.SetLineWidth(3)
        h.SetMaximum(ymax)
        h.SetMinimum(ymin)
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)

        if first:
            h.Draw("hist")
            first = False
        else:
            h.Draw("hist same")

        leg.AddEntry(h, label, "l")

    if vertical_line is not None:
        line_ymin = ymin if logy else 0.0
        line = ROOT.TLine(vertical_line, line_ymin, vertical_line, ymax * 0.92)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw("same")

    leg.Draw()
    draw_cms_label()
    c.SaveAs(outname + ".pdf")
    c.SaveAs(outname + ".png")
