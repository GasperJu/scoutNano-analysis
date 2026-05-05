# ******************************************** #
# Run 3 dijet scouting analysis:               #
# trigger efficiency studies                   #
# ******************************************** #

#!/usr/bin/env python3

#================
# Import modules
#================
from argparse import ArgumentParser
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)

#================
# Data info
#================
SAMPLES = {
    ##"histos_QCD-HT100to200_JECs.root": 25360000.0,
    #"histos_QCD-HT200to400_JECs.root": 1951000.0,
    #"histos_QCD-HT400to600_JECs.root": 96660.0,
    #"histos_QCD-HT600to800_JECs.root": 13684.0,
    #"histos_QCD-HT800to1000_JECs.root": 3047.0,
    #"histos_TT4Q_JECs.root": 762.1,
    "histosMC_TT4Q_newL1TrgLogic.root": 762.1,
    #"histos_Wto2Q_JECs.root": 16100.0,
}

BASE_PATH = "/eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff"


def getCanvas():
    d = ROOT.TCanvas("", "", 800, 700)
    d.SetRightMargin(0.15)
    d.SetLeftMargin(0.14)
    d.SetBottomMargin(0.15)
    return d

def AddText(setx=0.242, sety=0.91):
    tex = ROOT.TLatex(0.,0.,'Simulation Preliminary');
    tex.SetNDC();
    tex.SetX(setx);
    tex.SetY(sety);
    tex.SetTextFont(53);
    tex.SetTextSize(29);
    tex.SetLineWidth(2)
    return tex

def AddCMSText(setx=0.227, sety=0.91):
    texcms = ROOT.TLatex(0.,0., 'CMS');
    texcms.SetNDC();
    texcms.SetTextAlign(31);
    texcms.SetX(setx);
    texcms.SetY(sety);
    texcms.SetTextFont(63);
    texcms.SetLineWidth(2);
    texcms.SetTextSize(35);
    return texcms

def createLegend():
    legend = ROOT.TLegend(0.27, 0.15, 0.8, 0.3)
    legend.SetFillColor(0)
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    return legend

def SetStyle(h, color, marker=21):
    h.SetMarkerStyle(marker)
    h.SetMarkerColor(color)
    h.SetLineColor(color)
    return h

def get_sumw(fdir):
    h = fdir.Get("h_sumw")
    if not h:
        raise RuntimeError("h_sumw missing")
    return h.Integral()

def get_sumw2(fdir):
    h = fdir.Get("h_sumw2")
    if not h:
        raise RuntimeError("h_sumw2 missing")
    return h.Integral()

def clone_detach(h, name):
    out = h.Clone(name)
    out.SetDirectory(0)
    out.Sumw2()
    return out

def check_same_binning(href, hnew, label=""):
    if href.GetNbinsX() != hnew.GetNbinsX():
        raise RuntimeError(f"[{label}] nbins mismatch: {href.GetNbinsX()} vs {hnew.GetNbinsX()}")
    ax1, ax2 = href.GetXaxis(), hnew.GetXaxis()
    if ax1.GetXmin() != ax2.GetXmin() or ax1.GetXmax() != ax2.GetXmax():
        raise RuntimeError(f"[{label}] x-range mismatch: [{ax1.GetXmin()},{ax1.GetXmax()}] vs [{ax2.GetXmin()},{ax2.GetXmax()}]")
    # also catch variable binning differences
    for i in range(1, href.GetNbinsX() + 2):
        if ax1.GetBinLowEdge(i) != ax2.GetBinLowEdge(i):
            raise RuntimeError(f"[{label}] bin edge mismatch at edge {i}")


colors = {0: ROOT.kMagenta+2, #ROOT.kBlack,
          1: ROOT.kBlue,
          2: ROOT.kGreen+1,
          3: ROOT.kRed+1,
          4: ROOT.kOrange-3,
          5: ROOT.kAzure-3,
          6: ROOT.kTeal+3,
          }
          
def main(args):

    LUMI_PB = 1000.0 
    statOption = ROOT.TEfficiency.kFCP
    variables  = ["ht_inclusive", "pt_leading", "dijetMass"]
    triggers   = ["passed"] #HLT or DST (only the las tone now)

    for var in variables:

        # ----- Build canvas and legend -----
        c = getCanvas()
        leg = createLegend()

        # ----- Build numerator and denominator -----
        den_tot = None
        num_tot = {trg: None for trg in triggers}
    
        print("------------------------------------")
        for fname, xs in SAMPLES.items():
            f = ROOT.TFile.Open(f"{BASE_PATH}/{fname}")
            fdir = f.GetDirectory("InclusiveTrigNanoAOD")
            sumw = get_sumw(fdir)
            sumw2 = get_sumw2(fdir)
            scale = xs * LUMI_PB / sumw
            if sumw <= 0:
                raise RuntimeError(f"Invalid sumw={sumw} in {fname}")

            print(f"[{fname}] xs={xs} pb, sumw={sumw:.6g}, scale={scale:.6g}")

            h_den_in = fdir.Get(f"h_{var}_all")
            if not h_den_in:
                raise RuntimeError(f"Missing h_{var}_all in {fname}")

            den = clone_detach(h_den_in, f"den_{fname}")
            print(f"    den raw integral = {h_den_in.Integral():.6g}")
            den.Scale(scale)
            print(f"    den scaled integral = {den.Integral():.6g}")
            
            if den_tot is None:
                den_tot = den.Clone("den_tot")
                den_tot.SetDirectory(0)
            else:
                check_same_binning(den_tot, den, label=f"den {var}")
                den_tot.Add(den)
                
            for trg in triggers:
                h_num_in = fdir.Get(f"h_{var}_{trg}")
                if not h_num_in:
                    raise RuntimeError(f"Missing h_{var}_{trg} in {fname}")            

                num = clone_detach(h_num_in, f"num_{trg}_{fname}")
                print(f"    num({trg}) raw integral = {h_num_in.Integral():.6g}")
                num.Scale(scale)
                print(f"    num({trg}) scaled integral = {num.Integral():.6g}")
                
                if num_tot[trg] is None:
                    num_tot[trg] = num.Clone(f"num_tot_{trg}")
                    num_tot[trg].SetDirectory(0)
                else:
                    check_same_binning(num_tot[trg], num, label=f"num {var} {trg}")
                    num_tot[trg].Add(num)

            f.Close()
        print("------------------------------------")
            
        # Build TEfficiency from the TOTALS
        # ---------------------------------
        effs = {}
        for j, trg in enumerate(triggers):
            if not ROOT.TEfficiency.CheckConsistency(num_tot[trg], den_tot):
                # This typically means some bin has num > den (can happen due to rounding/overflow/underflow handling)
                raise RuntimeError(f"TEfficiency consistency check failed for {trg} (num>den in some bin?)")
            
            effs[trg] = ROOT.TEfficiency(num_tot[trg], den_tot)
            effs[trg].SetStatisticOption(ROOT.TEfficiency.kFCP)
            effs[trg] = SetStyle(effs[trg], colors[j])
                        
            if j == 0:
                effs[trg].Draw()
            else:
                effs[trg].Draw("same")                
            leg.AddEntry(effs[trg], trg.replace("passed", "DST_PFScouting_JetHT"), "lep") # "p" option to remove filled box

        c.Modified()
        c.Update()

        # ----- Styling stuff: axes ------
        graph = effs[triggers[0]].GetPaintedGraph()
        axis = graph.GetHistogram()

        axis.GetXaxis().SetTitleSize(0.05)
        axis.GetYaxis().SetTitleSize(0.05)
        axis.GetXaxis().SetLabelSize(0.045)
        axis.GetYaxis().SetLabelSize(0.045)
        axis.GetXaxis().SetTitleOffset(1.2)
        axis.GetYaxis().SetTitleOffset(1.25)
        axis.GetYaxis().SetRangeUser(-0.02,1.1)
        
        # ----- Add vertical line only for HT plots ------
        if var == "ht_inclusive":
            line = ROOT.TLine(280, -0.017, 280, 1.09)      # (x1, y1, x2, y2)
            line.SetLineStyle(2)                    
            line.SetLineWidth(2)
            line.SetLineColor(ROOT.kAzure+6)
            #line.Draw("same")
    
        # ----- Styling stuff: add texts ------    
        leg.Draw("same")

        tex_cms = AddCMSText()
        tex_cms.Draw("same")
        
        addText = AddText()
        addText.Draw("same")
        
        header = ROOT.TLatex()
        header.SetTextSize(0.045)
        header.DrawLatexNDC(0.6, 0.91, "2024 (13.6 TeV)")

        c.Update()
        c.Modified()
        for fs in args.formats:
            savename = f'/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/INCLUSIVE_TrgEff/mcTT4Q_newTriggerLogic_inclTrgEff_{var}{fs}'
            #savename = f'/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/INCLUSIVE_TrgEff/mcWeightedAll_JECs_inclTrgEff_{var}{fs}'
            c.SaveAs(savename)
            

if __name__ == "__main__":

    VERBOSE       = True
    YEAR          = "2024"
    FORMATS       = ['.png', '.pdf']

    parser = ArgumentParser(description="Derive the trigger scale factors")
    parser.add_argument("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true", help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))
    parser.add_argument("--year", dest="year", action="store", default=YEAR, help="Process year")
    parser.add_argument("--formats", dest="formats", default=FORMATS, action="store", help="Formats to save histograms")

    args = parser.parse_args()
    main(args)








    
