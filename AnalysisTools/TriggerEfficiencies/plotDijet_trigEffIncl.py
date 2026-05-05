#***********************************************
# Run 3 dijet scouting analysis:               #
# trigger efficiency studies                   #
# **********************************************

#!/usr/bin/env python3

#================
# Import modules
#================
from argparse import ArgumentParser
import ROOT
import math 
import time

def getCanvas():
    d = ROOT.TCanvas("", "", 800, 700)
    d.SetLeftMargin(0.12)
    d.SetRightMargin(0.15)
    d.SetLeftMargin(0.14)
    d.SetBottomMargin(0.15)
    return d

def AddText(setx=0.243, sety=0.91):
    tex = ROOT.TLatex(0.,0.,'Preliminary');
    tex.SetNDC();
    tex.SetX(setx);
    tex.SetY(sety);
    tex.SetTextFont(53);
    tex.SetTextSize(33);
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

def SetStyle(h, color, marker_style):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(marker_style)
    return h

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)

def SetStyle(h, COLOR):
    h.SetMarkerStyle(21)
    h.SetMarkerColor(COLOR)
    h.SetLineColor(COLOR)
    return h

colors = {0: ROOT.kAzure-3, #ROOT.kBlack,
          1: ROOT.kBlue,
          2: ROOT.kGreen+1,
          3: ROOT.kRed+1,
          4: ROOT.kOrange-3,
          5: ROOT.kMagenta+2,
          6: ROOT.kTeal+3,
          }
          

def main(args):

    f = ROOT.TFile(args.rfile, "READ")
    # Post Processor histDirName
    #fdir = f.GetDirectory("DijethtTrigAnalyzerNanoAOD")
    fdir = f.GetDirectory("InclusiveTrigNanoAOD")
    statOption = ROOT.TEfficiency.kFCP

    variables  = ["ht_inclusive", "pt_leading", "dijetMass"]
    triggers   = ["passed"] #HLT or DST (only the las tone now)

    for var in variables:
        c = getCanvas()

        # ----- Build legend -----
        leg = createLegend()

        # ----- Build numerator and denominator -----
        den = fdir.Get(f'h_{var}_all')
        nums = {}
        effs = {}

        for j, trg in enumerate(triggers):
            nums[trg] = fdir.Get(f'h_{var}_{trg}')
            effs[trg] = ROOT.TEfficiency(nums[trg], den)
            effs[trg].SetStatisticOption(statOption)
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
            line.SetLineStyle(2)                     # dashed
            line.SetLineWidth(2)
            line.SetLineColor(ROOT.kAzure+6)
            line.Draw("same")
    
        # ----- Styling stuff: add texts ------
    
        leg.Draw("same")

        tex_cms = AddCMSText()
        tex_cms.Draw("same")
        
        addText = AddText()
        addText.Draw("same")
        
        header = ROOT.TLatex()
        header.SetTextSize(0.045)
        header.DrawLatexNDC(0.581, 0.91, "2024H (13.6 TeV)")
        #header.DrawLatexNDC(0.57, 0.905, "2024I, #sqrt{s} = 13 TeV")

        c.Update()
        c.Modified()
        for fs in args.formats:
            #savename = f'/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/TEST_HT/inclTrigEffTightMuSelectionPlusMuCleaning_TriggEff_{var}{fs}'
            savename = f'/eos/user/e/elfontan/www/dijetAnaRun3/TRIGGER_EFF/INCLUSIVE_TrgEff/inclTrgEff_2024H_newL1TrgLogic_TriggEff_{var}{fs}'
            c.SaveAs(savename)
            
    '''
    # Plot 2D:
    c2D = getCanvas()

    den = fdir.Get('h_pfht_vs_pv_all')
    num = fdir.Get('h_pfht_vs_pv_passed')
    
    eff2D = ROOT.TEfficiency(num, den)
    eff2D.Draw("COLZ")
    c2D.Modified()
    c2D.Update()
    tex_cms = AddCMSText()
    tex_cms.Draw("same")
    private = AddText()
    private.Draw("same")
    header = ROOT.TLatex()
    header.SetTextSize(0.04)
    header.DrawLatexNDC(0.57, 0.905, "2023, #sqrt{s} = 13.6 TeV")
    c2D.Update()
    c2D.Modified()
    for fs in args.formats:
        c2D.SaveAs("plots/Eff2D_PFHTtrigger_PFHTvsPV%s" % (fs))
    '''
if __name__ == "__main__":

    VERBOSE       = True
    YEAR          = "2024"
    # ----- Histograms name
    TRGROOTFILE   = "/eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_inclTrgEffData_2024H_newL1TrgLogic.root"
    #TRGROOTFILE   = "/eos/cms/store/cmst3/user/elfontan/scoutAna/TriggerEff/histos_inclTrgEff_JECs_data2024H.root"
    FORMATS       = ['.png', '.pdf']

    parser = ArgumentParser(description="Derive the trigger scale factors")
    parser.add_argument("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true", help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))
    parser.add_argument("--rfile", dest="rfile", type=str, action="store", default=TRGROOTFILE, help="ROOT file containing the denominators and numerators [default: %s]" % (TRGROOTFILE))
    parser.add_argument("--year", dest="year", action="store", default=YEAR, help="Process year")
    parser.add_argument("--formats", dest="formats", default=FORMATS, action="store", help="Formats to save histograms")

    args = parser.parse_args()
    main(args)








    
