#!/usr/bin/env python3
'''
DESCRIPTION:


'''
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
    d.SetLeftMargin(0.13)
    return d

def AddPrivateWorkText(setx=0.21, sety=0.905):
    tex = ROOT.TLatex(0.,0.,'Simulation Preliminary'); #'Private Work'); #'HLT tutorial');
    tex.SetNDC();
    tex.SetX(setx);
    tex.SetY(sety);
    tex.SetTextFont(53);
    tex.SetTextSize(28);
    tex.SetLineWidth(2)
    return tex

def AddCMSText(setx=0.205, sety=0.905):
    texcms = ROOT.TLatex(0.,0., 'CMS');
    texcms.SetNDC();
    texcms.SetTextAlign(31);
    texcms.SetX(setx);
    texcms.SetY(sety);
    texcms.SetTextFont(63);
    texcms.SetLineWidth(2);
    texcms.SetTextSize(30);
    return texcms

def createLegend():
    legend = ROOT.TLegend(0.63,0.20, 0.78, 0.65) #(0.58, 0.20, 0.82, 0.35)
    legend.SetFillColor(0)
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextSize(0.025);
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

colors = {0: ROOT.kBlack,
          1: ROOT.kBlue,
          2: ROOT.kGreen+1,
          3: ROOT.kRed+1,
          4: ROOT.kOrange-3,
          5: ROOT.kMagenta+2,
          6: ROOT.kTeal+3,
          7: ROOT.kYellow+1,
          8: ROOT.kViolet-2,
          9: ROOT.kCyan+1,
          10: ROOT.kPink+2,
          11: ROOT.kSpring+3,
          12: ROOT.kAzure+1,
          13: ROOT.kGray+1,
          14: ROOT.kOrange+7,
          15: ROOT.kMagenta-4,
          16: ROOT.kTeal-2,
          17: ROOT.kBlue-7,
          18: ROOT.kRed-7,
          19: ROOT.kGreen-7,
          20: ROOT.kAzure+7,
          21: ROOT.kSpring-5,
          22: ROOT.kPink-7,
          23: ROOT.kViolet+6,
          24: ROOT.kCyan-3,
          25: ROOT.kYellow-7,
          26: ROOT.kOrange+3,
          27: ROOT.kMagenta+1,
          28: ROOT.kTeal+1,
          29: ROOT.kBlue+3,
          }
          
pathway_L1 = "DATA_test_L1/bkg_L1/"

def main(args):

    #Loop on mass number variable
    #mass = [300, 200, 	#100, 	50]
    #label = ["Partial_2024I", "TTto4Q_s_PDID", "TTto4Q_s_L1"] #"TTto4Q_s_UNPS"] 
    bkg = ["2024I", "QCD_100to200", "QCD_200to400", "QCD_400to600", "QCD_600to800", "QCD_800to1000", "QCD_1000to1200", "QCD_1200to1500", "QCD_1500to2000", "QCD_2000", "TTto4Q", "WWto4Q"]
    #["Partial_2024I", "QCD_100to200", "QCD_200to400", "QCD_400to600", "QCD_600to800", "QCD_800to1000", "QCD_1000to1200", "QCD_1200to1500", "QCD_1500to2000", "QCD_2000", "TTto4Q", "WWto4Q"]
     
    variables  = ["ht_inclusive", "pt_leading"] #["minv_1", "minv_2", "minv_3"] #, "pv"] After pass ref trigger
    triggers   = ["passed"] #, "passedEXC", "passedL1HTT280", "passedL1Jet", "passedL1HTDouble", "passedL1HTSingle"] #HLT or DST (only the las tone now)
    statOption = ROOT.TEfficiency.kFCP

    for var in variables:
        #print(args.rfile)
        #print(type(args.rfile))
        c = getCanvas()
        leg = createLegend()
	
        print(var)
        nums = {}
        effs = {}
        for k, (m, arg) in enumerate(zip(bkg, args.rfile)): #MASS OR BKG
            #print(arg)
            #print(type(arg))
            #print(l)
            #for j, trg in enumerate(triggers): #triggers
            f = ROOT.TFile(arg, "READ")
            fdir = f.GetDirectory("InclusiveTrigNanoAOD")#(str(m)+"_DijethtTrigAnalyzerNanoAOD_AN_Ortogonal") #_+str(m))
            # Pass Ref Trigger & Offline Selection
            den = fdir.Get(f'h_{var}_all')
                
            # Pass Ref Trigger, Signal Trigger & Offline Selection
            nums[m] = fdir.Get(f'h_{var}_passed')
            #num2s[trg] = pdir.Get(f'h_{var}_{trg}') 
            try:
                effs[m] = ROOT.TEfficiency(nums[m], den)
                effs[m].SetStatisticOption(statOption)
                effs[m] = SetStyle(effs[m], colors[k])
            except Exception as e:
                 continue
            #effs2[trg] = ROOT.TEfficiency(num2s[trg], pen)
            #effs2[trg].SetStatisticOption(statOption)
            #effs2[trg] = SetStyle(effs2[trg], colors[j+1])
            if k == 0:
                 effs[m].Draw()
                 #leg.AddEntry(effs[m], triggers[k].replace(str(trg), str(m)))
                 #continue
                 #effs2[trg].Draw("same")
            else:
                 effs[m].Draw("same")
            #leg.AddEntry(effs[trg], "OR) #"DST_PFScouting_JetHT_"    
            #leg.AddEntry(effs[m], triggers[j].replace("passed", str(m)))	#+str(m))) #+"_GEN"))
            leg.AddEntry(effs[m], triggers[0].replace("passed", str(m)))#"DST_PFScouting_JetHT_v"))
                      
        c.Modified()
        c.Update()
        #effs["pIassedL1"].GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.2)
        leg.Draw("same")
            
        # Styling stuff
        tex_cms = AddCMSText()
        tex_cms.Draw("same")
        
        private = AddPrivateWorkText()
        private.Draw("same")
        
        header = ROOT.TLatex()
        header.SetTextSize(0.04)
        header.DrawLatexNDC(0.63, 0.905, "2024 (13.6 TeV)") #0.57 0.905
        
        #c.Update()
        c.Modified()
        for fs in args.formats:
           savename = f'plots/TrgEffs_DijetHT_code_{var}{fs}' #DATA+BKG/TrgEffs_DijetHT_Minv_{var}{fs}'
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
    private = AddPrivateWorkText()
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
    # Post Processor histFileName
    TRGROOTFILE   = [ 
"histos_Data2024_inclTrigEff_wMuCleaning_GoldenJSON.root",
"histos_Data2024_inclTrigEff_wMuCleaning_GoldenJSON.root",
"histosQCD_HT_100to200_InclusiveTrigNanoAOD.root",
"histosQCD_HT_200to400_InclusiveTrigNanoAOD.root",
"histosQCD_HT_400to600_InclusiveTrigNanoAOD.root",
"histosQCD_HT_600to800_InclusiveTrigNanoAOD.root",
"histosQCD_HT_800to1000_InclusiveTrigNanoAOD.root",
"histosQCD_HT_1000to1200_InclusiveTrigNanoAOD.root",
"histosQCD_HT_1200to1500_InclusiveTrigNanoAOD.root",
"histosQCD_HT_1500to2000_InclusiveTrigNanoAOD.root",
"histosQCD_HT_2000_InclusiveTrigNanoAOD.root",
#"histosTTto4Qfull_InclusiveTrigNanoAOD.root",
"histosTTto4Qs_InclusiveTrigNanoAOD.root",
"histosWWto4Q_InclusiveTrigNanoAOD.root"]
#"histos_DijetHTTrigNanoAOD_AN_Ortogonal_Partial_2024I.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_Partial_2024I_PFID.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_only2replace.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_bothor.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_both.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_PFID.root"]
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_L1bitvalidation.root"]
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_L1bitvalidation2.root"]
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_only1or.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_only1.root",
#"DATA_test_L1/histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q_s_only2.root"
#]
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_100to200.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_200to400.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_400to600.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_600to800.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_800to1000.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_1000to1200.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_1200to1500.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_1500to2000.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_QCD_2000.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_TTto4Q.root",
#    pathway_L1+"histos_DijetHTTrigNanoAOD_AN_Ortogonal_WWto4Q.root"]
 
# SIGNAL
#["histos_DijetHTTrigNanoAOD_AN_Ortogonal_300.root", "histos_DijetHTTrigNanoAOD_AN_Ortogonal_200.root", #"histos_DijetHTTrigNanoAOD_AN_Ortogonal_100GEN.root", "histos_DijetHTTrigNanoAOD_AN_Ortogonal_50.root"]
    #"histos_DijetHTTrigNanoAOD_AN_Ortogonal.root"
    FORMATS       = ['.png', '.pdf']

    parser = ArgumentParser(description="Derive the trigger scale factors")
    parser.add_argument("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true", help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))
    parser.add_argument("--rfile", dest="rfile", type=str, action="store", default=TRGROOTFILE, help="ROOT file containing the denominators and numerators [default: %s]" % (TRGROOTFILE))
    parser.add_argument("--year", dest="year", action="store", default=YEAR, help="Process year")
    parser.add_argument("--formats", dest="formats", default=FORMATS, action="store", help="Formats to save histograms")

    args = parser.parse_args()
    main(args)
