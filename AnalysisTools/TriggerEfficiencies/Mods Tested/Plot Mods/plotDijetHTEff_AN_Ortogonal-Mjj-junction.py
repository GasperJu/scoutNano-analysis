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
from ctypes import c_double

def getCanvas():
    d = ROOT.TCanvas("", "", 800, 800)
    d.SetLeftMargin(0.12)
    d.SetRightMargin(0.15)
    d.SetLeftMargin(0.13)
    return d

def pad1Canvas():
    pad1 = ROOT.pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.032)  # Small gap between pads
    #pad1.SetLeftMargin(0.12)
    #pad1.SetRightMargin(0.15)
    #pad1.SetLeftMargin(0.13)
    return pad1

def pad2Canvas():
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.35)  # Space for x-axis labels
    #pad2.SetGridy()  # Horizontal grid line at y=1
    pad2.SetGrid(1,1)
    return pad2
'''
def AddPrivateWorkText(setx=0.185, sety=0.935, t = text): #0.65, 0.935
    tex = ROOT.TLatex(0.,0., text); # 'Simulation Premilinary'; 'Private Work'); #'HLT tutorial');
    tex.SetNDC();
    tex.SetX(setx);
    tex.SetY(sety);
    tex.SetTextFont(53);
    tex.SetTextSize(28);
    tex.SetLineWidth(2)
    return tex
'''

def AddCMSText(setx=0.175, sety=0.935):
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
    legend = ROOT.TLegend(0.63,0.50, 0.78, 0.65) #(0.58, 0.20, 0.82, 0.35)
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

def AddPrivateWorkText(setx=0.185, sety=0.935, text = 'Private Work'): #0.65, 0.935
    tex = ROOT.TLatex(0.,0., text); # 'Simulation Premilinary'; 'Private Work'); #'HLT tutorial');
    tex.SetNDC();
    tex.SetX(setx);
    tex.SetY(sety);
    tex.SetTextFont(53);
    tex.SetTextSize(28);
    tex.SetLineWidth(2)
    return tex

def ratio_print(eff_lower, number, var):
   
    eff_lower.SetLineColor(colors[number])
    eff_lower.SetMarkerColor(colors[number])
    eff_lower.SetMarkerStyle(21)
    eff_lower.SetMarkerSize(0.8)    
 
    # Style lower pad
    eff_lower.SetTitle("")
    eff_lower.GetYaxis().SetTitle("Data-MC")
    eff_lower.GetXaxis().SetTitle(var)
    eff_lower.GetYaxis().SetRangeUser(0.95, 1.05) 
    eff_lower.GetYaxis().SetNdivisions(505)
    eff_lower.GetYaxis().SetTitleSize(20)
    eff_lower.GetYaxis().SetTitleFont(43)
    eff_lower.GetYaxis().SetTitleOffset(1.5)
    eff_lower.GetYaxis().SetLabelFont(43)
    eff_lower.GetYaxis().SetLabelSize(15)
    eff_lower.GetXaxis().SetTitleSize(20)
    eff_lower.GetXaxis().SetTitleFont(43)
    eff_lower.GetXaxis().SetTitleOffset(0.8)
    eff_lower.GetXaxis().SetLabelFont(43)
    eff_lower.GetXaxis().SetLabelSize(15)

    return eff_lower

def create_ratio(graph1, graph2):

    # Create ratio TGraphAsymmErrors
    n_points = graph1.GetN()
    #print(n_points)
    ratio_graph = ROOT.TGraphAsymmErrors(n_points)

    for i in range(n_points):
        #x1, y1 = ROOT.Double(0), ROOT.Double(0)
        #x2, y2 = ROOT.Double(0), ROOT.Double(0)
        x = c_double()
        y1 = c_double()
        y2 = c_double()    

        graph1.GetPoint(i, x, y1)
        graph2.GetPoint(i, x, y2)
        if y2.value > 0 and y1.value > 0:
           ratio = y1.value / y2.value
        
           # Get errors
           err1_low = graph1.GetErrorYlow(i)
           err1_high = graph1.GetErrorYhigh(i)
           err2_low = graph2.GetErrorYlow(i)
           err2_high = graph2.GetErrorYhigh(i)
        
           # Error propagation for ratio
           err_low = ratio * ROOT.TMath.Sqrt((err1_low/y1.value)**2 + (err2_low/y2.value)**2)
           err_high = ratio * ROOT.TMath.Sqrt((err1_high/y1.value)**2 + (err2_high/y2.value)**2)
        
           ratio_graph.SetPoint(i, x.value, ratio)
           ratio_graph.SetPointError(i, graph1.GetErrorXlow(i), graph1.GetErrorXhigh(i),
                                   err_low, err_high)
        else:
           # Handle zero or negative values
           ratio_graph.SetPoint(i, x.value, 1.0)
           ratio_graph.SetPointError(i, 0, 0, 0, 0)
    
    # Style and draw ratio
    ratio_graph.SetMarkerStyle(20)
    ratio_graph.SetMarkerSize(0.8)

    # Style ratio plot
    ratio_graph.SetTitle("")
    ratio_graph.GetYaxis().SetTitle("Eff Data/Eff BKG")
    ratio_graph.GetYaxis().SetRangeUser(0.5, 0.8)
    ratio_graph.GetYaxis().SetNdivisions(505)
    ratio_graph.GetYaxis().SetTitleSize(20)
    ratio_graph.GetYaxis().SetTitleFont(43)
    ratio_graph.GetYaxis().SetTitleOffset(1.5)
    ratio_graph.GetYaxis().SetLabelFont(43)
    ratio_graph.GetYaxis().SetLabelSize(15)
    ratio_graph.GetXaxis().SetTitleSize(20)
    ratio_graph.GetXaxis().SetTitleFont(43)
    ratio_graph.GetXaxis().SetLabelFont(43)
    ratio_graph.GetXaxis().SetLabelSize(15)

    return ratio_graph

pathway_L1 = "DATA_test_L1/bkg_L1/"

def main(args):

    #Loop on mass number variable
    #mass = [300, 200, 	#100, 	50]
    bkg = ["Partial_2024I", "QCD_100to200", "QCD_200to400", "QCD_400to600", "QCD_600to800", "QCD_800to1000", "QCD_1000to1200", "QCD_1200to1500", "QCD_1500to2000", "QCD_2000", "TTto4Q", "WWto4Q"]
    ### ADD WEIGHT (xsec took from samples codes)
    weights = [1.0, 25360000.0, 1951000.0, 96660.0, 13684.0, 3047.0, 889.0, 388.1, 127.2, 26.5, 762.1, 16100.0]            #50.8]

    variables  = ["ht_inclusive", "pt_leading"]
    #["minv_1", "minv_2", "minv_3"] #, "pv"] After pass ref trigger
    triggers   = ["passed"] #HLT or DST (only the las tone now)
    statOption = ROOT.TEfficiency.kFCP
   
    c = getCanvas()
    #histog_all = ROOT.THStack("histog_all","")
    #histog_pass = ROOT.THStack("histog_pass","")
    for var in variables:
        #print(args.rfile)
        #print(type(args.rfile))
        pad1 = pad1Canvas()
        leg = createLegend()
        pad1.Draw()
        pad1.cd()	

        nums = {}
        effs = {}
        #histog_all = {}
        #histog_pass = {}
        for k, (m, arg) in enumerate(zip(bkg, args.rfile)): #MASS OR BKG
            #print(arg)
            #print(type(arg))
            
            f = ROOT.TFile(arg, "READ")
            fdir = f.GetDirectory("InclusiveTrigNanoAOD")
            #str(m)+"_DijethtTrigAnalyzerNanoAOD_AN_Ortogonal") #_+str(m))
            # Pass Ref Trigger & Offline Selection
            den = fdir.Get(f'h_{var}_all')
            #print(den.ls()) 
            #print({den.Class().GetName()})
            nums[m] = fdir.Get(f'h_{var}_passed')            

            # DATA
            if k == 0:
               print(m)
               #print({nums.Class()})
               #print({nums.ls()})
               #den = fdir.Get(f'h_{var}_all')
               nums[m] = fdir.Get(f'h_{var}_passed')
               effs[m] = ROOT.TEfficiency(nums[m], den)
               
               #PAD
               #print(effs[m])
               hist_data_x = effs[m].GetTotalHistogram()
               graph1 = effs[m].CreateGraph()

               effs[m].SetStatisticOption(statOption)
               effs[m] = SetStyle(effs[m], colors[k])
               effs[m].Draw()
               leg.AddEntry(effs[m], triggers[0].replace("passed", str("2024I")))
               #effs2[trg].Draw("same")
               #break
               f.Close()
            # BKG ADDED UP
            elif k == len(bkg)-1:
               print(m)
               print(k)
               #den =  histog_all #.Get(f'histog_pass')
               #nums[m] = histog_pass #.Get(f'histog_pass')
               effs[m] = ROOT.TEfficiency(histog_pass, histog_all)
               
               #PAD
               #graph1_temp = effs[m].GetTotalHistogram()
               graph2 = effs[m].CreateGraph()

               effs[m].SetStatisticOption(statOption)
               effs[m] = SetStyle(effs[m], colors[k])              
               effs[m].Draw("same")
               #leg.AddEntry(effs[trg], "OR) #"DST_PFScouting_JetHT_"
               leg.AddEntry(effs[m], triggers[0].replace("passed", str("Background")))  #+str(m))) #+"_GEN"))
               #break
               f.Close()
            #CLONE BKG 1
            elif k == 1:
               print(m)
               histog_all = den.Clone("histog_all")
               histog_all.Scale(weights[k])
               histog_all.SetDirectory(0)  # Detach from file
             
               #nums[m] = fdir.Get(f'h_{var}_passed')
               histog_pass = nums[m].Clone("histog_pass")
               histog_pass.Scale(weights[k])
               histog_pass.SetDirectory(0)  # Detach from file
               f.Close()
            # ADDING BKG
            else: #Probably not necessary
               #print(m)
               h_temp_all = den.Clone()
               h_temp_all.Scale(weights[k])
               #if k!= 0 or k != len(bkg)-1:
               histog_all.Add(h_temp_all)
               histog_all.SetDirectory(0)  # Detach from file

               h_temp_pass = nums[m].Clone()
               h_temp_pass.Scale(weights[k])
               #if k!= 0 or k != len(bkg)-1:
               histog_pass.Add(h_temp_pass)
               histog_pass.SetDirectory(0)  # Detach from file
               f.Close()
        
        #SAME X RANGE
        pad1.Update()
        xmax = graph1.GetXaxis().GetXmax()
        #x_title = graph1.GetXaxis().GetTitle()
        graph1.GetXaxis().SetLimits(0, xmax)
        graph1.GetXaxis().SetRangeUser(0, xmax)
        hist_data_x.GetYaxis().SetTitle('Efficiency')
        pad1.Modified()
        pad1.Update()
        
        #c.Update()
        # Lower pad for ratio
        c.cd()       
        pad2 = pad2Canvas()
        pad2.Draw()
        pad2.cd()
        #xmax = graph1_temp.GetXaxis().GetXmax()
        # Create ratio - convert TEfficiency to TGraphAsymmErrors
        #graph1 = effs[0].CreateGraph() #DATA
        #graph2 = effs[1].CreateGraph() #BKG
        pad_graph1 = ratio_print(graph1, 0, var)
        pad_graph2 = ratio_print(graph2, len(bkg)-1, var)
        # Style ratio plot
        pad_graph1.Draw("ap")
        pad_graph2.Draw("p same")
        
        x_title = hist_data_x.GetXaxis().GetTitle()
        #print(x_title)
        pad_graph1.GetXaxis().SetTitle(x_title)
        pad2.Modified()
        pad2.Update()

        # Set same x-axis range as upper pad
        #pad_graph1.GetXaxis().SetLimits(0, xmax)
        #pad_graph1.GetXaxis().SetRangeUser(0, xmax)
        #pad_graph1.GetXaxis().SetLimits(0, xmax)
        #pad_graph1.GetXaxis().SetRangeUser(0, xmax)
        
        # Reference line at 1
        line = ROOT.TLine(0, 1, xmax, 1)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kBlack)
        line.Draw("SAME")                         

        c.cd() 
        #pad1.GetYaxis().SetTitle('Efficiency')
        c.Update()
        c.Modified()
        #effs["passedL1"].GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.2)
        leg.Draw("same")

        # Styling stuff
        tex_cms = AddCMSText()
        tex_cms.Draw("same")

        private = AddPrivateWorkText(text = 'Preliminary')
        private.Draw("same")

        header = ROOT.TLatex()
        header.SetTextSize(0.04)
        header.DrawLatexNDC(0.65, 0.935, "2024 (13.6 TeV)") #0.57 0.905

        c.Update()
        c.Modified()
        for fs in args.formats:
           savename = f'plots/DATA+COMPBKG/TrgEffs_DijetHT_{var}_ratio_{fs}'
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
    YEAR          = "2024I"
    # Post Processor histFileName
    TRGROOTFILE   = [
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

# SIGNAL
#["histos_DijetHTTrigNanoAOD_AN_Ortogonal_300.root", "histos_DijetHTTrigNanoAOD_AN_Ortogonal_200.root", #"histos_DijetHTTrigNanoAOD_AN_Ortogonal_100GEN.root", "histos_DijetHTTrigNanoAOD_AN_Ortogonal_50.root"
    #"histos_DijetHTTrigNanoAOD_AN_Ortogonal.root"
    FORMATS       = ['.png', '.pdf']

    parser = ArgumentParser(description="Derive the trigger scale factors")
    parser.add_argument("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true", help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))
    parser.add_argument("--rfile", dest="rfile", type=str, action="store", default=TRGROOTFILE, help="ROOT file containing the denominators and numerators [default: %s]" % (TRGROOTFILE))
    parser.add_argument("--year", dest="year", action="store", default=YEAR, help="Process year")
    parser.add_argument("--formats", dest="formats", default=FORMATS, action="store", help="Formats to save histograms")

    args = parser.parse_args()
    main(args)
