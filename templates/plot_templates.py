import os, sys
import optparse
import ROOT
from tdrStyle import *
from binning import binning

sPlotsStore = 'plots'

print 'Welcome in plot_templates!'

# Considering or not decimals in bin boundaries
decimal = {
'mass4l': False,
'mass4l_zzfloating': False,
'njets_pt30_eta4p7': False,
'pT4l': False,
'pT4l_kL': False,
'rapidity4l': True,
'costhetaZ1': True,
'costhetaZ2': True,
'phi': True,
'phistar': True,
'costhetastar': True,
'massZ1': False,
'massZ2': False,
'pTj1': False,
'pTHj': False,
'mHj': False,
'pTj2': False,
'mjj': False,
'absdetajj': True,
'dphijj': True,
'pTHjj': False,
'TCjmax': False,
'TBjmax': False,
'D0m': True,
'Dcp': True,
'D0hp': True,
'Dint': True,
'DL1': True,
'DL1Zg': True,
'rapidity4l vs pT4l': True,
'njets_pt30_eta4p7 vs pT4l': False,
'pTj1 vs pTj2': False,
'pT4l vs pTHj': False,
'massZ1 vs massZ2': False,
'TCjmax vs pT4l': False
}

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='int',default=105,   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='int',default=140,   help='Upper bound for m4l')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

global opt, args, runAllSteps
parseOptions()

def setHistProperties(hist, lineWidth, lineStyle, lineColor, fillStyle=0, fillColor=0, xAxisTitle = "skip", yAxisTitle = "skip"):
    if not hist: return
    hist.SetLineWidth(lineWidth)
    hist.SetLineStyle(lineStyle)
    hist.SetLineColor(lineColor)
    hist.SetFillStyle(fillStyle)
    hist.SetFillColor(fillColor)
    hist.GetXaxis().SetNdivisions(510)
    hist.GetYaxis().SetNdivisions(510)
    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.24)
    if (xAxisTitle!="skip"): hist.GetXaxis().SetTitle(xAxisTitle)
    if (yAxisTitle!="skip"): hist.GetYaxis().SetTitle(yAxisTitle)

def cmsPreliminary(c, top):
    c.cd()

    CMSPrelim = ROOT.TLatex()
    CMSPrelim.SetNDC(ROOT.kTRUE)

    CMSPrelim.SetTextSize(0.8*c.GetTopMargin());
    CMSPrelim.SetTextFont(42)
    CMSPrelim.SetTextAlign(31) #align right
    CMSPrelim.DrawLatex(0.93, 0.96,"#bf{"+top+"}")

def setLegendProperties(leg, sHeader = "skip", fillStyle=0, fillColor=0):
    # sanity-check
    if not leg: return
    # titles
    if (sHeader!="skip"): leg.SetHeader(sHeader)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)

obs_bins, doubleDiff = binning(opt.OBSNAME)

#Prepare all strings for following steps
obsTag = opt.OBSNAME

if doubleDiff:
    N_BINS = len(obs_bins)
    for i in range(N_BINS):
        if not decimal[obsTag]: obs_bins[i] = [int(k) for k in obs_bins[i]]
    binRange = [str(obs_bins[i][0])+'_'+str(obs_bins[i][1])+'_'+str(obs_bins[i][2])+'_'+str(obs_bins[i][3]) for i in range(N_BINS)]
    binRangeLow = [str(obs_bins[i][0])+'_'+str(obs_bins[i][2]) for i in range(N_BINS)]
    binRangeHigh = [str(obs_bins[i][1])+'_'+str(obs_bins[i][3]) for i in range(N_BINS)]
    obsTag = obsTag.split(' vs ')[0]+'_'+obsTag.split(' vs ')[1]
    binRangeLeg = [str(obs_bins[i][0])+'<'+obsTag+'<'+str(obs_bins[i][1])+'/'+str(obs_bins[i][2])+'<'+obsTag+'<'+str(obs_bins[i][3]) for i in range(N_BINS)]
else:
    N_BINS = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)
    if not decimal[obsTag]: obs_bins = [int(k) for k in obs_bins]
    binRange = [str(obs_bins[i])+'_'+str(obs_bins[i+1]) for i in range(N_BINS)]
    binRangeLow = [str(obs_bins[i]) for i in range(N_BINS)]
    binRangeHigh = [str(obs_bins[i+1]) for i in range(N_BINS)]
    binRangeLeg = [str(int(obs_bins[i]))+'<'+obsTag+'<'+str(int(obs_bins[i+1])) for i in range(N_BINS)]

print(binRange)
print(binRangeLow)
print(binRangeHigh)
print(binRangeLeg)


# setup environment & canvas
setTDRStyle()
c1 = ROOT.TCanvas('c1',"myPlots",0,0,800,800)
c1.cd(1)
c1.SetLogy(0)
gStyle.SetOptStat('')
gStyle.SetPalette(1)
c1.GetPad(0).SetRightMargin(0.05)
c1.GetPad(0).SetLeftMargin(0.15)
c1.GetPad(0).SetTopMargin(0.05)
c1.GetPad(0).SetBottomMargin(0.15)

bkgName=['qqzz','ggzz','ZJetsCR']
if opt.YEAR == 'Full':
    year=['2016', '2017', '2018']
else:
    year=[opt.YEAR]
print(year)

for iYear in range(len(year)):
    sTemplateDirName = year[iYear]+'/'+obsTag
    fTemplateFile_2e2mu = {}
    fTemplateFile_4mu = {}
    fTemplateFile_4e = {}
    h1D_2e2mu = {}
    h1D_4mu = {}
    h1D_4e = {}

    for iBin in range(N_BINS):
        for iBkg in range(len(bkgName)):

            sTemplateFileName = "XSBackground_"+bkgName[iBkg]+"_2e2mu_"+obsTag+"_"+binRange[iBin]+".root"
            fTemplateFile_2e2mu[iBkg,iBin] = ROOT.TFile(sTemplateDirName+"/"+sTemplateFileName, "READ")
            h1D_2e2mu[iBkg,iBin] = ROOT.TH1D()
            h1D_2e2mu[iBkg,iBin] = fTemplateFile_2e2mu[iBkg,iBin].Get("m4l_"+obsTag+"_"+binRange[iBin])
            print(sTemplateFileName)
            print("m4l_"+obsTag+"_"+binRange[iBin])

            sTemplateFileName = "XSBackground_"+bkgName[iBkg]+"_4mu_"+obsTag+"_"+binRange[iBin]+".root"
            fTemplateFile_4mu[iBkg,iBin] = ROOT.TFile(sTemplateDirName+"/"+sTemplateFileName, "READ")
            h1D_4mu[iBkg,iBin] = ROOT.TH1D()
            h1D_4mu[iBkg,iBin] = fTemplateFile_4mu[iBkg,iBin].Get("m4l_"+obsTag+"_"+binRange[iBin])

            sTemplateFileName = "XSBackground_"+bkgName[iBkg]+"_4e_"+obsTag+"_"+binRange[iBin]+".root"
            fTemplateFile_4e[iBkg,iBin] = ROOT.TFile(sTemplateDirName+"/"+sTemplateFileName, "READ")
            h1D_4e[iBkg,iBin] = ROOT.TH1D()
            h1D_4e[iBkg,iBin] = fTemplateFile_4e[iBkg,iBin].Get("m4l_"+obsTag+"_"+binRange[iBin])

    # prepare dummy
    var_plotHigh = opt.UPPER_BOUND
    var_plotLow = opt.LOWER_BOUND
    var_nBins = 20
    varAxLabel = 'm_{4l} (GeV)'
    binWidth = (int(100*(var_plotHigh - var_plotLow)/var_nBins))/100.
    h1D_dummy = ROOT.TH1D("dummy", "dummy", var_nBins, var_plotLow, var_plotHigh)
    setHistProperties(h1D_dummy,1,1,ROOT.kBlue-7,0,0,varAxLabel,"Events/"+str(binWidth)+'(GeV)')

    # common properties
    lineWidth = 2
    leg_xl = 0.52
    leg_xr = 0.90
    leg_yb = 0.72
    leg_yt = 0.90

    # plot hists
    kBkg_qqZZ = 0
    kBkg_ggZZ = 1
    kBkg_ZJets = 2
    c1.cd()
    for iBin in range(N_BINS):

        ########## 2e2mu ##########
        #### qqZZZ + ggZZ +ZX ####
        h1D_dummy.SetMaximum(2.0*h1D_2e2mu[kBkg_qqZZ,iBin].GetMaximum())
        h1D_dummy.Draw()
        cmsPreliminary(c1, binRangeLeg[iBin]+"      2e2#mu      "+year[iYear])
        leg1 = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)
        setLegendProperties(leg1)
        setHistProperties(h1D_2e2mu[kBkg_qqZZ,iBin],lineWidth,1,ROOT.kBlack)
        h1D_2e2mu[kBkg_qqZZ,iBin].Draw("histsame")
        leg1.AddEntry(h1D_2e2mu[kBkg_qqZZ,iBin], "q#bar{q} #rightarrow ZZ","L")
        setHistProperties(h1D_2e2mu[kBkg_ggZZ,iBin],lineWidth,1,ROOT.kBlue-7)
        h1D_2e2mu[kBkg_ggZZ,iBin].Draw("histsame")
        leg1.AddEntry(h1D_2e2mu[kBkg_ggZZ,iBin], "gg #rightarrow ZZ","L");
        setHistProperties(h1D_2e2mu[kBkg_ZJets,iBin],lineWidth,1,ROOT.kRed-7)
        h1D_2e2mu[kBkg_ZJets,iBin].Draw("histsame")
        leg1.AddEntry(h1D_2e2mu[kBkg_ZJets,iBin], "Z + X","L")
        leg1.Draw()
        c1.SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_2e2mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".pdf")
        c1.SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_2e2mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".png")

        ########## 4mu ##########
        #### qqZZZ + ggZZ +ZX ####
        h1D_dummy.SetMaximum(2.0*h1D_4mu[kBkg_qqZZ,iBin].GetMaximum())
        h1D_dummy.Draw()
        cmsPreliminary(c1, binRangeLeg[iBin]+"      4#mu      "+year[iYear])
        leg1 = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)
        setLegendProperties(leg1)
        setHistProperties(h1D_4mu[kBkg_qqZZ,iBin],lineWidth,1,ROOT.kBlack)
        h1D_4mu[kBkg_qqZZ,iBin].Draw("histsame")
        leg1.AddEntry(h1D_4mu[kBkg_qqZZ,iBin], "q#bar{q} #rightarrow ZZ","L")
        setHistProperties(h1D_4mu[kBkg_ggZZ,iBin],lineWidth,1,ROOT.kBlue-7)
        h1D_4mu[kBkg_ggZZ,iBin].Draw("histsame")
        leg1.AddEntry(h1D_4mu[kBkg_ggZZ,iBin], "gg #rightarrow ZZ","L");
        setHistProperties(h1D_4mu[kBkg_ZJets,iBin],lineWidth,1,ROOT.kRed-7)
        h1D_4mu[kBkg_ZJets,iBin].Draw("histsame")
        leg1.AddEntry(h1D_4mu[kBkg_ZJets,iBin], "Z + X","L")
        leg1.Draw()
        c1.SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".pdf")
        c1.SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".png")

        ########## 4e ##########
        #### qqZZZ + ggZZ +ZX ####
        h1D_dummy.SetMaximum(2.0*h1D_4e[kBkg_qqZZ,iBin].GetMaximum())
        h1D_dummy.Draw()
        cmsPreliminary(c1, binRangeLeg[iBin]+"      4e      "+year[iYear])
        leg1 = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)
        setLegendProperties(leg1)
        setHistProperties(h1D_4e[kBkg_qqZZ,iBin],lineWidth,1,ROOT.kBlack)
        h1D_4e[kBkg_qqZZ,iBin].Draw("histsame")
        leg1.AddEntry(h1D_4e[kBkg_qqZZ,iBin], "q#bar{q} #rightarrow ZZ","L")
        setHistProperties(h1D_4e[kBkg_ggZZ,iBin],lineWidth,1,ROOT.kBlue-7)
        h1D_4e[kBkg_ggZZ,iBin].Draw("histsame")
        leg1.AddEntry(h1D_4e[kBkg_ggZZ,iBin], "gg #rightarrow ZZ","L");
        setHistProperties(h1D_4e[kBkg_ZJets,iBin],lineWidth,1,ROOT.kRed-7)
        h1D_4e[kBkg_ZJets,iBin].Draw("histsame")
        leg1.AddEntry(h1D_4e[kBkg_ZJets,iBin], "Z + X","L")
        leg1.Draw()
        c1.SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4e_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".pdf")
        c1.SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4e_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".png")
