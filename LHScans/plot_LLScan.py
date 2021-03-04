from ROOT import *
import numpy as np
import math
from functools import partial
import plotting as plot
import json
import argparse, optparse
import os.path, sys

NAMECOUNTER = 0

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)

### Define function for parsing options
def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "mass4l", "pT4l", "massZ2", "rapidity4l", "cosThetaStar", "nets_reco_pt30_eta4p7"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

def read(scan, param, files, ycut):
    goodfiles = [f for f in files if plot.TFileIsGood(f)]
    limit = plot.MakeTChain(goodfiles, 'limit')
    graph = plot.TGraphFromTree(limit, param, '2*deltaNLL', 'quantileExpected > -1.5')
    graph.SetName(scan)
    graph.Sort()
    plot.RemoveGraphXDuplicates(graph)
    plot.RemoveGraphYAbove(graph, ycut)
    #graph.Print()
    return graph


def Eval(obj, x, params):
    return obj.Eval(x[0])


def BuildScan(scan, param, files, color, yvals, ycut):
    graph = read(scan, param, files, ycut)
    bestfit = None

    # for i in xrange(graph.GetN()):
    #     if graph.GetY()[i] < 0.003:
    #         bestfit = graph.GetX()[i]
    #
    # if graph.GetY()[i] == 0:
    #     bestfit = graph.GetX()[i]
    #
    # if bestfit is None:
    #     for i in xrange(graph.GetN()):
    #         if graph.GetY()[i] < 0.02:
    #                 bestfit = graph.GetX()[i]
    #
    # if bestfit is None:
    #     for i in xrange(graph.GetN()):
    #         if graph.GetY()[i] < 0.25:
    #                 bestfit = graph.GetX()[i]

    graph.SetMarkerColor(color)
    spline = TSpline3("spline3", graph)
    global NAMECOUNTER
    func = TF1('splinefn'+str(NAMECOUNTER), partial(Eval, spline), graph.GetX()[0], graph.GetX()[graph.GetN() - 1], 1)
    bestfit = func.GetMinimumX() #AT
    print bestfit
    NAMECOUNTER += 1
    func.SetLineColor(color)
    func.SetLineWidth(3)
    assert(bestfit is not None)
    crossings = {}
    cross_1sig = None
    cross_2sig = None
    other_1sig = []
    other_2sig = []
    val = None
    val_2sig = None
    for yval in yvals:
        crossings[yval] = plot.FindCrossingsWithSpline(graph, func, yval)
        for cr in crossings[yval]:
            cr["contains_bf"] = cr["lo"] <= bestfit and cr["hi"] >= bestfit
    for cr in crossings[yvals[0]]:
        if cr['contains_bf']:
            val = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
            cross_1sig = cr
        else:
            other_1sig.append(cr)
    if len(yvals) > 1:
        for cr in crossings[yvals[1]]:
            if cr['contains_bf']:
                val_2sig = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
                cross_2sig = cr
            else:
                other_2sig.append(cr)
    else:
        val_2sig = (0., 0., 0.)
        cross_2sig = cross_1sig
    print val
    return {
        "graph"     : graph,
        "spline"    : spline,
        "func"      : func,
        "crossings" : crossings,
        "val"       : val,
        "val_2sig": val_2sig,
        "cross_1sig" : cross_1sig,
        "cross_2sig" : cross_2sig,
        "other_1sig" : other_1sig,
        "other_2sig" : other_2sig
    }

yvals = [1., 4.]


channel = ["Expected","Expected - no syst.","Observed","Observed - no syst."]

inputPath = '../combine_files/'

fileList = [ "higgsCombine_BIN_OBS.MultiDimFit.mH125.38.123456.root",
            "higgsCombine_BIN_OBS_NoSys_exp.MultiDimFit.mH125.38.123456.root",
       ]

if(opt.UNBLIND):
    fileList = [ "higgsCombine_BIN_OBS.MultiDimFit.mH125.38.123456.root",
                "higgsCombine_BIN_OBS_NoSys_exp.MultiDimFit.mH125.38.123456.root",
                    "higgsCombine_BIN_OBS.MultiDimFit.mH125.38.root",
                    "higgsCombine_BIN_OBS_NoSys.MultiDimFit.mH125.38.root"
                   ]

outString = "test"

titles = ["Expected","Expected - stat-only"]
idx_max = 1

if(opt.UNBLIND):
    print('Using real data')
    idx_max = 3
    titles = ["Expected","Expected - stat-only","Observed","Observed - stat-only"]
colors = [kRed, kRed, kBlack, kBlack]

resultsXS_data = {}
resultsXS_asimov = {}
if opt.OBSNAME=='mass4l':
    resultsXS_data_v2 = {}
    resultsXS_asimov_v2 = {}

year = opt.YEAR

if year == '2016':
    _lumi = '35.9'
elif year == '2017':
    _lumi = '41.5'
elif year == '2018':
    _lumi = '59.7'
else:
    _lumi = '137'

_poi    = 'SigmaBin'
obsName = opt.OBSNAME

doubleDiff = False
if(obsName == 'mass4l'): label = 'm_{4l}'
if(obsName == 'rapidity4l'): label = '|y_{H}|'
if(obsName == 'pT4l'): label = 'p_{T}^{H} (GeV)'
if(obsName == 'massZ1'): label = 'm_{Z1} (GeV)'
if(obsName == 'massZ2'): label = 'm_{Z2} (GeV)'
if(obsName == 'njets_pt30_eta2p5'): label = 'nJets, pT>30 GeV, |#eta|<2.5'
if(obsName == 'pTj1'): label = 'p_{T}^{(Lead. jet)} (GeV)'
if(obsName == 'mass4l'): label = 'm_{4\ell} (GeV)'
if(obsName == 'costhetastar'): label = 'cos(#theta^{*})'
if(obsName == 'costhetaZ1'): label = 'cos(#theta_{1})'
if(obsName == 'costhetaZ2'): label = 'cos(#theta_{2})'
if(obsName == 'phistar'): label = '#Phi^{#star}'
if(obsName == 'phi'): label = '#Phi'
if(obsName == 'massZ1 vs massZ2'):
    obsName_tmp = obsName.split(' vs ')
    obsName = obsName_tmp[0]+"_"+obsName_tmp[1]
    label = 'm_{Z1} (GeV)'
    label_2nd = 'm_{Z2} (GeV)'
    doubleDiff = True
if(obsName == 'njets_pt30_eta2p5 vs pTHj'):
    obsName_tmp = obsName.split(' vs ')
    obsName = obsName_tmp[0]+"_"+obsName_tmp[1]
    label = 'N_{jet}'
    label_2nd = 'p_{T}^{Hj} (GeV)'
    doubleDiff = True

sys.path.append('../inputs')
_temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
obs_bins = _temp.observableBins
sys.path.remove('../inputs')

nBins = len(obs_bins)
if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries
if obsName=='mass4l': nBins = nBins + 3 #in case of mass4l len(obs_bins)=1, we need to add +3 for cross section in the three different final states

for i in range(nBins):
    _bin = i
    _obs_bin = _poi+str(i)

    if obsName == 'mass4l':
        if _bin == 1:
            _obs_bin = 'r2e2muBin0'
        if _bin == 2:
            _obs_bin = 'r4muBin0'
        if _bin == 3:
            _obs_bin = 'r4eBin0'

    graphs = []
    grapherrs = []
    grapherrslow = []

    for ifile in range(len(fileList)):
        rfile = fileList[ifile].replace('OBS', _obs_bin)
        rfile = rfile.replace('BIN', obsName)
        graphs.append(TGraph())
        fname = inputPath+rfile
        inF = TFile.Open(fname,"READ")
        tree = inF.Get("limit")

        ipoint = 0
        for entry in tree :
            if((2*entry.deltaNLL<5)):
                if _bin == 0:
                    graphs[ifile].SetPoint(ipoint,entry.SigmaBin0,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                elif _bin == 1:
                    if obsName == 'mass4l':
                        graphs[ifile].SetPoint(ipoint,entry.r2e2muBin0,2.0*entry.deltaNLL)
                    else:
                        graphs[ifile].SetPoint(ipoint,entry.SigmaBin1,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                elif _bin == 2:
                    if obsName == 'mass4l':
                        graphs[ifile].SetPoint(ipoint,entry.r4muBin0,2.0*entry.deltaNLL)
                    else:
                        graphs[ifile].SetPoint(ipoint,entry.SigmaBin2,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                elif _bin == 3:
                    if obsName == 'mass4l':
                        graphs[ifile].SetPoint(ipoint,entry.r4eBin0,2.0*entry.deltaNLL)
                    else:
                        graphs[ifile].SetPoint(ipoint,entry.SigmaBin3,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                elif _bin == 4:
                    graphs[ifile].SetPoint(ipoint,entry.SigmaBin4,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                elif _bin == 5:
                    graphs[ifile].SetPoint(ipoint,entry.SigmaBin5,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                elif _bin == 6:
                    graphs[ifile].SetPoint(ipoint,entry.SigmaBin6,2.0*entry.deltaNLL)
                    ipoint = ipoint+1
                else:
                    graphs[ifile].SetPoint(ipoint,entry.SigmaBin7,2.0*entry.deltaNLL)
                    ipoint = ipoint+1

    c=TCanvas("c", "c", 1000, 800)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.08)
    c.cd()
    gStyle.SetOptTitle(0)

    graphs[0].SetLineColor(colors[0])
    graphs[0].SetLineWidth(3)
    graphs[0].Sort()
    graphs[0].SetTitle(titles[0])
    graphs[0].SetFillColor(colors[0])
    graphs[0].SetFillStyle(3000)

    mini = graphs[idx_max].GetXaxis().GetXmin()
    maxi = graphs[idx_max].GetXaxis().GetXmax()
    if mini > graphs[0].GetXaxis().GetXmin():
        mini = graphs[0].GetXaxis().GetXmin()
    if maxi < graphs[0].GetXaxis().GetXmax():
        maxi = graphs[0].GetXaxis().GetXmax()
    maxY = 8.

    graphs[0].Draw("AC")
    xtitle = "#sigma_{bin " + str(_bin) + "}"
    graphs[0].GetXaxis().SetTitle(xtitle)
    graphs[0].GetYaxis().SetTitle("-2#Delta ln L")
    graphs[0].GetYaxis().SetTitleOffset(0.9)
    graphs[0].GetXaxis().SetTitleOffset(0.85)
    graphs[0].GetXaxis().SetTitleSize(0.05)
    graphs[0].GetYaxis().SetTitleSize(0.05)

    graphs[0].GetYaxis().SetRangeUser(0,maxY)
    graphs[0].GetXaxis().SetRangeUser(mini,maxi)
    c.Modified()

    graphs[0].GetYaxis().SetLimits(0,maxY)
    graphs[0].GetXaxis().SetLimits(mini,maxi)
    c.Modified()
    c.Update()

    for ig in range(1,len(graphs)) :
        graphs[ig].SetLineColor(colors[ig])
        graphs[ig].SetFillColor(colors[ig])
        graphs[ig].SetFillStyle(3000)
        graphs[ig].SetLineWidth(3)
        if 'stat-only' in titles[ig]:
            graphs[ig].SetLineWidth(2)
            graphs[ig].SetLineStyle(9)
        graphs[ig].SetTitle(titles[ig])
        graphs[ig].Sort()

        graphs[ig].Draw("CSAME")

    lineone = TLine(mini,1,maxi,1)
    linetwo = TLine(mini,4.0,maxi,4.0)
    lineone.SetLineColor(15)#kGray+1)
    linetwo.SetLineColor(15)#kGray+1)
    lineone.Draw("SAME")
    linetwo.Draw("SAME")

    leg = TLegend(0.65,0.85,0.9,0.75)
    leg.SetLineColor(0)
    leg.SetLineStyle(0)
    leg.SetLineWidth(0)
    #leg.SetFillStyle(0)
    leg.SetShadowColor(10)
    leg.SetTextSize(0.030)
    leg.SetTextFont(42)
    for ip in range(0, len(graphs)):
        leg.AddEntry(graphs[ip], titles[ip], "l")

    leg.Draw("SAME")

    poi = _obs_bin
    fname = inputPath + "higgsCombine_"+obsName+"_"+poi+".MultiDimFit.mH125.38.123456.root"
    exp_scan = BuildScan('scan', poi, [fname], 2, yvals, 7.)
    exp_nom = exp_scan['val']
    exp_2sig = exp_scan['val_2sig']

    fname = inputPath + "higgsCombine_"+obsName+"_"+poi+"_NoSys_exp.MultiDimFit.mH125.38.123456.root"
    exp_scan_stat = BuildScan('scan', poi, [fname], 2, yvals, 7.)
    exp_nom_stat = exp_scan_stat['val']
    exp_2sig_stat = exp_scan_stat['val_2sig']

    exp_up_sys = np.sqrt(exp_nom[1]**2 - exp_nom_stat[1]**2)
    exp_do_sys = np.sqrt(abs(exp_nom[2])**2 - abs(exp_nom_stat[2])**2)

    if (opt.UNBLIND):
        fname = inputPath + "higgsCombine_"+obsName+"_"+poi+".MultiDimFit.mH125.38.root"
        obs_scan = BuildScan('scan', poi, [fname], 2, yvals, 7.)
        obs_nom = obs_scan['val']
        obs_2sig = obs_scan['val_2sig']

        fname = inputPath + "higgsCombine_"+obsName+"_"+poi+"_NoSys.MultiDimFit.mH125.38.root"
        obs_scan_stat = BuildScan('scan', poi, [fname], 2, yvals, 7.)
        obs_nom_stat = obs_scan_stat['val']
        obs_2sig_stat = obs_scan_stat['val_2sig']
        obs_up_sys = np.sqrt(obs_nom[1]**2 - obs_nom_stat[1]**2)
        obs_do_sys = np.sqrt(abs(obs_nom[2])**2 - abs(obs_nom_stat[2])**2)

    if(opt.UNBLIND):
        Text3 = TPaveText(0.18, 0.81,0.4,0.9,'brNDC')
    else:
    	Text3 = TPaveText(0.18, 0.76,0.4,0.84,'bfNDC')
    exp_fit = 'Exp. #sigma_{bin, %d} = %.2f^{#plus %.2f}_{#minus %.2f} (stat)^{#plus %.2f}_{#minus %.2f} (syst)' % (_bin, exp_nom[0], exp_nom_stat[1], abs(exp_nom_stat[2]), exp_up_sys, exp_do_sys)
    Text3.SetTextAlign(12);
    Text3.SetTextSize(0.036)
    Text3.AddText(exp_fit)
    Text3.SetFillStyle(0)
    Text3.SetLineStyle(0)
    Text3.SetBorderSize(0)
    Text3.Draw()

    if(opt.UNBLIND):
        Text4 = TPaveText(0.18, 0.71,0.4,0.8,'brNDC')
        obs_fit = 'Obs. #sigma_{bin, %d} = %.2f^{#plus %.2f}_{#minus %.2f} (stat)^{#plus %.2f}_{#minus %.2f} (syst)' % (_bin, obs_nom[0], obs_nom_stat[1], abs(obs_nom_stat[2]), obs_up_sys, obs_do_sys)
        Text4.SetTextAlign(12);
        Text4.SetTextSize(0.036)
        Text4.AddText(obs_fit)
        Text4.SetFillStyle(0)
        Text4.SetLineStyle(0)
        Text4.SetBorderSize(0)
        Text4.Draw()

    Text = TPaveText(0.58, 0.88,0.93,0.95,'brNDC')
    #Text.SetNDC()
    Text.SetTextAlign(31);
    Text.SetTextSize(0.03)
    leftText = "CMS"
    re = "#bf{%s fb^{-1} (13 TeV)}" %(_lumi)
    Text.AddText(re)
    Text.SetFillStyle(0)
    Text.SetLineStyle(0)
    Text.SetBorderSize(0)
    Text.Draw()

    Text2 = TPaveText(0.25, 0.88,0.15,0.95,'brNDC')
    Text2.SetTextAlign(31);
    Text2.SetTextSize(0.065)
    Text2.AddText(leftText)
    Text2.SetFillStyle(0)
    Text2.SetLineStyle(0)
    Text2.SetBorderSize(0)
    Text2.Draw()

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.04)
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    if not 'jet' in obsName or doubleDiff:
        if 'pTj1' in obsName:
            latex2.DrawLatex(0.55,0.65, str(obs_bins[_bin])+' < '+label+' < '+str(obs_bins[_bin+1]))
        elif 'mass4l' in obsName:
            latex2.DrawLatex(0.55,0.65, str(obs_bins[0])+' < '+label+' < '+str(obs_bins[1]))
        elif doubleDiff:
            latex2.DrawLatex(0.55,0.65, str(obs_bins[_bin][0])+' < '+label+' < '+str(obs_bins[_bin][1]))
            latex2.DrawLatex(0.55,0.60, str(obs_bins[_bin][2])+' < '+label_2nd+' < '+str(obs_bins[_bin][3]))
        else:
            latex2.DrawLatex(0.45,0.65, str(obs_bins[_bin])+' < '+label+' < '+str(obs_bins[_bin+1]))
    else:
        latex2.DrawLatex(0.45,0.65, str(_bin)+' jet(s)')
    latex2.DrawLatex(0.91,0.22, "#scale[0.7]{#color[12]{68% CL}}")
    latex2.DrawLatex(0.91,0.52, "#scale[0.7]{#color[12]{95% CL}}")

    graphs[0].GetYaxis().SetRangeUser(0,maxY)
    graphs[0].GetXaxis().SetRangeUser(mini,maxi)
    c.Modified()
    graphs[0].GetYaxis().SetLimits(0,maxY)
    graphs[0].GetXaxis().SetLimits(mini,maxi)
    c.Modified()
    c.Update()

    if(opt.UNBLIND):
        if obsName!='mass4l':
            resultsXS_data['SM_125_'+obsName+'_genbin'+str(i)] = {"uncerDn": -1.0*abs(obs_nom[2]), "uncerUp": obs_nom[1], "central": obs_nom[0]}
            resultsXS_data['SM_125_'+obsName+'_genbin'+str(i)+'_statOnly'] = {"uncerDn": -1.0*abs(obs_nom_stat[2]), "uncerUp": obs_nom_stat[1], "central": obs_nom[0]}
        else:
            if _bin==0:
                resultsXS_data['SM_125_'+obsName+'_genbin'+str(i)] = {"uncerDn": -1.0*abs(obs_nom[2]), "uncerUp": obs_nom[1], "central": obs_nom[0]}
                resultsXS_data['SM_125_'+obsName+'_genbin'+str(i)+'_statOnly'] = {"uncerDn": -1.0*abs(obs_nom_stat[2]), "uncerUp": obs_nom_stat[1], "central": obs_nom[0]}
            elif _bin==1:
                resultsXS_data_v2['SM_125_'+obsName+'_2e2mu_genbin0'] = {"uncerDn": -1.0*abs(obs_nom[2]), "uncerUp": obs_nom[1], "central": obs_nom[0]}
                resultsXS_data_v2['SM_125_'+obsName+'_2e2mu_genbin0_statOnly'] = {"uncerDn": -1.0*abs(obs_nom_stat[2]), "uncerUp": obs_nom_stat[1], "central": obs_nom[0]}
            elif _bin==2:
                resultsXS_data_v2['SM_125_'+obsName+'_4mu_genbin0'] = {"uncerDn": -1.0*abs(obs_nom[2]), "uncerUp": obs_nom[1], "central": obs_nom[0]}
                resultsXS_data_v2['SM_125_'+obsName+'_4mu_genbin0_statOnly'] = {"uncerDn": -1.0*abs(obs_nom_stat[2]), "uncerUp": obs_nom_stat[1], "central": obs_nom[0]}
            elif _bin==3:
                resultsXS_data_v2['SM_125_'+obsName+'_4e_genbin0'] = {"uncerDn": -1.0*abs(obs_nom[2]), "uncerUp": obs_nom[1], "central": obs_nom[0]}
                resultsXS_data_v2['SM_125_'+obsName+'_4e_genbin0_statOnly'] = {"uncerDn": -1.0*abs(obs_nom_stat[2]), "uncerUp": obs_nom_stat[1], "central": obs_nom[0]}

    if obsName!='mass4l':
        resultsXS_asimov['SM_125_'+obsName+'_genbin'+str(i)] = {"uncerDn": -1.0*abs(exp_nom[2]), "uncerUp": exp_nom[1], "central": exp_nom[0]}
        resultsXS_asimov['SM_125_'+obsName+'_genbin'+str(i)+'_statOnly'] = {"uncerDn": -1.0*abs(exp_nom_stat[2]), "uncerUp": exp_nom_stat[1], "central": exp_nom[0]}
    else:
        if _bin==0:
            resultsXS_asimov['SM_125_'+obsName+'_genbin'+str(i)] = {"uncerDn": -1.0*abs(exp_nom[2]), "uncerUp": exp_nom[1], "central": exp_nom[0]}
            resultsXS_asimov['SM_125_'+obsName+'_genbin'+str(i)+'_statOnly'] = {"uncerDn": -1.0*abs(exp_nom_stat[2]), "uncerUp": exp_nom_stat[1], "central": exp_nom[0]}
        elif _bin==1:
            resultsXS_asimov_v2['SM_125_'+obsName+'_2e2mu_genbin0'] = {"uncerDn": -1.0*abs(exp_nom[2]), "uncerUp": exp_nom[1], "central": exp_nom[0]}
            resultsXS_asimov_v2['SM_125_'+obsName+'_2e2mu_genbin0_statOnly'] = {"uncerDn": -1.0*abs(exp_nom_stat[2]), "uncerUp": exp_nom_stat[1], "central": exp_nom[0]}
        elif _bin==2:
            resultsXS_asimov_v2['SM_125_'+obsName+'_4mu_genbin0'] = {"uncerDn": -1.0*abs(exp_nom[2]), "uncerUp": exp_nom[1], "central": exp_nom[0]}
            resultsXS_asimov_v2['SM_125_'+obsName+'_4mu_genbin0_statOnly'] = {"uncerDn": -1.0*abs(exp_nom_stat[2]), "uncerUp": exp_nom_stat[1], "central": exp_nom[0]}
        elif _bin==3:
            resultsXS_asimov_v2['SM_125_'+obsName+'_4e_genbin0'] = {"uncerDn": -1.0*abs(exp_nom[2]), "uncerUp": exp_nom[1], "central": exp_nom[0]}
            resultsXS_asimov_v2['SM_125_'+obsName+'_4e_genbin0_statOnly'] = {"uncerDn": -1.0*abs(exp_nom_stat[2]), "uncerUp": exp_nom_stat[1], "central": exp_nom[0]}

    c.Update()
    c.SaveAs("plots/lhscan_compare_"+obsName+"_"+poi+".pdf")
    c.SaveAs("plots/lhscan_compare_"+obsName+"_"+poi+".png")

with open('resultsXS_LHScan_expected_'+obsName+'_v3.py', 'w') as f:
    f.write('resultsXS = '+str(resultsXS_asimov)+' \n')
if obsName=='mass4l':
    with open('resultsXS_LHScan_expected_'+obsName+'_v2.py', 'w') as f:
        f.write('resultsXS = '+str(resultsXS_asimov_v2)+' \n')

if(opt.UNBLIND):
    with open('resultsXS_LHScan_observed_'+obsName+'_v3.py', 'w') as f:
        f.write('resultsXS = '+str(resultsXS_data)+' \n')
    if obsName=='mass4l':
        with open('resultsXS_LHScan_observed_'+obsName+'_v2.py', 'w') as f:
            f.write('resultsXS = '+str(resultsXS_data_v2)+' \n')
# raw_input()
