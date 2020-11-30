import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *

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
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

from ROOT import *
from tdrStyle_lhs_compare import *
setTDRStyle()

resultsXS_asimov ={}
resultsXS_data ={}

obsName = opt.OBSNAME
if(obsName == 'rapidity4l'): label = '|y_{H}|'
if(obsName == 'pT4l'): label = 'p$_T^H$ (GeV)'
if(obsName == 'massZ1'): label = 'm$_{Z1}$ (GeV)'
if(obsName == 'massZ2'): label = 'm$_{Z2}$ (GeV)'
if(obsName == 'njets_pt30_eta2p5'): label = 'N$_{jet}$'
if(obsName == 'pTj1'): label = 'p$_T^{j1}$ (GeV)'
if(obsName == 'mass4l'): label = 'mass4l (GeV)'

bins = {'rapidity4l': [0, 0.15, 0.3, 0.6, 0.9, 1.2, 2.5], 'pT4l': [0, 10, 20, 30, 45, 80, 120, 200, 1300], 'njets_pt30_eta2p5': [0,1,2,3,4,20], 'pTj1': [0,30,55,95,200,1300], 'mass4l': [105,140]}
if obsName=="mass4l": obsbins = ['SigmaBin0','r2e2muBin0','r4muBin0','r4eBin0']
elif obsName=="pT4l": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5','SigmaBin6', 'SigmaBin7']
elif obsName=="massZ2": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5']
elif obsName=="massZ1": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5']
elif obsName=="rapidity4l": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5']
elif obsName=="njets_pt30_eta2p5": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4']
elif obsName=="pTj1": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4']
else: obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4']

for obsbin in obsbins:

    if (obsName=="cosTheta1" and obsbin=="0"): continue

    if(opt.UNBLIND):
        f_data = TFile("/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/combine_files/higgsCombine_"+obsName+"_"+obsbin+".MultiDimFit.mH125.38.root","READ")
        if (f_data==0): continue
    f_asimov = TFile("/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/combine_files/higgsCombine_"+obsName+"_"+obsbin+".MultiDimFit.mH125.38.123456.root","READ")
    if (f_asimov==0): continue

    if(opt.UNBLIND):
        limit_data = f_data.Get("limit")
        npoints_data = limit_data.GetEntries()
    limit_asimov = f_asimov.Get("limit")
    npoints_asimov = limit_asimov.GetEntries()

    if(opt.UNBLIND):
        sigma_data = []
        deltanll_data = []
        bestfit_data = 9999.0
    sigma_asimov = []
    deltanll_asimov = []
    bestfit_asimov = 9999.0
    index_print = 0 # It will be used when printing the bin boundaries.

    if(opt.UNBLIND):
        for point in range(0,npoints_data):
            limit_data.GetEntry(point)
            if (obsbin=="SigmaBin0"):
                index_print = 0
                if (point==0): bestfit_data=limit_data.SigmaBin0
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin0)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin1"):
                index_print = 1
                if (point==0): bestfit_data=limit_data.SigmaBin1
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin1)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin2"):
                index_print = 2
                if (point==0): bestfit_data=limit_data.SigmaBin2
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin2)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin3"):
                index_print = 3
                if (point==0): bestfit_data=limit_data.SigmaBin3
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin3)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin4"):
                index_print = 4
                if (point==0): bestfit_data=limit_data.SigmaBin4
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin4)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin5"):
                index_print = 5
                if (point==0): bestfit_data=limit_data.SigmaBin5
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin5)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin6"):
                index_print = 6
                if (point==0): bestfit_data=limit_data.SigmaBin6
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin6)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="SigmaBin7"):
                index_print = 7
                if (point==0): bestfit_data=limit_data.SigmaBin7
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.SigmaBin7)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="r2e2muBin0"):
                if (point==0): bestfit_data=limit_data.r2e2muBin0
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.r2e2muBin0)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="r4muBin0"):
                if (point==0): bestfit_data=limit_data.r4muBin0
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.r4muBin0)
                        deltanll_data.append(2.0*limit_data.deltaNLL)
            if (obsbin=="r4eBin0"):
                if (point==0): bestfit_data=limit_data.r4eBin0
                if (point>0):
                    if (limit_data.deltaNLL<2.5):
                        sigma_data.append(limit_data.r4eBin0)
                        deltanll_data.append(2.0*limit_data.deltaNLL)

            if point>0 and len(deltanll_data)>0:
                if deltanll_data[len(deltanll_data)-1]>5.0 and sigma_data[len(sigma_data)-1]>bestfit_data: break

    for point in range(0,npoints_asimov):
        limit_asimov.GetEntry(point)
        if (obsbin=="SigmaBin0"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin0
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin0)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin1"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin1
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin1)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin2"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin2
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin2)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin3"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin3
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin3)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin4"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin4
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin4)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin5"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin5
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin5)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin6"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin6
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin6)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="SigmaBin7"):
            if (point==0): bestfit_asimov=limit_asimov.SigmaBin7
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.SigmaBin7)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="r2e2muBin0"):
            if (point==0): bestfit_asimov=limit_asimov.r2e2muBin0
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.r2e2muBin0)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="r4muBin0"):
            if (point==0): bestfit_asimov=limit_asimov.r4muBin0
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.r4muBin0)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)
        if (obsbin=="r4eBin0"):
            if (point==0): bestfit_asimov=limit_asimov.r4eBin0
            if (point>0):
                if (limit_asimov.deltaNLL<2.5):
                    sigma_asimov.append(limit_asimov.r4eBin0)
                    deltanll_asimov.append(2.0*limit_asimov.deltaNLL)

        if point>0 and len(deltanll_asimov)>0:
            if deltanll_asimov[len(deltanll_asimov)-1]>5.0 and sigma_asimov[len(sigma_asimov)-1]>bestfit_asimov: break

    if(opt.UNBLIND):
        fstat_data = TFile("/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/combine_files/higgsCombine_"+obsName+"_"+obsbin+"_NoSys.MultiDimFit.mH125.38.root","READ")
        if (fstat_data==0): continue
    fstat_asimov = TFile("/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXSFWK/combine_files/higgsCombine_"+obsName+"_"+obsbin+"_NoSys_exp.MultiDimFit.mH125.38.123456.root","READ")
    if (fstat_asimov==0): continue

    if(opt.UNBLIND):
        limitstat_data = fstat_data.Get("limit")
        npointsstat_data = limitstat_data.GetEntries()
    limitstat_asimov = fstat_asimov.Get("limit")
    npointsstat_asimov = limitstat_asimov.GetEntries()

    if(opt.UNBLIND):
        sigmastat_data = []
        deltanllstat_data = []
        bestfitstat_data = 9999.0
    sigmastat_asimov = []
    deltanllstat_asimov = []
    bestfitstat_asimov = 9999.0

    if(opt.UNBLIND):
        for point in range(0,npointsstat_data):
            limitstat_data.GetEntry(point)
            if (obsbin=="SigmaBin0"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin0
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin0)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin1"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin1
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin1)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin2"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin2
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin2)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin3"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin3
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin3)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin4"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin4
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin4)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin5"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin5
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin5)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin6"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin6
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin6)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="SigmaBin7"):
                if (point==0): bestfitstat_data=limitstat_data.SigmaBin7
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.SigmaBin7)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="r2e2muBin0"):
                if (point==0): bestfitstat_data=limitstat_data.r2e2muBin0
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.r2e2muBin0)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="r4muBin0"):
                if (point==0): bestfitstat_data=limitstat_data.r4muBin0
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.r4muBin0)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)
            if (obsbin=="r4eBin0"):
                if (point==0): bestfitstat_data=limitstat_data.r4eBin0
                if (point>0):
                    if (limitstat_data.deltaNLL<2.5):
                        sigmastat_data.append(limitstat_data.r4eBin0)
                        deltanllstat_data.append(2.0*limitstat_data.deltaNLL)

            if point>0 and len(deltanllstat_data)>0:
                if deltanllstat_data[len(deltanllstat_data)-1]>5.0 and sigmastat_data[len(sigmastat_data)-1]>bestfitstat_data: break


    for point in range(0,npointsstat_asimov):
        limitstat_asimov.GetEntry(point)
        if (obsbin=="SigmaBin0"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin0
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin0)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin1"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin1
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin1)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin2"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin2
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin2)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin3"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin3
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin3)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin4"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin4
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin4)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin5"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin5
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin5)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin6"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin6
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin6)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="SigmaBin7"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.SigmaBin7
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.SigmaBin7)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="r2e2muBin0"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.r2e2muBin0
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.r2e2muBin0)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="r4muBin0"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.r4muBin0
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.r4muBin0)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)
        if (obsbin=="r4eBin0"):
            if (point==0): bestfitstat_asimov=limitstat_asimov.r4eBin0
            if (point>0):
                if (limitstat_asimov.deltaNLL<2.5):
                    sigmastat_asimov.append(limitstat_asimov.r4eBin0)
                    deltanllstat_asimov.append(2.0*limitstat_asimov.deltaNLL)

        if point>0 and len(deltanllstat_asimov)>0:
            if deltanllstat_asimov[len(deltanllstat_asimov)-1]>5.0 and sigmastat_asimov[len(sigmastat_asimov)-1]>bestfitstat_asimov: break


    print obsName,obsbin
    if(opt.UNBLIND):
        scan_data = TGraph(len(sigma_data),array('d',sigma_data),array('d',deltanll_data))
        scanstat_data = TGraph(len(sigmastat_data),array('d',sigmastat_data),array('d',deltanllstat_data))
    scan_asimov = TGraph(len(sigma_asimov),array('d',sigma_asimov),array('d',deltanll_asimov))
    scanstat_asimov = TGraph(len(sigmastat_asimov),array('d',sigmastat_asimov),array('d',deltanllstat_asimov))


    c = TCanvas("c","c",1000,800)
    gStyle.SetOptTitle(0)
    c.SetHighLightColor(2)
    c.SetBorderMode(0)
    c.SetBorderSize(2)
    c.SetLeftMargin(0.1454155)
    c.SetRightMargin(0.07378224)
    c.SetFrameBorderMode(0)
    c.SetFrameBorderMode(0)

    dummy = TH1D("dummy","dummy",1,0.0,sigmastat_asimov[len(sigmastat_asimov)-1]+0.2)
    dummy.SetMinimum(0.0)
    dummy.SetMaximum(8.0)
    dummy.SetLineColor(0)
    dummy.SetMarkerColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.GetYaxis().SetTitle("-2 #Delta lnL")
    dummy.GetXaxis().SetTitle("#sigma_{fid.} [fb]")
    dummy.Draw()

    mg = TMultiGraph()

    if(opt.UNBLIND):
        scan_data.SetLineWidth(2)
        scan_data.SetLineColor(1)
        # scan_data.Draw("same")
        mg.Add(scan_data)

    scan_asimov.SetLineWidth(2)
    scan_asimov.SetLineColor(2)
    # scan_asimov.Draw("same")
    mg.Add(scan_asimov)

    if(opt.UNBLIND):
        scanstat_data.SetLineWidth(2)
        scanstat_data.SetLineStyle(9)
        scanstat_data.SetLineColor(1)
        # scanstat_data.Draw("same")
        mg.Add(scanstat_data)

    scanstat_asimov.SetLineWidth(2)
    scanstat_asimov.SetLineStyle(9)
    scanstat_asimov.SetLineColor(2)
    # scanstat_asimov.Draw("same")
    mg.Add(scanstat_asimov)

    mg.Draw("sameapl")
    mg.SetMinimum(0)
    mg.SetMaximum(8)

    ## Fit
    # gStyle.SetOptFit(0)
    #
    # f1 = TF1("f1","pol2",0.0,5.0)
    # f1.SetLineColor(2)
    # f1.SetLineWidth(2)
    # scan_data.Fit("f1","N")
    # f1.Draw("same")
    #
    # f1stat = TF1("f1stat","pol8",0.0,5.0)
    # f1stat.SetLineColor(4)
    # f1stat.SetLineWidth(1)
    # # f1stat.SetLineStyle(3)
    # scanstat.Fit("f1stat","N")
    # # f1stat.Draw("same")
    #
    cl68 = TF1("cl68","1.0",0.0,5.0)
    cl68.SetLineStyle(3)
    cl68.SetLineColor(12)
    cl68.Draw("same")

    cl95 = TF1("cl95","3.84",0.0,5.0)
    cl95.SetLineStyle(3)
    cl95.SetLineColor(12)
    cl95.Draw("same")


    if(opt.UNBLIND):
        cl68up_data = 0.0
        cl68dn_data = 0.0
        cl95up_data = 0.0
        cl95dn_data = 0.0

    cl68up_asimov = 0.0
    cl68dn_asimov = 0.0
    cl95up_asimov = 0.0
    cl95dn_asimov = 0.0

    if(opt.UNBLIND):
        cl68upstat_data = 0.0
        cl68dnstat_data = 0.0
        cl95upstat_data = 0.0
        cl95dnstat_data = 0.0

    cl68upstat_asimov = 0.0
    cl68dnstat_asimov = 0.0
    cl95upstat_asimov = 0.0
    cl95dnstat_asimov = 0.0

    if(opt.UNBLIND):
        for i in range(0,100000):
            x = 0.+i/20000.
            #scanval = f1.Eval(x)
            scanval = scan_data.Eval(x)
            #if abs(scanval-3.84)<.001: print x,scanval
            if abs(scanval-1.0)<.003 and x<bestfit_data:
                cl68dn_data = round((bestfit_data-x),6)
            if abs(scanval-1.0)<.003 and x>bestfit_data:
                cl68up_data = round((x-bestfit_data),6)
            if abs(scanval-3.84)<.003 and x<bestfit_data:
                cl95dn_data = round((bestfit_data-x),6)
            if abs(scanval-3.84)<.003 and x>bestfit_data:
                cl95up_data = round((x-bestfit_data),6)

        if (cl68dn_data==0.0): cl68dn_data=round(bestfit_data,6)
        if (cl95dn_data==0.0): cl95dn_data=round(bestfit_data,6)

        print obsName,obsbin,round(bestfit_data,6),"+",cl68up_data,"-",cl68dn_data,"(68%)","+",cl95up_data,"-",cl95dn_data,"(95%)"

        for i in range(0,100000):
            x = 0.+i/20000.
            #scanval = f1stat.Eval(x)
            scanval = scanstat_data.Eval(x)
            #if abs(scanval-3.84)<.001: print x,scanval
            if abs(scanval-1.0)<.003 and x<bestfitstat_data:
                cl68dnstat_data = round((bestfitstat_data-x),6)
            if abs(scanval-1.0)<.003 and x>bestfitstat_data:
                cl68upstat_data = round((x-bestfitstat_data),6)
            if abs(scanval-3.84)<.003 and x<bestfitstat_data:
                cl95dnstat_data = round((bestfitstat_data-x),6)
            if abs(scanval-3.84)<.003 and x>bestfitstat_data:
                cl95upstat_data = round((x-bestfitstat_data),6)

        if (cl68dnstat_data==0.0): cl68dnstat_data=round(bestfitstat_data,6)
        if (cl95dnstat_data==0.0): cl95dnstat_data=round(bestfitstat_data,6)


        print obsName,obsbin,round(bestfitstat_data,6),"+",cl68upstat_data,"-",cl68dnstat_data,"(68%)","+",cl95upstat_data,"-",cl95dnstat_data,"(95%)"

    for i in range(0,100000):
        x = 0.+i/20000.
        #scanval = f1.Eval(x)
        scanval = scan_asimov.Eval(x)
        #if abs(scanval-3.84)<.001: print x,scanval
        if abs(scanval-1.0)<.003 and x<bestfit_asimov:
            cl68dn_asimov = round((bestfit_asimov-x),6)
        if abs(scanval-1.0)<.003 and x>bestfit_asimov:
            cl68up_asimov = round((x-bestfit_asimov),6)
        if abs(scanval-3.84)<.003 and x<bestfit_asimov:
            cl95dn_asimov = round((bestfit_asimov-x),6)
        if abs(scanval-3.84)<.003 and x>bestfit_asimov:
            cl95up_asimov = round((x-bestfit_asimov),6)

    if (cl68dn_asimov==0.0): cl68dn_asimov=round(bestfit_asimov,6)
    if (cl95dn_asimov==0.0): cl95dn_asimov=round(bestfit_asimov,6)

    print obsName,obsbin,round(bestfit_asimov,6),"+",cl68up_asimov,"-",cl68dn_asimov,"(68%)","+",cl95up_asimov,"-",cl95dn_asimov,"(95%)"

    for i in range(0,100000):
        x = 0.+i/20000.
        #scanval = f1stat.Eval(x)
        scanval = scanstat_asimov.Eval(x)
        #if abs(scanval-3.84)<.001: print x,scanval
        if abs(scanval-1.0)<.003 and x<bestfitstat_asimov:
            cl68dnstat_asimov = round((bestfitstat_asimov-x),6)
        if abs(scanval-1.0)<.003 and x>bestfitstat_asimov:
            cl68upstat_asimov = round((x-bestfitstat_asimov),6)
        if abs(scanval-3.84)<.003 and x<bestfitstat_asimov:
            cl95dnstat_asimov = round((bestfitstat_asimov-x),6)
        if abs(scanval-3.84)<.003 and x>bestfitstat_asimov:
            cl95upstat_asimov = round((x-bestfitstat_asimov),6)

    if (cl68dnstat_asimov==0.0): cl68dnstat_asimov=round(bestfitstat_asimov,6)
    if (cl95dnstat_asimov==0.0): cl95dnstat_asimov=round(bestfitstat_asimov,6)


    print obsName,obsbin,round(bestfitstat_asimov,6),"+",cl68upstat_asimov,"-",cl68dnstat_asimov,"(68%)","+",cl95upstat_asimov,"-",cl95dnstat_asimov,"(95%)"


    if(opt.UNBLIND):
        sysup_data = round(sqrt(max(0.0,cl68up_data**2-cl68upstat_data**2)),6)
        sysdn_data = round(sqrt(max(0.0,cl68dn_data**2-cl68dnstat_data**2)),6)
        print obsName,obsbin,round(bestfit_data,6),"+",cl68upstat_data,"-",cl68dnstat_data,"(stat.)","+",sysup_data,"-",sysdn_data,"(sys.)"
    sysup_asimov = round(sqrt(max(0.0,cl68up_asimov**2-cl68upstat_asimov**2)),6)
    sysdn_asimov = round(sqrt(max(0.0,cl68dn_asimov**2-cl68dnstat_asimov**2)),6)
    print obsName,obsbin,round(bestfit_asimov,6),"+",cl68upstat_asimov,"-",cl68dnstat_asimov,"(stat.)","+",sysup_asimov,"-",sysdn_asimov,"(sys.)"

    legend = TLegend(0.62,0.69,0.96,0.89)
    legend.SetBorderSize(1)
    legend.SetTextSize(0.03)
    legend.SetLineColor(0)
    legend.SetLineStyle(0)
    legend.SetLineWidth(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    if(opt.UNBLIND):
        legend.AddEntry(scan_data,"Observed", "l1")
        legend.AddEntry(scanstat_data,"Observed (stat only)","l1")
    legend.AddEntry(scan_asimov,"Expected","l1")
    legend.AddEntry(scanstat_asimov,"Expected (stat only)","l1")
    legend.Draw()


    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.04)
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    # print opt.LUMISCALE
    # if (not opt.LUMISCALE=="1.0"):
    #     lumi = round(137.1*float(opt.LUMISCALE),1)
    #     latex2.DrawLatex(0.87, 0.94,str(lumi)+" fb^{-1} (13 TeV)")
    # else:
    if(opt.YEAR=='2016'):
        latex2.DrawLatex(0.92, 0.95,"35.9 fb^{-1} (13 TeV)")
    if(opt.YEAR=='2017'):
        latex2.DrawLatex(0.92, 0.95,"41.5 fb^{-1} (13 TeV)")
    if(opt.YEAR=='2018'):
        latex2.DrawLatex(0.92, 0.95,"59.7 fb^{-1} (13 TeV)")
    if(opt.YEAR=='Full'):
        latex2.DrawLatex(0.92, 0.95,"137 fb^{-1} (13 TeV)")
    latex2.SetTextSize(0.8*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.15, 0.95, "CMS")
    latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    # latex2.DrawLatex(0.30, 0.95, "Preliminary")
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.45*c.GetTopMargin())
    # latex2.DrawLatex(0.50,0.85, obsName+" Bin"+obsbin)

    latex2.DrawLatex(0.20,0.87, str(bins[obsName][index_print])+' < '+label+' < '+str(bins[obsName][index_print+1]))
    # latex2.DrawLatex(0.30,0.78, "#sigma_{fid.} = "+str(round(bestfit,3))+" ^{+"+str(cl68up)+"}_{-"+str(cl68dn)+"} (68%) ^{+"+str(cl95up)+"}_{-"+str(cl95dn)+"} (95%)")
    if(opt.UNBLIND): latex2.DrawLatex(0.20,0.80, "#sigma_{fid.}^{obs} = "+str(round(bestfit_data,3))+" ^{+"+str(round(cl68upstat_data,3))+"}_{-"+str(round(cl68dnstat_data,3))+"} (stat.) ^{+"+str(round(sysup_data,3))+"}_{-"+str(round(sysdn_data,3))+"} (sys.)")
    latex2.DrawLatex(0.20,0.72, "#sigma_{fid.}^{exp} = "+str(round(bestfit_asimov,3))+" ^{+"+str(round(cl68upstat_asimov,3))+"}_{-"+str(round(cl68dnstat_asimov,3))+"} (stat.) ^{+"+str(round(sysup_asimov,3))+"}_{-"+str(round(sysdn_asimov,3))+"} (sys.)")
    latex2.DrawLatex(0.18,0.24, "#scale[0.7]{#color[12]{68% CL}}")
    latex2.DrawLatex(0.18,0.52, "#scale[0.7]{#color[12]{95% CL}}")


    # Expected
    if (obsName=="mass4l"):
        if (obsbin=="SigmaBin0"):
            resultsXS_asimov['SM_125_mass4l_genbin0'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_mass4l_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        else:
            resultsXS_asimov['SM_125_mass4l_'+obsbin.replace('r','').replace('Bin0','')+'_genbin0'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_mass4l_'+obsbin.replace('r','').replace('Bin0','')+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
    elif obsName=="pT4l":
        if (obsbin=="SigmaBin0"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin1"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin2"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin3"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin4"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin5"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin5'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin5_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin6"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin6'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin6_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin7"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin7'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin7_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
    elif obsName=="rapidity4l" or obsName=="massZ2" or obsName=="massZ1":
        if (obsbin=="SigmaBin0"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin1"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin2"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin3"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin4"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin5"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin5'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin5_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
    else:
        if (obsbin=="SigmaBin0"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin1"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin2"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin3"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}
        if (obsbin=="SigmaBin4"):
            resultsXS_asimov['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn_asimov, "uncerUp": cl68up_asimov, "central": bestfit_asimov}
            resultsXS_asimov['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat_asimov, "uncerUp": cl68upstat_asimov, "central": bestfit_asimov}

    if opt.UNBLIND:
        # Observed
        if (obsName=="mass4l"):
            if (obsbin=="SigmaBin0"):
                resultsXS_data['SM_125_mass4l_genbin0'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_mass4l_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            else:
                resultsXS_data['SM_125_mass4l_'+obsbin.replace('r','').replace('Bin0','')+'_genbin0'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_mass4l_'+obsbin.replace('r','').replace('Bin0','')+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
        elif obsName=="pT4l":
            if (obsbin=="SigmaBin0"):
                resultsXS_data['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin1"):
                resultsXS_data['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin2"):
                resultsXS_data['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin3"):
                resultsXS_data['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin4"):
                resultsXS_data['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin5"):
                resultsXS_data['SM_125_'+obsName+'_genbin5'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin5_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin6"):
                resultsXS_data['SM_125_'+obsName+'_genbin6'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin6_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin7"):
                resultsXS_data['SM_125_'+obsName+'_genbin7'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin7_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
        elif obsName=="rapidity4l" or obsName=="massZ2" or obsName=="massZ1":
            if (obsbin=="SigmaBin0"):
                resultsXS_data['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin1"):
                resultsXS_data['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin2"):
                resultsXS_data['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin3"):
                resultsXS_data['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin4"):
                resultsXS_data['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin5"):
                resultsXS_data['SM_125_'+obsName+'_genbin5'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin5_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
        else:
            if (obsbin=="SigmaBin0"):
                resultsXS_data['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin1"):
                resultsXS_data['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin2"):
                resultsXS_data['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin3"):
                resultsXS_data['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}
            if (obsbin=="SigmaBin4"):
                resultsXS_data['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn_data, "uncerUp": cl68up_data, "central": bestfit_data}
                resultsXS_data['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat_data, "uncerUp": cl68upstat_data, "central": bestfit_data}


    c.Update()
    c.SaveAs("plots/lhscan_compare_"+obsName+"_"+obsbin+".pdf")
    c.SaveAs("plots/lhscan_compare_"+obsName+"_"+obsbin+".png")


    if (obsName=="mass4l"):
        if (obsbin=="SigmaBin0"):
            with open('resultsXS_LHScan_expected_mass4l_v3.py', 'w') as f:
                f.write('resultsXS = '+str(resultsXS_asimov)+' \n')
            with open('resultsXS_LHScan_observed_mass4l_v3.py', 'w') as f:
                f.write('resultsXS = '+str(resultsXS_data)+' \n')
        else:
            with open('resultsXS_LHScan_expected_mass4l_v2.py', 'w') as f:
                f.write('resultsXS = '+str(resultsXS_asimov)+' \n')
            with open('resultsXS_LHScan_observed_mass4l_v2.py', 'w') as f:
                f.write('resultsXS = '+str(resultsXS_data)+' \n')
    else:
        with open('resultsXS_LHScan_expected_'+obsName+'_v3.py', 'w') as f:
            f.write('resultsXS = '+str(resultsXS_asimov)+' \n')
        with open('resultsXS_LHScan_observed_'+obsName+'_v3.py', 'w') as f:
            f.write('resultsXS = '+str(resultsXS_data)+' \n')
