import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *
from sample_shortnames import *

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
    parser.add_option('',   '--lumiscale', type='string', dest='LUMISCALE', default='1.0', help='Scale yields')
    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

from ROOT import *
from tdrStyle import *
setTDRStyle()

observables = [opt.OBSNAME]

resultsXS ={}

for obsName in observables:

    if obsName=="mass4l": obsbins = ['SigmaBin0','r2e2muBin0','r4muBin0','r4eBin0']
    elif obsName=="pT4l": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5','SigmaBin6']
    elif obsName=="massZ2": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5']
    elif obsName=="massZ1": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5']
    elif obsName=="rapidity4l": obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4','SigmaBin5']
    else: obsbins = ['SigmaBin0','SigmaBin1','SigmaBin2','SigmaBin3','SigmaBin4']

    for obsbin in obsbins:

        if (obsName=="cosTheta1" and obsbin=="0"): continue


        if (opt.UNBLIND): f = TFile("higgsCombine"+obsName+"_"+obsbin+"_data.MultiDimFit.mH125.root","READ")
        elif(not opt.UNBLIND): f = TFile("higgsCombine"+obsName+"_"+obsbin+".MultiDimFit.mH125.123456.root","READ")
        if (f==0): continue

        limit = f.Get("limit")
        npoints = limit.GetEntries()

        sigma = []
        deltanll = []
        bestfit = 9999.0

        for point in range(0,npoints):
            limit.GetEntry(point)
            if (obsbin=="SigmaBin0"):
                if (point==0): bestfit=limit.SigmaBin0
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin0)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="SigmaBin1"):
                if (point==0): bestfit=limit.SigmaBin1
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin1)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="SigmaBin2"):
                if (point==0): bestfit=limit.SigmaBin2
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin2)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="SigmaBin3"):
                if (point==0): bestfit=limit.SigmaBin3
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin3)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="SigmaBin4"):
                if (point==0): bestfit=limit.SigmaBin4
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin4)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="SigmaBin5"):
                if (point==0): bestfit=limit.SigmaBin5
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin5)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="SigmaBin6"):
                if (point==0): bestfit=limit.SigmaBin6
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.SigmaBin6)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="r2e2muBin0"):
                if (point==0): bestfit=limit.r2e2muBin0
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.r2e2muBin0)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="r4muBin0"):
                if (point==0): bestfit=limit.r4muBin0
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.r4muBin0)
                        deltanll.append(2.0*limit.deltaNLL)
            if (obsbin=="r4eBin0"):
                if (point==0): bestfit=limit.r4eBin0
                if (point>0):
                    if (limit.deltaNLL<2.5):
                        sigma.append(limit.r4eBin0)
                        deltanll.append(2.0*limit.deltaNLL)

            if point>0 and len(deltanll)>0:
                if deltanll[len(deltanll)-1]>5.0 and sigma[len(sigma)-1]>bestfit: break

        if (opt.UNBLIND): fstat = TFile("higgsCombine"+obsName+"_"+obsbin+"data.MultiDimFit.mH125.root","READ")
        elif(not opt.UNBLIND): fstat = TFile("higgsCombine"+obsName+"_"+obsbin+".MultiDimFit.mH125.123456.root","READ")
        if (fstat==0): continue

        limitstat = fstat.Get("limit")
        npointsstat = limitstat.GetEntries()

        sigmastat = []
        deltanllstat = []
        bestfitstat = 9999.0

        for point in range(0,npointsstat):
            limitstat.GetEntry(point)
            if (obsbin=="SigmaBin0"):
                if (point==0): bestfitstat=limitstat.SigmaBin0
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin0)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="SigmaBin1"):
                if (point==0): bestfitstat=limitstat.SigmaBin1
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin1)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="SigmaBin2"):
                if (point==0): bestfitstat=limitstat.SigmaBin2
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin2)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="SigmaBin3"):
                if (point==0): bestfitstat=limitstat.SigmaBin3
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin3)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="SigmaBin4"):
                if (point==0): bestfitstat=limitstat.SigmaBin4
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin4)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="SigmaBin5"):
                if (point==0): bestfitstat=limitstat.SigmaBin5
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin5)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="SigmaBin6"):
                if (point==0): bestfitstat=limitstat.SigmaBin6
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.SigmaBin6)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="r2e2muBin0"):
                if (point==0): bestfitstat=limitstat.r2e2muBin0
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.r2e2muBin0)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="r4muBin0"):
                if (point==0): bestfitstat=limitstat.r4muBin0
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.r4muBin0)
                        deltanllstat.append(2.0*limitstat.deltaNLL)
            if (obsbin=="r4eBin0"):
                if (point==0): bestfitstat=limitstat.r4eBin0
                if (point>0):
                    if (limitstat.deltaNLL<2.5):
                        sigmastat.append(limitstat.r4eBin0)
                        deltanllstat.append(2.0*limitstat.deltaNLL)

            if point>0 and len(deltanllstat)>0:
                if deltanllstat[len(deltanllstat)-1]>5.0 and sigmastat[len(sigmastat)-1]>bestfitstat: break


        print obsName,obsbin
        scan = TGraph(len(sigma),array('d',sigma),array('d',deltanll))
        scanstat = TGraph(len(sigmastat),array('d',sigmastat),array('d',deltanllstat))


        c = TCanvas("c","c",1000,800)

        dummy = TH1D("dummy","dummy",1,0.0,sigmastat[len(sigmastat)-1])
        dummy.SetMinimum(0.0)
        dummy.SetMaximum(5.0)
        dummy.SetLineColor(0)
        dummy.SetMarkerColor(0)
        dummy.SetLineWidth(0)
        dummy.SetMarkerSize(0)
        dummy.GetYaxis().SetTitle("-2 #Delta lnL")
        dummy.GetXaxis().SetTitle("#sigma_{fid.} [fb]")
        dummy.Draw()

        scan.SetLineWidth(2)
        scan.SetLineColor(2)
        scan.Draw("Lsame")

        scanstat.SetLineWidth(1)
        # scanstat.SetLineStyle(2)
        scanstat.SetLineColor(4)
        scanstat.Draw("Lsame")


        gStyle.SetOptFit(0)

        f1 = TF1("f1","pol8",0.0,5.0)
        f1.SetLineColor(2)
        f1.SetLineWidth(2)
        scan.Fit("f1","N")
        # f1.Draw("same")

        f1stat = TF1("f1stat","pol8",0.0,5.0)
        f1stat.SetLineColor(4)
        f1stat.SetLineWidth(1)
        # f1stat.SetLineStyle(3)
        scanstat.Fit("f1stat","N")
        # f1stat.Draw("same")

        cl68 = TF1("cl68","1.0",0.0,5.0)
        cl68.SetLineStyle(2)
        cl68.SetLineColor(1)
        cl68.Draw("same")

        cl95 = TF1("cl95","3.84",0.0,5.0)
        cl95.SetLineStyle(2)
        cl95.SetLineColor(1)
        cl95.Draw("same")

        cl68up = 0.0
        cl68dn = 0.0
        cl95up = 0.0
        cl95dn = 0.0

        cl68upstat = 0.0
        cl68dnstat = 0.0
        cl95upstat = 0.0
        cl95dnstat = 0.0

        for i in range(0,100000):
            x = 0.+i/20000.
            #scanval = f1.Eval(x)
            scanval = scan.Eval(x)
            #if abs(scanval-3.84)<.001: print x,scanval
            if abs(scanval-1.0)<.003 and x<bestfit:
                cl68dn = round((bestfit-x),6)
            if abs(scanval-1.0)<.003 and x>bestfit:
                cl68up = round((x-bestfit),6)
            if abs(scanval-3.84)<.003 and x<bestfit:
                cl95dn = round((bestfit-x),6)
            if abs(scanval-3.84)<.003 and x>bestfit:
                cl95up = round((x-bestfit),6)

        if (cl68dn==0.0): cl68dn=round(bestfit,6)
        if (cl95dn==0.0): cl95dn=round(bestfit,6)

        print obsName,obsbin,round(bestfit,6),"+",cl68up,"-",cl68dn,"(68%)","+",cl95up,"-",cl95dn,"(95%)"

        for i in range(0,100000):
            x = 0.+i/20000.
            #scanval = f1stat.Eval(x)
            scanval = scanstat.Eval(x)
            #if abs(scanval-3.84)<.001: print x,scanval
            if abs(scanval-1.0)<.003 and x<bestfitstat:
                cl68dnstat = round((bestfitstat-x),6)
            if abs(scanval-1.0)<.003 and x>bestfitstat:
                cl68upstat = round((x-bestfitstat),6)
            if abs(scanval-3.84)<.003 and x<bestfitstat:
                cl95dnstat = round((bestfitstat-x),6)
            if abs(scanval-3.84)<.003 and x>bestfitstat:
                cl95upstat = round((x-bestfitstat),6)

        if (cl68dnstat==0.0): cl68dnstat=round(bestfitstat,6)
        if (cl95dnstat==0.0): cl95dnstat=round(bestfitstat,6)


        print obsName,obsbin,round(bestfitstat,6),"+",cl68upstat,"-",cl68dnstat,"(68%)","+",cl95upstat,"-",cl95dnstat,"(95%)"

        sysup = round(sqrt(max(0.0,cl68up**2-cl68upstat**2)),6)
        sysdn = round(sqrt(max(0.0,cl68dn**2-cl68dnstat**2)),6)
        print obsName,obsbin,round(bestfit,6),"+",cl68upstat,"-",cl68dnstat,"(stat.)","+",sysup,"-",sysdn,"(sys.)"

        latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.5*c.GetTopMargin())
        latex2.SetTextFont(42)
        latex2.SetTextAlign(31) # align right
        print opt.LUMISCALE
        if (not opt.LUMISCALE=="1.0"):
            lumi = round(137.1*float(opt.LUMISCALE),1)
            latex2.DrawLatex(0.87, 0.94,str(lumi)+" fb^{-1} (13 TeV)")
        else:
            latex2.DrawLatex(0.95, 0.95,"137.1 fb^{-1} (13 TeV)")
        latex2.SetTextSize(0.8*c.GetTopMargin())
        latex2.SetTextFont(62)
        latex2.SetTextAlign(11) # align right
        # latex2.DrawLatex(0.19, 0.95, "CMS")
        latex2.SetTextSize(0.6*c.GetTopMargin())
        latex2.SetTextFont(52)
        latex2.SetTextAlign(11)
        # latex2.DrawLatex(0.30, 0.95, "Preliminary")
        latex2.SetTextFont(42)
        latex2.SetTextSize(0.45*c.GetTopMargin())
        latex2.DrawLatex(0.30,0.85, obsName+" Bin"+obsbin)
        latex2.DrawLatex(0.30,0.78, "#sigma_{fid.} = "+str(round(bestfit,3))+" ^{+"+str(cl68up)+"}_{-"+str(cl68dn)+"} (68%) ^{+"+str(cl95up)+"}_{-"+str(cl95dn)+"} (95%)")
        latex2.DrawLatex(0.37,0.78, "#sigma_{fid.} = "+str(round(bestfit,3))+" ^{+"+str(cl68upstat)+"}_{-"+str(cl68dnstat)+"} (stat.) ^{+"+str(sysup)+"}_{-"+str(sysdn)+"} (sys.)")


        if (obsName=="mass4l"):
            if (obsbin=="SigmaBin0"):
                resultsXS['SM_125_mass4l_genbin0'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_mass4l_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            else:
                resultsXS['SM_125_mass4l_'+obsbin.replace('r','').replace('Bin0','')+'_genbin0'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_mass4l_'+obsbin.replace('r','').replace('Bin0','')+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
#        elif obsName=="pT4l":
        elif obsName=="pT4l":
            if (obsbin=="SigmaBin0"):
                resultsXS['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin1"):
                resultsXS['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin2"):
                resultsXS['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin3"):
                resultsXS['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin4"):
                resultsXS['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin5"):
                resultsXS['SM_125_'+obsName+'_genbin5'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin5_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin6"):
                resultsXS['SM_125_'+obsName+'_genbin6'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin6_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
        elif obsName=="rapidity4l" or obsName=="massZ2" or obsName=="massZ1":
            if (obsbin=="SigmaBin0"):
                resultsXS['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin1"):
                resultsXS['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin2"):
                resultsXS['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin3"):
                resultsXS['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin4"):
                resultsXS['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin5"):
                resultsXS['SM_125_'+obsName+'_genbin5'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin5_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
        else:
            if (obsbin=="SigmaBin0"):
                resultsXS['SM_125_'+obsName+'_genbin0'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin0_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin1"):
                resultsXS['SM_125_'+obsName+'_genbin1'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin1_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin2"):
                resultsXS['SM_125_'+obsName+'_genbin2'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin2_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin3"):
                resultsXS['SM_125_'+obsName+'_genbin3'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin3_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}
            if (obsbin=="SigmaBin4"):
                resultsXS['SM_125_'+obsName+'_genbin4'] = {"uncerDn": -1.0*cl68dn, "uncerUp": cl68up, "central": bestfit}
                resultsXS['SM_125_'+obsName+'_genbin4_statOnly'] = {"uncerDn": -1.0*cl68dnstat, "uncerUp": cl68upstat, "central": bestfit}


        if (opt.UNBLIND): subdir = 'data'
        elif(not opt.UNBLIND): subdir = 'asimov'
        c.SaveAs("plots/"+obsName+"/"+subdir+"/lhscan_"+obsName+"_"+obsbin+".pdf")
        c.SaveAs("plots/"+obsName+"/"+subdir+"/lhscan_"+obsName+"_"+obsbin+".png")


        if (obsName=="mass4l"):
            if (obsbin=="SigmaBin0"):
                with open('resultsXS_LHScan_mass4l_v3.py', 'w') as f:
                    f.write('resultsXS = '+str(resultsXS)+' \n')
            else:
                with open('resultsXS_LHScan_mass4l_v2.py', 'w') as f:
                    f.write('resultsXS = '+str(resultsXS)+' \n')
        else:
            with open('resultsXS_LHScan_'+obsName+'_v3.py', 'w') as f:
                f.write('resultsXS = '+str(resultsXS)+' \n')
