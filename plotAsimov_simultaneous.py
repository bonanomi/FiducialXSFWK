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
    parser.add_option('-d', '--dir',    dest='SOURCEDIR',  type='string',default='./combine_files/', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--asimovModel',dest='ASIMOV',type='string',default='SM_125', help='Name of the asimov data mode')
    parser.add_option('',   '--asimovMass',dest='ASIMOVMASS',type='string',default='125.0', help='Asimov Mass')
    parser.add_option('',   '--unfoldModel',dest='UNFOLD',type='string',default='SM_125', help='Name of the unfolding model') 
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "inclusive", "pT", "eta", "Njets"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='Use results from fixed fraction fit, default is False')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)
    
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

def checkDir(folder_path):
    isdir = os.path.isdir(folder_path)
    if not isdir:
        print('Directory {} does not exist. Creating it.' .format(folder_path))
        os.mkdir(folder_path)
            
modelName = opt.UNFOLD
physicalModel = 'v3'
asimovDataModel = opt.ASIMOV
asimovPhysicalModel = 'v3'
obsName = opt.OBSNAME
observableBins = opt.OBSBINS.split('|')
observableBins.pop()
observableBins.pop(0)

def plotAsimov_sim(asimovDataModel, asimovPhysicalModel, modelName, physicalModel, obsName, fstate, observableBins, recobin):

    sourcedir = opt.SOURCEDIR
    theorymass = opt.THEORYMASS
    year = opt.YEAR

    if year == '2016':
        lumi = '35.9'
    elif year == '2017':
        lumi = '41.5'
    elif year == '2018':
        lumi = '59.7'
    else:
        lumi = '137'

    nBins = len(observableBins)-1
    channel = {"4mu":"1", "4e":"2", "2e2mu":"3", "4l":"2"} # 4l is dummy, won't be used

    # Load some libraries                                 
    ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
    ROOT.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so") #Print 0 in case of succesfull loading
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
    ROOT.gSystem.AddIncludePath("-Iinclude/")
    
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

    if(not opt.UNBLIND):
    	theorymass = theorymass + '.123456'
    fname = 'higgsCombine_'+obsName+'_SigmaBin0.MultiDimFit.mH'+theorymass+'.root'

    f_asimov = TFile(sourcedir + fname, "READ")

    if (not opt.UNBLIND):
        data = f_asimov.Get("toys/toy_asimov");
    w_asimov = f_asimov.Get("w")
    if (opt.UNBLIND):
        data = w_asimov.data("data_obs")
    w_asimov.loadSnapshot("clean")

    trueH_asimov = {}
    zjets_asimov = {}
    ggzz_asimov = {}
    fakeH_asimov = {}
    out_trueH_asimov = {}
    qqzz_asimov = {}
    n_trueH_asimov = {}
    n_trueH_asimov["4l"] = 0.0
    n_trueH_otherfid_asimov = {}
    n_trueH_otherfid_asimov["4l"] = 0.0
    n_zjets_asimov = {}
    n_zjets_asimov["4l"] = 0.0
    n_ggzz_asimov = {}
    n_ggzz_asimov["4l"] = 0.0
    n_fakeH_asimov = {}
    n_fakeH_asimov["4l"] = 0.0
    n_out_trueH_asimov = {}
    n_out_trueH_asimov["4l"] = 0.0
    n_qqzz_asimov = {}
    n_qqzz_asimov["4l"] = 0.0
    n_zz_asimov = {}
    n_zz_asimov["4l"] = 0.0
    
                                                                        
    fStates = ['4mu','4e','2e2mu']    
    for fState in fStates:
        for bin in range(nBins):
            trueH_asimov[fState+"Bin"+str(bin)] = w_asimov.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_trueH"+fState+"Bin"+str(bin))            
            print fState+"Bin"+str(bin) 
            print "n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_trueH"+fState+"Bin"+str(bin)
            print trueH_asimov[fState+"Bin"+str(bin)].getVal()

        zjets_asimov[fState] = w_asimov.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_bkg_zjets")
        ggzz_asimov[fState] = w_asimov.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_bkg_ggzz")
        fakeH_asimov[fState] = w_asimov.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_fakeH")
        out_trueH_asimov[fState] = w_asimov.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_out_trueH")
        qqzz_asimov[fState] = w_asimov.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_bkg_qqzz")
        n_trueH_otherfid_asimov[fState] = 0.0
        for bin in range(nBins):
            if (bin==recobin): n_trueH_asimov[fState] = trueH_asimov[fState+"Bin"+str(bin)].getVal()
            else: n_trueH_otherfid_asimov[fState] += trueH_asimov[fState+"Bin"+str(bin)].getVal()
        n_zjets_asimov[fState] = zjets_asimov[fState].getVal()
        n_ggzz_asimov[fState] = ggzz_asimov[fState].getVal() 
        n_fakeH_asimov[fState] = fakeH_asimov[fState].getVal()
        n_out_trueH_asimov[fState] = out_trueH_asimov[fState].getVal()
        n_qqzz_asimov[fState] = qqzz_asimov[fState].getVal()
        n_zz_asimov[fState] = n_ggzz_asimov[fState]+n_qqzz_asimov[fState]        
        n_trueH_asimov["4l"] += n_trueH_asimov[fState]
        n_trueH_otherfid_asimov["4l"] += n_trueH_otherfid_asimov[fState]
        n_zjets_asimov["4l"] += zjets_asimov[fState].getVal()
        n_ggzz_asimov["4l"] += ggzz_asimov[fState].getVal()
        n_fakeH_asimov["4l"] += fakeH_asimov[fState].getVal()
        n_out_trueH_asimov["4l"] += out_trueH_asimov[fState].getVal()
        n_qqzz_asimov["4l"] += qqzz_asimov[fState].getVal()
        n_zz_asimov["4l"] += n_ggzz_asimov[fState]+n_qqzz_asimov[fState]                                                                
        

    f_modelfit = TFile(sourcedir + fname, "READ")
    w_modelfit = f_modelfit.Get("w")    
    sim = w_modelfit.pdf("model_s")
    #sim.Print("v")
    if (fstate=="4l"): pdfi = sim.getPdf("ch1_ch1") # dummy won't be used
    else: pdfi = sim.getPdf("ch"+channel[fstate]+"_ch"+str(recobin+1))     
    CMS_zz4l_mass = w_modelfit.var("CMS_zz4l_mass")
    w_modelfit.loadSnapshot("MultiDimFit")
    #pdfi.Print("v")

    trueH_modelfit = {}
    zjets_modelfit = {}
    ggzz_modelfit = {}
    fakeH_modelfit = {}
    out_trueH_modelfit = {}
    qqzz_modelfit = {}
    n_trueH_modelfit = {}
    n_trueH_modelfit["4l"] = 0.0
    n_trueH_otherfid_modelfit = {}
    n_trueH_otherfid_modelfit["4l"] = 0.0
    n_zjets_modelfit = {}
    n_zjets_modelfit["4l"] = 0.0
    n_ggzz_modelfit = {}
    n_ggzz_modelfit["4l"] = 0.0
    n_fakeH_modelfit = {}
    n_fakeH_modelfit["4l"] = 0.0
    n_out_trueH_modelfit = {}
    n_out_trueH_modelfit["4l"] = 0.0
    n_qqzz_modelfit = {}
    n_qqzz_modelfit["4l"] = 0.0
    n_zz_modelfit = {}
    n_zz_modelfit["4l"] = 0.0
        
    for fState in fStates:
        for bin in range(nBins):
            trueH_modelfit[fState+"Bin"+str(bin)] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_trueH"+fState+"Bin"+str(bin))
            trueH_modelfit[fState+"Bin"+str(bin)] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_trueH"+fState+"Bin"+str(bin))            
            print fState+"Bin"+str(bin) 
            print "n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_trueH"+fState+"Bin"+str(bin)
            print trueH_modelfit[fState+"Bin"+str(bin)].getVal()

        zjets_modelfit[fState] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_bkg_zjets")
        ggzz_modelfit[fState] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_bkg_ggzz")
        fakeH_modelfit[fState] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_fakeH")
        out_trueH_modelfit[fState] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_out_trueH")
        qqzz_modelfit[fState] = w_modelfit.function("n_exp_final_binch"+channel[fState]+"_ch"+str(recobin+1)+"_proc_bkg_qqzz")
        n_trueH_otherfid_modelfit[fState] = 0.0
        for bin in range(nBins):
            if (bin==recobin): n_trueH_modelfit[fState] = trueH_modelfit[fState+"Bin"+str(bin)].getVal()
            else: n_trueH_otherfid_modelfit[fState] += trueH_modelfit[fState+"Bin"+str(bin)].getVal()
        n_zjets_modelfit[fState] = zjets_modelfit[fState].getVal()
        n_ggzz_modelfit[fState] = ggzz_modelfit[fState].getVal()
        n_fakeH_modelfit[fState] = fakeH_modelfit[fState].getVal()
        n_out_trueH_modelfit[fState] = out_trueH_modelfit[fState].getVal()
        n_qqzz_modelfit[fState] = qqzz_modelfit[fState].getVal()
        n_zz_modelfit[fState] = n_ggzz_modelfit[fState]+n_qqzz_modelfit[fState]
        n_trueH_modelfit["4l"] += n_trueH_modelfit[fState]
        n_trueH_otherfid_modelfit["4l"] += n_trueH_otherfid_modelfit[fState]
        n_zjets_modelfit["4l"] += zjets_modelfit[fState].getVal()
        n_ggzz_modelfit["4l"] += ggzz_modelfit[fState].getVal()
        n_fakeH_modelfit["4l"] += fakeH_modelfit[fState].getVal()
        n_out_trueH_modelfit["4l"] += out_trueH_modelfit[fState].getVal()
        n_qqzz_modelfit["4l"] += qqzz_modelfit[fState].getVal()
        n_zz_modelfit["4l"] += n_ggzz_modelfit[fState]+n_qqzz_modelfit[fState]


    CMS_channel = w.cat("CMS_channel")
    mass = w.var("CMS_zz4l_mass").frame(RooFit.Bins(15))
    #mass = w.var("CMS_zz4l_mass").frame(RooFit.Bins(45))

    if (fstate=="4l"):

        datacut = ''
        for fState in fStates:
            datacut += "CMS_channel==CMS_channel::ch"+channel[fState]+"_ch"+str(recobin+1)+" || "
        datacut = datacut.rstrip(" || ")
        data = data.reduce(RooFit.Cut(datacut))        
        data.plotOn(mass)
        sim.plotOn(mass,RooFit.LineColor(kRed), RooFit.ProjWData(data,True))

        comp_otherfid = ''
        for bin in range(nBins):            
            if bin==recobin: continue
            for fState in fStates:
                comp_otherfid += "shapeSig_trueH"+fState+"Bin"+str(bin)+"_ch"+channel[fState]+"_ch"+str(recobin+1)+","
        comp_otherfid = comp_otherfid.rstrip(',')

        comp_out = ''
        comp_fake = ''
        comp_zz = ''
        comp_zx = ''
        for fState in fStates:
            comp_out += "shapeBkg_out_trueH_ch"+channel[fState]+"_ch"+str(recobin+1)+","
            comp_fake += "shapeBkg_fakeH_ch"+channel[fState]+"_ch"+str(recobin+1)+","
            comp_zz += "shapeBkg_bkg_ggzz_ch"+channel[fState]+"_ch"+str(recobin+1)+",shapeBkg_bkg_qqzz_ch"+channel[fState]+"_ch"+str(recobin+1)+","
            comp_zx += "shapeBkg_bkg_zjets_ch"+channel[fState]+"_ch"+str(recobin+1)+","
        comp_out = comp_out.rstrip(',')
        comp_fake = comp_fake.rstrip(',')
        comp_zz = comp_zz.rstrip(',')
        comp_zx = comp_zx.rstrip(',')
        sim.plotOn(mass, RooFit.LineColor(kGray+2), RooFit.LineStyle(2), RooFit.Components(comp_zx+","+comp_zz+","+comp_fake+","+comp_otherfid+","+comp_out), RooFit.ProjWData(data,True))
        sim.plotOn(mass, RooFit.LineColor(kRed), RooFit.LineStyle(2), RooFit.Components(comp_zx+","+comp_zz+","+comp_fake+","+comp_otherfid), RooFit.ProjWData(data,True))
        sim.plotOn(mass, RooFit.LineColor(kOrange), RooFit.Components(comp_zx+","+comp_zz+","+comp_fake), RooFit.ProjWData(data,True))
        sim.plotOn(mass, RooFit.LineColor(kAzure-3), RooFit.Components(comp_zx+","+comp_zz), RooFit.ProjWData(data,True))
        sim.plotOn(mass, RooFit.LineColor(kGreen+3), RooFit.Components(comp_zx), RooFit.ProjWData(data,True))
        data.plotOn(mass)
    else:
        sbin = "ch"+channel[fstate]+"_ch"+str(recobin+1)
        data = data.reduce(RooFit.Cut("CMS_channel==CMS_channel::"+sbin)) 
        data.plotOn(mass)
        comp_otherfid = ''
        for bin in range(nBins):
            if bin==recobin: continue
            comp_otherfid += "shapeSig_trueH"+fstate+"Bin"+str(bin)+"_"+sbin+","
        comp_otherfid = comp_otherfid.rstrip(',')
        pdfi.plotOn(mass, RooFit.Components("pdf_binch"+channel[fstate]+"_ch"+str(recobin+1)),RooFit.LineColor(kRed), RooFit.Slice(CMS_channel,sbin),RooFit.ProjWData(RooArgSet(CMS_channel),data,True))
        pdfi.plotOn(mass, RooFit.LineColor(kGray+2), RooFit.LineStyle(2), RooFit.Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+sbin+",shapeBkg_fakeH_"+sbin+",shapeBkg_out_trueH_"+sbin+","+comp_otherfid), RooFit.Slice(CMS_channel,sbin),RooFit.ProjWData(RooArgSet(CMS_channel),data,True))
        pdfi.plotOn(mass, RooFit.LineColor(kRed), RooFit.LineStyle(2), RooFit.Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+sbin+",shapeBkg_fakeH_"+sbin+","+comp_otherfid), RooFit.Slice(CMS_channel,sbin),RooFit.ProjWData(RooArgSet(CMS_channel),data,True))
        pdfi.plotOn(mass, RooFit.LineColor(kOrange), RooFit.Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+sbin+",shapeBkg_fakeH_"+sbin), RooFit.Slice(CMS_channel,sbin),RooFit.ProjWData(RooArgSet(CMS_channel),data,True))
        pdfi.plotOn(mass, RooFit.LineColor(kAzure-3), RooFit.Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+sbin), RooFit.Slice(CMS_channel,sbin),RooFit.ProjWData(RooArgSet(CMS_channel),data,True))
        pdfi.plotOn(mass, RooFit.LineColor(kGreen+3), RooFit.Components("shapeBkg_bkg_zjets_"+sbin), RooFit.Slice(CMS_channel,sbin),RooFit.ProjWData(RooArgSet(CMS_channel),data,True))
        data.plotOn(mass)
        
    gStyle.SetOptStat(0)

    c = TCanvas("c","c",1000,800)
    c.cd()

    #dummy = TH1D("","",1,105.6,140.6)
    dummy = TH1D("","",1,105.0,140.0)
    dummy.SetBinContent(1,2)
    dummy.SetFillColor(0)
    dummy.SetLineColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.SetMarkerColor(0) 
    dummy.GetYaxis().SetTitle("Events / (2.33 GeV)")
    dummy.GetXaxis().SetTitle("m_{"+fstate.replace("mu","#mu")+"} [GeV]")
    if (opt.UNBLIND):
        dummy.SetMaximum(max(3.0*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        if (obsName=="massZ2" and recobin==0): dummy.SetMaximum(max(6.0*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),3.5))
    else:
        dummy.SetMaximum(max(1.5*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        if (obsName=="massZ2" and recobin==0): dummy.SetMaximum(max(6.0*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),3.5))
        #dummy.SetMaximum(0.5*max(n_trueH_asimov[fstate],n_zz_asimov[fstate],2.5))
    dummy.SetMinimum(0.0)
    dummy.Draw()

    dummy_data = TH1D()
    dummy_data.SetMarkerColor(kBlack)
    dummy_data.SetMarkerStyle(20)
    dummy_fid = TH1D()
    dummy_fid.SetLineColor(kRed)
    dummy_fid.SetLineWidth(2)
    dummy_other = TH1D()
    dummy_other.SetLineColor(kRed)
    dummy_other.SetLineWidth(2)
    dummy_other.SetLineStyle(2)
    dummy_out = TH1D()
    dummy_out.SetLineColor(kGray+2)
    dummy_out.SetLineWidth(2)
    dummy_out.SetLineStyle(2)
    dummy_fake = TH1D()
    dummy_fake.SetLineColor(kOrange)
    dummy_fake.SetLineWidth(2)
    dummy_zz = TH1D()
    dummy_zz.SetLineColor(kAzure-3)
    dummy_zz.SetLineWidth(2)
    dummy_zx = TH1D()
    dummy_zx.SetLineColor(kGreen+3)
    dummy_zx.SetLineWidth(2)


    legend = TLegend(.20,.41,.53,.89)
    if(not opt.UNBLIND):
        legend.AddEntry(dummy_data,"Asimov Data (SM m(H) = "+opt.ASIMOVMASS+" GeV)","ep")
    else:
        legend.AddEntry(dummy_data,"Data","ep")   
    legend.AddEntry(dummy_fid,"N_{fid.}^{fit} = %.2f (exp. = %.2f)"%(n_trueH_modelfit[fstate],n_trueH_asimov[fstate]), "l")
    legend.AddEntry(dummy_other,"N_{other fid.}^{fit} = %.2f (exp = %.2f)"%(n_trueH_otherfid_modelfit[fstate],n_trueH_otherfid_asimov[fstate]), "l")
    legend.AddEntry(dummy_out, "N_{out}^{fit} = %.2f (exp. = %.2f)"%(n_out_trueH_modelfit[fstate],n_out_trueH_asimov[fstate]), "l")
    legend.AddEntry(dummy_fake, "N_{wrong}^{fit} = %.2f (exp. = %.2f)"%(n_fakeH_modelfit[fstate],n_fakeH_asimov[fstate]), "l")
    legend.AddEntry(dummy_zz, "N_{ZZ}^{fit} = %.2f (exp. = %.2f)"%(n_zz_modelfit[fstate],n_zz_asimov[fstate]), "l")
    legend.AddEntry(dummy_zx, "N_{Z+X}^{fit} = %.2f (exp. = %.2f)"%(n_zjets_modelfit[fstate],n_zjets_asimov[fstate]), "l")

    #legend.AddEntry(dummy_fid,"Fiducial Signal", "l")
    #legend.AddEntry(dummy_other,"Other Bin Fiducial Signal", "l")
    #legend.AddEntry(dummy_out, "Non-fiducial Signal", "l")
    #legend.AddEntry(dummy_fake, "Non-resonant Signal", "l")
    #legend.AddEntry(dummy_zz, "ZZ", "l")
    #legend.AddEntry(dummy_zx, "Z+X", "l")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw()

    mass.Draw("same")

    if (obsName=="pT4l"):
        label="p_{T}^{H}"
        unit="GeV"
    elif (obsName=="massZ2"):
        label = "m(Z_{2})"
        unit = "GeV"
    elif (obsName=="massZ1"):
        label = "m(Z_{1})"
        unit = "GeV"
    elif (obsName=="nJets" or obsName=="njets_reco_pt30_eta4p7"):
        label = "N(jets) |#eta|<4.7"
        unit = ""
    elif (obsName=="njets_pt30_eta2p5"):
        label = "N(jets) |#eta|<2.5"
        unit = ""
    elif (obsName=='pt_leadingjet_pt30_eta4p7'):
        label = "p_{T}(jet)"
        unit = "GeV"
    elif ((obsName=='pt_leadingjet_pt30_eta2p5') | (obsName=='pT1j')):
        label = "p_{T}(jet) |#eta|<2.5"
        unit = "GeV"
    elif (obsName=='absdeltarapidity_hleadingjet_pt30_eta4p7'):
        label = "|y(H)-y(jet)|"
        unit = ""
    elif (obsName=='absdeltarapidity_hleadingjet_pt30_eta2p5'):
        label = "|y(H)-y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=='absrapidity_leadingjet_pt30_eta4p7'):
        label = "|y(jet)|"
        unit = ""
    elif (obsName=='absrapidity_leadingjet_pt30_eta2p5'):
        label = "|y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=="rapidity4l"):
        label = "|y^{H}|"
        unit = ""
    elif (obsName=="cosThetaStar"):
        label = "|cos#theta*|"
        unit = ""
    elif (obsName=="cosTheta1"):
        label = "|cos#theta_{1}|"
        unit = ""
    elif (obsName=="cosTheta2"):
        label = "|cos#theta_{2}|"
        unit = ""
    elif (obsName=="Phi"):
        label = "|#Phi|"
        unit = ""
    elif (obsName=="Phi1"):
        label = "|#Phi_{1}|"
        unit = ""
    elif (obsName=="mass4l"):
        label = "inclusive"
        unit = ""
    else:
        label = obsName
        unit = ""
                
    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.87, 0.95,lumi+" fb^{-1} (#sqrt{s} = 13 TeV)")
    latex2.SetTextSize(0.8*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.19, 0.95, "CMS")
    latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    latex2.DrawLatex(0.30, 0.95, "Preliminary")
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.45*c.GetTopMargin())
    #latex2.DrawLatex(0.20,0.85, observableBins[recobin]+" "+unit+" < "+label+" < "+observableBins[recobin+1]+" "+unit+"    Unfolding model: "+modelName.replace("_"," ")+" GeV")                                                                                                            
    latex2.DrawLatex(0.65,0.85, observableBins[recobin]+" "+unit+" < "+label+" < "+observableBins[recobin+1]+" "+unit)

    checkDir("plots")
    checkDir("plots/"+obsName)
 
    if (not opt.UNBLIND):
        checkDir("plots/"+obsName+"/asimov")
        c.SaveAs("plots/"+obsName+"/asimov/asimovdata_"+asimovDataModel+"_"+year+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".pdf")
        c.SaveAs("plots/"+obsName+"/asimov/asimovdata_"+asimovDataModel+"_"+year+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".png")
    else:
        checkDir("plots/"+obsName+"/data")
        c.SaveAs("plots/"+obsName+"/data/data_unfoldwith_"+modelName+"_"+year+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".pdf")
        c.SaveAs("plots/"+obsName+"/data/data_unfoldwith_"+modelName+"_"+year+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".png")


fStates = ["4e","4mu","2e2mu","4l"]
for fState in fStates:
    for recobin in range(len(observableBins)-1):
        plotAsimov_sim(asimovDataModel, asimovPhysicalModel, modelName, physicalModel, obsName, fState, observableBins, recobin)
