import sys, os, string, re, pwd,  ast, optparse, shlex, time
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
    parser.add_option('',   '--asimovMass',dest='ASIMOVMASS',type='string',default='125.38', help='Asimov Mass')
    parser.add_option('',   '--unfoldModel',dest='UNFOLD',type='string',default='SM_125', help='Name of the unfolding model')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "inclusive", "pT", "eta", "Njets"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='Use results from fixed fraction fit, default is False')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='int',default=105.0,   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='int',default=140.0,   help='Upper bound for m4l')
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

def generateName(_year, _fStateNumber, _recobin, _fState, _bin, _physicalModel, _observableBins, _obsName):
    _years = {"1":"2016", "2":"2017", "3":"2018"}
    if _physicalModel != 'v3':
        # In case of mass4l we have only one bin, so the third 'chan' part is not included in the name of the function
        if (_obsName!='mass4l' and _obsName!='mass4l_zzfloating'):
            binName = "ch"+_year+"_ch"+_fStateNumber+"_ch"+str(_recobin+1)
            procName = _fState+"Bin"+str(_bin)
        else:
            binName = "ch"+_year+"_ch"+_fStateNumber
            procName = _fState+"Bin"+str(_bin)
    else:
        _obsName_v3 = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta4p7': 'NJ'}
        if _obsName not in _obsName_v3:
            _obsName_v3[_obsName] = _obsName

        if '_' in obsName and not 'floating' in obsName and not 'kL' in obsName and obsName != 'njets_pt30_eta4p7':
            _recobin_final = str(_observableBins[_recobin][0]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin][1]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin][2]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin][3]).replace('.', 'p').replace('-','m')
            _genbin_final = str(_observableBins[_bin][0]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin][1]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin][2]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin][3]).replace('.', 'p').replace('-','m')
        else:
            _recobin_final = str(_observableBins[_recobin]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin+1]).replace('.', 'p').replace('-','m')
            if int(_observableBins[_recobin+1]) > 1000:
                _recobin_final = 'GT'+str(int(_observableBins[_recobin]))

            _genbin_final = str(_observableBins[_bin]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin+1]).replace('.', 'p').replace('-','m')
            if int(_observableBins[_bin+1]) > 1000:
                _genbin_final = 'GT'+str(int(_observableBins[_bin]))


        binName = "hzz_" + _obsName_v3[_obsName] + "_" + _recobin_final + "_cat" + _fState + "_" + _years[_year]
        procName = _obsName_v3[_obsName] + "_" + _genbin_final

    return binName, procName



def plotAsimov_sim(modelName, physicalModel, obsName, fstate, observableBins, recobin):
    sourcedir = opt.SOURCEDIR
    theorymass = opt.THEORYMASS
    year = opt.YEAR
    if year == '2016':
        lumi = '36.3'
        years = [""]
    elif year == '2017':
        lumi = '41.5'
        years = [""]
    elif year == '2018':
        lumi = '59.7'
        years = [""]
    else:
        lumi = '138'
        years = ["1", "2", "3"]

    # nBins = len(observableBins)
    # if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries

    channel = {"4mu":"1", "4e":"2", "2e2mu":"3", "4l":"2"} # 4l is dummy, won't be used
    run = {"2016":"1", "2017":"2", "2018":"3", "Full": "2"}
    SignalNames = {"v3":"smH_", "v2":"trueH", "v4":"trueH", "kLambda": "trueH"}
    CombNames = {"v3":"nonResH", "v2":"fakeH", "v4":"fakeH", "kLambda":"fakeH"}
    OutNames = {"v3":"OutsideAcceptance", "v2":"out_trueH", "v4":"out_trueH", "kLambda":"out_trueH"}

    # Load some libraries
    ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
    ROOT.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
    ROOT.gSystem.AddIncludePath("-Iinclude/")
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

    if(not opt.UNBLIND):
        theorymass = theorymass + '.123456'
    if physicalModel == 'v3':
        fname = 'higgsCombine_'+obsName+'_r_smH_0.MultiDimFit.mH'+theorymass+'.root'
    elif physicalModel == 'kLambda':
        fname = 'higgsCombine_'+obsName+'.MultiDimFit.mH'+theorymass+'.root'
    else:
        fname = 'higgsCombine_'+obsName+'_r2e2muBin0.MultiDimFit.mH'+theorymass+'.root'
    print 'combine file: ', fname

    f_asimov = TFile(sourcedir + fname, "READ")

    if (not opt.UNBLIND):
        data = f_asimov.Get("toys/toy_asimov");
    w_asimov = f_asimov.Get("w")
    if (opt.UNBLIND):
        data = w_asimov.data("data_obs")
    w_asimov.loadSnapshot("clean")

    ##### ------------------------ Normalisation values for Asimov dataset ------------------------ #####
    trueH_asimov = {}
    zjets_asimov = {}
    ggzz_asimov = {}
    fakeH_asimov = {}
    out_trueH_asimov = {}
    qqzz_asimov = {}
    n_trueH_asimov = {}
    for year in years:
        n_trueH_asimov["4l_"+year] = 0.0
    n_trueH_otherfid_asimov = {}
    for year in years:
        n_trueH_otherfid_asimov["4l_"+year] = 0.0
    n_zjets_asimov = {}
    for year in years:
        n_zjets_asimov["4l_"+year] = 0.0
    n_ggzz_asimov = {}
    for year in years:
        n_ggzz_asimov["4l_"+year] = 0.0
    n_fakeH_asimov = {}
    for year in years:
        n_fakeH_asimov["4l_"+year] = 0.0
    n_out_trueH_asimov = {}
    for year in years:
        n_out_trueH_asimov["4l_"+year] = 0.0
    n_qqzz_asimov = {}
    for year in years:
        n_qqzz_asimov["4l_"+year] = 0.0
    n_zz_asimov = {}
    for year in years:
        n_zz_asimov["4l_"+year] = 0.0

    fStates = ['4mu','4e','2e2mu']
    for year in years:
        for fState in fStates:
            for bin in range(nBins):
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                trueH_asimov[fState+"_"+year+"Bin"+str(bin)] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)
                # print (fState+"_"+year+"Bin"+str(bin))
                # print ("n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)
                # print (trueH_asimov[fState+"_"+year+"Bin"+str(bin)].getVal())
                # print ()

            zjets_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_bkg_zjets")
            ggzz_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_bkg_ggzz")
            fakeH_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_"+CombNames[physicalModel])
            out_trueH_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_"+OutNames[physicalModel])
            qqzz_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_bkg_qqzz")

            n_trueH_otherfid_asimov[fState+"_"+year] = 0.0
            for bin in range(nBins):
                if (bin==recobin): n_trueH_asimov[fState+"_"+year] = trueH_asimov[fState+"_"+year+"Bin"+str(bin)].getVal()
                else: n_trueH_otherfid_asimov[fState+"_"+year] += trueH_asimov[fState+"_"+year+"Bin"+str(bin)].getVal()
            n_zjets_asimov[fState+"_"+year] = zjets_asimov[fState+"_"+year].getVal()
            n_ggzz_asimov[fState+"_"+year] = ggzz_asimov[fState+"_"+year].getVal()
            n_fakeH_asimov[fState+"_"+year] = fakeH_asimov[fState+"_"+year].getVal()
            n_out_trueH_asimov[fState+"_"+year] = out_trueH_asimov[fState+"_"+year].getVal()
            n_qqzz_asimov[fState+"_"+year] = qqzz_asimov[fState+"_"+year].getVal()
            n_zz_asimov[fState+"_"+year] = n_ggzz_asimov[fState+"_"+year]+n_qqzz_asimov[fState+"_"+year]
            n_trueH_asimov["4l_"+year] += n_trueH_asimov[fState+"_"+year]
            n_trueH_otherfid_asimov["4l_"+year] += n_trueH_otherfid_asimov[fState+"_"+year]
            n_zjets_asimov["4l_"+year] += zjets_asimov[fState+"_"+year].getVal()
            n_ggzz_asimov["4l_"+year] += ggzz_asimov[fState+"_"+year].getVal()
            n_fakeH_asimov["4l_"+year] += fakeH_asimov[fState+"_"+year].getVal()
            n_out_trueH_asimov["4l_"+year] += out_trueH_asimov[fState+"_"+year].getVal()
            n_qqzz_asimov["4l_"+year] += qqzz_asimov[fState+"_"+year].getVal()
            n_zz_asimov["4l_"+year] += n_ggzz_asimov[fState+"_"+year]+n_qqzz_asimov[fState+"_"+year]

    # Sum over the years
    for fState in fStates:
        n_trueH_asimov[fState] = 0.0
        n_trueH_otherfid_asimov[fState] = 0.0
        n_fakeH_asimov[fState] = 0.0
        n_out_trueH_asimov[fState] = 0.0
        n_zz_asimov[fState] = 0.0
        n_qqzz_asimov[fState] = 0.0
        n_zjets_asimov[fState] = 0.0
        for year in ["1", "2", "3"]:
            n_trueH_asimov[fState] += n_trueH_asimov[fState+"_"+year]
            n_trueH_otherfid_asimov[fState] += n_trueH_otherfid_asimov[fState+"_"+year]
            n_fakeH_asimov[fState] += n_fakeH_asimov[fState+"_"+year]
            n_out_trueH_asimov[fState] += n_out_trueH_asimov[fState+"_"+year]
            n_zz_asimov[fState] += n_zz_asimov[fState+"_"+year]
            n_qqzz_asimov[fState] += n_qqzz_asimov[fState+"_"+year]
            n_zjets_asimov[fState] += n_zjets_asimov[fState+"_"+year]

    # Sum over the final states
    n_trueH_asimov["4l"] = 0.0
    n_trueH_otherfid_asimov["4l"] = 0.0
    n_fakeH_asimov["4l"] = 0.0
    n_out_trueH_asimov["4l"] = 0.0
    n_zjets_asimov["4l"] = 0.0
    n_zz_asimov["4l"] = 0.0
    n_qqzz_asimov["4l"] = 0.0
    for year in ["1", "2", "3"]:
        n_trueH_asimov["4l"] += n_trueH_asimov["4l_"+year]
        n_trueH_otherfid_asimov["4l"] += n_trueH_otherfid_asimov["4l_"+year]
        n_fakeH_asimov["4l"] += n_fakeH_asimov["4l_"+year]
        n_out_trueH_asimov["4l"] += n_out_trueH_asimov["4l_"+year]
        n_zjets_asimov["4l"] += n_zjets_asimov["4l_"+year]
        n_qqzz_asimov["4l"] += n_zz_asimov["4l_"+year]
        n_zz_asimov["4l"] += n_zz_asimov["4l_"+year]

    ##### ------------------------ Normalisation values for modelfit ------------------------ #####
    f_modelfit = TFile(sourcedir + fname, "READ")
    w_modelfit = f_modelfit.Get("w")
    sim = w_modelfit.pdf("model_s")
    CMS_zz4l_mass = w_modelfit.var("CMS_zz4l_mass")
    w_modelfit.loadSnapshot("MultiDimFit")

    trueH_modelfit = {}
    zjets_modelfit = {}
    ggzz_modelfit = {}
    fakeH_modelfit = {}
    out_trueH_modelfit = {}
    qqzz_modelfit = {}
    n_trueH_modelfit = {}
    for year in years:
        n_trueH_modelfit["4l_"+year] = 0.0
    n_trueH_otherfid_modelfit = {}
    for year in years:
        n_trueH_otherfid_modelfit["4l_"+year] = 0.0
    n_zjets_modelfit = {}
    for year in years:
        n_zjets_modelfit["4l_"+year] = 0.0
    n_ggzz_modelfit = {}
    for year in years:
        n_ggzz_modelfit["4l_"+year] = 0.0
    n_fakeH_modelfit = {}
    for year in years:
        n_fakeH_modelfit["4l_"+year] = 0.0
    n_out_trueH_modelfit = {}
    for year in years:
        n_out_trueH_modelfit["4l_"+year] = 0.0
    n_qqzz_modelfit = {}
    for year in years:
        n_qqzz_modelfit["4l_"+year] = 0.0
    n_zz_modelfit = {}
    for year in years:
        n_zz_modelfit["4l_"+year] = 0.0

    fStates = ['4mu','4e','2e2mu']
    for year in ["1", "2", "3"]:
        for fState in fStates:
            for bin in range(nBins):
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                trueH_modelfit[fState+"_"+year+"Bin"+str(bin)] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)
                # print (fState+"_"+year+"Bin"+str(bin))
                # print ("n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)
                # print (trueH_modelfit[fState+"_"+year+"Bin"+str(bin)].getVal())

            zjets_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_bkg_zjets")
            ggzz_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_bkg_ggzz")
            fakeH_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_"+CombNames[physicalModel])
            out_trueH_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_"+OutNames[physicalModel])
            qqzz_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_bkg_qqzz")

            n_trueH_otherfid_modelfit[fState+"_"+year] = 0.0
            for bin in range(nBins):
                if (bin==recobin): n_trueH_modelfit[fState+"_"+year] = trueH_modelfit[fState+"_"+year+"Bin"+str(bin)].getVal()
                else: n_trueH_otherfid_modelfit[fState+"_"+year] += trueH_modelfit[fState+"_"+year+"Bin"+str(bin)].getVal()
            n_zjets_modelfit[fState+"_"+year] = zjets_modelfit[fState+"_"+year].getVal()
            n_ggzz_modelfit[fState+"_"+year] = ggzz_modelfit[fState+"_"+year].getVal()
            n_fakeH_modelfit[fState+"_"+year] = fakeH_modelfit[fState+"_"+year].getVal()
            n_out_trueH_modelfit[fState+"_"+year] = out_trueH_modelfit[fState+"_"+year].getVal()
            n_qqzz_modelfit[fState+"_"+year] = qqzz_modelfit[fState+"_"+year].getVal()
            n_zz_modelfit[fState+"_"+year] = n_ggzz_modelfit[fState+"_"+year]+n_qqzz_modelfit[fState+"_"+year]
            n_trueH_modelfit["4l_"+year] += n_trueH_modelfit[fState+"_"+year]
            n_trueH_otherfid_modelfit["4l_"+year] += n_trueH_otherfid_modelfit[fState+"_"+year]
            n_zjets_modelfit["4l_"+year] += zjets_modelfit[fState+"_"+year].getVal()
            n_ggzz_modelfit["4l_"+year] += ggzz_modelfit[fState+"_"+year].getVal()
            n_fakeH_modelfit["4l_"+year] += fakeH_modelfit[fState+"_"+year].getVal()
            n_out_trueH_modelfit["4l_"+year] += out_trueH_modelfit[fState+"_"+year].getVal()
            n_qqzz_modelfit["4l_"+year] += qqzz_modelfit[fState+"_"+year].getVal()
            n_zz_modelfit["4l_"+year] += n_ggzz_modelfit[fState+"_"+year]+n_qqzz_modelfit[fState+"_"+year]

    # Sum over the years
    for fState in fStates:
        n_trueH_modelfit[fState] = 0.0
        n_trueH_otherfid_modelfit[fState] = 0.0
        n_fakeH_modelfit[fState] = 0.0
        n_out_trueH_modelfit[fState] = 0.0
        n_zz_modelfit[fState] = 0.0
        n_qqzz_modelfit[fState] = 0.0
        n_zjets_modelfit[fState] = 0.0
        for year in ["1", "2", "3"]:
            n_trueH_modelfit[fState] += n_trueH_modelfit[fState+"_"+year]
            n_trueH_otherfid_modelfit[fState] += n_trueH_otherfid_modelfit[fState+"_"+year]
            n_fakeH_modelfit[fState] += n_fakeH_modelfit[fState+"_"+year]
            n_out_trueH_modelfit[fState] += n_out_trueH_modelfit[fState+"_"+year]
            n_zz_modelfit[fState] += n_zz_modelfit[fState+"_"+year]
            n_qqzz_modelfit[fState] += n_qqzz_modelfit[fState+"_"+year]
            n_zjets_modelfit[fState] += n_zjets_modelfit[fState+"_"+year]

    # Sum over the final states
    n_trueH_modelfit["4l"] = 0.0
    n_trueH_otherfid_modelfit["4l"] = 0.0
    n_fakeH_modelfit["4l"] = 0.0
    n_out_trueH_modelfit["4l"] = 0.0
    n_zjets_modelfit["4l"] = 0.0
    n_zz_modelfit["4l"] = 0.0
    n_qqzz_modelfit["4l"] = 0.0
    for year in ["1", "2", "3"]:
        n_trueH_modelfit["4l"] += n_trueH_modelfit["4l_"+year]
        n_trueH_otherfid_modelfit["4l"] += n_trueH_otherfid_modelfit["4l_"+year]
        n_fakeH_modelfit["4l"] += n_fakeH_modelfit["4l_"+year]
        n_out_trueH_modelfit["4l"] += n_out_trueH_modelfit["4l_"+year]
        n_zjets_modelfit["4l"] += n_zjets_modelfit["4l_"+year]
        n_qqzz_modelfit["4l"] += n_zz_modelfit["4l_"+year]
        n_zz_modelfit["4l"] += n_zz_modelfit["4l_"+year]


    ##### ------------------------ Data ------------------------ #####
    CMS_channel = w.cat("CMS_channel")
    mass = w.var("CMS_zz4l_mass").frame(RooFit.Bins(30))

    if (fstate=="4l"):
        datacut = ''
        for year in ["1", "2", "3"]:
            for fState in fStates:
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                if(obsName!='mass4l' and obsName!='mass4l_zzfloating'):
                    datacut += "CMS_channel==CMS_channel::"+bin_name+" || "
                else:
                    datacut += "CMS_channel==CMS_channel::"+bin_name+" || "
        datacut = datacut.rstrip(" || ")
        data = data.reduce(RooFit.Cut(datacut))
        data.plotOn(mass)
        sim.plotOn(mass,RooFit.LineColor(kOrange-3), RooFit.ProjWData(data,True))
    else:
        datacut = ''
        for year in ["1", "2", "3"]:
            bin_name, process_name = generateName(year, channel[fstate], recobin, fstate, bin, physicalModel, observableBins, obsName)
            if(obsName!='mass4l' and obsName!='mass4l_zzfloating'):
                datacut += "CMS_channel==CMS_channel::"+bin_name+" || "
            else:
                datacut += "CMS_channel==CMS_channel::"+bin_name+" || "
        datacut = datacut.rstrip(" || ")
        data = data.reduce(RooFit.Cut(datacut))
        data.plotOn(mass)
        sim.plotOn(mass,RooFit.LineColor(kOrange-3), RooFit.ProjWData(data,True))

    ##### ------------------------ Shapes ------------------------ #####
    if (fstate!="4l"):
        comp_otherfid = ''
        for bin in range(nBins):
            if bin==recobin: continue
            for year in ["1", "2", "3"]:
                bin_name, process_name = generateName(year, channel[fstate], recobin, fstate, bin, physicalModel, observableBins, obsName)
                comp_otherfid += "shapeSig_"+SignalNames[physicalModel]+process_name+"_"+bin_name+","
        comp_otherfid = comp_otherfid.rstrip(',')

        comp_out = ''
        comp_fake = ''
        comp_zz = ''
        comp_zx = ''
        for year in ["1", "2", "3"]:
            bin_name, process_name = generateName(year, channel[fstate], recobin, fstate, 0, physicalModel, observableBins, obsName)
            comp_out += "shapeBkg_"+OutNames[physicalModel]+"_"+bin_name+","
            comp_fake += "shapeBkg_"+CombNames[physicalModel]+"_"+bin_name+","
            comp_zz += "shapeBkg_bkg_ggzz_"+bin_name+",shapeBkg_bkg_qqzz_"+bin_name+","
            comp_zx += "shapeBkg_bkg_zjets_"+bin_name+","
        comp_out = comp_out.rstrip(',')
        comp_fake = comp_fake.rstrip(',')
        comp_zz = comp_zz.rstrip(',')
        comp_zx = comp_zx.rstrip(',')

    else: #if fstate!=4 a loop over the threee final states is necessary
        comp_otherfid = ''
        for bin in range(nBins):
            if bin==recobin: continue
            for year in ["1", "2", "3"]:
                for fState in fStates:
                    bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                    comp_otherfid += "shapeSig_"+SignalNames[physicalModel]+process_name+"_"+bin_name+","
        comp_otherfid = comp_otherfid.rstrip(',')

        comp_out = ''
        comp_fake = ''
        comp_zz = ''
        comp_zx = ''
        for year in ["1", "2", "3"]:
            for fState in fStates:
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, 0, physicalModel, observableBins, obsName)
                comp_out += "shapeBkg_"+OutNames[physicalModel]+"_"+bin_name+","
                comp_fake += "shapeBkg_"+CombNames[physicalModel]+"_"+bin_name+","
                comp_zz += "shapeBkg_bkg_ggzz_"+bin_name+",shapeBkg_bkg_qqzz_"+bin_name+","
                comp_zx += "shapeBkg_bkg_zjets_"+bin_name+","
        comp_out = comp_out.rstrip(',')
        comp_fake = comp_fake.rstrip(',')
        comp_zz = comp_zz.rstrip(',')
        comp_zx = comp_zx.rstrip(',')

    sim.plotOn(mass, RooFit.LineColor(kGreen+2), RooFit.Components(comp_zx+","+comp_zz+","+comp_fake+","+comp_otherfid+","+comp_out), RooFit.ProjWData(data,True))
    sim.plotOn(mass, RooFit.LineColor(kOrange-3), RooFit.LineStyle(2), RooFit.Components(comp_zx+","+comp_zz+","+comp_fake+","+comp_otherfid), RooFit.ProjWData(data,True))
    sim.plotOn(mass, RooFit.LineColor(kAzure-3), RooFit.Components(comp_zx+","+comp_zz+","+comp_fake), RooFit.ProjWData(data,True))
    sim.plotOn(mass, RooFit.LineColor(kViolet), RooFit.Components(comp_zx+","+comp_zz), RooFit.ProjWData(data,True))
    sim.plotOn(mass, RooFit.LineColor(kViolet+2), RooFit.Components(comp_zx), RooFit.ProjWData(data,True))
    data.plotOn(mass)


    ##### ------------------------ Plot ------------------------ #####
    gStyle.SetOptStat(0)

    c = TCanvas("c","c",1000,800)
    c.cd()

    #dummy = TH1D("","",1,105.6,160.6)
    dummy = TH1D("","",1,opt.LOWER_BOUND,opt.UPPER_BOUND)
    dummy.SetBinContent(1,2)
    dummy.SetFillColor(0)
    dummy.SetLineColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.SetMarkerColor(0)
    dummy.GetYaxis().SetTitle("Events / (1.83 GeV)")
    dummy.GetXaxis().SetTitle("m_{"+fstate.replace("mu","#mu")+"} [GeV]")
    if (opt.UNBLIND):
        dummy.SetMaximum(max(0.8*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
    else:
        # if fstate=='4e': dummy.SetMaximum(max(1*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        # elif fstate=='4l': dummy.SetMaximum(max(0.2*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        # else: dummy.SetMaximum(max(0.5*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        # if (obsName=="massZ2" and recobin==0): dummy.SetMaximum(max(3.0*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),3.5))
        dummy.SetMaximum(max(1*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        if (obsName=="massZ2" and recobin==0): dummy.SetMaximum(max(3.0*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),3.5))
        #dummy.SetMaximum(0.5*max(n_trueH_asimov[fstate],n_zz_asimov[fstate],2.5))
    dummy.SetMinimum(0.0)
    dummy.Draw()

    dummy_data = TH1D()
    dummy_data.SetMarkerColor(kBlack)
    dummy_data.SetMarkerStyle(20)
    dummy_fid = TH1D()
    dummy_fid.SetLineColor(kOrange-3)
    dummy_fid.SetLineWidth(2)
    dummy_other = TH1D()
    dummy_other.SetLineColor(kOrange-3)
    dummy_other.SetLineWidth(2)
    dummy_other.SetLineStyle(2)
    dummy_out = TH1D()
    dummy_out.SetLineColor(kGreen+2)
    dummy_out.SetLineWidth(2)
    dummy_fake = TH1D()
    dummy_fake.SetLineColor(kAzure-3)
    dummy_fake.SetLineWidth(2)
    dummy_zz = TH1D()
    dummy_zz.SetLineColor(kViolet)
    dummy_zz.SetLineWidth(2)
    dummy_zx = TH1D()
    dummy_zx.SetLineColor(kViolet+2)
    dummy_zx.SetLineWidth(2)

    # if opt.UPPER_BOUND == 160: legend = TLegend(.60,.41,.93,.89)
    # else: legend = TLegend(.20,.41,.53,.89)
    legend = TLegend(.60,.60,.93,.89)
    # legend.AddEntry(dummy_data,"Asimov Data (SM m(H) = "+opt.ASIMOVMASS+" GeV)","ep")
    if (not opt.UNBLIND):
       legend.AddEntry(dummy_data,"Asimov Data (SM m(H) = "+opt.ASIMOVMASS+" GeV)","ep")
    else:
       legend.AddEntry(dummy_data,"Data","ep")
    # legend.AddEntry(dummy_fid,"N_{fid.}^{fit} = %.2f (exp. = %.2f)"%(n_trueH_modelfit[fstate],n_trueH_asimov[fstate]), "l")
    # legend.AddEntry(dummy_other,"N_{other fid.}^{fit} = %.2f (exp = %.2f)"%(n_trueH_otherfid_modelfit[fstate],n_trueH_otherfid_asimov[fstate]), "l")
    # legend.AddEntry(dummy_out, "N_{out}^{fit} = %.2f (exp. = %.2f)"%(n_out_trueH_modelfit[fstate],n_out_trueH_asimov[fstate]), "l")
    # legend.AddEntry(dummy_fake, "N_{wrong}^{fit} = %.2f (exp. = %.2f)"%(n_fakeH_modelfit[fstate],n_fakeH_asimov[fstate]), "l")
    # legend.AddEntry(dummy_zz, "N_{ZZ}^{fit} = %.2f (exp. = %.2f)"%(n_zz_modelfit[fstate],n_zz_asimov[fstate]), "l")
    # legend.AddEntry(dummy_zx, "N_{Z+X}^{fit} = %.2f (exp. = %.2f)"%(n_zjets_modelfit[fstate],n_zjets_asimov[fstate]), "l")

    legend.AddEntry(dummy_fid,"N_{fid}", "l")
    legend.AddEntry(dummy_out, "N_{nonfid}", "l")
    legend.AddEntry(dummy_fake, "N_{nonres}", "l")
    legend.AddEntry(dummy_zz, "N_{ZZ}", "l")
    legend.AddEntry(dummy_zx, "N_{ZX}", "l")

    # print('qqZZ', fstate, n_qqzz_modelfit[fstate])
    # print('ggZZ', fstate, n_ggzz_modelfit[fstate])

    #legend.SetTextSize(0.03)
    #if (not opt.UNBLIND):
    #    legend.AddEntry(dummy_data,"Asimov Data (SM m(H) = "+opt.ASIMOVMASS+" GeV)","ep")
    #else:
    #    legend.AddEntry(dummy_data,"Data","ep")
    #legend.AddEntry(dummy_fid,"Fiducial Signal", "l")
    #legend.AddEntry(dummy_other,"Other Bin Fiducial Signal", "l")
    #legend.AddEntry(dummy_out, "Non-fiducial Signal", "l")
    #legend.AddEntry(dummy_fake, "Non-resonant Signal", "l")
    #legend.AddEntry(dummy_zz, "ZZ", "l")
    #legend.AddEntry(dummy_zx, "Z+X", "l")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetNColumns(2);
    legend.SetTextSize(0.045);
    legend.Draw()

    mass.Draw("same")

    if (obsName=="pT4l"):
        label="p_{T}^{H}"
        unit="GeV"
    elif(obsName=="pT4lj"):
        label="p_{T}^{Hj}"
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
    elif (obsName=="njets_pt30_eta4p7"):
        label = "N(jets) |#eta|<4.7"
        unit = ""
    elif (obsName=='pt_leadingjet_pt30_eta4p7'):
        label = "p_{T}(jet)"
        unit = "GeV"
    elif (obsName=='pt_leadingjet_pt30_eta2p5'):
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
    elif (obsName=="costhetastar"):
        label = "|cos#theta*|"
        unit = ""
    elif (obsName=="costhetaZ1"):
        label = "|cos#theta_{1}|"
        unit = ""
    elif (obsName=="cosTheta2"):
        label = "|cos#theta_{2}|"
        unit = ""
    elif (obsName=="phi"):
        label = "|#Phi|"
        unit = ""
    elif (obsName=="Phi"):
        label = "|#Phi^{#star}|"
        unit = ""
    elif (obsName=="mass4l") or (obsName=='mass4l_zzfloating'):
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
#    latex2.DrawLatex(0.87, 0.95,"35.9 fb^{-1} at #sqrt{s} = 13 TeV")
#    latex2.DrawLatex(0.87, 0.95,"41.4 fb^{-1} at #sqrt{s} = 13 TeV") # 2017
#    latex2.DrawLatex(0.87, 0.95,"59.7 fb^{-1} at #sqrt{s} = 13 TeV") # 2017
#    latex2.DrawLatex(0.87, 0.95,"136.0 fb^{-1} at #sqrt{s} = 13 TeV") # 2017
    latex2.DrawLatex(0.95, 0.95,lumi+" fb^{-1} (#sqrt{s} = 13 TeV)") # 2017
    latex2.SetTextSize(0.8*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.16, 0.94, "CMS")
    latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    # latex2.DrawLatex(0.30, 0.95, "Preliminary")
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.45*c.GetTopMargin())
    # if (obsName!='mass4l' and obsName!='mass4l_zzfloating'): latex2.DrawLatex(0.65,0.85, str(observableBins[recobin])+" "+unit+" < "+label+" < "+str(observableBins[recobin+1])+" "+unit)

    checkDir("plots")
    checkDir("plots/"+obsName)

    if (not opt.UNBLIND):
        checkDir("plots/"+obsName+"/asimov")
        checkDir("plots/"+obsName+"/asimov/model")
        c.SaveAs("plots/"+obsName+"/asimov/model/asimovdata_"+physicalModel+"_"+opt.YEAR+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".pdf")
        c.SaveAs("plots/"+obsName+"/asimov/model/asimovdata_"+physicalModel+"_"+opt.YEAR+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".png")
        c.SaveAs("plots/"+obsName+"/asimov/model/asimovdata_"+physicalModel+"_"+opt.YEAR+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".root")
    else:
        checkDir("plots/"+obsName+"/data")
        checkDir("plots/"+obsName+"/data/model")
        c.SaveAs("plots/"+obsName+"/data/model/data_unfoldwith_"+modelName+"_"+physicalModel+"_"+opt.YEAR+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".pdf")
        c.SaveAs("plots/"+obsName+"/data/model/data_unfoldwith_"+modelName+"_"+physicalModel+"_"+opt.YEAR+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".png")
        c.SaveAs("plots/"+obsName+"/data/model/data_unfoldwith_"+modelName+"_"+physicalModel+"_"+opt.YEAR+"_"+obsName+'_'+fstate+"_recobin"+str(recobin)+".root")




######### ---------------------------------------- #########
######### ----------------- Main ----------------- #########
######### ---------------------------------------- #########
modelName = opt.UNFOLD

if 'vs' in opt.OBSNAME:
    obsName_tmp = opt.OBSNAME.split(' vs ')
    obsName = obsName_tmp[0]+'_'+obsName_tmp[1]
    doubleDiff = True
else:
    obsName = opt.OBSNAME
    doubleDiff = False

sys.path.append("inputs")
_temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
observableBins = _temp.observableBins
print observableBins
sys.path.remove("inputs")

if obsName.startswith("mass4l"):
    PhysicalModels = ['v3']
elif obsName == 'D0m' or obsName == 'Dcp' or obsName == 'D0hp' or obsName == 'Dint' or obsName == 'DL1' or obsName == 'DL1Zg' or obsName == 'costhetaZ1' or obsName == 'costhetaZ2'or obsName == 'costhetastar' or obsName == 'phi' or obsName == 'phistar' or obsName == 'massZ1' or obsName == 'massZ2':
    PhysicalModels = ['v4','v3']
elif 'kL' in obsName:
    PhysicalModels = ['kLambda']
else:
    PhysicalModels = ['v3']


nBins = len(observableBins)
if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries
print nBins
fStates = ["4e","4mu","2e2mu","4l"]
for fState in fStates:
    for recobin in range(nBins):
        for physicalModel in PhysicalModels:
            plotAsimov_sim(opt.UNFOLD, physicalModel, obsName, fState, observableBins, recobin)
