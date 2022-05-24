# this script is called once for each reco bin (obsBin)
# in each reco bin there are (nBins) signals (one for each gen bin)

import ROOT
import os,sys,subprocess
from math import trunc

sys.path.append('../helperstuff')

from paths import path

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
'rapidity4l_pT4l': True,
'njets_pt30_eta4p7 vs pT4l': False,
'pTj1_pTj2': False,
'pT4l_pTHj': False,
'massZ1_massZ2': False,
'TCjmax_pT4l': False
}


sys.path.append('../../inputs/')
sys.path.append('../../templates/')

def createXSworkspace(obsName, channel, nBins, obsBin, observableBins, addfakeH, modelName, physicalModel, year, JES, doubleDiff, lowerBound, upperBound, rawObsName):
    print '\n'
    print 'Creating WorkSpace', year

    if not doubleDiff:
        obsBin_low = observableBins[obsBin]
        obsBin_high = observableBins[obsBin+1]
        obs_bin_lowest = observableBins[0]
        obs_bin_highest = observableBins[len(observableBins)-1]
    elif doubleDiff:
        obsBin_low = observableBins[obsBin][0]
        obsBin_high = observableBins[obsBin][1]
        obs_bin_lowest = min(x[0] for x in observableBins.values())
        obs_bin_highest = max(x[1] for x in observableBins.values())
        #second variables
        obsBin_2nd_low = observableBins[obsBin][2]
        obsBin_2nd_high = observableBins[obsBin][3]
        obs_bin_2nd_lowest = min(x[2] for x in observableBins.values())
        obs_bin_2nd_highest = max(x[3] for x in observableBins.values())

    recobin = "recobin"+str(obsBin)

    print recobin
    JES = False
    doJES = False

    # Load some libraries
    ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
    ROOT.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so") #Print 0 in case of succesfull loading
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
    ROOT.gSystem.AddIncludePath("-Iinclude/")

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    # from inputs_sig imports some coefficients
    _temp = __import__('inputs_sig_extrap_'+obsName+'_'+year, globals(), locals(), ['acc','eff','inc_wrongfrac','binfrac_wrongfrac','outinratio'], -1)
    acc = _temp.acc
    eff = _temp.eff
    outinratio = _temp.outinratio
    # lambdajesup = _temp.lambdajesup
    # lambdajesdn = _temp.lambdajesdn
    inc_wrongfrac = _temp.inc_wrongfrac
    binfrac_wrongfrac = _temp.binfrac_wrongfrac
    #number_fake = _temp.number_fake

    # import h4l xs br
    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs4l_br'], -1)
    higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br

    _temp = __import__('inputs_bkg_'+obsName+'_'+year, globals(), locals(), ['fractionsBackground'], -1)
    fractionsBackground = _temp.fractionsBackground

    _temp = __import__('observables', globals(), locals(), ['observables'], -1)
    observables = _temp.observables

    if (not obsName=="mass4l"):
        obsName_help = observables[rawObsName]['obs_reco']
        if doubleDiff: obsName_2nd_help = observables[rawObsName]['obs_reco_2nd']

        observable = ROOT.RooRealVar(obsName_help,obsName_help,float(obs_bin_lowest),float(obs_bin_highest))
        if doubleDiff: observable_2nd = ROOT.RooRealVar(obsName_2nd_help,obsName_2nd_help,float(obs_bin_2nd_lowest),float(obs_bin_2nd_highest)) #ATdouble

        observable.Print()
        if doubleDiff: observable_2nd.Print()

    if doubleDiff:
        _recobin = str(observableBins[obsBin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin][3]).replace('.', 'p').replace('-','m')
    else:
        _recobin = str(observableBins[obsBin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin+1]).replace('.', 'p').replace('-','m')
        if int(observableBins[obsBin+1]) > 1000:
            _recobin = 'GT'+str(int(observableBins[obsBin]))

    _obsName = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta2p5': 'NJ'}
    if obsName not in _obsName:
        _obsName[obsName] = obsName

    # Parameters of doubleCB signal
    m = ROOT.RooRealVar("CMS_zz4l_mass", "CMS_zz4l_mass", lowerBound, upperBound)
    #MH = ROOT.RooRealVar("MH","MH",125.7,109.55,1000.05)
    MH = ROOT.RooRealVar("MH","MH", 125.38, lowerBound, upperBound)

    if( physicalModel=='v3'):
        comb_name = 'nonResH'
        sig_name = 'smH'
    else:
        comb_name = 'fakeH'
        sig_name = 'trueH'

    if(year == '2018'):
        CMS_zz4l_mean_m_sig_2018 = ROOT.RooRealVar("CMS_zz4l_mean_m_sig_2018","CMS_zz4l_mean_m_sig_2018",-10,10)
        CMS_zz4l_mean_e_sig_2018 = ROOT.RooRealVar("CMS_zz4l_mean_e_sig_2018","CMS_zz4l_mean_e_sig_2018",-10,10)
        CMS_zz4l_sigma_m_sig_2018 = ROOT.RooRealVar("CMS_zz4l_sigma_m_sig_2018","CMS_zz4l_sigma_m_sig_2018",-10,10)
        CMS_zz4l_sigma_e_sig_2018 = ROOT.RooRealVar("CMS_zz4l_sigma_e_sig_2018","CMS_zz4l_sigma_e_sig_2018",-10,10)
        lumi = ROOT.RooRealVar("lumi_132018","lumi_132018", 59.7)
        if (channel=='2e2mu'):
            CMS_zz4l_mean_m_err_3_2018 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_3_2018","CMS_zz4l_mean_m_err_3_2018",0.0004,0.0004,0.0004)
            CMS_zz4l_mean_e_err_3_2018 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_3_2018","CMS_zz4l_mean_e_err_3_2018",0.003,0.003,0.003)
            CMS_zz4l_n_sig_3_2018 = ROOT.RooRealVar("CMS_zz4l_n_sig_3_2018","CMS_zz4l_n_sig_3_2018",-10,10)
            CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2018", "(124.512837+0.99356*(@0-125)) + (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2018,CMS_zz4l_mean_e_sig_2018,CMS_zz4l_mean_m_err_3_2018,CMS_zz4l_mean_e_err_3_2018))
            CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018","CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018))
            CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2018", "(1.706305+0.0096*(@0-125))*(TMath::Sqrt((1+@1)*(1+@2)))",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2018,CMS_zz4l_sigma_e_sig_2018))
            CMS_zz4l_alpha_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2018","(1.00105+0.00414*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2018","(3.280255+(-0.06493)*(MH-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_3_2018))
            CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2018","(1.784923+(-0.00265)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2018","(3.590384+0.0087*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha_3_centralValue_2e2murecobin2018, CMS_zz4l_n_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018, CMS_zz4l_n2_3_centralValue_2e2murecobin2018)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha_3_centralValue_2e2murecobin2018, CMS_zz4l_n_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018, CMS_zz4l_n2_3_centralValue_2e2murecobin2018)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha_3_centralValue_2e2murecobin2018, CMS_zz4l_n_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018, CMS_zz4l_n2_3_centralValue_2e2murecobin2018)
        if (channel=='4e'):
            CMS_zz4l_mean_e_err_2_2018 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_2_2018","CMS_zz4l_mean_e_err_2_2018",0.003,0.003,0.003)
            CMS_zz4l_n_sig_2_2018 = ROOT.RooRealVar("CMS_zz4l_n_sig_2_2018","CMS_zz4l_n_sig_2_2018",-10,10)
            CMS_zz4l_mean_sig_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2018","CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2018", "(124.120149+0.99301*(@0-125)) + @0*@1*@2", ROOT.RooArgList(MH,CMS_zz4l_mean_e_sig_2018,CMS_zz4l_mean_e_err_2_2018))
            CMS_zz4l_mean_sig_NoConv_2_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_2_4erecobin2018","CMS_zz4l_mean_sig_NoConv_2_4erecobin2018","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_2_centralValue_4erecobin2018))
            CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2018","CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2018", "(2.242479+0.00793*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_e_sig_2018))
            CMS_zz4l_alpha_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2018","CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2018","(1.000023+0.00224*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n_2_centralValue_4e" + recobin + "2018","CMS_zz4l_n_2_centralValue_4e" + recobin + "2018","(3.84493+(-0.0734)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_2_2018))
            CMS_zz4l_alpha2_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2018","CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2018","(1.924843+(-0.01309)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n2_2_centralValue_4e" + recobin + "2018","CMS_zz4l_n2_2_centralValue_4e" + recobin + "2018","(3.670917+0.09786*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2018, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018, CMS_zz4l_alpha_2_centralValue_4erecobin2018, CMS_zz4l_n_2_centralValue_4erecobin2018, CMS_zz4l_alpha2_2_centralValue_4erecobin2018, CMS_zz4l_n2_2_centralValue_4erecobin2018)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2018, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018, CMS_zz4l_alpha_2_centralValue_4erecobin2018, CMS_zz4l_n_2_centralValue_4erecobin2018, CMS_zz4l_alpha2_2_centralValue_4erecobin2018, CMS_zz4l_n2_2_centralValue_4erecobin2018)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2018, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018, CMS_zz4l_alpha_2_centralValue_4erecobin2018, CMS_zz4l_n_2_centralValue_4erecobin2018, CMS_zz4l_alpha2_2_centralValue_4erecobin2018, CMS_zz4l_n2_2_centralValue_4erecobin2018)
        if (channel=='4mu'):
            CMS_zz4l_mean_m_err_1_2018 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_1_2018","CMS_zz4l_mean_m_err_1_2018",0.0004,0.0004,0.0004)
            CMS_zz4l_n_sig_1_2018 = ROOT.RooRealVar("CMS_zz4l_n_sig_1_2018","CMS_zz4l_n_sig_1_2018",-10,10)
            CMS_zz4l_mean_sig_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2018", "(124.842384+0.99991*(@0-125)) + @0*@1*@2",ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2018,CMS_zz4l_mean_m_err_1_2018))
            CMS_zz4l_mean_sig_NoConv_1_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_1_4murecobin2018","CMS_zz4l_mean_sig_NoConv_1_4murecobin2018","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_1_centralValue_4murecobin2018))
            CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2018", "(1.160456+0.00942*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2018))
            CMS_zz4l_alpha_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2018","(1.265095+(-0.0023)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_n_1_centralValue_4mu" + recobin + "2018","(2.048659+0.00083*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_1_2018))
            CMS_zz4l_alpha2_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2018","(1.94597+(-0.00139)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2018","(2.529541+0.02094*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2018, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018, CMS_zz4l_alpha_1_centralValue_4murecobin2018, CMS_zz4l_n_1_centralValue_4murecobin2018, CMS_zz4l_alpha2_1_centralValue_4murecobin2018, CMS_zz4l_n2_1_centralValue_4murecobin2018)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2018, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018, CMS_zz4l_alpha_1_centralValue_4murecobin2018, CMS_zz4l_n_1_centralValue_4murecobin2018, CMS_zz4l_alpha2_1_centralValue_4murecobin2018, CMS_zz4l_n2_1_centralValue_4murecobin2018)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2018, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018, CMS_zz4l_alpha_1_centralValue_4murecobin2018, CMS_zz4l_n_1_centralValue_4murecobin2018, CMS_zz4l_alpha2_1_centralValue_4murecobin2018, CMS_zz4l_n2_1_centralValue_4murecobin2018)

        # Wrong signal combination events
        if (channel=='4mu'):
            p1_12018 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_12018", "p1_12018","137.63080525+0.20718*(@0-125)",ROOT.RooArgList(MH))
            p2_12018 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_12018","p2_12018","17.73939875+(-0.0113475)*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_12018, p2_12018)
        if (channel=='4e'):
            p1_22018 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_22018","p1_22018","136.9123585+-0.020945*(@0-125)",ROOT.RooArgList(MH))
            p2_22018 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_22018","p2_22018","16.10426275+0.080675*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_22018, p2_22018)
        if (channel=='2e2mu'):
            p1_32018 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_32018","p1_32018","139.11565075+0.2334875*(@0-125)",ROOT.RooArgList(MH))
            p2_32018 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_32018","p2_32018","17.39454325+0.27387*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_32018, p2_32018)

    elif(year == '2017'):
        CMS_zz4l_mean_m_sig_2017 = ROOT.RooRealVar("CMS_zz4l_mean_m_sig_2017","CMS_zz4l_mean_m_sig_2017",-10,10)
        CMS_zz4l_mean_e_sig_2017 = ROOT.RooRealVar("CMS_zz4l_mean_e_sig_2017","CMS_zz4l_mean_e_sig_2017",-10,10)
        CMS_zz4l_sigma_m_sig_2017 = ROOT.RooRealVar("CMS_zz4l_sigma_m_sig_2017","CMS_zz4l_sigma_m_sig_2017",-10,10)
        CMS_zz4l_sigma_e_sig_2017 = ROOT.RooRealVar("CMS_zz4l_sigma_e_sig_2017","CMS_zz4l_sigma_e_sig_2017",-10,10)
        lumi = ROOT.RooRealVar("lumi_132017","lumi_132017", 41.5)
        if (channel=='2e2mu'):
            CMS_zz4l_mean_m_err_3_2017 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_3_2017","CMS_zz4l_mean_m_err_3_2017",0.0004,0.0004,0.0004)
            CMS_zz4l_mean_e_err_3_2017 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_3_2017","CMS_zz4l_mean_e_err_3_2017",0.003,0.003,0.003)
            CMS_zz4l_n_sig_3_2017 = ROOT.RooRealVar("CMS_zz4l_n_sig_3_2017","CMS_zz4l_n_sig_3_2017",-10,10)
            CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2017", "(124.498988+0.99545*(@0-125)) + (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2017,CMS_zz4l_mean_e_sig_2017,CMS_zz4l_mean_m_err_3_2017,CMS_zz4l_mean_e_err_3_2017))
            CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017","CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017))
            CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2017", "(1.679491+0.0043*(@0-125))*(TMath::Sqrt((1+@1)*(1+@2)))",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2017,CMS_zz4l_sigma_e_sig_2017))
            CMS_zz4l_alpha_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2017","(1.000075+(-0.00322)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2017","(3.081563+(-0.01503)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_3_2017))
            CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2017","(1.770038+(-0.01497)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2017","(3.645579+0.04511*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha_3_centralValue_2e2murecobin2017, CMS_zz4l_n_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017, CMS_zz4l_n2_3_centralValue_2e2murecobin2017)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha_3_centralValue_2e2murecobin2017, CMS_zz4l_n_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017, CMS_zz4l_n2_3_centralValue_2e2murecobin2017)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha_3_centralValue_2e2murecobin2017, CMS_zz4l_n_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017, CMS_zz4l_n2_3_centralValue_2e2murecobin2017)

        if (channel=='4e'):
            CMS_zz4l_mean_e_err_2_2017 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_2_2017","CMS_zz4l_mean_e_err_2_2017",0.003,0.003,0.003)
            CMS_zz4l_n_sig_2_2017 = ROOT.RooRealVar("CMS_zz4l_n_sig_2_2017","CMS_zz4l_n_sig_2_2017",-10,10)
            CMS_zz4l_mean_sig_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2017","CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2017", "(124.031616+0.99862*(@0-125)) + @0*@1*@2", ROOT.RooArgList(MH,CMS_zz4l_mean_e_sig_2017,CMS_zz4l_mean_e_err_2_2017))
            CMS_zz4l_mean_sig_NoConv_2_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_2_4erecobin2017","CMS_zz4l_mean_sig_NoConv_2_4erecobin2017","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_2_centralValue_4erecobin2017))
            CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2017","CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2017", "(2.243124+0.00452*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_e_sig_2017))
            CMS_zz4l_alpha_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2017","CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2017","(1.0+(-0.01389)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n_2_centralValue_4e" + recobin + "2017","CMS_zz4l_n_2_centralValue_4e" + recobin + "2017","(3.675354+0.09118*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_2_2017))
            CMS_zz4l_alpha2_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2017","CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2017","(1.927716+(-0.01581)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n2_2_centralValue_4e" + recobin + "2017","CMS_zz4l_n2_2_centralValue_4e" + recobin + "2017","(3.907514+0.09642*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2017, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017, CMS_zz4l_alpha_2_centralValue_4erecobin2017, CMS_zz4l_n_2_centralValue_4erecobin2017, CMS_zz4l_alpha2_2_centralValue_4erecobin2017, CMS_zz4l_n2_2_centralValue_4erecobin2017)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2017, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017, CMS_zz4l_alpha_2_centralValue_4erecobin2017, CMS_zz4l_n_2_centralValue_4erecobin2017, CMS_zz4l_alpha2_2_centralValue_4erecobin2017, CMS_zz4l_n2_2_centralValue_4erecobin2017)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2017, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017, CMS_zz4l_alpha_2_centralValue_4erecobin2017, CMS_zz4l_n_2_centralValue_4erecobin2017, CMS_zz4l_alpha2_2_centralValue_4erecobin2017, CMS_zz4l_n2_2_centralValue_4erecobin2017)

        if (channel=='4mu'):
            CMS_zz4l_mean_m_err_1_2017 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_1_2017","CMS_zz4l_mean_m_err_1_2017",0.0004,0.0004,0.0004)
            CMS_zz4l_n_sig_1_2017 = ROOT.RooRealVar("CMS_zz4l_n_sig_1_2017","CMS_zz4l_n_sig_1_2017",-10,10)
            CMS_zz4l_mean_sig_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2017", "(124.84515+0.9947*(@0-125)) + @0*@1*@2",ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2017,CMS_zz4l_mean_m_err_1_2017))
            CMS_zz4l_mean_sig_NoConv_1_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_1_4murecobin2017","CMS_zz4l_mean_sig_NoConv_1_4murecobin2017","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_1_centralValue_4murecobin2017))
            CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2017", "(1.16062+0.00777*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2017))
            CMS_zz4l_alpha_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2017","(1.2673+0.00332*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_n_1_centralValue_4mu" + recobin + "2017","(2.043382+(-0.01556)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_1_2017))
            CMS_zz4l_alpha2_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2017","(1.970747+(-0.00004)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2017","(2.454519+0.0031*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2017, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017, CMS_zz4l_alpha_1_centralValue_4murecobin2017, CMS_zz4l_n_1_centralValue_4murecobin2017, CMS_zz4l_alpha2_1_centralValue_4murecobin2017, CMS_zz4l_n2_1_centralValue_4murecobin2017)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2017, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017, CMS_zz4l_alpha_1_centralValue_4murecobin2017, CMS_zz4l_n_1_centralValue_4murecobin2017, CMS_zz4l_alpha2_1_centralValue_4murecobin2017, CMS_zz4l_n2_1_centralValue_4murecobin2017)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2017, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017, CMS_zz4l_alpha_1_centralValue_4murecobin2017, CMS_zz4l_n_1_centralValue_4murecobin2017, CMS_zz4l_alpha2_1_centralValue_4murecobin2017, CMS_zz4l_n2_1_centralValue_4murecobin2017)

        # Wrong signal combination events
        if (channel=='4mu'):
            p1_12017 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_12017", "p1_12017","137.5872+-0.0187475*(@0-125)",ROOT.RooArgList(MH))
            p2_12017 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_12017","p2_12017","18.05642225+-0.182495*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_12017, p2_12017)
        if (channel=='4e'):
            p1_22017 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_22017","p1_22017","133.98211375+0.790945*(@0-125)",ROOT.RooArgList(MH))
            p2_22017 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_22017","p2_22017","16.82234175+-0.2156775*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_22017, p2_22017)
        if (channel=='2e2mu'):
            p1_32017 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_32017","p1_32017","139.4166885+0.149575*(@0-125)",ROOT.RooArgList(MH))
            p2_32017 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_32017","p2_32017","17.0934935+0.2527325*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_32017, p2_32017)

    elif(year == '2016'):
        CMS_zz4l_mean_m_sig_2016 = ROOT.RooRealVar("CMS_zz4l_mean_m_sig_2016","CMS_zz4l_mean_m_sig_2016",-10,10)
        CMS_zz4l_mean_e_sig_2016 = ROOT.RooRealVar("CMS_zz4l_mean_e_sig_2016","CMS_zz4l_mean_e_sig_2016",-10,10)
        CMS_zz4l_sigma_m_sig_2016 = ROOT.RooRealVar("CMS_zz4l_sigma_m_sig_2016","CMS_zz4l_sigma_m_sig_2016",-10,10)
        CMS_zz4l_sigma_e_sig_2016 = ROOT.RooRealVar("CMS_zz4l_sigma_e_sig_2016","CMS_zz4l_sigma_e_sig_2016",-10,10)
        lumi = ROOT.RooRealVar("lumi_132016","lumi_132016", 35.9)
        if (channel=='2e2mu'):
            CMS_zz4l_mean_m_err_3_2016 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_3_2016","CMS_zz4l_mean_m_err_3_2016",0.0004,0.0004,0.0004)
            CMS_zz4l_mean_e_err_3_2016 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_3_2016","CMS_zz4l_mean_e_err_3_2016",0.003,0.003,0.003)
            CMS_zz4l_n_sig_3_2016 = ROOT.RooRealVar("CMS_zz4l_n_sig_3_2016","CMS_zz4l_n_sig_3_2016",-10,10)
            CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2016", "(124.510306+1.00214*(@0-125)) + (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2016,CMS_zz4l_mean_e_sig_2016,CMS_zz4l_mean_m_err_3_2016,CMS_zz4l_mean_e_err_3_2016))
            CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016","CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016))
            CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2016", "(1.656768+0.00585*(MH-125))*(TMath::Sqrt((1+@1)*(1+@2)))",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2016,CMS_zz4l_sigma_e_sig_2016))
            CMS_zz4l_alpha_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2016","(1.000125+(-0.0131)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2016","(2.959593+0.05254*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_3_2016))
            CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2016","(1.733022+0.00459*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2016","(3.801617+(-0.02044)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha_3_centralValue_2e2murecobin2016, CMS_zz4l_n_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016, CMS_zz4l_n2_3_centralValue_2e2murecobin2016)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha_3_centralValue_2e2murecobin2016, CMS_zz4l_n_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016, CMS_zz4l_n2_3_centralValue_2e2murecobin2016)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha_3_centralValue_2e2murecobin2016, CMS_zz4l_n_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016, CMS_zz4l_n2_3_centralValue_2e2murecobin2016)
        if (channel=='4e'):
            CMS_zz4l_mean_e_err_2_2016 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_2_2016","CMS_zz4l_mean_e_err_2_2016",0.003,0.003,0.003)
            CMS_zz4l_n_sig_2_2016 = ROOT.RooRealVar("CMS_zz4l_n_sig_2_2016","CMS_zz4l_n_sig_2_2016",-10,10)
            CMS_zz4l_mean_sig_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2016","CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2016", "(123.997243+1.00921*(@0-125)) + @0*@1*@2", ROOT.RooArgList(MH,CMS_zz4l_mean_e_sig_2016,CMS_zz4l_mean_e_err_2_2016))
            CMS_zz4l_mean_sig_NoConv_2_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_2_4erecobin2016","CMS_zz4l_mean_sig_NoConv_2_4erecobin2016","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_2_centralValue_4erecobin2016))
            CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2016","CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2016", "(2.236925+(-0.00431)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_e_sig_2016))
            CMS_zz4l_alpha_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2016","CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2016","(1.000197+(-0.02628)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n_2_centralValue_4e" + recobin + "2016","CMS_zz4l_n_2_centralValue_4e" + recobin + "2016","(3.445092+0.18446*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_2_2016))
            CMS_zz4l_alpha2_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2016","CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2016","(1.915296+(-0.02106)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n2_2_centralValue_4e" + recobin + "2016","CMS_zz4l_n2_2_centralValue_4e" + recobin + "2016","(3.924614+0.0656*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2016, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016, CMS_zz4l_alpha_2_centralValue_4erecobin2016, CMS_zz4l_n_2_centralValue_4erecobin2016, CMS_zz4l_alpha2_2_centralValue_4erecobin2016, CMS_zz4l_n2_2_centralValue_4erecobin2016)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2016, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016, CMS_zz4l_alpha_2_centralValue_4erecobin2016, CMS_zz4l_n_2_centralValue_4erecobin2016, CMS_zz4l_alpha2_2_centralValue_4erecobin2016, CMS_zz4l_n2_2_centralValue_4erecobin2016)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2016, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016, CMS_zz4l_alpha_2_centralValue_4erecobin2016, CMS_zz4l_n_2_centralValue_4erecobin2016, CMS_zz4l_alpha2_2_centralValue_4erecobin2016, CMS_zz4l_n2_2_centralValue_4erecobin2016)
        if (channel=='4mu'):
            CMS_zz4l_mean_m_err_1_2016 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_1_2016","CMS_zz4l_mean_m_err_1_2016",0.0004,0.0004,0.0004)
            CMS_zz4l_n_sig_1_2016 = ROOT.RooRealVar("CMS_zz4l_n_sig_1_2016","CMS_zz4l_n_sig_1_2016",-10,10)
            CMS_zz4l_mean_sig_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2016", "(124.844749+0.99461*(@0-125)) + @0*@1*@2",ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2016,CMS_zz4l_mean_m_err_1_2016))
            CMS_zz4l_mean_sig_NoConv_1_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_1_4murecobin2016","CMS_zz4l_mean_sig_NoConv_1_4murecobin2016","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_1_centralValue_4murecobin2016))
            CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2016", "(1.154947+0.01331*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2016))
            CMS_zz4l_alpha_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2016","(1.2625+0.00496*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_n_1_centralValue_4mu" + recobin + "2016","(2.045535+(-0.01482)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_1_2016))
            CMS_zz4l_alpha2_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2016","(1.870032+0.01029*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2016","(2.681714+(-0.01777)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB(sig_name+'_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2016, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016, CMS_zz4l_alpha_1_centralValue_4murecobin2016, CMS_zz4l_n_1_centralValue_4murecobin2016, CMS_zz4l_alpha2_1_centralValue_4murecobin2016, CMS_zz4l_n2_1_centralValue_4murecobin2016)
            ggH = ROOT.RooDoubleCB('ggH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2016, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016, CMS_zz4l_alpha_1_centralValue_4murecobin2016, CMS_zz4l_n_1_centralValue_4murecobin2016, CMS_zz4l_alpha2_1_centralValue_4murecobin2016, CMS_zz4l_n2_1_centralValue_4murecobin2016)
            xH = ROOT.RooDoubleCB('xH_' + _obsName[obsName], "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2016, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016, CMS_zz4l_alpha_1_centralValue_4murecobin2016, CMS_zz4l_n_1_centralValue_4murecobin2016, CMS_zz4l_alpha2_1_centralValue_4murecobin2016, CMS_zz4l_n2_1_centralValue_4murecobin2016)

        # Wrong signal combination events
        if (channel=='4mu'):
            p1_12016 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_12016", "p1_12016","136.7719185+0.3965175*(@0-125)",ROOT.RooArgList(MH))
            p2_12016 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_12016","p2_12016","17.86584175+-0.0833375*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_12016, p2_12016)
        if (channel=='4e'):
            p1_22016 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_22016","p1_22016","135.6637315+1.07089*(@0-125)",ROOT.RooArgList(MH))
            p2_22016 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_22016","p2_22016","15.976523+0.3065*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_22016, p2_22016)
        if (channel=='2e2mu'):
            p1_32016 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p1_32016","p1_32016","139.007959+0.569655*(@0-125)",ROOT.RooArgList(MH))
            p2_32016 = ROOT.RooFormulaVar("CMS_"+comb_name+"_p2_32016","p2_32016","18.213078+0.378485*(@0-125)",ROOT.RooArgList(MH))
            nonResH = ROOT.RooLandau(comb_name, "landau", m, p1_32016, p2_32016)

    # Out of acceptance events
    # (same shape as in acceptance shape)
    OutsideAcceptance = trueH.Clone()
    if( physicalModel=='v3'):
        OutsideAcceptance.SetName("OutsideAcceptance")
    else:
        OutsideAcceptance.SetName("out_trueH")

    # Coefficients for wrong signal combination events
    if (addfakeH):
        inc_wrongfrac_ggH=inc_wrongfrac["ggH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
        inc_wrongfrac_qqH=inc_wrongfrac["VBFH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
        inc_wrongfrac_WH=inc_wrongfrac["WH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
        inc_wrongfrac_ZH=inc_wrongfrac["ZH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
        inc_wrongfrac_ttH=inc_wrongfrac["ttH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    else:
        inc_wrongfrac_ggH=0.0
        inc_wrongfrac_qqH=0.0
        inc_wrongfrac_WH=0.0
        inc_wrongfrac_ZH=0.0
        inc_wrongfrac_ttH=0.0

    binfrac_wrongfrac_ggH=binfrac_wrongfrac["ggH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    binfrac_wrongfrac_qqH=binfrac_wrongfrac["VBFH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    binfrac_wrongfrac_WH=binfrac_wrongfrac["WH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    binfrac_wrongfrac_ZH=binfrac_wrongfrac["ZH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    binfrac_wrongfrac_ttH=binfrac_wrongfrac["ttH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]

    if (channel=='4e'):
        n_nonResH = (0.24*inc_wrongfrac_WH*binfrac_wrongfrac_WH+0.20*inc_wrongfrac_ZH*binfrac_wrongfrac_ZH+0.10*inc_wrongfrac_ttH*binfrac_wrongfrac_ttH)
    if (channel=='4mu'):
        n_nonResH = (0.45*inc_wrongfrac_WH*binfrac_wrongfrac_WH+0.38*inc_wrongfrac_ZH*binfrac_wrongfrac_ZH+0.20*inc_wrongfrac_ttH*binfrac_wrongfrac_ttH)
    if (channel=='2e2mu'):
        n_nonResH = (0.57*inc_wrongfrac_WH*binfrac_wrongfrac_WH+0.51*inc_wrongfrac_ZH*binfrac_wrongfrac_ZH+0.25*inc_wrongfrac_ttH*binfrac_wrongfrac_ttH)

    #numberFake_WH = number_fake["WH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    #numberFake_ZH = number_fake["ZH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    #numberFake_ttH = number_fake["ttH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]

    #n_nonResH = numberFake_WH + numberFake_ZH + numberFake_ttH

    n_nonResH_var = ROOT.RooRealVar("n_"+comb_name+"_var_"+recobin+"_"+channel+"_"+year,"n_"+comb_name+"_var_"+recobin+"_"+channel+"_"+year,n_nonResH)
    nonResH_norm = ROOT.RooFormulaVar(comb_name+"_norm","@0",ROOT.RooArgList(n_nonResH_var))


    # signal shape in different recobin
    trueH_shape = {}; ggH_shape = {}; xH_shape = {}
    fideff = {}
    fideff_ggH = {}; fideff_xH = {}
    fideff_var = {}; fideff_ggH_var = {}; fideff_xH_var = {}
    trueH_norm = {}; GGH_norm = {}; XH_norm = {}

    rBin_channel = {}

    trueH_norm_final = {}
    fracBin = {}
    fracSM4eBin = {}; fracGGH4eBin = {}; fracXH4eBin = {}
    fracSM4muBin = {}; fracGGH4muBin = {}; fracXH4muBin = {}
    K1Bin = {}; K1GGHBin = {}; K1XHBin = {}
    K2Bin = {}; K2GGHBin = {}; K2XHBin = {}
    SigmaBin = {}; SigmaBin_ggH = {}; SigmaBin_xH = {}
    SigmaHBin = {}; SigmaHBin_ggH = {}; SigmaHBin_xH = {}

    for genbin in range(nBins):

        # Set name of the process according to conventions for combination
        if not doubleDiff:
            if observableBins[genbin+1] > 1000:
                _binName = 'GT'+str(int(observableBins[genbin]))
            else:
                _binName = str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m')
        else:
            _binName = str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p').replace('-','m')


        if( physicalModel=='v3'):
            processName = sig_name+'_' + _obsName[obsName]
            processName = processName+'_'+_binName
            trueH_shape[genbin] = trueH.Clone();
            trueH_shape[genbin].SetName(processName)
        else:
            processName = "trueH"+channel+"Bin"+str(genbin)
            trueH_shape[genbin] = trueH.Clone();
            trueH_shape[genbin].SetName(processName)

        fideff[genbin] = eff[modelName+"_"+channel+"_"+obsName+"_genbin"+str(genbin)+"_"+recobin]
        print "fideff[genbin]", fideff[genbin]
        print "model name is ", modelName

        _effName = "eff_hzz_"+sig_name+"_"+_obsName[obsName]+"_13TeV_"+_binName
        _effName = _effName + '_hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel

        fideff_var[genbin] = ROOT.RooRealVar(_effName, _effName, fideff[genbin]);

    for genbin in range(nBins):

        # Set name of the process according to conventions for combination
        if not doubleDiff:
            if observableBins[genbin+1] > 1000:
                _binName = 'GT'+str(int(observableBins[genbin]))
            else:
                _binName = str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m')
        else:
            _binName = str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p').replace('-','m')

        if( physicalModel=='v3'):
            processName = sig_name+'_' + _obsName[obsName]
            processName = processName+'_'+_binName
        else:
            processName = "trueH"+channel+"Bin"+str(genbin)

        fidxs = {}; fidxs_ggH = {}; fidxs_xH = {}
        for fState in ['4e','4mu', '2e2mu']:
            fidxs[fState] = 0; fidxs_ggH[fState] = 0; fidxs_xH[fState] = 0

            fidxs[fState] += higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs[fState] += higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs[fState] += higgs_xs['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs[fState] += higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs[fState] += higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]

            fidxs_ggH[fState] += higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs_xH[fState] += higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs_xH[fState] += higgs_xs['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs_xH[fState] += higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fidxs_xH[fState] += higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]

        ggHName = 'ggH_' + _obsName[obsName]
        ggHName = ggHName+'_'+_binName
        ggH_shape[genbin] = ggH.Clone();
        ggH_shape[genbin].SetName(ggHName)
        fideff_ggH[genbin] = eff["ggH125_"+channel+"_"+obsName+"_genbin"+str(genbin)+"_"+recobin]

        xHName = 'xH_' + _obsName[obsName]
        xHName = xHName+'_'+_binName
        xH_shape[genbin] = xH.Clone();
        xH_shape[genbin].SetName(xHName)

        xheff = 0.0
        sumxsec = 0.0
        for prodMode in ['VBFH', 'ZH', 'WH', 'ttH']:
            _prodMode = prodMode
            if prodMode == 'VBFH':
                _prodMode = 'VBF'
            xheff += higgs_xs[_prodMode+'_125.0']*eff[prodMode+"125_"+channel+"_"+obsName+"_genbin"+str(genbin)+"_"+recobin]*acc[prodMode+'125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            sumxsec += higgs_xs[_prodMode+'_125.0']
        fideff_xH[genbin] = xheff/sumxsec
        xheff = 0.0

        _effName_ggH = "eff_hzz_GGH_"+_obsName[obsName]+"_13TeV_"+_binName
        _effName_ggH = _effName_ggH + '_hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel
        _effName_xH = "eff_hzz_XH_"+_obsName[obsName]+"_13TeV_"+_binName
        _effName_xH = _effName_xH + '_hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel

        fideff_ggH_var[genbin] = ROOT.RooRealVar(_effName_ggH, _effName_ggH, fideff_ggH[genbin]);
        fideff_xH_var[genbin] = ROOT.RooRealVar(_effName_xH, _effName_xH, fideff_xH[genbin]);

        if (physicalModel=='v2'):

            # In these models the xsec is left floting, being directly the POI
            trueH_norm[genbin] = ROOT.RooFormulaVar(processName+"_norm","@0*@1", ROOT.RooArgList(fideff_var[genbin], lumi) );
            rBin_channel[str(genbin)] = ROOT.RooRealVar("r"+channel+"Bin"+str(genbin),"r"+channel+"Bin"+str(genbin), 1.0, 0.0, 10.0)
            rBin_channel[str(genbin)].setConstant(True)
            trueH_norm_final[genbin] = ROOT.RooFormulaVar(processName+"_"+recobin+"_final","@0*@1*@2", ROOT.RooArgList(rBin_channel[str(genbin)], fideff_var[genbin],lumi))

        elif (physicalModel=='v4'):

            # In these models the xsec is left floting, being directly the POI
            trueH_norm[genbin] = ROOT.RooFormulaVar(processName+"_norm","@0*@1", ROOT.RooArgList(fideff_var[genbin], lumi) );
            if channel == '2e2mu':
                rBin_channel['2e2mu'+str(genbin)] = ROOT.RooRealVar("r"+channel+"Bin"+str(genbin),"r"+channel+"Bin"+str(genbin), 1.0, 0.0, 10.0)
                trueH_norm_final[genbin] = ROOT.RooFormulaVar(processName+"_"+recobin+"_final","@0*@1*@2", ROOT.RooArgList(rBin_channel['2e2mu'+str(genbin)], fideff_var[genbin],lumi))
            else: #4e+4mu
                fidxs = {}
                for fState in ['4e','4mu']:
                    fidxs[fState] = 0
                    fidxs[fState] += higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fidxs[fState] += higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fidxs[fState] += higgs_xs['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fidxs[fState] += higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fidxs[fState] += higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                fidxs['4l'] = fidxs['4e'] + fidxs['4mu']
                fracSM4eBin[str(genbin)] = ROOT.RooRealVar('fracSM4eBin'+str(genbin), 'fracSM4eBin'+str(genbin), fidxs['4e']/fidxs['4l'])
                fracSM4eBin[str(genbin)].setConstant(True)
                fracSM4muBin[str(genbin)] = ROOT.RooRealVar('fracSM4muBin'+str(genbin), 'fracSM4muBin'+str(genbin), fidxs['4mu']/fidxs['4l'])
                fracSM4muBin[str(genbin)].setConstant(True)

                rBin_channel[str(genbin)] = ROOT.RooRealVar("r4lBin"+str(genbin),"r4lBin"+str(genbin), fidxs['4l'] , 0.0, 10.0)

                rBin_channel['4e'+str(genbin)] = ROOT.RooFormulaVar("r4eBin"+str(genbin),"(@0*@1)", ROOT.RooArgList(rBin_channel[str(genbin)], fracSM4eBin[str(genbin)]))
                rBin_channel['4mu'+str(genbin)] = ROOT.RooFormulaVar("r4muBin"+str(genbin),"(@0*@1)", ROOT.RooArgList(rBin_channel[str(genbin)], fracSM4muBin[str(genbin)]))
                trueH_norm_final[genbin] = ROOT.RooFormulaVar(processName+"_"+recobin+"_final","@0*@1*@2", ROOT.RooArgList(rBin_channel[channel+str(genbin)], fideff_var[genbin],lumi))

        elif (physicalModel=='kLambda'):
            muBin = {}
            C1_map = {}
            tot_xs = {}
            fidxs_fl = {}
            C1_ggH = {}; mu_ggH = {}; scale_ggH = {}; fidxs_ggH = {}
            C1_VH = {}; mu_VH = {}; scale_VH = {}; fidxs_WH = {}; fidxs_ZH = {}
            C1_ttH = {}; mu_ttH = {}; scale_ttH = {}; fidxs_ttH = {}
            C1_VBFH = {}; mu_VBFH = {}; scale_VBFH = {}; fidxs_VBFH = {}
            C1_HZZ = {}; mu_BR = {}
            C1_tot = {}
            dZH = {}

            kappa_lambda = ROOT.RooRealVar("kappa_lambda", "kappa_lambda", 1.0, -10.0, 20.0)

            C1i_ttH = [0.0530525859571, 0.0472618825815, 0.0392337055167, 0.0278818345971, 0.0141882242091]
            C1i_VH = [0.0165863149378, 0.012328663897, 0.00774755197694, 0.0034957241269, 0.00024199147094]
            for genbin in range(nBins):
                C1_ggH[str(genbin)] = ROOT.RooRealVar("C1_ggH_"+str(genbin), "C1_ggH_"+str(genbin), 0.0066, 0.0066, 0.0066)

                C1_ttH[str(genbin)] = ROOT.RooRealVar("C1_ttH_"+str(genbin), "C1_ttH_"+str(genbin), C1i_ttH[genbin], C1i_ttH[genbin], C1i_ttH[genbin])

                C1_VH[str(genbin)] = ROOT.RooRealVar("C1_VH_"+str(genbin), "C1_VH_"+str(genbin), C1i_VH[genbin], C1i_VH[genbin], C1i_VH[genbin])

                C1_VBFH[str(genbin)] = ROOT.RooRealVar("C1_VBFH_"+str(genbin), "C1_VBFH_"+str(genbin), 0.0063, 0.0063, 0.0063)

                # hzz4l
                C1_HZZ[str(genbin)] = ROOT.RooRealVar("C1_HZZ_"+str(genbin), "C1_HZZ_"+str(genbin), 0.0083, 0.0083, 0.0083)

                C1_tot[str(genbin)] = ROOT.RooRealVar("C1_tot_"+str(genbin), "C1_tot", 2.5e-3, 2.5e-3, 2.5e-3)

                #Define dZH constant variable
                dZH[str(genbin)] = ROOT.RooRealVar("dZH_"+str(genbin), "dZH_"+str(genbin), -1.536e-3, -1.536e-3, -1.536e-3)

                mu_ggH[str(genbin)] = ROOT.RooFormulaVar("XSscal_ggH_"+str(genbin), "XSscal_ggH_"+str(genbin), "(1+@0*@1+@2)/((1-(@0*@0-1)*@2)*(1+@1+@2))", ROOT.RooArgList(kappa_lambda, C1_ggH[str(genbin)],dZH[str(genbin)]))
                mu_VBFH[str(genbin)] = ROOT.RooFormulaVar("XSscal_VBFH_"+str(genbin), "XSscal_VBFH_"+str(genbin), "(1+@0*@1+@2)/((1-(@0*@0-1)*@2)*(1+@1+@2))", ROOT.RooArgList(kappa_lambda, C1_VBFH[str(genbin)],dZH[str(genbin)]))
                mu_VH [str(genbin)]= ROOT.RooFormulaVar("XSscal_VH_"+str(genbin), "XSscal_VH_"+str(genbin), "(1+@0*@1+@2)/((1-(@0*@0-1)*@2)*(1+@1+@2))", ROOT.RooArgList(kappa_lambda, C1_VH[str(genbin)],dZH[str(genbin)]))
                mu_ttH[str(genbin)] = ROOT.RooFormulaVar("XSscal_ttH_"+str(genbin), "XSscal_ttH_"+str(genbin), "(1+@0*@1+@2)/((1-(@0*@0-1)*@2)*(1+@1+@2))", ROOT.RooArgList(kappa_lambda, C1_ttH[str(genbin)],dZH[str(genbin)]))

                mu_BR[str(genbin)] = ROOT.RooFormulaVar("BRscal_hzz_"+str(genbin), "BRscal_hzz_"+str(genbin), "(1+(((@0-1)*(@1-@2))/(1+(@0-1)*@2)))", ROOT.RooArgList(kappa_lambda, C1_HZZ[str(genbin)], C1_tot[str(genbin)]))

                scale_ggH[str(genbin)] = ROOT.RooFormulaVar("scale_ggH_"+str(genbin), "scale_ggH_"+str(genbin), "@0*@1", ROOT.RooArgList(mu_ggH[str(genbin)],mu_BR[str(genbin)]))
                scale_VBFH[str(genbin)] = ROOT.RooFormulaVar("scale_VBFH_"+str(genbin), "scale_VBFH_"+str(genbin), "@0*@1", ROOT.RooArgList(mu_VBFH[str(genbin)],mu_BR[str(genbin)]))
                scale_VH[str(genbin)] = ROOT.RooFormulaVar("scale_VH_"+str(genbin), "scale_VH_"+str(genbin),"@0*@1", ROOT.RooArgList(mu_VH[str(genbin)],mu_BR[str(genbin)]))
                scale_ttH[str(genbin)] = ROOT.RooFormulaVar("scale_ttH_"+str(genbin), "scale_ttH_"+str(genbin),"@0*@1", ROOT.RooArgList(mu_ttH[str(genbin)],mu_BR[str(genbin)]))

                fidxs_ggH[str(genbin)] = ROOT.RooRealVar("fidxs_ggH_"+str(genbin), "fidxs_ggH_"+str(genbin), 1.0*higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                fidxs_ggH[str(genbin)].setConstant(True)
                fidxs_VBFH[str(genbin)] = ROOT.RooRealVar("fidxs_VBFH_"+str(genbin), "fidxs_VBFH_"+str(genbin), 1.000*higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                fidxs_VBFH[str(genbin)].setConstant(True)
                fidxs_WH[str(genbin)] = ROOT.RooRealVar("fidxs_WH_"+str(genbin), "fidxs_WH_"+str(genbin), 1.000*higgs_xs['WH_125.0']*higgs4l_br['125.0_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                fidxs_WH[str(genbin)].setConstant(True)
                fidxs_ZH[str(genbin)] = ROOT.RooRealVar("fidxs_ZH_"+str(genbin), "fidxs_ZH_"+str(genbin), 1.000*higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                fidxs_ZH[str(genbin)].setConstant(True)
                fidxs_ttH[str(genbin)] = ROOT.RooRealVar("fidxs_ttH_"+str(genbin), "fidxs_ttH_"+str(genbin), 1.000*higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                fidxs_ttH[str(genbin)].setConstant(True)

                fidxs_fl[str(genbin)] = ROOT.RooFormulaVar("fidxs_fl_"+str(genbin), "fidxs_fl_"+str(genbin), "(@0*@5+@1*@6+(@2+@3)*@7+@4*@8)", ROOT.RooArgList(fidxs_ggH[str(genbin)], fidxs_VBFH[str(genbin)], fidxs_WH[str(genbin)], fidxs_ZH[str(genbin)], fidxs_ttH[str(genbin)], scale_ggH[str(genbin)], scale_VBFH[str(genbin)], scale_VH[str(genbin)], scale_ttH[str(genbin)]))

                print('ggH Bin', genbin, higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                print('VBF Bin', genbin, higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                print('WH Bin', genbin, higgs_xs['WH_125.0']*higgs4l_br['125.0_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                print('ZH Bin', genbin, higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])
                print('ttH Bin', genbin, higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)])

                # Set name of the process according to conventions for combination
                if observableBins[genbin+1] > 1000:
                  _binName = 'GT'+str(int(observableBins[genbin]))
                else:
                  _binName = str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m')

                # Need to be redefined here for proper import in ws
                if( physicalModel=='v3'):
                    processName = sig_name+'_' + _obsName[obsName]
                    processName = processName+'_'+_binName
                else:
                    processName = "trueH"+channel+"Bin"+str(genbin)

                # For klambda we need to define ad hoc these because of the C1 scalings
                trueH_norm_final[genbin] = ROOT.RooFormulaVar(processName+"_"+recobin+"_final","@0*@1*@2", ROOT.RooArgList(fideff_var[genbin], lumi, fidxs_fl[str(genbin)]));
                trueH_norm[genbin] = ROOT.RooFormulaVar(processName+"_norm","@0*@1*@2", ROOT.RooArgList(fideff_var[genbin], lumi, fidxs_fl[str(genbin)]));

        else:

            fidxs['4l'] = fidxs['4e'] + fidxs['4mu'] + fidxs['2e2mu']
            fidxs_ggH['4l'] = fidxs_ggH['4e'] + fidxs_ggH['4mu'] + fidxs_ggH['2e2mu']
            fidxs_xH['4l'] = fidxs_xH['4e'] + fidxs_xH['4mu'] + fidxs_xH['2e2mu']

            fracSM4eBin[str(genbin)] = ROOT.RooRealVar('fracSM4eBin'+str(genbin), 'fracSM4eBin'+str(genbin), fidxs['4e']/fidxs['4l'])
            fracSM4eBin[str(genbin)].setConstant(True)
            fracSM4muBin[str(genbin)] = ROOT.RooRealVar('fracSM4muBin'+str(genbin), 'fracSM4muBin'+str(genbin), fidxs['4mu']/fidxs['4l'])
            fracSM4muBin[str(genbin)].setConstant(True)

            fracGGH4eBin[str(genbin)] = ROOT.RooRealVar('fracGGH4eBin'+str(genbin), 'fracGGH4eBin'+str(genbin), fidxs_ggH['4e']/fidxs_ggH['4l'])
            fracGGH4eBin[str(genbin)].setConstant(True)
            fracGGH4muBin[str(genbin)] = ROOT.RooRealVar('fracGGH4muBin'+str(genbin), 'fracGGH4muBin'+str(genbin), fidxs_ggH['4mu']/fidxs_ggH['4l'])
            fracGGH4muBin[str(genbin)].setConstant(True)

            fracXH4eBin[str(genbin)] = ROOT.RooRealVar('fracXH4eBin'+str(genbin), 'fracXH4eBin'+str(genbin), fidxs_xH['4e']/fidxs_xH['4l'])
            fracXH4eBin[str(genbin)].setConstant(True)
            fracXH4muBin[str(genbin)] = ROOT.RooRealVar('fracXH4muBin'+str(genbin), 'fracXH4muBin'+str(genbin), fidxs_xH['4mu']/fidxs_xH['4l'])
            fracXH4muBin[str(genbin)].setConstant(True)

            K1Bin[str(genbin)] = ROOT.RooRealVar('K1Bin'+str(genbin), 'K1Bin'+str(genbin), 1.0, 0.0,  1.0/fracSM4eBin[str(genbin)].getVal())
            K2Bin[str(genbin)] = ROOT.RooRealVar('K2Bin'+str(genbin), 'K2Bin'+str(genbin), 1.0, 0.0, (1.0-fracSM4eBin[str(genbin)].getVal())/fracSM4muBin[str(genbin)].getVal())

            K1GGHBin[str(genbin)] = ROOT.RooRealVar('K1GGHBin'+str(genbin), 'K1GGHBin'+str(genbin), 1.0, 0.0,  1.0/fracGGH4eBin[str(genbin)].getVal())
            K2GGHBin[str(genbin)] = ROOT.RooRealVar('K2GGHBin'+str(genbin), 'K2GGHBin'+str(genbin), 1.0, 0.0, (1.0-fracGGH4eBin[str(genbin)].getVal())/fracGGH4muBin[str(genbin)].getVal())

            K1XHBin[str(genbin)] = ROOT.RooRealVar('K1XHBin'+str(genbin), 'K1XHBin'+str(genbin), 1.0, 0.0,  1.0/fracXH4eBin[str(genbin)].getVal())
            K2XHBin[str(genbin)] = ROOT.RooRealVar('K2XHBin'+str(genbin), 'K2XHBin'+str(genbin), 1.0, 0.0, (1.0-fracXH4eBin[str(genbin)].getVal())/fracXH4muBin[str(genbin)].getVal())

            if not doubleDiff:
                SigmaBin[str(genbin)] = ROOT.RooRealVar('xs_hzz_'+sig_name+'_'+_obsName[obsName]+'_'+str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m'), 'xs_hzz_'+sig_name+'_'+_obsName[obsName]+'_'+str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m'), fidxs['4l'], 0.0, 10.0)
            else:
                SigmaBin[str(genbin)] = ROOT.RooRealVar('xs_hzz_'+sig_name+'_'+_obsName[obsName]+'_'+str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p').replace('-','m'), 'xs_hzz_'+sig_name+'_'+_obsName[obsName]+'_'+str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1])+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]), fidxs['4l'], 0.0, 10.0)
            SigmaBin[str(genbin)].setConstant(True)
            SigmaHBin['4e'+str(genbin)] = ROOT.RooFormulaVar("Sigma4eBin"+str(genbin),"(@0*@1*@2)", ROOT.RooArgList(SigmaBin[str(genbin)], fracSM4eBin[str(genbin)], K1Bin[str(genbin)]))
            SigmaHBin['4mu'+str(genbin)] = ROOT.RooFormulaVar("Sigma4muBin"+str(genbin),"(@0*(1.0-@1*@2)*@3*@4/(1.0-@1))", ROOT.RooArgList(SigmaBin[str(genbin)], fracSM4eBin[str(genbin)], K1Bin[str(genbin)], K2Bin[str(genbin)], fracSM4muBin[str(genbin)]))
            SigmaHBin['2e2mu'+str(genbin)] = ROOT.RooFormulaVar("Sigma2e2muBin"+str(genbin),"(@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1)))", ROOT.RooArgList(SigmaBin[str(genbin)], fracSM4eBin[str(genbin)], K1Bin[str(genbin)], K2Bin[str(genbin)], fracSM4muBin[str(genbin)]))

            if not doubleDiff:
                SigmaBin_ggH[str(genbin)] = ROOT.RooRealVar('xs_hzz_ggH_'+_obsName[obsName]+'_'+str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m'), 'xs_hzz_ggH_'+_obsName[obsName]+'_'+str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m'), fidxs_ggH['4l'], 0.0, 10.0)
            else:
                SigmaBin_ggH[str(genbin)] = ROOT.RooRealVar('xs_hzz_ggH_'+_obsName[obsName]+'_'+str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p').replace('-','m'), 'xs_hzz_ggH_'+_obsName[obsName]+'_'+str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p'), fidxs_ggH['4l'], 0.0, 10.0)
            SigmaBin_ggH[str(genbin)].setConstant(True)
            SigmaHBin_ggH['4e'+str(genbin)] = ROOT.RooFormulaVar("SigmaGGH4eBin"+str(genbin),"(@0*@1*@2)", ROOT.RooArgList(SigmaBin_ggH[str(genbin)], fracGGH4eBin[str(genbin)], K1GGHBin[str(genbin)]))
            SigmaHBin_ggH['4mu'+str(genbin)] = ROOT.RooFormulaVar("SigmaGGH4muBin"+str(genbin),"(@0*(1.0-@1*@2)*@3*@4/(1.0-@1))", ROOT.RooArgList(SigmaBin_ggH[str(genbin)], fracGGH4eBin[str(genbin)], K1GGHBin[str(genbin)], K2GGHBin[str(genbin)], fracGGH4muBin[str(genbin)]))
            SigmaHBin_ggH['2e2mu'+str(genbin)] = ROOT.RooFormulaVar("SigmaGGH2e2muBin"+str(genbin),"(@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1)))", ROOT.RooArgList(SigmaBin_ggH[str(genbin)], fracGGH4eBin[str(genbin)], K1GGHBin[str(genbin)], K2GGHBin[str(genbin)], fracGGH4muBin[str(genbin)]))

            if not doubleDiff:
                SigmaBin_xH[str(genbin)] = ROOT.RooRealVar('xs_hzz_xH_'+_obsName[obsName]+'_'+str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m'), 'xs_hzz_xH_'+_obsName[obsName]+'_'+str(observableBins[genbin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin+1]).replace('.', 'p').replace('-','m'), fidxs_xH['4l'], 0.0, 10.0)
            else:
                SigmaBin_xH[str(genbin)] = ROOT.RooRealVar('xs_hzz_xH_'+_obsName[obsName]+'_'+str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p').replace('-','m'), 'xs_hzz_xH_'+_obsName[obsName]+'_'+str(observableBins[genbin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[genbin][3]).replace('.', 'p').replace('-','m'), fidxs_xH['4l'], 0.0, 10.0)
            SigmaBin_xH[str(genbin)].setConstant(True)
            SigmaHBin_xH['4e'+str(genbin)] = ROOT.RooFormulaVar("SigmaXH4eBin"+str(genbin),"(@0*@1*@2)", ROOT.RooArgList(SigmaBin_xH[str(genbin)], fracXH4eBin[str(genbin)], K1XHBin[str(genbin)]))
            SigmaHBin_xH['4mu'+str(genbin)] = ROOT.RooFormulaVar("SigmaXH4muBin"+str(genbin),"(@0*(1.0-@1*@2)*@3*@4/(1.0-@1))", ROOT.RooArgList(SigmaBin_xH[str(genbin)], fracXH4eBin[str(genbin)], K1XHBin[str(genbin)], K2XHBin[str(genbin)], fracXH4muBin[str(genbin)]))
            SigmaHBin_xH['2e2mu'+str(genbin)] = ROOT.RooFormulaVar("SigmaXH2e2muBin"+str(genbin),"(@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1)))", ROOT.RooArgList(SigmaBin_xH[str(genbin)], fracXH4eBin[str(genbin)], K1XHBin[str(genbin)], K2XHBin[str(genbin)], fracXH4muBin[str(genbin)]))

            # Here we define the scalings used for the Combination (and v3 measurements)
            # The POIs are now signal strenght modifiers, hence xsecSM is constant
            trueH_norm_final[genbin] = ROOT.RooFormulaVar(processName+"_"+recobin+"_final","@0*@1*@2" ,ROOT.RooArgList(SigmaHBin[channel+str(genbin)],fideff_var[genbin],lumi))
            trueH_norm[genbin] = ROOT.RooFormulaVar(processName+"_norm","@0*@1*@2", ROOT.RooArgList(SigmaHBin[channel+str(genbin)],fideff_var[genbin], lumi) );
            GGH_norm[genbin] = ROOT.RooFormulaVar(ggHName+"_norm","@0*@1*@2", ROOT.RooArgList(SigmaHBin_ggH[channel+str(genbin)],fideff_ggH_var[genbin], lumi) );
            XH_norm[genbin] = ROOT.RooFormulaVar(xHName+"_norm","@0*@1*@2", ROOT.RooArgList(SigmaHBin_xH[channel+str(genbin)],fideff_xH_var[genbin], lumi) );

    outin = outinratio[modelName+"_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    # print "outin",obsBin,outin
    outin_var = ROOT.RooRealVar("outfracBin_"+recobin+"_"+channel+year,"outfracBin_"+recobin+"_"+channel+year, outin);
    outin_var.setConstant(True)
    OutsideAcceptance_norm_args = ROOT.RooArgList(outin_var)
    OutsideAcceptance_norm_func = "@0*("
    for i in range(nBins):
        OutsideAcceptance_norm_args.add(trueH_norm_final[i])
        OutsideAcceptance_norm_func = OutsideAcceptance_norm_func+"@"+str(i+1)+"+"
    OutsideAcceptance_norm_func = OutsideAcceptance_norm_func.replace(str(nBins)+"+",str(nBins)+")")
    if( physicalModel=='v3'):
        OutsideAcceptance_norm = ROOT.RooFormulaVar("OutsideAcceptance_norm",OutsideAcceptance_norm_func,OutsideAcceptance_norm_args)
    else:
        OutsideAcceptance_norm = ROOT.RooFormulaVar("out_trueH_norm",OutsideAcceptance_norm_func,OutsideAcceptance_norm_args)

    frac_qqzz = fractionsBackground['qqzz_'+channel+'_'+obsName+'_'+recobin]
    frac_qqzz_var  = ROOT.RooRealVar("frac_qqzz_"+recobin+"_"+channel+"_"+year,"frac_qqzz_"+recobin+"_"+channel+"_"+year, frac_qqzz);

    frac_ggzz = fractionsBackground['ggzz_'+channel+'_'+obsName+'_'+recobin]
    frac_ggzz_var = ROOT.RooRealVar("frac_ggzz_"+recobin+"_"+channel+"_"+year,"frac_ggzz_"+recobin+"_"+channel+"_"+year, frac_ggzz);

    frac_zjets = fractionsBackground['ZJetsCR_'+channel+'_'+obsName+'_'+recobin]
    frac_zjets_var = ROOT.RooRealVar("frac_zjet_"+recobin+"_"+channel+"_"+year,"frac_zjet_"+recobin+"_"+channel+"_"+year, frac_zjets);

    print obsBin,"frac_qqzz",frac_qqzz,"frac_ggzz",frac_ggzz,"frac_zjets",frac_zjets

    os.chdir('../../templates/'+year+"/"+obsName+"/")

    if doubleDiff and decimal[obsName]:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high)+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high)+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high)+".root"
    elif doubleDiff:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high))+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high))+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high))+".root"
    elif decimal[obsName]:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
    elif not doubleDiff:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+".root"
    if obsName == 'rapidity4l':
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"

    # if (not obsName=="mass4l"):
    #     template_zjetsName = "/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/templates/"+year+"/"+obsName+"/XSBackground_ZJetsCR_AllChans_"+obsName+"_"+obsBin_low+"_"+obsBin_high+".root"
    # else:
    #     template_zjetsName = "/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/templates/"+year+"/"+obsName+"/XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+obsBin_low+"_"+obsBin_high+".root"

    qqzzTempFile = ROOT.TFile(template_qqzzName,"READ")
    if decimal[obsName] and doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high))
    elif decimal[obsName]: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    elif not doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high)))
    elif doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    # if obsName == 'rapidity4l': qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    # elif doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    print qqzzTempFile
    print qqzzTemplate.GetName()
    print 'qqZZ bins',qqzzTemplate.GetNbinsX(),qqzzTemplate.GetBinLowEdge(1),qqzzTemplate.GetBinLowEdge(qqzzTemplate.GetNbinsX()+1)

    ggzzTempFile = ROOT.TFile(template_ggzzName,"READ")
    if decimal[obsName] and doubleDiff: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high))
    elif decimal[obsName]: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    elif not doubleDiff: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high)))
    elif doubleDiff: ggzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    # if obsName == 'rapidity4l': ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    # elif doubleDiff: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    print 'ggZZ bins',ggzzTemplate.GetNbinsX(),ggzzTemplate.GetBinLowEdge(1),ggzzTemplate.GetBinLowEdge(ggzzTemplate.GetNbinsX()+1)

    zjetsTempFile = ROOT.TFile(template_zjetsName,"READ")
    if decimal[obsName] and doubleDiff: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high))
    elif decimal[obsName]: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    elif not doubleDiff: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high)))
    elif doubleDiff: zjetsTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    # if obsName == 'rapidity4l': zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    # elif doubleDiff: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    print 'zjets bins',zjetsTemplate.GetNbinsX(),zjetsTemplate.GetBinLowEdge(1),zjetsTemplate.GetBinLowEdge(zjetsTemplate.GetNbinsX()+1)

    os.chdir('../../../datacard/datacard_'+year)

    binscale = 3
    qqzzTemplateNew = ROOT.TH1F("qqzzTemplateNew","qqzzTemplateNew",binscale*qqzzTemplate.GetNbinsX(),qqzzTemplate.GetBinLowEdge(1),qqzzTemplate.GetBinLowEdge(qqzzTemplate.GetNbinsX()+1))
    for i in range(1,qqzzTemplate.GetNbinsX()+1):
        for j in range(binscale):
            qqzzTemplateNew.SetBinContent((i-1)*binscale+j+1,qqzzTemplate.GetBinContent(i)/binscale)
    ggzzTemplateNew = ROOT.TH1F("ggzzTemplateNew","ggzzTemplateNew",binscale*ggzzTemplate.GetNbinsX(),ggzzTemplate.GetBinLowEdge(1),ggzzTemplate.GetBinLowEdge(ggzzTemplate.GetNbinsX()+1))
    for i in range(1,ggzzTemplate.GetNbinsX()+1):
        for j in range(binscale):
            ggzzTemplateNew.SetBinContent((i-1)*binscale+j+1,ggzzTemplate.GetBinContent(i)/binscale)
    zjetsTemplateNew = ROOT.TH1F("zjetsTemplateNew","zjetsTemplateNew",binscale*zjetsTemplate.GetNbinsX(),zjetsTemplate.GetBinLowEdge(1),zjetsTemplate.GetBinLowEdge(zjetsTemplate.GetNbinsX()+1))
    for i in range(1,zjetsTemplate.GetNbinsX()+1):
        for j in range(binscale):
            zjetsTemplateNew.SetBinContent((i-1)*binscale+j+1,zjetsTemplate.GetBinContent(i)/binscale)

    qqzzTemplateName = "qqzz_"+channel+recobin+year
    ggzzTemplateName = "ggzz_"+channel+recobin+year
    zjetsTemplateName = "zjets_"+channel+recobin+year

    qqzzTempDataHist = ROOT.RooDataHist(qqzzTemplateName,qqzzTemplateName,ROOT.RooArgList(m),qqzzTemplateNew)
    ggzzTempDataHist = ROOT.RooDataHist(ggzzTemplateName,ggzzTemplateName,ROOT.RooArgList(m),ggzzTemplateNew)
    zjetsTempDataHist = ROOT.RooDataHist(zjetsTemplateName,zjetsTemplateName,ROOT.RooArgList(m),zjetsTemplateNew)

    qqzzTemplatePdf = ROOT.RooHistPdf("qqzz","qqzz",ROOT.RooArgSet(m),qqzzTempDataHist)
    ggzzTemplatePdf = ROOT.RooHistPdf("ggzz","ggzz",ROOT.RooArgSet(m),ggzzTempDataHist)
    zjetsTemplatePdf = ROOT.RooHistPdf("zjets","zjets",ROOT.RooArgSet(m),zjetsTempDataHist)

    # bkg fractions in reco bin; implemented in terms of fractions

    #if( not (obsName=='nJets' or ("jet" in obsName) ) or (not doJES)) :
    qqzz_norm = ROOT.RooFormulaVar("bkg_qqzz_norm", "@0", ROOT.RooArgList(frac_qqzz_var) )
    ggzz_norm = ROOT.RooFormulaVar("bkg_ggzz_norm", "@0", ROOT.RooArgList(frac_ggzz_var) )
    zjets_norm = ROOT.RooFormulaVar("bkg_zjets_norm", "@0", ROOT.RooArgList(frac_zjets_var) )
    #else :
    #    qqzz_norm = ROOT.RooFormulaVar("bkg_qqzz_norm", "@0*(1-@1)", ROOT.RooArgList(frac_qqzz_var, JES_qqzz_rfv) )
    #    ggzz_norm = ROOT.RooFormulaVar("bkg_ggzz_norm", "@0*(1-@1)", ROOT.RooArgList(frac_ggzz_var, JES_ggzz_rfv) )
    #    zjets_norm = ROOT.RooFormulaVar("bkg_zjets_norm", "@0*(1-@1)", ROOT.RooArgList(frac_zjets_var, JES_zjets_rfv) )

    # Data
    # if not os.path.isfile('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/reducedTree_'+year+'.root'):
    #     os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS')
    #     os.system('c++ -o  skim_data_tree skim_data_tree.cpp `root-config --cflags --glibs`')
    #     os.system('./skim_data_tree '+year)
    #     os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard_'+year)
    data_obs_file = ROOT.TFile(path['eos_path']+'Data/reducedTree_AllData_'+year+'.root')
    data_obs_tree = data_obs_file.Get('SR')

    print obsName,obsBin_low,obsBin_high
    chan = ROOT.RooRealVar("chan", "chan", 0, 3)
    # if (obsName == "nJets"): obsName = "njets_reco_pt30_eta4p7"
    if (channel=='4mu'):
        if (obsName.startswith("mass4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && chan == 1)")
        # elif (obsName.startswith("rapidity4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>105.0 && CMS_zz4l_mass<160.0 && abs("+obsName_help+")>="+obsBin_low+" && abs("+obsName_help+")<"+obsBin_high+" && chan == 1)")
        elif not doubleDiff and not obsName.startswith("mass4l"): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && chan == 1)")
        elif doubleDiff: data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,observable_2nd,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && "+obsName_2nd_help+">="+str(obsBin_2nd_low)+" && "+obsName_2nd_help+"<"+str(obsBin_2nd_high)+" && chan == 1)")
        print data_obs.numEntries()
    if (channel=='4e'):
        if (obsName.startswith("mass4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && chan == 2)")
        # elif (obsName.startswith("rapidity4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>105.0 && CMS_zz4l_mass<160.0 && abs("+obsName_help+")>="+obsBin_low+" && abs("+obsName_help+")<"+obsBin_high+" && chan == 2)")
        elif not doubleDiff and not obsName.startswith("mass4l"): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && chan == 2)")
        elif doubleDiff: data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,observable_2nd,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && "+obsName_2nd_help+">="+str(obsBin_2nd_low)+" && "+obsName_2nd_help+"<"+str(obsBin_2nd_high)+" && chan == 2)")
        print data_obs.numEntries()
    if (channel=='2e2mu'):
        if (obsName.startswith("mass4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && chan == 3)")
        # elif (obsName.startswith("rapidity4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>105.0 && CMS_zz4l_mass<160.0 && abs("+obsName_help+")>="+obsBin_low+" && abs("+obsName_help+")<"+obsBin_high+" && chan == 3)")
        elif not doubleDiff and not obsName.startswith("mass4l"): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && chan == 3)")
        elif doubleDiff: data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,observable_2nd,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && "+obsName_2nd_help+">="+str(obsBin_2nd_low)+" && "+obsName_2nd_help+"<"+str(obsBin_2nd_high)+" && chan == 3)")
        print data_obs.numEntries()
    data_obs_file.Close()

    # Create workspace
    wout = ROOT.RooWorkspace("w","w")
    if (channel=='2e2mu'):
        if (year == '2016'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_3_2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_3_2016,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2017'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_3_2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_3_2017,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2018'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_3_2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_3_2018,ROOT.RooFit.RecycleConflictNodes())
    if (channel=='4e'):
        if (year == '2016'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_2_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_2_2016,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2017'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_2_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_2_2017,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2018'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_2_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_2_2018,ROOT.RooFit.RecycleConflictNodes())

    if (channel=='4mu'):
        if (year == '2016'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_1_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_1_2016,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2017'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_1_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_1_2017,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2018'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_1_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_1_2018,ROOT.RooFit.RecycleConflictNodes())

    for genbin in range(nBins):
        getattr(wout,'import')(trueH_shape[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        getattr(wout,'import')(trueH_norm[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        if(physicalModel!='v3'): continue
        getattr(wout,'import')(ggH_shape[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        getattr(wout,'import')(xH_shape[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        getattr(wout,'import')(GGH_norm[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        getattr(wout,'import')(XH_norm[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

    getattr(wout,'import')(OutsideAcceptance,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
    getattr(wout,'import')(OutsideAcceptance_norm,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

    getattr(wout,'import')(nonResH,ROOT.RooFit.Silence())
    getattr(wout,'import')(nonResH_norm,ROOT.RooFit.Silence())

    #print "trueH norm: ",n_trueH,"nonResH norm:",n_nonResH
    qqzzTemplatePdf.SetName("bkg_qqzz")
    qqzzTemplatePdf.Print("v")
    getattr(wout,'import')(qqzzTemplatePdf,ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
    getattr(wout,'import')(qqzz_norm,ROOT.RooFit.Silence())

    ggzzTemplatePdf.SetName("bkg_ggzz")
    ggzzTemplatePdf.Print("v")
    getattr(wout,'import')(ggzzTemplatePdf,ROOT.RooFit.RecycleConflictNodes())
    getattr(wout,'import')(ggzz_norm,ROOT.RooFit.Silence())

    zjetsTemplatePdf.SetName("bkg_zjets")
    zjetsTemplatePdf.Print("v")
    getattr(wout,'import')(zjetsTemplatePdf, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
    getattr(wout,'import')(zjets_norm,ROOT.RooFit.Silence())

    ## data
    # getattr(wout,'import')(data_obs.reduce(ROOT.RooArgSet(m)),ROOT.RooFit.Silence()) #AT https://root-forum.cern.ch/t/rooworkspace-import-roofit-silence-does-not-work-when-importing-datasets/32591
    getattr(wout,'import')(data_obs.reduce(ROOT.RooArgSet(m)))

    if (addfakeH):
        fout = ROOT.TFile('hzz4l_'+channel+'S_13TeV_xs_'+modelName+'_'+_obsName[obsName]+'_'+physicalModel+'.Databin'+str(obsBin)+'.root','RECREATE')
    else:
        fout = ROOT.TFile('hzz4l_'+channel+'S_13TeV_xs_'+modelName+'_'+_obsName[obsName]+'_'+physicalModel+'.Databin'+str(obsBin)+'.NoFakeH.root','RECREATE')

    print "write ws to fout"
    fout.WriteTObject(wout)
    fout.Close()
    os.chdir('../../fit')
    return data_obs.numEntries()
