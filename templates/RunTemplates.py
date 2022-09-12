import matplotlib
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import pandas as pd
import uproot3 as uproot
from math import sqrt, log
import itertools
import optparse
import math
import ROOT
import json
# from tdrStyle import *
import random #-*-*-*-*-*-*-*-*-*-*-*-* Temporary - since we do not have discriminantsa in data yet, we perform a random generation -*-*-*-*-*-*-*-*-*-*-*-*
from binning import binning
from paths import path

# sys.path.append('../inputs/')
# from observables import observables

print 'Welcome in RunTemplates!'

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='int',default=105,   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='int',default=140,   help='Upper bound for m4l')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
    #     parser.error('Bin boundaries not specified for differential measurement. Exiting...')
    #     sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

def checkDir(folder_path):
    isdir = os.path.isdir(folder_path)
    if not isdir:
        print('Directory {} does not exist. Creating it.' .format(folder_path))
        os.mkdir(folder_path)

# ------------------------------- FUNCTIONS TO GENERATE DATAFRAME FOR ggZZ AND qqZZ ----------------------------------------------------
# Weights for histogram
def weight(df, xsec, gen, lumi, additional = None):
    weight = (lumi * 1000 * xsec * df.overallEventWeight * df.SFcorr * df.L1prefiringWeight) / gen #Common structure
    if additional == 'ggH':
        weight *= df.ggH_NNLOPS_weight
    elif additional == 'qqzz':
        weight *= df.KFactor_EW_qqZZ*df.KFactor_QCD_qqZZ_M
    elif additional == 'ggzz':
        weight *= df.KFactor_QCD_ggZZ_Nominal
    df['weight'] = weight
    return df

# Uproot to generate pandas
def prepareTrees(year):
    d_bkg = {}

    for bkg in bkgs:
        fname = eos_path + 'MC_samples_UL/%s_MELA' %year
        # if year == 2016:
        #     fname += '_CorrectBTag'
        # if (year == 2018) & (bkg == 'ZZTo4lext'):
        #     bkg += '1'
        fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'
        d_bkg[bkg] = uproot.open(fname)[key]

    return d_bkg

# Calculate cross sections
def xsecs(year):
    xsec_bkg = {}
    d_bkg = prepareTrees(year)

    for bkg in bkgs:
        # if (year == 2018) & (bkg == 'ZZTo4lext'):
        #     bkg += '1'
        xsec_bkg[bkg] = d_bkg[bkg].pandas.df('xsec').xsec[0]

    return xsec_bkg

# Get the "number" of MC events to divide the weights
def generators(year):
    gen_bkg = {}
    for bkg in bkgs:
        fname = eos_path + 'MC_samples_UL/%s_MELA' %year
        # if year == 2016:
        #     fname += '_CorrectBTag'
        # if (year == 2018) & (bkg == 'ZZTo4lext'):
        #     bkg += '1'
        fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("Counters")
        gen_bkg[bkg] = hCounters.GetBinContent(40)

    return gen_bkg

# Jets variables
def add_njets(pt,eta):
    n = 0
    for i in range(len(pt)):
        if pt[i]>30 and abs(eta[i])<4.7: n=n+1
    return n
def add_leadjet(pt,eta):
    _pTj1 = 0.0
    if len(pt) == 0:
        return _pTj1
    else:
        for i in range(len(pt)):
            if (pt[i]>30 and abs(eta[i])<4.7 and pt[i] > _pTj1): _pTj1 = pt[i]
        return _pTj1

# Rapidity
def rapidity(p, eta):
    return np.abs(np.log((np.sqrt(125*125 + p*p*np.cosh(eta)*np.cosh(eta))+p*np.sinh(eta))/np.sqrt(125*125+p*p)))
def add_rapidity(df):
    df['ZZy'] = rapidity(df['ZZPt'], df['ZZEta'])
    return df

# Define the final state
def add_fin_state(i, j):
    if abs(i) == 121 and abs(j) == 121:
        fin = '4e'
    elif abs(i) == 169 and abs(j) == 169:
        fin = '4mu'
    elif (abs(i) == 121 and abs(j) == 169) or (abs(i) == 169 and abs(j) == 121):
        fin = '2e2mu'
    else: print('Problem with add_fin_state')
    return fin

# Set up data frames
def dataframes(year):
    if year == '2016pre':
        lumi = 19.52
    elif year == '2016post':
        lumi = 16.81
    elif year == '2017':
        lumi = 41.48
    else:
        lumi = 59.83
    d_df_bkg = {}
    d_bkg = prepareTrees(year)
    gen_bkg = generators(year)
    xsec_bkg = xsecs(year)
    for bkg in bkgs:
        # if (year == 2018) & (bkg == 'ZZTo4lext'):
        #     bkg += '1'
        b_bkg = ['ZZMass', 'ZZPt', 'Z1Mass', 'Z2Mass', 'Z1Flav', 'Z2Flav', 'ZZEta', 'LepPt',
                 'overallEventWeight', 'L1prefiringWeight', 'JetPt', 'JetEta',
                 'costhetastar', 'helcosthetaZ1','helcosthetaZ2','helphi','phistarZ1',
                 'pTHj', 'TCjmax', 'TBjmax', 'mjj', 'pTj1', 'pTj2', 'mHj', 'mHjj', 'pTHjj',
                 'njets_pt30_eta4p7', 'absdetajj', 'SFcorr', 'dphijj',
                 'Dcp', 'D0m', 'D0hp', 'Dint', 'DL1', 'DL1int', 'DL1Zg', 'DL1Zgint']
        if (bkg == 'ZZTo4l'):
            b_bkg.append('KFactor_EW_qqZZ'); b_bkg.append('KFactor_QCD_qqZZ_M')
        else:
            b_bkg.append('KFactor_QCD_ggZZ_Nominal')
        gen = gen_bkg[bkg]
        xsec = xsec_bkg[bkg]
        df = d_bkg[bkg].pandas.df(b_bkg, flatten = False)
        df['FinState'] = [add_fin_state(i, j) for i,j in zip(df.Z1Flav, df.Z2Flav)]
        # df['njets_pt30_eta2p5'] = [add_njets(i,j) for i,j in zip(df['JetPt'],df['JetEta'])]
        #df['pTj1'] = [add_leadjet(i,j) for i,j in zip(df['JetPt'],df['JetEta'])]
        df = add_rapidity(df)
        if (bkg != 'ZZTo4l'):
            d_df_bkg[bkg] = weight(df, xsec, gen, lumi, 'ggzz')
        else:
            d_df_bkg[bkg] = weight(df, xsec, gen, lumi, 'qqzz')
    print('Background df created, %s' %year)
    return d_df_bkg

# Sort production modes in view of the histogram (VBF, ggH, Others, qqZZ, ggZZ)
def skim_df(year):
    d_df_bkg = dataframes(year)
    d_skim_bkg = {}
    frames = []
    for bkg in bkgs:
        # if (year == 2018) & (bkg == 'ZZTo4lext'):
        #     bkg += '1'
        if (bkg == 'ZZTo4l'):
            d_skim_bkg['qqzz'] = d_df_bkg[bkg]
        else:
            frames.append(d_df_bkg[bkg])
    d_skim_bkg['ggzz'] = pd.concat(frames)
    print('%s skimmed df created' %year)
    return d_skim_bkg

# ------------------------------- FUNCTIONS TO GENERATE DATAFRAME FOR ZX ----------------------------------------------------
def FindFinalState(z1_flav, z2_flav):
    if(z1_flav == -121):
        if(z2_flav == +121): return 0 # 4e
        if(z2_flav == +169): return 2 # 2e2mu
    if(z1_flav == -169):
        if(z2_flav == +121): return 3 # 2mu2e
        if(z2_flav == +169): return 1 # 4mu

def GetFakeRate(lep_Pt, lep_eta, lep_ID):
    if(lep_Pt >= 80.):
        my_lep_Pt = 79.
    else:
        my_lep_Pt = lep_Pt
    my_lep_ID = abs(lep_ID)
    if((my_lep_Pt > 5) & (my_lep_Pt <= 7)): bin = 0
    if((my_lep_Pt >  7) & (my_lep_Pt <= 10)): bin = 1
    if((my_lep_Pt > 10) & (my_lep_Pt <= 20)): bin = 2
    if((my_lep_Pt > 20) & (my_lep_Pt <= 30)): bin = 3
    if((my_lep_Pt > 30) & (my_lep_Pt <= 40)): bin = 4
    if((my_lep_Pt > 40) & (my_lep_Pt <= 50)): bin = 5
    if((my_lep_Pt > 50) & (my_lep_Pt <= 80)): bin = 6
    if(abs(my_lep_ID) == 11): bin = bin-1 # There is no [5, 7] bin in the electron fake rate
    if(my_lep_ID == 11):
        if(abs(lep_eta) < 1.479): return g_FR_e_EB.GetY()[bin]
        else: return g_FR_e_EE.GetY()[bin]
    if(my_lep_ID == 13):
        if(abs(lep_eta) < 1.2): return g_FR_mu_EB.GetY()[bin]
        else: return g_FR_mu_EE.GetY()[bin]

# Open Fake Rates files
def openFR(year):
    fnameFR = eos_path_FR + 'FRfiles_UL/FakeRates_SS_%i.root' %year
    file = uproot.open(fnameFR)
    # Retrieve FR from TGraphErrors
    input_file_FR = ROOT.TFile(fnameFR)
    g_FR_mu_EB = input_file_FR.Get("FR_SS_muon_EB")
    g_FR_mu_EE = input_file_FR.Get("FR_SS_muon_EE")
    g_FR_e_EB  = input_file_FR.Get("FR_SS_electron_EB")
    g_FR_e_EE  = input_file_FR.Get("FR_SS_electron_EE")
    return g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE

# Find final state
def findFSZX(df):
    df['FinState'] = [FindFinalState(x,y) for x,y in zip(df['Z1Flav'], df['Z2Flav'])]
    return df

# Define combination coefficients
def comb(year):
    if year == 2016:
        cb_SS = np.array([
            1.175,   # 4e
            0.975,   # 4mu
            1.052,    # 2e2mu
            1.148,    # 2mu2e
        ])
    elif year == 2017:
        cb_SS = np.array([
            1.094,   # 4e
            0.948,   # 4mu
            0.930,    # 2e2mu
            1.139,    # 2mu2e
        ])
    else:
        cb_SS = np.array([
            1.157,   # 4e
            0.974,   # 4mu
            0.930,    # 2e2mu
            1.143,    # 2mu2e
        ])
    return cb_SS

# Define ration OppositeSign/SameSign
def ratio(year):
    if year == 2016:
        fs_ROS_SS = np.array([
            1.0039,   # 4e
            0.999103,  # 4mu
            1.0332,   # 2e2mu
            1.00216,  # 2mu2e
            ])
    elif year == 2017:
        fs_ROS_SS = np.array([
            0.990314,   # 4e
            1.02903,  # 4mu
            1.0262,   # 2e2mu
            1.00154,  # 2mu2e
            ])
    else:
        fs_ROS_SS = np.array([
            1.00322,   # 4e
            1.0187,  # 4mu
            1.04216,   # 2e2mu
            0.996253,  # 2mu2e
            ])
    return fs_ROS_SS

# Calculate yield for Z+X (data in CRZLL control region are scaled in signal region through yields)
def ZXYield(df, year):
    cb_SS = comb(year)
    fs_ROS_SS = ratio(year)
    vec = df.to_numpy()
    Yield = np.zeros(len(vec), float)
    for i in range(len(vec)):
        finSt  = vec[i][len(branches_ZX)] #Final state information is in the last column which is added afterthe last column of branches_ZX
        lepPt  = vec[i][5]
        lepEta = vec[i][4]
        lepID  = vec[i][3]
        Yield[i] = cb_SS[finSt] * fs_ROS_SS[finSt] * GetFakeRate(lepPt[2], lepEta[2], lepID[2]) * GetFakeRate(lepPt[3], lepEta[3], lepID[3])
    return Yield

def doZX(year):
    keyZX = 'CRZLL'
    data = eos_path + 'Data_UL/reducedTree_AllData_'+str(year)+'.root'
    ttreeZX = uproot.open(data)[keyZX]
    dfZX = ttreeZX.pandas.df(branches_ZX, flatten = False)
    dfZX = dfZX[dfZX.Z2Flav > 0] #Keep just same-sign events
    dfZX = findFSZX(dfZX)
    dfZX['yield_SR'] = ZXYield(dfZX, year)
    return dfZX

# ------------------------------- FUNCTIONS FOR TEMPLATES ----------------------------------------------------
def smoothAndNormaliseTemplate(h1d, norm):
    #smooth
    h1d.Smooth(10000)
    #norm + floor + norm
    #normaliseHist(h1d, norm)
    fillEmptyBinsHist(h1d,.01/(h1d.GetNbinsX()))
    normaliseHist(h1d, norm)

def normaliseHist(h1d, norm):
    if (h1d.Integral() > 0): #return -1
	h1d.Scale(norm/h1d.Integral())
    else:
        return -1

def fillEmptyBinsHist(h1d, floor):
    nXbins=h1d.GetNbinsX()
    for i in range(1, nXbins+1): h1d.SetBinContent(i, h1d.GetBinContent(i)+floor)

def doTemplates(df_irr, df_red, binning, var, var_string, var_2nd='None'):
    for year in years:
        checkDir(str(year))
        checkDir(str(year)+"/"+var_string)
        fractionBkg = {}
        nBins = len(obs_bins)
        if not doubleDiff: nBins = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)
        # qqzz and ggzz
        for bkg in ['qqzz', 'ggzz']:
            for f in ['2e2mu', '4e', '4mu']:
                df = df_irr[year][bkg][(df_irr[year][bkg].FinState == f) & (df_irr[year][bkg].ZZMass >= opt.LOWER_BOUND) & (df_irr[year][bkg].ZZMass <= opt.UPPER_BOUND)].copy()
                len_tot = df['weight'].sum() # Total number of bkg b events in final state f
                yield_bkg[year,bkg,f] = len_tot
                print(year, bkg, f, len_tot)
                for i in range(nBins):
                    if not doubleDiff:
                        bin_low = binning[i]
                        bin_high = binning[i+1]
                    else:
                        bin_low = binning[i][0]
                        bin_high = binning[i][1]
                        bin_low_2nd = binning[i][2]
                        bin_high_2nd = binning[i][3]
                    sel_bin_low = df_irr[year][bkg][var] >= bin_low
                    sel_bin_high = df_irr[year][bkg][var] < bin_high
                    if doubleDiff:
                        sel_bin_2nd_low = df_irr[year][bkg][var_2nd] >= bin_low_2nd
                        sel_bin_2nd_high = df_irr[year][bkg][var_2nd] < bin_high_2nd
                    sel_bin_mass_low = df_irr[year][bkg].ZZMass >= opt.LOWER_BOUND
                    sel_bin_mass_high = df_irr[year][bkg].ZZMass <= opt.UPPER_BOUND
                    sel_Z2_mass = df_irr[year][bkg].Z2Mass < 60 ## Uncomment below to cut mZ2 at 60 GeV, hence removing non-reso evts
                    sel_fstate = df_irr[year][bkg]['FinState'] == f

                    sel = sel_bin_low & sel_bin_high & sel_bin_mass_low & sel_bin_mass_high & sel_fstate
                    if doubleDiff: sel &= sel_bin_2nd_low & sel_bin_2nd_high

                    if 'zzfloating' in obs_name:
                        df_2016_qqzz = df_irr[2016]['qqzz'][(df_irr[2016]['qqzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2016]['qqzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2016]['qqzz'][var] >= bin_low) & (df_irr[2016]['qqzz'][var] < bin_high)].copy()
                        df_2017_qqzz = df_irr[2017]['qqzz'][(df_irr[2017]['qqzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2017]['qqzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2017]['qqzz'][var] >= bin_low) & (df_irr[2017]['qqzz'][var] < bin_high)].copy()
                        df_2018_qqzz = df_irr[2018]['qqzz'][(df_irr[2018]['qqzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2018]['qqzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2018]['qqzz'][var] >= bin_low) & (df_irr[2018]['qqzz'][var] < bin_high)].copy()
                        df_2016_ggzz = df_irr[2016]['ggzz'][(df_irr[2016]['ggzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2016]['ggzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2016]['ggzz'][var] >= bin_low) & (df_irr[2016]['ggzz'][var] < bin_high)].copy()
                        df_2017_ggzz = df_irr[2017]['ggzz'][(df_irr[2017]['ggzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2017]['ggzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2017]['ggzz'][var] >= bin_low) & (df_irr[2017]['ggzz'][var] < bin_high)].copy()
                        df_2018_ggzz = df_irr[2018]['ggzz'][(df_irr[2018]['ggzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2018]['ggzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2018]['ggzz'][var] >= bin_low) & (df_irr[2018]['ggzz'][var] < bin_high)].copy()
                        df = pd.concat([df_2016_qqzz.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M']),df_2017_qqzz.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M']),
                                        df_2018_qqzz.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M']),df_2016_ggzz.drop(columns=['KFactor_QCD_ggZZ_Nominal']),
                                        df_2017_ggzz.drop(columns=['KFactor_QCD_ggZZ_Nominal']),df_2018_ggzz.drop(columns=['KFactor_QCD_ggZZ_Nominal'])])
                        # In case of zzfloating len_tot is overwritten (previous definition at the beginning of for loops)
                        len_tot = df['weight'].sum() # Total number of bkg b events in all final states and across years
                        yield_bkg['ZZ_'+str(i)] = len_tot
                        #### 2e2mu ####
                        df_2016_qqzz = df_irr[2016]['qqzz'][(df_irr[2016]['qqzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2016]['qqzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2016]['qqzz'][var] >= bin_low) & (df_irr[2016]['qqzz'][var] < bin_high) & (df_irr[2016]['qqzz']['FinState'] == f)].copy()
                        df_2017_qqzz = df_irr[2017]['qqzz'][(df_irr[2017]['qqzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2017]['qqzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2017]['qqzz'][var] >= bin_low) & (df_irr[2017]['qqzz'][var] < bin_high) & (df_irr[2017]['qqzz']['FinState'] == f)].copy()
                        df_2018_qqzz = df_irr[2018]['qqzz'][(df_irr[2018]['qqzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2018]['qqzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2018]['qqzz'][var] >= bin_low) & (df_irr[2018]['qqzz'][var] < bin_high) & (df_irr[2018]['qqzz']['FinState'] == f)].copy()
                        df_2016_ggzz = df_irr[2016]['ggzz'][(df_irr[2016]['ggzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2016]['ggzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2016]['ggzz'][var] >= bin_low) & (df_irr[2016]['ggzz'][var] < bin_high) & (df_irr[2016]['ggzz']['FinState'] == f)].copy()
                        df_2017_ggzz = df_irr[2017]['ggzz'][(df_irr[2017]['ggzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2017]['ggzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2017]['ggzz'][var] >= bin_low) & (df_irr[2017]['ggzz'][var] < bin_high) & (df_irr[2017]['ggzz']['FinState'] == f)].copy()
                        df_2018_ggzz = df_irr[2018]['ggzz'][(df_irr[2018]['ggzz'].ZZMass >= opt.LOWER_BOUND) & (df_irr[2018]['ggzz'].ZZMass <= opt.UPPER_BOUND) & (df_irr[2018]['ggzz'][var] >= bin_low) & (df_irr[2018]['ggzz'][var] < bin_high) & (df_irr[2018]['ggzz']['FinState'] == f)].copy()
                        df = pd.concat([df_2016_qqzz.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M']),df_2017_qqzz.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M']),
                                        df_2018_qqzz.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M']),df_2016_ggzz.drop(columns=['KFactor_QCD_ggZZ_Nominal']),
                                        df_2017_ggzz.drop(columns=['KFactor_QCD_ggZZ_Nominal']),df_2018_ggzz.drop(columns=['KFactor_QCD_ggZZ_Nominal'])])
                        # In case of zzfloating len_tot is overwritten (previous definition at the beginning of for loops)
                        # len_tot = df['weight'].sum() # Total number of bkg b events in all final states and across years
                        yield_bkg['ZZ_'+f] = df['weight'].sum()

                    df = df_irr[year][bkg][sel].copy()
                    len_bin = df['weight'].sum() # Number of bkg events in bin i
                    if(len_tot <= 0):
                        fractionBkg[bkg+'_'+f+'_'+var_string+'_recobin'+str(i)] = 0.0
                    else:
                        fractionBkg[bkg+'_'+f+'_'+var_string+'_recobin'+str(i)] = float(len_bin/len_tot)
                    if 'zzfloating' in obs_name:
                        fractionBkg[bkg+'_'+f+'_'+var_string+'_recobin'+str(i)+'_v2'] = float(len_bin/yield_bkg['ZZ_'+f])
                    # ------
                    sel = sel_bin_low & sel_bin_high & sel_fstate
                    if doubleDiff: sel &= sel_bin_2nd_low & sel_bin_2nd_high
                    df = df_irr[year][bkg][sel].copy()
                    mass4l = df['ZZMass'].to_numpy()
                    mass4l = np.asarray(mass4l).astype('float')
                    w = df['weight'].to_numpy()
                    w = np.asarray(w).astype('float')
                    # ------

                    if((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | ('deta' in obs_name) | acFlag):
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                    elif doubleDiff and 'rapidity' in var_string:
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+"_"+str(bin_low_2nd)+"_"+str(bin_high_2nd), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+'_'+str(bin_low_2nd)+"_"+str(bin_high_2nd), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                    elif doubleDiff:
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+"_"+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+'_'+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                    else:
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)

                    print (histo.GetName())
                    histo.FillN(len(mass4l), mass4l, w)
                    smoothAndNormaliseTemplate(histo, 1)

                    if ((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | ('deta' in obs_name) | acFlag):
                        outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_"+bkg+"_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+".root", "RECREATE")
                    elif doubleDiff and 'rapidity' in var_string:
                        outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_"+bkg+"_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+"_"+str(bin_low_2nd)+"_"+str(bin_high_2nd)+".root", "RECREATE")
                    elif doubleDiff:
                        outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_"+bkg+"_"+f+"_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+"_"+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd))+".root", "RECREATE")
                    else:
                        outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_"+bkg+"_"+f+"_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+".root", "RECREATE")
                    outFile.cd()

                    histo.Write()
                    outFile.Close()
                    histo.Delete()
        # ZX for different final states
        for f in ['2e2mu', '4e', '4mu']:
            if(f == '4e'):
                sel_f_state_zx = df_red[year]['FinState'] == 0
            elif(f == '4mu'):
                sel_f_state_zx = df_red[year]['FinState'] == 1
            elif(f == '2e2mu'):
                sel_f_state_zx = (df_red[year]['FinState'] == 2) | (df_red[year]['FinState'] == 3)
            df = df_red[year][(sel_f_state_zx) & (df_red[year].ZZMass >= opt.LOWER_BOUND) & (df_red[year].ZZMass <=opt.UPPER_BOUND)].copy()
            df_inclusive = df.copy()
            len_tot = df['yield_SR'].sum() # Total number of bkg events in final state f
            yield_bkg[year,'ZX',f] = len_tot
            print(year, f, len_tot)
            for i in range(nBins):
                if not doubleDiff:
                    bin_low = binning[i]
                    bin_high = binning[i+1]
                else:
                    bin_low = binning[i][0]
                    bin_high = binning[i][1]
                    bin_low_2nd = binning[i][2]
                    bin_high_2nd = binning[i][3]
                sel_bin_low = df_red[year][var] >= bin_low
                sel_bin_high = df_red[year][var] < bin_high
                if doubleDiff:
                    sel_bin_2nd_low = df_red[year][var_2nd] >= bin_low_2nd
                    sel_bin_2nd_high = df_red[year][var_2nd] < bin_high_2nd
                sel_bin_mass_low = df_red[year]['ZZMass'] >= opt.LOWER_BOUND
                sel_bin_mass_high = df_red[year]['ZZMass'] <= opt.UPPER_BOUND

                sel_Z2_mass = df_red[year]['Z2Mass'] < 60 ## Uncomment below to cut mZ2 at 60 GeV, hence removing non-reso evts
                sel = sel_bin_low & sel_bin_high & sel_f_state_zx & sel_bin_mass_low & sel_bin_mass_high #& sel_Z2_mass
                if doubleDiff: sel &= sel_bin_2nd_low & sel_bin_2nd_high

                df = df_red[year][sel].copy()
                len_bin = df['yield_SR'].sum() # Number of bkg events in bin i
                fractionBkg['ZJetsCR_'+f+'_'+var_string+'_recobin'+str(i)] = float(len_bin/len_tot)
                # ------
                if(len_bin <= 0): df = df_inclusive
                mass4l = df['ZZMass'].to_numpy()
                mass4l = np.asarray(mass4l).astype('float')
                w = df['yield_SR'].to_numpy()
                w = np.asarray(w).astype('float')
                # ------
                if((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | ('deta' in obs_name) | acFlag):
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                elif doubleDiff and 'rapidity' in var_string:
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+"_"+str(bin_low_2nd)+"_"+str(bin_high_2nd), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+str(bin_low_2nd)+"_"+str(bin_high_2nd), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                elif doubleDiff:
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+"_"+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                else:
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                histo.FillN(len(mass4l), mass4l, w)
                smoothAndNormaliseTemplate(histo, 1)
                if((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | ('deta' in obs_name) | acFlag):
                    outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_ZJetsCR_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+".root", "RECREATE")
                elif doubleDiff and 'rapidity' in var_string:
                    outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_ZJetsCR_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+"_"+str(bin_low_2nd)+"_"+str(bin_high_2nd)+".root", "RECREATE")
                elif doubleDiff:
                    outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_ZJetsCR_"+f+"_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+"_"+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd))+".root", "RECREATE")
                else:
                    outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_ZJetsCR_"+f+"_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+".root", "RECREATE")
                outFile.cd()
                histo.Write()
                outFile.Close()
                histo.Delete()
        with open('../inputs/inputs_bkg_'+var_string+'_'+str(year)+'.py', 'w') as f:
            f.write('observableBins = '+json.dumps(binning)+';\n')
            f.write('fractionsBackground = '+json.dumps(fractionBkg))

# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------

# General settings
bkgs = ['ZZTo4l', 'ggTo2e2mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701',
        'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701']
eos_path_FR = path['eos_path_FR']
eos_path = path['eos_path']
key = 'candTree'

if (opt.YEAR == '2016'):
    years_MC = ['2016pre', '2016post']
    years = [2016]
if (opt.YEAR == '2017'):
    years_MC = ['2017']
    years = [2017]
if (opt.YEAR == '2018'):
    years_MC = ['2018']
    years = [2018]
if (opt.YEAR == 'Full'):
    years_MC = ['2016pre', '2016post', '2017', '2018']
    years = [2016,2017,2018]

obs_bins, doubleDiff = binning(opt.OBSNAME)

obs_name = opt.OBSNAME
if obs_name == 'D0m' or obs_name == 'D0hp' or obs_name == 'Dcp' or obs_name == 'Dint' or obs_name == 'DL1' or obs_name == 'DL1Zg': acFlag = True
else: acFlag = False
print acFlag

_temp = __import__('observables', globals(), locals(), ['observables'], -1)
observables = _temp.observables

if doubleDiff:
    obs_reco_2nd = observables[obs_name]['obs_reco_2nd']
obs_reco = observables[obs_name]['obs_reco']

if doubleDiff:
    obs_name = opt.OBSNAME.split(' vs ')[0]
    obs_name_2nd = opt.OBSNAME.split(' vs ')[1]
    obs_name = obs_name + '_' + obs_name_2nd

print 'Following observables extracted from dictionary: RECO = ',obs_reco
if doubleDiff:
    print 'It is a double-differential measurement: RECO_2nd = ',obs_reco_2nd

# Generate pandas for ggZZ and qqZZ
d_bkg_tmp = {}
for year in years_MC:
    bkg = skim_df(year)
    d_bkg_tmp[year] = bkg

# # Create pandas with int as indeces and 2016post+2016pre
d_bkg = {}
if (opt.YEAR == '2016' or opt.YEAR == 'Full'):
    d_bkg_2016 = {}
    d_bkg_2016['qqzz'] = pd.concat([d_bkg_tmp['2016post']['qqzz'], d_bkg_tmp['2016pre']['qqzz']])
    d_bkg_2016['ggzz'] = pd.concat([d_bkg_tmp['2016post']['ggzz'], d_bkg_tmp['2016pre']['ggzz']])
    d_bkg[2016] = d_bkg_2016
if (opt.YEAR == '2017' or opt.YEAR == 'Full'):
    d_bkg[2017] = d_bkg_tmp['2017']
if (opt.YEAR == '2018' or opt.YEAR == 'Full'):
    d_bkg[2018] = d_bkg_tmp['2018']


# Generate pandas for ZX
branches_ZX = ['ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt', 'Z1Mass', 'Z2Mass', 'ZZPt',
               'ZZEta', 'JetPt', 'JetEta', 'costhetastar', 'helcosthetaZ1','helcosthetaZ2','helphi','phistarZ1',
               'pTHj', 'TCjmax', 'TBjmax', 'mjj', 'pTj1', 'pTj2', 'mHj', 'mHjj', 'pTHjj', 'njets_pt30_eta4p7',
               'Dcp', 'D0m', 'D0hp', 'Dint', 'DL1', 'DL1int', 'DL1Zg', 'DL1Zgint','absdetajj', 'dphijj']
dfZX={}
for year in years:
    g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE = openFR(year)
    dfZX[year] = doZX(year)
    # dfZX[year]['njets_pt30_eta2p5'] = [add_njets(i,j) for i,j in zip(dfZX[year]['JetPt'],dfZX[year]['JetEta'])]
    #dfZX[year]['pTj1'] = [add_leadjet(i,j) for i,j in zip(dfZX[year]['JetPt'],dfZX[year]['JetEta'])]
    dfZX[year] = add_rapidity(dfZX[year])

    print(year,'done')

yield_bkg = {}
if not doubleDiff:doTemplates(d_bkg, dfZX, obs_bins, obs_reco, obs_name)
else: doTemplates(d_bkg, dfZX, obs_bins, obs_reco, obs_name, obs_reco_2nd)

#Write file with expected background yields
with open('../inputs/inputs_bkgTemplate_'+obs_name+'.py', 'w') as f:
    f.write('observableBins = '+str(obs_bins)+';\n')
    f.write('expected_yield = '+str(yield_bkg))
