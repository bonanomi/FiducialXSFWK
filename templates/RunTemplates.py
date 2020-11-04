import matplotlib
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import uproot
from math import sqrt, log
import itertools
import math
import ROOT
import json
from tdrStyle import *

# ------------------------------- FUNCTIONS TO GENERATE DATAFRAME FOR ggZZ AND qqZZ ----------------------------------------------------
# Weights for histogram
def weight(df, xsec, gen, lumi, additional = None):
    weight = (lumi * 1000 * xsec * df.overallEventWeight * df.L1prefiringWeight) / gen #Common structure
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
        fname = eos_path + 'MC_%i' %year
        if year == 2016:
            fname += '_CorrectBTag'
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        fname += '/%s/ZZ4lAnalysis.root' %bkg
        d_bkg[bkg] = uproot.open(fname)[key]

    return d_bkg

# Calculate cross sections
def xsecs(year):
    xsec_bkg = {}
    d_bkg = prepareTrees(year)

    for bkg in bkgs:
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        xsec_bkg[bkg] = d_bkg[bkg].pandas.df('xsec').xsec[0]

    return xsec_bkg

# Get the "number" of MC events to divide the weights
def generators(year):
    gen_bkg = {}
    for bkg in bkgs:
        fname = eos_path + 'MC_%i' %year
        if year == 2016:
            fname += '_CorrectBTag'
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        fname += '/%s/ZZ4lAnalysis.root' %bkg
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("ZZTree/Counters")
        gen_bkg[bkg] = hCounters.GetBinContent(40)

    return gen_bkg

# Finding the leading lepton
def add_lead_lep(df):
    df['Pt_l1'] = [max(i) for i in df.LepPt]
    return df

# Rapidity
def rapidity(p, eta):
    return np.log((np.sqrt(125*125 + p*p*np.cosh(eta)*np.cosh(eta))+p*np.sinh(eta))/np.sqrt(125*125+p*p))
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
    if year == 2016:
        lumi = 35.9
    elif year == 2017:
        lumi = 41.5
    else:
        lumi = 59.7
    d_df_bkg = {}
    d_bkg = prepareTrees(year)
    gen_bkg = generators(year)
    xsec_bkg = xsecs(year)
    for bkg in bkgs:
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        b_bkg = ['ZZMass', 'ZZPt', 'Z1Mass', 'Z2Mass', 'Z1Flav', 'Z2Flav', 'ZZEta', 'LepPt', 'overallEventWeight', 'L1prefiringWeight']
        if (bkg == 'ZZTo4lext') | (bkg == 'ZZTo4lext1'):
            b_bkg.append('KFactor_EW_qqZZ'); b_bkg.append('KFactor_QCD_qqZZ_M')
        else:
            b_bkg.append('KFactor_QCD_ggZZ_Nominal')
        gen = gen_bkg[bkg]
        xsec = xsec_bkg[bkg]
        df = d_bkg[bkg].pandas.df(b_bkg, flatten = False)
#         df = add_lead_lep(df)
        df['FinState'] = [add_fin_state(i, j) for i,j in zip(df.Z1Flav, df.Z2Flav)]
        df = add_rapidity(df)
        if (bkg != 'ZZTo4lext') & (bkg != 'ZZTo4lext1'):
            d_df_bkg[bkg] = weight(df, xsec, gen, lumi, 'ggzz')
        else:
            d_df_bkg[bkg] = weight(df, xsec, gen, lumi, 'qqzz')
    print('Background df created, %i' %year)
    return d_df_bkg

# Sort production modes in view of the histogram (VBF, ggH, Others, qqZZ, ggZZ)
def skim_df(year):
    d_df_bkg = dataframes(year)
    d_skim_bkg = {}
    frames = []
    for bkg in bkgs:
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        if (bkg == 'ZZTo4lext') | (bkg == 'ZZTo4lext1'):
            d_skim_bkg['qqzz'] = d_df_bkg[bkg]
        else:
            frames.append(d_df_bkg[bkg])
    d_skim_bkg['ggzz'] = pd.concat(frames)
    print('%i skimmed df created' %year)
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
    fnameFR = eos_path + 'FRfiles/FakeRates_SS_%i.root' %year
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
            1.23628,   # 4e
            0.95433,   # 4mu
            1.0726,    # 2e2mu
            1.0726,    # 2mu2e
        ])
    elif year == 2017:
        cb_SS = np.array([
            1.1934,   # 4e
            0.99669,   # 4mu
            1.0569,    # 2e2mu
            1.0569,    # 2mu2e
        ])
    else:
        cb_SS = np.array([
            1.2087,   # 4e
            0.9878,   # 4mu
            1.0552,    # 2e2mu
            1.0552,    # 2mu2e
        ])
    return cb_SS

# Define ration OppositeSign/SameSign
def ratio(year):
    if year == 2016:
        fs_ROS_SS = np.array([
            1.00245,   # 4e
            0.998863,  # 4mu
            1.03338,   # 2e2mu
            0.998852,  # 2mu2e
            ])
    elif year == 2017:
        fs_ROS_SS = np.array([
            1.01198,   # 4e
            1.03949,  # 4mu
            1.013128,   # 2e2mu
            1.00257,  # 2mu2e
            ])
    else:
        fs_ROS_SS = np.array([
            1.00568,   # 4e
            1.02926,  # 4mu
            1.03226,   # 2e2mu
            1.00432,  # 2mu2e
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
    keyZX = 'CRZLLTree/candTree'
    data = eos_path + 'Data_%i/AllData/ZZ4lAnalysis.root' %year
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
    normaliseHist(h1d, norm)
    fillEmptyBinsHist(h1d,.001/(h1d.GetNbinsX()))
    normaliseHist(h1d, norm)

def normaliseHist(h1d, norm):
    if (h1d.Integral() == 0): return -1
    h1d.Scale(norm/h1d.Integral())

def fillEmptyBinsHist(h1d, floor):
    nXbins=h1d.GetNbinsX()
    for i in range(1, nXbins+1): h1d.SetBinContent(i, h1d.GetBinContent(i)+floor)

def doTemplates(df_irr, df_red, binning, var, var_string):
    for year in years:
        fractionBkg = {}
        # qqzz and ggzz
        for bkg in ['qqzz', 'ggzz']:
            for f in ['2e2mu', '4e', '4mu']:
                df = df_irr[year][bkg][(df_irr[year][bkg].FinState == f) & (df_irr[year][bkg].ZZMass >= 105) & (df_irr[year][bkg].ZZMass <= 140)].copy()
                len_tot = df['weight'].sum() # Total number of bkg b events in final state f
                for i in range(len(binning)-1):
                    bin_low = binning[i]
                    bin_high = binning[i+1]
                    sel_bin_low = df_irr[year][bkg][var] >= bin_low
                    sel_bin_high = df_irr[year][bkg][var] < bin_high
                    sel_bin_mass_low = df_irr[year][bkg].ZZMass >= 105
                    sel_bin_mass_high = df_irr[year][bkg].ZZMass <= 140
                    sel_fstate = df_irr[year][bkg]['FinState'] == f
                    sel = sel_bin_low & sel_bin_high & sel_bin_mass_low & sel_bin_mass_high & sel_fstate
                    df = df_irr[year][bkg][sel].copy()
                    len_bin = df['weight'].sum() # Number of bkg events in bin i
                    fractionBkg[bkg+'_'+f+'_'+var_string+'_recobin'+str(i)] = float(len_bin/len_tot)
                    # ------
                    sel = sel_bin_low & sel_bin_high & sel_fstate
                    df = df_irr[year][bkg][sel].copy()
                    mass4l = df['ZZMass'].to_numpy()
                    mass4l = np.asarray(mass4l).astype('float')
                    w = df['weight'].to_numpy()
                    w = np.asarray(w).astype('float')
                    # ------
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), 20, 105, 140)
                    histo.FillN(len(mass4l), mass4l, w)
                    smoothAndNormaliseTemplate(histo, 1)
                    outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_"+bkg+"_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+".root", "RECREATE")
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
            df = df_red[year][(sel_f_state_zx) & (df_red[year].ZZMass >= 105) & (df_red[year].ZZMass <=140)].copy()
            len_tot = df['yield_SR'].sum() # Total number of bkg events in final state f
            for i in range(len(binning)-1):
                bin_low = binning[i]
                bin_high = binning[i+1]
                sel_bin_low = df_red[year][var] >= bin_low
                sel_bin_high = df_red[year][var] < bin_high
                sel_bin_mass_low = df_red[year]['ZZMass'] >= 105
                sel_bin_mass_high = df_red[year]['ZZMass'] <= 140
                sel = sel_bin_low & sel_bin_high & sel_f_state_zx & sel_bin_mass_low & sel_bin_mass_high
                df = df_red[year][sel].copy()
                len_bin = df['yield_SR'].sum() # Number of bkg events in bin i
                fractionBkg['ZJetsCR_'+f+'_'+var_string+'_recobin'+str(i)] = float(len_bin/len_tot)
                # ------
                mass4l = df['ZZMass'].to_numpy()
                mass4l = np.asarray(mass4l).astype('float')
                w = df['yield_SR'].to_numpy()
                w = np.asarray(w).astype('float')
                # ------
                histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), 20, 105, 140)
                histo.FillN(len(mass4l), mass4l, w)
                smoothAndNormaliseTemplate(histo, 1)
                outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_ZJetsCR_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+".root", "RECREATE")
                outFile.cd()
                histo.Write()
                outFile.Close()
                histo.Delete()
        # with open('../inputs/inputs_bkg_'+var_string+'_'+str(year)+'.py', 'w') as f:
        #     f.write('observableBins = '+json.dumps(binning)+';\n')
        #     f.write('fractionsBackground = '+json.dumps(fractionBkg))

# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------

# General settings
bkgs = ['ZZTo4lext', 'ggTo2e2mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701',
        'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701']
eos_path = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/'
key = 'ZZTree/candTree'
years = [2016, 2017, 2018]

# Generate pandas for ggZZ and qqZZ
d_bkg = {}
for year in years:
    bkg = skim_df(year)
    d_bkg[year] = bkg

# Generate pandas for ZX
branches_ZX = ['ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt', 'Z1Mass', 'Z2Mass', 'ZZPt', 'ZZEta']
dfZX={}
for year in years:
    g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE = openFR(year)
    dfZX[year] = doZX(year)
    dfZX[year] = add_lead_lep(dfZX[year])
    dfZX[year] = add_rapidity(dfZX[year])
    print(year,'done')

doTemplates(d_bkg, dfZX, [0,15,30,45,80,120,200,1300], 'ZZPt', 'pT4l')
doTemplates(d_bkg, dfZX, [0.0,0.15,0.3,0.6,0.9,1.2,2.5], 'ZZy', 'rapidity4l')
