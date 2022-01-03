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

sys.path.append('../inputs/')
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

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


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
        fname = eos_path + 'MC_samples/%i_MELA' %year
        # if year == 2016:
        #     fname += '_CorrectBTag'
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'
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
        fname = eos_path + 'MC_samples/%i_MELA' %year
        # if year == 2016:
        #     fname += '_CorrectBTag'
        if (year == 2018) & (bkg == 'ZZTo4lext'):
            bkg += '1'
        fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("Counters")
        gen_bkg[bkg] = hCounters.GetBinContent(40)

    return gen_bkg

# Jets variables
def add_njets(pt,eta):
    n = 0
    for i in range(len(pt)):
        if pt[i]>30 and abs(eta[i])<2.5: n=n+1
    return n
def add_leadjet(pt,eta):
    _pTj1 = 0.0
    if len(pt) == 0:
        return _pTj1
    else:
        for i in range(len(pt)):
            if (pt[i]>30 and abs(eta[i])<2.5 and pt[i] > _pTj1): _pTj1 = pt[i]
        return _pTj1

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
        b_bkg = ['ZZMass', 'ZZPt', 'Z1Mass', 'Z2Mass', 'Z1Flav', 'Z2Flav', 'ZZEta', 'LepPt',
                 'overallEventWeight', 'L1prefiringWeight', 'JetPt', 'JetEta',
                 'costhetastar', 'helcosthetaZ1','helcosthetaZ2','helphi','phistarZ1',
                 'pTHj', 'TCjmax', 'TBjmax', 'mjj', 'pTj1', 'pTj2', 'mHj', 'mHjj', 'pTHjj',
                 'njets_pt30_eta4p7',
                 'Dcp', 'D0m', 'D0hp', 'Dint', 'DL1', 'DL1int', 'DL1Zg', 'DL1Zgint']
        if (bkg == 'ZZTo4lext') | (bkg == 'ZZTo4lext1'):
            b_bkg.append('KFactor_EW_qqZZ'); b_bkg.append('KFactor_QCD_qqZZ_M')
        else:
            b_bkg.append('KFactor_QCD_ggZZ_Nominal')
        gen = gen_bkg[bkg]
        xsec = xsec_bkg[bkg]
        df = d_bkg[bkg].pandas.df(b_bkg, flatten = False)
        df['FinState'] = [add_fin_state(i, j) for i,j in zip(df.Z1Flav, df.Z2Flav)]
        df['njets_pt30_eta2p5'] = [add_njets(i,j) for i,j in zip(df['JetPt'],df['JetEta'])]
        #df['pTj1'] = [add_leadjet(i,j) for i,j in zip(df['JetPt'],df['JetEta'])]
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
    fnameFR = eos_path_FR + 'FRfiles/newData_FakeRates_SS_%i.root' %year
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
    keyZX = 'CRZLL'
    data = eos_path + 'Data/reducedTree_AllData_'+str(year)+'.root'
    ttreeZX = uproot.open(data)[keyZX]
    dfZX = ttreeZX.pandas.df(branches_ZX, flatten = False)
    dfZX = dfZX[dfZX.Z2Flav > 0] #Keep just same-sign events
    dfZX = findFSZX(dfZX)
    dfZX['yield_SR'] = ZXYield(dfZX, year)
    return dfZX

# ------------------------------- FUNCTIONS FOR TEMPLATES ----------------------------------------------------
def smoothAndNormaliseTemplate(h1d, norm):
    #smooth
    h1d.Smooth()#10000)
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
        checkDir(str(year)+"/"+var_string)
        fractionBkg = {}
        nBins = len(obs_bins)
        if not doubleDiff: nBins = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)
        # qqzz and ggzz
        for bkg in ['qqzz', 'ggzz']:
            for f in ['2e2mu', '4e', '4mu']:
                #df = df_irr[year][bkg][(df_irr[year][bkg].FinState == f) & (df_irr[year][bkg].Z2Mass < 60)  & (df_irr[year][bkg].ZZMass >= 105) & (df_irr[year][bkg].ZZMass <= 160)].copy()
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

                    df = df_irr[year][bkg][sel].copy()
                    len_bin = df['weight'].sum() # Number of bkg events in bin i
                    if(len_tot <= 0):
                        fractionBkg[bkg+'_'+f+'_'+var_string+'_recobin'+str(i)] = 0.0
                    else:
                        fractionBkg[bkg+'_'+f+'_'+var_string+'_recobin'+str(i)] = float(len_bin/len_tot)
                    # ------
                    sel = sel_bin_low & sel_bin_high & sel_fstate
                    if doubleDiff: sel &= sel_bin_2nd_low & sel_bin_2nd_high
                    df = df_irr[year][bkg][sel].copy()
                    mass4l = df['ZZMass'].to_numpy()
                    mass4l = np.asarray(mass4l).astype('float')
                    w = df['weight'].to_numpy()
                    w = np.asarray(w).astype('float')
                    # ------

                    if((obs_name == 'rapidity4l') | acFlag):
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                    elif doubleDiff:
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+"_"+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                    else:
                        histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)

                    print (histo.GetName())
                    histo.FillN(len(mass4l), mass4l, w)
                    smoothAndNormaliseTemplate(histo, 1)

                    if ((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | acFlag):
                        outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_"+bkg+"_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+".root", "RECREATE")
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
            #df = df_red[year][(sel_f_state_zx) & (df_red[year].Z2Mass < 60) & (df_red[year].ZZMass >= 105) & (df_red[year].ZZMass <=160)].copy()
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
                if((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | acFlag):
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), "m4l_"+var_string+"_"+str(bin_low)+"_"+str(bin_high), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                elif doubleDiff:
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+"_"+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high))+str(int(bin_low_2nd))+"_"+str(int(bin_high_2nd)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                else:
                    histo = ROOT.TH1D("m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), "m4l_"+var_string+"_"+str(int(bin_low))+"_"+str(int(bin_high)), 20, opt.LOWER_BOUND, opt.UPPER_BOUND)
                histo.FillN(len(mass4l), mass4l, w)
                smoothAndNormaliseTemplate(histo, 1)
                if((obs_name == 'rapidity4l') | ('cos' in obs_name) | ('phi' in obs_name) | acFlag):
                    outFile = ROOT.TFile.Open(str(year)+"/"+var_string+"/XSBackground_ZJetsCR_"+f+"_"+var_string+"_"+str(bin_low)+"_"+str(bin_high)+".root", "RECREATE")
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
bkgs = ['ZZTo4lext', 'ggTo2e2mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701',
        'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701']
eos_path_FR = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/'
eos_path = '/eos/user/a/atarabin/'
key = 'candTree'
# years = [2016, 2017, 2018]

if (opt.YEAR == '2016'): years = [2016]
if (opt.YEAR == '2017'): years = [2017]
if (opt.YEAR == '2018'): years = [2018]
if (opt.YEAR == 'Full'): years = [2016,2017,2018]

if not 'vs' in opt.OBSBINS: #It is not a double-differential analysis
    obs_bins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
    obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
    doubleDiff = False
    print 'It is a single-differential measurement, binning', obs_bins
else: #It is a double-differential analysis
    doubleDiff = True
    # The structure of obs_bins is:
    # index of the dictionary is the number of the bin
    # [obs_bins_low, obs_bins_high, obs_bins_low_2nd, obs_bins_high_2nd]
    # The first two entries are the lower and upper bound of the first variable
    # The second two entries are the lower and upper bound of the second variable
    if opt.OBSBINS.count('vs')==1 and opt.OBSBINS.count('/')>=1: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|'
        obs_bins_tmp = opt.OBSBINS.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|']
        obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
        obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
        obs_bins_tmp = obs_bins_tmp[1].split(' / ') #['|0|10|20|45|90|250|', '|0|10|20|80|250|', '|0|20|90|250|', '|0|25|250|']
        obs_bins_2nd = {}
        for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_2nd {0: ['0', '10', '20', '45', '90', '250'], 1: ['0', '10', '20', '80', '250'], 2: ['0', '20', '90', '250'], 3: ['0', '25', '250']}
            obs_bins_2nd[i] = obs_bins_tmp[i].split('|')[1:len(obs_bins_tmp[i].split('|'))-1]
            obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
        obs_bins = {}
        k = 0 #Bin index
        for i in range(len(obs_bins_1st)-1):
            for j in range(len(obs_bins_2nd[i])-1):
                obs_bins[k] = []
                obs_bins[k].append(obs_bins_1st[i])
                obs_bins[k].append(obs_bins_1st[i+1])
                obs_bins[k].append(obs_bins_2nd[i][j])
                obs_bins[k].append(obs_bins_2nd[i][j+1])
                k +=1
    elif opt.OBSBINS.count('vs')>1 and opt.OBSBINS.count('/')>1: #Situation like this one '|50|80| vs |10|30| / |50|80| vs |30|60| / |80|110| vs |10|25| / |80|110| vs |25|30|'
        obs_bins_tmp = opt.OBSBINS.split(' / ') #['|50|80| vs |10|30|', '|50|80| vs |30|60|', '|80|110| vs |10|25|', '|80|110| vs |25|30|']
        obs_bins_1st={}
        obs_bins_2nd={}
        obs_bins={}
        for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_1st {0: ['50', '80'], 1: ['50', '80'], 2: ['80', '110'], 3: ['80', '110']} and obs_bins_2nd {0: ['10', '30'], 1: ['30', '60'], 2: ['10', '25'], 3: ['25', '30']}
            obs_bins_tmp_bis = obs_bins_tmp[i].split(' vs ')
            obs_bins_1st[i] = obs_bins_tmp_bis[0].split('|')[1:len(obs_bins_tmp_bis[0].split('|'))-1]
            obs_bins_1st[i] = [float(j) for j in obs_bins_1st[i]] #Convert a list of str to a list of float
            obs_bins_2nd[i] = obs_bins_tmp_bis[1].split('|')[1:len(obs_bins_tmp_bis[1].split('|'))-1]
            obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
            obs_bins[i] = []
            obs_bins[i].append(obs_bins_1st[i][0])
            obs_bins[i].append(obs_bins_1st[i][1])
            obs_bins[i].append(obs_bins_2nd[i][0])
            obs_bins[i].append(obs_bins_2nd[i][1])
    elif opt.OBSBINS.count('vs')==1 and opt.OBSBINS.count('/')==0: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250|'
        obs_bins_tmp = opt.OBSBINS.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250|']
        obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
        obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
        obs_bins_2nd = obs_bins_tmp[1].split('|')[1:len(obs_bins_tmp[1].split('|'))-1] #['0', '10', '20', '45', '90', '250']
        obs_bins_2nd = [float(i) for i in obs_bins_2nd] #Convert a list of str to a list of float
        obs_bins = {}
        k = 0 #Bin index
        for i in range(len(obs_bins_1st)-1):
            for j in range(len(obs_bins_2nd)-1):
                obs_bins[k] = []
                obs_bins[k].append(obs_bins_1st[i])
                obs_bins[k].append(obs_bins_1st[i+1])
                obs_bins[k].append(obs_bins_2nd[j])
                obs_bins[k].append(obs_bins_2nd[j+1])
                k +=1
    else:
        print 'Problem in the definition of the binning'
        quit()
    print 'It is a double-differential measurement, binning for the 1st variable', obs_bins_1st, 'and for the 2nd variable', obs_bins_2nd
    print obs_bins

obs_name = opt.OBSNAME
if obs_name == 'D0m' or obs_name == 'D0hp' or obs_name == 'Dcp' or obs_name == 'Dint' or obs_name == 'DL1': acFlag = True
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
d_bkg = {}
for year in years:
    bkg = skim_df(year)
    d_bkg[year] = bkg

# Generate pandas for ZX
branches_ZX = ['ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt', 'Z1Mass', 'Z2Mass', 'ZZPt',
               'ZZEta', 'JetPt', 'JetEta', 'costhetastar', 'helcosthetaZ1','helcosthetaZ2','helphi','phistarZ1',
               'pTHj', 'TCjmax', 'TBjmax', 'mjj', 'pTj1', 'pTj2', 'mHj', 'mHjj', 'pTHjj', 'njets_pt30_eta4p7',
               'Dcp', 'D0m', 'D0hp', 'Dint', 'DL1', 'DL1int', 'DL1Zg', 'DL1Zgint']
dfZX={}
for year in years:
    g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE = openFR(year)
    dfZX[year] = doZX(year)
    dfZX[year]['njets_pt30_eta2p5'] = [add_njets(i,j) for i,j in zip(dfZX[year]['JetPt'],dfZX[year]['JetEta'])]
    #dfZX[year]['pTj1'] = [add_leadjet(i,j) for i,j in zip(dfZX[year]['JetPt'],dfZX[year]['JetEta'])]
    dfZX[year] = add_rapidity(dfZX[year])

    # #-*-*-*-*-*-*-*-*-*-*-*-* Temporary - since we do not still have discriminantsa in data, we perform a random generation -*-*-*-*-*-*-*-*-*-*-*-*
    # dfZX[year]['Dcp'] = np.random.choice(d_bkg[year]['qqzz']['Dcp'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['D0m'] = np.random.choice(d_bkg[year]['qqzz']['D0m'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['D0hp'] = np.random.choice(d_bkg[year]['qqzz']['D0hp'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['Dint'] = np.random.choice(d_bkg[year]['qqzz']['Dint'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['DL1'] = np.random.choice(d_bkg[year]['qqzz']['DL1'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['DL1int'] = np.random.choice(d_bkg[year]['qqzz']['DL1int'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['DL1Zg'] = np.random.choice(d_bkg[year]['qqzz']['DL1Zg'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # dfZX[year]['DL1Zgint'] = np.random.choice(d_bkg[year]['qqzz']['DL1Zgint'].tolist(), len(dfZX[year]['ZZMass']), replace=False)
    # #-*-*-*-*-*-*-*-*-*-*-*-* Temporary - since we do not still have discriminantsa in data, we perform a random generation -*-*-*-*-*-*-*-*-*-*-*-*

    print(year,'done')

yield_bkg = {}
if not doubleDiff:doTemplates(d_bkg, dfZX, obs_bins, obs_reco, obs_name)
else: doTemplates(d_bkg, dfZX, obs_bins, obs_reco, obs_name, obs_reco_2nd)

#Write file with expected background yields
with open('../inputs/inputs_bkgTemplate_'+obs_name+'.py', 'w') as f:
    f.write('observableBins = '+str(obs_bins)+';\n')
    f.write('expected_yield = '+str(yield_bkg))
