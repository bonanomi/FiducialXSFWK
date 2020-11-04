import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import uproot
from math import sqrt, log
import sys
import itertools
import math
import ROOT
import json


# ------------------------------- FUNCTIONS TO GENERATE DATAFRAME FOR ggZZ AND qqZZ ----------------------------------------------------
# Gen-Reco lepton matching
def GenMatching(lep_pt, lep_eta, lep_phi, lep_id, genlep_pt, genlep_eta, genlep_phi, genlep_id):

    lep_genindex = [-1, -1, -1, -1]
    index = [0,1,2,3]
    for i in range(4): #Reco
        minDr = 9999.0

        reco = ROOT.TLorentzVector()
        gen = ROOT.TLorentzVector()

        if abs(lep_id[i]) == 11: mass = 0.0005109989461
        if abs(lep_id[i]) == 13: mass = 0.1056583745
        if abs(lep_id[i]) == 15: mass = 1.77686
        reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],mass)

        for j in index: #Gen

            if (genlep_id[j] != lep_id[i]): continue

            if abs(genlep_id[j]) == 11: mass = 0.0005109989461
            if abs(genlep_id[j]) == 13: mass = 0.1056583745
            if abs(genlep_id[j]) == 15: mass = 1.77686
            gen.SetPtEtaPhiM(genlep_pt[j],genlep_eta[j],genlep_phi[j],mass)
            thisDr = ROOT.TLorentzVector.DeltaR(reco, gen)

            if thisDr<minDr and thisDr<0.5:
                lep_genindex[i]=j
                minDr=thisDr

        if lep_genindex[i] != -1:index.remove(lep_genindex[i])

    return lep_genindex


# Gen-Reco lepton matching for failed events -> [-2,-2,-2,-2] default
def GenMatching_bis(genlep_pt, genlep_eta, genlep_phi, genlep_id):
    lep_genindex = [-2, -2, -2, -2]
    return lep_genindex


# Lepton mass
def mass_lep(flavour):
    if abs(float(flavour)) == 11: return 0.0005109989461
    elif abs(float(flavour)) == 13: return 0.1056583745
    elif abs(float(flavour)) == 15: return 1.77686
    elif abs(float(flavour)) == 0: return 0
    else: print('Houston, we have a problem', sapore)


# Lepton tetravector
def tetra(pt, eta, phi, ident):
    t0 = ROOT.TLorentzVector()
    t1 = ROOT.TLorentzVector()
    t2 = ROOT.TLorentzVector()
    t3 = ROOT.TLorentzVector()
    mass = [mass_lep(i) for i in ident]

    t0.SetPtEtaPhiM(pt[0],eta[0],phi[0],mass[0])
    t1.SetPtEtaPhiM(pt[1],eta[1],phi[1],mass[1])
    t2.SetPtEtaPhiM(pt[2],eta[2],phi[2],mass[2])
    t3.SetPtEtaPhiM(pt[3],eta[3],phi[3],mass[3])

    tet = [t0, t1, t2, t3]

    return tet


# Weights for histogram
def weight(df, fail, xsec, gen, lumi, additional = None):
    #Coefficient to calculate weights for histograms
    coeff = (lumi * 1000 * xsec) / gen
    #Gen
    weight_gen = df.genHEPMCweight * df.PUWeight
    weight_histo_gen = weight_gen * coeff
    #Reco
    if(fail == False):
        weight_reco = (weight_gen * df.dataMCWeight * df.L1prefiringWeight * df.trigEffWeight)
        weight_histo_reco = weight_reco * coeff
    elif(fail == True):
        weight_reco = 0
        weight_histo_reco = weight_reco * coeff
    #Columns in pandas
    df['weight_gen'] = weight_gen #Powheg
    df['weight_reco'] = weight_reco #Powheg
    df['weight_histo_gen'] = weight_histo_gen #Powheg
    df['weight_histo_reco'] = weight_histo_reco #Powheg
    if additional == 'ggH':
        weight_gen_NNLOPS = weight_gen * df.ggH_NNLOPS_weight
        weight_reco_NNLOPS = weight_reco * df.ggH_NNLOPS_weight
        weight_histo_gen_NNLOPS = weight_histo_gen * df.ggH_NNLOPS_weight
        weight_histo_reco_NNLOPS = weight_histo_reco * df.ggH_NNLOPS_weight
        df['weight_gen_NNLOPS'] = weight_gen_NNLOPS #NNLOPS (only ggH)
        df['weight_reco_NNLOPS'] = weight_reco_NNLOPS #NNLOPS (only ggH)
        df['weight_histo_gen_NNLOPS'] = weight_histo_gen_NNLOPS #NNLOPS (only ggH)
        df['weight_histo_reco_NNLOPS'] = weight_histo_reco_NNLOPS #NNLOPS (only ggH)
    return df


# Uproot to generate pandas
def prepareTrees(year):
    d_sig = {}
    d_sig_failed = {}
    for signal in signals_original:
        fname = eos_path_sig + '%i' %year
        fname += '/'+signal+'/reducedTree_MC_'+str(year)+'_'+signal+'.root'
        d_sig[signal] = uproot.open(fname)[key]
        d_sig_failed[signal] = uproot.open(fname)[key_failed]

    return d_sig, d_sig_failed


# Order Pt leptons in magnitude
def add_order_lep(df):
    sortdf=pd.DataFrame(np.sort(df[['GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt']].values))
    df['GenLep1PtOrder'] = sortdf.iloc[:,-1]
    df['GenLep2PtOrder'] = sortdf.iloc[:,-2]
    df['GenLep3PtOrder'] = sortdf.iloc[:,-3]
    df['GenLep4PtOrder'] = sortdf.iloc[:,-4]
    return df


# Calculate cross sections
def xsecs(year):
    xsec_sig = {}
    d_sig, d_sig_failed = prepareTrees(year)
    for signal in signals_original:
        xsec_sig[signal] = d_sig[signal].pandas.df('xsec').xsec[0]
    return xsec_sig



# RAPIDITY
def rapidity(p, eta):
    return np.log((np.sqrt(125*125 + p*p*np.cosh(eta)*np.cosh(eta))+p*np.sinh(eta))/np.sqrt(125*125+p*p))
def add_rapidity(df):
    df['ZZy'] = rapidity(df['ZZPt'], df['ZZEta'])
    return df


def add_fin_state(i, j):
    if abs(i) == 121 and abs(j) == 121:
        fin = '4e'
    elif abs(i) == 169 and abs(j) == 169:
        fin = '4mu'
    elif (abs(i) == 121 and abs(j) == 169) or (abs(i) == 169 and abs(j) == 121):
        fin = '2e2mu'
    elif (abs(i) == 225 and abs(j) == 169) or (abs(i) == 169 and abs(j) == 225):
        fin = '2tau2mu'
    elif (abs(i) == 225 and abs(j) == 121) or (abs(i) == 121 and abs(j) == 225):
        fin = '2tau2e'
    elif abs(i) == 225 and abs(j) == 225:
        fin = '4tau'
    elif abs(i) == 0 and abs(j) == 0:
        fin = 'other'
    return fin


# Get the "number" of MC events to divide the weights
def generators(year):
    gen_sig = {}
    for signal in signals_original:
        fname = eos_path_sig + '%i' %year
        fname += '/'+signal+'/reducedTree_MC_'+str(year)+'_'+signal+'.root'
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("Counters")
        gen_sig[signal] = hCounters.GetBinContent(40)

    return gen_sig


# Set up data frames
def dataframes(year):
    if year == 2016:
        lumi = 35.9
    elif year == 2017:
        lumi = 41.5
    elif year == 2018:
        lumi = 59.7

    d_df_sig = {}
    d_df_sig_failed = {}
    d_sig, d_sig_failed = prepareTrees(year)
    gen_sig = generators(year)
    xsec_sig = xsecs(year)

    # PASSED
    for signal in signals_original:
        print('Processing', signal, year)
        b_sig = ['ZZMass', 'ZZPt', 'Z1Mass', 'Z2Mass', 'ZZEta', 'Z1Flav', 'Z2Flav',
                 'GenHMass', 'GenHPt', 'GenHRapidity', 'GenZ1Flav', 'GenZ2Flav','GenZ1Mass', 'GenZ2Mass',
                 'GenZ1Pt', 'GenZ2Pt',
                 'GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt',
                 'GenLep1Eta', 'GenLep2Eta', 'GenLep3Eta', 'GenLep4Eta',
                 'GenLep1Phi', 'GenLep2Phi', 'GenLep3Phi', 'GenLep4Phi',
                 'GenLep1Id', 'GenLep2Id', 'GenLep3Id', 'GenLep4Id',
                 'GenLepPtSorted', 'GenLepEtaSorted', 'GenLepPhiSorted', 'GenLepIdSorted',
                 'GenAssocLep1Id', 'GenAssocLep2Id', 'GenAssocLep1Pt', 'GenAssocLep2Pt', 'passedFiducialSelection',
                 'LepPt', 'LepEta', 'LepPhi', 'LepLepId',
                 'PUWeight', 'genHEPMCweight', 'overallEventWeight', 'L1prefiringWeight',
                 'dataMCWeight', 'trigEffWeight']
        if signal == 'ggH125':
            b_sig.append('ggH_NNLOPS_weight') #Additional entry for the weight in case of ggH
        gen = gen_sig[signal]
        xsec = xsec_sig[signal]
        df = d_sig[signal].pandas.df(b_sig, flatten = False)
        df['gen'] = gen
        df['xsec'] = xsec
        df = add_rapidity(df)
        df['FinState_reco'] = [add_fin_state(i, j) for i,j in zip(df.Z1Flav, df.Z2Flav)]
        df['FinState_gen'] = [add_fin_state(i, j) for i,j in zip(df.GenZ1Flav, df.GenZ2Flav)]
#         df = add_order_lep(df)
        df['GenLepPt'] = df[['GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt']].values.tolist()
        df['GenLepEta'] = df[['GenLep1Eta', 'GenLep2Eta', 'GenLep3Eta', 'GenLep4Eta']].values.tolist()
        df['GenLepPhi'] = df[['GenLep1Phi', 'GenLep2Phi', 'GenLep3Phi', 'GenLep4Phi']].values.tolist()
        df['GenLepId'] = df[['GenLep1Id', 'GenLep2Id', 'GenLep3Id', 'GenLep4Id']].values.tolist()
#         df['GenLepPtOrder'] = df[['GenLep1PtOrder', 'GenLep2PtOrder', 'GenLep3PtOrder',
#                                   'GenLep4PtOrder']].values.tolist()
#         df.drop(columns = ['GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt', 'GenLep1Eta',
#                            'GenLep2Eta', 'GenLep3Eta', 'GenLep4Eta', 'GenLep1Phi', 'GenLep2Phi', 'GenLep3Phi',
#                            'GenLep4Phi', 'GenLep1Id', 'GenLep2Id', 'GenLep3Id', 'GenLep4Id',
#                            'GenLep1PtOrder', 'GenLep2PtOrder', 'GenLep3PtOrder', 'GenLep4PtOrder'])
        df['lep_genindex'] = [GenMatching(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7])
                              for row in df[['LepPt', 'LepEta', 'LepPhi', 'LepLepId', 'GenLepPt', 'GenLepEta',
                                             'GenLepPhi', 'GenLepId']].values]
#         df['GenTetravector'] = [tetra(row[0],row[1],row[2],row[3]) for row in df[['GenLepPt', 'GenLepEta',
#                                              'GenLepPhi', 'GenLepId']].values]
        if signal != 'ggH125':
            d_df_sig[signal] = weight(df, False, xsec, gen, lumi)
        else:
            d_df_sig[signal] = weight(df, False, xsec, gen, lumi, 'ggH')
            d_df_sig[signal] = d_df_sig[signal].drop(columns=['ggH_NNLOPS_weight'])
    print('SIGNAL df CREATED, %i' %year)

    # FAILED
    for signal in signals_original:
        print('Processing', signal, year)
        b_sig = ['GenHMass', 'GenHPt', 'GenHRapidity', 'GenZ1Flav', 'GenZ2Flav','GenZ1Mass', 'GenZ2Mass',
                 'GenZ1Pt', 'GenZ2Pt',
                 'GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt',
                 'GenLep1Eta', 'GenLep2Eta', 'GenLep3Eta', 'GenLep4Eta',
                 'GenLep1Phi', 'GenLep2Phi', 'GenLep3Phi', 'GenLep4Phi',
                 'GenLep1Id', 'GenLep2Id', 'GenLep3Id', 'GenLep4Id',
                 'GenLepPtSorted', 'GenLepEtaSorted', 'GenLepPhiSorted', 'GenLepIdSorted',
                 'GenAssocLep1Id', 'GenAssocLep2Id', 'GenAssocLep1Pt', 'GenAssocLep2Pt', 'passedFiducialSelection',
                 'PUWeight', 'genHEPMCweight']
        if signal == 'ggH125':
            b_sig.append('ggH_NNLOPS_weight') #Additional entry for the weight in case of ggH
        gen = gen_sig[signal]
        xsec = xsec_sig[signal]
        df = d_sig_failed[signal].pandas.df(b_sig, flatten = False)
        df['overallEventWeight'] = -1
        df['L1prefiringWeight'] = -1
        df['dataMCWeight'] = -1
        df['trigEffWeight'] = -1
        df['gen'] = gen
        df['xsec'] = xsec
        df['ZZy'] = -1 #Negative ZZy for failed events (it is useful when creating fiducial pandas)
        df['FinState_reco'] = 'fail'
        df['FinState_gen'] = [add_fin_state(i, j) for i,j in zip(df.GenZ1Flav, df.GenZ2Flav)]
#         df = add_order_lep(df)
        df['GenLepPt'] = df[['GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt']].values.tolist()
        df['GenLepEta'] = df[['GenLep1Eta', 'GenLep2Eta', 'GenLep3Eta', 'GenLep4Eta']].values.tolist()
        df['GenLepPhi'] = df[['GenLep1Phi', 'GenLep2Phi', 'GenLep3Phi', 'GenLep4Phi']].values.tolist()
        df['GenLepId'] = df[['GenLep1Id', 'GenLep2Id', 'GenLep3Id', 'GenLep4Id']].values.tolist()
#         df['GenLepPtOrder'] = df[['GenLep1PtOrder', 'GenLep2PtOrder', 'GenLep3PtOrder',
#                                   'GenLep4PtOrder']].values.tolist()
#         df.drop(columns = ['GenLep1Pt', 'GenLep2Pt', 'GenLep3Pt', 'GenLep4Pt', 'GenLep1Eta',
#                            'GenLep2Eta', 'GenLep3Eta', 'GenLep4Eta', 'GenLep1Phi', 'GenLep2Phi', 'GenLep3Phi',
#                            'GenLep4Phi', 'GenLep1Id', 'GenLep2Id', 'GenLep3Id', 'GenLep4Id',
#                            'GenLep1PtOrder', 'GenLep2PtOrder', 'GenLep3PtOrder', 'GenLep4PtOrder'])
        df['lep_genindex'] = [GenMatching_bis(row[0],row[1],row[2],row[3])
                              for row in df[['GenLepPt', 'GenLepEta', 'GenLepPhi', 'GenLepId']].values]
#         df['GenTetravector'] = [tetra(row[0],row[1],row[2],row[3]) for row in df[['GenLepPt', 'GenLepEta',
#                                              'GenLepPhi', 'GenLepId']].values]
        df.insert(0, 'Z2Mass', -1) #Negative Z1Mass for failed events (it is useful when creating fiducial pandas)
        df.insert(0, 'Z1Mass', -1) #Negative Z2Mass for failed events (it is useful when creating fiducial pandas)
        df.insert(0, 'ZZPt', -1) #Negative ZZPt for failed events (it is useful when creating fiducial pandas)
        df.insert(0, 'ZZMass', -1) #Negative ZZMass for failed events (it is useful when creating fiducial pandas)
        if signal != 'ggH125':
            d_df_sig_failed[signal] = weight(df, True, xsec, gen, lumi)
        else:
            d_df_sig_failed[signal] = weight(df, True, xsec, gen, lumi, 'ggH')
            d_df_sig_failed[signal] = d_df_sig_failed[signal].drop(columns=['ggH_NNLOPS_weight'])
    print('SIGNAL FAILED df CREATED, %i' %year)

    return d_df_sig, d_df_sig_failed

# Merge WplusH125 and WminusH125
def skim_df(year):
    d_df_sig, d_df_sig_failed = dataframes(year)

    d_skim_sig = {}
    d_skim_sig_failed = {}

    frames = []
    for signal in signals_original:
        if (signal == 'WplusH125') or (signal == 'WminusH125'):
            frames.append(d_df_sig[signal])
        else:
            d_skim_sig[signal] = d_df_sig[signal]
    d_skim_sig['WH125'] = pd.concat(frames)

    frames = []
    for signal in signals_original:
        if (signal == 'WplusH125') or (signal == 'WminusH125'):
            frames.append(d_df_sig_failed[signal])
        else:
            d_skim_sig_failed[signal] = d_df_sig_failed[signal]
    d_skim_sig_failed['WH125'] = pd.concat(frames)
    print('%i SKIMMED df CREATED' %year)

    return d_skim_sig, d_skim_sig_failed

# ------------------------------- FUNCTIONS TO CALCULATE COEFFICIENTS ----------------------------------------------------
def getCoeff(channel, m4l_bins, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, year):
    #RecoBin limits I'm considering
    obs_reco_low = obs_bins[recobin]
    obs_reco_high = obs_bins[recobin+1]

    #GenBin limits I'm considering
    obs_gen_low = obs_bins[genbin]
    obs_gen_high = obs_bins[genbin+1]

    #Extrimities of gen area
    obs_gen_lowest = obs_bins[0]
    obs_gen_highest = obs_bins[len(obs_bins)-1]

    for signal in signals:
        processBin = signal+'_'+channel+'_'+obs_name+'_genbin'+str(genbin)+'_recobin'+str(recobin)

        # Selections
        cutm4l_reco = (d_sig_tot[year][signal]['ZZMass'] >= m4l_low) & (d_sig_tot[year][signal]['ZZMass'] <= m4l_high)
        cutobs_reco = (abs(d_sig_tot[year][signal][obs_reco]) >= obs_reco_low) & (abs(d_sig_tot[year][signal][obs_reco]) < obs_reco_high)
        PassedFullSelection = d_sig_tot[year][signal]['FinState_reco'] != 'fail'
#         cuth4l_reco = (d_sig_tot[year][signal].lep_genindex.apply(lambda x: any(item for item in [-1] if item not in x))) | ((d_sig_tot[year][signal].lep_genindex.apply(lambda x: any(item for item in [-1] if item in x))) & (d_sig_tot[year][signal].GenAssocLep1Id == 0))
#         cutnoth4l_reco = (d_sig_tot[year][signal].lep_genindex.apply(lambda x: any(item for item in [-1] if item in x))) & (d_sig_tot[year][signal].GenAssocLep1Id != 0)
        cuth4l_reco = d_sig_tot[year][signal].lep_genindex.apply(lambda x: any(item for item in [-1] if item not in x))
        cutnoth4l_reco = d_sig_tot[year][signal].lep_genindex.apply(lambda x: any(item for item in [-1] if item in x))
        PassedFiducialSelection = d_sig_tot[year][signal]['passedFiducialSelection'] == True
        NotPassedFiducialSelection = d_sig_tot[year][signal]['passedFiducialSelection'] != True
        cutchan_gen = d_sig_tot[year][signal]['FinState_gen'] == channel
        cutchan_reco = d_sig_tot[year][signal]['FinState_reco'] == channel #AT
        cutobs_gen = (abs(d_sig_tot[year][signal][obs_gen]) >= obs_gen_low) & (abs(d_sig_tot[year][signal][obs_gen]) < obs_gen_high)
        cutobs_gen_otherfid = ((abs(d_sig_tot[year][signal][obs_gen]) >= obs_gen_lowest) & (abs(d_sig_tot[year][signal][obs_gen]) < obs_gen_low)) | ((abs(d_sig_tot[year][signal][obs_gen]) >= obs_gen_high) & (abs(d_sig_tot[year][signal][obs_gen]) <= obs_gen_highest))

        # --------------- acceptance ---------------
        acc_num = d_sig_tot[year][signal][PassedFiducialSelection & cutobs_gen & cutchan_gen]['weight_histo_gen'].sum()
        acc_den = d_sig_tot[year][signal][cutchan_gen]['weight_histo_gen'].sum()
        acceptance[processBin] = acc_num/acc_den
        #Error
        acc_num = d_sig_tot[year][signal][PassedFiducialSelection & cutobs_gen & cutchan_gen]['weight_gen'].sum()
        acc_den = d_sig_tot[year][signal][cutchan_gen]['weight_gen'].sum()
        acc = acc_num/acc_den
        err_acceptance[processBin] = sqrt((acc*(1-acc))/acc_den)

        # --------------- EffRecoToFid ---------------
        eff_num = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection & cuth4l_reco &
                                       PassedFiducialSelection & cutchan_gen & cutobs_gen]['weight_histo_reco'].sum()
        eff_den = d_sig_tot[year][signal][PassedFiducialSelection & cutobs_gen & cutchan_gen]['weight_histo_gen'].sum()
        effrecotofid[processBin] = eff_num/eff_den
        if effrecotofid[processBin] == 0: effrecotofid[processBin] = 1e-06
        #Error
        eff_num = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection & cuth4l_reco &
                                       PassedFiducialSelection & cutchan_gen & cutobs_gen]['weight_reco'].sum()
        eff_den = d_sig_tot[year][signal][PassedFiducialSelection & cutobs_gen & cutchan_gen]['weight_gen'].sum()
        eff = eff_num/eff_den
        if (eff*(1-eff))/eff_den > 0:
            err_effrecotofid[processBin] = sqrt((eff*(1-eff))/eff_den)
        else:
            err_effrecotofid[processBin] = 0.00001

        # --------------- outinratio ---------------
        oir_num = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection &
                                       cuth4l_reco & NotPassedFiducialSelection & cutchan_reco]['weight_histo_reco'].sum()
        oir_den_1 = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection &
                                         cuth4l_reco & PassedFiducialSelection &
                                         cutchan_gen & cutobs_gen]['weight_histo_reco'].sum()
        oir_den_2 = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection &
                                         cuth4l_reco & PassedFiducialSelection &
                                         cutchan_gen & cutobs_gen_otherfid]['weight_histo_reco'].sum()
        outinratio[processBin] = oir_num/(oir_den_1+oir_den_2)
        #Error
        oir_num = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection &
                                       cuth4l_reco & NotPassedFiducialSelection & cutchan_reco]['weight_reco'].sum()
        oir_den_1 = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection &
                                         cuth4l_reco & PassedFiducialSelection &
                                         cutchan_gen & cutobs_gen]['weight_reco'].sum()
        oir_den_2 = d_sig_tot[year][signal][cutm4l_reco & cutobs_reco & PassedFullSelection &
                                         cuth4l_reco & PassedFiducialSelection &
                                         cutchan_gen & cutobs_gen_otherfid]['weight_reco'].sum()
        out = oir_num/(oir_den_1+oir_den_2)
        if (out*(1-out))/(oir_den_1+oir_den_2) > 0:
            err_outinratio[processBin] = sqrt((out*(1-out))/(oir_den_1+oir_den_2))
        else:
            err_outinratio[processBin] = 0.00001

        # --------------- wrongfrac ---------------
        wf_num = d_sig_tot[year][signal][PassedFullSelection & cutm4l_reco & cutnoth4l_reco]['weight_histo_reco'].sum()
        wf_den = d_sig_tot[year][signal][PassedFullSelection & cutm4l_reco]['weight_histo_reco'].sum()
        wrongfrac[processBin] = wf_num/wf_den

        # --------------- binfrac_wrongfrac ---------------
        binwf_num = d_sig_tot[year][signal][PassedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco]['weight_histo_reco'].sum()
        binwf_den = d_sig_tot[year][signal][PassedFullSelection & cutm4l_reco & cutnoth4l_reco]['weight_histo_reco'].sum()
        if binwf_den == 0: binwf_den = 1e-06 #To avoid the division by zero
        binfrac_wrongfrac[processBin] = binwf_num/binwf_den

        # --------------- numberFake ---------------
#         numberFake[processBin] = d_sig_tot[year][signal][PassedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco]['weight'].sum()
        numberFake[processBin] = d_sig_tot[year][signal][PassedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco & cutchan_reco]['weight_histo_reco'].sum()


def doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins):
    chans = ['4e', '4mu', '2e2mu']
    m4l_bins = 35
    m4l_low = 105.0
    m4l_high = 140.0

    for year in years:
        for chan in chans:
            for recobin in range(len(obs_bins)-1):
                for genbin in range(len(obs_bins)-1):
                    getCoeff(chan, m4l_bins, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, year)

        with open('../inputs/inputs_sig_'+obs_name+'_'+str(year)+'.py', 'w') as f:
            f.write('observableBins = '+str(obs_bins)+';\n')
            f.write('acc = '+str(acceptance)+' \n')
            f.write('err_acc = '+str(err_acceptance)+' \n')
            f.write('eff = '+str(effrecotofid)+' \n')
            f.write('err_eff = '+str(err_effrecotofid)+' \n')
            f.write('outinratio = '+str(outinratio)+' \n')
            f.write('err_outinratio = '+str(err_outinratio)+' \n')
            f.write('inc_wrongfrac = '+str(wrongfrac)+' \n')
            f.write('binfrac_wrongfrac = '+str(binfrac_wrongfrac)+' \n')
            f.write('number_fake = '+str(numberFake))

# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------
signals_original = ['VBFH125', 'ggH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']
signals = ['ggH125', 'VBFH125', 'WH125', 'ZH125', 'ttH125']
eos_path_sig = '/eos/user/a/atarabin/MC_samples/'
key = 'candTree'
key_failed = 'candTree_failed'
# years = [2016, 2017, 2018]
years = [2017]


# Generate dataframes
d_sig = {}
d_sig_failed = {}
for year in years:
    sig, sig_failed = skim_df(year)
    d_sig[year] = sig
    d_sig_failed[year] = sig_failed

# Create dataframe with all the events
d_sig_tot = {}
for year in years:
    d_sup = {}
    for signal in signals:
        print(year, signal)
        d_sig_sup = d_sig[year][signal].drop(columns = ['Z1Flav', 'Z2Flav', 'ZZEta', 'LepPt', 'LepPhi', 'LepEta',
                                                       'LepLepId'])
        d_sup[signal] = pd.concat([d_sig_sup, d_sig_failed[year][signal]], ignore_index=True)
    d_sig_tot[year] = d_sup

# Create dataframe FullRun2
d_sig_full = {}
for signal in signals:
    frame = [d_sig_tot[year][signal] for year in years]
    d_sig_full[signal] = pd.concat(frame, ignore_index=True)

obs_bins = [0, 0.15, 0.3, 0.6, 0.9, 1.2, 2.5]
obs_reco = 'ZZy'
obs_gen = 'GenHRapidity'
obs_name = 'rapidity4l'

wrongfrac = {}
binfrac_wrongfrac = {}
binfrac_outfrac = {}
outinratio = {}
err_outinratio = {}
effrecotofid = {}
err_effrecotofid = {}
acceptance = {}
err_acceptance = {}
numberFake = {}

doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins)
