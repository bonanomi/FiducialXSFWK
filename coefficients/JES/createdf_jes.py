import numpy as np
import pandas as pd
import uproot
from math import sqrt, log
import sys,os
import optparse
import itertools
import math
import json
import ROOT

jesNames = ['Total', 'Abs', 'Abs_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavQCD', 'HF', 'HF_year', 'RelBal', 'RelSample_year']
signals_original = ['ggH125', 'VBFH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']
eos_path_sig = '/eos/user/a/atarabin/MC_samples/'
key = 'candTree'


# ------------------------------- FUNCTIONS TO GENERATE DATAFRAMES ----------------------------------------------------
# Weights for histogram
def weight(df, fail, xsec, gen, lumi, additional = None):
    #Coefficient to calculate weights for histograms
    coeff = (lumi * 1000 * xsec) / gen
    #Reco
    if(fail == False):
        weight_reco = (df.overallEventWeight * df.L1prefiringWeight)
        weight_histo_reco = weight_reco * coeff
    elif(fail == True):
        weight_reco = 0
        weight_histo_reco = weight_reco * coeff
    #Columns in pandas
    df['weight_reco'] = weight_reco #Powheg
    df['weight_histo_reco'] = weight_histo_reco #Powheg
    if additional == 'ggH':
        weight_reco_NNLOPS = weight_reco * df.ggH_NNLOPS_weight
        weight_histo_reco_NNLOPS = weight_histo_reco * df.ggH_NNLOPS_weight
        df['weight_reco_NNLOPS'] = weight_reco_NNLOPS #NNLOPS (only ggH)
        df['weight_histo_reco_NNLOPS'] = weight_histo_reco_NNLOPS #NNLOPS (only ggH)
    else:
        df['weight_reco_NNLOPS'] = -1
        df['weight_histo_reco_NNLOPS'] = -1
    return df

# Uproot to generate pandas
def prepareTrees(year):
    d_sig = {}
    for signal in signals_original:
        # if(opt.AC==False):
        fname = eos_path_sig + '%i_MELA' %year
        # else:
	        # fname = eos_path_sig + 'AC%i' %year
        fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
        d_sig[signal] = uproot.open(fname)[key]

    return d_sig


# Calculate cross sections
def xsecs(year):
    xsec_sig = {}
    d_sig = prepareTrees(year)
    for signal in signals_original:
        xsec_sig[signal] = d_sig[signal].pandas.df('xsec').xsec[0]
    return xsec_sig


def add_fin_state_reco(i, j):
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
        # if(opt.AC==False):
        fname = eos_path_sig + '%i_MELA' %year
        # else:
	        # fname = eos_path_sig + 'AC%i' %year
        fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("Counters")
        gen_sig[signal] = hCounters.GetBinContent(40)
    return gen_sig

def add_leadjet(pt,eta,phi,mass):
	_pTj1 = 0.0
	_mj1 = 0.0
	_etaj1 = 0.0
	_phij1 = 0.0
	index = -1
	_j = ROOT.TLorentzVector()
	if len(pt) == 0:
		_j.SetPtEtaPhiM(0,0,0,0)
	else:
	    for i in range(len(pt)):
	        if (pt[i]>30 and abs(eta[i])<2.5 and pt[i] > _pTj1):
	        	_pTj1 = pt[i]
	        	index = i
	    _j.SetPtEtaPhiM(_pTj1,eta[index],phi[index],mass[index])
	return _j


def add_subleadjet(pt,eta,phi,mass,leadJet):
    _pTj2 = 0.0
    _mj2 = 0.0
    _etaj2 = 0.0
    _phij2 = 0.0
    _j = ROOT.TLorentzVector()
    for i in range(len(pt)):
        if (pt[i]>30 and abs(eta[i])<2.5 and pt[i] > _pTj2 and leadJet.Pt()-pt[i] > 0.00001):
            _pTj2 = pt[i]
            _mj2 = mass[i]
            _etaj2 = eta[i]
            _phij2 = phi[i]
    _j.SetPtEtaPhiM(_pTj2,_etaj2,_phij2,_mj2)
    return _j


def count_jets(pt,eta,phi,mass):
    n = 0
    for i in range(len(pt)):
        if pt[i]>30 and abs(eta[i])<2.5: n = n + 1
    return n

def varHiggsOneJets_jes(jet,Hmass,Heta,Hphi,Hpt):
    Higgs = ROOT.TLorentzVector()
    Higgs.SetPtEtaPhiM(Hpt,Heta,Hphi,Hmass)
    return (Higgs+jet).Pt()

def varHiggsTwoJets_jes(jet1,jet2,Hmass,Heta,Hphi,Hpt):
    Higgs = ROOT.TLorentzVector()
    Higgs.SetPtEtaPhiM(Hpt,Heta,Hphi,Hmass)
    return (Higgs+jet1+jet2).Pt()

def tetra_Higgs(mass,eta,phi,pt):
    h = ROOT.TLorentzVector()
    h.SetPtEtaPhiM(pt,eta,phi,mass)
    return h


def createDataframe(d_sig,fail,gen,xsec,signal,lumi,obs_reco,obs_reco_2nd='None'):
    b_sig = ['EventNumber', 'passedFiducialSelection_bbf', 'PUWeight', 'genHEPMCweight',
             'ZZMass', 'ZZPt','ZZEta', 'ZZPhi', 'Z1Flav', 'Z2Flav', 'JetPt', 'JetMass', 'JetEta', 'JetPhi',
             'overallEventWeight', 'L1prefiringWeight','dataMCWeight', 'trigEffWeight',
             'pTj1', 'Mj1', 'ETAj1', 'PHIj1',
             'pTj2', 'Mj2', 'ETAj2', 'PHIj2']
    if signal == 'ggH125': b_sig.append('ggH_NNLOPS_weight') #Additional entry for the weight in case of ggH
    for i in jesNames:
        b_sig.extend(['JetPt_JESUp_'+i,'JetPt_JESDown_'+i])
    if obs_reco != 'pTj1' and obs_reco != 'pTj2': b_sig.append(obs_reco)
    if (obs_reco_2nd!='None' and obs_reco_2nd != 'pTj1' and obs_reco_2nd != 'pTj2'): b_sig.append(obs_reco_2nd)

    df = d_sig.pandas.df(b_sig, flatten = False)
    df['gen'] = gen
    df['xsec'] = xsec
    # Leading jets
    for i in jesNames:
        df['j1_jesup_'+i] = [add_leadjet(row[0],row[1],row[2],row[3]) for row in df[['JetPt_JESUp_'+i,'JetEta','JetPhi','JetMass']].values]
        df['j1_jesdn_'+i] = [add_leadjet(row[0],row[1],row[2],row[3]) for row in df[['JetPt_JESDown_'+i,'JetEta','JetPhi','JetMass']].values]
    # Subleading jets
    for i in jesNames:
        df['j2_jesup_'+i] = [add_subleadjet(row[0],row[1],row[2],row[3],row[4]) for row in df[['JetPt_JESUp_'+i,'JetEta','JetPhi','JetMass','j1_jesup_'+i]].values]
        df['j2_jesdn_'+i] = [add_subleadjet(row[0],row[1],row[2],row[3],row[4]) for row in df[['JetPt_JESDown_'+i,'JetEta','JetPhi','JetMass','j1_jesdn_'+i]].values]
    df['FinState_reco'] = [add_fin_state_reco(i, j) for i,j in zip(df.Z1Flav, df.Z2Flav)]
    # Calculus of up and down variations for observables different from leading jet
    df['Higgs'] = [tetra_Higgs(row[0],row[1],row[2],row[3]) for row in df[['ZZMass', 'ZZEta', 'ZZPhi', 'ZZPt']].values]
    for i in jesNames:
        if obs_reco == 'pTj1' or obs_reco_2nd == 'pTj1':
            df['pTj1_jesup_'+i] = [x.Pt() for x in df['j1_jesup_'+i]]
            df['pTj1_jesdn_'+i] = [x.Pt() for x in df['j1_jesdn_'+i]]
        if obs_reco == 'pTj2' or obs_reco_2nd == 'pTj2':
            df['pTj2_jesup_'+i] = [x.Pt() for x in df['j2_jesup_'+i]]
            df['pTj2_jesdn_'+i] = [x.Pt() for x in df['j2_jesdn_'+i]]
        if 'njets' in obs_reco or 'njets' in obs_reco_2nd:
            df['njets_pt30_eta2p5_jesup_'+i] = [count_jets(row[0],row[1],row[2],row[3]) for row in df[['JetPt_JESUp_'+i,'JetEta','JetPhi','JetMass']].values]
            df['njets_pt30_eta2p5_jesdn_'+i] = [count_jets(row[0],row[1],row[2],row[3]) for row in df[['JetPt_JESDown_'+i,'JetEta','JetPhi','JetMass']].values]
        if obs_reco == 'mjj' or obs_reco_2nd == 'mjj':
            df['mjj_jesup_'+i] = [(j1+j2).M() for j1,j2 in zip(df['j1_jesup_'+i],df['j2_jesup_'+i])]
            df['mjj_jesdn_'+i] = [(j1+j2).M() for j1,j2 in zip(df['j1_jesdn_'+i],df['j2_jesdn_'+i])]
        if obs_reco == 'pTHj' or obs_reco_2nd == 'pTHj':
            df['pTHj_jesup_'+i] = [(H+j1).Pt() for H,j1 in zip(df['Higgs'],df['j1_jesup_'+i])]
            df['pTHj_jesdn_'+i] = [(H+j1).Pt() for H,j1 in zip(df['Higgs'],df['j1_jesdn_'+i])]
        if obs_reco == 'pTHjj' or obs_reco_2nd == 'pTHjj':
            df['pTHjj_jesup_'+i] = [(row[0]+row[1]+row[2]).Pt() for row in df[['Higgs','j1_jesup_'+i,'j2_jesup_'+i]].values]
            df['pTHjj_jesdn_'+i] = [(row[0]+row[1]+row[2]).Pt() for row in df[['Higgs','j1_jesdn_'+i,'j2_jesdn_'+i]].values]
        if obs_reco == 'mHj' or obs_reco_2nd == 'mHj':
            df['mHj_jesup_'+i] = [(H+j1).M() for H,j1 in zip(df['Higgs'],df['j1_jesup_'+i])]
            df['mHj_jesdn_'+i] = [(H+j1).M() for H,j1 in zip(df['Higgs'],df['j1_jesdn_'+i])]
        if obs_reco == 'mHjj' or obs_reco_2nd == 'mHjj':
            df['mHjj_jesup_'+i] = [(row[0]+row[1]+row[2]).M() for row in df[['Higgs','j1_jesup_'+i,'j2_jesup_'+i]].values]
            df['mHjj_jesdn_'+i] = [(row[0]+row[1]+row[2]).M() for row in df[['Higgs','j1_jesdn_'+i,'j2_jesdn_'+i]].values]

    if signal != 'ggH125':
        df = weight(df, fail, xsec, gen, lumi)
    else:
        df = weight(df, fail, xsec, gen, lumi, 'ggH')
        df = df.drop(columns=['ggH_NNLOPS_weight'])

    return df


# Set up data frames
def dataframes(year, doubleDiff, obs_reco, obs_reco_2nd):
    if year == 2016:
        lumi = 35.9
    elif year == 2017:
        lumi = 41.5
    elif year == 2018:
        lumi = 59.7
    d_df_sig = {}
    d_sig = prepareTrees(year)
    gen_sig = generators(year)
    xsec_sig = xsecs(year)
    for signal in signals_original:
        print ('Processing', signal, year)
        if doubleDiff:
            d_df_sig[signal] = createDataframe(d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi,obs_reco,obs_reco_2nd)
        else:
            d_df_sig[signal] = createDataframe(d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi,obs_reco)
        print ('Signal created')
    return d_df_sig


# Merge WplusH125 and WminusH125
def skim_df(year, doubleDiff, obs_reco, obs_reco_2nd = ''):
    d_df_sig = dataframes(year, doubleDiff, obs_reco, obs_reco_2nd)
    d_skim_sig = {}
    frames = []
    for signal in signals_original:
        if (signal == 'WplusH125') or (signal == 'WminusH125'):
            frames.append(d_df_sig[signal])
        else:
            d_skim_sig[signal] = d_df_sig[signal]
    d_skim_sig['WH125'] = pd.concat(frames)
    print ('%i SKIMMED df CREATED' %year)
    return d_skim_sig
