import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import uproot
from math import sqrt, log
import sys,os
import optparse
import itertools
import math
import ROOT
import json

print 'Welcome in RunCoefficients!'

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--verbose', action='store_true', dest='VERBOSE', default=False, help='print values')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()


# ------------------------------- FUNCTIONS TO GENERATE DATAFRAMES ----------------------------------------------------
# Weights for histogram
def weight(df, fail, xsec, gen, lumi, additional = None):
    #Coefficient to calculate weights for histograms
    coeff = (lumi * 1000 * xsec) / gen
    #Gen
    weight_gen = df.genHEPMCweight * df.PUWeight
    weight_histo_gen = weight_gen * coeff
    #Reco
    if(fail == False):
        weight_reco = (df.overallEventWeight * df.L1prefiringWeight)
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
        fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
        d_sig[signal] = uproot.open(fname)[key]
        d_sig_failed[signal] = uproot.open(fname)[key_failed]

    return d_sig, d_sig_failed


# Calculate cross sections
def xsecs(year):
    xsec_sig = {}
    d_sig, d_sig_failed = prepareTrees(year)
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


def add_fin_state_gen(lepId, Hindex, number):
    if (Hindex[0]==99) | (Hindex[1]==99) | (Hindex[2]==99) | (Hindex[3]==99):
        return 'other'
    if (abs(lepId[Hindex[0]])==11) & (abs(lepId[Hindex[2]])==11):
        fin = '4e'
    elif (abs(lepId[Hindex[0]])==13) & (abs(lepId[Hindex[2]])==13):
        fin = '4mu'
    elif ((abs(lepId[Hindex[0]])==11) & (abs(lepId[Hindex[2]])==13)) | ((abs(lepId[Hindex[0]])==13) & (abs(lepId[Hindex[2]])==11)):
        fin = '2e2mu'
    else:
        fin = 'other'
    return fin


def add_fin_state_gen_out(ZdauId,event):
    if (abs(ZdauId[0])==11) and (abs(ZdauId[1])==11):
        fin = '4e'
    elif (abs(ZdauId[0])==13) and (abs(ZdauId[1])==13):
        fin = '4mu'
    elif ((abs(ZdauId[0])==11) and (abs(ZdauId[1])==13)) or ((abs(ZdauId[0])==13) and (abs(ZdauId[1])==11)):
        fin = '2e2mu'
    else:
        fin = 'other'
    return fin


def add_fin_state_gen_out_ZH(ZdauId,momId):
    if ((abs(ZdauId[0])==11) and (abs(ZdauId[1])==11) and (momId[0]==25) and (momId[1]==25)) or ((abs(ZdauId[0])==11) and (abs(ZdauId[2])==11) and (momId[0]==25) and (momId[2]==25)) or ((abs(ZdauId[1])==11) and (abs(ZdauId[2])==11) and (momId[1]==25) and (momId[2]==25)):
        fin = '4e'
    elif ((abs(ZdauId[0])==13) and (abs(ZdauId[1])==13) and (momId[0]==25) and (momId[1]==25)) or ((abs(ZdauId[0])==13) and (abs(ZdauId[2])==13)&(momId[0]==25) and (momId[2]==25)) or ((abs(ZdauId[1])==13) and (abs(ZdauId[2])==13) and (momId[1]==25) and (momId[2]==25)):
        fin = '4mu'
    elif (momId[0]==25 and (ZdauId[0]==11 or ZdauId[0]==13) and momId[1]==25 and (ZdauId[1]==11 or ZdauId[1]==13) and (ZdauId[0]!=ZdauId[1])) or (momId[0]==25 and (ZdauId[0]==11 or ZdauId[0]==13) and momId[2]==25 and (ZdauId[2]==11 or ZdauId[2]==13) and (ZdauId[0]!=ZdauId[2])) or (momId[1]==25 and (ZdauId[1]==11 or ZdauId[1]==13) and momId[2]==25 and (ZdauId[2]==11 or ZdauId[2]==13) and (ZdauId[1]!=ZdauId[2])):
        fin = '2e2mu'
    else:
        fin = 'other'
    return fin


def add_cuth4l_gen(momMomId,Hindex):
    if (Hindex[0]==99) | (Hindex[1]==99) | (Hindex[2]==99) | (Hindex[3]==99):
        return False
    if momMomId[Hindex[0]]==25 and momMomId[Hindex[1]]==25 and momMomId[Hindex[2]]==25 and momMomId[Hindex[3]]==25:
        return True
    else:
        return False


def add_cuth4l_reco(Hindex,genIndex,momMomId,momId):
    if (Hindex[0]==99) | (Hindex[1]==99) | (Hindex[2]==99) | (Hindex[3]==99):
        return False
    if ((genIndex[Hindex[0]]>-0.5)*momMomId[max(0,genIndex[Hindex[0]])]==25) and ((genIndex[Hindex[0]]>-0.5)*momId[max(0,genIndex[Hindex[0]])]==23) and ((genIndex[Hindex[1]]>-0.5)*momMomId[max(0,genIndex[Hindex[1]])]==25) and ((genIndex[Hindex[1]]>-0.5)*momId[max(0,genIndex[Hindex[1]])]==23) and ((genIndex[Hindex[2]]>-0.5)*momMomId[max(0,genIndex[Hindex[2]])]==25) and ((genIndex[Hindex[2]]>-0.5)*momId[max(0,genIndex[Hindex[2]])]==23) and ((genIndex[Hindex[3]]>-0.5)*momMomId[max(0,genIndex[Hindex[3]])]==25) and ((genIndex[Hindex[3]]>-0.5)*momId[max(0,genIndex[Hindex[3]])]==23):
        return True
    else:
        return False


# Get the "number" of MC events to divide the weights
def generators(year):
    gen_sig = {}
    for signal in signals_original:
        fname = eos_path_sig + '%i' %year
        fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("Counters")
        gen_sig[signal] = hCounters.GetBinContent(40)
    return gen_sig


def createDataframe(d_sig,fail,gen,xsec,signal,lumi):
    b_sig = ['EventNumber','GENmass4l', 'GENpT4l', 'GENrapidity4l', 'GENeta4l',
             'GENlep_id', 'GENlep_MomId', 'GENlep_MomMomId', 'GENlep_Hindex',
             'GENZ_DaughtersId', 'GENZ_MomId', 'passedFiducialSelection_bbf',
             'PUWeight', 'genHEPMCweight','GENnjets_pt30_eta2p5',
             'GenCleanedJetPt', 'GenCleanedJetEta', 'GENpTj1', 'GENmassZ2', 'GENmassZ1',
             'GENcosThetaStar', 'GENcosTheta1','GENcosTheta2','GENPhi','GENPhi1',
             'GENpTHj']
    if signal == 'ggH125': b_sig.append('ggH_NNLOPS_weight') #Additional entry for the weight in case of ggH
    if not fail: b_sig.extend(['ZZMass', 'ZZPt', 'ZZy', 'Z1Mass', 'Z2Mass', 'ZZEta', 'Z1Flav', 'Z2Flav',
                          'lep_genindex', 'lep_Hindex', 'overallEventWeight', 'L1prefiringWeight','dataMCWeight', 'trigEffWeight', 'njets_pt30_eta2p5',
                          'njets_pt30_eta2p5_jesup', 'njets_pt30_eta2p5_jesdn', 'pTj1',
                          'costhetastar', 'helcosthetaZ1','helcosthetaZ2','helphi','phistarZ1',
                          'pTHj']) #Additioanl entries for passing events
    df = d_sig.pandas.df(b_sig, flatten = False)
    if fail: #Negative branches for failed events (it is useful when creating fiducial pandas)
        df['ZZMass'] = -1
        df['ZZPt'] = -1
        df['ZZy'] = -1
        df['Z1Mass'] = -1
        df['Z2Mass'] = -1
        df['ZZEta'] = -1
        df['Z1Flav'] = -1
        df['Z2Flav'] = -1
        df['lep_genindex'] = -1
        df['lep_Hindex'] = -1
        df['overallEventWeight'] = -1
        df['L1prefiringWeight'] = -1
        df['dataMCWeight'] = -1
        df['trigEffWeight'] = -1
        df['njets_pt30_eta2p5'] = -1
        df['njets_pt30_eta2p5_jesup'] = -1
        df['njets_pt30_eta2p5_jesdn'] = -1
        df['pTj1'] = -1
        df['costhetastar'] = -1
        df['helcosthetaZ1'] = -1
        df['helcosthetaZ2'] = -1
        df['helphi'] = -1
        df['phistarZ1'] = -1
        df['pTHj'] = -1
    df['gen'] = gen
    df['xsec'] = xsec
    if not fail:
        df['FinState_reco'] = [add_fin_state_reco(i, j) for i,j in zip(df.Z1Flav, df.Z2Flav)]
    elif fail:
        df['FinState_reco'] = 'fail'
    df['FinState_gen'] = [add_fin_state_gen(row[0],row[1],row[2]) for row in df[['GENlep_id', 'GENlep_Hindex', 'EventNumber']].values]
    if signal != 'ZH125':
        df['FinState_gen_out'] = [add_fin_state_gen_out(i,j) for i,j in zip(df.GENZ_DaughtersId,df.EventNumber)]
    else:
        df['FinState_gen_out'] = [add_fin_state_gen_out_ZH(i,j) for i,j in zip(df.GENZ_DaughtersId,df.GENZ_MomId)]
    df['cuth4l_gen'] = [add_cuth4l_gen(i,j) for i,j in zip(df.GENlep_MomMomId,df.GENlep_Hindex)]
    if not fail:
        df['cuth4l_reco'] = [add_cuth4l_reco(row[0],row[1],row[2],row[3]) for row in df[['lep_Hindex','lep_genindex','GENlep_MomMomId','GENlep_MomId']].values]
    elif fail:
        df['cuth4l_reco'] = False
    if signal != 'ggH125':
        df = weight(df, fail, xsec, gen, lumi)
    else:
        df = weight(df, fail, xsec, gen, lumi, 'ggH')
        df = df.drop(columns=['ggH_NNLOPS_weight'])

    return df


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
    for signal in signals_original:
        print 'Processing', signal, year
        d_df_sig[signal] = createDataframe(d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi)
        print 'Signal created'
        d_df_sig_failed[signal] = createDataframe(d_sig_failed[signal],True,gen_sig[signal],xsec_sig[signal],signal,lumi)
        print 'Signal failed created'
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
    print '%i SKIMMED df CREATED' %year
    return d_skim_sig, d_skim_sig_failed

# ------------------------------- FUNCTIONS TO CALCULATE COEFFICIENTS ----------------------------------------------------
def getCoeff(channel, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, type, year, obs_reco_2nd = 'None', obs_gen_2nd = 'None', obs_name_2nd = 'None'):
    if not doubleDiff:
        #RecoBin limits I'm considering
        obs_reco_low = obs_bins[recobin]
        obs_reco_high = obs_bins[recobin+1]
        #GenBin limits I'm considering
        obs_gen_low = obs_bins[genbin]
        obs_gen_high = obs_bins[genbin+1]
        #Extrimities of gen area
        obs_gen_lowest = obs_bins[0]
        obs_gen_highest = obs_bins[len(obs_bins)-1]
    elif doubleDiff:
        obs_reco_low = obs_bins[recobin][0]
        obs_reco_high = obs_bins[recobin][1]
        obs_gen_low = obs_bins[genbin][0]
        obs_gen_high = obs_bins[genbin][1]
        obs_gen_lowest = min(x[0] for x in obs_bins.values())
        obs_gen_highest = max(x[1] for x in obs_bins.values())
        #Second variable
        obs_reco_2nd_low = obs_bins[recobin][2]
        obs_reco_2nd_high = obs_bins[recobin][3]
        obs_gen_2nd_low = obs_bins[genbin][2]
        obs_gen_2nd_high = obs_bins[genbin][3]
        obs_gen_2nd_lowest = min(x[2] for x in obs_bins.values())
        obs_gen_2nd_highest = max(x[3] for x in obs_bins.values())

    for signal in signals:
        if type=='std':
            datafr = d_sig_tot[year][signal]
            genweight = 'weight_gen'
            recoweight = 'weight_reco'
        elif type=='full':
            datafr = d_sig_full[signal]
            genweight = 'weight_gen'
            recoweight = 'weight_reco'
        elif type=='fullNNLOPS' and signal=='ggH125':
            datafr = d_sig_full[signal]
            genweight = 'weight_gen_NNLOPS'
            recoweight = 'weight_reco_NNLOPS'
        elif type=='fullNNLOPS' and signal!='ggH125': # In case of fullNNLOPS we are interested in ggH125 only
            continue

        if doubleDiff:
            processBin = signal+'_'+channel+'_'+obs_name+'_'+obs_name_2nd+'_genbin'+str(genbin)+'_recobin'+str(recobin)
        else:
            processBin = signal+'_'+channel+'_'+obs_name+'_genbin'+str(genbin)+'_recobin'+str(recobin)

        if type=='fullNNLOPS':
            processBin = signal+'_NNLOPS_'+channel+'_'+obs_name+'_genbin'+str(genbin)+'_recobin'+str(recobin)
        if type=='fullNNLOPS' and doubleDiff:
            processBin = signal+'_NNLOPS_'+channel+'_'+obs_name+'_'+obs_name_2nd+'_genbin'+str(genbin)+'_recobin'+str(recobin)

        # Selections
        cutobs_reco = (abs(datafr[obs_reco]) >= obs_reco_low) & (abs(datafr[obs_reco]) < obs_reco_high)
        #cutobs_reco &= (datafr['Z2Mass'] < 60)
        cutobs_gen = (abs(datafr[obs_gen]) >= obs_gen_low) & (abs(datafr[obs_gen]) < obs_gen_high)
        if doubleDiff:
            cutobs_reco &= (abs(datafr[obs_reco_2nd]) >= obs_reco_2nd_low) & (abs(datafr[obs_reco_2nd]) < obs_reco_2nd_high)
            cutobs_gen &= (abs(datafr[obs_gen_2nd]) >= obs_gen_2nd_low) & (abs(datafr[obs_gen_2nd]) < obs_gen_2nd_high)
        #cutobs_gen &= (datafr['GENmassZ2'] < 60)

        if 'jet' in obs_name:
            cutobs_reco_jesup = (datafr[obs_reco+'_jesup'] >= obs_reco_low) & (datafr[obs_reco+'_jesup'] < obs_reco_high)
            cutobs_reco_jesdn = (datafr[obs_reco+'_jesdn'] >= obs_reco_low) & (datafr[obs_reco+'_jesdn'] < obs_reco_high)
        if 'jet' in obs_name_2nd: # If it is not doubleDiff obs_reco_2nd = 'None', this if is valid only for doubleDiff
            if 'cutobs_reco_jesup' in locals(): # If there is already cutobs_reco_jes* we need to add additional conditions (&=)
                cutobs_reco_jesup &= (datafr[obs_reco_2nd+'_jesup'] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd+'_jesup'] < obs_reco_2nd_high)
                cutobs_reco_jesdn &= (datafr[obs_reco_2nd+'_jesdn'] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd+'_jesdn'] < obs_reco_2nd_high)
            else: # Otherwise it is a first declaration (=)
                cutobs_reco_jesup = (datafr[obs_reco_2nd+'_jesup'] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd+'_jesup'] < obs_reco_2nd_high)
                cutobs_reco_jesdn = (datafr[obs_reco_2nd+'_jesdn'] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd+'_jesdn'] < obs_reco_2nd_high)
        cutobs_gen_otherfid = ((abs(datafr[obs_gen]) >= obs_gen_lowest) & (abs(datafr[obs_gen]) < obs_gen_low)) | ((abs(datafr[obs_gen]) >= obs_gen_high) & (abs(datafr[obs_gen]) <= obs_gen_highest))
        if doubleDiff:
            cutobs_gen_otherfid |= ((abs(datafr[obs_gen_2nd]) >= obs_gen_2nd_lowest) & (abs(datafr[obs_gen_2nd]) < obs_gen_2nd_low)) | ((abs(datafr[obs_gen_2nd]) >= obs_gen_2nd_high) & (abs(datafr[obs_gen_2nd]) <= obs_gen_2nd_highest))
        cutm4l_gen = (datafr['GENmass4l'] > m4l_low) & (datafr['GENmass4l'] < m4l_high)
        cutnotm4l_gen = (datafr['GENmass4l'] <= m4l_low) | (datafr['GENmass4l'] >= m4l_high)
        cuth4l_gen = datafr['cuth4l_gen'] == True
        cutnoth4l_gen = datafr['cuth4l_gen'] == False
        cuth4l_reco = datafr['cuth4l_reco'] == True
        cutnoth4l_reco = datafr['cuth4l_reco'] == False
        passedFullSelection = datafr['FinState_reco'] != 'fail'
        passedFiducialSelection = datafr['passedFiducialSelection_bbf'] == True
        notPassedFiducialSelection = datafr['passedFiducialSelection_bbf'] == False
        if channel != '4l':
            cutm4l_reco = (datafr['ZZMass'] > m4l_low) & (datafr['ZZMass'] < m4l_high) & (datafr['FinState_reco'] == channel)
            cutchan_gen = datafr['FinState_gen'] == channel
            cutchan_gen_out = datafr['FinState_gen_out'] == channel
        else:
            cutm4l_reco = (datafr['ZZMass'] > m4l_low) & (datafr['ZZMass'] < m4l_high)
            cutchan_gen = (datafr['FinState_gen'] == '2e2mu') | (datafr['FinState_gen'] == '4e') | (datafr['FinState_gen'] == '4mu')
            cutchan_gen_out = (datafr['FinState_gen_out'] == '2e2mu') | (datafr['FinState_gen_out'] == '4e') | (datafr['FinState_gen_out'] == '4mu')


        # --------------- acceptance ---------------
        acc_num = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen][genweight].sum()
        acc_den = datafr[cutchan_gen_out][genweight].sum()
        if acc_den>0:
            acceptance[processBin] = acc_num/acc_den
            err_acceptance[processBin] = sqrt((acceptance[processBin]*(1-acceptance[processBin]))/acc_den)
        else:
            acceptance[processBin] = -1.0
            err_acceptance[processBin] = -1.0


        if type=='fullNNLOPS': continue # In case of fullNNLOPS we are interested in acceptance only

        # --------------- EffRecoToFid ---------------
        eff_num = datafr[cutm4l_reco & cutobs_reco & passedFullSelection & cuth4l_reco &
                                       passedFiducialSelection & cuth4l_gen & cutm4l_gen & cutchan_gen & cutobs_gen][recoweight].sum()
        eff_den = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen][genweight].sum()
        if eff_den>10:
            effrecotofid[processBin] = eff_num/eff_den
            if effrecotofid[processBin] == 0: effrecotofid[processBin] = 1e-06
            if (effrecotofid[processBin]*(1-effrecotofid[processBin]))/eff_den > 0:
                err_effrecotofid[processBin] = sqrt((effrecotofid[processBin]*(1-effrecotofid[processBin]))/eff_den)
            else:
                err_effrecotofid[processBin] = 1e-06
        else:
            effrecotofid[processBin] = -1.0
            err_effrecotofid[processBin] = -1.0

        # --------------- outinratio ---------------
        oir_num = datafr[cutm4l_reco & cutobs_reco & passedFullSelection &
                                       cuth4l_reco & cutchan_gen_out & (notPassedFiducialSelection | cutnoth4l_gen | cutnotm4l_gen)][recoweight].sum()
        oir_den_1 = datafr[cutm4l_reco & cutobs_reco & passedFullSelection &
                                         cuth4l_reco & passedFiducialSelection & cuth4l_gen & cutm4l_gen &
                                         cutchan_gen & cutobs_gen][recoweight].sum()
        oir_den_2 = datafr[cutm4l_reco & cutobs_reco & passedFullSelection &
                                         cuth4l_reco & passedFiducialSelection & cuth4l_gen & cutm4l_gen &
                                         cutchan_gen & cutobs_gen_otherfid][recoweight].sum()
        if oir_den_1+oir_den_2>0:
            outinratio[processBin] = oir_num/(oir_den_1+oir_den_2)
            if (outinratio[processBin]*(1-outinratio[processBin]))/(oir_den_1+oir_den_2) > 0:
                err_outinratio[processBin] = sqrt((outinratio[processBin]*(1-outinratio[processBin]))/(oir_den_1+oir_den_2))
            else:
                err_outinratio[processBin] = 1e-06
        else:
            outinratio[processBin] = 0.0
            err_outinratio[processBin] = 0.0

        # --------------- wrongfrac ---------------
        wf_num = datafr[passedFullSelection & cutm4l_reco & cutnoth4l_reco][recoweight].sum()
        wf_den = datafr[passedFullSelection & cutm4l_reco][recoweight].sum()
        if wf_den>0:
            wrongfrac[processBin] = wf_num/wf_den
        else:
            wrongfrac[processBin] = -1.0

        # --------------- binfrac_wrongfrac ---------------
        binwf_num = datafr[passedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco][recoweight].sum()
        binwf_den = datafr[passedFullSelection & cutm4l_reco & cutnoth4l_reco][recoweight].sum()
        if binwf_den>0:
            binfrac_wrongfrac[processBin] = binwf_num/binwf_den
        else:
            binfrac_wrongfrac[processBin] = -1.0

        # --------------- numberFake ---------------
#         numberFake[processBin] = datafr[PassedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco & cutchan_reco]['weight_histo_reco'].sum()
        numberFake[processBin] = -1

        if 'jet' in obs_name:
            # --------------- lambdajesup ---------------
            lambdajesup_num_1 = datafr[passedFullSelection & cutm4l_reco & cutobs_reco_jesup & cuth4l_reco][recoweight].sum()
            lambdajesup_num_2 = lambdajesup_den = datafr[passedFullSelection & cutm4l_reco & cutobs_reco & cuth4l_reco][recoweight].sum()
            if lambdajesup_den>0:
                lambdajesup[processBin] = (lambdajesup_num_1 - lambdajesup_num_2)/lambdajesup_den
            else:
                lambdajesup[processBin] = 0.0
            # --------------- lambdajesdn ---------------
            lambdajesdn_num_1 = datafr[passedFullSelection & cutm4l_reco & cutobs_reco_jesdn & cuth4l_reco][recoweight].sum()
            lambdajesdn_num_2 = lambdajesdn_den = datafr[passedFullSelection & cutm4l_reco & cutobs_reco & cuth4l_reco][recoweight].sum()
            if lambdajesdn_den>0:
                lambdajesdn[processBin] = (lambdajesdn_num_1 - lambdajesdn_num_2)/lambdajesdn_den
            else:
                lambdajesdn[processBin] = 0.0
        else:
            lambdajesup[processBin] = 0.0
            lambdajesdn[processBin] = 0.0

        if opt.VERBOSE:
            print processBin,'acc',round(acceptance[processBin],4),'eff',round(effrecotofid[processBin],4),'outinratio',round(outinratio[processBin],4), '\n'



def doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, type, obs_reco_2nd = 'None', obs_gen_2nd = 'None', obs_name_2nd = 'None',):
    if obs_reco != 'ZZMass':
        chans = ['4e', '4mu', '2e2mu']
    else:
        chans = ['4l', '4e', '4mu', '2e2mu']
    m4l_low = 105.0
    m4l_high = 140.0

    nBins = len(obs_bins)
    if not doubleDiff: nBins = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)

    if type=='std':
        for year in years:
            for chan in chans:
                for recobin in range(nBins):
                    for genbin in range(nBins):
                        getCoeff(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, type, year, obs_reco_2nd, obs_gen_2nd, obs_name_2nd)
            # Write dictionaries
            if doubleDiff: obs_name_dic = obs_name+'_'+obs_name_2nd
            else: obs_name_dic = obs_name
            if (os.path.exists('../inputs/inputs_sig_'+obs_name_dic+'_'+str(year)+'_ORIG.py')):
                os.system('rm ../inputs/inputs_sig_'+obs_name_dic+'_'+str(year)+'_ORIG.py')
            with open('../inputs/inputs_sig_'+obs_name_dic+'_'+str(year)+'.py', 'w') as f:
                f.write('observableBins = '+str(obs_bins)+';\n')
                f.write('acc = '+str(acceptance)+' \n')
                f.write('err_acc = '+str(err_acceptance)+' \n')
                f.write('eff = '+str(effrecotofid)+' \n')
                f.write('err_eff = '+str(err_effrecotofid)+' \n')
                f.write('outinratio = '+str(outinratio)+' \n')
                f.write('err_outinratio = '+str(err_outinratio)+' \n')
                f.write('inc_wrongfrac = '+str(wrongfrac)+' \n')
                f.write('binfrac_wrongfrac = '+str(binfrac_wrongfrac)+' \n')
                f.write('number_fake = '+str(numberFake)+' \n')
                f.write('lambdajesup = '+str(lambdajesup)+' \n')
                f.write('lambdajesdn = '+str(lambdajesup))

    elif type=='full' or type=='fullNNLOPS':
        for chan in chans:
            for recobin in range(nBins):
                for genbin in range(nBins):
                    getCoeff(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, type, 'None', obs_reco_2nd, obs_gen_2nd, obs_name_2nd)

        # Write dictionaries
        if doubleDiff: obs_name_dic = obs_name+'_'+obs_name_2nd
        else: obs_name_dic = obs_name
        if type=='full':
            if (os.path.exists('../inputs/inputs_sig_'+obs_name_dic+'_Full_ORIG.py')):
                os.system('rm ../inputs/inputs_sig_'+obs_name_dic+'_Full_ORIG.py')
            with open('../inputs/inputs_sig_'+obs_name_dic+'_Full.py', 'w') as f:
                f.write('observableBins = '+str(obs_bins)+';\n')
                f.write('acc = '+str(acceptance)+' \n')
                f.write('err_acc = '+str(err_acceptance)+' \n')
                f.write('eff = '+str(effrecotofid)+' \n')
                f.write('err_eff = '+str(err_effrecotofid)+' \n')
                f.write('outinratio = '+str(outinratio)+' \n')
                f.write('err_outinratio = '+str(err_outinratio)+' \n')
                f.write('inc_wrongfrac = '+str(wrongfrac)+' \n')
                f.write('binfrac_wrongfrac = '+str(binfrac_wrongfrac)+' \n')
                f.write('number_fake = '+str(numberFake)+' \n')
                f.write('lambdajesup = '+str(lambdajesup)+' \n')
                f.write('lambdajesdn = '+str(lambdajesup))
        elif type=='fullNNLOPS':
            if(opt.YEAR == 'Full'):
                with open('../inputs/inputs_sig_'+obs_name_dic+'_NNLOPS_Full.py', 'w') as f:
                    f.write('observableBins = '+str(obs_bins)+';\n')
                    f.write('acc = '+str(acceptance)+' \n')
                    f.write('err_acc = '+str(err_acceptance)+' \n')
            else:
                with open('../inputs/inputs_sig_'+obs_name_dic+'_NNLOPS_'+opt.YEAR+'.py', 'w') as f:
                    f.write('observableBins = '+str(obs_bins)+';\n')
                    f.write('acc = '+str(acceptance)+' \n')
                    f.write('err_acc = '+str(err_acceptance)+' \n')

# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------
signals_original = ['VBFH125', 'ggH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']
signals = ['ggH125', 'VBFH125', 'WH125', 'ZH125', 'ttH125']
eos_path_sig = '/eos/user/a/atarabin/MC_samples/'
key = 'candTree'
key_failed = 'candTree_failed'

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

if doubleDiff:
    obs_name = opt.OBSNAME.split(' vs ')[0]
    obs_name_2nd = opt.OBSNAME.split(' vs ')[1]
else:
    obs_name = opt.OBSNAME

if(opt.OBSNAME == 'rapidity4l'):
    obs_reco = 'ZZy'
    obs_gen = 'GENrapidity4l'
elif(opt.OBSNAME == 'pT4l'):
    obs_reco = 'ZZPt'
    obs_gen = 'GENpT4l'
elif(opt.OBSNAME == 'massZ1'):
    obs_reco = 'Z1Mass'
    obs_gen = 'GENmassZ1'
elif(opt.OBSNAME == 'massZ2'):
    obs_reco = 'Z2Mass'
    obs_gen = 'GENmassZ2'
elif(opt.OBSNAME == 'mass4l'):
    obs_reco = 'ZZMass'
    obs_gen = 'GENmass4l'
elif(opt.OBSNAME == "njets_pt30_eta2p5"):
    obs_reco = "njets_pt30_eta2p5"
    obs_gen = "GENnjets_pt30_eta2p5"
elif(opt.OBSNAME == 'pTj1'):
    obs_reco = 'pTj1'
    obs_gen = 'GENpTj1'
elif(opt.OBSNAME == 'mass4l'):
    obs_reco = 'ZZMass'
    obs_gen = 'GENmass4l'
elif(opt.OBSNAME == 'costhetastar'):
    obs_reco = 'costhetastar'
    obs_gen = 'GENcosThetaStar'
elif(opt.OBSNAME == 'costhetaZ1'):
    obs_reco = 'helcosthetaZ1'
    obs_gen  = 'GENcosTheta1'
elif(opt.OBSNAME == 'costhetaZ2'):
    obs_reco = 'helcosthetaZ2'
    obs_gen  = 'GENcosTheta2'
elif(opt.OBSNAME == 'phi'):
    obs_reco = 'helphi'
    obs_gen  = 'GENPhi'
elif(opt.OBSNAME == 'phistar'):
    obs_reco = 'phistarZ1'
    obs_gen  = 'GENPhi1'
elif(opt.OBSNAME == 'njets_pt30_eta2p5 vs pT4l'):
    obs_reco = 'njets_pt30_eta2p5'
    obs_reco_2nd = 'ZZPt'
    obs_gen = 'GENnjets_pt30_eta2p5'
    obs_gen_2nd = 'GENpT4l'
elif(opt.OBSNAME == 'massZ1 vs massZ2'):
    obs_reco = 'Z1Mass'
    obs_reco_2nd = 'Z2Mass'
    obs_gen = 'GENmassZ1'
    obs_gen_2nd = 'GENmassZ2'
elif(opt.OBSNAME == 'njets_pt30_eta2p5 vs pTHj'):
    obs_reco = 'njets_pt30_eta2p5'
    obs_reco_2nd = 'pTHj'
    obs_gen = 'GENnjets_pt30_eta2p5'
    obs_gen_2nd = 'GENpTHj'


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
        print year, signal
        d_sup[signal] = pd.concat([d_sig[year][signal], d_sig_failed[year][signal]], ignore_index=True, sort=True)
    d_sig_tot[year] = d_sup


# Create dataframe FullRun2
if(opt.YEAR == 'Full'):
    d_sig_full = {}
    for signal in signals:
        frame = [d_sig_tot[year][signal] for year in years]
        d_sig_full[signal] = pd.concat(frame, ignore_index=True, sort=True)
else: # If I work with one year only, the FullRun2 df coincides with d_sig_tot (it is useful when fullNNLOPS is calculated)
    d_sig_full = d_sig_tot[int(opt.YEAR)]
print 'Dataframes created successfully'

print 'Coeff std'
wrongfrac = {}
binfrac_wrongfrac = {}
binfrac_outfrac = {}
outinratio = {}
err_outinratio = {}
effrecotofid = {}
err_effrecotofid = {}
acceptance = {}
err_acceptance = {}
lambdajesup = {}
lambdajesdn = {}
numberFake = {}
if doubleDiff:
    doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, 'std', obs_reco_2nd, obs_gen_2nd, obs_name_2nd)
else:
    doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, 'std')

if (opt.YEAR == 'Full'):
    print 'Coeff full'
    wrongfrac = {}
    binfrac_wrongfrac = {}
    binfrac_outfrac = {}
    outinratio = {}
    effrecotofid = {}
    err_effrecotofid = {}
    acceptance = {}
    err_acceptance = {}
    numberFake = {}
    lambdajesup = {}
    lambdajesdn = {}
    if doubleDiff:
        doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, 'full', obs_reco_2nd, obs_gen_2nd, obs_name_2nd)
    else:
        doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, 'full')

print 'Coeff fullNNLOPS'
acceptance = {}
err_acceptance = {}
if doubleDiff:
    doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, 'fullNNLOPS', obs_reco_2nd, obs_gen_2nd, obs_name_2nd)
else:
    doGetCoeff(obs_reco, obs_gen, obs_name, obs_bins, 'fullNNLOPS')
