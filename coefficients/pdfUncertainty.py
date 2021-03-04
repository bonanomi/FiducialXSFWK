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

print 'Welcome in pdfUncertainty!'

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

    # if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
    #     parser.error('Bin boundaries not specified for differential measurement. Exiting...')
    #     sys.exit()


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
    #Specific for pdfUncertainty
    df['weight_gen_NNLO_0'] = weight_gen * df.LHEweight_QCDscale_muR1_muF1
    df['weight_gen_NNLO_1'] = weight_gen * df.LHEweight_QCDscale_muR1_muF2
    df['weight_gen_NNLO_2'] = weight_gen * df.LHEweight_QCDscale_muR1_muF0p5
    df['weight_gen_NNLO_3'] = weight_gen * df.LHEweight_QCDscale_muR2_muF1
    df['weight_gen_NNLO_4'] = weight_gen * df.LHEweight_QCDscale_muR2_muF2
    df['weight_gen_NNLO_5'] = weight_gen * df.LHEweight_QCDscale_muR2_muF0p5
    df['weight_gen_NNLO_6'] = weight_gen * df.LHEweight_QCDscale_muR0p5_muF1
    df['weight_gen_NNLO_7'] = weight_gen * df.LHEweight_QCDscale_muR0p5_muF2
    df['weight_gen_NNLO_8'] = weight_gen * df.LHEweight_QCDscale_muR0p5_muF0p5
    df['weight_gen_PDFup'] = weight_gen * df.LHEweight_PDFVariation_Up
    df['weight_gen_PDFdn'] = weight_gen * df.LHEweight_PDFVariation_Dn
    if additional == 'ggH':
        weight_gen_NNLOPS = weight_gen * df.ggH_NNLOPS_weight
        weight_reco_NNLOPS = weight_reco * df.ggH_NNLOPS_weight
        weight_histo_gen_NNLOPS = weight_histo_gen * df.ggH_NNLOPS_weight
        weight_histo_reco_NNLOPS = weight_histo_reco * df.ggH_NNLOPS_weight
        df['weight_gen_NNLOPS'] = weight_gen_NNLOPS #NNLOPS (only ggH)
        df['weight_reco_NNLOPS'] = weight_reco_NNLOPS #NNLOPS (only ggH)
        df['weight_histo_gen_NNLOPS'] = weight_histo_gen_NNLOPS #NNLOPS (only ggH)
        df['weight_histo_reco_NNLOPS'] = weight_histo_reco_NNLOPS #NNLOPS (only ggH)
        #Specific for pdfUncertainty
        #In case of NNLOPS I should use always nnloWeights[0], for coding purposes I set all weight_gen_NNLOPS_NNLO_* equal to weight_gen_NNLOPS_NNLO_0
        df['weight_gen_NNLOPS_NNLO_0'] = weight_gen_NNLOPS * df.LHEweight_QCDscale_muR1_muF1
        df['weight_gen_NNLOPS_NNLO_1'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_2'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_3'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_4'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_5'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_6'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_7'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_NNLO_8'] = df['weight_gen_NNLOPS_NNLO_0']
        df['weight_gen_NNLOPS_PDFup'] = weight_gen_NNLOPS * df.LHEweight_PDFVariation_Up
        df['weight_gen_NNLOPS_PDFdn'] = weight_gen_NNLOPS * df.LHEweight_PDFVariation_Dn
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
             'GenCleanedJetPt', 'GenCleanedJetEta', 'GENpTj1',
             'LHEweight_QCDscale_muR1_muF1','LHEweight_QCDscale_muR1_muF2','LHEweight_QCDscale_muR1_muF0p5',
             'LHEweight_QCDscale_muR2_muF1','LHEweight_QCDscale_muR2_muF2','LHEweight_QCDscale_muR2_muF0p5',
             'LHEweight_QCDscale_muR0p5_muF1','LHEweight_QCDscale_muR0p5_muF2','LHEweight_QCDscale_muR0p5_muF0p5',
             'LHEweight_PDFVariation_Up', 'LHEweight_PDFVariation_Dn',
             'GENmassZ1', 'GENmassZ2',
             'GENcosThetaStar','GENcosTheta1','GENcosTheta2','GENPhi','GENPhi1',
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
    # d_skim_sig['WH125'] = pd.concat(frames)
    frames = []
    for signal in signals_original:
        if (signal == 'WplusH125') or (signal == 'WminusH125'):
            frames.append(d_df_sig_failed[signal])
        else:
            d_skim_sig_failed[signal] = d_df_sig_failed[signal]
    # d_skim_sig_failed['WH125'] = pd.concat(frames)
    print '%i SKIMMED df CREATED' %year
    return d_skim_sig, d_skim_sig_failed


# ------------------------------- FUNCTIONS TO CALCULATE COEFFICIENTS ----------------------------------------------------
def getPdfUncert(channel, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin, obs_name, type, year='None'):
    if not doubleDiff:
        #GenBin limits I'm considering
        obs_gen_low = obs_bins[genbin]
        obs_gen_high = obs_bins[genbin+1]
        #Extrimities of gen area
        obs_gen_lowest = obs_bins[0]
        obs_gen_highest = obs_bins[lenObsBins-1]
    else:
        obs_gen_low = obs_bins[genbin][0]
        obs_gen_high = obs_bins[genbin][1]
        obs_gen_lowest = min(x[0] for x in obs_bins.values())
        obs_gen_highest = max(x[1] for x in obs_bins.values())
        #Second variable
        obs_gen_2nd_low = obs_bins[genbin][2]
        obs_gen_2nd_high = obs_bins[genbin][3]
        obs_gen_2nd_lowest = min(x[2] for x in obs_bins.values())
        obs_gen_2nd_highest = max(x[3] for x in obs_bins.values())

    for signal in signals:
        datafr = d_sig_tot[year][signal]
        if type=='std':
            nnlops = ''
            processBin = signal+'_'+channel+'_'+obs_name+'_genbin'+str(genbin)
        elif type=='NNLOPS':
            nnlops = '_NNLOPS'
            processBin = signal+'_NNLOPS_'+channel+'_'+obs_name+'_genbin'+str(genbin)

        # Selections
        cutm4l_gen = (datafr['GENmass4l'] > m4l_low) & (datafr['GENmass4l'] < m4l_high)
        if obs_reco.startswith('njets') and not doubleDiff:
            cutobs_gen = abs(datafr[obs_gen]) >= obs_gen_low
        else:
            cutobs_gen = (abs(datafr[obs_gen]) >= obs_gen_low) & (abs(datafr[obs_gen]) < obs_gen_high)
            if doubleDiff:
                cutobs_gen &= (abs(datafr[obs_gen_2nd]) >= obs_gen_2nd_low) & (abs(datafr[obs_gen_2nd]) < obs_gen_2nd_high)
        if channel != '4l':
            cutm4l_reco = (datafr['ZZMass'] > m4l_low) & (datafr['ZZMass'] < m4l_high) & (datafr['FinState_reco'] == channel)
            cutchan_gen = datafr['FinState_gen'] == channel
            cutchan_gen_out = datafr['FinState_gen_out'] == channel
        else:
            cutm4l_reco = (datafr['ZZMass'] > m4l_low) & (datafr['ZZMass'] < m4l_high)
            cutchan_gen = (datafr['FinState_gen'] == '2e2mu') | (datafr['FinState_gen'] == '4e') | (datafr['FinState_gen'] == '4mu')
            cutchan_gen_out = (datafr['FinState_gen_out'] == '2e2mu') | (datafr['FinState_gen_out'] == '4e') | (datafr['FinState_gen_out'] == '4mu')
        cuth4l_gen = datafr['cuth4l_gen'] == True
        cutnoth4l_gen = datafr['cuth4l_gen'] == False
        passedFullSelection = datafr['FinState_reco'] != 'fail'
        passedFiducialSelection = datafr['passedFiducialSelection_bbf'] == True
        notPassedFiducialSelection = datafr['passedFiducialSelection_bbf'] == False

        coefficients = {}
        coefficients['fs'] = datafr[cutchan_gen_out]['weight_gen'+nnlops+'_NNLO_0'].sum()
        for i in range(0,9):
            if i==5 or i==7: continue
            coefficients['fs'+str(i)] = datafr[cutchan_gen_out]['weight_gen'+nnlops+'_NNLO_'+str(i)].sum()
            # coefficients['fid'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen'+nnlops+'_NNLO_'+str(i)].sum()
            coefficients['fid'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen_NNLO_'+str(i)].sum()
            coefficients['fid'+str(i)] = coefficients['fid'+str(i)] * (1/coefficients['fs']) #This stands for "Scale" in the original code
            # coefficients['fidraw'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen'+nnlops+'_NNLO_'+str(i)].sum()
            coefficients['fidraw'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen_NNLO_'+str(i)].sum()
            coefficients['fidraw'+str(i)] = coefficients['fidraw'+str(i)] * (1/coefficients['fs'+str(i)]) #This stands for "Scale" in the original code
            coefficients['fs'+str(i)] = coefficients['fs'+str(i)] * (1/coefficients['fs'+str(i)]) #This stands for "Scale" in the original code

        coefficients['fidPDF_up'] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen'+nnlops+'_PDFup'].sum()
        coefficients['fidPDF_up'] = coefficients['fidPDF_up'] * (1/coefficients['fs']) #This stands for "Scale" in the original code
        coefficients['fidPDF_dn'] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen'+nnlops+'_PDFdn'].sum()
        coefficients['fidPDF_dn'] = coefficients['fidPDF_dn'] * (1/coefficients['fs']) #This stands for "Scale" in the original code
        fsintegral = coefficients['fs']
        coefficients['fs'] = coefficients['fs'] * (1/coefficients['fs'])

        if coefficients['fs']>0:
            if verbose: print coefficients['fs'],coefficients['fid0']
            acceptance[processBin] = coefficients['fid0']/coefficients['fs']
            qcderrup=1.0; qcderrdn=1.0;
            accerrup=1.0; accerrdn=1.0;
            if verbose: print processBin,coefficients['fid0']
            for i in range(0,9):
                if i==5 or i==7: continue

                ratio = coefficients['fid'+str(i)]/coefficients['fid0']
                if verbose: print i, 'ratio', ratio
                if ratio>qcderrup: qcderrup = coefficients['fid'+str(i)]/coefficients['fid0']
                if ratio<qcderrdn: qcderrdn = coefficients['fid'+str(i)]/coefficients['fid0']

                acci = coefficients['fidraw'+str(i)]/coefficients['fs'+str(i)]
                if verbose: print i, 'acc', acci
                if verbose: print coefficients['fidraw'+str(i)], coefficients['fs'+str(i)]
                if acci/acceptance[processBin]>accerrup: accerrup=acci/acceptance[processBin]
                if acci/acceptance[processBin]<accerrdn: accerrdn=acci/acceptance[processBin]
            qcdUncert[processBin]={'uncerDn':abs(qcderrdn-1.0),'uncerUp':abs(qcderrup-1.0)}
            pdferr_up = coefficients['fidPDF_up']/coefficients['fid0']
            pdferr_dn = coefficients['fidPDF_dn']/coefficients['fid0']
            pdfUncert[processBin] = {"uncerDn":abs(pdferr_dn-1.0),"uncerUp":abs(pdferr_up-1.0)}

            if verbose: print(processBin,acceptance[processBin],0.0,qcderrup,qcderrdn,pdferr_up,pdferr_dn)

# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------
signals_original = signals = ['ggH125']
eos_path_sig = '/eos/user/a/atarabin/MC_samples/'
key = 'candTree'
key_failed = 'candTree_failed'
verbose = False

if (opt.YEAR == '2016'): years = [2016]
elif (opt.YEAR == '2017'): years = [2017]
elif (opt.YEAR == '2018'): years = [2018]
elif (opt.YEAR == 'Full'): years = [2016,2017,2018]

obs_name = opt.OBSNAME
doubleDiff=False
if(obs_name == 'rapidity4l'):
    obs_reco = 'ZZy'
    obs_gen = 'GENrapidity4l'
elif(obs_name == 'pT4l'):
    obs_reco = 'ZZPt'
    obs_gen = 'GENpT4l'
elif(obs_name == 'massZ1'):
    obs_reco = 'Z1Mass'
    obs_gen = 'GENmassZ1'
elif(obs_name == 'massZ2'):
    obs_reco = 'Z2Mass'
    obs_gen = 'GENmassZ2'
elif(obs_name == 'mass4l'):
    obs_reco = 'ZZMass'
    obs_gen = 'GENmass4l'
elif(obs_name == "njets_pt30_eta2p5"):
    obs_reco = "njets_pt30_eta2p5"
    obs_gen = "GENnjets_pt30_eta2p5"
elif(obs_name == 'pTj1'):
    obs_reco = 'pTj1'
    obs_gen = 'GENpTj1'
elif(obs_name == 'mass4l'):
    obs_reco = 'ZZMass'
    obs_gen = 'GENmass4l'
elif(obs_name == 'costhetastar'):
    obs_reco = 'costhetastar'
    obs_gen = 'GENcosThetaStar'
elif(obs_name == 'costhetaZ1'):
    obs_reco = 'helcosthetaZ1'
    obs_gen  = 'GENcosTheta1'
elif(obs_name == 'costhetaZ2'):
    obs_reco = 'helcosthetaZ2'
    obs_gen  = 'GENcosTheta2'
elif(obs_name == 'phi'):
    obs_reco = 'helphi'
    obs_gen  = 'GENPhi'
elif(obs_name == 'phistar'):
    obs_reco = 'phistarZ1'
    obs_gen  = 'GENPhi1'
elif(obs_name == 'massZ1 vs massZ2'):
    doubleDiff = True
    obs_name = 'massZ1_massZ2'
    obs_reco = 'Z1Mass'
    obs_reco_2nd = 'Z2Mass'
    obs_gen = 'GENmassZ1'
    obs_gen_2nd = 'GENmassZ2'
elif(obs_name == 'njets_pt30_eta2p5 vs pT4l'):
    doubleDiff = True
    obs_name = 'njets_pt30_eta2p5_pT4l'
    obs_reco = 'njets_pt30_eta2p5'
    obs_reco_2nd = 'ZZPt'
    obs_gen = 'GENnjets_pt30_eta2p5'
    obs_gen_2nd = 'GENpT4l'
elif(obs_name == 'njets_pt30_eta2p5 vs pTHj'):
    doubleDiff = True
    obs_name = 'njets_pt30_eta2p5_pTHj'
    obs_reco = 'njets_pt30_eta2p5'
    obs_reco_2nd = 'pTHj'
    obs_gen = 'GENnjets_pt30_eta2p5'
    obs_gen_2nd = 'GENpTHj'
else:
    print 'Missing varibale'
    sys.exit()


sys.path.append('../inputs')
_temp = __import__('inputs_sig_'+obs_name+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
obs_bins = _temp.observableBins
sys.path.remove('../inputs')

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
        d_sup[signal] = pd.concat([d_sig[year][signal], d_sig_failed[year][signal]], ignore_index=True)
    d_sig_tot[year] = d_sup


# Create dataframe FullRun2
if(opt.YEAR == 'Full'):
    d_sig_full = {}
    for signal in signals:
        frame = [d_sig_tot[year][signal] for year in years]
        d_sig_full[signal] = pd.concat(frame, ignore_index=True)
else: # If I work with one year only, the FullRun2 df coincides with d_sig_tot (it is useful when fullNNLOPS is calculated)
    d_sig_full = d_sig_tot[int(opt.YEAR)]


chans = ['4e', '4mu', '2e2mu', '4l']
m4l_low = 105.0
m4l_high = 140.0
acceptance = {}
qcdUncert = {}
pdfUncert = {}
lenObsBins = len(obs_bins)
if doubleDiff: lenObsBins = lenObsBins+1
for chan in chans:
    for genbin in range(lenObsBins-1):
        getPdfUncert(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin, obs_name, 'std', year)
        getPdfUncert(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin, obs_name, 'NNLOPS', year)
if (obs_reco.startswith("njets")) and not doubleDiff:
    for chan in chans:
        for genbin in range(lenObsBins-2): # last bin is >=3
            for signal in signals:
                processBin = signal+'_'+chan+'_'+obs_reco+'_genbin'+str(genbin)
                processBinPlus1 = signal+'_'+chan+'_'+obs_reco+'_genbin'+str(genbin+1)
                acceptance[processBin] = acceptance[processBin]-acceptance[processBinPlus1]
                qcdUncert[processBin]['uncerUp'] = sqrt(qcdUncert[processBin]['uncerUp']*qcdUncert[processBin]['uncerUp']+qcdUncert[processBinPlus1]['uncerUp']*qcdUncert[processBinPlus1]['uncerUp'])
                qcdUncert[processBin]['uncerDn'] = sqrt(qcdUncert[processBin]['uncerDn']*qcdUncert[processBin]['uncerDn']+qcdUncert[processBinPlus1]['uncerDn']*qcdUncert[processBinPlus1]['uncerDn'])

                processBin = signal+'_NNLOPS_'+chan+'_'+obs_name+'_genbin'+str(genbin)
                processBinPlus1 = signal+'_NNLOPS_'+chan+'_'+obs_name+'_genbin'+str(genbin+1)
                acceptance[processBin] = acceptance[processBin]-acceptance[processBinPlus1]
                qcdUncert[processBin]['uncerUp'] = sqrt(qcdUncert[processBin]['uncerUp']*qcdUncert[processBin]['uncerUp']+qcdUncert[processBinPlus1]['uncerUp']*qcdUncert[processBinPlus1]['uncerUp'])
                qcdUncert[processBin]['uncerDn'] = sqrt(qcdUncert[processBin]['uncerDn']*qcdUncert[processBin]['uncerDn']+qcdUncert[processBinPlus1]['uncerDn']*qcdUncert[processBinPlus1]['uncerDn'])
with open('../inputs/accUnc_'+obs_name+'.py', 'w') as f:
    f.write('acc = '+str(acceptance)+';\n')
    f.write('qcdUncert = '+str(qcdUncert)+' \n')
    f.write('pdfUncert = '+str(pdfUncert)+' \n')
