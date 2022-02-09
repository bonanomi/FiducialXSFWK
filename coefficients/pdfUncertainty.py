import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import uproot3 as uproot
from math import sqrt, log
import sys,os
import optparse
import itertools
import math
import ROOT
import json

from paths import path

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
    parser.add_option('',   '--nnlops', action='store_true', dest='NNLOPS', default=False, help='Calculate uncert for ggH_NNLOPS sample')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='int',default=105.0,   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='int',default=140.0,   help='Upper bound for m4l')
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
    # coeff = (lumi * 1000 * xsec) / gen
    #Gen
    weight_gen = df.genHEPMCweight # * df.PUWeight //Not sure if PUWeight have to be included
    #Specific for pdfUncertainty
    df['qcdWeight_0'] = weight_gen * df.LHEweight_QCDscale_muR1_muF1
    df['qcdWeight_1'] = weight_gen * df.LHEweight_QCDscale_muR1_muF2
    df['qcdWeight_2'] = weight_gen * df.LHEweight_QCDscale_muR1_muF0p5
    df['qcdWeight_3'] = weight_gen * df.LHEweight_QCDscale_muR2_muF1
    df['qcdWeight_4'] = weight_gen * df.LHEweight_QCDscale_muR2_muF2
    df['qcdWeight_5'] = weight_gen * df.LHEweight_QCDscale_muR2_muF0p5
    df['qcdWeight_6'] = weight_gen * df.LHEweight_QCDscale_muR0p5_muF1
    df['qcdWeight_7'] = weight_gen * df.LHEweight_QCDscale_muR0p5_muF2
    df['qcdWeight_8'] = weight_gen * df.LHEweight_QCDscale_muR0p5_muF0p5
    df['PDFupWeight'] = (weight_gen * df.LHEweight_QCDscale_muR1_muF1) * df.LHEweight_PDFVariation_Up / abs(weight_gen * df.LHEweight_QCDscale_muR1_muF1) # qcdWeights[0]*pdfENVup/abs(qcdWeights[0])
    df['PDFdnWeight'] = (weight_gen * df.LHEweight_QCDscale_muR1_muF1) * df.LHEweight_PDFVariation_Dn / abs(weight_gen * df.LHEweight_QCDscale_muR1_muF1) # qcdWeights[0]*pdfENVdn/abs(qcdWeights[0])
    if opt.NNLOPS:
        #In case of NNLOPS I should always use nnloWeights[0], for coding purposes I set all weight_gen_NNLOPS_NNLO_* equal to weight_gen_NNLOPS_NNLO_0
        for i in range(0,27):
            df['nnlopsWeight_'+str(i)] = [x[i] for x in df['nnlopsWeight']]
        ref = df['nnlopsWeight_0'] * weight_gen
        for i in range(0,27):
            df['nnlopsWeight_'+str(i)] = df['nnlopsWeight_'+str(i)] * weight_gen / ref
        # In case of NNLOPS re-write the PDF weights
        df['PDFupWeight'] = df['nnlopsWeight_0'] * df.LHEweight_PDFVariation_Up / abs(weight_gen * df.LHEweight_QCDscale_muR1_muF1) # nnloWeights[0]*pdfENVup/abs(qcdWeights[0])
        df['PDFdnWeight'] = df['nnlopsWeight_0'] * df.LHEweight_PDFVariation_Dn / abs(weight_gen * df.LHEweight_QCDscale_muR1_muF1) # nnloWeights[0]*pdfENVdn/abs(qcdWeights[0])
    return df

# Uproot to generate pandas
def prepareTrees(year):
    d_sig = {}
    d_sig_failed = {}
    for signal in signals_original:
        # if opt.AC_HYP: fname = eos_path_sig + 'AC%i_MELA' %year
        # else: fname = eos_path_sig + '%i_MELA' %year
        fname = eos_path_sig + '%i_MELA' %year
        if not opt.NNLOPS: fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
        else: fname += '/'+signal+'_NNLOPS/'+signal+'_NNLOPS_reducedTree_MC_'+str(year)+'.root'
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
        # if opt.AC_HYP: fname = eos_path_sig + 'AC%i_MELA' %year
        # else: fname = eos_path_sig + '%i_MELA' %year
        fname = eos_path_sig + '%i_MELA' %year
        if not opt.NNLOPS: fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
        else: fname += '/'+signal+'_NNLOPS/'+signal+'_NNLOPS_reducedTree_MC_'+str(year)+'.root'
        input_file = ROOT.TFile(fname)
        hCounters = input_file.Get("Counters")
        gen_sig[signal] = hCounters.GetBinContent(40)
    return gen_sig

def createDataframe(d_sig,fail,gen,xsec,signal,lumi,obs_reco,obs_gen,obs_reco_2nd='None',obs_gen_2nd='None'):
    b_sig = ['EventNumber','GENmass4l', 'GENlep_id', 'GENlep_MomId',
             'GENlep_MomMomId', 'GENlep_Hindex', 'GENZ_DaughtersId',
             'GENZ_MomId', 'passedFiducialSelection_bbf', 'PUWeight', 'genHEPMCweight',
             'LHEweight_QCDscale_muR1_muF1','LHEweight_QCDscale_muR1_muF2','LHEweight_QCDscale_muR1_muF0p5',
             'LHEweight_QCDscale_muR2_muF1','LHEweight_QCDscale_muR2_muF2','LHEweight_QCDscale_muR2_muF0p5',
             'LHEweight_QCDscale_muR0p5_muF1','LHEweight_QCDscale_muR0p5_muF2','LHEweight_QCDscale_muR0p5_muF0p5',
             'LHEweight_PDFVariation_Up', 'LHEweight_PDFVariation_Dn']
    if opt.NNLOPS: b_sig.append('nnlopsWeight')
    if (obs_gen != 'GENmass4l'): b_sig.append(obs_gen)
    if (obs_gen_2nd!='None'): b_sig.append(obs_gen_2nd)
    if signal == 'ggH125' and not opt.NNLOPS: b_sig.append('ggH_NNLOPS_weight') #Additional entry for the weight in case of ggH
    if not fail:
        b_sig.extend(['ZZMass', 'Z1Flav', 'Z2Flav', 'lep_genindex', 'lep_Hindex', 'overallEventWeight',
                      'L1prefiringWeight','dataMCWeight', 'trigEffWeight'])
        if (obs_reco!='ZZMass'): b_sig.append(obs_reco)
        if (obs_reco_2nd!='None'): b_sig.append(obs_reco_2nd)

    df = d_sig.pandas.df(b_sig, flatten = False)
    if fail: #Negative branches for failed events (it is useful when creating fiducial pandas)
        df['ZZMass'] = -1
        df['Z1Flav'] = -1
        df['Z2Flav'] = -1
        df['lep_genindex'] = -1
        df['lep_Hindex'] = -1
        df['overallEventWeight'] = -1
        df['L1prefiringWeight'] = -1
        df['dataMCWeight'] = -1
        df['trigEffWeight'] = -1
        if (obs_reco != 'ZZMass'): df[obs_reco] = -1
        if (obs_reco_2nd!='None'): df[obs_reco_2nd] = -1
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
        # df = df.drop(columns=['ggH_NNLOPS_weight'])

    return df

# Set up data frames
def dataframes(year, doubleDiff):
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
        if doubleDiff:
            d_df_sig[signal] = createDataframe(d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi,obs_reco,obs_gen,obs_reco_2nd,obs_gen_2nd)
        else:
            d_df_sig[signal] = createDataframe(d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi,obs_reco,obs_gen)
        print 'Signal created'
        if doubleDiff:
            d_df_sig_failed[signal] = createDataframe(d_sig_failed[signal],True,gen_sig[signal],xsec_sig[signal],signal,lumi,obs_reco,obs_gen,obs_reco_2nd,obs_gen_2nd)
        else:
            d_df_sig_failed[signal] = createDataframe(d_sig_failed[signal],True,gen_sig[signal],xsec_sig[signal],signal,lumi,obs_reco,obs_gen)
        print 'Signal failed created'
    return d_df_sig, d_df_sig_failed


# Merge WplusH125 and WminusH125
def skim_df(year, doubleDiff):
    d_df_sig, d_df_sig_failed = dataframes(year, doubleDiff)
    d_skim_sig = {}
    d_skim_sig_failed = {}
    frames = []
    for signal in signals_original:
        if (signal == 'WplusH125') or (signal == 'WminusH125'):
            frames.append(d_df_sig[signal])
        else:
            d_skim_sig[signal] = d_df_sig[signal]

    frames = []
    for signal in signals_original:
        if (signal == 'WplusH125') or (signal == 'WminusH125'):
            frames.append(d_df_sig_failed[signal])
        else:
            d_skim_sig_failed[signal] = d_df_sig_failed[signal]

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
            cutobs_gen = datafr[obs_gen] >= obs_gen_low
        else:
            cutobs_gen = (datafr[obs_gen] >= obs_gen_low) & (datafr[obs_gen] < obs_gen_high)
            if (obs_name=='Dcp'): cutobs_gen = (datafr[obs_gen] >= obs_gen_low) & (datafr[obs_gen] < obs_gen_high)
            if doubleDiff:
                cutobs_gen &= (datafr[obs_gen_2nd] >= obs_gen_2nd_low) & (datafr[obs_gen_2nd] < obs_gen_2nd_high)
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
        if not opt.NNLOPS: coefficients['fs'] = datafr[cutchan_gen_out]['qcdWeight_0'].sum()
        else: coefficients['fs'] = datafr[cutchan_gen_out]['nnlopsWeight_0'].sum()
        if type=='std':
            for i in range(0,9):
                maxWeight = 9
                coefficients['fs'+str(i)] = datafr[cutchan_gen_out]['qcdWeight_'+str(i)].sum()
                # coefficients['fid'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen'+nnlops+'_NNLO_'+str(i)].sum()
                coefficients['fid'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['qcdWeight_'+str(i)].sum()
                coefficients['fid'+str(i)] = coefficients['fid'+str(i)] * (1/coefficients['fs']) #This stands for "Scale" in the original code
                # coefficients['fidraw'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen'+nnlops+'_NNLO_'+str(i)].sum()
                # coefficients['fidraw'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['weight_gen_NNLO_'+str(i)].sum()
                # coefficients['fidraw'+str(i)] = coefficients['fidraw'+str(i)] * (1/coefficients['fs'+str(i)]) #This stands for "Scale" in the original code
                coefficients['fs'+str(i)] = coefficients['fs'+str(i)] * (1/coefficients['fs'+str(i)]) #This stands for "Scale" in the original code
        else:
            for i in range(0,27):
                maxWeight = 27
                if (i==5 or i==7 or i==11 or i==14 or i==15 or i==16 or i==17 or i==19 or i==21 or i==22 or i==23 or i==25): continue
                coefficients['fs'+str(i)] = datafr[cutchan_gen_out]['nnlopsWeight_'+str(i)].sum()
                coefficients['fid'+str(i)] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['nnlopsWeight_'+str(i)].sum()
                coefficients['fid'+str(i)] = coefficients['fid'+str(i)] * (1/coefficients['fs']) #This stands for "Scale" in the original code
                coefficients['fs'+str(i)] = coefficients['fs'+str(i)] * (1/coefficients['fs'+str(i)]) #This stands for "Scale" in the original code

        coefficients['fidPDF_up'] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['PDFupWeight'].sum()
        coefficients['fidPDF_up'] = coefficients['fidPDF_up'] * (1/coefficients['fs']) #This stands for "Scale" in the original code
        coefficients['fidPDF_dn'] = datafr[passedFiducialSelection & cutm4l_gen & cutobs_gen & cutchan_gen & cuth4l_gen]['PDFdnWeight'].sum()
        coefficients['fidPDF_dn'] = coefficients['fidPDF_dn'] * (1/coefficients['fs']) #This stands for "Scale" in the original code
        fsintegral = coefficients['fs']
        coefficients['fs'] = coefficients['fs'] * (1/coefficients['fs'])

        if coefficients['fs']>0:
            if verbose: print coefficients['fs'],coefficients['fid0']
            # acceptance[processBin] = coefficients['fid0']/coefficients['fs']
            qcderrup=1.0; qcderrdn=1.0;
            # accerrup=1.0; accerrdn=1.0;
            if verbose: print processBin,coefficients['fid0']
            for i in range(0,maxWeight):
                if type=='NNLOPS' and (i==5 or i==7 or i==11 or i==14 or i==15 or i==16 or i==17 or i==19 or i==21 or i==22 or i==23 or i==25): continue
                ratio = coefficients['fid'+str(i)]/coefficients['fid0']
                # print(channel, genbin, i, ratio, coefficients['fid'+str(i)], coefficients['fid0'])
                if verbose: print i, 'ratio', ratio
                if ratio>qcderrup: qcderrup = coefficients['fid'+str(i)]/coefficients['fid0']
                if ratio<qcderrdn: qcderrdn = coefficients['fid'+str(i)]/coefficients['fid0']
            # print('')
                # acci = coefficients['fidraw'+str(i)]/coefficients['fs'+str(i)]
                # if verbose: print i, 'acc', acci
                # if verbose: print coefficients['fidraw'+str(i)], coefficients['fs'+str(i)]
                # if acci/acceptance[processBin]>accerrup: accerrup=acci/acceptance[processBin]
                # if acci/acceptance[processBin]<accerrdn: accerrdn=acci/acceptance[processBin]
            qcdUncert[processBin]={'uncerDn':abs(qcderrdn-1.0),'uncerUp':abs(qcderrup-1.0)}
            pdferr_up = coefficients['fidPDF_up']/coefficients['fid0']
            pdferr_dn = coefficients['fidPDF_dn']/coefficients['fid0']
            pdfUncert[processBin] = {"uncerDn":abs(pdferr_dn-1.0),"uncerUp":abs(pdferr_up-1.0)}

            if verbose: print(processBin,acceptance[processBin],0.0,qcderrup,qcderrdn,pdferr_up,pdferr_dn)

# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------
# if opt.AC_HYP: signals_original = signals = ['ggH'+opt.AC_HYP+'_M125']
# else: signals_original = signals = ['ggH125']
signals_original = signals = ['ggH125']
eos_path_sig = path['eos_path_sig']
key = 'candTree'
key_failed = 'candTree_failed'
verbose = False

if (opt.YEAR == '2016'): years = [2016]
elif (opt.YEAR == '2017'): years = [2017]
elif (opt.YEAR == '2018'): years = [2018]
elif (opt.YEAR == 'Full'): years = [2016,2017,2018]

obs_name = opt.OBSNAME
doubleDiff=False
if ' vs ' in obs_name: doubleDiff=True

if doubleDiff:
    obs_name = opt.OBSNAME.split(' vs ')[0]
    obs_name_2nd = opt.OBSNAME.split(' vs ')[1]
    obs_name_2d = opt.OBSNAME
    obs_name = obs_name + '_' + obs_name_2nd
else:
    obs_name = opt.OBSNAME

# sys.path.append('../inputs/')
_temp = __import__('observables', globals(), locals(), ['observables'], -1)
observables = _temp.observables
# sys.path.remove('../inputs/')

if doubleDiff:
    obs_reco = observables[obs_name_2d]['obs_reco']
    obs_reco_2nd = observables[obs_name_2d]['obs_reco_2nd']
    obs_gen = observables[obs_name_2d]['obs_gen']
    obs_gen_2nd = observables[obs_name_2d]['obs_gen_2nd']
else:
    obs_reco = observables[obs_name]['obs_reco']
    obs_gen = observables[obs_name]['obs_gen']

print 'Following observables extracted from dictionary: RECO = ',obs_reco,' GEN = ',obs_gen
if doubleDiff:
    print 'It is a double-differential measurement: RECO_2nd = ',obs_reco_2nd,' GEN_2nd = ',obs_gen_2nd

sys.path.append('../inputs')
_temp = __import__('inputs_sig_'+obs_name+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
obs_bins = _temp.observableBins
sys.path.remove('../inputs')

# Generate dataframes
d_sig = {}
d_sig_failed = {}
for year in years:
    sig, sig_failed = skim_df(year, doubleDiff)
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
m4l_low = opt.LOWER_BOUND
m4l_high = opt.UPPER_BOUND
acceptance = {}
qcdUncert = {}
pdfUncert = {}
lenObsBins = len(obs_bins)
if doubleDiff: lenObsBins = lenObsBins+1
for chan in chans:
    for genbin in range(lenObsBins-1):
        if not opt.NNLOPS: getPdfUncert(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin, obs_name, 'std', year)
        else: getPdfUncert(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, genbin, obs_name, 'NNLOPS', year)
if (obs_reco.startswith("njets")) and not doubleDiff:
    for chan in chans:
        for genbin in range(lenObsBins-2): # last bin is >=3
            for signal in signals:
                if not opt.NNLOPS:
                    processBin = signal+'_'+chan+'_'+obs_reco+'_genbin'+str(genbin)
                    processBinPlus1 = signal+'_'+chan+'_'+obs_reco+'_genbin'+str(genbin+1)
                else:
                    processBin = signal+'_NNLOPS_'+chan+'_'+obs_name+'_genbin'+str(genbin)
                    processBinPlus1 = signal+'_NNLOPS_'+chan+'_'+obs_name+'_genbin'+str(genbin+1)
                acceptance[processBin] = acceptance[processBin]-acceptance[processBinPlus1]
                qcdUncert[processBin]['uncerUp'] = sqrt(qcdUncert[processBin]['uncerUp']*qcdUncert[processBin]['uncerUp']+qcdUncert[processBinPlus1]['uncerUp']*qcdUncert[processBinPlus1]['uncerUp'])
                qcdUncert[processBin]['uncerDn'] = sqrt(qcdUncert[processBin]['uncerDn']*qcdUncert[processBin]['uncerDn']+qcdUncert[processBinPlus1]['uncerDn']*qcdUncert[processBinPlus1]['uncerDn'])
                # acceptance[processBin] = acceptance[processBin]-acceptance[processBinPlus1]
                # qcdUncert[processBin]['uncerUp'] = sqrt(qcdUncert[processBin]['uncerUp']*qcdUncert[processBin]['uncerUp']+qcdUncert[processBinPlus1]['uncerUp']*qcdUncert[processBinPlus1]['uncerUp'])
                # qcdUncert[processBin]['uncerDn'] = sqrt(qcdUncert[processBin]['uncerDn']*qcdUncert[processBin]['uncerDn']+qcdUncert[processBinPlus1]['uncerDn']*qcdUncert[processBinPlus1]['uncerDn'])
if opt.NNLOPS: obs_name = obs_name+'_NNLOPS'
# elif opt.AC_HYP: obs_name = obs_name+'_'+opt.AC_HYP
with open('../inputs/accUnc_'+obs_name+'.py', 'w') as f:
    f.write('acc = '+str(acceptance)+';\n')
    f.write('qcdUncert = '+str(qcdUncert)+' \n')
    f.write('pdfUncert = '+str(pdfUncert)+' \n')
