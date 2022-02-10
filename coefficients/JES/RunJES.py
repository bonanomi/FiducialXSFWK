import numpy
import os
import pandas as pd
import matplotlib as mpl
# if os.environ.get('DISPLAY','') == '':
#     print('no display found. Using non-interactive Agg backend')
#     mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import optparse, os, sys
import json
import numpy
import ROOT
from binning import binning
from createdf_jes import skim_df
from tabulate import tabulate
from zx import zx

print 'Welcome in RunJES!'

jesNames = ['Total', 'Abs', 'Abs_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavQCD', 'HF', 'HF_year', 'RelBal', 'RelSample_year']

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
    parser.add_option('',   '--AC', action='store_true', dest='AC', default=False, help='AC samples')
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


def computeRatio(nominal, upVar, dnVar):
    if nominal == 0:
        return '-'
    else:
        dn_ratio = round(dnVar/nominal,3)
        up_ratio = round(upVar/nominal,3)
        if up_ratio==dn_ratio and up_ratio==1.000:
             return '-'
        elif up_ratio == dn_ratio:
            return str(dn_ratio)
        else:
            return str(dn_ratio)+'/'+str(up_ratio)

# ------------------------------- FUNCTIONS TO CALCULATE JES BOTH IN 1D AND 2D CASE -------------------------------
def getJes(channel, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, obs_name, year, obs_reco_2nd = 'None', obs_gen_2nd = 'None', obs_name_2nd = 'None'):
    if not doubleDiff:
        #RecoBin limits I'm considering
        obs_reco_low = obs_bins[recobin]
        obs_reco_high = obs_bins[recobin+1]
    elif doubleDiff:
        #First variable
        obs_reco_low = obs_bins[recobin][0]
        obs_reco_high = obs_bins[recobin][1]
        #Second variable
        obs_reco_2nd_low = obs_bins[recobin][2]
        obs_reco_2nd_high = obs_bins[recobin][3]

    datafr = d_sig[year]
    datafr_qqzz = d_bkg[year]['qqzz']
    datafr_ggzz = d_bkg[year]['ggzz']
    datafr_zx = d_ZX[year]

    # ------------- ZX computation (at this stage it is done only the inclusive JES) -------------
    if doubleDiff:
        processBin = '_'+channel+'_'+str(year)+'_'+obs_reco+'_'+obs_reco_2nd+'_recobin'+str(recobin)
    else:
        processBin = '_'+channel+'_'+str(year)+'_'+obs_reco+'_recobin'+str(recobin)

    cutobs_reco_zx = (datafr_zx[obs_reco] >= obs_reco_low) & (datafr_zx[obs_reco] < obs_reco_high)
    cutobs_reco_jesup_zx = (datafr_zx[obs_reco+'_jesup'] >= obs_reco_low) & (datafr_zx[obs_reco+'_jesup'] < obs_reco_high)
    cutobs_reco_jesdn_zx = (datafr_zx[obs_reco+'_jesdn'] >= obs_reco_low) & (datafr_zx[obs_reco+'_jesdn'] < obs_reco_high)
    cutchan_reco_zx = (datafr_zx['FinState_reco'] == channel)
    cutm4l_reco_zx = (datafr_zx['ZZMass'] > m4l_low) & (datafr_zx['ZZMass'] < m4l_high) & (datafr_zx['FinState_reco'] == channel)
    if doubleDiff:
        cutobs_reco_zx &= (datafr_zx[obs_reco_2nd] >= obs_reco_2nd_low) & (datafr_zx[obs_reco_2nd] < obs_reco_2nd_high)
        cutobs_reco_jesup_zx &= (datafr_zx[obs_reco_2nd+'_jesup'] >= obs_reco_2nd_low) & (datafr_zx[obs_reco_2nd+'_jesup'] < obs_reco_2nd_high)
        cutobs_reco_jesdn_zx &= (datafr_zx[obs_reco_2nd+'_jesdn'] >= obs_reco_2nd_low) & (datafr_zx[obs_reco_2nd+'_jesdn'] < obs_reco_2nd_high)

    evts['ZX'+processBin] = datafr_zx[cutm4l_reco_zx & cutobs_reco_zx & cutchan_reco_zx]['yield_SR'].sum()
    evts['ZX_jesup'+processBin] = datafr_zx[cutm4l_reco_zx & cutobs_reco_jesup_zx & cutchan_reco_zx]['yield_SR'].sum()
    evts['ZX_jesdn'+processBin] = datafr_zx[cutm4l_reco_zx & cutobs_reco_jesdn_zx & cutchan_reco_zx]['yield_SR'].sum()
    evts_noWeight['ZX'+processBin] = len(datafr_zx[cutm4l_reco_zx & cutobs_reco_zx & cutchan_reco_zx])
    evts_noWeight['ZX_jesup'+processBin] = len(datafr_zx[cutm4l_reco_zx & cutobs_reco_jesup_zx & cutchan_reco_zx])
    evts_noWeight['ZX_jesdn'+processBin] = len(datafr_zx[cutm4l_reco_zx & cutobs_reco_jesdn_zx & cutchan_reco_zx])

    ratio['ZX'+processBin] = computeRatio(evts['ZX'+processBin], evts['ZX_jesup'+processBin], evts['ZX_jesdn'+processBin])

    # ------------- signal, qqZZ, and ggZZ computation -------------
    for i in jesNames:
        if doubleDiff:
            processBin = '_'+i+'_'+channel+'_'+str(year)+'_'+obs_reco+'_'+obs_reco_2nd+'_recobin'+str(recobin)
        else:
            processBin = '_'+i+'_'+channel+'_'+str(year)+'_'+obs_reco+'_recobin'+str(recobin)

        cutobs_reco = (datafr[obs_reco] >= obs_reco_low) & (datafr[obs_reco] < obs_reco_high)
        cutobs_reco_qqzz = (datafr_qqzz[obs_reco] >= obs_reco_low) & (datafr_qqzz[obs_reco] < obs_reco_high)
        cutobs_reco_ggzz = (datafr_ggzz[obs_reco] >= obs_reco_low) & (datafr_ggzz[obs_reco] < obs_reco_high)
        cutobs_reco_jesup = (datafr[obs_reco+'_jesup_'+i] >= obs_reco_low) & (datafr[obs_reco+'_jesup_'+i] < obs_reco_high)
        cutobs_reco_jesdn = (datafr[obs_reco+'_jesdn_'+i] >= obs_reco_low) & (datafr[obs_reco+'_jesdn_'+i] < obs_reco_high)
        cutobs_reco_jesup_qqzz = (datafr_qqzz[obs_reco+'_jesup_'+i] >= obs_reco_low) & (datafr_qqzz[obs_reco+'_jesup_'+i] < obs_reco_high)
        cutobs_reco_jesdn_qqzz = (datafr_qqzz[obs_reco+'_jesdn_'+i] >= obs_reco_low) & (datafr_qqzz[obs_reco+'_jesdn_'+i] < obs_reco_high)
        cutobs_reco_jesup_ggzz = (datafr_ggzz[obs_reco+'_jesup_'+i] >= obs_reco_low) & (datafr_ggzz[obs_reco+'_jesup_'+i] < obs_reco_high)
        cutobs_reco_jesdn_ggzz = (datafr_ggzz[obs_reco+'_jesdn_'+i] >= obs_reco_low) & (datafr_ggzz[obs_reco+'_jesdn_'+i] < obs_reco_high)

        if doubleDiff:
            cutobs_reco &= (datafr[obs_reco_2nd] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd] < obs_reco_2nd_high)
            cutobs_reco_jesup &= (datafr[obs_reco_2nd+'_jesup_'+i] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd+'_jesup_'+i] < obs_reco_2nd_high)
            cutobs_reco_jesdn &= (datafr[obs_reco_2nd+'_jesdn_'+i] >= obs_reco_2nd_low) & (datafr[obs_reco_2nd+'_jesdn_'+i] < obs_reco_2nd_high)
            cutobs_reco_jesup_qqzz &= (datafr_qqzz[obs_reco_2nd+'_jesup_'+i] >= obs_reco_2nd_low) & (datafr_qqzz[obs_reco_2nd+'_jesup_'+i] < obs_reco_2nd_high)
            cutobs_reco_jesdn_qqzz &= (datafr_qqzz[obs_reco_2nd+'_jesdn_'+i] >= obs_reco_2nd_low) & (datafr_qqzz[obs_reco_2nd+'_jesdn_'+i] < obs_reco_2nd_high)
            cutobs_reco_jesup_ggzz &= (datafr_ggzz[obs_reco_2nd+'_jesup_'+i] >= obs_reco_2nd_low) & (datafr_ggzz[obs_reco_2nd+'_jesup_'+i] < obs_reco_2nd_high)
            cutobs_reco_jesdn_ggzz &= (datafr_ggzz[obs_reco_2nd+'_jesdn_'+i] >= obs_reco_2nd_low) & (datafr_ggzz[obs_reco_2nd+'_jesdn_'+i] < obs_reco_2nd_high)

        cutm4l_reco = (datafr['ZZMass'] > m4l_low) & (datafr['ZZMass'] < m4l_high) & (datafr['FinState_reco'] == channel)
        cutm4l_reco_qqzz = (datafr_qqzz['ZZMass'] > m4l_low) & (datafr_qqzz['ZZMass'] < m4l_high) & (datafr_qqzz['FinState_reco'] == channel)
        cutm4l_reco_ggzz = (datafr_ggzz['ZZMass'] > m4l_low) & (datafr_ggzz['ZZMass'] < m4l_high) & (datafr_ggzz['FinState_reco'] == channel)

        # --------------- signal events ---------------
        evts['signal'+processBin] = datafr[cutm4l_reco & cutobs_reco]['weight_histo_reco'].sum()
        evts['signal_jesup'+processBin] = datafr[cutm4l_reco & cutobs_reco_jesup]['weight_histo_reco'].sum()
        evts['signal_jesdn'+processBin] = datafr[cutm4l_reco & cutobs_reco_jesdn]['weight_histo_reco'].sum()
        evts_noWeight['signal'+processBin] = len(datafr[cutm4l_reco & cutobs_reco]['weight_histo_reco'])
        evts_noWeight['signal_jesup'+processBin] = len(datafr[cutm4l_reco & cutobs_reco_jesup])
        evts_noWeight['signal_jesdn'+processBin] = len(datafr[cutm4l_reco & cutobs_reco_jesdn])

        ratio['signal'+processBin] = computeRatio(evts['signal'+processBin], evts['signal_jesup'+processBin], evts['signal_jesdn'+processBin])

        # --------------- qqzz events ---------------
        evts['qqzz'+processBin] = datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_qqzz]['weight_histo_reco'].sum()
        evts['qqzz_jesup'+processBin] = datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_jesup_qqzz]['weight_histo_reco'].sum()
        evts['qqzz_jesdn'+processBin] = datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_jesdn_qqzz]['weight_histo_reco'].sum()
        evts_noWeight['qqzz'+processBin] = len(datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_qqzz]['weight_histo_reco'])
        evts_noWeight['qqzz_jesup'+processBin] = len(datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_jesup_qqzz])
        evts_noWeight['qqzz_jesdn'+processBin] = len(datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_jesdn_qqzz])

        ratio['qqzz'+processBin] = computeRatio(evts['qqzz'+processBin], evts['qqzz_jesup'+processBin], evts['qqzz_jesdn'+processBin])

        # --------------- ggzz events ---------------
        evts['ggzz'+processBin] = datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_ggzz]['weight_histo_reco'].sum()
        evts['ggzz_jesup'+processBin] = datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_jesup_ggzz]['weight_histo_reco'].sum()
        evts['ggzz_jesdn'+processBin] = datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_jesdn_ggzz]['weight_histo_reco'].sum()
        evts_noWeight['ggzz'+processBin] = len(datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_ggzz])
        evts_noWeight['ggzz_jesup'+processBin] = len(datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_jesup_ggzz])
        evts_noWeight['ggzz_jesdn'+processBin] = len(datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_jesdn_ggzz])

        ratio['ggzz'+processBin] = computeRatio(evts['ggzz'+processBin], evts['ggzz_jesup'+processBin], evts['ggzz_jesdn'+processBin])


def doGetJes(obs_reco, obs_gen, obs_name, obs_bins, obs_reco_2nd = 'None', obs_gen_2nd = 'None', obs_name_2nd = 'None',):
    chans = ['4e', '4mu', '2e2mu']
    m4l_low = 105
    m4l_high = 140

    nBins = len(obs_bins)
    if not doubleDiff: nBins = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)
    for year in years:
        for chan in chans:
            for recobin in range(nBins):
                getJes(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, obs_name, year, obs_reco_2nd, obs_gen_2nd, obs_name_2nd)


        with open('JESNP_'+obsname_out+'_'+str(year)+'.py', 'w') as f:
                f.write('obsbins = ' + str(obs_bins) + '\n')
                f.write('JESNP = ' + str(ratio) + '\n')

        with open('JESNP_evts_'+obsname_out+'_'+str(year)+'.py', 'w') as f:
                f.write('obsbins = ' + str(obs_bins) + '\n')
                f.write('evts = ' + str(evts) + '\n')


## ---------------------------- Main ----------------------------
sys.path.append('../../inputs/')
from observables import observables
_temp = __import__('observables', globals(), locals(), ['observables'])
observables = _temp.observables
sys.path.remove('../../inputs/')

obsname = opt.OBSNAME
obsname_out = obsname
doubleDiff = False

if ' vs ' in obsname:
    doubleDiff = True
    obsname_1st = opt.OBSNAME.split(' vs ')[0]
    obsname_2nd = opt.OBSNAME.split(' vs ')[1]
    obsname_out = obsname_1st + '_' + obsname_2nd

obs_reco = observables[obsname]['obs_reco']
obs_gen = observables[obsname]['obs_gen']
if doubleDiff:
    obs_reco_2nd = observables[obsname]['obs_reco_2nd']
    obs_gen_2nd = observables[obsname]['obs_gen_2nd']

# year    = opt.YEAR
if (opt.YEAR == '2016'): years = [2016]
if (opt.YEAR == '2017'): years = [2017]
if (opt.YEAR == '2018'): years = [2018]
if (opt.YEAR == 'Full'): years = [2016,2017,2018]


m4l_low = opt.LOWER_BOUND
m4l_high = opt.UPPER_BOUND

print (obsname, years, m4l_low, m4l_high, obsname_out)

obs_bins, doubleDiff = binning(opt.OBSNAME)

# Generate dataframes
d_sig = {}
d_bkg = {}
for year in years:
    if doubleDiff: sig, bkg = skim_df(year, doubleDiff, obs_reco, obs_reco_2nd)
    else: sig, bkg = skim_df(year, doubleDiff, obs_reco)
    d_sig[year] = sig
    d_sig[year] = pd.concat([d_sig[year]['ggH125'], d_sig[year]['VBFH125'], d_sig[year]['WH125'], d_sig[year]['ZH125'], d_sig[year]['ttH125']])
    d_bkg[year] = bkg

d_ZX = {}
d_ZX = zx()


ratio = {} # Dict with ratio of jesup and jesdown variations wrt nominal value
evts = {} # Dict with number of weighted events for each process
evts_noWeight = {} # Dict with number of events for each process
if not doubleDiff: doGetJes(obs_reco, obs_gen, obsname, obs_bins)
else: doGetJes(obs_reco, obs_gen, obsname, obs_bins, obs_reco_2nd, obsname_2nd)

with open('JESNP_'+obsname_out+'.py', 'w') as f:
        f.write('obsbins = ' + str(obs_bins) + '\n')
        f.write('JESNP = ' + str(ratio) + '\n')

with open('JESNP_evts_'+obsname_out+'.py', 'w') as f:
        f.write('obsbins = ' + str(obs_bins) + '\n')
        f.write('evts = ' + str(evts) + '\n')
        f.write('evts_noWeight = ' + str(evts_noWeight) + '\n')