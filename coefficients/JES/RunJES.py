import numpy
import os
import pandas as pd
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import optparse, os, sys
import json
import numpy
import ROOT
from binning import binning
from createdf_jes import skim_df

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

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()


# ------------------------------- FUNCTIONS TO CALCULATE JES BOTH IN 1D AND 2D CASE -------------------------------
# -------------- 1D --------------
def computeJES(obsname, obs_bins, year, fState, m4l_low, m4l_high):

    if fState != '4l':
        sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high)) & (d_sig[year].FinState_reco == fState)
    else:
        sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high))

    for i in jesNames:
        fig = plt.figure(figsize=(6,5))
        # Upper frame
        frame1 = fig.add_axes((.1, .35, .8, .6))
        plt.title('CMS', weight = 'bold', loc = 'left', fontsize = 15)
        plt.title(fState, loc = 'right', fontsize = 14)
        n, bins, _= plt.hist([d_sig[year][sel][obsname],d_sig[year][sel][obsname+'_jesup_'+i],d_sig[year][sel][obsname+'_jesdn_'+i]],
                             bins = obs_bins,
                             weights = [d_sig[year][sel]['weight_reco'],d_sig[year][sel]['weight_reco'],d_sig[year][sel]['weight_reco']],
                             histtype='step',
                             color = ['blue', 'orange', 'green'],
                             label = ['nominal', 'jesup_'+i, 'jesdn_'+i])
        if obsname == 'pTj1': plt.xlim(0,250)
        ratio_up = n[1]/n[0]
        ratio_dn = n[2]/n[0]
        # Fill JES dictionary
        for b in range(len(bins)-1):
            if abs(ratio_dn[b] - ratio_up[b]) < 0.001: # If they are equal, only one value in datacards
                jesNP[year,fState,'recobin_'+str(b),i] = '%.3f' %(ratio_dn[b])
            else:
                jesNP[year,fState,'recobin_'+str(b),i] = '%.3f/%.3f' %(ratio_dn[b], ratio_up[b])
        plt.legend()
        frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
        # Lower (ratio) frame
        frame2=fig.add_axes((.1,.1,.8,.2))
        for index in range(len(bins)-1):
            plt.hlines(ratio_up[index], bins[index], bins[index+1], 'orange')
            plt.hlines(ratio_dn[index], bins[index], bins[index+1], 'green')
        plt.axhline(1,color='blue',ls='--')
        if obsname == 'pTj1': plt.xlim(0,250)
        plt.savefig('plots/'+obs_reco+'/jesNP_'+str(year)+'_'+fState+'_'+i+'.png',dpi=400, bbox_inches='tight')
        # plt.show()
        plt.show()
        plt.close(fig)

# -------------- 2D --------------
def computeJES_2D(obsname, obsname_2nd, obs_bins, year, fState, m4l_low, m4l_high):

    if fState != '4l':
        sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high)) & (d_sig[year].FinState_reco == fState)
    else:
        sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high))

    for i in jesNames:
        fig = plt.figure(figsize=(6,5))
        # Upper frame
        frame1 = fig.add_axes((.1, .35, .8, .6))
        plt.title('CMS', weight = 'bold', loc = 'left', fontsize = 15)
        plt.title(fState, loc = 'right', fontsize = 14)
        for index in obs_bins:
            sel_nominal = sel & (d_sig[year][obsname] >= obs_bins[index][0]) & (d_sig[year][obsname] < obs_bins[index][1]) & (d_sig[year][obsname_2nd] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd] < obs_bins[index][3])
            plt.hlines(d_sig[year][sel_nominal]['weight_reco'].sum(), index, index+1, 'blue', label='nominal')

            sel_up = sel & (d_sig[year][obsname+'_jesup_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesup_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesup_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesup_'+i] < obs_bins[index][3])
            plt.hlines(d_sig[year][sel_up]['weight_reco'].sum(), index, index+1, 'orange', label='jesup_'+i)

            sel_dn = sel & (d_sig[year][obsname+'_jesdn_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesdn_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] < obs_bins[index][3])
            plt.hlines(d_sig[year][sel_dn]['weight_reco'].sum(), index, index+1, 'green', label='jesdn_'+i)

        # We have to set the legend entries manually
        legend_elements = [Line2D([0], [0], color='blue', label='nominal'),
                           Line2D([0], [0], color='orange', label='jesup_'+i),
                           Line2D([0], [0], color='green', label='jesdown_'+i)]
        plt.legend(handles=legend_elements)
        plt.gca().get_xaxis().set_major_locator(ticker.MultipleLocator(1))
        plt.xlim(0,len(obs_bins))
        frame1.set_xticklabels([])
        # Lower (ratio) frame
        frame2=fig.add_axes((.1,.1,.8,.2))
        plt.axhline(1,color='blue',ls='--')
        ratio_up = []
        ratio_dn = []
        for index in obs_bins:
            sel_nominal = sel & (d_sig[year][obsname] >= obs_bins[index][0]) & (d_sig[year][obsname] < obs_bins[index][1]) & (d_sig[year][obsname_2nd] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd] < obs_bins[index][3])

            sel_up = sel & (d_sig[year][obsname+'_jesup_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesup_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesup_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesup_'+i] < obs_bins[index][3])
            plt.hlines(d_sig[year][sel_up]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum(), index, index+1, 'orange', label='jesup_'+i)
            ratio_up.append(d_sig[year][sel_up]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum())

            sel_dn = sel & (d_sig[year][obsname+'_jesdn_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesdn_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] < obs_bins[index][3])
            ratio_dn.append(d_sig[year][sel_dn]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum())
            plt.hlines(d_sig[year][sel_dn]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum(), index, index+1, 'green', label='jesdn_'+i)

        plt.xlim(0,len(obs_bins))
        xticks = [n+0.5 for n in range(len(obs_bins))]
        xticks_label = ['Bin'+str(n+1) for n in range(len(obs_bins))]
        plt.xticks(xticks, xticks_label)
        plt.savefig('plots/'+obsname+'_'+obsname_2nd+'/jesNP_'+str(year)+'_'+fState+'_'+i+'.png',dpi=400, bbox_inches='tight')
        plt.show()

        for b in range(len(obs_bins)):
            if abs(ratio_dn[b] - ratio_up[b]) < 0.001: # If they are equal, only one value in datacards
                jesNP[year,fState,'recobin_'+str(b),i] = '%.3f' %(ratio_dn[b])
            else:
                jesNP[year,fState,'recobin_'+str(b),i] = '%.3f/%.3f' %(ratio_dn[b], ratio_up[b])


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

if doubleDiff:
    obs_reco_2nd = observables[obsname]['obs_reco_2nd']

obs_reco = observables[obsname]['obs_reco']


# year    = opt.YEAR
if (opt.YEAR == '2016'): years = [2016]
if (opt.YEAR == '2017'): years = [2017]
if (opt.YEAR == '2018'): years = [2018]
if (opt.YEAR == 'Full'): years = [2016,2017,2018]


m4l_low = opt.LOWER_BOUND
m4l_high = opt.UPPER_BOUND

print (obsname, years, m4l_low, m4l_high, obsname_out)

obs_bins, doubleDiff = binning(opt.OBSBINS)


# Generate dataframes
d_sig = {}
for year in years:
    if doubleDiff: sig = skim_df(year, doubleDiff, obs_reco, obs_reco_2nd)
    else: sig = skim_df(year, doubleDiff, obs_reco)
    d_sig[year] = sig
    d_sig[year] = pd.concat([d_sig[year]['ggH125'], d_sig[year]['VBFH125'], d_sig[year]['WH125'], d_sig[year]['ZH125'], d_sig[year]['ttH125']])


# Check if the folder for plots exist
if not os.path.exists('plots/'+obsname_out):
    os.makedirs('plots/'+obsname_out)

jesNP = {}
if not doubleDiff:
    for year in years:
        for fState in ['2e2mu', '4mu' ,'4e','4l']:
            # year = str(year)
            computeJES(obs_reco, obs_bins, year, fState, m4l_low, m4l_high)
else:
    for year in years:
        for fState in ['2e2mu', '4mu' ,'4e','4l']:
            computeJES_2D(obs_reco, obs_reco_2nd, obs_bins, year, fState, m4l_low, m4l_high)

with open('../../inputs/JESNP_'+obsname_out+'.py', 'w') as f:
        f.write('obsbins = ' + str(obs_bins) + '\n')
        f.write('JESNP = ' + str(jesNP) + '\n')
