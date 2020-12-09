import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from math import sqrt, log, trunc
import sys
import optparse
import itertools
import math
import json
import copy

sys.path.append('../inputs/')

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

def tickBin(binning, name): # Create labels for each bin
    tick = []
    if doubleDiff: # Double differential
        if name == 'njets_pt30_eta2p5_pT4l': # Do not print range but just the number of jets
            for i in range(len(binning)):
                if trunc(binning[i][0])==3: # Last bin-> >=3
                    tick.append('$\geq$'+str(trunc(binning[i][0]))+'/'+str(trunc(binning[i][2]))+'-'+str(trunc(binning[i][3])))
                else:
                    tick.append(str(trunc(binning[i][0]))+'/'+str(trunc(binning[i][2]))+'-'+str(trunc(binning[i][3])))
        else:
            for i in range(len(binning)):
                tick.append(str(trunc(binning[i][0]))+'-'+str(trunc(binning[i][1]))+'/'+str(trunc(binning[i][2]))+'-'+str(trunc(binning[i][3])))
    else: # Single differential
        for i in range(len(binning)-1):
            tick.append(str(binning[i])+'-'+str(binning[i+1]))
    return tick

def matrix(obs_bins, obs_name, label):
    chans = ['4e', '4mu', '2e2mu']

    if doubleDiff:
        obs_bins_label = [i for i in range(len(obs_bins))]
        obs_bins_label_medium = [i+0.5 for i in range(len(obs_bins))]
        bin_max = len(obs_bins)
    else:
        obs_bins_label = [i for i in range(len(obs_bins))]
        obs_bins_label_medium = [i+0.5 for i in range(len(obs_bins)-1)]
        bin_max = len(obs_bins)-1

    tickLabel = tickBin(obs_bins, obs_name)
    print(tickLabel)

    for year in years:
        _temp = __import__('inputs_sig_'+obs_name+'_'+str(year), globals(), locals(), ['eff', 'err_eff']) # Open file with coefficients
        eff = _temp.eff
        err_eff = _temp.err_eff
        for signal in signals:
            for fState in ['4e', '4mu', '2e2mu']:
                if not doubleDiff:
                    eps = np.zeros((len(obs_bins_label)-1, len(obs_bins_label)-1)) # Efficienty matrix
                    eps_err = np.zeros((len(obs_bins_label)-1, len(obs_bins_label)-1)) # Efficienty error matrix
                    for genbin in range(len(obs_bins_label)-1):
                        for recobin in range(len(obs_bins_label)-1):
                            processBin = signal+'_'+fState+'_'+obs_name+'_genbin'+str(genbin)+'_recobin'+str(recobin)
                            eps[recobin, genbin] = eff.get(processBin)
                            eps_err[recobin, genbin] = err_eff.get(processBin)
                elif doubleDiff:
                    eps = np.zeros((len(obs_bins_label), len(obs_bins_label))) # Efficienty matrix
                    eps_err = np.zeros((len(obs_bins_label), len(obs_bins_label))) # Efficienty error matrix
                    for genbin in range(len(obs_bins_label)):
                        for recobin in range(len(obs_bins_label)):
                            processBin = signal+'_'+fState+'_'+obs_name+'_genbin'+str(genbin)+'_recobin'+str(recobin)
                            eps[recobin, genbin] = eff.get(processBin)
                            eps_err[recobin, genbin] = err_eff.get(processBin)
                print(eps)

                # The following two lines set white color for bin with zero efficiency
                my_cmap = copy.copy(matplotlib.cm.get_cmap("rainbow"))
                my_cmap.set_under('w')

                # Define pcolormesh
                if doubleDiff: plt.pcolormesh(eps, cmap = my_cmap, alpha=0.6)
                else: plt.pcolormesh(obs_bins_label, obs_bins_label, eps, cmap = my_cmap, alpha=0.6)
                # Fill the matrix
                for (i, j), z in np.ndenumerate(eps):
                    if(z >= 0.01):
                        if not doubleDiff:
                            if (eps_err[i,j] < 0.002):
                                if not doubleDiff:
                                    plt.text(j+0.5, i+0.5, '{:.4f}'.format(z), ha='center', va='baseline', size = 'medium')
                                    plt.text(j+0.5, i+0.5,'\n\n$\pm${:.4f}'.format(eps_err[i,j]), ha='center', va='center_baseline', size = 'x-small')
                            else:
                                plt.text(j+0.5, i+0.5, '{:.3f}'.format(z), ha='center', va='baseline', size = 'medium')
                                plt.text(j+0.5, i+0.5,'\n\n$\pm${:.3f}'.format(eps_err[i,j]), ha='center', va='center_baseline', size = 'x-small')
                        elif doubleDiff:
                            plt.text(j+0.5, i+0.5, '{:.2f}'.format(z), ha='center', va='center', size = 'x-small')
                plt.clim(vmin=0.01, vmax=1) # Range colorbar
                plt.colorbar(label = 'EFFICIENCY (> 0.01)')
                plt.xticks(obs_bins_label_medium, tickLabel)
                for tick in plt.gca().get_xaxis().get_major_ticks(): # Remove central tick
                    tick.tick1line.set_markersize(0)
                if not doubleDiff: plt.setp(plt.gca().get_xticklabels(), ha="right",rotation = 25, rotation_mode="anchor") # Rotate xticks
                elif doubleDiff: plt.setp(plt.gca().get_xticklabels(), ha="right",rotation = 40, rotation_mode="anchor") # Rotate xticks
                plt.yticks(obs_bins_label_medium, tickLabel)
                for tick in plt.gca().get_yaxis().get_major_ticks(): # Remove central tick
                    tick.tick1line.set_markersize(0)
                plt.xlim(0,bin_max)
                plt.ylim(0,bin_max)
                if not doubleDiff: _fontsize = 14
                elif doubleDiff: _fontsize = 13
                if obs_name != 'rapidity4l' and not doubleDiff:
                    plt.xlabel(label+' (fid) [GeV]', fontsize=_fontsize)
                    plt.ylabel(label+' (reco) [GeV]', fontsize=_fontsize)
                else:
                    plt.xlabel(label+' (fid)', fontsize=_fontsize)
                    plt.ylabel(label+' (reco)', fontsize=_fontsize)
                plt.rcParams['xtick.labelsize']=13
                plt.rcParams['ytick.labelsize']=13

                plt.gca().get_xaxis().set_minor_locator(ticker.MultipleLocator(1))
                plt.gca().get_yaxis().set_minor_locator(ticker.MultipleLocator(1))
                for tick in plt.gca().get_xaxis().get_minor_ticks():
                    tick.tick1line.set_markersize(5)
                for tick in plt.gca().get_yaxis().get_minor_ticks():
                    tick.tick1line.set_markersize(5)
                # Grey lines for sub-matrices in case of double-differential measurements
                if type(obs_bins) is dict: # If binning is a dictionary it is a double differential analysis
                    if obs_name == 'njets_pt30_eta2p5_pT4l':
                        plt.axhline(y=4, color='grey', linewidth=0.5, ls='dotted')
                        plt.axhline(y=8, color='grey', linewidth=0.5, ls='dotted')
                        plt.axhline(y=12, color='grey', linewidth=0.5, ls='dotted')
                        plt.axvline(x=4, color='grey', linewidth=0.5, ls='dotted')
                        plt.axvline(x=8, color='grey', linewidth=0.5, ls='dotted')
                        plt.axvline(x=12, color='grey', linewidth=0.5, ls='dotted')

                #Final steps
                plt.title('%s - %s - %s' %(year, signal, fState), loc = 'left', fontweight = 'bold')
                plt.savefig('matrix_eff/%s/eff_%s_%s_%s_%s.png' %(year, year, obs_name, signal, fState), bbox_inches='tight')
                plt.savefig('matrix_eff/%s/eff_%s_%s_%s_%s.pdf' %(year, year, obs_name, signal, fState), bbox_inches='tight')
                plt.tight_layout()
                plt.close()


def nonFid(obs_bins, obs_name, label):
    chans = ['4e', '4mu', '2e2mu']

    if doubleDiff:
        obs_bins_label_x = [i for i in range(len(obs_bins))]
        obs_bins_label_x_medium = [i+0.5 for i in range(len(obs_bins))]
        # bin_max = len(obs_bins)
    else:
        obs_bins_label_x = [i for i in range(len(obs_bins))]
        obs_bins_label_x_medium = [i+0.5 for i in range(len(obs_bins)-1)]
        # bin_max = len(obs_bins)-1
    obs_bins_label_y = [0,1,2,3,4,5] # Number of row and columns
    obs_bins_label_y_medium = [0.5,1.5,2.5,3.5,4.5] # Medium point of each row and column
    bin_max = len(obs_bins)-1


    tickLabel = tickBin(obs_bins, obs_name)
    print(tickLabel)

    for year in years:
        _temp = __import__('inputs_sig_'+obs_name+'_'+str(year), globals(), locals(), ['outinratio', 'err_outinratio']) # Open file with coefficients
        outinratio = _temp.outinratio
        err_outinratio = _temp.err_outinratio
        for fState in ['4e', '4mu', '2e2mu']:
            eps = np.zeros((len(obs_bins_label_y)-1, len(obs_bins_label_x)-1)) # Efficienty matrix
            eps_err = np.zeros((len(obs_bins_label_y)-1, len(obs_bins_label_x)-1)) # Efficienty error matrix
            index_signal = 0
            for signal in signals:
                for recobin in range(len(obs_bins_label_x)-1):
                    processBin = signal+'_'+fState+'_'+obs_name+'_genbin0_recobin'+str(recobin)
                    eps[index_signal, recobin] = outinratio.get(processBin)
                    eps_err[index_signal, recobin] = err_outinratio.get(processBin)
                index_signal +=1
            print(eps)

            # The following two lines set white color for bin with zero efficiency
            my_cmap = copy.copy(matplotlib.cm.get_cmap("rainbow"))
            my_cmap.set_under('w')

            plt.pcolormesh(obs_bins_label_x, obs_bins_label_y, eps, cmap = my_cmap, alpha=0.6)
            for (i, j), z in np.ndenumerate(eps):
                if(z >= 0.00001):
                    if not doubleDiff:
                        if (eps_err[i,j] < 0.001):
                            plt.text(j+0.5, i+0.5, '{:.4f}'.format(z), ha='center', va='baseline', size = 'medium')
                            plt.text(j+0.5, i+0.5,'\n\n$\pm${:.4f}'.format(eps_err[i,j]), ha='center', va='center_baseline', size = 'x-small')
                        else:
                            plt.text(j+0.5, i+0.5, '{:.3f}'.format(z), ha='center', va='baseline', size = 'medium')
                            plt.text(j+0.5, i+0.5,'\n\n$\pm${:.3f}'.format(eps_err[i,j]), ha='center', va='center_baseline', size = 'x-small')
                    elif doubleDiff:
                            plt.text(j+0.5, i+0.5, '{:.3f}'.format(z), ha='center', va='center', size = 'small')
            plt.clim(vmin=0.0001, vmax=0.5) # Range colorbar
            plt.colorbar(label = 'FRACTION (> 0.00001)')
            plt.xticks(obs_bins_label_x_medium, tickLabel)
            for tick in plt.gca().get_xaxis().get_major_ticks(): # Remove central tick
                tick.tick1line.set_markersize(0)
            if not doubleDiff: plt.setp(plt.gca().get_xticklabels(), ha="right",rotation = 25, rotation_mode="anchor") # Rotate xticks
            elif doubleDiff: plt.setp(plt.gca().get_xticklabels(), ha="right",rotation = 40, rotation_mode="anchor") # Rotate xticks
            plt.yticks(obs_bins_label_y_medium, ['ggH', 'VBF', 'WH', 'ZH', 'ttH'])
            for tick in plt.gca().get_yaxis().get_major_ticks(): # Remove central tick
                tick.tick1line.set_markersize(0)
            plt.xlim(0,bin_max)
            plt.ylim(0,5)
            if not doubleDiff: _fontsize = 14
            elif doubleDiff: _fontsize = 13
            if obs_name != 'rapidity4l' and not doubleDiff:
                plt.xlabel(label+' (fid) [GeV]', fontsize=_fontsize)
            else:
                plt.xlabel(label+' (fid)', fontsize=_fontsize)
            plt.rcParams['xtick.labelsize']=13
            plt.rcParams['ytick.labelsize']=13
            # Dimension and resolution of the figure
            if doubleDiff:
                plt.rcParams['figure.figsize']= [10,4.8]
                plt.rcParams['figure.dpi']= 600.

            plt.gca().get_xaxis().set_minor_locator(ticker.MultipleLocator(1))
            plt.gca().get_yaxis().set_minor_locator(ticker.MultipleLocator(1))
            for tick in plt.gca().get_xaxis().get_minor_ticks():
                tick.tick1line.set_markersize(5)
            for tick in plt.gca().get_yaxis().get_minor_ticks():
                tick.tick1line.set_markersize(5)
            #Final steps
            plt.title('%s - %s' %(year, fState), loc = 'left', fontweight = 'bold')
            plt.savefig('matrix_nonfid/%s/nonFid_%s_%s_%s.png' %(year, year, obs_name, fState), bbox_inches='tight')
            plt.savefig('matrix_nonfid/%s/nonFid_%s_%s_%s.pdf' %(year, year, obs_name, fState), bbox_inches='tight')
            plt.tight_layout()
            plt.close()


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

# obs_bins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
# obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
obs_name = opt.OBSNAME
if(obs_name == 'rapidity4l'):
    label = '|y$_H$|'
elif(obs_name == 'pT4l'):
    label = 'p$_T^H$ (GeV)'
elif(obs_name == 'massZ1'):
    label = 'm$_{Z1}$ (GeV)'
elif(obs_name == 'massZ2'):
    label = 'm$_{Z2}$ (GeV)'
elif (obs_name == "njets_pt30_eta2p5"):
    label = 'N$_{jet}$'
elif(obs_name == 'njets_pt30_eta2p5 vs pT4l'):
    obs_name = 'njets_pt30_eta2p5_pT4l' #Change name of obs_name
    label = 'N$_{jet}$/p$_T^H$(GeV)'
elif(obs_name == 'massZ1 vs massZ2'):
    obs_name = 'massZ1_massZ2' #Change name of obs_name
    label = 'm$_{Z1}$(GeV)/m$_{Z2}$(GeV)'

_temp = __import__('inputs_sig_'+obs_name+'_'+opt.YEAR, globals(), locals(), ['observableBins']) # Open file to retrieve the binning
obs_bins = _temp.observableBins

doubleDiff = False
if type(obs_bins) is dict: doubleDiff = True # If binning is a dictionary it is a double differential analysis

matrix(obs_bins, obs_name, label)
nonFid(obs_bins, obs_name, label)
sys.path.remove('../inputs/')
