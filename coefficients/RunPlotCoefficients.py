import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from math import sqrt, log
import sys
import itertools
import math
import json

def tickBin(binning): # Create labels for each bin
    tick = []
    for i in range(len(binning)-1):
        tick.append(str(binning[i])+'-'+str(binning[i+1]))
    return tick

def matrix(obs_bins, obs_reco, obs_gen, obs_name, label):
    chans = ['4e', '4mu', '2e2mu']
    sys.path.append('../inputs/')

    if len(obs_bins)-1 == 6:
        obs_bins_label = [0,1,2,3,4,5,6] # Number of row and columns
        obs_bins_label_medium = [0.5,1.5,2.5,3.5,4.5,5.5] # Medium point of each row and column
        bin_max = 6
    elif len(obs_bins)-1 == 7:
        obs_bins_label = [0,1,2,3,4,5,6,7] # Number of row and columns
        obs_bins_label_medium = [0.5,1.5,2.5,3.5,4.5,5.5,6.5] # Medium point of each row and column
        bin_max = 7

    tickLabel = tickBin(obs_bins)
    print(tickLabel)

    for year in years:
        _temp = __import__('inputs_sig_'+obs_name+'_'+str(year), globals(), locals(), ['eff', 'err_eff']) # Open file with coefficients
        eff = _temp.eff
        err_eff = _temp.err_eff
        for signal in signals:
            for fState in ['4e', '4mu', '2e2mu']:
                eps = np.zeros((len(obs_bins_label)-1, len(obs_bins_label)-1)) # Efficienty matrix
                eps_err = np.zeros((len(obs_bins_label)-1, len(obs_bins_label)-1)) # Efficienty erro matrix
                for genbin in range(len(obs_bins_label)-1):
                    for recobin in range(len(obs_bins_label)-1):
                        processBin = signal+'_'+fState+'_'+obs_name+'_genbin'+str(genbin)+'_recobin'+str(recobin)
                        eps[recobin, genbin] = eff.get(processBin)
                        eps_err[recobin, genbin] = err_eff.get(processBin)
                print(eps)

                # The following two lines set white color for bin with zero efficiency
                my_cmap = matplotlib.cm.get_cmap('rainbow')
                my_cmap.set_under('w')

                plt.pcolormesh(obs_bins_label, obs_bins_label, eps, cmap = my_cmap, alpha=0.6)
                for (i, j), z in np.ndenumerate(eps):
                    if(z >= 0.01):
                        if (eps_err[i,j] < 0.002):
                            plt.text(j+0.5, i+0.5, '{:.4f}'.format(z), ha='center', va='baseline', size = 'medium')
                            plt.text(j+0.5, i+0.5,'\n\n$\pm${:.4f}'.format(eps_err[i,j]), ha='center', va='center_baseline', size = 'x-small')
                        else:
                            plt.text(j+0.5, i+0.5, '{:.3f}'.format(z), ha='center', va='baseline', size = 'medium')
                            plt.text(j+0.5, i+0.5,'\n\n$\pm${:.3f}'.format(eps_err[i,j]), ha='center', va='center_baseline', size = 'x-small')
                plt.clim(vmin=0.01, vmax=1) # Range colorbar
                plt.colorbar(label = 'EFFICIENCY (> 0.01)')
                plt.xticks(obs_bins_label_medium, tickLabel)
                for tick in plt.gca().get_xaxis().get_major_ticks(): # Remove central tick
                    tick.tick1line.set_markersize(0)
                plt.setp(plt.gca().get_xticklabels(), ha="right",rotation = 25, rotation_mode="anchor") # Rotate xticks
                plt.yticks(obs_bins_label_medium, tickLabel)
                for tick in plt.gca().get_yaxis().get_major_ticks(): # Remove central tick
                    tick.tick1line.set_markersize(0)
                plt.xlim(0,bin_max)
                plt.ylim(0,bin_max)
                if obs_name != 'rapidity4l':
                    plt.xlabel(label+' (fid) [GeV]', fontsize=14)
                    plt.ylabel(label+' (reco) [GeV]', fontsize=14)
                else:
                    plt.xlabel(label+' (fid)', fontsize=14)
                    plt.ylabel(label+' (reco)', fontsize=14)
                plt.rcParams['xtick.labelsize']=13
                plt.rcParams['ytick.labelsize']=13

                plt.gca().get_xaxis().set_minor_locator(ticker.MultipleLocator(1))
                plt.gca().get_yaxis().set_minor_locator(ticker.MultipleLocator(1))
                for tick in plt.gca().get_xaxis().get_minor_ticks():
                    tick.tick1line.set_markersize(5)
                for tick in plt.gca().get_yaxis().get_minor_ticks():
                    tick.tick1line.set_markersize(5)
                #Final steps
                plt.title('%s - %s - %s' %(year, signal, fState), loc = 'left', fontweight = 'bold')
                plt.savefig('matrix_eff/%s/eff_%s_%s_%s_%s.png' %(year, year, obs_name, signal, fState), bbox_inches='tight')
                plt.savefig('matrix_eff/%s/eff_%s_%s_%s_%s.pdf' %(year, year, obs_name, signal, fState), bbox_inches='tight')
                plt.tight_layout()
                # plt.clf()
                # plt.cla()
                plt.close()

    sys.path.remove('../inputs/')

signals_original = ['VBFH125', 'ggH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']
signals = ['ggH125', 'VBFH125', 'WH125', 'ZH125', 'ttH125']
eos_path_sig = '/eos/user/a/atarabin/MC_samples/'
key = 'candTree'
key_failed = 'candTree_failed'
# years = [2016, 2017, 2018]
years = [2017]

obs_bins = [0, 0.15, 0.3, 0.6, 0.9, 1.2, 2.5]
obs_reco = 'ZZy'
obs_gen = 'GenHRapidity'
obs_name = 'rapidity4l'
label = '|y$_H$|'

matrix(obs_bins, obs_reco, obs_gen, obs_name, label)
