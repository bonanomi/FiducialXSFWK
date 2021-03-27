## Small tweak to use Helvetica font
import matplotlib.font_manager as font_manager
import matplotlib as mpl
# font_dirs = ['/eos/home-m/mbonanom/SWAN_projects/LHE/fonts/', ]
# font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
# font_list = font_manager.createFontList(font_files)
# font_manager.fontManager.ttflist.extend(font_list)

import uproot
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from os import listdir
from os.path import isfile, join
import optparse, sys

url = 'https://gist.githubusercontent.com/bonanomi/d14780f7562cb2a22fdd753a9d4459d4/raw/c77dcb028433c73b26f3ad84ef61f54abf13236e/MyMPLStyle'
plt.style.use(url)

sys.path.append('./inputs')
sys.path.append('./LHScans')
from plotUtils import *

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-d', '--dir',    dest='SOURCEDIR',  type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--unfoldModel',dest='UNFOLD',type='string',default='SM_125', help='Name of the unfolding model for central value')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "inclusive", "pT", "eta", "Njets"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='Use results from fixed fraction fit, default is False')
    parser.add_option('',   '--logScale', action='store_true', dest='LOGSCALE', default=False, help='Use log scale for differential plot, default is False')
    parser.add_option('',   '--setLog', action='store_true', dest='SETLOG', default=False, help='set plot to log scale y, default is False')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--lumiscale', type='string', dest='LUMISCALE', default='1.0', help='Scale yields')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

obs_bins   = parseBins(opt.OBSBINS)
obs_name   = opt.OBSNAME
obs_name, doubleDiff   = getObsName(obs_name)
print(obs_name)
if doubleDiff: 
    obs_label  = opt.OBSNAME.split(' vs ')[1]
else: 
    obs_label = opt.OBSNAME
logScale   = opt.LOGSCALE
theoryMass = opt.THEORYMASS
model      = opt.UNFOLD
if (opt.FIXFRAC): floatfix = '_fixfrac'
else: floatfix = ''

if (opt.UNBLIND): _temp = __import__('resultsXS_LHScan_observed_'+obs_name+'_v3'+floatfix, globals(), locals(), ['resultsXS'], 0)
else: _temp = __import__('resultsXS_LHScan_expected_'+obs_name+'_v3'+floatfix, globals(), locals(), ['resultsXS'], 0)
resultsXS = _temp.resultsXS

def getPoints(resultsXS, obs_name, obs_bins):
    toPlot = []
    toPlot_sys = []
    bin_centers = []
    bin_widths = []
    if doubleDiff:
        for _bin in obs_bins:
            xsec = resultsXS[model+'_'+obs_name+'_genbin%i' %_bin]['central']
            up = resultsXS[model+'_'+obs_name+'_genbin%i' %_bin]['uncerUp']
            dn = resultsXS[model+'_'+obs_name+'_genbin%i' %_bin]['uncerDn']
            upstat = resultsXS[model+'_'+obs_name+'_genbin%i_statOnly' %_bin]['uncerUp']
            dnstat = resultsXS[model+'_'+obs_name+'_genbin%i_statOnly' %_bin]['uncerDn']
            upsys = np.sqrt(up**2 - upstat**2)
            dnsys = np.sqrt(dn**2 - dnstat**2)
            width = (obs_bins[_bin][3] - obs_bins[_bin][2]) # bin center	        
            #if width > 100:
            #    width = 60
            if width == 800: # first bin in njets measurements
                width = 1
            if (_bin == 0) and ('njets' in obs_name): width = 1
            center = obs_bins[_bin][2]+0.5*width
            bin_centers.append(center)
            bin_widths.append(width)
            toPlot.append([xsec/width, up/width, dn/width, width])
            toPlot_sys.append([xsec/width, upsys/width, dnsys/width, width])
    else:
        for idx, _bin in enumerate(range(len(obs_bins)-1)):
            xsec = resultsXS[model+'_'+obs_name+'_genbin%i' %_bin]['central']
            up = resultsXS[model+'_'+obs_name+'_genbin%i' %_bin]['uncerUp']
            dn = resultsXS[model+'_'+obs_name+'_genbin%i' %_bin]['uncerDn']
            upstat = resultsXS[model+'_'+obs_name+'_genbin%i_statOnly' %_bin]['uncerUp']
            dnstat = resultsXS[model+'_'+obs_name+'_genbin%i_statOnly' %_bin]['uncerDn']
            upsys = np.sqrt(up**2 - upstat**2)
            dnsys = np.sqrt(dn**2 - dnstat**2)
            width = (obs_bins[idx+1] - obs_bins[idx]) # bin center
            
            if width > 100:
                width = 60
            center = obs_bins[idx]+0.5*width
            bin_centers.append(center)
            bin_widths.append(width)
            toPlot.append([xsec/width, up/width, dn/width, width])
            toPlot_sys.append([xsec/width, upsys/width, dnsys/width, width])

    return toPlot, toPlot_sys, bin_centers, bin_widths

def makePlot(theory, data):
    ggH_powheg, ggH_minloHJ, XH = theory 
    toPlot, toPlot_sys, bin_centers, bin_widths = data

    i = 0
    xpos = []
    ratio_data = []
    _widths = []
    fig = plt.figure(figsize=(6,5))
    frame1=fig.add_axes((.1,.3,.8,.6))
    plt.title('CMS', weight = 'bold', loc = 'left', fontsize = 15)
    for p, psys, th, xh, w, width in zip(toPlot, toPlot_sys, ggH_minloHJ, XH, bin_centers, bin_widths):
        if doubleDiff: w = i
        plt.errorbar(w, p[0], yerr=np.array([[abs(p[2]), p[1]]]).T, marker = 'o', markersize = 4., capsize = 0, color = 'k')
        plt.errorbar(w, psys[0], yerr=np.array([[abs(psys[2]), psys[1]]]).T, marker = 'None', linewidth = 3, capsize = 0, color = 'r')
        
        # Error bars on th pred
    #     plt.errorbar(w-0.25*w, psys[0], yerr=np.array([[abs(psys[2]), psys[1]]]).T, marker = 'None', linewidth = 4, color = 'b', alpha = 0.6)
    #     plt.errorbar(w+0.25*w, psys[0], yerr=np.array([[abs(psys[2]), psys[1]]]).T, marker = 'None', linewidth = 4, color = 'darkorange', alpha = 0.6)

        xpos.append(w)
        _widths.append(width)
        ratio_data.append([(p[0]/((th+xh)/width)), abs(p[2])/((th+xh)/width), p[1]/((th+xh)/width), abs(psys[2])/((th+xh)/width), psys[1]/((th+xh)/width)])

        i+=1    

    _widths = np.array(_widths)
    ## Make Differential also TH pred

    if not doubleDiff:
    	ggH_powheg = (np.array(XH)+np.array(ggH_powheg))/np.array(_widths)
    	ggH_minloHJ = (np.array(XH)+np.array(ggH_minloHJ))/np.array(_widths)
    	XH = np.array(XH)/np.array(_widths)
    else:
    	ggH_powheg = (np.array(XH)+np.array(ggH_powheg))/np.array(_widths)
    	ggH_minloHJ = (np.array(XH)+np.array(ggH_minloHJ))/np.array(_widths)
    	XH = np.array(XH)/np.array(_widths)

    plt.step(xpos, ggH_powheg, where = 'mid', label = 'POWHEG')
    ## just to close up the line
    plt.hlines(ggH_powheg[0], np.array(xpos[0]) - 0.5, xpos[0], color = 'tab:blue')
    plt.hlines(ggH_powheg[-1:], xpos[-1:], np.array(xpos[-1:])+0.5, color = 'tab:blue')

    ratio_powheg = np.array(ggH_powheg)/np.array(ggH_minloHJ)
    ratio_nnlops = np.array(ggH_minloHJ)/np.array(ggH_minloHJ)

    plt.step(xpos, ggH_minloHJ, where = 'mid', color = 'darkorange', label = 'NNLOPS')
    ## just to close up the line
    plt.hlines(ggH_minloHJ[0], np.array(xpos[0]) - 0.5, xpos[0], color = 'darkorange')
    plt.hlines(ggH_minloHJ[-1:], xpos[-1:], np.array(xpos[-1:]) + 0.5, color = 'darkorange')

    # Add XH histo
    if doubleDiff: 
        xh_width = 0.8
    else:
        xh_width = np.array(bin_centers) - np.array(bin_centers)*0.5
	    
    plt.fill_between(np.array(xpos), XH, step = 'mid', hatch = 'xxx', facecolor = 'none', linewidth = 0.0, label = 'XH', edgecolor = 'darkgreen', alpha = 0.5)
    plt.fill_between([np.array(xpos)[0]-0.5,np.array(xpos)[0]], XH[0], step = 'mid', hatch = 'xxx', facecolor = 'none', linewidth = 0.0, edgecolor = 'darkgreen', alpha = 0.5)
    plt.fill_between([np.array(xpos)[-1],np.array(xpos)[-1]+0.5], XH[-1], step = 'mid', hatch = 'xxx', facecolor = 'none', linewidth = 0.0, edgecolor = 'darkgreen', alpha = 0.5)
    
    #xticks = np.arange(0, len(xpos))
    xticks = np.arange(-0.5, len(xpos)-0.5, 0.5)
    xticks = [xt for xt in xticks if xt%1 != 0]
    if doubleDiff:
        labels = []
        for thebin in obs_bins:
            thelabel = str([obs_bins[thebin][2], obs_bins[thebin][3]])
            thelabel = thelabel[1:-1]
            if (thebin == 0) and ('njets' in obs_name) and ('pT4l' not in obs_name): thelabel = r'$N_{jets}$ = 0'
            if '1000.0' in thelabel: thelabel.replace('1000.0', 'inf')
            labels.append(thelabel)
    if not doubleDiff: plt.xticks(obs_bins[:-1])

    plt.legend(frameon=False, fontsize = 14, loc = 1)
    plt.tick_params(axis='x', which = 'minor', length = 0)

    if logScale: plt.yscale('log')

    plt.ylabel(r'd$\sigma_{\mathrm{fid}}$/d%s (fb)' %obs_name, ha='right', y=1.0)
    if doubleDiff: plt.xlim(np.array(xpos)[0]-0.5, np.array(xpos)[-1]+0.5)

    frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
    frame2=fig.add_axes((.1,.08,.8,.2)) 

    plt.step(xpos, ratio_powheg, where = 'mid', label = 'POWHEG')
    # ## just to close up the line
    plt.hlines(ratio_powheg[0], np.array(xpos[0]) - 0.5, xpos[0], color = 'tab:blue')
    plt.hlines(ratio_powheg[-1], xpos[-1], np.array(xpos[-1]) + 0.5, color = 'tab:blue')

    plt.step(xpos, ratio_nnlops, where = 'post', label = 'NNLOPS')
    ## just to close up the line
    plt.hlines(ratio_nnlops[0], np.array(xpos[0]) - 0.5, xpos[0], color = 'darkorange')
    plt.hlines(ratio_nnlops[-1], xpos[-1], np.array(xpos[-1]) + 0.5, color = 'darkorange')

    for idx, ratio in zip(xpos, ratio_data):
        plt.errorbar(idx, ratio[0], yerr=np.array([[ratio[1], ratio[2]]]).T, marker = 'o', markersize = 4., capsize = 0, color = 'k')
        plt.errorbar(idx, ratio[0], yerr=np.array([[ratio[3], ratio[4]]]).T, marker = 'None', linewidth = 2, capsize = 0, color = 'r')
    # plt.ylim(0, 2)

    # Add to reset axis aspect
    plt.fill_between(np.array(xpos), XH, step = 'mid', hatch = 'xxx', facecolor = 'none', alpha = 0)

    #labels = ['Bin %i' %(i+1) for i in xticks]
    if not doubleDiff:
        plt.xticks()
    #     plt.xticks(obs_bins[:-1])
    else:
        plt.xticks(np.array(xticks), labels, fontsize = 10, rotation = 30)
    plt.tick_params(axis='x', which = 'minor', length = 0)
    plt.ylabel('Ratio to NNLOPS')
    plt.ylim(0, 2)
    if doubleDiff: plt.xlim(np.array(xpos)[0]-0.5, np.array(xpos)[-1]+0.5)
    plt.xlabel(obs_label+' [GeV]', ha='right', x=1.0)
    plt.savefig('plot_%s' %obs_name, bbox_inches = 'tight')

def plotXS(obs_name, obs_bins):

    _temp = __import__('inputs_sig_'+obs_name+'_2017', globals(), locals(), ['acc'], 0)
    acc = _temp.acc

    if(opt.YEAR=='Full'):
        _temp = __import__('inputs_sig_'+obs_name+'_NNLOPS_Full', globals(), locals(), ['acc'], 0)
        acc_NNLOPS = _temp.acc
    else:
        _temp = __import__('inputs_sig_'+obs_name+'_NNLOPS_'+opt.YEAR, globals(), locals(), ['acc'], 0)
        acc_NNLOPS = _temp.acc

    ggH_powheg, ggH_minloHJ, XH = computeXSTH(acc, acc_NNLOPS, theoryMass, obs_name, obs_bins, doubleDiff)
    theory = [ggH_powheg, ggH_minloHJ, XH]

    ## TODO unc_nnlops, unc_powheg, unc_XH = computeUNCTH(acc, acc_NNLOPS, theoryMass, obs_name, obs_bins)

    toPlot, toPlot_sys, bin_centers, bin_widths  = getPoints(resultsXS, obs_name, obs_bins)
    data = [toPlot, toPlot_sys, bin_centers, bin_widths]
    makePlot(theory, data)

plotXS(obs_name, obs_bins)
