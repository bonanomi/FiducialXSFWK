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
from tabulate import tabulate

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


# # ------------------------------- FUNCTIONS TO CALCULATE JES BOTH IN 1D AND 2D CASE -------------------------------
# # -------------- 1D --------------
# def computeJES(obsname, obs_bins, year, fState, m4l_low, m4l_high):
#
#     if fState != '4l':
#         sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high)) & (d_sig[year].FinState_reco == fState)
#     else:
#         sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high))
#
#     for i in jesNames:
#         fig = plt.figure(figsize=(6,5))
#         # Upper frame
#         frame1 = fig.add_axes((.1, .35, .8, .6))
#         plt.title('CMS', weight = 'bold', loc = 'left', fontsize = 15)
#         plt.title(fState, loc = 'right', fontsize = 14)
#         n, bins, _= plt.hist([d_sig[year][sel][obsname],d_sig[year][sel][obsname+'_jesup_'+i],d_sig[year][sel][obsname+'_jesdn_'+i]],
#                              bins = obs_bins,
#                              weights = [d_sig[year][sel]['weight_reco'],d_sig[year][sel]['weight_reco'],d_sig[year][sel]['weight_reco']],
#                              histtype='step',
#                              color = ['blue', 'orange', 'green'],
#                              label = ['nominal', 'jesup_'+i, 'jesdn_'+i])
#         if obsname == 'pTj1': plt.xlim(0,250)
#         ratio_up = n[1]/n[0]
#         ratio_dn = n[2]/n[0]
#         # Fill JES dictionary
#         for b in range(len(bins)-1):
#             if abs(ratio_dn[b] - ratio_up[b]) < 0.001: # If they are equal, only one value in datacards
#                 jesNP[year,fState,'recobin_'+str(b),i] = '%.3f' %(ratio_dn[b])
#             else:
#                 jesNP[year,fState,'recobin_'+str(b),i] = '%.3f/%.3f' %(ratio_dn[b], ratio_up[b])
#         plt.legend()
#         frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
#         # Lower (ratio) frame
#         frame2=fig.add_axes((.1,.1,.8,.2))
#         for index in range(len(bins)-1):
#             plt.hlines(ratio_up[index], bins[index], bins[index+1], 'orange')
#             plt.hlines(ratio_dn[index], bins[index], bins[index+1], 'green')
#         plt.axhline(1,color='blue',ls='--')
#         if obsname == 'pTj1': plt.xlim(0,250)
#         plt.savefig('plots/'+obs_reco+'/jesNP_'+str(year)+'_'+fState+'_'+i+'.png',dpi=400, bbox_inches='tight')
#         # plt.show()
#         plt.show()
#         plt.close(fig)
#
# # -------------- 2D --------------
# def computeJES_2D(obsname, obsname_2nd, obs_bins, year, fState, m4l_low, m4l_high):
#
#     if fState != '4l':
#         sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high)) & (d_sig[year].FinState_reco == fState)
#     else:
#         sel = (d_sig[year].ZZMass >= int(m4l_low)) & (d_sig[year].ZZMass <= int(m4l_high))
#
#     for i in jesNames:
#         fig = plt.figure(figsize=(6,5))
#         # Upper frame
#         frame1 = fig.add_axes((.1, .35, .8, .6))
#         plt.title('CMS', weight = 'bold', loc = 'left', fontsize = 15)
#         plt.title(fState, loc = 'right', fontsize = 14)
#         for index in obs_bins:
#             sel_nominal = sel & (d_sig[year][obsname] >= obs_bins[index][0]) & (d_sig[year][obsname] < obs_bins[index][1]) & (d_sig[year][obsname_2nd] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd] < obs_bins[index][3])
#             plt.hlines(d_sig[year][sel_nominal]['weight_reco'].sum(), index, index+1, 'blue', label='nominal')
#
#             sel_up = sel & (d_sig[year][obsname+'_jesup_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesup_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesup_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesup_'+i] < obs_bins[index][3])
#             plt.hlines(d_sig[year][sel_up]['weight_reco'].sum(), index, index+1, 'orange', label='jesup_'+i)
#
#             sel_dn = sel & (d_sig[year][obsname+'_jesdn_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesdn_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] < obs_bins[index][3])
#             plt.hlines(d_sig[year][sel_dn]['weight_reco'].sum(), index, index+1, 'green', label='jesdn_'+i)
#
#         # We have to set the legend entries manually
#         legend_elements = [Line2D([0], [0], color='blue', label='nominal'),
#                            Line2D([0], [0], color='orange', label='jesup_'+i),
#                            Line2D([0], [0], color='green', label='jesdown_'+i)]
#         plt.legend(handles=legend_elements)
#         plt.gca().get_xaxis().set_major_locator(ticker.MultipleLocator(1))
#         plt.xlim(0,len(obs_bins))
#         frame1.set_xticklabels([])
#         # Lower (ratio) frame
#         frame2=fig.add_axes((.1,.1,.8,.2))
#         plt.axhline(1,color='blue',ls='--')
#         ratio_up = []
#         ratio_dn = []
#         for index in obs_bins:
#             sel_nominal = sel & (d_sig[year][obsname] >= obs_bins[index][0]) & (d_sig[year][obsname] < obs_bins[index][1]) & (d_sig[year][obsname_2nd] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd] < obs_bins[index][3])
#
#             sel_up = sel & (d_sig[year][obsname+'_jesup_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesup_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesup_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesup_'+i] < obs_bins[index][3])
#             plt.hlines(d_sig[year][sel_up]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum(), index, index+1, 'orange', label='jesup_'+i)
#             ratio_up.append(d_sig[year][sel_up]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum())
#
#             sel_dn = sel & (d_sig[year][obsname+'_jesdn_'+i] >= obs_bins[index][0]) & (d_sig[year][obsname+'_jesdn_'+i] < obs_bins[index][1]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] >= obs_bins[index][2]) & (d_sig[year][obsname_2nd+'_jesdn_'+i] < obs_bins[index][3])
#             ratio_dn.append(d_sig[year][sel_dn]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum())
#             plt.hlines(d_sig[year][sel_dn]['weight_reco'].sum()/d_sig[year][sel_nominal]['weight_reco'].sum(), index, index+1, 'green', label='jesdn_'+i)
#
#         plt.xlim(0,len(obs_bins))
#         xticks = [n+0.5 for n in range(len(obs_bins))]
#         xticks_label = ['Bin'+str(n+1) for n in range(len(obs_bins))]
#         plt.xticks(xticks, xticks_label)
#         plt.savefig('plots/'+obsname+'_'+obsname_2nd+'/jesNP_'+str(year)+'_'+fState+'_'+i+'.png',dpi=400, bbox_inches='tight')
#         plt.show()
#
#         for b in range(len(obs_bins)):
#             if abs(ratio_dn[b] - ratio_up[b]) < 0.001: # If they are equal, only one value in datacards
#                 jesNP[year,fState,'recobin_'+str(b),i] = '%.3f' %(ratio_dn[b])
#             else:
#                 jesNP[year,fState,'recobin_'+str(b),i] = '%.3f/%.3f' %(ratio_dn[b], ratio_up[b])

# ------------------------------- FUNCTIONS TO CALCULATE JES BOTH IN 1D AND 2D CASE -------------------------------
def getJes(channel, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, year, obs_reco_2nd = 'None', obs_gen_2nd = 'None', obs_name_2nd = 'None'):
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

    datafr = d_sig[year]
    datafr_qqzz = d_bkg[year]['qqzz']
    datafr_ggzz = d_bkg[year]['ggzz']

    for i in jesNames:
        if doubleDiff:
            processBin = '_'+i+'_'+channel+'_'+obs_reco+'_'+obs_reco_2nd+'_genbin'+str(genbin)+'_recobin'+str(recobin)
        else:
            processBin = '_'+i+'_'+channel+'_'+obs_reco+'_genbin'+str(genbin)+'_recobin'+str(recobin)

        cutobs_reco = (abs(datafr[obs_reco]) >= obs_reco_low) & (abs(datafr[obs_reco]) < obs_reco_high)
        cutobs_reco_qqzz = (abs(datafr_qqzz[obs_reco]) >= obs_reco_low) & (abs(datafr_qqzz[obs_reco]) < obs_reco_high)
        cutobs_reco_ggzz = (abs(datafr_ggzz[obs_reco]) >= obs_reco_low) & (abs(datafr_ggzz[obs_reco]) < obs_reco_high)
        cutobs_reco_jesup = (abs(datafr[obs_reco+'_jesup_'+i]) >= obs_reco_low) & (abs(datafr[obs_reco+'_jesup_'+i]) < obs_reco_high)
        cutobs_reco_jesdn = (abs(datafr[obs_reco+'_jesdn_'+i]) >= obs_reco_low) & (abs(datafr[obs_reco+'_jesdn_'+i]) < obs_reco_high)
        cutobs_reco_jesup_qqzz = (abs(datafr_qqzz[obs_reco+'_jesup_'+i]) >= obs_reco_low) & (abs(datafr_qqzz[obs_reco+'_jesup_'+i]) < obs_reco_high)
        cutobs_reco_jesdn_qqzz = (abs(datafr_qqzz[obs_reco+'_jesdn_'+i]) >= obs_reco_low) & (abs(datafr_qqzz[obs_reco+'_jesdn_'+i]) < obs_reco_high)
        cutobs_reco_jesup_ggzz = (abs(datafr_ggzz[obs_reco+'_jesup_'+i]) >= obs_reco_low) & (abs(datafr_ggzz[obs_reco+'_jesup_'+i]) < obs_reco_high)
        cutobs_reco_jesdn_ggzz = (abs(datafr_ggzz[obs_reco+'_jesdn_'+i]) >= obs_reco_low) & (abs(datafr_ggzz[obs_reco+'_jesdn_'+i]) < obs_reco_high)

        cutobs_gen = (abs(datafr[obs_gen]) >= obs_gen_low) & (abs(datafr[obs_gen]) < obs_gen_high)
        if doubleDiff:
            cutobs_reco &= (abs(datafr[obs_reco_2nd]) >= obs_reco_2nd_low) & (abs(datafr[obs_reco_2nd]) < obs_reco_2nd_high)
            cutobs_reco_jesup &= (abs(datafr[obs_reco_2nd+'_jesup_'+i]) >= obs_reco_2nd_low) & (abs(datafr[obs_reco_2nd+'_jesup_'+i]) < obs_reco_2nd_high)
            cutobs_reco_jesdn &= (abs(datafr[obs_reco_2nd+'_jesdn_'+i]) >= obs_reco_2nd_low) & (abs(datafr[obs_reco_2nd+'_jesdn_'+i]) < obs_reco_2nd_high)
            cutobs_gen &= (abs(datafr[obs_gen_2nd]) >= obs_gen_2nd_low) & (abs(datafr[obs_gen_2nd]) < obs_gen_2nd_high)

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
            cutm4l_reco_qqzz = (datafr_qqzz['ZZMass'] > m4l_low) & (datafr_qqzz['ZZMass'] < m4l_high)
            cutm4l_reco_ggzz = (datafr_ggzz['ZZMass'] > m4l_low) & (datafr_ggzz['ZZMass'] < m4l_high)
            cutchan_reco_qqzz = datafr_qqzz['FinState_reco'] == channel
            cutchan_reco_ggzz = datafr_ggzz['FinState_reco'] == channel
            cutchan_gen = datafr['FinState_gen'] == channel
            cutchan_gen_out = datafr['FinState_gen_out'] == channel
        else:
            cutm4l_reco = (datafr['ZZMass'] > m4l_low) & (datafr['ZZMass'] < m4l_high)
            cutm4l_reco_qqzz = (datafr_qqzz['ZZMass'] > m4l_low) & (datafr_qqzz['ZZMass'] < m4l_high)
            cutm4l_reco_ggzz = (datafr_ggzz['ZZMass'] > m4l_low) & (datafr_ggzz['ZZMass'] < m4l_high)
            cutchan_gen = (datafr['FinState_gen'] == '2e2mu') | (datafr['FinState_gen'] == '4e') | (datafr['FinState_gen'] == '4mu')
            cutchan_reco_qqzz = (datafr_qqzz['FinState_reco'] == '2e2mu') | (datafr_qqzz['FinState_reco'] == '4e') | (datafr_qqzz['FinState_reco'] == '4mu')
            cutchan_reco_ggzz = (datafr_ggzz['FinState_reco'] == '2e2mu') | (datafr_ggzz['FinState_reco'] == '4e') | (datafr_ggzz['FinState_reco'] == '4mu')
            cutchan_gen_out = (datafr['FinState_gen_out'] == '2e2mu') | (datafr['FinState_gen_out'] == '4e') | (datafr['FinState_gen_out'] == '4mu')


        # --------------- Fiducial events ---------------
        evts['fiducial'+processBin] = datafr[cutm4l_reco & cutobs_reco & passedFullSelection & cuth4l_reco &
                                       passedFiducialSelection & cuth4l_gen & cutm4l_gen & cutchan_gen & cutobs_gen]['weight_reco'].sum()
        evts['fiducial_jesup'+processBin] = datafr[cutm4l_reco & cutobs_reco_jesup & passedFullSelection & cuth4l_reco &
                                       passedFiducialSelection & cuth4l_gen & cutm4l_gen & cutchan_gen & cutobs_gen]['weight_reco'].sum()
        evts['fiducial_jesdn'+processBin] = datafr[cutm4l_reco & cutobs_reco_jesdn & passedFullSelection & cuth4l_reco &
                                       passedFiducialSelection & cuth4l_gen & cutm4l_gen & cutchan_gen & cutobs_gen]['weight_reco'].sum()

        if evts['fiducial'+processBin] == 0:
            ratio['fiducial'+processBin] = '-'
        else:
            dn_ratio = round(evts['fiducial_jesdn'+processBin]/evts['fiducial'+processBin],3)
            # if dn_ratio == 1.000: dn_ratio = '-'
            up_ratio = round(evts['fiducial_jesup'+processBin]/evts['fiducial'+processBin],3)
            # if up_ratio == 1.000: up_ratio = '-'

            if up_ratio==dn_ratio and up_ratio==1.000: ratio['fiducial'+processBin] = '-'
            elif up_ratio == dn_ratio: ratio['fiducial'+processBin] = str(dn_ratio)
            else : ratio['fiducial'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

        if genbin!=recobin: continue #Non-fiducial events, non-resonant, and bkgs events do not depend on the genbin

        # --------------- Non fiducial events ---------------
        evts['nonFiducial'+processBin] = datafr[cutm4l_reco & cutobs_reco & passedFullSelection &
                                       cuth4l_reco & cutchan_gen_out & (notPassedFiducialSelection | cutnoth4l_gen | cutnotm4l_gen)]['weight_reco'].sum()
        evts['nonFiducial_jesup'+processBin] = datafr[cutm4l_reco & cutobs_reco_jesup & passedFullSelection &
                                       cuth4l_reco & cutchan_gen_out & (notPassedFiducialSelection | cutnoth4l_gen | cutnotm4l_gen)]['weight_reco'].sum()
        evts['nonFiducial_jesdn'+processBin] = datafr[cutm4l_reco & cutobs_reco_jesdn & passedFullSelection &
                                       cuth4l_reco & cutchan_gen_out & (notPassedFiducialSelection | cutnoth4l_gen | cutnotm4l_gen)]['weight_reco'].sum()

        if evts['nonFiducial'+processBin] == 0:
            ratio['nonFiducial'+processBin] = '-'
        else:
            dn_ratio = round(evts['nonFiducial_jesdn'+processBin]/evts['nonFiducial'+processBin],3)
            # if dn_ratio == 1.000: dn_ratio = '-'
            up_ratio = round(evts['nonFiducial_jesup'+processBin]/evts['nonFiducial'+processBin],3)
            # if up_ratio == 1.000: up_ratio = '-'
            ratio['nonFiducial'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

            if up_ratio==dn_ratio and up_ratio==1.000: ratio['nonFiducial'+processBin] = '-'
            elif up_ratio == dn_ratio: ratio['nonFiducial'+processBin] = str(dn_ratio)
            else : ratio['nonFiducial'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

        # --------------- Non resonant events ---------------
        evts['nonResonant'+processBin] = datafr[passedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco]['weight_reco'].sum()
        evts['nonResonant_jesup'+processBin] = datafr[passedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco_jesup]['weight_reco'].sum()
        evts['nonResonant_jesdn'+processBin] = datafr[passedFullSelection & cutm4l_reco & cutnoth4l_reco & cutobs_reco_jesdn]['weight_reco'].sum()

        if evts['nonResonant'+processBin] == 0:
            ratio['nonResonant'+processBin] = '-'
        else:
            dn_ratio = round(evts['nonResonant_jesdn'+processBin]/evts['nonResonant'+processBin],3)
            # if dn_ratio == 1.000: dn_ratio = '-'
            up_ratio = round(evts['nonResonant_jesup'+processBin]/evts['nonResonant'+processBin],3)
            # if up_ratio == 1.000: up_ratio = '-'
            ratio['nonResonant'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

            if up_ratio==dn_ratio and up_ratio==1.000: ratio['nonResonant'+processBin] = '-'
            elif up_ratio == dn_ratio: ratio['nonResonant'+processBin] = str(dn_ratio)
            else : ratio['nonResonant'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

        # --------------- qqzz events ---------------
        evts['qqzz'+processBin] = datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_qqzz & cutchan_reco_qqzz]['weight_reco'].sum()
        evts['qqzz_jesup'+processBin] = datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_jesup_qqzz & cutchan_reco_qqzz]['weight_reco'].sum()
        evts['qqzz_jesdn'+processBin] = datafr_qqzz[cutm4l_reco_qqzz & cutobs_reco_jesdn_qqzz & cutchan_reco_qqzz]['weight_reco'].sum()

        if evts['qqzz'+processBin] == 0:
            ratio['qqzz'+processBin] = '-'
        else:
            dn_ratio = round(evts['qqzz_jesdn'+processBin]/evts['qqzz'+processBin],3)
            # if dn_ratio == 1.000: dn_ratio = '-'
            up_ratio = round(evts['qqzz_jesup'+processBin]/evts['qqzz'+processBin],3)
            # if up_ratio == 1.000: up_ratio = '-'
            ratio['qqzz'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

            if up_ratio==dn_ratio and up_ratio==1.000: ratio['qqzz'+processBin] = '-'
            elif up_ratio == dn_ratio: ratio['qqzz'+processBin] = str(dn_ratio)
            else : ratio['qqzz'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

        # --------------- ggzz events ---------------
        evts['ggzz'+processBin] = datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_ggzz & cutchan_reco_ggzz]['weight_reco'].sum()
        evts['ggzz_jesup'+processBin] = datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_jesup_ggzz & cutchan_reco_ggzz]['weight_reco'].sum()
        evts['ggzz_jesdn'+processBin] = datafr_ggzz[cutm4l_reco_ggzz & cutobs_reco_jesdn_ggzz & cutchan_reco_ggzz]['weight_reco'].sum()

        if evts['ggzz'+processBin] == 0:
            ratio['ggzz'+processBin] = '-'
        else:
            dn_ratio = round(evts['ggzz_jesdn'+processBin]/evts['ggzz'+processBin],3)
            # if dn_ratio == 1.000: dn_ratio = '-'
            up_ratio = round(evts['ggzz_jesup'+processBin]/evts['ggzz'+processBin],3)
            # if up_ratio == 1.000: up_ratio = '-'
            ratio['ggzz'+processBin] = str(dn_ratio)+'/'+str(up_ratio)

            if up_ratio==dn_ratio and up_ratio==1.000: ratio['ggzz'+processBin] = '-'
            elif up_ratio == dn_ratio: ratio['ggzz'+processBin] = str(dn_ratio)
            else : ratio['ggzz'+processBin] = str(dn_ratio)+'/'+str(up_ratio)


def doGetJes(obs_reco, obs_gen, obs_name, obs_bins, obs_reco_2nd = 'None', obs_gen_2nd = 'None', obs_name_2nd = 'None',):
    if obs_reco != 'ZZMass':
        chans = ['4e', '4mu', '2e2mu']
    else:
        chans = ['4l', '4e', '4mu', '2e2mu']
    m4l_low = 105
    m4l_high = 140

    nBins = len(obs_bins)
    if not doubleDiff: nBins = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)
    for year in years:
        for chan in chans:
            for recobin in range(nBins):
                for genbin in range(nBins):
                    getJes(chan, m4l_low, m4l_high, obs_reco, obs_gen, obs_bins, recobin, genbin, obs_name, year, obs_reco_2nd, obs_gen_2nd, obs_name_2nd)

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

obs_bins, doubleDiff = binning(opt.OBSBINS)


# Generate dataframes
d_sig = {}
d_bkg = {}
for year in years:
    if doubleDiff: sig, bkg = skim_df(year, doubleDiff, obs_reco, obs_reco_2nd)
    else: sig, bkg = skim_df(year, doubleDiff, obs_reco)
    d_sig[year] = sig
    d_sig[year] = pd.concat([d_sig[year]['ggH125'], d_sig[year]['VBFH125'], d_sig[year]['WH125'], d_sig[year]['ZH125'], d_sig[year]['ttH125']])
    d_bkg[year] = bkg


ratio = {} # Dict with ratio of jesup and jesdown variations wrt nominal value
evts = {} # Dict with number of events for each process
if not doubleDiff: doGetJes(obs_reco, obs_gen, obsname, obs_bins)
else: doGetJes(obs_reco, obs_gen, obsname, obs_bins, obs_reco_2nd, obsname_2nd)

with open('JESNP_'+obsname_out+'_'+str(year)+'.py', 'w') as f:
        f.write('obsbins = ' + str(obs_bins) + '\n')
        f.write('JESNP = ' + str(ratio) + '\n')

with open('JESNP_evts_'+obsname_out+'_'+str(year)+'.py', 'w') as f:
        f.write('obsbins = ' + str(obs_bins) + '\n')
        f.write('evts = ' + str(evts) + '\n')

'''
# Check if the folder for tables exist
if not os.path.exists('tables/'+obsname_out):
    os.makedirs('tables/'+obsname_out)
# Tables with numerical values
tables = {}
for index, jesName in enumerate(jesNames):
    for fState in ['2e2mu', '4e', '4mu']:
        print(jesName,fState)
        # Indices to print nonFiducial and nonResonant contribution only once per bin
        nonIndex = [(index*3)+k+1 for k in range(0,len(ratio)-1,36*len(obs_bins))] + [(index*3)+k+2 for k in range(0,len(ratio)-1,36*len(obs_bins))]
        nonIndex.sort()

        table = []
        nominal_incl = 0
        up_incl = 0
        dn_incl = 0
        for i in range(0,len(ratio)-1):
            if not (jesName in list(evts.keys())[i*3] and fState in list(evts.keys())[i*3] and jesName in list(evts.keys())[i*3] and not jesName+'_year' in list(evts.keys())[i*3]): continue
            if 'fiducial' in list(evts.keys())[i*3]:
                nominal_incl += list(evts.values())[i*3]
                up_incl += list(evts.values())[i*3+1]
                dn_incl += list(evts.values())[i*3+2]
                table.append([list(evts.keys())[i*3], list(evts.values())[i*3], list(evts.values())[i*3+1], list(evts.values())[i*3+2], list(ratio.values())[i]])
                if 'genbin'+str(len(obs_bins)-1) in list(evts.keys())[i*3]:
                    table.append(['inclusive', nominal_incl, up_incl, dn_incl, str(round(dn_incl/nominal_incl,3))+'/'+str(round(up_incl/nominal_incl,3))])
                    nominal_incl = 0
                    up_incl = 0
                    dn_incl = 0
                    if not 'recobin'+str(len(obs_bins)-1) in list(evts.keys())[i*3]: table.append([])
            if i in nonIndex:
                nominal_incl += list(evts.values())[i*3]
                up_incl += list(evts.values())[i*3+1]
                dn_incl += list(evts.values())[i*3+2]
                table.append([list(evts.keys())[i*3], list(evts.values())[i*3], list(evts.values())[i*3+1], list(evts.values())[i*3+2], list(ratio.values())[i]])
        tables[jesName,fState] = tabulate(table, headers=['bin', 'nominal', 'jes_up', 'jes_dn', 'ratio(dn/up)'], numalign="right", tablefmt="github")
        # print(tables[jesName,fState],'\n')

        with open('tables/JESNP_table_'+obsname_out+'_'+str(year)+'_'+jesName+'_'+fState+'.py', 'w') as f:
            f.write(tables[jesName,fState])

'''
# jesNP = {}
# if not doubleDiff:
#     for year in years:
#         for fState in ['2e2mu', '4mu' ,'4e','4l']:
#             # year = str(year)
#             computeJES(obs_reco, obs_bins, year, fState, m4l_low, m4l_high)
# else:
#     for year in years:
#         for fState in ['2e2mu', '4mu' ,'4e','4l']:
#             computeJES_2D(obs_reco, obs_reco_2nd, obs_bins, year, fState, m4l_low, m4l_high)
#
# with open('../../inputs/JESNP_'+obsname_out+'.py', 'w') as f:
#         f.write('obsbins = ' + str(obs_bins) + '\n')
#         f.write('JESNP = ' + str(jesNP) + '\n')
