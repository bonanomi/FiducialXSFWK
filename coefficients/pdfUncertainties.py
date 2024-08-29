#!/usr/bin/env python
# coding: utf-8

import ROOT

import sys
import uproot 
import awkward as ak
from collections import defaultdict

import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.append('/eos/user/m/mbonanom/run3_trees/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/FiducialXSFWK/inputs')


def get_zh_chan_cuts(channel, tree):

    genz_idx0 = abs(tree['GENZ_DaughtersId'][:,0])
    genz_idx1 = abs(tree['GENZ_DaughtersId'][:,1])
    genz_idx2 = abs(tree['GENZ_DaughtersId'][:,2])
    
    mom_id0 = tree['GENZ_MomId'][:,0]
    mom_id1 = tree['GENZ_MomId'][:,1]
    mom_id2 = tree['GENZ_MomId'][:,2]

    if channel == "4l":
        sel_l0l1 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (genz_idx0 != genz_idx1))
        sel_l0l2 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx0 != genz_idx2))
        sel_l1l2 = ((mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx1 != genz_idx2))
        cutchan_gen_out_4e = ((genz_idx0 == 11) & (genz_idx1 == 11) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 11) & (genz_idx2 == 11) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 11) & (genz_idx2 == 11) & (mom_id1 == 25) & (mom_id2 == 25))
        cutchan_gen_out_4mu = ((genz_idx0 == 13) & (genz_idx1 == 13) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 13) & (genz_idx2 == 13) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 13) & (genz_idx2 == 13) & (mom_id1 == 25) & (mom_id2 == 25))
        cutchan_gen_out_2e2mu = sel_l0l1 | sel_l0l2 | sel_l1l2
        cutchan_gen_out = cutchan_gen_out_4e | cutchan_gen_out_4mu | cutchan_gen_out_2e2mu
    elif channel == "4e":
        cutchan_gen_out = ((genz_idx0 == 11) & (genz_idx1 == 11) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 11) & (genz_idx2 == 11) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 11) & (genz_idx2 == 11) & (mom_id1 == 25) & (mom_id2 == 25))
    elif channel == "4mu":
        cutchan_gen_out = ((genz_idx0 == 13) & (genz_idx1 == 13) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 13) & (genz_idx2 == 13) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 13) & (genz_idx2 == 13) & (mom_id1 == 25) & (mom_id2 == 25))
    elif channel == "2e2mu":
        sel_l0l1 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (genz_idx0 != genz_idx1))
        sel_l0l2 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx0 != genz_idx2))
        sel_l1l2 = ((mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx1 != genz_idx2))
        cutchan_gen_out = sel_l0l1 | sel_l0l2 | sel_l1l2
        
    return cutchan_gen_out

def get_chan_cuts(channel, tree):

    genz_idx0 = abs(tree['GENZ_DaughtersId'][:,0])
    genz_idx1 = abs(tree['GENZ_DaughtersId'][:,1])

    if channel == "4l":
        cutchan_gen_out = (((genz_idx0 == 11) | (genz_idx0 == 13)) & ((genz_idx1 == 11) | (genz_idx1 == 13)))
    elif channel == "4e":
        cutchan_gen_out = (genz_idx0 == 11) & (genz_idx1 == 11)
    elif channel == "4mu":
        cutchan_gen_out = (genz_idx0 == 13) & (genz_idx1 == 13)
    elif channel == "2e2mu":
        cutchan_gen_out = (((genz_idx0 == 11) & (genz_idx1 == 13)) | ((genz_idx0 == 13) & (genz_idx1 == 11)))

    return cutchan_gen_out

def m4l_cuts(tree, obs_gen, obs_gen_low, obs_gen_high):
    cutm4l_gen = (tree['GENmass4l'] > m4l_low) & (tree['GENmass4l'] < m4l_high)
    
    if 'rapidity4l' in obs_gen:
        obs_val = abs(tree[obs_gen])
    else:
        obs_val = tree[obs_gen]

    cutobs_gen = (obs_val >= obs_gen_low) & (obs_val < obs_gen_high)

    return cutm4l_gen, cutobs_gen

def get_qcd_weights(qcd_weights):
    w_nom = qcd_weights[:, 4]
    w_scale = {}
    for i in range(9):
        if ((i == 2) | (i == 4) | (i == 6)): continue
        w_scale[i] = qcd_weights[:, i]

    return w_nom, w_scale

def get_pdf_weights(pdf_weights):
    w_nom = pdf_weights[:, 0]
    w_scale = {}
    for i in range(1,101):
        w_scale[i] = pdf_weights[:, i]
        
    w_as = {
            0: pdf_weights[:, -1],
            1: pdf_weights[:, -2]
           }

    return w_nom, w_scale, w_as
    
def build_histos(tree, sel, w_gen, w_nom, w_scale):
    h_nom = ak.sum(w_gen[sel] * w_nom[sel])
    h_scale = {}

    for i in w_scale:
        h_scale[i] = ak.sum(w_gen[sel] * w_scale[i][sel])

    return h_nom, h_scale

def get_scale_unc(process, channel, tree, scale, observable, nnlops = False):
    obs_gen, obs_gen_low, obs_gen_high = observable
    cutm4l_gen, cutobs_gen = m4l_cuts(tree, obs_gen, obs_gen_low, obs_gen_high)
    if process == "ZH125":
        cutchan_gen = get_zh_chan_cuts(channel, tree)
    else:
        cutchan_gen = get_chan_cuts(channel, tree)
        
    full_sel = cutm4l_gen & cutobs_gen & cutchan_gen & tree['passedFiducial'] == 1
    
    w_gen = tree['genHEPMCweight']
    
    if nnlops:
        w_nnlops = tree['ggH_NNLOPS_weight']
        w_gen = w_gen * w_nnlops
        
    if scale == "qcd":
        w_name = "LHEScaleWeight"
        scale_weights = tree[w_name]
        w_nom, w_scale = get_qcd_weights(scale_weights)
        w_var = w_scale
    elif ((scale == "pdf") or (scale == "as")):
        w_name = "LHEPdfWeight"
        scale_weights = tree[w_name]
        w_nom, w_scale, w_as = get_pdf_weights(scale_weights)
        if scale == "pdf":
            w_var = w_scale
        else:
            w_var = w_as
    else:
        assert '... Unsupported type of uncertainty! Supported types are "qcd", "pdf", or "as" ...'
        
    h_nom, h_scale = build_histos(tree, full_sel, w_gen, w_nom, w_var)
    h_tot, h_tot_scale = build_histos(tree, cutchan_gen, w_gen, w_nom, w_var)
    
    acc_nom = h_nom/h_tot

    scale_vars = {}
    for i in h_scale:
        acc_scale = h_scale[i]/h_tot_scale[i]
        scale_vars[i] = acc_scale/acc_nom
        
    if scale == "pdf":
        scale_vars = np.array(list(scale_vars.values()))*acc_nom
        scale_unc = np.sqrt(np.sum(np.square(scale_vars - acc_nom)))
        unc_up = (acc_nom + scale_unc)/acc_nom
        unc_dn = (acc_nom - scale_unc)/acc_nom
    else:
        unc_up = max(scale_vars.values())
        unc_dn = min(scale_vars.values())
        
    return unc_up, unc_dn

def load_tree(fname):
    with uproot.open(f'{fname}') as f:
        tree = f['ZZTree/candTree'].arrays()
        tree_failed = f['ZZTree/candTree_failed'].arrays()

    tree_tot = ak.concatenate([tree, tree_failed])
    return tree_tot

def get_th_xsec(process, obs_gen, suffix):
    # TODO: Improve
    obs_name_dict = {"GENmass4l": "mass4l", "GENpT4l": "PTH", "GENrapidity4l": "YH"}
    obs_name = obs_name_dict[obs_gen]
    
    if (('ZH' not in process) and ('W' not in process)):
        th_xs = __import__(f'fidXS_{suffix}{obs_name}_{process.split("125")[0]}', globals(), locals(), vars)
    else:
        th_xs = __import__(f'fidXS_{suffix}{obs_name}_VH', globals(), locals(), vars)
    return th_xs
        
def get_unc_dict(obs_gen, tree, th_xs, channel, bins, nnlops):
    unc_var_dn = {}
    unc_var_up = {}
    for var in ["qcd", "pdf", "as"]:
        unc_bin_dn = {}
        unc_bin_up = {}
        for idx in range(len(bins)-1):
            obs_gen_low = bins[idx]
            obs_gen_high = bins[idx+1]
            if "mass4l" in obs_gen:
                obs_gen_low = 105
                obs_gen_high = 160
            observable = [obs_gen, obs_gen_low, obs_gen_high]

            var_up, var_dn = get_scale_unc(process, channel, tree, var, observable, nnlops)

            up = 1 + (var_up-1) + YR4_UNC[process][var]['up']
            # Correlate
            # np.sqrt((var_up-1)**2 + YR4_UNC[process][var]['up']**2)
            dn = 1 - ((1-var_dn) + YR4_UNC[process][var]['dn'])
            # Correlate
            # np.sqrt((1-var_dn)**2 + YR4_UNC[process][var]['dn']**2)

            print(f'({var}) {obs_gen}, bin_{idx} ({channel}): {up*th_xs.fidXS[idx]}/{dn*th_xs.fidXS[idx]}')

            unc_bin_up[idx] = up*th_xs.fidXS[idx]
            unc_bin_dn[idx] = dn*th_xs.fidXS[idx]
        
        unc_var_up[var] = unc_bin_up
        unc_var_dn[var] = unc_bin_dn
    
    return unc_var_up, unc_var_dn

def get_uncerainties(obs_gen, process, channel, bins, nnlops):
    unc = {}
    fname = f"/eos/user/m/mbonanom/Run3RedTrees/lheWeights/2022EE/{process}_MC_2022EE_skimmed.root"
    
    if (process == "ggH125") and nnlops:
        suffix = "NNLOPS_"
    else:
        suffix = ""

    tree = load_tree(fname)
    
    th_xs = get_th_xsec(process, obs_gen, suffix)
    
    unc_var_up, unc_var_dn = get_unc_dict(obs_gen, tree, th_xs, channel, bins, nnlops)

    unc[process] = (unc_var_dn, unc_var_up)
    
    return unc

def save_uncertainties(process, obs_gen, nnlops):
    obs_name = obs_gen.split('GEN')[1]
    pname = process.split("125")[0]
    if (process == "ggH125") and nnlops:
        suffix = "NNLOPS_"
    else:
        suffix = ""

    th_xs = get_th_xsec(process, obs_gen, suffix)
    with open(f'fidXS_{suffix}{obs_name}_{pname}.py', mode = 'w') as f:
        f.write(f'Boundaries = {bins}\n')
        f.write(f'fidXS = {th_xs.fidXS}\n')
        f.write(f'fidXS_scale_up = {list(unc[process][1]["qcd"].values())}\n')
        f.write(f'fidXS_scale_dn = {list(unc[process][0]["qcd"].values())}\n')
        f.write(f'fidXS_pdf_up = {list(unc[process][1]["pdf"].values())}\n')
        f.write(f'fidXS_pdf_dn = {list(unc[process][0]["pdf"].values())}\n')
        f.write(f'fidXS_alpha_up = {list(unc[process][1]["as"].values())}\n')
        f.write(f'fidXS_alpha_dn = {list(unc[process][0]["as"].values())}\n')

YR4_UNC = {'ggH125': {'qcd': {'up': 0.077, 'dn': 0.088},
                     'pdf': {'up': 0.018, 'dn': 0.018},
                     'as':{'up': 0.025, 'dn': 0.025}},
           
           'VBFH125': {'qcd': {'up': 0.004, 'dn': 0.003},
                     'pdf': {'up': 0.021, 'dn': 0.021},
                     'as':{'up': 0.005, 'dn': 0.005}},
           
           'ZH125': {'qcd': {'up': 0.038, 'dn': 0.031},
                     'pdf': {'up': 0.013, 'dn': 0.013},
                     'as':{'up': 0.009, 'dn': 0.009}},
           
           'ttH125': {'qcd': {'up': 0.058, 'dn': 0.092},
                     'pdf': {'up': 0.03, 'dn': 0.03},
                     'as':{'up': 0.02, 'dn': 0.02}},
          }



if __name__ == '__main__':
    # TODO: Improve and generalize this part
    # global obs
    m4l_low = 105
    m4l_high = 160

    obs_gen = 'GENmass4l'
    bins = [0, 2, 4, 6, 8]

    channel = '4l'
    process = 'ggH125'

    for process in ["ggH125"]:#, "VBFH125", "ttH125", "ZH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, True)
        save_uncertainties(process, obs_gen, True)

    obs_gen = 'GENrapidity4l'
    bins = [0.0, 0.15, 0.30, 0.6, 0.9, 2.5]

    for process in ["ggH125", "VBFH125", "ttH125", "ZH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, False)
        save_uncertainties(process, obs_gen, False)
        
    for process in ["ggH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, True)
        save_uncertainties(process, obs_gen, True)

    obs_gen = 'GENpT4l'
    bins = [0, 30, 80, 10000]
    for process in ["VBFH125", "ttH125", "ZH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, False)
        save_uncertainties(process, obs_gen, False)