import ROOT
import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json

from higgs_xsbr_13TeV import *
from binning import binning

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-d', '--dir',      dest='SOURCEDIR',type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--doHIG', action='store_true', dest='DOHIG', default=False, help='use HIG 19 001 acceptances')
    parser.add_option('',   '--nnlops', action='store_true', dest='NNLOPS', default=False, help='nnlops prediction')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()



# parse the arguments and options
global opt, args
parseOptions()

def exp_xsec():
    # prepare the set of bin boundaries to run over, only 1 bin in case of the inclusive measurement
    observableBins, doubleDiff = binning(opt.OBSNAME)
    if doubleDiff:
        obs_name = opt.OBSNAME.split(' vs ')[0]
        obs_name_2nd = opt.OBSNAME.split(' vs ')[1]
        obs_name_2d = opt.OBSNAME
        doubleDiff = True
    else:
        obs_name = opt.OBSNAME
        doubleDiff = False
    ## Run for the given observable
    # obsName = opt.OBSNAME
    if doubleDiff: obsName = obs_name+'_'+obs_name_2nd
    else: obsName = obs_name
    _th_MH = opt.THEORYMASS
    print 'Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n'
    print 'Theory xsec and BR at MH = '+_th_MH

    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs_xs_136TeV','higgs4l_br'], -1)
    if(opt.YEAR=='Run3'):
        higgs_xs = _temp.higgs_xs_136TeV
    else:
        higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br
    fname = 'inputs_sig_'+obsName+'_'+opt.YEAR
    if opt.DOHIG: fname = fname + '_HIG19001'
    _temp = __import__(fname, globals(), locals(), ['acc'], -1)
    acc = _temp.acc
    if opt.NNLOPS:
        fname = 'inputs_sig_'+obsName+'_NNLOPS_'+opt.YEAR
        _temp = __import__(fname, globals(), locals(), ['acc'], -1)
        acc_ggh = _temp.acc
    else:
        acc_ggh = acc

    if opt.NNLOPS:
        suffix = 'NNLOPS_'
    else:
        suffix = ''

    XH = []
    XH_ggh = []
    XH_vbf = []
    XH_zh = []
    XH_wh = []
    XH_ttH = []
    nBins = len(observableBins)
    if not doubleDiff: nBins = nBins - 1
    xs = {}
    xs_ggh = {}
    xs_vbf = {}
    xs_zh = {}
    xs_wh = {}
    xs_vh = {}
    xs_tth = {}
    xs_xh = {}

    for obsBin in range(nBins):
        XH.append(0.0)
        XH_ggh.append(0.0)
        XH_vbf.append(0.0)
        XH_zh.append(0.0)
        XH_wh.append(0.0)
        XH_ttH.append(0.0)

        for channel in ['4e','4mu','2e2mu']:
            print(channel)
            print(acc_ggh['ggH125_'+suffix+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])
            xxs_ggh = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_ggh['ggH125_'+suffix+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_vbf = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_zh = higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_tth = higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_wh = higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

            XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_ggh['ggH125_'+suffix+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH[obsBin]+=XH_fs

            XH_ggh[obsBin]+=xxs_ggh
            XH_vbf[obsBin]+=xxs_vbf
            XH_zh[obsBin]+=xxs_zh
            XH_wh[obsBin]+=xxs_wh
            XH_ttH[obsBin]+=xxs_tth

        _obsxsec = XH[obsBin]

        print '\n'
        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' = ', _obsxsec
        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' (ggh) = ', XH_ggh[obsBin]
        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' (vbf) = ', XH_vbf[obsBin]
        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' (zh) = ', XH_zh[obsBin]
        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' (wh) = ', XH_wh[obsBin]
        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' (tth) = ', XH_ttH[obsBin]
        print '(xcheck : ', XH_ggh[obsBin]+XH_vbf[obsBin]+XH_zh[obsBin]+XH_wh[obsBin]+XH_ttH[obsBin],')'
        print '(XH : ', XH_vbf[obsBin]+XH_zh[obsBin]+XH_wh[obsBin]+XH_ttH[obsBin],')'
        print '\n\n'
        xs['SigmaBin'+str(obsBin)] = _obsxsec
        xs_ggh['SigmaBin'+str(obsBin)] = XH_ggh[obsBin]
        xs_vbf['SigmaBin'+str(obsBin)] = XH_vbf[obsBin]
        xs_vh['SigmaBin'+str(obsBin)] = XH_zh[obsBin]+XH_wh[obsBin]
        xs_wh['SigmaBin'+str(obsBin)] = XH_wh[obsBin]
        xs_tth['SigmaBin'+str(obsBin)] = XH_ttH[obsBin]
        
        xs_xh['SigmaBin'+str(obsBin)] = XH_vbf[obsBin]+XH_zh[obsBin]+XH_wh[obsBin]+XH_ttH[obsBin]

    with open('../inputs/xsec_'+obsName+'.py', 'w') as f:
        f.write('xsec = '+str(xs)+' \n')


    if obsName == 'pT4l':
        obsFull = 'PTH'

    if obsName == 'rapidity4l':
        obsFull = 'YH'
    else:
        obsFull = obsName

    with open('../inputs/fidXS_'+suffix+obsFull+'_ggH.py', 'w') as f:
        f.write('Boundaries = '+str(observableBins)+'\n')
        f.write('fidXS = '+str(xs_ggh.values())+'\n')
        f.write('fidXS_scale_up = '+str(xs_ggh.values())+'\n')
        f.write('fidXS_scale_dn = '+str(xs_ggh.values())+'\n')
        f.write('fidXS_pdf_up = '+str(xs_ggh.values())+'\n')
        f.write('fidXS_pdf_dn = '+str(xs_ggh.values())+'\n')
        f.write('fidXS_alpha_up = '+str(xs_ggh.values())+'\n')
        f.write('fidXS_alpha_dn = '+str(xs_ggh.values())+'\n')

    with open('../inputs/fidXS_'+obsFull+'_VBFH.py', 'w') as f:
        f.write('Boundaries = '+str(observableBins)+'\n')
        f.write('fidXS = '+str(xs_vbf.values())+'\n')
        f.write('fidXS_scale_up = '+str(xs_vbf.values())+'\n')
        f.write('fidXS_scale_dn = '+str(xs_vbf.values())+'\n')
        f.write('fidXS_pdf_up = '+str(xs_vbf.values())+'\n')
        f.write('fidXS_pdf_dn = '+str(xs_vbf.values())+'\n')
        f.write('fidXS_alpha_up = '+str(xs_vbf.values())+'\n')
        f.write('fidXS_alpha_dn = '+str(xs_vbf.values())+'\n')

    with open('../inputs/fidXS_'+obsFull+'_VH.py', 'w') as f:
        f.write('Boundaries = '+str(observableBins)+'\n')
        f.write('fidXS = '+str(xs_vh.values())+'\n')
        f.write('fidXS_scale_up = '+str(xs_vh.values())+'\n')
        f.write('fidXS_scale_dn = '+str(xs_vh.values())+'\n')
        f.write('fidXS_pdf_up = '+str(xs_vh.values())+'\n')
        f.write('fidXS_pdf_dn = '+str(xs_vh.values())+'\n')
        f.write('fidXS_alpha_up = '+str(xs_vh.values())+'\n')
        f.write('fidXS_alpha_dn = '+str(xs_vh.values())+'\n')

    with open('../inputs/fidXS_'+obsFull+'_ttH.py', 'w') as f:
        f.write('Boundaries = '+str(observableBins)+'\n')
        f.write('fidXS = '+str(xs_tth.values())+'\n')
        f.write('fidXS_scale_up = '+str(xs_tth.values())+'\n')
        f.write('fidXS_scale_dn = '+str(xs_tth.values())+'\n')
        f.write('fidXS_pdf_up = '+str(xs_tth.values())+'\n')
        f.write('fidXS_pdf_dn = '+str(xs_tth.values())+'\n')
        f.write('fidXS_alpha_up = '+str(xs_tth.values())+'\n')
        f.write('fidXS_alpha_dn = '+str(xs_tth.values())+'\n')

    with open('../inputs/fidXS_'+obsFull+'_xH.py', 'w') as f:
        f.write('Boundaries = '+str(observableBins)+'\n')
        f.write('fidXS = '+str(xs_xh.values())+'\n')
        f.write('fidXS_scale_up = '+str(xs_xh.values())+'\n')
        f.write('fidXS_scale_dn = '+str(xs_xh.values())+'\n')
        f.write('fidXS_pdf_up = '+str(xs_xh.values())+'\n')
        f.write('fidXS_pdf_dn = '+str(xs_xh.values())+'\n')
        f.write('fidXS_alpha_up = '+str(xs_xh.values())+'\n')
        f.write('fidXS_alpha_dn = '+str(xs_xh.values())+'\n')
exp_xsec()