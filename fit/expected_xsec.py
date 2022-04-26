import ROOT
import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json

sys.path.append('../inputs/')

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

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()



# parse the arguments and options
global opt, args
parseOptions()

def exp_xsec():
    # prepare the set of bin boundaries to run over, only 1 bin in case of the inclusive measurement
    observableBins, doubleDiff = binning(opt.OBSNAME)
    ## Run for the given observable
    obsName = opt.OBSNAME
    _th_MH = opt.THEORYMASS
    print 'Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n'
    print 'Theory xsec and BR at MH = '+_th_MH

    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs4l_br'], -1)
    higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br
    fname = 'inputs_sig_extrap_'+obsName+'_'+opt.YEAR
    if opt.DOHIG: fname = fname + '_HIG19001'
    _temp = __import__(fname, globals(), locals(), ['acc'], -1)
    acc = _temp.acc
    XH = []
    nBins = len(observableBins)
    xs = {}
    for obsBin in range(nBins-1):
        XH.append(0.0)
        # if('mass4l' not in obsName):
        for channel in ['4e','4mu','2e2mu']:
            XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH[obsBin]+=XH_fs
        # else:
        #     XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+'4l']*acc['ggH125_4l_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
        #     XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+'4l']*acc['VBFH125_4l_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
        #     XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+'4l']*acc['WH125_4l_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
        #     XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+'4l']*acc['ZH125_4l_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
        #     XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+'4l']*acc['ttH125_4l_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
        #     XH[obsBin]+=XH_fs

        _obsxsec = XH[obsBin]

        print 'Bin ', obsBin, '\t SigmaBin', obsBin, ' = ', _obsxsec
        xs['SigmaBin'+str(obsBin)] = _obsxsec

    with open('../inputs/xsec_'+obsName+'.py', 'w') as f:
        f.write('xsec = '+str(xs)+' \n')

exp_xsec()
