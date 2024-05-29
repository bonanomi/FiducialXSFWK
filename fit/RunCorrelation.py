import ROOT
import sys, os, pwd
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json
from collections import OrderedDict as od
from collections import defaultdict
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


sys.path.append('../inputs/')
from higgs_xsbr_13TeV import *

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
            + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-d', '--dir',      dest='SOURCEDIR',type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--asimovModelName',dest='ASIMOVMODEL',type='string',default='SM_125', help='Name of the Asimov Model')
    parser.add_option('',   '--asimovMass',dest='ASIMOVMASS',type='string',default='125.0', help='Asimov Mass')
    parser.add_option('',   '--ModelNames',dest='MODELNAMES',type='string',default='SM_125',help='Names of models for unfolding, separated by | . Default is "SM_125"')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--fixMass',  dest='FIXMASS',  type='string',default='125.0',   help='Fix mass, default is a string "125.09" or can be changed to another string, e.g."125.6" or "False"')
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--ZZfloating',action='store_true', dest='ZZ',default=False, help='Let ZZ normalisation to float')
    # Unblind option
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    # Calculate Systematic Uncertainties
    # parser.add_option('',   '--calcSys', action='store_true', dest='SYS', default=False, help='Calculate Systematic Uncertainties (in addition to stat+sys)')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
parseOptions()

# Define function for processing of os command
def processCmd(cmd, quiet = 0):
    output = '\n'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=-1)
    for line in iter(p.stdout.readline, ''):
        output=output+str(line)
        print(line,)
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError("%r failed, exit status: %d" % (cmd, p.returncode))
    return output

def RunCombineCorrelation():
    _th_MH = opt.THEORYMASS

    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs_xs_136TeV','higgs4l_br'])
    if(opt.YEAR=="Run3"):
        higgs_xs = _temp.higgs_xs_136TeV
    else:
        higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br


    os.chdir('../combine_files/')
    # print 'Current directory: combine_files'

    for physicalModel in PhysicalModels:
        if physicalModel == 'v2': # In this case implemented for mass4l only
            cmd = 'combine -n _'+obsName+'_'+physicalModel+' -M MultiDimFit ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v2.root -m 125.38 --freezeParameters MH --floatOtherPOIs=1 --saveWorkspace --setParameterRanges r4eBin0=0.0,2.5:r4muBin0=0.0,2.5:r2e2muBin0=0.0,2.5 --redefineSignalPOI r4eBin0,r4muBin0,r2e2muBin0 --algo=singles --cminDefaultMinimizerStrategy 0 --saveInactivePOI=1 --robustHesse 1 --robustHesseSave 1'
            if not opt.UNBLIND:
                cmd += ' -t -1 --setParameters '
                for channel in ['4e', '4mu', '2e2mu']:
                    fidxs = 0
                    fidxs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                    fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                    fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                    fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                    fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                    cmd += 'r'+channel+'Bin0='+str(round(fidxs,4))+','
                cmd = cmd[:-1]
            print(cmd, '\n')
            output = processCmd(cmd)
            # cmds.append(cmd)

        if physicalModel == 'v4':
            # ----- 2e2mu -----
            cmd = 'combine -n _'+obsName+'_'+physicalModel+' -M MultiDimFit ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v4.root -m 125.38 --freezeParameters MH --floatOtherPOIs=1 --saveWorkspace --algo=singles --cminDefaultMinimizerStrategy 0 --saveInactivePOI=1 --robustHesse 1 --robustHesseSave 1 --setParameterRanges '

            for obsBin in range(nBins):
                cmd += 'r2e2muBin'+str(obsBin)+'=0.0,2.5:r4lBin'+str(obsBin)+'=0.0,2.5:'

            cmd = cmd[:-1]
            cmd += ' --redefineSignalPOI '

            for obsBin in range(nBins):
                cmd += 'r2e2muBin'+str(obsBin)+',r4lBin'+str(obsBin)+','
            cmd = cmd[:-1]

            if not opt.UNBLIND:
                cmd += ' -t -1 --setParameters '
                for obsBin in range(nBins):
                    fidxs = 0
                    fidxs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['ggH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['VBFH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['WH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['ZH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['ttH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    if(not opt.UNBLIND): cmd += 'r2e2muBin'+str(obsBin)+'='+str(round(fidxs,4))+','

                    fidxs = 0
                    # 4e
                    fidxs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4e']*acc['ggH125_4e_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4e']*acc['VBFH125_4e_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4e']*acc['WH125_4e_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4e']*acc['ZH125_4e_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4e']*acc['ttH125_4e_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    # 4mu
                    fidxs += higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4mu']*acc['ggH125_4mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4mu']*acc['VBFH125_4mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4mu']*acc['WH125_4mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4mu']*acc['ZH125_4mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_4mu']*acc['ttH125_4mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    if(not opt.UNBLIND): cmd += 'r4lBin'+str(obsBin)+'='+str(round(fidxs,4))+','
                cmd = cmd[:-1]

            print(cmd, '\n')
            output = processCmd(cmd)
            # cmds.append(cmd)


        elif physicalModel == 'v3':
            _obsName = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta4p7': 'NJ'}
            if obsName not in _obsName:
                _obsName[obsName] = obsName
            fitName = _obsName[obsName]

            cmd = 'combine -n _'+obsName+'_'+physicalModel+' -M MultiDimFit ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125.38 --freezeParameters MH --floatOtherPOIs=1 --saveWorkspace --algo=singles --cminDefaultMinimizerStrategy 0 --robustHesse 1 --robustHesseSave 1 --setParameterRanges '
            for obsBin in range(nBins):
                cmd += 'r_smH_'+fitName+'_'+str(obsBin)+'=0.0,5.0:'
            cmd = cmd[:-1]
            cmd += ' --redefineSignalPOI '
            for obsBin in range(nBins):
                cmd += 'r_smH_'+fitName+'_'+str(obsBin)+','
            cmd = cmd[:-1]

            if not opt.UNBLIND:
                cmd += ' -t -1 --setParameters '
                XH = []
                for obsBin in range(nBins):
                    # XH.append(0.0)
                    # for channel in ['4e','4mu','2e2mu']:
                    #     XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    #     XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    #     XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    #     XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    #     XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    #     XH[obsBin]+=XH_fs
                    #
                    # _obsxsec = XH[obsBin]

                    cmd += 'r_smH_'+fitName+'_'+str(obsBin)+'=1,'
                cmd = cmd[:-1]
            print(cmd, '\n')
            output = processCmd(cmd)
            # cmds.append(cmd)

        # processCmd('rm ../combine_files/robustHesse_'+obsName+'_'+physicalModel+'.root')
        # processCmd('mv robustHesse_'+obsName+'_'+physicalModel+'.root ../combine_files/.')


if 'vs' in opt.OBSNAME:
    obsName_tmp = opt.OBSNAME.split(' vs ')
    obsName = obsName_tmp[0]+'_'+obsName_tmp[1]
    doubleDiff = True
else:
    obsName = opt.OBSNAME
    doubleDiff = False

DataModelName = 'SM_125'
if obsName.startswith("mass4l"):
    PhysicalModels = ['v2','v3']
elif obsName == 'D0m' or obsName == 'Dcp' or obsName == 'D0hp' or obsName == 'Dint' or obsName == 'DL1' or obsName == 'DL1Zg' or obsName == 'costhetaZ1' or obsName == 'costhetaZ2'or obsName == 'costhetastar' or obsName == 'phi' or obsName == 'phistar' or obsName == 'massZ1' or obsName == 'massZ2':
    PhysicalModels = ['v3','v4']
else:
    PhysicalModels = ['v3']

# prepare the set of bin boundaries to run over, it is retrieved from inputs file
# _temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins', 'acc'], -1)
_temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins'])
observableBins = _temp.observableBins
_temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['acc'])
acc = _temp.acc
# print 'Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n'
# print 'Theory xsec and BR at MH = '+_th_MH
# print 'Current directory: python'

nBins = len(observableBins)
if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries

RunCombineCorrelation()

sys.path.remove('../inputs/')
