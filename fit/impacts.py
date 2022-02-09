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
from createXSworkspace import createXSworkspace
from createDatacard import createDatacard

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
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='fix the fractions of 4e and 4mu when extracting the results, default is False')
    parser.add_option('',   '--physicsModel',dest='PHYSICSMODEL',type='string',default='v3',help='In case of mass4l specify explicitly v2, physicsModel to calculate impacts plots for r2e2mu,r4e,r4mu')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')


    # Unblind option
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
    #     parser.error('Bin boundaries not specified for differential measurement. Exiting...')
    #     sys.exit()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

def checkDir(folder_path):
    isdir = os.path.isdir(folder_path)
    if not isdir:
        print('Directory {} does not exist. Creating it.' .format(folder_path))
        os.mkdir(folder_path)

# Define function for processing of os command
def processCmd(cmd, quiet = 0):
    output = '\n'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=-1)
    for line in iter(p.stdout.readline, ''):
        output=output+str(line)
        print line,
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError("%r failed, exit status: %d" % (cmd, p.returncode))
    return output

def impactPlots():

    sys.path.append('../inputs/')
    _temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
    observableBins = _temp.observableBins
    sys.path.remove('../inputs/')
    ## Run for the given observable
    print 'NP impacts calculation - '+obsName+' - bin boundaries: ', observableBins, '\n'

    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs4l_br'], -1)
    higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br

    _temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['acc'], -1)
    acc = _temp.acc


    # Impact plot
    checkDir('../impacts/')
    os.chdir('../impacts/')
    print 'Current directory: impacts'
    nBins = len(observableBins)
    if '_' in obsName and obsName!='mass4l_zzfloating': nBins = len(observableBins)+1
    #nBins = len(observableBins)
    tmp_xs = {}
    tmp_xs_sm = {}
    xsec = []
    _th_MH = opt.THEORYMASS
    if opt.PHYSICSMODEL=='v3':
        cmd_XSEC =''
        for obsBin in range(nBins-1):
            for channel in ['4e','4mu','2e2mu']:
                fidxs_sm = 0
                fidxs_sm += higgs_xs['ggH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += higgs_xs['VBF_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += higgs_xs['WH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += higgs_xs['ZH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += higgs_xs['ttH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                fidxs = 0
                fidxs += higgs_xs['ggH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['VBF_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['WH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ZH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ttH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                tmp_xs_sm[channel+'_genbin'+str(obsBin)] = fidxs_sm
                tmp_xs[channel+'_genbin'+str(obsBin)] = fidxs
            cmd_XSEC += 'SigmaBin'+str(obsBin)+'='+str(tmp_xs['2e2mu_genbin'+str(obsBin)]+tmp_xs['4e_genbin'+str(obsBin)]+tmp_xs['4mu_genbin'+str(obsBin)])+','
        cmd_XSEC = cmd_XSEC[:-1]

        cmd_BR = ''
        for obsBin in range(nBins-1):
            fidxs4e = tmp_xs['4e_genbin'+str(obsBin)]
            fidxs4mu = tmp_xs['4mu_genbin'+str(obsBin)]
            fidxs2e2mu = tmp_xs['2e2mu_genbin'+str(obsBin)]
            frac4e = fidxs4e/(fidxs4e+fidxs4mu+fidxs2e2mu)
            frac4mu = fidxs4mu/(fidxs4e+fidxs4mu+fidxs2e2mu)
            fidxs4e_sm = tmp_xs_sm['4e_genbin'+str(obsBin)]
            fidxs4mu_sm = tmp_xs_sm['4mu_genbin'+str(obsBin)]
            fidxs2e2mu_sm = tmp_xs_sm['2e2mu_genbin'+str(obsBin)]
            frac4e_sm = fidxs4e_sm/(fidxs4e_sm+fidxs4mu_sm+fidxs2e2mu_sm)
            frac4mu_sm = fidxs4mu_sm/(fidxs4e_sm+fidxs4mu_sm+fidxs2e2mu_sm)
            K1 = frac4e/frac4e_sm
            K2 = frac4mu/frac4mu_sm * (1.0-frac4e_sm)/(1.0-frac4e)

            cmd_BR += 'K1Bin'+str(obsBin)+'='+str(K1)+',K2Bin'+str(obsBin)+'='+str(K2)+','
        print(cmd_BR)

        cmd_sigma = ''
        for obsBin in range(nBins-1):
            cmd_sigma += 'SigmaBin'+str(obsBin)+','
        cmd_sigma = cmd_sigma[:-1]
        print(cmd_sigma)

    elif opt.PHYSICSMODEL=='v2': #This model is used only for mass4l
        for channel in ['4e','4mu','2e2mu']:
            fidxs_sm = 0
            fidxs_sm += higgs_xs['ggH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs_sm += higgs_xs['VBF_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs_sm += higgs_xs['WH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs_sm += higgs_xs['ZH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs_sm += higgs_xs['ttH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin0_recobin0']

            fidxs = 0
            fidxs += higgs_xs['ggH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs += higgs_xs['VBF_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs += higgs_xs['WH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs += higgs_xs['ZH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            fidxs += higgs_xs['ttH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin0_recobin0']

            tmp_xs_sm[channel+'_genbin0'] = fidxs_sm
            tmp_xs[channel+'_genbin0'] = fidxs
        cmd_XSEC = 'r2e2muBin0='+str(tmp_xs['2e2mu_genbin0'])+',r4muBin0='+str(tmp_xs['4mu_genbin0'])+',r4eBin0='+str(tmp_xs['4e_genbin0'])

        cmd_sigma = 'r2e2muBin0,r4eBin0,r4muBin0'
        print(cmd_sigma)

    elif opt.PHYSICSMODEL=='v4':
        cmd_XSEC =''
        for obsBin in range(nBins-1):
            for channel in ['4e','4mu','2e2mu']:
                fidxs = 0
                fidxs += higgs_xs['ggH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['VBF_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['WH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ZH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ttH_'+_th_MH]*higgs4l_br[_th_MH+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                tmp_xs[channel+'_genbin'+str(obsBin)] = fidxs
            cmd_XSEC += 'r2e2muBin'+str(obsBin)+'='+str(tmp_xs['2e2mu_genbin'+str(obsBin)])+',r4lBin'+str(obsBin)+'='+str(tmp_xs['4e_genbin'+str(obsBin)]+tmp_xs['4mu_genbin'+str(obsBin)])+','
        cmd_XSEC = cmd_XSEC[:-1]

        cmd_BR = ''

        cmd_sigma = ''
        for obsBin in range(nBins-1):
            cmd_sigma += 'r2e2muBin'+str(obsBin)+',r4lBin'+str(obsBin)+','
        cmd_sigma = cmd_sigma[:-1]
        print(cmd_sigma)

    if (obsName.startswith("mass4l")): max_sigma = '5'
    # else: max_sigma = '2.5'
    else: max_sigma = '5'

    ### First step (Files from asimov and data have the same name)
    cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_'+opt.PHYSICSMODEL+'.root -m 125.38 --cminDefaultMinimizerStrategy 0 --doInitialFit --robustFit 1 --redefineSignalPOIs '
    for obsBin in range(nBins-1):
        if opt.PHYSICSMODEL=='v3':
            cmd += 'SigmaBin' + str(obsBin) + ','
        elif opt.PHYSICSMODEL=='v2':
            cmd += 'r2e2muBin' + str(obsBin) + ',r4eBin' + str(obsBin) +',r4muBin' + str(obsBin) +','
        elif opt.PHYSICSMODEL=='v4':
            cmd += 'r2e2muBin' + str(obsBin) + ',r4lBin' + str(obsBin) +','
    cmd = cmd[:-1]
    cmd += ' --setParameterRanges MH=125.38,125.38'
    for obsBin in range(nBins-1):
        if opt.PHYSICSMODEL=='v3':
            cmd += ':SigmaBin' + str(obsBin) + '=0,'+max_sigma
        elif opt.PHYSICSMODEL=='v2':
            cmd += ':r2e2muBin' + str(obsBin) + '=0,'+max_sigma+':r4muBin' + str(obsBin) + '=0,'+max_sigma+':r4eBin' + str(obsBin) + '=0,'+max_sigma
        elif opt.PHYSICSMODEL=='v4':
            cmd += ':r2e2muBin' + str(obsBin) + '=0,'+max_sigma+':r4lBin' + str(obsBin) + '=0,'+max_sigma
    if (not opt.UNBLIND):
        if opt.PHYSICSMODEL=='v3':
            cmd = cmd + ' -t -1 --setParameters MH=125.38,' + cmd_BR[:-1] + ',' + cmd_XSEC
        elif opt.PHYSICSMODEL=='v2':
            cmd = cmd + ' -t -1 --setParameters MH=125.38,' + cmd_XSEC
        elif opt.PHYSICSMODEL=='v4':
            cmd = cmd + ' -t -1 --setParameters MH=125.38,' + cmd_XSEC
    print '---------------------------'
    print cmd, '\n'
    print '---------------------------'
    cmds.append(cmd)
    output = processCmd(cmd)

    ### Second step (Files from asimov and data have the same name)
    cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_'+opt.PHYSICSMODEL+'.root -m 125.38 --cminDefaultMinimizerStrategy 0 --doFits --robustFit 1 --parallel 10 --redefineSignalPOIs '
    for obsBin in range(nBins-1):
        if opt.PHYSICSMODEL=='v3':
            cmd += 'SigmaBin' + str(obsBin) + ','
        elif opt.PHYSICSMODEL=='v2':
            cmd += 'r2e2muBin' + str(obsBin) + ',r4eBin' + str(obsBin) +',r4muBin' + str(obsBin) +','
        elif opt.PHYSICSMODEL=='v4':
            cmd += 'r2e2muBin' + str(obsBin) + ',r4lBin' + str(obsBin) +','
    cmd = cmd[:-1]
    cmd += ' --setParameterRanges MH=125.38,125.38'
    for obsBin in range(nBins-1):
        if opt.PHYSICSMODEL=='v3':
            cmd += ':SigmaBin' + str(obsBin) + '=0,'+max_sigma
        elif opt.PHYSICSMODEL=='v2':
            cmd += ':r2e2muBin' + str(obsBin) + '=0,'+max_sigma+':r4muBin' + str(obsBin) + '=0,'+max_sigma+':r4eBin' + str(obsBin) + '=0,'+max_sigma
        elif opt.PHYSICSMODEL=='v4':
            cmd += ':r2e2muBin' + str(obsBin) + '=0,'+max_sigma+':r4lBin' + str(obsBin) + '=0,'+max_sigma
    if (not opt.UNBLIND):
        if opt.PHYSICSMODEL=='v3':
            cmd = cmd + ' -t -1 --setParameters MH=125.38,' + cmd_BR[:-1] + ',' + cmd_XSEC
        elif opt.PHYSICSMODEL=='v2':
            cmd = cmd + ' -t -1 --setParameters MH=125.38,' + cmd_XSEC
        elif opt.PHYSICSMODEL=='v4':
            cmd = cmd + ' -t -1 --setParameters MH=125.38,' + cmd_XSEC
    print '---------------------------'
    print cmd, '\n'
    print '---------------------------'
    cmds.append(cmd)
    output = processCmd(cmd)

    ### Third step
    if opt.PHYSICSMODEL=='v3':
        for obsBin in range(nBins-1):
            cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125.38 -t -1 --setParameters MH=125.38,'+cmd_BR[:-1]+','+cmd_XSEC+' --redefineSignalPOIs '+cmd_sigma
            cmd += ' -o impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_'
            if (not opt.UNBLIND):
                cmd = cmd + 'asimov.json'
            elif (opt.UNBLIND):
                cmd = cmd + 'data.json'
            print '---------------------------'
            print cmd, '\n'
            print '---------------------------'
            cmds.append(cmd)
            output = processCmd(cmd)
            # plot
            cmd = 'plotImpacts.py -i impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_'
            if (not opt.UNBLIND): cmd = cmd + 'asimov.json -o impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_asimov --POI SigmaBin'+str(obsBin)
            elif (opt.UNBLIND): cmd = cmd + 'data.json -o impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_data --POI SigmaBin'+str(obsBin)
            print '---------------------------'
            print cmd, '\n'
            print '---------------------------'
            cmds.append(cmd)
            output = processCmd(cmd)

    elif opt.PHYSICSMODEL=='v2':
        for obsBin in ['2e2muBin0','4eBin0','4muBin0']:
            cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v2.root -m 125.38 -t -1 --setParameters MH=125.38,'+cmd_XSEC+' --redefineSignalPOIs '+cmd_sigma
            cmd += ' -o impacts_v2_'+obsName+'_r'+str(obsBin)+'_'
            if (not opt.UNBLIND):
                cmd = cmd + 'asimov.json'
            elif (opt.UNBLIND):
                cmd = cmd + 'data.json'
            print '---------------------------'
            print cmd, '\n'
            print '---------------------------'
            cmds.append(cmd)
            output = processCmd(cmd)
            # plot
            cmd = 'plotImpacts.py -i impacts_v2_'+obsName+'_r'+str(obsBin)+'_'
            if (not opt.UNBLIND): cmd = cmd + 'asimov.json -o impacts_v2_'+obsName+'_r'+str(obsBin)+'_asimov --POI r'+str(obsBin)
            elif (opt.UNBLIND): cmd = cmd + 'data.json -o impacts_v2_'+obsName+'_r'+str(obsBin)+'_data --POI r'+str(obsBin)
            print '---------------------------'
            print cmd, '\n'
            print '---------------------------'
            cmds.append(cmd)
            output = processCmd(cmd)

    elif opt.PHYSICSMODEL=='v4':
        for nBin in range(nBins-1):
            for obsBin in ['2e2muBin'+str(nBin),'4lBin'+str(nBin)]:
                cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v4.root -m 125.38 -t -1 --cminDefaultMinimizerStrategy 0 --setParameters MH=125.38,'+cmd_XSEC+' --redefineSignalPOIs '+cmd_sigma
                cmd += ' -o impacts_v4_'+obsName+'_r'+str(obsBin)+'_'
                if (not opt.UNBLIND):
                    cmd = cmd + 'asimov.json'
                elif (opt.UNBLIND):
                    cmd = cmd + 'data.json'
                print '---------------------------'
                print cmd, '\n'
                print '---------------------------'
                cmds.append(cmd)
                output = processCmd(cmd)
                # plot
                cmd = 'plotImpacts.py -i impacts_v4_'+obsName+'_r'+str(obsBin)+'_'
                if (not opt.UNBLIND): cmd = cmd + 'asimov.json -o impacts_v4_'+obsName+'_r'+str(obsBin)+'_asimov --POI r'+str(obsBin)
                elif (opt.UNBLIND): cmd = cmd + 'data.json -o impacts_v4_'+obsName+'_r'+str(obsBin)+'_data --POI r'+str(obsBin)
                print '---------------------------'
                print cmd, '\n'
                print '---------------------------'
                cmds.append(cmd)
                output = processCmd(cmd)


# ----------------- Main -----------------
cmds = [] #List of all cmds
if 'vs' in opt.OBSNAME:
    obsName_tmp = opt.OBSNAME.split(' vs ')
    obsName = obsName_tmp[0]+'_'+obsName_tmp[1]
else:
    obsName = opt.OBSNAME
impactPlots()
if (os.path.exists('commands_impacts_'+obsName+'_'+opt.PHYSICSMODEL+'.py')):
    os.system('rm commands_impacts_'+obsName+'_'+opt.PHYSICSMODEL+'.py')
with open('commands_impacts_'+obsName+'_'+opt.PHYSICSMODEL+'.py', 'w') as f:
    for i in cmds:
        f.write(str(i)+' \n')
        f.write('\n')

print "Impacts plots done."
