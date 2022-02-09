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
    parser.add_option('',   '--fixMass',  dest='FIXMASS',  type='string',default='125.0',   help='Fix mass, default is a string "125.09" or can be changed to another string, e.g."125.6" or "False"')
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='fix the fractions of 4e and 4mu when extracting the results, default is False')
    # action options - "only"
    # parser.add_option('',   '--effOnly',       action='store_true', dest='effOnly',       default=False, help='Extract the eff. factors only, default is False')
    # parser.add_option('',   '--templatesOnly', action='store_true', dest='templatesOnly', default=False, help='Prepare the bkg shapes and fractions only, default is False')
    # parser.add_option('',   '--uncertOnly',    action='store_true', dest='uncertOnly',    default=False, help='Extract the uncertanties only, default is False')
    # parser.add_option('',   '--resultsOnly',   action='store_true', dest='resultsOnly',   default=False, help='Run the measurement only, default is False')
    # parser.add_option('',   '--finalplotsOnly',action='store_true', dest='finalplotsOnly',default=False, help='Make the final plots only, default is False')
    parser.add_option('',   '--impactsOnly',action='store_true', dest='impactsOnly',default=False, help='Make the impacts plots only, default is False')
    parser.add_option('',   '--combineOnly',action='store_true', dest='combineOnly',default=False, help='Run the measurement only, default is False')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='int',default=105.0,   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='int',default=140.0,   help='Upper bound for m4l')
    parser.add_option('',   '--ZZfloating',action='store_true', dest='ZZ',default=False, help='Let ZZ normalisation to float')
    # Unblind option
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    # Calculate Systematic Uncertainties
    # parser.add_option('',   '--calcSys', action='store_true', dest='SYS', default=False, help='Calculate Systematic Uncertainties (in addition to stat+sys)')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # prepare the global flag if all the step should be run
    runAllSteps = not(opt.combineOnly or opt.impactsOnly)

    # if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
    #     parser.error('Bin boundaries not specified for differential measurement. Exiting...')
    #     sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

if (opt.YEAR == '2016'): years = ['2016']
if (opt.YEAR == '2017'): years = ['2017']
if (opt.YEAR == '2018'): years = ['2018']
if (opt.YEAR == 'Full'): years = ['2016','2017','2018']


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


### Produce datacards for given obs and bin, for all final states
def produceDatacards(obsName, observableBins, ModelName, physicalmodel):
    print '\n'
    print '[Producing workspace/datacards for obsName '+obsName+', bins '+str(observableBins)+']'
    fStates = ['2e2mu','4mu','4e']
    nBins = len(observableBins)
    if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries
    # if (('pTj1' in obsName) | ('pTHj' in obsName) | ('mHj' in obsName) | ('pTj2' in obsName) | ('mjj' in obsName) | ('absdetajj' in obsName) | ('dphijj' in obsName) | ('pTHjj' in obsName)  | ('TCjmax' in obsName) | ('TBjmax' in obsName) | ('njets_pt30_eta4p7' in obsName)):
    #     JES = True
    # else:
    #     JES = False
    # # JES = False
    os.chdir('../datacard/datacard_'+years[0])
    for year in years:
        os.chdir('../datacard_'+year)
        print 'Current diretory: datacard_'+year
        for fState in fStates:
            if (not obsName.startswith("mass4l")):
                for obsBin in range(nBins):
                    ndata = createXSworkspace(obsName,fState, nBins, obsBin, observableBins, False, True, ModelName, physicalmodel, year, JES, doubleDiff, opt.LOWER_BOUND, opt.UPPER_BOUND, opt.OBSNAME)
                    createDatacard(obsName, fState, nBins, obsBin, observableBins, physicalmodel, year, ndata, JES, opt.LOWER_BOUND, opt.UPPER_BOUND, opt.YEAR)
                    os.chdir('../datacard/datacard_'+year)
            else:
                ndata = createXSworkspace(obsName,fState, nBins, 0, observableBins, False, True, ModelName, physicalmodel, year, JES, doubleDiff, opt.LOWER_BOUND, opt.UPPER_BOUND, opt.OBSNAME)
                createDatacard(obsName, fState, nBins, 0, observableBins, physicalmodel, year, ndata, JES, opt.LOWER_BOUND, opt.UPPER_BOUND, opt.YEAR)
                os.chdir('../datacard/datacard_'+year)

                # if obsName=='mass4l': os.system("cp xs_125.0_1bin/hzz4l_"+fState+"S_13TeV_xs_inclusive_bin0.txt xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")
                # if obsName=='mass4lREFIT': os.system("cp xs_125.0_1bin/hzz4l_"+fState+"S_13TeV_xs_inclusiveREFIT_bin0.txt xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")
                # os.system("sed -i 's~observation [0-9]*~observation "+str(ndata)+"~g' xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")
                # os.system("sed -i 's~_xs.Databin0~_xs_"+ModelName+"_"+obsName+"_"+PhysicalModel+".Databin0~g' xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")


def runFiducialXS():
    # variable for double-differential measurements and obsName
    # global doubleDiff
    # if 'vs' in opt.OBSNAME:
    #     obsName_tmp = opt.OBSNAME.split(' vs ')
    #     obsName = obsName_tmp[0]+'_'+obsName_tmp[1]
    #     doubleDiff = True
    # else:
    #     obsName = opt.OBSNAME
    #     doubleDiff = False
    _th_MH = opt.THEORYMASS
    # prepare the set of bin boundaries to run over, it is retrieved from inputs file
    _temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
    observableBins = _temp.observableBins
    print 'Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n'
    print 'Theory xsec and BR at MH = '+_th_MH
    print 'Current directory: python'

    # _fit_dir = os.getcwd()

    ## addConstrainedModel
    years_bis = years
    if(opt.YEAR == 'Full'):
        years_bis.append('Full')
    for year in years_bis:
        if not os.path.exists('../inputs/inputs_sig_'+obsName+'_'+year+'_ORIG.py'):
            cmd = 'python addConstrainedModel.py -l -q -b --obsName="'+obsName+'" --year="'+year+'"'
            if doubleDiff: cmd += ' --doubleDiff'
            print cmd
            output = processCmd(cmd)
            cmds.append(cmd)
            print output
            print 'addConstrainedModel DONE'
        elif os.path.exists('../inputs/inputs_sig_'+obsName+'_'+year+'_ORIG.py'):
            print 'addConstrainedModel '+year+' already done'
    if 'Full' in years: years.remove('Full')
    # "__import__" to SetParameters when running the expected measurement
    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs4l_br'], -1)
    higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br
    _temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['acc'], -1)
    acc = _temp.acc

    DataModelName = 'SM_125'
    if obsName.startswith("mass4l"):
        PhysicalModels = ['v2','v3']
    elif obsName == 'D0m' or obsName == 'Dcp' or obsName == 'D0hp' or obsName == 'Dint' or obsName == 'DL1' or obsName == 'DL1Zg' or obsName == 'costhetaZ1' or obsName == 'costhetaZ2'or obsName == 'costhetastar' or obsName == 'phi' or obsName == 'phistar' or obsName == 'massZ1' or obsName == 'massZ2':
        PhysicalModels = ['v3','v4']
    else:
        PhysicalModels = ['v3']

    for physicalModel in PhysicalModels:
        produceDatacards(obsName, observableBins, DataModelName, physicalModel)

        # combination of bins (if there is just one bin, it is essentially a change of name from _bin0_ to _bin_)
        fStates = ['2e2mu','4mu','4e']
        nBins = len(observableBins)
        if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries

        os.chdir(_fit_dir)
        for year in years:
            #We are already in datacard dir at this point
            os.chdir('../datacard/datacard_'+year)
            print 'Current directory: datacard_'+year
            for fState in fStates:
                if(nBins>0):
                    cmd = 'combineCards.py '
                    for obsBin in range(nBins):
                        cmd = cmd + 'hzz4l_'+fState+'S_13TeV_xs_'+obsName+'_bin'+str(obsBin)+'_'+physicalModel+'.txt '
                    cmd = cmd + '> hzz4l_'+fState+'S_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt'
                    print cmd, '\n'
                    processCmd(cmd,1)
                    cmds.append(cmd)
                else:
                    print 'There is a problem during the combination over bins'

            # combine 3 final states
            cmd = 'combineCards.py hzz4l_4muS_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt hzz4l_4eS_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt hzz4l_2e2muS_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
            cmds.append(cmd)

            os.chdir(_fit_dir)

        # Combine 3 years
    	# we go back from datacard_Y to datacard folder
        os.chdir('../datacard/')
        print 'Current directory: datacard'
        if (opt.YEAR == 'Full'):
            cmd = 'combineCards.py datacard_2016/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt datacard_2017/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt datacard_2018/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
            cmds.append(cmd)
            if obsName == 'mass4l_zzfloating': # Remove bkg theo nuisances in case of zz floating
                cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_2016 CMS_hzz2e2mu_Zjets_2017 CMS_hzz2e2mu_Zjets_2018 CMS_hzz4e_Zjets_2016 CMS_hzz4e_Zjets_2017 CMS_hzz4e_Zjets_2018 CMS_hzz4mu_Zjets_2016 CMS_hzz4mu_Zjets_2017 CMS_hzz4mu_Zjets_2018 lumi_13TeV_2016_uncorrelated lumi_13TeV_2017_uncorrelated lumi_13TeV_2018_uncorrelated lumi_13TeV_correlated_16_17_18 lumi_13TeV_correlated_17_18 norm_fakeH CMS_zz4l_sigma_e_sig_2017 CMS_zz4l_sigma_e_sig_2016 CMS_zz4l_sigma_m_sig_2018 CMS_zz4l_sigma_m_sig_2017 CMS_zz4l_sigma_m_sig_2016 CMS_zz4l_n_sig_3_2016 CMS_zz4l_n_sig_3_2017 CMS_zz4l_mean_e_sig_2016 CMS_zz4l_mean_e_sig_2017 CMS_zz4l_n_sig_3_2018 CMS_zz4l_mean_m_sig_2018 CMS_zz4l_mean_m_sig_2016 CMS_zz4l_mean_m_sig_2017 CMS_zz4l_sigma_e_sig_2018 CMS_zz4l_mean_e_sig_2018'
            else:
                cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_2016 CMS_hzz2e2mu_Zjets_2017 CMS_hzz2e2mu_Zjets_2018 CMS_hzz4e_Zjets_2016 CMS_hzz4e_Zjets_2017 CMS_hzz4e_Zjets_2018 CMS_hzz4mu_Zjets_2016 CMS_hzz4mu_Zjets_2017 CMS_hzz4mu_Zjets_2018 QCDscale_VV QCDscale_ggVV kfactor_ggzz lumi_13TeV_2016_uncorrelated lumi_13TeV_2017_uncorrelated lumi_13TeV_2018_uncorrelated lumi_13TeV_correlated_16_17_18 lumi_13TeV_correlated_17_18 norm_fakeH pdf_gg pdf_qqbar CMS_zz4l_sigma_e_sig_2017 CMS_zz4l_sigma_e_sig_2016 CMS_zz4l_sigma_m_sig_2018 CMS_zz4l_sigma_m_sig_2017 CMS_zz4l_sigma_m_sig_2016 CMS_zz4l_n_sig_3_2016 CMS_zz4l_n_sig_3_2017 CMS_zz4l_mean_e_sig_2016 CMS_zz4l_mean_e_sig_2017 CMS_zz4l_n_sig_3_2018 CMS_zz4l_mean_m_sig_2018 CMS_zz4l_mean_m_sig_2016 CMS_zz4l_mean_m_sig_2017 CMS_zz4l_sigma_e_sig_2018 CMS_zz4l_mean_e_sig_2018'
            if JES:
                cmd += ' CMS_scale_j_Abs CMS_scale_j_Abs_2016 CMS_scale_j_BBEC1 CMS_scale_j_BBEC1_2016 CMS_scale_j_EC2 CMS_scale_j_EC2_2016 CMS_scale_j_FlavQCD CMS_scale_j_HF CMS_scale_j_HF_2016 CMS_scale_j_RelBal CMS_scale_j_RelSample_2016 CMS_scale_j_Abs CMS_scale_j_Abs_2017 CMS_scale_j_BBEC1 CMS_scale_j_BBEC1_2017 CMS_scale_j_EC2 CMS_scale_j_EC2_2017 CMS_scale_j_FlavQCD CMS_scale_j_HF CMS_scale_j_HF_2017 CMS_scale_j_RelBal CMS_scale_j_RelSample_2017 CMS_scale_j_Abs CMS_scale_j_Abs_2018 CMS_scale_j_BBEC1 CMS_scale_j_BBEC1_2018 CMS_scale_j_EC2 CMS_scale_j_EC2_2018 CMS_scale_j_FlavQCD CMS_scale_j_HF CMS_scale_j_HF_2018 CMS_scale_j_RelBal CMS_scale_j_RelSample_2018 CMS_scale_j_ZX'
            cmd += '" >> hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
            cmds.append(cmd)
        else:
            cmd = 'cp datacard_'+str(opt.YEAR)+'/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
            cmds.append(cmd)
            cmd = "sed -i 's|hzz4l|datacard_"+str(opt.YEAR)+"/hzz4l|g' hzz4l_all_13TeV_xs_"+obsName+"_bin_"+physicalModel+".txt" # Specify the right pattern to datacards (Before it was not necessary because there was a further combination)
            print cmd, '\n'
            processCmd(cmd,1)
            cmds.append(cmd)
            if obsName == 'mass4l_zzfloating':
                cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_'+str(opt.YEAR)+' CMS_hzz4e_Zjets_'+str(opt.YEAR)+' CMS_hzz4mu_Zjets_'+str(opt.YEAR)+' lumi_13TeV_'+str(opt.YEAR)+' norm_fakeH CMS_zz4l_sigma_e_sig_'+str(opt.YEAR)+' CMS_zz4l_sigma_m_sig_'+str(opt.YEAR)+' CMS_zz4l_n_sig_3_'+str(opt.YEAR)+' CMS_zz4l_mean_e_sig_'+str(opt.YEAR)+' CMS_zz4l_mean_m_sig_'+str(opt.YEAR)+' CMS_zz4l_sigma_e_sig_'+str(opt.YEAR)
            else:
                cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_'+str(opt.YEAR)+' CMS_hzz4e_Zjets_'+str(opt.YEAR)+' CMS_hzz4mu_Zjets_'+str(opt.YEAR)+' QCDscale_VV QCDscale_ggVV kfactor_ggzz lumi_13TeV_'+str(opt.YEAR)+' norm_fakeH pdf_gg pdf_qqbar CMS_zz4l_sigma_e_sig_'+str(opt.YEAR)+' CMS_zz4l_sigma_m_sig_'+str(opt.YEAR)+' CMS_zz4l_n_sig_3_'+str(opt.YEAR)+' CMS_zz4l_mean_e_sig_'+str(opt.YEAR)+' CMS_zz4l_mean_m_sig_'+str(opt.YEAR)+' CMS_zz4l_sigma_e_sig_'+str(opt.YEAR)
            if JES:
                cmd += ' CMS_scale_j_Abs CMS_scale_j_Abs_'+str(opt.YEAR)+' CMS_scale_j_BBEC1 CMS_scale_j_BBEC1_'+str(opt.YEAR)+' CMS_scale_j_EC2 CMS_scale_j_EC2_'+str(opt.YEAR)+' CMS_scale_j_FlavQCD CMS_scale_j_HF CMS_scale_j_HF_'+str(opt.YEAR)+' CMS_scale_j_RelBal CMS_scale_j_RelSample_'+str(opt.YEAR)+' CMS_scale_j_ZX'
            cmd += '" >> hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt'
            processCmd(cmd,1)
            cmds.append(cmd)
            print cmd, '\n'

        # text-to-workspace
        if (physicalModel=="v3"):
            cmd = 'text2workspace.py hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt -P HiggsAnalysis.CombinedLimit.HZZ4L_Fiducial:differentialFiducialV3 --PO higgsMassRange=115,135 --PO nBin='+str(nBins)+' -o hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.root'
            print cmd, '\n'
            processCmd(cmd)
            cmds.append(cmd)
        elif (physicalModel=="v4"):
            cmd = 'text2workspace.py hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt -P HiggsAnalysis.CombinedLimit.HZZ4L_Fiducial_v2:differentialFiducialV4 --PO higgsMassRange=115,135 --PO nBin='+str(nBins)+' -o hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.root'
            print cmd, '\n'
            processCmd(cmd)
            cmds.append(cmd)
        elif (physicalModel=="v2"):
            cmd = 'text2workspace.py hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.txt -P HiggsAnalysis.CombinedLimit.HZZ4L_Fiducial_v2:differentialFiducialV2 --PO higgsMassRange=115,135 --PO nBin='+str(nBins)+' -o hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.root'
            print cmd, '\n'
            processCmd(cmd)
            cmds.append(cmd)


        # The workspace got from text2workspace changes name from hzz4l_ to SM_125 and it is transferred to the combine_files directory
        cmd = 'cp hzz4l_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.root ../combine_files/'+DataModelName+'_all_13TeV_xs_'+obsName+'_bin_'+physicalModel+'.root'
        print cmd, '\n'
        processCmd(cmd,1)
        cmds.append(cmd)

    	# From datacard directory to combine_files, to store fit results
        os.chdir('../combine_files/')
        print 'Current directory: combine_files'
        # nBins = len(observableBins)
        if physicalModel == 'v2': # In this case implemented for mass4l only
            for channel in ['4e', '4mu', '2e2mu']:
                cmd = 'combine -n _'+obsName+'_r'+channel+'Bin0 -M MultiDimFit SM_125_all_13TeV_xs_'+obsName+'_bin_v2.root -m 125.38 --freezeParameters MH -P r'+channel+'Bin0 --floatOtherPOIs=1 --saveWorkspace --setParameterRanges r'+channel+'Bin0=0.0,2.5 --redefineSignalPOI r'+channel+'Bin0 --algo=grid --points=100 --cminDefaultMinimizerStrategy 0 --saveInactivePOI=1'

                fidxs = 0
                fidxs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys --setParameters r'+channel+'Bin0='+str(round(fidxs,4))
                print cmd, '\n'
                output = processCmd(cmd)
                cmds.append(cmd)
                # Stat-only
                cmd = 'combine -n _'+obsName+'_r'+channel+'Bin0_NoSys'
                if(not opt.UNBLIND): cmd = cmd + '_exp'
                cmd = cmd + ' -M MultiDimFit higgsCombine_'+obsName+'_r'+channel+'Bin0.MultiDimFit.mH125.38'
                if(not opt.UNBLIND): cmd = cmd + '.123456'
                cmd = cmd + '.root -w w --snapshotName "MultiDimFit" -m 125.38 -P r'+channel+'Bin0 --floatOtherPOIs=1 --saveWorkspace --setParameterRanges SigmaBin0=0.0,2.5 --redefineSignalPOI r'+channel+'Bin0 --algo=grid --points=100 --cminDefaultMinimizerStrategy 0 --freezeNuisanceGroups nuis'
                if (opt.YEAR == 'Full'): cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
                else: cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_1'+str(opt.YEAR)+',CMS_fakeH_p3_1'+str(opt.YEAR)+',CMS_fakeH_p1_2'+str(opt.YEAR)+',CMS_fakeH_p3_2'+str(opt.YEAR)+',CMS_fakeH_p1_3'+str(opt.YEAR)+',CMS_fakeH_p3_3'+str(opt.YEAR)
                if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys --setParameters r'+channel+'Bin0='+str(round(fidxs,4))
                print cmd+'\n'
                output = processCmd(cmd)
                cmds.append(cmd)

        if physicalModel == 'v4':
            for obsBin in range(nBins):
                # ----- 2e2mu -----
                cmd = 'combine -n _'+obsName+'_r2e2muBin'+str(obsBin)+' -M MultiDimFit SM_125_all_13TeV_xs_'+obsName+'_bin_v4.root -m 125.38 --freezeParameters MH -P r2e2muBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges r2e2muBin'+str(obsBin)+'=0.0,2.5 --redefineSignalPOI r2e2muBin'+str(obsBin)+' --algo=grid --points=100 --cminDefaultMinimizerStrategy 0 --saveInactivePOI=1'

                fidxs = 0
                fidxs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['ggH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['VBFH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['WH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['ZH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_2e2mu']*acc['ttH125_2e2mu_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys --setParameters r2e2muBin'+str(obsBin)+'='+str(round(fidxs,4))
                print cmd, '\n'
                output = processCmd(cmd)
                cmds.append(cmd)
                # Stat-only
                cmd = 'combine -n _'+obsName+'_r2e2muBin'+str(obsBin)+'_NoSys'
                if(not opt.UNBLIND): cmd = cmd + '_exp'
                cmd = cmd + ' -M MultiDimFit higgsCombine_'+obsName+'_r2e2muBin'+str(obsBin)+'.MultiDimFit.mH125.38'
                if(not opt.UNBLIND): cmd = cmd + '.123456'
                cmd = cmd + '.root -w w --snapshotName "MultiDimFit" -m 125.38 -P r2e2muBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges r2e2muBin'+str(obsBin)+'=0.0,2.5 --redefineSignalPOI r2e2muBin'+str(obsBin)+' --algo=grid --points=100 --cminDefaultMinimizerStrategy 0 --freezeNuisanceGroups nuis'
                if (opt.YEAR == 'Full'): cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
                else: cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_1'+str(opt.YEAR)+',CMS_fakeH_p3_1'+str(opt.YEAR)+',CMS_fakeH_p1_2'+str(opt.YEAR)+',CMS_fakeH_p3_2'+str(opt.YEAR)+',CMS_fakeH_p1_3'+str(opt.YEAR)+',CMS_fakeH_p3_3'+str(opt.YEAR)
                if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys --setParameters r2e2muBin'+str(obsBin)+'='+str(round(fidxs,4))
                print cmd+'\n'
                output = processCmd(cmd)
                cmds.append(cmd)

                # ----- 4e+4mu = 4l -----
                cmd = 'combine -n _'+obsName+'_r4lBin'+str(obsBin)+' -M MultiDimFit SM_125_all_13TeV_xs_'+obsName+'_bin_v4.root -m 125.38 --freezeParameters MH -P r4lBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges r4lBin'+str(obsBin)+'=0.0,2.5 --redefineSignalPOI r4lBin'+str(obsBin)+' --algo=grid --points=100 --cminDefaultMinimizerStrategy 0'

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
                if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys --setParameters r4lBin'+str(obsBin)+'='+str(round(fidxs,4))
                print cmd, '\n'
                output = processCmd(cmd)
                cmds.append(cmd)
                # Stat-only
                cmd = 'combine -n _'+obsName+'_r4lBin'+str(obsBin)+'_NoSys'
                if(not opt.UNBLIND): cmd = cmd + '_exp'
                cmd = cmd + ' -M MultiDimFit higgsCombine_'+obsName+'_r4lBin'+str(obsBin)+'.MultiDimFit.mH125.38'
                if(not opt.UNBLIND): cmd = cmd + '.123456'
                cmd = cmd + '.root -w w --snapshotName "MultiDimFit" -m 125.38 -P r4lBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges r4lBin'+str(obsBin)+'=0.0,2.5 --redefineSignalPOI r4lBin'+str(obsBin)+' --algo=grid --points=100 --cminDefaultMinimizerStrategy 0 --freezeNuisanceGroups nuis'
                if (opt.YEAR == 'Full'): cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
                else: cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_1'+str(opt.YEAR)+',CMS_fakeH_p3_1'+str(opt.YEAR)+',CMS_fakeH_p1_2'+str(opt.YEAR)+',CMS_fakeH_p3_2'+str(opt.YEAR)+',CMS_fakeH_p1_3'+str(opt.YEAR)+',CMS_fakeH_p3_3'+str(opt.YEAR)
                if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys --setParameters r4lBin'+str(obsBin)+'='+str(round(fidxs,4))
                print cmd+'\n'
                output = processCmd(cmd)
                cmds.append(cmd)


        elif physicalModel == 'v3':
            XH = []
            tmp_xs = {}
            tmp_xs_sm = {}
            for channel in ['4e','4mu','2e2mu']:
                for obsBin in range(nBins):
                    fidxs_sm = 0
                    fidxs_sm += higgs_xs['ggH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs_sm += higgs_xs['VBF_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs_sm += higgs_xs['WH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs_sm += higgs_xs['ZH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs_sm += higgs_xs['ttH_'+'125.0']*higgs4l_br['125.0'+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                    fidxs = 0

                    fidxs += higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    fidxs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                    # fidxs = fidxs_sm

                    tmp_xs_sm[channel+'_genbin'+str(obsBin)] = fidxs_sm
                    tmp_xs[channel+'_genbin'+str(obsBin)] = fidxs

            cmd_BR = ""
            for obsBin in range(nBins):
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

            for obsBin in range(nBins):
                XH.append(0.0)
                for channel in ['4e','4mu','2e2mu']:
                    XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH[obsBin]+=XH_fs

                _obsxsec = XH[obsBin]
                if obsName.startswith("mass4l"): max_range = '5.0'
                else: max_range = '2.5'
                ## The inclusive xsec for 2j phase space is about 2.49 fb, hence enlarge fit range
                if ('jj' in obsName) and (obsBin == 0): max_range = '5.0'
                if ('njet' in obsName) and ('pTj' in obsName) and (obsBin == 0): max_range = '5.0'
                cmd = 'combine -n _'+obsName+'_SigmaBin'+str(obsBin)+' -M MultiDimFit SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125.38 --freezeParameters MH -P SigmaBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges SigmaBin'+str(obsBin)+'=0.0,'+max_range+' --redefineSignalPOI SigmaBin'+str(obsBin)+' --algo=grid --points=100 --cminDefaultMinimizerStrategy 0'
                if(not opt.UNBLIND):
                    cmd = cmd + ' -t -1 --saveToys --setParameters SigmaBin'+str(obsBin)+'='+str(round(_obsxsec,4))
                    if opt.FIXFRAC:
                        cmd = cmd+','+cmd_BR
                if(opt.UNBLIND and opt.FIXFRAC):
                    cmd = cmd+' --setParameters '+cmd_BR
                print cmd, '\n'
                output = processCmd(cmd)
                cmds.append(cmd)
            # Stat-only
            XH = []
            for obsBin in range(nBins):
                XH.append(0.0)
                for channel in ['4e','4mu','2e2mu']:
                    XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH[obsBin]+=XH_fs
                _obsxsec = XH[obsBin]
                cmd = 'combine -n _'+obsName+'_SigmaBin'+str(obsBin)+'_NoSys'
                if(not opt.UNBLIND): cmd = cmd + '_exp'
                cmd = cmd + ' -M MultiDimFit higgsCombine_'+obsName+'_SigmaBin'+str(obsBin)+'.MultiDimFit.mH125.38'
                if(not opt.UNBLIND): cmd = cmd + '.123456'
                cmd = cmd + '.root -w w --snapshotName "MultiDimFit" -m 125.38 -P SigmaBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges SigmaBin0=0.0,'+max_range+' --redefineSignalPOI SigmaBin'+str(obsBin)+' --algo=grid --points=100 --cminDefaultMinimizerStrategy 0 --freezeNuisanceGroups nuis'
                if (opt.YEAR == 'Full'): cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
                else: cmd = cmd + ' --freezeParameters MH,CMS_fakeH_p1_1'+str(opt.YEAR)+',CMS_fakeH_p3_1'+str(opt.YEAR)+',CMS_fakeH_p1_2'+str(opt.YEAR)+',CMS_fakeH_p3_2'+str(opt.YEAR)+',CMS_fakeH_p1_3'+str(opt.YEAR)+',CMS_fakeH_p3_3'+str(opt.YEAR)
                if(not opt.UNBLIND):
                    cmd = cmd + ' -t -1 --saveToys --setParameters SigmaBin'+str(obsBin)+'='+str(round(_obsxsec,4))
                    if(opt.FIXFRAC):
                        cmd = cmd+','+cmd_BR
                if(opt.UNBLIND and opt.FIXFRAC):
                    cmd = cmd+' --setParameters '+cmd_BR
                print cmd+'\n'
                output = processCmd(cmd)
                cmds.append(cmd)

# ----------------- Main -----------------
_fit_dir = os.getcwd()
cmds = [] #List of all cmds
global doubleDiff
if 'vs' in opt.OBSNAME:
    obsName_tmp = opt.OBSNAME.split(' vs ')
    obsName = obsName_tmp[0]+'_'+obsName_tmp[1]
    doubleDiff = True
else:
    obsName = opt.OBSNAME
    doubleDiff = False
if (('pTj1' in obsName) | ('pTHj' in obsName) | ('mHj' in obsName) | ('pTj2' in obsName) | ('mjj' in obsName) | ('absdetajj' in obsName) | ('dphijj' in obsName) | ('pTHjj' in obsName)  | ('TCjmax' in obsName) | ('TBjmax' in obsName) | ('njets_pt30_eta4p7' in obsName)):
    JES = True
else:
    JES = False
runFiducialXS()
os.chdir(_fit_dir)
if (os.path.exists('commands_'+obsName+'.py')):
    os.system('rm commands_'+obsName+'.py')
with open('commands_'+obsName+'.py', 'w') as f:
    for i in cmds:
        f.write(str(i)+' \n')
        f.write('\n')
print "all modules successfully compiled"
sys.path.remove('../inputs/')
