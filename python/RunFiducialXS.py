import ROOT
import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json

sys.path.append('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/inputs')
from higgs_xsbr_13TeV import *
from createXSworkspace import createXSworkspace
from createDatacard import createDatacard

years = [2016, 2017, 2018]

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
    parser.add_option('',   '--fixMass',  dest='FIXMASS',  type='string',default='125.0',   help='Fix mass, default is a string "125.09" or can be changed to another string, e.g."125.6" or "False"')
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='fix the fractions of 4e and 4mu when extracting the results, default is False')
    # action options - "only"
    parser.add_option('',   '--effOnly',       action='store_true', dest='effOnly',       default=False, help='Extract the eff. factors only, default is False')
    parser.add_option('',   '--templatesOnly', action='store_true', dest='templatesOnly', default=False, help='Prepare the bkg shapes and fractions only, default is False')
    parser.add_option('',   '--uncertOnly',    action='store_true', dest='uncertOnly',    default=False, help='Extract the uncertanties only, default is False')
    parser.add_option('',   '--resultsOnly',   action='store_true', dest='resultsOnly',   default=False, help='Run the measurement only, default is False')
    parser.add_option('',   '--finalplotsOnly',action='store_true', dest='finalplotsOnly',default=False, help='Make the final plots only, default is False')
    # Unblind option
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    # Calculate Systematic Uncertainties
    parser.add_option('',   '--calcSys', action='store_true', dest='SYS', default=False, help='Calculate Systematic Uncertainties (in addition to stat+sys)')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # prepare the global flag if all the step should be run
    runAllSteps = not(opt.effOnly or opt.templatesOnly or opt.uncertOnly or opt.resultsOnly or opt.finalplotsOnly)

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


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


# parse the fit results from the MultiDim fit output "resultLog", for the bin and final state designated by "rTag"
def parseXSResults(resultLog, rTag):
    try:
        fXS_c = float(resultLog.split(rTag)[1].split(' (68%)')[0].strip().split(" ")[0])
        fXS_d = float('-'+resultLog.split(rTag)[1].split(' (68%)')[0].strip().split(" -")[1].split("/+")[0])
        fXS_u = float('+'+resultLog.split(rTag)[1].split(' (68%)')[0].strip().split(" -")[1].split("/+")[1])
        fXS = {'central':fXS_c, 'uncerDn':fXS_d, 'uncerUp':fXS_u}
        return fXS
    except IndexError:
        print "Parsing Failed!!! Inserting dummy values!!! check log!!!"
        fXS = {'central':-1.0, 'uncerDn':0.0, 'uncerUp':0.0}
        return fXS


### Create the asimov dataset and return fit results
def createAsimov(obsName, observableBins, ModelName, resultsXS, PhysicalModel):
    print '[Creating Asimov dataset for '+obsName+' using '+ModelName+']'

    # combination of bins (if there is just one bin, it is essentially a change of name from _bin0_ to _bin_)
    fStates = ['2e2mu','4mu','4e']
    nBins = len(observableBins)
    for year in years:
        os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard/datacard_'+str(year)+'/xs_125.0')
        for fState in fStates:
            if(nBins>1):
                cmd = 'combineCards.py '
                for obsBin in range(nBins-1):
                    cmd = cmd + 'hzz4l_'+fState+'S_13TeV_xs_'+obsName+'_bin'+str(obsBin)+'_'+PhysicalModel+'.txt '
                cmd = cmd + '> hzz4l_'+fState+'S_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
                print cmd, '\n'
                processCmd(cmd,1)
            else:
                print 'There is a problem during the combination over bins'

        # combine 3 final states
        cmd = 'combineCards.py hzz4l_4muS_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt hzz4l_4eS_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt hzz4l_2e2muS_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
        print cmd, '\n'
        processCmd(cmd,1)

    # Combine 3 years
    os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard')
    cmd = 'combineCards.py datacard_2016/xs_125.0/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt datacard_2017/xs_125.0/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt datacard_2018/xs_125.0/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
    print cmd, '\n'
    processCmd(cmd,1)

    cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_2016 CMS_hzz2e2mu_Zjets_2017 CMS_hzz2e2mu_Zjets_2018 CMS_hzz4e_Zjets_2016 CMS_hzz4e_Zjets_2017 CMS_hzz4e_Zjets_2018 CMS_hzz4mu_Zjets_2016 CMS_hzz4mu_Zjets_2017 CMS_hzz4mu_Zjets_2018 QCDscale_VV QCDscale_ggVV kfactor_ggzz lumi_13TeV_2016 lumi_13TeV_2017 lumi_13TeV_2018 norm_fakeH pdf_gg pdf_qqbar CMS_zz4l_sigma_e_sig_2017 CMS_zz4l_sigma_e_sig_2016 CMS_zz4l_sigma_m_sig_2018 CMS_zz4l_sigma_m_sig_2017 CMS_zz4l_sigma_m_sig_2016 CMS_zz4l_n_sig_3_2016 CMS_zz4l_n_sig_3_2017 CMS_zz4l_mean_e_sig_2016 CMS_zz4l_mean_e_sig_2017 CMS_zz4l_n_sig_3_2018 CMS_zz4l_mean_m_sig_2018 CMS_zz4l_mean_m_sig_2016 CMS_zz4l_mean_m_sig_2017 CMS_zz4l_sigma_e_sig_2018 CMS_zz4l_mean_e_sig_2018" >> hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
    processCmd(cmd,1)

    # text-to-workspace
    if (PhysicalModel=="v2"):
        cmd = 'text2workspace.py hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt -P HiggsAnalysis.CombinedLimit.HZZ4L_Fiducial:differentialFiducialV3 --PO higgsMassRange=115,135 --PO nBin='+str(nBins-1)+' -o hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root'
        print cmd, '\n'
        processCmd(cmd)

    # The workspace got from text2workspace changes name from hzz4l_ to SM_125 and it is transferred to the parent directory
    cmd = 'cp hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root ../'+ModelName+'_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root'
    print cmd, '\n'
    processCmd(cmd,1)
    os.chdir('..')

    # import acc factors from the Full Run2 coefficients file
    _temp = __import__('inputs_sig_'+obsName+'_Full', globals(), locals(), ['acc'], -1)
    acc = _temp.acc

    # Run the Combine
    if (PhysicalModel=="v2"):
        cmd =  'combine -n '+obsName+' -M MultiDimFit  '+ModelName+'_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root -m '+opt.ASIMOVMASS+' --setParameters '
        for fState in fStates:
            nBins = len(observableBins)
            for obsBin in range(nBins-1):
                fidxs = 0
                fidxs += higgs_xs['ggH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'ggH', fState, obsBin, fidxs
                fidxs += higgs_xs['VBF_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'VBF', fState, obsBin, fidxs
                fidxs += higgs_xs['WH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'WH', fState, obsBin, fidxs
                fidxs += higgs_xs['ZH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'ZH', fState, obsBin, fidxs
                fidxs += higgs_xs['ttH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'ttH', fState, obsBin, fidxs
                cmd = cmd + 'r'+fState+'Bin'+str(obsBin)+'='+str(fidxs)+','
        cmd =  cmd+ 'MH='+opt.ASIMOVMASS
        for fState in fStates:
            nBins = len(observableBins)
            for obsBin in range(nBins-1):
                cmd = cmd + ' -P r'+fState+'Bin'+str(obsBin)
        if (opt.FIXMASS=="False"):
            cmd = cmd + ' -P MH '
        else:
            cmd = cmd + ' --floatOtherPOIs=0'
        cmd = cmd +' -t -1 --saveWorkspace --saveToys'
        #cmd += ' --X-rtd TMCSO_PseudoAsimov=1000000'
        #cmd += ' --freezeNuisances r4muBin0,r4eBin0,r2e2muBin0'
        print cmd, '\n'
        output = processCmd(cmd)
        processCmd('mv higgsCombine'+obsName+'.MultiDimFit.mH'+opt.FIXMASS.rstrip('.0')+'.123456.root '+ModelName+'_all_'+obsName+'_13TeV_Asimov_'+PhysicalModel+'.root',1)
        #cmd = cmd.replace(' --freezeNuisances r4muBin0,r4eBin0,r2e2muBin0','')
        #cmd = cmd.replace(' --X-rtd TMCSO_PseudoAsimov=1000000','')
        cmd = cmd + ' --algo=singles --cl=0.68 --robustFit=1' # AT Cancellato RobustFit
        print cmd, '\n'
        output = processCmd(cmd)

    # parse the results for all the bins and the given final state
    tmp_resultsXS = {}
    for fState in fStates:
        for obsBin in range(len(observableBins)-1):
            binTag = str(obsBin)
            tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag] = parseXSResults(output,'r'+fState+'Bin'+str(obsBin)+' :')

    # merge the results for 3 final states, for the given bins
    for obsBin in range(len(observableBins)-1):
        binTag = str(obsBin)
        resultsXS['AsimovData_'+obsName+'_genbin'+binTag] = {'central':0.0, 'uncerDn':0.0, 'uncerUp':0.0}
        for fState in fStates:
            resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag] = {'central':0.0, 'uncerDn':0.0, 'uncerUp':0.0}
        tmp_central = 0.0
        tmp_uncerDn = 0.0
        tmp_uncerUp = 0.0
        for fState in fStates:
            tmp_central += tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag]['central']
            tmp_uncerDn += tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag]['uncerDn']**2
            tmp_uncerUp += tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag]['uncerUp']**2
            resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag]['central'] = float("{0:.5f}".format(tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag]['central']))
            resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag]['uncerDn'] = -float("{0:.5f}".format(tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag]['uncerDn']))
            resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag]['uncerUp'] = +float("{0:.5f}".format(tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag]['uncerUp']))
        resultsXS['AsimovData_'+obsName+'_genbin'+binTag]['central'] = float("{0:.5f}".format(tmp_central))
        resultsXS['AsimovData_'+obsName+'_genbin'+binTag]['uncerDn'] = -float("{0:.5f}".format(tmp_uncerDn**0.5))
        resultsXS['AsimovData_'+obsName+'_genbin'+binTag]['uncerUp'] = +float("{0:.5f}".format(tmp_uncerUp**0.5))
        # resultsXS['AsimovData_'+obsName+'_genbin'+binTag]['uncerDn'] = float("{0:.5f}".format(tmp_uncerDn**0.5))
        # resultsXS['AsimovData_'+obsName+'_genbin'+binTag]['uncerUp'] = float("{0:.5f}".format(tmp_uncerUp**0.5))

    # run combine with no systematics for stat only uncertainty
    if (opt.SYS):
        cmd = cmd + ' --freezeNuisanceGroups nuis --freezeParameters CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
        print cmd, '\n'
        output = processCmd(cmd)
        for fState in fStates:
            for obsBin in range(len(observableBins)-1):
                binTag = str(obsBin)
                tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly'] = parseXSResults(output,'r'+fState+'Bin'+str(obsBin)+' :')

        # merge the results for 3 final states, for the given bins
        for obsBin in range(len(observableBins)-1):
            binTag = str(obsBin)
            resultsXS['AsimovData_'+obsName+'_genbin'+binTag+'_statOnly'] = {'central':0.0, 'uncerDn':0.0, 'uncerUp':0.0}
            for fState in fStates:
                resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag+'_statOnly'] = {'central':0.0, 'uncerDn':0.0, 'uncerUp':0.0}
            tmp_central = 0.0
            tmp_uncerDn = 0.0
            tmp_uncerUp = 0.0
            for fState in fStates:
                tmp_central += tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly']['central']
                tmp_uncerDn += tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly']['uncerDn']**2
                tmp_uncerUp += tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly']['uncerUp']**2
                resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag+'_statOnly']['central'] = float("{0:.5f}".format(tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly']['central']))
                resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag+'_statOnly']['uncerDn'] = -float("{0:.5f}".format(tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly']['uncerDn']))
                resultsXS['AsimovData_'+obsName+'_'+fState+'_genbin'+binTag+'_statOnly']['uncerUp'] = +float("{0:.5f}".format(tmp_resultsXS[ModelName+'_'+fState+'_'+obsName+'_genbin'+binTag+'_statOnly']['uncerUp']))
            resultsXS['AsimovData_'+obsName+'_genbin'+binTag+'_statOnly']['central'] = float("{0:.5f}".format(tmp_central))
            resultsXS['AsimovData_'+obsName+'_genbin'+binTag+'_statOnly']['uncerDn'] = -float("{0:.5f}".format(tmp_uncerDn**0.5))
            resultsXS['AsimovData_'+obsName+'_genbin'+binTag+'_statOnly']['uncerUp'] = +float("{0:.5f}".format(tmp_uncerUp**0.5))

    return resultsXS


### Produce datacards for given obs and bin, for all final states
def produceDatacards(obsName, observableBins, ModelName, PhysicalModel):
    print '\n'
    print '[Producing workspace/datacards for obsName '+obsName+', bins '+str(observableBins)+']'
    fStates = ['2e2mu','4mu','4e']
    nBins = len(observableBins)
    for year in years:
        os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard/datacard_'+str(year))
        for fState in fStates:
            if (not obsName.startswith("mass4l")):
                for obsBin in range(nBins-1):
                    ndata = createXSworkspace(obsName,fState, nBins, obsBin, observableBins, False, True, ModelName, PhysicalModel, year)
                    createDatacard(obsName, fState, nBins, obsBin, observableBins, PhysicalModel, year, ndata)
            else:
                ndata = createXSworkspace(obsName,fState, nBins, 0, observableBins, False, True, ModelName, PhysicalModel, year)
                if obsName=='mass4l': os.system("cp xs_125.0_1bin/hzz4l_"+fState+"S_13TeV_xs_inclusive_bin0.txt xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")
                if obsName=='mass4lREFIT': os.system("cp xs_125.0_1bin/hzz4l_"+fState+"S_13TeV_xs_inclusiveREFIT_bin0.txt xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")
                os.system("sed -i 's~observation [0-9]*~observation "+str(ndata)+"~g' xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")
                os.system("sed -i 's~_xs.Databin0~_xs_"+ModelName+"_"+obsName+"_"+PhysicalModel+".Databin0~g' xs_125.0/hzz4l_"+fState+"S_13TeV_xs_"+obsName+"_bin0_"+PhysicalModel+".txt")



def runFiducialXS():

    # parse the arguments and options
    global opt, args, runAllSteps
    parseOptions()

    # prepare the set of bin boundaries to run over, only 1 bin in case of the inclusive measurement
    observableBins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
    ## Run for the given observable
    obsName = opt.OBSNAME
    print 'Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n'

    ## Create the Asimov dataset
    if(runAllSteps):
        for year in ['2016', '2017', '2018', 'Full']:
            if not os.path.exists('../inputs/inputs_sig_'+obsName+'_'+year+'_ORIG.py'):
                cmd = 'python addConstrainedModel.py -l -q -b --obsName="'+opt.OBSNAME+'" --obsBins="'+opt.OBSBINS+'" --year="'+year+'"'
                print cmd
                output = processCmd(cmd)
                print output
            elif os.path.exists('../inputs/inputs_sig_'+obsName+'_'+year+'_ORIG.py'):
                print 'addConstrainedModel '+year+' already done'

        DataModelName = 'SM_125'
        PhysicalModel = 'v3'
        produceDatacards(obsName, observableBins, DataModelName, PhysicalModel)


        # combination of bins (if there is just one bin, it is essentially a change of name from _bin0_ to _bin_)
        fStates = ['2e2mu','4mu','4e']
        nBins = len(observableBins)
        for year in years:
            os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard/datacard_'+str(year))
            for fState in fStates:
                if(nBins>1):
                    cmd = 'combineCards.py '
                    for obsBin in range(nBins-1):
                        cmd = cmd + 'hzz4l_'+fState+'S_13TeV_xs_'+obsName+'_bin'+str(obsBin)+'_'+PhysicalModel+'.txt '
                    cmd = cmd + '> hzz4l_'+fState+'S_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
                    print cmd, '\n'
                    processCmd(cmd,1)
                else:
                    print 'There is a problem during the combination over bins'

            # combine 3 final states
            cmd = 'combineCards.py hzz4l_4muS_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt hzz4l_4eS_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt hzz4l_2e2muS_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)

        # Combine 3 years
        os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard')
        cmd = 'combineCards.py datacard_2016/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt datacard_2017/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt datacard_2018/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
        print cmd, '\n'
        processCmd(cmd,1)

        cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_2016 CMS_hzz2e2mu_Zjets_2017 CMS_hzz2e2mu_Zjets_2018 CMS_hzz4e_Zjets_2016 CMS_hzz4e_Zjets_2017 CMS_hzz4e_Zjets_2018 CMS_hzz4mu_Zjets_2016 CMS_hzz4mu_Zjets_2017 CMS_hzz4mu_Zjets_2018 QCDscale_VV QCDscale_ggVV kfactor_ggzz lumi_13TeV_2016 lumi_13TeV_2017 lumi_13TeV_2018 norm_fakeH pdf_gg pdf_qqbar CMS_zz4l_sigma_e_sig_2017 CMS_zz4l_sigma_e_sig_2016 CMS_zz4l_sigma_m_sig_2018 CMS_zz4l_sigma_m_sig_2017 CMS_zz4l_sigma_m_sig_2016 CMS_zz4l_n_sig_3_2016 CMS_zz4l_n_sig_3_2017 CMS_zz4l_mean_e_sig_2016 CMS_zz4l_mean_e_sig_2017 CMS_zz4l_n_sig_3_2018 CMS_zz4l_mean_m_sig_2018 CMS_zz4l_mean_m_sig_2016 CMS_zz4l_mean_m_sig_2017 CMS_zz4l_sigma_e_sig_2018 CMS_zz4l_mean_e_sig_2018" >> hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
        processCmd(cmd,1)

        # text-to-workspace
        if (PhysicalModel=="v3"):
            cmd = 'text2workspace.py hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt -P HiggsAnalysis.CombinedLimit.HZZ4L_Fiducial:differentialFiducialV3 --PO higgsMassRange=115,135 --PO nBin='+str(nBins-1)+' -o hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root'
            print cmd, '\n'
            processCmd(cmd)

        # The workspace got from text2workspace changes name from hzz4l_ to SM_125 and it is transferred to the parent directory
        cmd = 'cp hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root ../'+DataModelName+'_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root'
        print cmd, '\n'
        processCmd(cmd,1)
        os.chdir('..')

        nBins = len(observableBins)
        for obsBin in range(nBins-1):
            cmd = 'combine -n _'+obsName+'_SigmaBin'+str(obsBin)+' -M MultiDimFit SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125.0 --setParameters MH=125.0 -P SigmaBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges MH=125.0,125.0:SigmaBin'+str(obsBin)+'=0.0,2.5 --redefineSignalPOI SigmaBin0 --algo=grid --points=150'
            if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys'
            print cmd, '\n'
            output = processCmd(cmd)
        # Stat-only
        for obsBin in range(nBins-1):
            cmd = 'combine -n _'+obsName+'_SigmaBin'+str(obsBin)+'_NoSys -M MultiDimFit higgsCombine_'+obsName+'_SigmaBin'+str(obsBin)+'.MultiDimFit.mH125'
            if(not opt.UNBLIND): cmd = cmd + '.123456'
            cmd = cmd + '.root -w w --snapshotName "MultiDimFit" -m 125.0 --setParameters MH=125.0 -P SigmaBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges MH=125.0,125.0:SigmaBin0=0.0,2.5 --redefineSignalPOI SigmaBin'+str(obsBin)+' --algo=grid --points=150 --freezeNuisanceGroups nuis --freezeParameters CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
            if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys'
            print cmd, '\n'
            output = processCmd(cmd)


        return

        cmd =  'combine -n '+obsName+' -M MultiDimFit  '+ModelName+'_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root -m '+opt.ASIMOVMASS+' --setParameters '
        for fState in fStates:
            nBins = len(observableBins)
            for obsBin in range(nBins-1):
                fidxs = 0
                fidxs += higgs_xs['ggH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'ggH', fState, obsBin, fidxs
                fidxs += higgs_xs['VBF_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'VBF', fState, obsBin, fidxs
                fidxs += higgs_xs['WH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'WH', fState, obsBin, fidxs
                fidxs += higgs_xs['ZH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'ZH', fState, obsBin, fidxs
                fidxs += higgs_xs['ttH_'+opt.ASIMOVMASS]*higgs4l_br[opt.ASIMOVMASS+'_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                # print 'ttH', fState, obsBin, fidxs
                cmd = cmd + 'r'+fState+'Bin'+str(obsBin)+'='+str(fidxs)+','
        cmd =  cmd+ 'MH='+opt.ASIMOVMASS
        for fState in fStates:
            nBins = len(observableBins)
            for obsBin in range(nBins-1):
                cmd = cmd + ' -P r'+fState+'Bin'+str(obsBin)
        if (opt.FIXMASS=="False"):
            cmd = cmd + ' -P MH '
        else:
            cmd = cmd + ' --floatOtherPOIs=0'
        cmd = cmd +' -t -1 --saveWorkspace --saveToys'
        #cmd += ' --X-rtd TMCSO_PseudoAsimov=1000000'
        #cmd += ' --freezeNuisances r4muBin0,r4eBin0,r2e2muBin0'
        print cmd, '\n'
        output = processCmd(cmd)
        processCmd('mv higgsCombine'+obsName+'.MultiDimFit.mH'+opt.FIXMASS.rstrip('.0')+'.123456.root '+ModelName+'_all_'+obsName+'_13TeV_Asimov_'+PhysicalModel+'.root',1)
        #cmd = cmd.replace(' --freezeNuisances r4muBin0,r4eBin0,r2e2muBin0','')
        #cmd = cmd.replace(' --X-rtd TMCSO_PseudoAsimov=1000000','')
        cmd = cmd + ' --algo=singles --cl=0.68 --robustFit=1' # AT Cancellato RobustFit
        print cmd, '\n'
        output = processCmd(cmd)




















    #     resultsXS = {}
    #     resultsXS = createAsimov(obsName, observableBins, asimovDataModelName, resultsXS, asimovPhysicalModel)
    #     print 'resultsXS: \n', resultsXS
    #     # plot the asimov predictions for data, signal, and backround in differential bins
    #     if (not obsName.startswith("mass4l")):
    #         cmd = 'python plotDifferentialBins.py -l -q -b --obsName="'+obsName+'" --obsBins="'+opt.OBSBINS+'" --asimovModel="'+asimovDataModelName+'"'
    #         if (opt.UNBLIND): cmd = cmd + ' --unblind'
    #         print cmd, '\n'
    #         output = processCmd(cmd)
    #         print output
    #
    #
    # ## Extract the results
    # # Use constrained SM
    # #ModelNames = ['SM_125','SMup_125','SMdn_125'] #AT
    # ModelNames = ['SM_125']
    # print "ModelNames", ModelNames
    #
    # if (obsName.startswith("mass4l")): PhysicalModels = ["v2","v3"]
    # else: PhysicalModels = ["v3"]
    #
    # if(runAllSteps or opt.resultsOnly):
    #     for PhysicalModel in PhysicalModels:
    #         for ModelName in ModelNames:
    #             produceDatacards(obsName, observableBins, ModelName, PhysicalModel)
    #             resultsXS = extractResults(obsName, observableBins, ModelName, PhysicalModel, asimovDataModelName, asimovPhysicalModel, resultsXS)
    #             print "resultsXS: \n", resultsXS
    #             # plot the fit results
    #             if (not obsName.startswith("mass4l")):
    #                 cmd = 'python plotAsimov_simultaneous.py -l -q -b --obsName="'+obsName+'" --obsBins="'+opt.OBSBINS+'" --asimovModel="'+asimovDataModelName+'" --unfoldModel="'+ModelName+'"'# +' --lumiscale=str(opt.LUMISCALE)'
    #                 if (opt.UNBLIND): cmd = cmd + ' --unblind'
    #                 print cmd, '\n'
    #                 output = processCmd(cmd)
    #                 print output
    #             # elif (PhysicalModel=="v2"):
    #             #     cmd = 'python plotAsimov_inclusive.py -l -q -b --obsName="'+obsName+'" --obsBins="'+opt.OBSBINS+'" --asimovModel="'+asimovDataModelName+'" --unfoldModel="'+ModelName+'"' #+' --lumiscale=str(opt.LUMISCALE)'
    #             #     if (opt.UNBLIND): cmd = cmd + ' --unblind'
    #             #     print cmd, '\n'
    #             #     output = processCmd(cmd)
    #             #     print output
    #
    #         # Calculate model dependance uncertainties
    #         modelIndependenceUncert = addModelIndependenceUncert(obsName, observableBins, resultsXS, asimovDataModelName,PhysicalModel)
    #         print "modelIndependenceUncert: \n", modelIndependenceUncert
    #         if (opt.FIXFRAC): floatfix = '_fixfrac'
    #         else: floatfix = ''
    #         with open('resultsXS_'+obsName+'_'+PhysicalModel+floatfix+'.py', 'w') as f:
    #             f.write('ModelNames = '+json.dumps(ModelNames)+';\n')
    #             f.write('asimovDataModelName = '+json.dumps(asimovDataModelName)+';\n')
    #             f.write('resultsXS = '+json.dumps(resultsXS)+';\n')
    #             f.write('modelIndUncert = '+json.dumps(modelIndependenceUncert))
    #
    #
    # # Make final differential plots
    # if(runAllSteps or finalplotsOnly):
    #     os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS')
    #     if (opt.UNBLIND):
    #         # In files doLScan*.sh change toy_asimov with data_obs
    #         cmd = "sed -i 's~-D toy_asimov~-D data_obs~g' doLScan*.sh"
    #         output = processCmd(cmd)
    #     else:
    #         cmd = "sed -i 's~-D data_obs~-D toy_asimov~g' doLScan*.sh"
    #         output = processCmd(cmd)
    #
    #     if (obsName=="mass4l"):
    #         cmd = './doLScan_mass4l.sh'
    #     if (obsName=="pT4l"):
    #         cmd = './doLScan_pT4l.sh'
    #     if (obsName=="njets_pt30_eta4p7"):
    #         cmd = './doLScan_njets.sh'
    #     if (obsName=="njets_pt30_eta2p5"):
    #         cmd = './doLScan_njets2p5.sh'
    #     if (obsName=="pt_leadingjet_pt30_eta4p7"):
    #         cmd = './doLScan_ptjet1.sh'
    #     if (obsName=="pt_leadingjet_pt30_eta2p5"):
    #         cmd = './doLScan_ptjet12p5.sh'
    #     if (obsName=="rapidity4l"):
    #         cmd = './doLScan_rapidity4l.sh'
    #     if (obsName=="cosThetaStar"):
    #         cmd = './doLScan_cosThetaStar.sh'
    #     if (obsName=="massZ2"):
    #         cmd = './doLScan_massZ2.sh'
    #     if (obsName=="massZ1"):
    #         cmd = './doLScan_massZ1.sh'
    #     output = processCmd(cmd)
    #
    #     if (opt.UNBLIND):
    #         cmd = 'python plotLHScans.py -l -q -b --obsName='+obsName+' --unblind'
    #         output = processCmd(cmd)
    #     elif (not opt.UNBLIND):
    #         cmd = 'python plotLHScans.py -l -q -b --obsName='+obsName
    #         output = processCmd(cmd)
    #
    #     for ModelName in ModelNames:
    #         if (not opt.FIXMASS=="False"):
    #             cmd = 'python producePlots.py -l -q -b --obsName="'+obsName+'" --obsBins="'+opt.OBSBINS+'" --unfoldModel="'+ModelName+'" --theoryMass="'+opt.FIXMASS+'"'
    #         else:
    #             cmd = 'python producePlots.py -l -q -b --obsName="'+obsName+'" --obsBins="'+opt.OBSBINS+'" --unfoldModel="'+ModelName+'" --theoryMass="125.0"'
    #         if (opt.FIXFRAC): cmd = cmd + ' --fixFrac'
    #         if (opt.UNBLIND): cmd = cmd + ' --unblind'
    #         print cmd, '\n'
    #         output = processCmd(cmd)
    #         print output
    #         cmd = cmd + ' --setLog'
    #         print cmd, '\n'
    #         output = processCmd(cmd)
    #         print output

# ----------------- Main -----------------
runFiducialXS()
print "all modules successfully compiled"
