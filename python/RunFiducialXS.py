import ROOT
import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json

sys.path.append('../inputs')
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
    parser.add_option('',   '--combineOnly',action='store_true', dest='impactsOnly',default=False, help='Run the measurement only, default is False')
    # Unblind option
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    # Calculate Systematic Uncertainties
    # parser.add_option('',   '--calcSys', action='store_true', dest='SYS', default=False, help='Calculate Systematic Uncertainties (in addition to stat+sys)')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # prepare the global flag if all the step should be run
    runAllSteps = not(opt.combineOnly or opt.impactsOnly)

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


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
def produceDatacards(obsName, observableBins, ModelName, PhysicalModel):
    print '\n'
    print '[Producing workspace/datacards for obsName '+obsName+', bins '+str(observableBins)+']'
    fStates = ['2e2mu','4mu','4e']
    nBins = len(observableBins)
    for year in years:
        os.chdir('../datacard/datacard_'+year)
        print 'Current diretory: datacard_'+year
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

    # prepare the set of bin boundaries to run over, only 1 bin in case of the inclusive measurement
    observableBins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
    ## Run for the given observable
    obsName = opt.OBSNAME
    print 'Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n'
    print 'Current directory: python'

    ## addConstrainedModel
    if(runAllSteps or opt.combineOnly):
        years_bis = years
        if(opt.YEAR == 'Full'):
            years_bis.append('Full')
        for year in years_bis:
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
            os.chdir('/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard/datacard_'+year)
            print 'Current directory: datacard_'+year
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
        os.chdir('/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard')
        print 'Current directory: datacard'
        if (opt.YEAR == 'Full'):
            cmd = 'combineCards.py datacard_2016/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt datacard_2017/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt datacard_2018/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt > hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
            cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_2016 CMS_hzz2e2mu_Zjets_2017 CMS_hzz2e2mu_Zjets_2018 CMS_hzz4e_Zjets_2016 CMS_hzz4e_Zjets_2017 CMS_hzz4e_Zjets_2018 CMS_hzz4mu_Zjets_2016 CMS_hzz4mu_Zjets_2017 CMS_hzz4mu_Zjets_2018 QCDscale_VV QCDscale_ggVV kfactor_ggzz lumi_13TeV_2016 lumi_13TeV_2017 lumi_13TeV_2018 norm_fakeH pdf_gg pdf_qqbar CMS_zz4l_sigma_e_sig_2017 CMS_zz4l_sigma_e_sig_2016 CMS_zz4l_sigma_m_sig_2018 CMS_zz4l_sigma_m_sig_2017 CMS_zz4l_sigma_m_sig_2016 CMS_zz4l_n_sig_3_2016 CMS_zz4l_n_sig_3_2017 CMS_zz4l_mean_e_sig_2016 CMS_zz4l_mean_e_sig_2017 CMS_zz4l_n_sig_3_2018 CMS_zz4l_mean_m_sig_2018 CMS_zz4l_mean_m_sig_2016 CMS_zz4l_mean_m_sig_2017 CMS_zz4l_sigma_e_sig_2018 CMS_zz4l_mean_e_sig_2018" >> hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
        else:
            cmd = 'cp datacard_'+str(opt.YEAR)+'/hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
            print cmd, '\n'
            processCmd(cmd,1)
            cmd = "sed -i 's|hzz4l|datacard_"+str(opt.YEAR)+"/hzz4l|g' hzz4l_all_13TeV_xs_"+obsName+"_bin_"+PhysicalModel+".txt" # Specify the right pattern to datacards (Before it was not necessary because there was a further combination)
            print cmd, '\n'
            processCmd(cmd,1)
            cmd = 'echo "nuis group = CMS_eff_e CMS_eff_m CMS_hzz2e2mu_Zjets_'+str(opt.YEAR)+' CMS_hzz4e_Zjets_'+str(opt.YEAR)+' CMS_hzz4mu_Zjets_'+str(opt.YEAR)+' QCDscale_VV QCDscale_ggVV kfactor_ggzz lumi_13TeV_'+str(opt.YEAR)+' norm_fakeH pdf_gg pdf_qqbar CMS_zz4l_sigma_e_sig_'+str(opt.YEAR)+' CMS_zz4l_sigma_m_sig_'+str(opt.YEAR)+' CMS_zz4l_n_sig_3_'+str(opt.YEAR)+' CMS_zz4l_mean_e_sig_'+str(opt.YEAR)+' CMS_zz4l_mean_m_sig_'+str(opt.YEAR)+' CMS_zz4l_sigma_e_sig_'+str(opt.YEAR)+'" >> hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt'
            processCmd(cmd,1)
            print cmd, '\n'

        # text-to-workspace
        if (PhysicalModel=="v3"):
            cmd = 'text2workspace.py hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.txt -P HiggsAnalysis.CombinedLimit.HZZ4L_Fiducial:differentialFiducialV3 --PO higgsMassRange=115,135 --PO nBin='+str(nBins-1)+' -o hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root'
            print cmd, '\n'
            processCmd(cmd)

        # The workspace got from text2workspace changes name from hzz4l_ to SM_125 and it is transferred to the combine_files directory
        cmd = 'cp hzz4l_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root ../combine_files/'+DataModelName+'_all_13TeV_xs_'+obsName+'_bin_'+PhysicalModel+'.root'
        print cmd, '\n'
        processCmd(cmd,1)
        os.chdir('..')

        os.chdir('/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/combine_files/')
        print 'Current directory: combine_files'
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
            cmd = cmd + '.root -w w --snapshotName "MultiDimFit" -m 125.0 --setParameters MH=125.0 -P SigmaBin'+str(obsBin)+' --floatOtherPOIs=1 --saveWorkspace --setParameterRanges MH=125.0,125.0:SigmaBin0=0.0,2.5 --redefineSignalPOI SigmaBin'+str(obsBin)+' --algo=grid --points=150 --freezeNuisanceGroups nuis'
            if (opt.YEAR == 'Full'): cmd = cmd + '--freezeParameters CMS_fakeH_p1_12018,CMS_fakeH_p3_12018,CMS_fakeH_p1_22018,CMS_fakeH_p3_22018,CMS_fakeH_p1_32018,CMS_fakeH_p3_32018,CMS_fakeH_p1_12017,CMS_fakeH_p3_12017,CMS_fakeH_p1_22017,CMS_fakeH_p3_22017,CMS_fakeH_p1_32017,CMS_fakeH_p3_32017,CMS_fakeH_p1_12016,CMS_fakeH_p3_12016,CMS_fakeH_p1_22016,CMS_fakeH_p3_22016,CMS_fakeH_p1_32016,CMS_fakeH_p3_32016'
            else: cmd = cmd + ' --freezeParameters CMS_fakeH_p1_1'+str(opt.YEAR)+',CMS_fakeH_p3_1'+str(opt.YEAR)+',CMS_fakeH_p1_2'+str(opt.YEAR)+',CMS_fakeH_p3_2'+str(opt.YEAR)+',CMS_fakeH_p1_3'+str(opt.YEAR)+',CMS_fakeH_p3_3'+str(opt.YEAR)
            if(not opt.UNBLIND): cmd = cmd + ' -t -1 --saveToys'
            print cmd+'\n'
            output = processCmd(cmd)

    # Impact plot
    if(runAllSteps or opt.impactsOnly):
        os.chdir('/afs/cern.ch/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/impacts/')
        print 'Current directory: impacts'
        nBins = len(observableBins)
        # First step (Files from asimov and data have the same name)
        cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125 --setParameters MH=125 --setParameterRanges MH=125,125 --doInitialFit --robustFit 1'
        if (not opt.UNBLIND): cmd = cmd + ' -t -1'
        print cmd, '\n'
        output = processCmd(cmd)
        # Second step (Files from asimov and data have the same name)
        cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125 --setParameters MH=125 --setParameterRanges MH=125,125 --doFits --parallel 4'
        if (not opt.UNBLIND): cmd = cmd + ' -t -1'
        print cmd, '\n'
        output = processCmd(cmd)
        for obsBin in range(nBins-1):
            # Third step
            cmd = 'combineTool.py -M Impacts -d ../combine_files/SM_125_all_13TeV_xs_'+obsName+'_bin_v3.root -m 125 --setParameters MH=125 --setParameterRanges MH=125,125:SigmaBin'+str(obsBin)+'=0.0,2.5 -o impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_'
            if (not opt.UNBLIND): cmd = cmd + 'asimov.json -t -1'
            elif (opt.UNBLIND): cmd = cmd + 'data.json -t -1'
            print cmd, '\n'
            output = processCmd(cmd)
            # plot
            cmd = 'plotImpacts.py -i impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_'
            if (not opt.UNBLIND): cmd = cmd + 'asimov.json -o impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_asimov --POI SigmaBin'+str(obsBin)
            elif (opt.UNBLIND): cmd = cmd + 'data.json -o impacts_v3_'+obsName+'_SigmaBin'+str(obsBin)+'_data --POI SigmaBin'+str(obsBin)
            print cmd, '\n'
            output = processCmd(cmd)

# ----------------- Main -----------------
runFiducialXS()
print "all modules successfully compiled"
