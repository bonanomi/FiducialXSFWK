import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)

### Define function for parsing options
def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-d', '--dir',    dest='SOURCEDIR',  type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--modelName',dest='MODELNAME',type='string',default='SM', help='Name of the Higgs production or spin-parity model, default is "SM", supported: "SM", "ggH", "VBF", "WH", "ZH", "ttH", "exotic","all"')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "mass4l", "pT4l", "massZ2", "rapidity4l", "cosThetaStar", "nets_reco_pt30_eta4p7"')
    parser.add_option('',   '--year',dest='YEAR',    type='string',default='',   help='year')
    parser.add_option('',   '--finalPlot',dest='FINALPLOT',    type='string',default='',   help='If making final plot, it changes sys.path to avoid conflicts')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('-f', '--doFit', action="store_true", dest='DOFIT', default=False, help='doFit, default false')
    parser.add_option('', '--doubleDiff', action="store_true", dest='DOUBLEDIFF', default=False, help='double-differential measurements, default false')
    parser.add_option('-p', '--doPlots', action="store_true", dest='DOPLOTS', default=False, help='doPlots, default false')
    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

sys.path.append('../inputs')
obsName = opt.OBSNAME

# observableBins = opt.OBSBINS
# observableBins = observableBins.split('|')
# observableBins.pop()
# observableBins.pop(0)

_temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['observableBins','acc','eff','err_eff','outinratio','err_outinratio','inc_wrongfrac','binfrac_wrongfrac','number_fake','lambdajesup','lambdajesdn'], -1)
observableBins = _temp.observableBins
if not opt.DOUBLEDIFF: lenObsBins = len(observableBins)-1
elif opt.DOUBLEDIFF: lenObsBins = len(observableBins)
acc = _temp.acc
eff = _temp.eff
err_eff = _temp.err_eff
outinratio = _temp.outinratio
err_outinratio = _temp.err_outinratio
inc_wrongfrac = _temp.inc_wrongfrac
binfrac_wrongfrac = _temp.binfrac_wrongfrac
number_fake = _temp.number_fake
lambdajesup = _temp.lambdajesup
lambdajesdn = _temp.lambdajesdn

_temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs4l_br','higgs_xs'], -1)
higgs4l_br = _temp.higgs4l_br
higgs_xs = _temp.higgs_xs
higgs4l_br['125.0_4l'] = higgs4l_br['125.0_2e2mu']+higgs4l_br['125.0_4e']+higgs4l_br['125.0_4mu']

# if (obsName=="mass4l"): fStates = ['4e','4mu','2e2mu','4l'] #AT Mi sembra che 4l non venga usato
if (obsName=="mass4l"): fStates = ['4e','4mu','2e2mu']
else: fStates = ['4e','4mu','2e2mu']

for fState in fStates:


    upfactors_diag = {}
    dnfactors_diag = {}

    for recobin in range(lenObsBins):

        fout_ggH = max(outinratio['ggH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)],0.0)
        # fout_VBF = max(outinratio['VBFH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)],0.0)
        # fout_WH = max(outinratio['WH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)],0.0)
        # fout_ZH = max(outinratio['ZH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)],0.0)
        # fout_ttH = max(outinratio['ttH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)],0.0)

        ggHxs_allgen=0.0
        VBFxs_allgen=0.0
        WHxs_allgen=0.0
        ZHxs_allgen=0.0
        ttHxs_allgen=0.0

        for genbin in range(lenObsBins):
            ggHxs_allgen += higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # VBFxs_allgen += higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # WHxs_allgen += higgs_xs['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # ZHxs_allgen += higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # ttHxs_allgen += higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]

        outin_SM = ggHxs_allgen*fout_ggH/(ggHxs_allgen)#+VBFxs_allgen+WHxs_allgen+ZHxs_allgen+ttHxs_allgen)
        # outin_SM += VBFxs_allgen*fout_VBF/(ggHxs_allgen+VBFxs_allgen+WHxs_allgen+ZHxs_allgen+ttHxs_allgen)
        # outin_SM += WHxs_allgen*fout_WH/(ggHxs_allgen+VBFxs_allgen+WHxs_allgen+ZHxs_allgen+ttHxs_allgen)
        # outin_SM += ZHxs_allgen*fout_ZH/(ggHxs_allgen+VBFxs_allgen+WHxs_allgen+ZHxs_allgen+ttHxs_allgen)
        # outin_SM += ttHxs_allgen*fout_ttH/(ggHxs_allgen+VBFxs_allgen+WHxs_allgen+ZHxs_allgen+ttHxs_allgen)

        for genbin in range(lenObsBins):

            ggHxs = higgs_xs['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            #ggHxs = acc['ggH_HRes_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]

            # VBFxs = higgs_xs['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # WHxs = higgs_xs['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # ZHxs = higgs_xs['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            # ttHxs = higgs_xs['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]


            effsm = ggHxs/(ggHxs)*max(eff['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)

            # effsm = ggHxs/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)*max(eff['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
            # #effsm = ggHxs/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)*max(eff['ggH_HRes_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
            # effsm += VBFxs/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)*max(eff['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
            # effsm += WHxs/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)*max(eff['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
            # effsm += ZHxs/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)*max(eff['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
            # effsm += ttHxs/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)*max(eff['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)

            err_effsm = err_eff['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]

            # err_effsm = sqrt(err_eff['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_eff['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_eff['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_eff['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2)# +
            #             #err_eff['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2)

            err_outinratiosm = err_outinratio['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]

            # err_outinratiosm = sqrt(err_outinratio['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_outinratio['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_outinratio['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_outinratio['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2 +
            #             err_outinratio['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)]**2)


            err_eff['SM_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = err_effsm
            err_outinratio['SM_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = err_outinratiosm
            eff['SM_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = effsm
            outinratio['SM_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = outin_SM

            effUp = effsm
            effDn = effsm
            outinUp = outin_SM
            outinDn = outin_SM

            upfactors = {}
            dnfactors = {}

            # Any pure model (closure test)
            #for f_ggH in [0.00001,0.5,1.0,2.0]:
            #    for f_VBF in [0.00001,0.8,1.0,1.52]:
            #        for f_WH in [0.00001,0.5,1.5,2.0]:
            #            for f_ZH in [0.00001,0.5,1.0,1.5,2.0]:
            #                for f_ttH in [0.00001,0.5,1.5,2.0]:

            #for f_ggH in [0.0,1.0]:
            #    for f_VBF in [0.0,1.0]:
            #        for f_WH in [0.0,1.0]:
            #            for f_ZH in [0.0,1.0]:
            #                for f_ttH in [0.0,1.0]:


            # 10% unc. 0 jet bin, 40% 3 jet bin
            #for f_ggH in [0.2,0.5,1.0,1.5,2.0]:
            #    for f_VBF in [0.0,0.5,1.0,1.5,2.0,5.0,10.0]:
            #        for f_WH in [0.0,0.5,1.0,1.5,2.0,5.0,10.0]:
            #            for f_ZH in [0.0,0.5,1.0,1.5,2.0,5.0,10.0]:
            #                for f_ttH in [0.0,0.5,1.0,1.5,2.0,5.0,10.0]:

            # 95 % C.L. from combination paper
            for f_ggH in [0.5,1.0,1.2,1.5]:
                for f_VBF in [0.4,1.0,1.5,2.0]:
                    for f_WH in [0.2,0.5,1.0,1.5,1.8]:
                        for f_ZH in [0.2,0.5,1.0,1.5,1.8]:
                            for f_ttH in [0.5,1.0,1.5,2.0,5.0]:


                                #if (obsName.startswith("njets")): f_ttH = 0.0
                                if ("jet" in obsName): f_ttH = 0.0
                                if (f_ggH+f_VBF+f_WH+f_ZH+f_ttH <0.00001): continue
                                tmp_eff = f_ggH*ggHxs/(f_ggH*ggHxs)*max(eff['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)

                                # tmp_eff = f_ggH*ggHxs/(f_ggH*ggHxs+f_VBF*VBFxs+f_WH*WHxs+f_ZH*ZHxs+f_ttH*ttHxs)*max(eff['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
                                # #tmp_eff = f_ggH*ggHxs/(f_ggH*ggHxs+f_VBF*VBFxs+f_WH*WHxs+f_ZH*ZHxs+f_ttH*ttHxs)*max(eff['ggH_HRes_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
                                # tmp_eff += f_VBF*VBFxs/(f_ggH*ggHxs+f_VBF*VBFxs+f_WH*WHxs+f_ZH*ZHxs+f_ttH*ttHxs)*max(eff['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
                                # tmp_eff += f_WH*WHxs/(f_ggH*ggHxs+f_VBF*VBFxs+f_WH*WHxs+f_ZH*ZHxs+f_ttH*ttHxs)*max(eff['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
                                # tmp_eff += f_ZH*ZHxs/(f_ggH*ggHxs+f_VBF*VBFxs+f_WH*WHxs+f_ZH*ZHxs+f_ttH*ttHxs)*max(eff['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)
                                # tmp_eff += f_ttH*ttHxs/(f_ggH*ggHxs+f_VBF*VBFxs+f_WH*WHxs+f_ZH*ZHxs+f_ttH*ttHxs)*max(eff['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)],0.0)

                                tmp_outin = f_ggH*ggHxs_allgen*fout_ggH/(f_ggH*ggHxs_allgen)    
                                # tmp_outin = f_ggH*ggHxs_allgen*fout_ggH/(f_ggH*ggHxs_allgen+f_VBF*VBFxs_allgen+f_WH*WHxs_allgen+f_ZH*ZHxs_allgen+f_ttH*ttHxs_allgen)
                                # tmp_outin += f_VBF*VBFxs_allgen*fout_VBF/(f_ggH*ggHxs_allgen+f_VBF*VBFxs_allgen+f_WH*WHxs_allgen+f_ZH*ZHxs_allgen+f_ttH*ttHxs_allgen)
                                # tmp_outin += f_WH*WHxs_allgen*fout_WH/(f_ggH*ggHxs_allgen+f_VBF*VBFxs_allgen+f_WH*WHxs_allgen+f_ZH*ZHxs_allgen+f_ttH*ttHxs_allgen)
                                # tmp_outin += f_ZH*ZHxs_allgen*fout_ZH/(f_ggH*ggHxs_allgen+f_VBF*VBFxs_allgen+f_WH*WHxs_allgen+f_ZH*ZHxs_allgen+f_ttH*ttHxs_allgen)
                                # tmp_outin += f_ttH*ttHxs_allgen*fout_ttH/(f_ggH*ggHxs_allgen+f_VBF*VBFxs_allgen+f_WH*WHxs_allgen+f_ZH*ZHxs_allgen+f_ttH*ttHxs_allgen)

                                if (tmp_eff*effUp>0.0):
                                    if ( (1-tmp_outin)/tmp_eff > (1-outinUp)/effUp ):
                                        effUp = tmp_eff
                                        outinUp = tmp_outin
                                        upfactors['ggH']=f_ggH
                                        # upfactors['VBF']=f_VBF
                                        # upfactors['WH']=f_WH
                                        # upfactors['ZH']=f_ZH
                                        # upfactors['ttH']=f_ttH
                                        if recobin==0  and genbin==0:
                                            upfactors_diag['ggH']=f_ggH
                                            # upfactors_diag['VBF']=f_VBF
                                            # upfactors_diag['WH']=f_WH
                                            # upfactors_diag['ZH']=f_ZH
                                            # upfactors_diag['ttH']=f_ttH

                                if (tmp_eff*effDn>0.0):
                                    if ( (1-tmp_outin)/tmp_eff < (1-outinDn)/effDn ):
                                        effDn = tmp_eff
                                        outinDn = tmp_outin
                                        dnfactors['ggH']=f_ggH
                                        # dnfactors['VBF']=f_VBF
                                        # dnfactors['WH']=f_WH
                                        # dnfactors['ZH']=f_ZH
                                        # dnfactors['ttH']=f_ttH
                                        if recobin==genbin:
                                            dnfactors_diag['ggH']=f_ggH
                                            # dnfactors_diag['VBF']=f_VBF
                                            # dnfactors_diag['WH']=f_WH
                                            # dnfactors_diag['ZH']=f_ZH
                                            # dnfactors_diag['ttH']=f_ttH




            eff['SMup_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = effUp
            eff['SMdn_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = effDn
            outinratio['SMup_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = outinUp
            outinratio['SMdn_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)] = outinDn


            #print fState,obsName,'genbin',str(genbin),'recobin',str(recobin),'effsm',effsm,'effUp',effUp,'effDn',effDn,'outinSM',outin_SM,'outinUp',outinUp,'outinDn',outinDn
            #if (effsm>0.0):
            #    print '(1-fout)/eff SM:',(1.0-outin_SM)/effsm,'(1-fout)/eff SMup',(1.0-outinUp)/effUp,'(1-fout)/eff SMdn',(1.0-outinDn)/effDn
            #print 'dnfactors',upfactors,'dnfactors',dnfactors


    # for recobin in range(lenObsBins):

    #     jesSM = ggHxs*(1.0+lambdajesup['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)])/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)
    #     #jesSM = ggHxs*(1.0+lambdajesup['ggH_HRes_125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)])/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)
    #     jesSM += VBFxs*(1.0+lambdajesup['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)])/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)
    #     jesSM += WHxs*(1.0+lambdajesup['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)])/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)
    #     jesSM += ZHxs*(1.0+lambdajesup['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)])/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)
    #     jesSM += ttHxs*(1.0+lambdajesup['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(recobin)])/(ggHxs+VBFxs+WHxs+ZHxs+ttHxs)
    #     jesSM = jesSM-1.0
    #     lambdajesup['SM_125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)] = jesSM


    #     jesSMup = upfactors_diag['ggH']*ggHxs*(1.0+lambdajesup['ggH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(upfactors_diag['ggH']*ggHxs+upfactors_diag['VBF']*VBFxs+upfactors_diag['WH']*WHxs+upfactors_diag['ZH']*ZHxs+upfactors_diag['ttH']*ttHxs)
    #     #jesSMup = upfactors_diag['ggH']*ggHxs*(1.0+lambdajesup['ggH_HRes_125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(upfactors_diag['ggH']*ggHxs+upfactors_diag['VBF']*VBFxs+upfactors_diag['WH']*WHxs+upfactors_diag['ZH']*ZHxs+upfactors_diag['ttH']*ttHxs)
    #     jesSMup += upfactors_diag['VBF']*VBFxs*(1.0+lambdajesup['VBFH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(upfactors_diag['ggH']*ggHxs+upfactors_diag['VBF']*VBFxs+upfactors_diag['WH']*WHxs+upfactors_diag['ZH']*ZHxs+upfactors_diag['ttH']*ttHxs)
    #     jesSMup += upfactors_diag['WH']*WHxs*(1.0+lambdajesup['WH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(upfactors_diag['ggH']*ggHxs+upfactors_diag['VBF']*VBFxs+upfactors_diag['WH']*WHxs+upfactors_diag['ZH']*ZHxs+upfactors_diag['ttH']*ttHxs)
    #     jesSMup += upfactors_diag['ZH']*ZHxs*(1.0+lambdajesup['ZH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(upfactors_diag['ggH']*ggHxs+upfactors_diag['VBF']*VBFxs+upfactors_diag['WH']*WHxs+upfactors_diag['ZH']*ZHxs+upfactors_diag['ttH']*ttHxs)
    #     jesSMup += upfactors_diag['ttH']*ttHxs*(1.0+lambdajesup['ttH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(upfactors_diag['ggH']*ggHxs+upfactors_diag['VBF']*VBFxs+upfactors_diag['WH']*WHxs+upfactors_diag['ZH']*ZHxs+upfactors_diag['ttH']*ttHxs)
    #     jesSMup = jesSMup-1.0
    #     lambdajesup['SMup_125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)] = jesSMup


    #     #print upfactors_diag
    #     #print dnfactors_diag

    #     #jesSMdn = dnfactors_diag['ggH']*ggHxs*(1.0+lambdajesup['ggH_HRes_125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(dnfactors_diag['ggH']*ggHxs+dnfactors_diag['VBF']*VBFxs+dnfactors_diag['WH']*WHxs+dnfactors_diag['ZH']*ZHxs+dnfactors_diag['ttH']*ttHxs)
    #     jesSMdn = dnfactors_diag['ggH']*ggHxs*(1.0+lambdajesup['ggH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(dnfactors_diag['ggH']*ggHxs+dnfactors_diag['VBF']*VBFxs+dnfactors_diag['WH']*WHxs+dnfactors_diag['ZH']*ZHxs+dnfactors_diag['ttH']*ttHxs)
    #     jesSMdn += dnfactors_diag['VBF']*VBFxs*(1.0+lambdajesup['VBFH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(dnfactors_diag['ggH']*ggHxs+dnfactors_diag['VBF']*VBFxs+dnfactors_diag['WH']*WHxs+dnfactors_diag['ZH']*ZHxs+dnfactors_diag['ttH']*ttHxs)
    #     jesSMdn += dnfactors_diag['WH']*WHxs*(1.0+lambdajesup['WH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(dnfactors_diag['ggH']*ggHxs+dnfactors_diag['VBF']*VBFxs+dnfactors_diag['WH']*WHxs+dnfactors_diag['ZH']*ZHxs+dnfactors_diag['ttH']*ttHxs)
    #     jesSMdn += dnfactors_diag['ZH']*ZHxs*(1.0+lambdajesup['ZH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(dnfactors_diag['ggH']*ggHxs+dnfactors_diag['VBF']*VBFxs+dnfactors_diag['WH']*WHxs+dnfactors_diag['ZH']*ZHxs+dnfactors_diag['ttH']*ttHxs)
    #     jesSMdn += dnfactors_diag['ttH']*ttHxs*(1.0+lambdajesup['ttH125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)])/(dnfactors_diag['ggH']*ggHxs+dnfactors_diag['VBF']*VBFxs+dnfactors_diag['WH']*WHxs+dnfactors_diag['ZH']*ZHxs+dnfactors_diag['ttH']*ttHxs)
    #     jesSMdn = jesSMdn-1.0
    #     lambdajesup['SMdn_125_'+fState+'_'+obsName+'_genbin'+str(recobin)+'_recobin'+str(recobin)] = jesSMdn


        # print fState,obsName,'recobin',str(recobin),'jesSM',jesSM,'jesSMup',jesSMup,'jesSMdn',jesSMdn


os.system('cp ../inputs/inputs_sig_'+opt.OBSNAME+'_'+opt.YEAR+'.py ../inputs/inputs_sig_'+opt.OBSNAME+'_'+opt.YEAR+'_ORIG.py')

with open('../inputs/inputs_sig_'+opt.OBSNAME+'_'+opt.YEAR+'.py', 'w') as f:
    f.write('observableBins = '+str(observableBins)+' \n')
    f.write('acc = '+str(acc)+' \n')
    f.write('eff = '+str(eff)+' \n')
    f.write('err_eff = '+str(err_eff)+' \n')
    f.write('outinratio = '+str(outinratio)+' \n')
    f.write('err_outinratio = '+str(err_outinratio)+' \n')
    f.write('inc_wrongfrac = '+str(inc_wrongfrac)+' \n')
    f.write('binfrac_wrongfrac = '+str(binfrac_wrongfrac)+' \n')
    f.write('number_fake = '+str(number_fake)+' \n')
    # f.write('lambdajesup = '+str(lambdajesup)+' \n')
    # f.write('lambdajesdn = '+str(lambdajesdn))
