import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *
#from sample_shortnames import *
def checkDir(folder_path):
    isdir = os.path.isdir(folder_path)
    if not isdir:
        print('Directory {} does not exist. Creating it.' .format(folder_path))
        os.mkdir(folder_path)

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
    parser.add_option('',   '--unfoldModel',dest='UNFOLD',type='string',default='SM_125', help='Name of the unfolding model for central value')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='',   help='Name of the observalbe, supported: "inclusive", "pT", "eta", "Njets"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='Use results from fixed fraction fit, default is False')
    parser.add_option('',   '--setLog', action='store_true', dest='SETLOG', default=False, help='set plot to log scale y, default is False')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--lumiscale', type='string', dest='LUMISCALE', default='1.0', help='Scale yields')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
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

if (not os.path.exists("plots")):
    os.system("mkdir plots")

from ROOT import *

from tdrStyle import *
setTDRStyle()

datamodel = opt.UNFOLD

sys.path.append('./inputs')
sys.path.append('./LHScans')

def plotXS(obsName, obs_bins, obs_bins_boundaries = False):

    _temp = __import__('inputs_sig_'+obsName+'_'+opt.YEAR, globals(), locals(), ['acc'], -1)
    acc = _temp.acc
    print 'inputs_sig_'+obsName+'_'+opt.YEAR
    # eff = _temp.eff
    # outinratio = _temp.outinratio
    if(opt.YEAR=='Full'):
        _temp = __import__('inputs_sig_'+obsName+'_NNLOPS_Full', globals(), locals(), ['acc'], -1)
        acc_NNLOPS = _temp.acc
    else:
        _temp = __import__('inputs_sig_'+obsName+'_NNLOPS_'+opt.YEAR, globals(), locals(), ['acc'], -1)
        acc_NNLOPS = _temp.acc

    if(acFlag):
        _temp = __import__('inputs_sig_ACggH_'+acSample+'_'+obsName+'_Full', globals(), locals(), ['acc'], -1)
        acc_AC = _temp.acc
        if(acFlagBis):
            _temp = __import__('inputs_sig_ACggH_'+acSampleBis+'_'+obsName+'_Full', globals(), locals(), ['acc'], -1)
            acc_ACbis = _temp.acc

    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs4l_br'], -1)
    higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br
    if (opt.FIXFRAC): floatfix = '_fixfrac'
    else: floatfix = ''
    if (obsName.startswith("mass4l")):
        if opt.UNBLIND: _temp = __import__('resultsXS_LHScan_observed_'+obsName+'_v3'+floatfix, globals(), locals(), ['resultsXS'], -1)
        else: _temp = __import__('resultsXS_LHScan_expected_'+obsName+'_v3'+floatfix, globals(), locals(), ['resultsXS'], -1)
        # modelNames = _temp.modelNames
        # asimovDataModelName = _temp.asimovDataModelName
        resultsXS = _temp.resultsXS
        # modelIndUncert = _temp.modelIndUncert
        if opt.UNBLIND: _temp = __import__('resultsXS_LHScan_observed_'+obsName+'_v2'+floatfix, globals(), locals(), ['resultsXS'], -1)
        else: _temp = __import__('resultsXS_LHScan_expected_'+obsName+'_v2'+floatfix, globals(), locals(), ['resultsXS'], -1)
        # modelNames_v2 = _temp.modelNames
        # asimovDataModelName_v2 = _temp.asimovDataModelName
        resultsXS_v2 = _temp.resultsXS
        # modelIndUncert_v2 = _temp.modelIndUncert
    else:
        # _temp = __import__('resultsXS_'+obsName+'_v3'+floatfix, globals(), locals(), ['modelNames', 'asimovDataModelName', 'resultsXS', 'modelIndUncert'], -1)
        # modelNames = _temp.modelNames
        # asimovDataModelName = _temp.asimovDataModelName
        # modelIndUncert = _temp.modelIndUncert
        if (opt.UNBLIND): _temp = __import__('resultsXS_LHScan_observed_'+obsName+'_v3'+floatfix, globals(), locals(), ['resultsXS'], -1)
        else: _temp = __import__('resultsXS_LHScan_expected_'+obsName+'_v3'+floatfix, globals(), locals(), ['resultsXS'], -1)
        resultsXS = _temp.resultsXS


    # Calculation of coefficients to scale ggH acceptances
    if(acFlag):
        accCorrFact = {}
        for pmode in ['VBFH125', 'WH125', 'ZH125', 'ttH125']:
            for fState in ['2e2mu', '4e', '4mu']:
                for obsBin in range(len(obs_bins)-1):
                    # if opt.TUNEAC:
                    #     accCorrFact[pmode+'_'+fState+'_genbin'+str(obsBin)] = acc[pmode+'_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0'] / acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(obsBin)+'_recobin0']
                    #     print pmode+'_'+fState+'_genbin'+str(obsBin), accCorrFact[pmode+'_'+fState+'_genbin'+str(obsBin)]
                    # else:
                    accCorrFact[pmode+'_'+fState+'_genbin'+str(obsBin)] = 1


    #print _temp
    #print resultsXS

    # acc_ggH_powheg = {}
    pdfunc_ggH_powheg = {}
    qcdunc_ggH_powheg = {}
    _temp = __import__('accUnc_'+obsName, globals(), locals(), ['pdfUncert','qcdUncert'], -1)
    # acc_ggH_powheg = _temp.acc #AT Non viene mai usato
    pdfunc_ggH_powheg = _temp.pdfUncert
    qcdunc_ggH_powheg = _temp.qcdUncert
    if not acFlag:
        _temp = __import__('accUnc_'+obsName+'_NNLOPS', globals(), locals(), ['pdfUncert','qcdUncert'], -1)
        pdfunc_ggH_nnlops = _temp.pdfUncert
        qcdunc_ggH_nnlops = _temp.qcdUncert
    else:
        _temp = __import__('accUnc_'+obsName+'_NNLOPS', globals(), locals(), ['pdfUncert','qcdUncert'], -1)
        pdfunc_ggH_nnlops = _temp.pdfUncert
        qcdunc_ggH_nnlops = _temp.qcdUncert

    # cross sections
    ggH_powheg = []
    ggH_powheg_unc_hi = []
    ggH_powheg_unc_lo = []
    ggH_minloHJ = []
    ggH_minloHJ_unc_hi = []
    ggH_minloHJ_unc_lo = []
    # ggH_HRes = []
    # ggH_HRes_unc_hi = []
    # ggH_HRes_unc_lo = []
    # NNLO theory unc
    ggH_powheg_NNLOunc_hi = []
    ggH_powheg_NNLOunc_lo = []
    ggH_minloHJ_NNLOunc_hi = []
    ggH_minloHJ_NNLOunc_lo = []
    # ggH_HRes_NNLOunc_hi = []
    # ggH_HRes_NNLOunc_lo = []
    # NLO theory unc
    ggH_powheg_NLOunc_hi = []
    ggH_powheg_NLOunc_lo = []
    ggH_minloHJ_NLOunc_hi = []
    ggH_minloHJ_NLOunc_lo = []
    # AC predictions
    if(acFlag):
        ggH_AC = []
        ggH_AC_unc_hi = []
        ggH_AC_unc_lo = []
        ggH_AC_NNLOunc_hi = []
        ggH_AC_NNLOunc_lo = []
        ggH_AC_NLOunc_hi = []
        ggH_AC_NLOunc_lo = []
        XH_AC = []
        XH_AC_unc = []
        if(acFlagBis):
            ggH_ACbis = []
            ggH_ACbis_unc_hi = []
            ggH_ACbis_unc_lo = []
            ggH_ACbis_NNLOunc_hi = []
            ggH_ACbis_NNLOunc_lo = []
            ggH_ACbis_NLOunc_hi = []
            ggH_ACbis_NLOunc_lo = []
            XH_ACbis = []
            XH_ACbis_unc = []
    # XH unc
    XH = []
    XH_unc = []
    # Data
    data = []
    data_hi = []
    data_lo = []
    data_hi2 = []
    data_lo2 = []
    asimovdata = []
    # Systematic unc.
    systematics_hi = []
    systematics_lo = []
    systematics_hi2 = []
    systematics_lo2 = []
    modeldep_hi = []
    modeldep_lo = []
    data_hi_allunc = []
    data_lo_allunc = []
    stat_hi = []
    stat_lo = []

    #process ggH qqH WH ZH ttH bkg_qqzz bkg_ggzz bkg_zjets
    #pdf_gg lnN 1.0720 - - - 1.0780 - 1.0710 -
    #pdf_qqbar lnN - 1.0270 1.0350 1.0350 - 1.0342 - -
    #pdf_hzz4l_accept lnN 1.02 1.02 1.02 1.02 1.02 - - -
    #QCDscale_ggH lnN 1.0750 - - - - - - -
    #QCDscale_qqH lnN - 1.0020 - - - - - -
    #QCDscale_VH lnN - - 1.0040 1.0155 - - - -
    #QCDscale_ttH lnN - - - - 1.0655 - - -
    #QCDscale_ggVV lnN - - - - - - 1.2435 -
    #BRhiggs_hzz4l lnN 1.02 1.02 1.02 1.02 1.02 - - -
    # unc_theory_ggH_hi = sqrt(0.072**2+0.075**2+0.02**2+0.02**2)
    # unc_theory_ggH_lo = sqrt(0.078**2+0.069**2+0.02**2+0.02**2)
    # unc_theory_XH_hi  = sqrt(0.027**2+0.02**2+0.002**2+0.02**2)
    # unc_theory_XH_lo  = unc_theory_XH_hi
    # unc_VBF = sqrt(0.027**2+0.02**2+0.002**2+0.02**2)
    # unc_WH = sqrt(0.035**2+0.02**2+0.004**2+0.02**2)
    # unc_ZH = sqrt(0.035**2+0.02**2+0.0155**2+0.02**2)
    # unc_ttH =  sqrt(0.078**2+0.02**2+0.0655**2+0.02**2)

    unc_acc = 0.02
    unc_br = 0.02

    #unc_pdf_ggH_hi = 0.075
    unc_pdf_ggH_hi = 0.032
    unc_pdf_ggH_lo = 0.032
    unc_pdf_VBF = 0.021
    unc_pdf_WH = 0.023
    unc_pdf_ZH = 0.025
    unc_pdf_ttH = 0.081

    #unc_qcd_ggH_hi = 0.072
    unc_qcd_ggH_hi = 0.039
    unc_qcd_ggH_lo = 0.039
    unc_qcd_VBF = 0.002
    unc_qcd_WH = 0.01
    unc_qcd_ZH = 0.031
    unc_qcd_ttH = 0.0655

    nBins=len(obs_bins)
    for obsBin in range(nBins-1):

        # theory cross sections
        ggH_powheg.append(0.0)
        ggH_powheg_unc_hi.append(0.0)
        ggH_powheg_unc_lo.append(0.0)
        ggH_minloHJ.append(0.0)
        ggH_minloHJ_unc_hi.append(0.0)
        ggH_minloHJ_unc_lo.append(0.0)
        # ggH_HRes.append(0.0)
        # ggH_HRes_unc_hi.append(0.0)
        # ggH_HRes_unc_lo.append(0.0)
        # NNLO theory unc
        ggH_powheg_NNLOunc_hi.append(0.0)
        ggH_powheg_NNLOunc_lo.append(0.0)
        ggH_minloHJ_NNLOunc_hi.append(0.0)
        ggH_minloHJ_NNLOunc_lo.append(0.0)
        # ggH_HRes_NNLOunc_hi.append(0.0)
        # ggH_HRes_NNLOunc_lo.append(0.0)
        # NLO theory unc
        ggH_powheg_NLOunc_hi.append(0.0)
        ggH_powheg_NLOunc_lo.append(0.0)
        ggH_minloHJ_NLOunc_hi.append(0.0)
        ggH_minloHJ_NLOunc_lo.append(0.0)
        # XH
        XH.append(0.0)
        XH_unc.append(0.0)
        # AC
        if(acFlag):
            ggH_AC.append(0.0)
            ggH_AC_unc_hi.append(0.0)
            ggH_AC_unc_lo.append(0.0)
            ggH_AC_NNLOunc_hi.append(0.0)
            ggH_AC_NNLOunc_lo.append(0.0)
            ggH_AC_NLOunc_hi.append(0.0)
            ggH_AC_NLOunc_lo.append(0.0)
            XH_AC.append(0.0)
            XH_AC_unc.append(0.0)
            if(acFlagBis):
                ggH_ACbis.append(0.0)
                ggH_ACbis_unc_hi.append(0.0)
                ggH_ACbis_unc_lo.append(0.0)
                ggH_ACbis_NNLOunc_hi.append(0.0)
                ggH_ACbis_NNLOunc_lo.append(0.0)
                ggH_ACbis_NLOunc_hi.append(0.0)
                ggH_ACbis_NLOunc_lo.append(0.0)
                XH_ACbis.append(0.0)
                XH_ACbis_unc.append(0.0)
        # Data
        data.append(0.0)
        data_hi.append(0.0)
        data_lo.append(0.0)
        data_hi2.append(0.0)
        data_lo2.append(0.0)
        asimovdata.append(0.0)
        # Systematic unc
        modeldep_hi.append(0.0)
        modeldep_lo.append(0.0)
        systematics_hi.append(0.0)
        systematics_lo.append(0.0)
        systematics_hi2.append(0.0)
        systematics_lo2.append(0.0)
        data_hi_allunc.append(0.0)
        data_lo_allunc.append(0.0)
        stat_hi.append(0.0)
        stat_lo.append(0.0)

        for channel in ['4e','4mu','2e2mu']:
            XH_fs = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

            XH[obsBin]+=XH_fs
            #XH_unc[obsBin]+= unc_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            #XH_unc[obsBin]+= unc_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            #XH_unc[obsBin]+= unc_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            #XH_unc[obsBin]+= unc_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

            # branching ratio uncertainty
            XH_unc_fs = (unc_br*XH_fs)**2
            # acceptance uncertainty
            XH_unc_fs += (unc_acc*XH_fs)**2

            # qcd scale
            XH_qcdunc_fs = (unc_qcd_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            XH_qcdunc_fs += (unc_qcd_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            XH_qcdunc_fs += (unc_qcd_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            XH_qcdunc_fs += (unc_qcd_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            XH_unc_fs += XH_qcdunc_fs

            # pdf
            XH_qqpdfunc_fs = (unc_pdf_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                              +unc_pdf_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                              +unc_pdf_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            XH_unc_fs += XH_qqpdfunc_fs

            # add pdf uncertainty for ttH to total XH uncertainty
            XH_unc_fs += (unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            # total XH uncertainty
            XH_unc[obsBin]+=sqrt(XH_unc_fs)

            # AC
            if(acFlag):
                XH_AC_fs = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['VBFH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                XH_AC_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['WH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                XH_AC_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ZH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                XH_AC_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ttH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                XH_AC[obsBin]+=XH_AC_fs

                # branching ratio uncertainty
                XH_AC_unc_fs = (unc_br*XH_AC_fs)**2
                # acceptance uncertainty
                XH_AC_unc_fs += (unc_acc*XH_AC_fs)**2

                XH_AC_qcdunc_fs = (unc_qcd_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['VBFH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                XH_AC_qcdunc_fs += (unc_qcd_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['WH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                XH_AC_qcdunc_fs += (unc_qcd_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ZH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                XH_AC_qcdunc_fs += (unc_qcd_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ttH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                XH_AC_unc_fs += XH_qcdunc_fs

                # pdf
                XH_AC_qqpdfunc_fs = (unc_pdf_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]**accCorrFact['VBFH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                  +unc_pdf_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['WH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                  +unc_pdf_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ZH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                XH_AC_unc_fs += XH_qqpdfunc_fs

                # add pdf uncertainty for ttH to total XH uncertainty
                XH_AC_unc_fs += (unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ttH125_'+fState+'_genbin'+str(obsBin)]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

                # total XH uncertainty
                XH_AC_unc[obsBin]+=sqrt(XH_AC_unc_fs)

                if(acFlagBis):
                    XH_ACbis_fs = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['VBFH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_ACbis_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['WH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_ACbis_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ZH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                    XH_ACbis_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ttH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                    XH_ACbis[obsBin]+=XH_ACbis_fs

                    # branching ratio uncertainty
                    XH_ACbis_unc_fs = (unc_br*XH_ACbis_fs)**2
                    # acceptance uncertainty
                    XH_ACbis_unc_fs += (unc_acc*XH_ACbis_fs)**2

                    XH_ACbis_qcdunc_fs = (unc_qcd_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['VBFH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    XH_ACbis_qcdunc_fs += (unc_qcd_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['WH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    XH_ACbis_qcdunc_fs += (unc_qcd_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ZH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    XH_ACbis_qcdunc_fs += (unc_qcd_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ttH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    XH_ACbis_unc_fs += XH_qcdunc_fs

                    # pdf
                    XH_ACbis_qqpdfunc_fs = (unc_pdf_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]**accCorrFact['VBFH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                      +unc_pdf_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['WH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                      +unc_pdf_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ZH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    XH_ACbis_unc_fs += XH_qqpdfunc_fs

                    # add pdf uncertainty for ttH to total XH uncertainty
                    XH_ACbis_unc_fs += (unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*accCorrFact['ttH125_'+fState+'_genbin'+str(obsBin)]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

                    # total XH uncertainty
                    XH_ACbis_unc[obsBin]+=sqrt(XH_ACbis_unc_fs)


            # ggH cross sections
            ggH_xsBR = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]
            #print "ggH_xsBR",ggH_xsBR
            #print "ggH_xsBR_emutau",higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_emutau']

            ggH_powheg[obsBin]+=ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            #ggH_minloHJ[obsBin]+=ggH_xsBR*acc_ggH_powheg['ggH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin'+str(obsBin)]
            ggH_minloHJ[obsBin]+=ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            # AC
            if(acFlag):
                ggH_AC[obsBin]+=ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                if(acFlagBis):
                    ggH_ACbis[obsBin]+=ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

            # for total uncertainty, correlate br and acc uncertainties across all channels (XH+ggH)
            total_NNLOunc_fs_powheg_hi =  (unc_br*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            total_NNLOunc_fs_powheg_lo =  (unc_br*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            total_NNLOunc_fs_powheg_hi +=  (unc_acc*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            total_NNLOunc_fs_powheg_lo +=  (unc_acc*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2

            total_NNLOunc_fs_minloHJ_hi =  (unc_br*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            total_NNLOunc_fs_minloHJ_lo =  (unc_br*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            total_NNLOunc_fs_minloHJ_hi +=  (unc_acc*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            total_NNLOunc_fs_minloHJ_lo +=  (unc_acc*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2

            # NLO and NNLO are the same at this point
            total_NLOunc_fs_powheg_hi = total_NNLOunc_fs_powheg_hi
            total_NLOunc_fs_powheg_lo = total_NNLOunc_fs_powheg_lo
            if acFlag:
                total_NLOunc_fs_AC_hi =  (unc_br*(XH_AC_fs+ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                total_NLOunc_fs_AC_lo =  (unc_br*(XH_AC_fs+ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                total_NLOunc_fs_AC_hi +=  (unc_acc*(XH_AC_fs+ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                total_NLOunc_fs_AC_lo +=  (unc_acc*(XH_AC_fs+ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                if acFlagBis:
                    total_NLOunc_fs_ACbis_hi =  (unc_br*(XH_ACbis_fs+ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                    total_NLOunc_fs_ACbis_lo =  (unc_br*(XH_ACbis_fs+ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                    total_NLOunc_fs_ACbis_hi +=  (unc_acc*(XH_ACbis_fs+ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
                    total_NLOunc_fs_ACbis_lo +=  (unc_acc*(XH_ACbis_fs+ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]))**2
            #total_NLOunc_fs_powheg_hi = 0.0
            #total_NLOunc_fs_powheg_lo = 0.0

            total_NLOunc_fs_minloHJ_hi = total_NNLOunc_fs_minloHJ_hi
            total_NLOunc_fs_minloHJ_lo = total_NNLOunc_fs_minloHJ_lo


            # add ggH qcd uncertainties (uncorrelated with anything else)
            #NNLO

            total_NNLOunc_fs_powheg_hi += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_powheg_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            total_NNLOunc_fs_minloHJ_hi += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qcdunc_fs
            if (obsName.startswith("mass4l")):
                total_NNLOunc_fs_minloHJ_hi += (unc_qcd_ggH_hi
                                                *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                total_NNLOunc_fs_minloHJ_lo += (unc_qcd_ggH_lo
                                                *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            else:
                total_NNLOunc_fs_minloHJ_hi += (qcdunc_ggH_nnlops["ggH125_NNLOPS_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                                *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                total_NNLOunc_fs_minloHJ_lo += (qcdunc_ggH_nnlops["ggH125_NNLOPS_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                                *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            #NLO
            total_NLOunc_fs_powheg_hi += XH_qcdunc_fs
            total_NLOunc_fs_powheg_lo += XH_qcdunc_fs
            total_NLOunc_fs_powheg_hi += (qcdunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NLOunc_fs_powheg_lo += (qcdunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            if acFlag:
                total_NLOunc_fs_AC_hi += XH_AC_qcdunc_fs
                total_NLOunc_fs_AC_hi += XH_AC_qcdunc_fs
                total_NLOunc_fs_AC_hi += (qcdunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                              *ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                total_NLOunc_fs_AC_lo += (qcdunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                              *ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

                if acFlagBis:
                    total_NLOunc_fs_ACbis_hi += XH_ACbis_qcdunc_fs
                    total_NLOunc_fs_ACbis_hi += XH_ACbis_qcdunc_fs
                    total_NLOunc_fs_ACbis_hi += (qcdunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                                  *ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    total_NLOunc_fs_ACbis_lo += (qcdunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                                  *ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2


            print channel,total_NLOunc_fs_powheg_hi

            total_NLOunc_fs_minloHJ_hi += XH_qcdunc_fs
            total_NLOunc_fs_minloHJ_lo += XH_qcdunc_fs
            total_NLOunc_fs_minloHJ_hi += (qcdunc_ggH_nnlops["ggH125_NNLOPS_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                           *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NLOunc_fs_minloHJ_lo += (qcdunc_ggH_nnlops["ggH125_NNLOPS_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                           *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            # add pdf unc, anti correlate ggH and ttH
            #NNLO
            if (obsName.startswith("mass4l")):
                obsUnc_pdf_ggH_hi = unc_pdf_ggH_hi
                obsUnc_pdf_ggH_lo = unc_pdf_ggH_lo
            else:
                obsUnc_pdf_ggH_hi = pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                obsUnc_pdf_ggH_lo = pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']

            total_NNLOunc_fs_powheg_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_hi += (obsUnc_pdf_ggH_hi*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_powheg_lo += (obsUnc_pdf_ggH_lo*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2


            total_NNLOunc_fs_minloHJ_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_hi += (obsUnc_pdf_ggH_hi*ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_minloHJ_lo += (obsUnc_pdf_ggH_lo*ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            #NLO
            total_NLOunc_fs_powheg_hi += XH_qqpdfunc_fs
            total_NLOunc_fs_powheg_lo += XH_qqpdfunc_fs
            total_NLOunc_fs_powheg_hi += (pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NLOunc_fs_powheg_lo += (pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            if acFlag:
                total_NLOunc_fs_AC_hi += XH_AC_qqpdfunc_fs
                total_NLOunc_fs_AC_lo += XH_AC_qqpdfunc_fs
                total_NLOunc_fs_AC_hi += (pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                              *ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                              -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                total_NLOunc_fs_AC_lo += (pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                              *ggH_xsBR*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                              -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_AC['ggH'+acSample+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                if acFlagBis:
                    total_NLOunc_fs_ACbis_hi += XH_ACbis_qqpdfunc_fs
                    total_NLOunc_fs_ACbis_lo += XH_ACbis_qqpdfunc_fs
                    total_NLOunc_fs_ACbis_hi += (pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                                  *ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                                  -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
                    total_NLOunc_fs_ACbis_lo += (pdfunc_ggH_powheg["ggH125_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                                  *ggH_xsBR*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                                  -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_ACbis['ggH'+acSampleBis+'_M125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2


            total_NLOunc_fs_minloHJ_hi += XH_qqpdfunc_fs
            total_NLOunc_fs_minloHJ_lo += XH_qqpdfunc_fs
            total_NLOunc_fs_minloHJ_hi += (pdfunc_ggH_nnlops["ggH125_NNLOPS_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerUp']
                                          *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NLOunc_fs_minloHJ_lo += (pdfunc_ggH_nnlops["ggH125_NNLOPS_"+channel+"_"+obsName.replace('_reco','_gen')+"_genbin"+str(obsBin)]['uncerDn']
                                          *ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                          -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            # finally total uncertainty (different final states are correlated)
            # NNLO
            ggH_powheg_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_powheg_hi)
            ggH_powheg_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_powheg_lo)
            ggH_minloHJ_NNLOunc_hi[obsBin]+=sqrt(total_NNLOunc_fs_minloHJ_hi)
            ggH_minloHJ_NNLOunc_lo[obsBin]+=sqrt(total_NNLOunc_fs_minloHJ_lo)
            # NLO
            ggH_powheg_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_powheg_hi)
            ggH_powheg_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_powheg_lo)
            ggH_minloHJ_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_minloHJ_hi)
            ggH_minloHJ_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_minloHJ_lo)
            if(acFlag):
                ggH_AC_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_AC_hi)
                ggH_AC_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_AC_lo)
                if(acFlagBis):
                    ggH_ACbis_NLOunc_hi[obsBin]+=sqrt(total_NLOunc_fs_ACbis_hi)
                    ggH_ACbis_NLOunc_lo[obsBin]+=sqrt(total_NLOunc_fs_ACbis_lo)

        ggH_powheg[obsBin]+=XH[obsBin]
        if(acFlag):
            ggH_AC[obsBin]+=XH_AC[obsBin]
            if(acFlagBis): ggH_ACbis[obsBin]+=XH_ACbis[obsBin]
        ggH_minloHJ[obsBin]+=XH[obsBin]

        if (opt.UNBLIND):
            data[obsBin] = resultsXS[datamodel+"_"+obsName+"_genbin"+str(obsBin)]["central"]
            data_hi[obsBin] = resultsXS[datamodel+"_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            data_lo[obsBin] = -1.0*resultsXS[datamodel+"_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]
            #data_hi[obsBin] = resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            #data_lo[obsBin] = -1.0*resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]
        else:
            data[obsBin] = resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["central"]
            data_hi[obsBin] = resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            data_lo[obsBin] = -1.0*resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]

        if (opt.UNBLIND):
            # modeldep_hi[obsBin] = modelIndUncert["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            # modeldep_lo[obsBin] = -1.0*modelIndUncert["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]
            systematics_hi[obsBin] = sqrt(max(0.0,resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            systematics_lo[obsBin] = sqrt(max(0.0,resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))
            #systematics_hi[obsBin] = sqrt(max(0.0,resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            #systematics_lo[obsBin] = sqrt(max(0.0,resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_LHScan["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))
        else:
            # modeldep_hi[obsBin] = modelIndUncert["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]
            # modeldep_lo[obsBin] = -1.0*modelIndUncert["AsimovData_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]
            systematics_hi[obsBin] = sqrt(max(0.0,resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2))
            systematics_lo[obsBin] = sqrt(max(0.0,resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2))
            stat_hi[obsBin] = resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]
            stat_lo[obsBin] = resultsXS["SM_125_"+obsName+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]

        data_hi_allunc[obsBin] = sqrt(data_hi[obsBin]**2+modeldep_hi[obsBin]**2)
        data_lo_allunc[obsBin] = sqrt(data_lo[obsBin]**2+modeldep_lo[obsBin]**2)

    if (obsName.startswith("mass4l")):
        for channel in ['2e2mu','4mu','4e']:

            # theory cross sections
            ggH_powheg.append(0.0)
            ggH_powheg_unc_hi.append(0.0)
            ggH_powheg_unc_lo.append(0.0)
            ggH_minloHJ.append(0.0)
            ggH_minloHJ_unc_hi.append(0.0)
            ggH_minloHJ_unc_lo.append(0.0)
            #ggH_HRes.append(0.0)
            #ggH_HRes_unc_hi.append(0.0)
            #ggH_HRes_unc_lo.append(0.0)
            # NNLO theory unc
            ggH_powheg_NNLOunc_hi.append(0.0)
            ggH_powheg_NNLOunc_lo.append(0.0)
            ggH_minloHJ_NNLOunc_hi.append(0.0)
            ggH_minloHJ_NNLOunc_lo.append(0.0)
            #ggH_HRes_NNLOunc_hi.append(0.0)
            #ggH_HRes_NNLOunc_lo.append(0.0)
            # NLO theory unc
            ggH_powheg_NLOunc_hi.append(0.0)
            ggH_powheg_NLOunc_lo.append(0.0)
            ggH_minloHJ_NLOunc_hi.append(0.0)
            ggH_minloHJ_NLOunc_lo.append(0.0)
            # XH
            XH.append(0.0)
            XH_unc.append(0.0)
            # Data
            data.append(0.0)
            data_hi.append(0.0)
            data_lo.append(0.0)
            # Systematic unc
            modeldep_hi.append(0.0)
            modeldep_lo.append(0.0)
            systematics_hi.append(0.0)
            systematics_lo.append(0.0)
            data_hi_allunc.append(0.0)
            data_lo_allunc.append(0.0)

            if (channel=='2e2mu'): bin = 1
            if (channel=='4mu'): bin = 2
            if (channel=='4e'): bin = 3

            obsBin=0

            print obsBin,acc
            XH_fs = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin0_recobin0']
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin0_recobin0']

            XH[bin]+=XH_fs

            #XH_unc[bin]+= unc_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBF_powheg_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            #XH_unc[bin]+= unc_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            #XH_unc[bin]+= unc_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']
            #XH_unc[bin]+= unc_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH_powheg_JHUgen_125_'+channel+'_'+obsName+'_genbin0_recobin0']

            # branching ratio uncertainty
            XH_unc_fs = (unc_br*XH_fs)**2
            # acceptance uncertainty
            XH_unc_fs += (unc_acc*XH_fs)**2

            # qcd scale
            XH_qcdunc_fs = (unc_qcd_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_qcdunc_fs += (unc_qcd_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_unc_fs += XH_qcdunc_fs

            # pdf
            XH_qqpdfunc_fs = (unc_pdf_VBF*higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                              +unc_pdf_WH*higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin0_recobin0']
                              +unc_pdf_ZH*higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2
            XH_unc_fs += XH_qqpdfunc_fs

            # add pdf uncertainty for ttH to total XH uncertainty
            XH_unc_fs += (unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin0_recobin0'])**2

            # total XH uncertainty
            XH_unc[bin]+=sqrt(XH_unc_fs)

            # ggH cross sections
            ggH_xsBR = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]

            ggH_powheg[bin]+=ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            ggH_minloHJ[bin]+=ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

             # for total uncertainty, correlate br and acc uncertainties across all channels (XH+ggH)
            total_NNLOunc_fs_powheg_hi =  (unc_br*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_powheg_lo =  (unc_br*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_powheg_hi +=  (unc_acc*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_powheg_lo +=  (unc_acc*(XH_fs+ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin0_recobin0']))**2

            total_NNLOunc_fs_minloHJ_hi = (unc_br*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_minloHJ_lo = (unc_br*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_minloHJ_hi += (unc_acc*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin0_recobin0']))**2
            total_NNLOunc_fs_minloHJ_lo += (unc_acc*(XH_fs+ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin0_recobin0']))**2

            # add ggH qcd uncertainties (uncorrelated with anything else)
            #NNLO
            total_NNLOunc_fs_powheg_hi += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qcdunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_powheg_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            total_NNLOunc_fs_minloHJ_hi += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qcdunc_fs
            total_NNLOunc_fs_minloHJ_hi += (unc_qcd_ggH_hi*ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_minloHJ_lo += (unc_qcd_ggH_lo*ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            # add pdf unc, anti correlate ggH and ttH
            #NNLO
            total_NNLOunc_fs_powheg_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_powheg_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_powheg_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            total_NNLOunc_fs_minloHJ_hi += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_lo += XH_qqpdfunc_fs
            total_NNLOunc_fs_minloHJ_hi += (unc_pdf_ggH_hi*ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2
            total_NNLOunc_fs_minloHJ_lo += (unc_pdf_ggH_lo*ggH_xsBR*acc_NNLOPS['ggH125_NNLOPS_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                                           -unc_pdf_ttH*higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])**2

            # finally total uncertainty (different final states are correlated)
            # NNLO
            ggH_powheg_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_powheg_hi)
            ggH_powheg_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_powheg_lo)
            ggH_minloHJ_NNLOunc_hi[bin]+=sqrt(total_NNLOunc_fs_minloHJ_hi)
            ggH_minloHJ_NNLOunc_lo[bin]+=sqrt(total_NNLOunc_fs_minloHJ_lo)

            ggH_powheg[bin]+=XH[bin]
            ggH_minloHJ[bin]+=XH[bin]

            data[bin] = resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["central"]
            data_hi[bin] = resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
            data_lo[bin] = -1.0*resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]

            if (opt.UNBLIND):
                # modeldep_hi[bin] = modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                # modeldep_lo[bin] = -1.0*modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                systematics_hi[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                systematics_lo[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)
            else:
                #modeldep_hi[bin] = modelIndUncert_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                #modeldep_lo[bin] = -1.0*modelIndUncert_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                #systematics_hi[bin] = sqrt(resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                #systematics_lo[bin] = sqrt(resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)
                if (obsName.startswith("mass4l")):
                    # modeldep_hi[bin] = modelIndUncert_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                    # modeldep_lo[bin] = -1.0*modelIndUncert_v2["AsimovData_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                    systematics_hi[bin] = sqrt(resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                    systematics_lo[bin] = sqrt(resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)
                else:
                    # modeldep_hi[bin] = modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]
                    # modeldep_lo[bin] = -1.0*modelIndUncert_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]
                    #systematics_hi[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                    #systematics_lo[bin] = sqrt(resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)
                    systematics_hi[bin] = sqrt(resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerUp"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerUp"]**2)
                    systematics_lo[bin] = sqrt(resultsXS_v2[datamodel+"_"+obsName+"_"+channel+"_genbin"+str(obsBin)]["uncerDn"]**2-resultsXS_v2["SM_125_"+obsName+"_"+channel+"_genbin"+str(obsBin)+'_statOnly']["uncerDn"]**2)

            data_hi_allunc[bin] = sqrt(data_hi[bin]**2+modeldep_hi[bin]**2)
            data_lo_allunc[bin] = sqrt(data_lo[bin]**2+modeldep_lo[bin]**2)



    print 'data',data
    sumdata = 0.0
    for i in range(len(data)):
        sumdata+=data[i]
    print obsName,'sum data',sumdata
    print 'data_hi',data_hi
    print 'data_lo',data_lo
    #print 'ggH HRes + XH',ggH_HRes
    #print 'NNLO ggH HRes + XH',ggH_HRes_NNLOunc_hi
    #print 'NNLO ggH HRes + XH',ggH_HRes_NNLOunc_lo
    print 'ggH_powheg',ggH_powheg
    if acFlag:
        print 'ggH_AC',ggH_AC
        if acFlagBis:
            print 'ggH_ACbis',ggH_ACbis
    print 'NLO ggH_powheg_hi',ggH_powheg_NLOunc_hi
    print 'NLO ggH_powheg_lo',ggH_powheg_NLOunc_lo
    print 'NNLO ggH_powheg_hi',ggH_powheg_NNLOunc_hi
    print 'NNLO ggH_powheg_lo',ggH_powheg_NNLOunc_lo
    print 'OLD ggH_powheg_hi',ggH_powheg_unc_hi
    print 'OLD ggH_powheg_lo',ggH_powheg_unc_lo
    print 'ggH_minloHJ',ggH_minloHJ
    print 'ggH_minloHJ_NNLOunc_hi',ggH_minloHJ_NNLOunc_hi
    print 'ggH_minloHJ_NNLOunc_lo',ggH_minloHJ_NNLOunc_lo
    print 'XH',XH
    print 'XH_unc',XH_unc
    if acFlag:
        print 'XH_AC',XH_AC
        print 'XH_AC_unc',XH_AC_unc
        print 'ggH_AC_NLOunc_hi', ggH_AC_NLOunc_hi
        print 'ggH_AC_NLOunc_lo', ggH_AC_NLOunc_lo
        if acFlagBis:
            print 'XH_ACbis',XH_ACbis
            print 'XH_ACbis_unc',XH_ACbis_unc
            print 'ggH_ACbis_NLOunc_hi', ggH_ACbis_NLOunc_hi
            print 'ggH_ACbis_NLOunc_lo', ggH_ACbis_NLOunc_lo
    print 'modedlep_hi',modeldep_hi
    print 'modeldep_lo',modeldep_lo
    print 'systematics_hi',systematics_hi
    print 'systematics_lo',systematics_lo
    print 'stat_hi',stat_hi
    print 'stat_lo',stat_lo
    if (obsName.startswith("mass4l")):
        offset = 0.0

        a_observable  = array('d',[0.5+i for i in range(0,4)])
        v_observable    = TVectorD(len(a_observable),a_observable)
        a_dobservable = array('d',[0.5 for i in range(0,4)])
        v_dobservable = TVectorD(len(a_dobservable),a_dobservable)

        a_zeros = array('d',[0.0 for i in range(0,4)])
        v_zeros = TVectorD(len(a_zeros),a_zeros)
        a_twos = array('d',[0.2*a_dobservable[i] for i in range(0,4)])
        v_twos = TVectorD(len(a_twos),a_twos)

        a_observable_1  = array('d',[(0.5*(float(obs_bins[i])+float(obs_bins[i+1]))+min(offset,0.25*(float(obs_bins[i+1])-float(obs_bins[i])))) for i in range(len(obs_bins)-1)])
        v_observable_1  = TVectorD(len(a_observable_1),a_observable_1)
        a_dobservable_1 = array('d',[0.125*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable_1 = TVectorD(len(a_dobservable_1),a_dobservable_1)

        a_observable_2  = array('d',[(0.5*(float(obs_bins[i])+float(obs_bins[i+1]))-min(offset,0.25*(float(obs_bins[i+1])-float(obs_bins[i])))) for i in range(len(obs_bins)-1)])
        v_observable_2  = TVectorD(len(a_observable_2),a_observable_2)
        a_dobservable_2 = array('d',[0.125*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable_2 = TVectorD(len(a_dobservable_2),a_dobservable_2)

        a_ggH_powheg = array('d',[ggH_powheg[i] for i in range(len(ggH_powheg))])
        v_ggH_powheg = TVectorD(len(a_ggH_powheg),a_ggH_powheg)
        a_ggH_powheg_unc_hi =  array('d',[ggH_powheg_NNLOunc_hi[i] for i in range(len(ggH_powheg))])
        a_ggH_powheg_unc_lo =  array('d',[ggH_powheg_NNLOunc_lo[i] for i in range(len(ggH_powheg))])
        v_ggH_powheg_unc_hi = TVectorD(len(a_ggH_powheg_unc_hi),a_ggH_powheg_unc_hi)
        v_ggH_powheg_unc_lo = TVectorD(len(a_ggH_powheg_unc_lo),a_ggH_powheg_unc_lo)

        print 'a_ggH_powheg_hi',a_ggH_powheg_unc_hi
        print 'a_ggH_powheg_lo',a_ggH_powheg_unc_lo

        a_ggH_minloHJ = array('d',[ggH_minloHJ[i] for i in range(len(ggH_minloHJ))])
        v_ggH_minloHJ = TVectorD(len(a_ggH_minloHJ),a_ggH_minloHJ)
        a_ggH_minloHJ_unc_hi =  array('d',[ggH_minloHJ_NNLOunc_hi[i] for i in range(len(ggH_minloHJ_unc_hi))])
        a_ggH_minloHJ_unc_lo =  array('d',[ggH_minloHJ_NNLOunc_lo[i] for i in range(len(ggH_minloHJ_unc_lo))])
        v_ggH_minloHJ_unc_hi = TVectorD(len(a_ggH_minloHJ_unc_hi),a_ggH_minloHJ_unc_hi)
        v_ggH_minloHJ_unc_lo = TVectorD(len(a_ggH_minloHJ_unc_lo),a_ggH_minloHJ_unc_lo)

        print 'a_ggH_minloHJ',a_ggH_minloHJ
        print 'a_ggH_minloHJ_hi',a_ggH_minloHJ_unc_hi
        print 'a_ggH_minloHJ_lo',a_ggH_minloHJ_unc_lo


        '''
        a_ggH_HRes = array('d',[ggH_HRes[i] for i in range(len(ggH_HRes))])
        v_ggH_HRes = TVectorD(len(a_ggH_HRes),a_ggH_HRes)
        a_ggH_HRes_unc_hi =  array('d',[ggH_HRes_NNLOunc_hi[i] for i in range(len(ggH_HRes_unc_hi))])
        a_ggH_HRes_unc_lo =  array('d',[ggH_HRes_NNLOunc_lo[i] for i in range(len(ggH_HRes_unc_lo))])
        v_ggH_HRes_unc_hi = TVectorD(len(a_ggH_HRes_unc_hi),a_ggH_HRes_unc_hi)
        v_ggH_HRes_unc_lo = TVectorD(len(a_ggH_HRes_unc_lo),a_ggH_HRes_unc_lo)
        print 'a_ggH_HRes',a_ggH_HRes
        print 'a_ggH_HRes_hi',a_ggH_HRes_unc_hi
        print 'a_ggH_HRes_lo',a_ggH_HRes_unc_lo
        '''

        a_XH = array('d',[XH[0],XH[1],XH[2],XH[3]])
        v_XH = TVectorD(len(a_XH),a_XH)

        a_XH_hi = array('d',[XH_unc[i] for i in range(len(a_XH))])
        v_XH_hi = TVectorD(len(a_XH_hi),a_XH_hi)

        a_XH_lo = array('d',[XH_unc[i] for i in range(len(a_XH))])
        v_XH_lo = TVectorD(len(a_XH_lo),a_XH_lo)

        a_data = array('d',[data[0],data[1],data[2],data[3]])
        v_data = TVectorD(len(a_data),a_data)
        a_data_hi = array('d',[data_hi[0],data_hi[1],data_hi[2],data_hi[3]])
        v_data_hi = TVectorD(len(a_data_hi),a_data_hi)
        a_data_lo = array('d',[data_lo[0],data_lo[1],data_lo[2],data_lo[3]])
        v_data_lo = TVectorD(len(a_data_lo),a_data_lo)

        a_systematics_hi = array('d',[systematics_hi[0],systematics_hi[1],systematics_hi[2],systematics_hi[3]])
        v_systematics_hi = TVectorD(len(a_systematics_hi),a_systematics_hi)
        a_systematics_lo = array('d',[systematics_lo[0],systematics_lo[1],systematics_lo[2],systematics_lo[3]])
        v_systematics_lo = TVectorD(len(a_systematics_lo),a_systematics_lo)

        a_modeldep_hi = array('d',[modeldep_hi[0],modeldep_hi[1],modeldep_hi[2],modeldep_hi[3]])
        v_modeldep_hi = TVectorD(len(a_modeldep_hi),a_modeldep_hi)
        a_modeldep_lo = array('d',[modeldep_lo[0],modeldep_lo[1],modeldep_lo[2],modeldep_lo[3]])
        v_modeldep_lo = TVectorD(len(a_modeldep_lo),a_modeldep_lo)

        v_data_hi_allunc = TVectorD(len(data_hi_allunc), array('d',[data_hi_allunc[i] for i in range(len(data_hi_allunc))]))
        v_data_lo_allunc = TVectorD(len(data_lo_allunc), array('d',[data_lo_allunc[i] for i in range(len(data_lo_allunc))]))

    elif doubleDiff:
        offset = 20.0
        a_observable  = array('d',[0.5*(float(obs_bins[i])+float(obs_bins[i+1])) for i in range(len(obs_bins)-1)])
        v_observable  = TVectorD(len(a_observable),a_observable)
        a_dobservable = array('d',[0.5*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable = TVectorD(len(a_dobservable),a_dobservable)

        a_observable_1  = array('d',[(0.5*(float(obs_bins[i])+float(obs_bins[i+1]))+min(offset,0.25*(float(obs_bins[i+1])-float(obs_bins[i])))) for i in range(len(obs_bins)-1)])
        v_observable_1  = TVectorD(len(a_observable_1),a_observable_1)
        a_dobservable_1 = array('d',[0.125*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable_1 = TVectorD(len(a_dobservable_1),a_dobservable_1)

        a_observable_2  = array('d',[(0.5*(float(obs_bins[i])+float(obs_bins[i+1]))-min(offset,0.25*(float(obs_bins[i+1])-float(obs_bins[i])))) for i in range(len(obs_bins)-1)])
        v_observable_2  = TVectorD(len(a_observable_2),a_observable_2)
        a_dobservable_2 = array('d',[0.125*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable_2 = TVectorD(len(a_dobservable_2),a_dobservable_2)

        a_zeros = array('d',[0.0 for i in range(len(obs_bins)-1)])
        v_zeros = TVectorD(len(a_zeros),a_zeros)
        a_twos = array('d',[0.015*(float(obs_bins[len(obs_bins)-1])-float(obs_bins[0])) for i in range(len(obs_bins)-1)])
        v_twos = TVectorD(len(a_twos),a_twos)

        a_ggH_powheg = array('d',[ggH_powheg[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_powheg))])
        v_ggH_powheg = TVectorD(len(a_ggH_powheg),a_ggH_powheg)
        a_ggH_powheg_unc_hi =  array('d',[ggH_powheg_NLOunc_hi[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_powheg_unc_hi))])
        a_ggH_powheg_unc_lo =  array('d',[ggH_powheg_NLOunc_lo[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_powheg_unc_lo))])
        v_ggH_powheg_unc_hi = TVectorD(len(a_ggH_powheg_unc_hi),a_ggH_powheg_unc_hi)
        v_ggH_powheg_unc_lo = TVectorD(len(a_ggH_powheg_unc_lo),a_ggH_powheg_unc_lo)

        a_ggH_minloHJ = array('d',[ggH_minloHJ[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_minloHJ))])
        v_ggH_minloHJ = TVectorD(len(a_ggH_minloHJ),a_ggH_minloHJ)
        a_ggH_minloHJ_unc_hi =  array('d',[ggH_minloHJ_NNLOunc_hi[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_minloHJ_unc_hi))])
        a_ggH_minloHJ_unc_lo =  array('d',[ggH_minloHJ_NNLOunc_lo[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_minloHJ_unc_lo))])
        v_ggH_minloHJ_unc_hi = TVectorD(len(a_ggH_minloHJ_unc_hi),a_ggH_minloHJ_unc_hi)
        v_ggH_minloHJ_unc_lo = TVectorD(len(a_ggH_minloHJ_unc_lo),a_ggH_minloHJ_unc_lo)

        '''
        a_ggH_HRes = array('d',[ggH_HRes[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_HRes))])
        v_ggH_HRes = TVectorD(len(a_ggH_HRes),a_ggH_HRes)
        a_ggH_HRes_unc_hi =  array('d',[ggH_HRes_NNLOunc_hi[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_HRes_unc_hi))])
        a_ggH_HRes_unc_lo =  array('d',[ggH_HRes_NNLOunc_lo[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(ggH_HRes_unc_lo))])
        v_ggH_HRes_unc_hi = TVectorD(len(a_ggH_HRes_unc_hi),a_ggH_HRes_unc_hi)
        v_ggH_HRes_unc_lo = TVectorD(len(a_ggH_HRes_unc_lo),a_ggH_HRes_unc_lo)
        '''

        a_XH = array('d',[XH[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(XH))])
        v_XH = TVectorD(len(a_XH),a_XH)

        a_XH_hi = array('d',[XH_unc[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(XH_unc))])
        a_XH_lo = array('d',[XH_unc[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(XH_unc))])
        v_XH_hi = TVectorD(len(a_XH_hi),a_XH_hi)
        v_XH_lo = TVectorD(len(a_XH_lo),a_XH_lo)

        a_data = array('d',[data[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(data))])
        v_data = TVectorD(len(a_data),a_data)
        a_data_hi = array('d',[data_hi[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(data_hi))])
        v_data_hi = TVectorD(len(a_data_hi),a_data_hi)
        a_data_lo = array('d',[data_lo[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(data_lo))])
        v_data_lo = TVectorD(len(a_data_lo),a_data_lo)

        a_data_hi2 = array('d',[data_hi2[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(data_hi2))])
        v_data_hi2 = TVectorD(len(a_data_hi2),a_data_hi2)
        a_data_lo2 = array('d',[data_lo2[i]/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(data_lo2))])
        v_data_lo2 = TVectorD(len(a_data_lo2),a_data_lo2)
        a_systematics_hi = array('d',[(systematics_hi[i])/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(systematics_hi))])
        v_systematics_hi = TVectorD(len(a_systematics_hi),a_systematics_hi)
        a_systematics_lo = array('d',[(systematics_lo[i])/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(systematics_lo))])
        v_systematics_lo = TVectorD(len(a_systematics_lo),a_systematics_lo)

        a_systematics_hi2 = array('d',[(systematics_hi2[i])/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(systematics_hi2))])
        v_systematics_hi2 = TVectorD(len(a_systematics_hi2),a_systematics_hi2)
        a_systematics_lo2 = array('d',[(systematics_lo2[i])/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(systematics_lo2))])
        v_systematics_lo2 = TVectorD(len(a_systematics_lo2),a_systematics_lo2)

        a_modeldep_hi = array('d',[(modeldep_hi[i])/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(modeldep_hi))])
        v_modeldep_hi = TVectorD(len(a_modeldep_hi),a_modeldep_hi)
        a_modeldep_lo = array('d',[(modeldep_lo[i])/((float(obs_bins_boundaries[i][1])-float(obs_bins_boundaries[i][0]))*(float(obs_bins_boundaries[i][3])-float(obs_bins_boundaries[i][2]))) for i in range(len(modeldep_lo))])
        v_modeldep_lo = TVectorD(len(a_modeldep_lo),a_modeldep_lo)

        v_data_hi_allunc = TVectorD(len(data_hi_allunc), array('d',[data_hi_allunc[i] for i in range(len(data_hi_allunc))]))
        v_data_lo_allunc = TVectorD(len(data_lo_allunc), array('d',[data_lo_allunc[i] for i in range(len(data_lo_allunc))]))

    else:
        if (obsName=="pT4l"): offset=20.0
        elif (obsName=="massZ2"): offset=20.0
        elif (obsName=="massZ1"): offset=20.0
        elif (obsName=="rapidity4l"): offset=20.0
        elif (obsName=="pTj1jet_pt30_eta2p5"): offset=30.0
        elif (obsName=="pTj1"): offset=30.0
        elif (obsName=="njets_pt30_eta2p5"): offset=999.0
        else: offset = 20.0
        a_observable  = array('d',[0.5*(float(obs_bins[i])+float(obs_bins[i+1])) for i in range(len(obs_bins)-1)])
        v_observable  = TVectorD(len(a_observable),a_observable)
        a_dobservable = array('d',[0.5*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable = TVectorD(len(a_dobservable),a_dobservable)

        a_observable_1  = array('d',[(0.5*(float(obs_bins[i])+float(obs_bins[i+1]))+min(offset,0.25*(float(obs_bins[i+1])-float(obs_bins[i])))) for i in range(len(obs_bins)-1)])
        v_observable_1  = TVectorD(len(a_observable_1),a_observable_1)
        a_dobservable_1 = array('d',[0.125*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable_1 = TVectorD(len(a_dobservable_1),a_dobservable_1)

        a_observable_2  = array('d',[(0.5*(float(obs_bins[i])+float(obs_bins[i+1]))-min(offset,0.25*(float(obs_bins[i+1])-float(obs_bins[i])))) for i in range(len(obs_bins)-1)])
        v_observable_2  = TVectorD(len(a_observable_2),a_observable_2)
        a_dobservable_2 = array('d',[0.125*(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(obs_bins)-1)])
        v_dobservable_2 = TVectorD(len(a_dobservable_2),a_dobservable_2)

        a_zeros = array('d',[0.0 for i in range(len(obs_bins)-1)])
        v_zeros = TVectorD(len(a_zeros),a_zeros)
        a_twos = array('d',[0.015*(float(obs_bins[len(obs_bins)-1])-float(obs_bins[0])) for i in range(len(obs_bins)-1)])
        v_twos = TVectorD(len(a_twos),a_twos)

        a_ggH_powheg = array('d',[ggH_powheg[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg))])
        v_ggH_powheg = TVectorD(len(a_ggH_powheg),a_ggH_powheg)
        a_ggH_powheg_unc_hi =  array('d',[ggH_powheg_NLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_hi))])
        a_ggH_powheg_unc_lo =  array('d',[ggH_powheg_NLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_powheg_unc_lo))])
        v_ggH_powheg_unc_hi = TVectorD(len(a_ggH_powheg_unc_hi),a_ggH_powheg_unc_hi)
        v_ggH_powheg_unc_lo = TVectorD(len(a_ggH_powheg_unc_lo),a_ggH_powheg_unc_lo)

        if acFlag:
            a_ggH_AC = array('d',[ggH_AC[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_AC))])
            v_ggH_AC = TVectorD(len(a_ggH_AC),a_ggH_AC)
            a_ggH_AC_unc_hi =  array('d',[ggH_AC_NLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_AC_unc_hi))])
            a_ggH_AC_unc_lo =  array('d',[ggH_AC_NLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_AC_unc_lo))])
            v_ggH_AC_unc_hi = TVectorD(len(a_ggH_AC_unc_hi),a_ggH_AC_unc_hi)
            v_ggH_AC_unc_lo = TVectorD(len(a_ggH_AC_unc_lo),a_ggH_AC_unc_lo)
            if acFlagBis:
                a_ggH_ACbis = array('d',[ggH_ACbis[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_ACbis))])
                v_ggH_ACbis = TVectorD(len(a_ggH_ACbis),a_ggH_ACbis)
                a_ggH_ACbis_unc_hi =  array('d',[ggH_ACbis_NLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_ACbis_unc_hi))])
                a_ggH_ACbis_unc_lo =  array('d',[ggH_ACbis_NLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_ACbis_unc_lo))])
                v_ggH_ACbis_unc_hi = TVectorD(len(a_ggH_ACbis_unc_hi),a_ggH_ACbis_unc_hi)
                v_ggH_ACbis_unc_lo = TVectorD(len(a_ggH_ACbis_unc_lo),a_ggH_ACbis_unc_lo)

        a_ggH_minloHJ = array('d',[ggH_minloHJ[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ))])
        v_ggH_minloHJ = TVectorD(len(a_ggH_minloHJ),a_ggH_minloHJ)
        a_ggH_minloHJ_unc_hi =  array('d',[ggH_minloHJ_NNLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ_unc_hi))])
        a_ggH_minloHJ_unc_lo =  array('d',[ggH_minloHJ_NNLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_minloHJ_unc_lo))])
        v_ggH_minloHJ_unc_hi = TVectorD(len(a_ggH_minloHJ_unc_hi),a_ggH_minloHJ_unc_hi)
        v_ggH_minloHJ_unc_lo = TVectorD(len(a_ggH_minloHJ_unc_lo),a_ggH_minloHJ_unc_lo)

        '''
        a_ggH_HRes = array('d',[ggH_HRes[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes))])
        v_ggH_HRes = TVectorD(len(a_ggH_HRes),a_ggH_HRes)
        a_ggH_HRes_unc_hi =  array('d',[ggH_HRes_NNLOunc_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes_unc_hi))])
        a_ggH_HRes_unc_lo =  array('d',[ggH_HRes_NNLOunc_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(ggH_HRes_unc_lo))])
        v_ggH_HRes_unc_hi = TVectorD(len(a_ggH_HRes_unc_hi),a_ggH_HRes_unc_hi)
        v_ggH_HRes_unc_lo = TVectorD(len(a_ggH_HRes_unc_lo),a_ggH_HRes_unc_lo)
        '''

        a_XH = array('d',[XH[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH))])
        v_XH = TVectorD(len(a_XH),a_XH)

        a_XH_hi = array('d',[XH_unc[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH_unc))])
        a_XH_lo = array('d',[XH_unc[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(XH_unc))])
        v_XH_hi = TVectorD(len(a_XH_hi),a_XH_hi)
        v_XH_lo = TVectorD(len(a_XH_lo),a_XH_lo)

        a_data = array('d',[data[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data))])
        v_data = TVectorD(len(a_data),a_data)
        a_data_hi = array('d',[data_hi[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data_hi))])
        v_data_hi = TVectorD(len(a_data_hi),a_data_hi)
        a_data_lo = array('d',[data_lo[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data_lo))])
        v_data_lo = TVectorD(len(a_data_lo),a_data_lo)

        a_data_hi2 = array('d',[data_hi2[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data_hi2))])
        v_data_hi2 = TVectorD(len(a_data_hi2),a_data_hi2)
        a_data_lo2 = array('d',[data_lo2[i]/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(data_lo2))])
        v_data_lo2 = TVectorD(len(a_data_lo2),a_data_lo2)
        a_systematics_hi = array('d',[(systematics_hi[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(systematics_hi))])
        v_systematics_hi = TVectorD(len(a_systematics_hi),a_systematics_hi)
        a_systematics_lo = array('d',[(systematics_lo[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(systematics_lo))])
        v_systematics_lo = TVectorD(len(a_systematics_lo),a_systematics_lo)

        a_systematics_hi2 = array('d',[(systematics_hi2[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(systematics_hi2))])
        v_systematics_hi2 = TVectorD(len(a_systematics_hi2),a_systematics_hi2)
        a_systematics_lo2 = array('d',[(systematics_lo2[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(systematics_lo2))])
        v_systematics_lo2 = TVectorD(len(a_systematics_lo2),a_systematics_lo2)

        a_modeldep_hi = array('d',[(modeldep_hi[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(modeldep_hi))])
        v_modeldep_hi = TVectorD(len(a_modeldep_hi),a_modeldep_hi)
        a_modeldep_lo = array('d',[(modeldep_lo[i])/(float(obs_bins[i+1])-float(obs_bins[i])) for i in range(len(modeldep_lo))])
        v_modeldep_lo = TVectorD(len(a_modeldep_lo),a_modeldep_lo)

        v_data_hi_allunc = TVectorD(len(data_hi_allunc), array('d',[data_hi_allunc[i] for i in range(len(data_hi_allunc))]))
        v_data_lo_allunc = TVectorD(len(data_lo_allunc), array('d',[data_lo_allunc[i] for i in range(len(data_lo_allunc))]))


    g_ggH_powheg = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
    g_ggH_powheg.SetFillStyle(3254);
    g_ggH_powheg.SetFillColor(ROOT.kAzure+2)
    g_ggH_powheg.SetLineColor(ROOT.kAzure+2)
    g_ggH_powheg.SetLineWidth(2)
    g_ggH_powheg.SetMarkerColor(ROOT.kAzure+2)

    g_ggH_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
    g_ggH_powhegBorder.SetFillStyle(0)
    g_ggH_powhegBorder.SetFillColor(ROOT.kAzure+2)
    g_ggH_powhegBorder.SetLineColor(ROOT.kAzure+2)
    g_ggH_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

    g_ggH_powhege0 = TGraphAsymmErrors(v_observable,v_ggH_powheg,v_dobservable,v_dobservable,v_zeros,v_zeros)
    g_ggH_powhege0.SetFillStyle(3254);
    g_ggH_powhege0.SetFillColor(ROOT.kAzure+2)
    g_ggH_powhege0.SetLineColor(ROOT.kAzure+2)
    g_ggH_powhege0.SetMarkerColor(ROOT.kAzure+2)

    if (obsName.startswith("mass4l")):

        g_ggH_powheg = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
        g_ggH_powheg.SetFillStyle(3254);
        g_ggH_powheg.SetFillColor(ROOT.kAzure+2)
        g_ggH_powheg.SetLineColor(ROOT.kAzure+2)
        g_ggH_powheg.SetLineWidth(2)
        g_ggH_powheg.SetMarkerColor(ROOT.kAzure+2)

        g_ggH_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
        g_ggH_powhegBorder.SetFillStyle(0)
        g_ggH_powhegBorder.SetFillColor(ROOT.kAzure+2)
        g_ggH_powhegBorder.SetLineColor(ROOT.kAzure+2)
        g_ggH_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        h_ggH_powheg = TH1D("h_ggH_powheg","h_ggH_powheg",4, 0, 4)
        for i in range(4):
            h_ggH_powheg.SetBinContent(i+1,v_ggH_powheg[i])
    elif doubleDiff:
        g_ggH_powheg = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
        g_ggH_powheg.SetFillStyle(3254);
        g_ggH_powheg.SetFillColor(ROOT.kAzure+2)
        g_ggH_powheg.SetLineColor(ROOT.kAzure+2)
        g_ggH_powheg.SetLineWidth(2)
        g_ggH_powheg.SetMarkerColor(ROOT.kAzure+2)

        g_ggH_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
        g_ggH_powhegBorder.SetFillStyle(0)
        g_ggH_powhegBorder.SetFillColor(ROOT.kAzure+2)
        g_ggH_powhegBorder.SetLineColor(ROOT.kAzure+2)
        g_ggH_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        h_ggH_powheg = TH1D("h_ggH_powheg","h_ggH_powheg",4, 0, 4)
        for i in range(4):
            h_ggH_powheg.SetBinContent(i+1,v_ggH_powheg[i])
    else:
        g_ggH_powheg = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
        g_ggH_powheg.SetFillStyle(3254);
        g_ggH_powheg.SetFillColor(ROOT.kAzure+2)
        g_ggH_powheg.SetLineColor(ROOT.kAzure+2)
        g_ggH_powheg.SetLineWidth(2)
        g_ggH_powheg.SetMarkerColor(ROOT.kAzure+2)

        g_ggH_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ggH_powheg,v_dobservable_1,v_dobservable_1,v_ggH_powheg_unc_lo,v_ggH_powheg_unc_hi)
        g_ggH_powhegBorder.SetFillStyle(0)
        g_ggH_powhegBorder.SetFillColor(ROOT.kAzure+2)
        g_ggH_powhegBorder.SetLineColor(ROOT.kAzure+2)
        g_ggH_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        if acFlag:
            g_ggH_AC = TGraphAsymmErrors(v_observable_2,v_ggH_AC,v_dobservable_2,v_dobservable_2,v_ggH_AC_unc_lo,v_ggH_AC_unc_hi)
            g_ggH_AC.SetFillStyle(3254);
            g_ggH_AC.SetFillColor(ROOT.kRed+2)
            g_ggH_AC.SetLineColor(ROOT.kRed+2)
            g_ggH_AC.SetLineWidth(2)
            g_ggH_AC.SetMarkerColor(ROOT.kRed+2)

            g_ggH_ACBorder = TGraphAsymmErrors(v_observable_2,v_ggH_AC,v_dobservable_2,v_dobservable_2,v_ggH_AC_unc_lo,v_ggH_AC_unc_hi)
            g_ggH_ACBorder.SetFillStyle(0)
            g_ggH_ACBorder.SetFillColor(ROOT.kRed+2)
            g_ggH_ACBorder.SetLineColor(ROOT.kRed+2)
            g_ggH_ACBorder.SetMarkerColor(ROOT.kRed+2)

            if acFlagBis:
                g_ggH_ACbis = TGraphAsymmErrors(v_observable_2,v_ggH_ACbis,v_dobservable_2,v_dobservable_2,v_ggH_ACbis_unc_lo,v_ggH_ACbis_unc_hi)
                g_ggH_ACbis.SetFillStyle(3254);
                g_ggH_ACbis.SetFillColor(ROOT.kViolet+2)
                g_ggH_ACbis.SetLineColor(ROOT.kViolet+2)
                g_ggH_ACbis.SetLineWidth(2)
                g_ggH_ACbis.SetMarkerColor(ROOT.kViolet+2)

                g_ggH_ACbisBorder = TGraphAsymmErrors(v_observable_2,v_ggH_ACbis,v_dobservable_2,v_dobservable_2,v_ggH_ACbis_unc_lo,v_ggH_ACbis_unc_hi)
                g_ggH_ACbisBorder.SetFillStyle(0)
                g_ggH_ACbisBorder.SetFillColor(ROOT.kViolet+2)
                g_ggH_ACbisBorder.SetLineColor(ROOT.kViolet+2)
                g_ggH_ACbisBorder.SetMarkerColor(ROOT.kViolet+2)


        if ("jet" in obsName and (not obsName.startswith("njets"))) :
            h_ggH_powheg = TH1D("h_ggH_powheg","h_ggH_powheg",nBins-2, array('d',[float(obs_bins[i]) for i in range(1,len(obs_bins))]) )
            for i in range(1,nBins-1):
                h_ggH_powheg.SetBinContent(i,v_ggH_powheg[i])
        else:
            h_ggH_powheg = TH1D("h_ggH_powheg","h_ggH_powheg",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            for i in range(nBins-1):
                h_ggH_powheg.SetBinContent(i+1,v_ggH_powheg[i])
            if acFlag:
                h_ggH_AC = TH1D("h_ggH_AC","h_ggH_AC",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
                for i in range(nBins-1):
                    h_ggH_AC.SetBinContent(i+1,v_ggH_AC[i])

                if acFlagBis:
                    h_ggH_ACbis = TH1D("h_ggH_ACbis","h_ggH_ACbis",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
                    for i in range(nBins-1):
                        h_ggH_ACbis.SetBinContent(i+1,v_ggH_ACbis[i])

    h_ggH_powheg.SetLineColor(ROOT.kAzure+2)
    h_ggH_powheg.SetLineWidth(2)

    if acFlag:
        h_ggH_AC.SetLineColor(ROOT.kRed+2)
        h_ggH_AC.SetLineWidth(2)
        if acFlagBis:
            h_ggH_ACbis.SetLineColor(ROOT.kViolet+2)
            h_ggH_ACbis.SetLineWidth(2)


    g_ggH_minloHJ = TGraphAsymmErrors(v_observable_2,v_ggH_minloHJ,v_dobservable_2,v_dobservable_2,v_ggH_minloHJ_unc_lo,v_ggH_minloHJ_unc_hi)
    g_ggH_minloHJ.SetFillStyle(3245);
    g_ggH_minloHJ.SetFillColor(ROOT.kOrange+2)
    g_ggH_minloHJ.SetLineColor(ROOT.kOrange+2)
    g_ggH_minloHJ.SetLineWidth(2)
    g_ggH_minloHJ.SetMarkerColor(ROOT.kOrange+2)

    g_ggH_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ggH_minloHJ,v_dobservable_2,v_dobservable_2,v_ggH_minloHJ_unc_lo,v_ggH_minloHJ_unc_hi)
    g_ggH_minloHJBorder.SetFillStyle(0)
    g_ggH_minloHJBorder.SetFillColor(ROOT.kOrange+2)
    g_ggH_minloHJBorder.SetLineColor(ROOT.kOrange+2)
    g_ggH_minloHJBorder.SetMarkerColor(ROOT.kOrange+2)

    g_ggH_minloHJe0 = TGraphAsymmErrors(v_observable,v_ggH_minloHJ,v_dobservable,v_dobservable,v_zeros,v_zeros)
    g_ggH_minloHJe0.SetFillStyle(3245);
    g_ggH_minloHJe0.SetFillColor(ROOT.kOrange+2)
    g_ggH_minloHJe0.SetLineColor(ROOT.kOrange+2)
    g_ggH_minloHJe0.SetMarkerColor(ROOT.kOrange+2)

    if (obsName.startswith("mass4l")):
        h_ggH_minloHJ = TH1D("h_ggH_minloHJ","h_ggH_minloHJ",4, 0, 4)
        for i in range(4):
            h_ggH_minloHJ.SetBinContent(i+1,v_ggH_minloHJ[i])
    elif doubleDiff:
        h_ggH_minloHJ = TH1D("h_ggH_minloHJ","h_ggH_minloHJ",nBins-1, 0, nBins-1)
        for i in range(nBins-1):
            h_ggH_minloHJ.SetBinContent(i+1,v_ggH_minloHJ[i])

    g_ggH_minloHJ = TGraphAsymmErrors(v_observable_2,v_ggH_minloHJ,v_dobservable_2,v_dobservable_2,v_ggH_minloHJ_unc_lo,v_ggH_minloHJ_unc_hi)
    g_ggH_minloHJ.SetFillStyle(3245);
    g_ggH_minloHJ.SetFillColor(ROOT.kOrange+2)
    g_ggH_minloHJ.SetLineColor(ROOT.kOrange+2)
    g_ggH_minloHJ.SetLineWidth(2)
    g_ggH_minloHJ.SetMarkerColor(ROOT.kOrange+2)

    g_ggH_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ggH_minloHJ,v_dobservable_2,v_dobservable_2,v_ggH_minloHJ_unc_lo,v_ggH_minloHJ_unc_hi)
    g_ggH_minloHJBorder.SetFillStyle(0)
    g_ggH_minloHJBorder.SetFillColor(ROOT.kOrange+2)
    g_ggH_minloHJBorder.SetLineColor(ROOT.kOrange+2)
    g_ggH_minloHJBorder.SetMarkerColor(ROOT.kOrange+2)

    if (not obsName.startswith("mass4l") and not doubleDiff):
        if ("jet" in obsName and (not obsName.startswith("njets"))):
            h_ggH_minloHJ = TH1D("h_ggH_minloHJ","h_ggH_minloHJ",nBins-2, array('d',[float(obs_bins[i]) for i in range(1,len(obs_bins))]) )
            for i in range(1,nBins-1):
                h_ggH_minloHJ.SetBinContent(i,v_ggH_minloHJ[i])
        else:
            h_ggH_minloHJ = TH1D("h_ggH_minloHJ","h_ggH_minloHJ",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            for i in range(nBins-1):
                h_ggH_minloHJ.SetBinContent(i+1,v_ggH_minloHJ[i])

    h_ggH_minloHJ.SetLineColor(ROOT.kOrange+2)
    h_ggH_minloHJ.SetLineWidth(2)



    '''
    g_ggH_HRes = TGraphAsymmErrors(v_observable,v_ggH_HRes,v_dobservable,v_dobservable,v_ggH_HRes_unc_lo,v_ggH_HRes_unc_hi)
    g_ggH_HRes.SetFillStyle(3254);
    g_ggH_HRes.SetFillColor(ROOT.kAzure+2)
    g_ggH_HRes.SetLineColor(ROOT.kAzure+2)
    g_ggH_HRes.SetLineWidth(2)
    g_ggH_HRes.SetMarkerColor(ROOT.kAzure+2)
    '''

    g_XH = TGraphAsymmErrors(v_observable,v_XH,v_dobservable,v_dobservable,v_XH_lo,v_XH_hi)
    g_XH.SetFillColor(ROOT.kGreen+3)
    g_XH.SetLineColor(ROOT.kGreen+3)

    g_data = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_data_lo,v_data_hi)
    g_data.SetMarkerColor(ROOT.kBlack)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.SetLineWidth(2)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1.4)

    g_data_e0 = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_zeros,v_zeros)
    g_data_e0.SetMarkerColor(ROOT.kBlack)
    g_data_e0.SetLineColor(ROOT.kBlack)
    g_data_e0.SetLineWidth(2)
    g_data_e0.SetMarkerStyle(20)
    g_data_e0.SetMarkerSize(1.4)

    # Ratio
    if not acFlag:
        v_ratio_data = TVectorD(len(data), array('d',[data[i]/ggH_minloHJ[i] for i in range(len(data))]))
        v_ratio_data_hi = TVectorD(len(data), array('d',[data_hi[i]/ggH_minloHJ[i] for i in range(len(data))]))
        v_ratio_data_lo = TVectorD(len(data), array('d',[data_lo[i]/ggH_minloHJ[i] for i in range(len(data))]))

        v_ratio_sys_hi = TVectorD(len(data), array('d',[systematics_hi[i]/ggH_minloHJ[i] for i in range(len(systematics_hi))]))
        v_ratio_sys_lo = TVectorD(len(data), array('d',[systematics_lo[i]/ggH_minloHJ[i] for i in range(len(systematics_lo))]))

        v_ratio_minloHJ = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ[i]/ggH_minloHJ[i] for i in range(len(ggH_minloHJ))]))
        v_ratio_minloHJ_hi = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ_NNLOunc_hi[i]/ggH_minloHJ[i] for i in range(len(ggH_minloHJ))]))
        v_ratio_minloHJ_lo = TVectorD(len(ggH_minloHJ), array('d',[ggH_minloHJ_NNLOunc_lo[i]/ggH_minloHJ[i] for i in range(len(ggH_minloHJ))]))

        v_ratio_powheg = TVectorD(len(ggH_powheg), array('d',[ggH_powheg[i]/ggH_minloHJ[i] for i in range(len(ggH_powheg))]))
        v_ratio_powheg_hi = TVectorD(len(ggH_powheg), array('d',[ggH_powheg_NLOunc_hi[i]/ggH_minloHJ[i] for i in range(len(ggH_powheg))]))
        v_ratio_powheg_lo = TVectorD(len(ggH_powheg), array('d',[ggH_powheg_NLOunc_lo[i]/ggH_minloHJ[i] for i in range(len(ggH_powheg))]))

        g_ratio_data = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_ratio_data_lo,v_ratio_data_hi)
        g_ratio_data.SetMarkerColor(ROOT.kBlack)
        g_ratio_data.SetLineColor(ROOT.kBlack)
        g_ratio_data.SetLineWidth(2)
        g_ratio_data.SetMarkerStyle(20)
        g_ratio_data.SetMarkerSize(1.4)

        g_ratio_datae0 = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_zeros,v_zeros)
        g_ratio_datae0.SetMarkerColor(ROOT.kBlack)
        g_ratio_datae0.SetLineColor(ROOT.kBlack)
        g_ratio_datae0.SetLineWidth(2)
        g_ratio_datae0.SetMarkerStyle(20)
        g_ratio_datae0.SetMarkerSize(1.4)

        g_ratio_sys = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_ratio_sys_lo,v_ratio_sys_hi)
        g_ratio_sys.SetMarkerColor(ROOT.kRed)
        g_ratio_sys.SetLineColor(ROOT.kRed)
        g_ratio_sys.SetFillColor(ROOT.kRed)
        g_ratio_sys.SetLineWidth(5)

        # g_ratio_powheg = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        # g_ratio_powheg.SetFillStyle(3254);
        # g_ratio_powheg.SetFillColor(ROOT.kAzure+2)
        # g_ratio_powheg.SetLineColor(ROOT.kAzure+2)
        # g_ratio_powheg.SetLineWidth(2)
        # g_ratio_powheg.SetMarkerColor(ROOT.kAzure+2)
        #
        # g_ratio_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        # g_ratio_powhegBorder.SetFillStyle(0);
        # g_ratio_powhegBorder.SetFillColor(ROOT.kAzure+2)
        # g_ratio_powhegBorder.SetLineColor(ROOT.kAzure+2)
        # g_ratio_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        g_ratio_powhege0 = TGraphAsymmErrors(v_observable,v_ratio_powheg,v_dobservable,v_dobservable,v_zeros,v_zeros)
        g_ratio_powhege0.SetFillStyle(3254);
        g_ratio_powhege0.SetFillColor(ROOT.kAzure+2)
        g_ratio_powhege0.SetLineColor(ROOT.kAzure+2)
        g_ratio_powhege0.SetLineWidth(2)
        g_ratio_powhege0.SetMarkerColor(ROOT.kAzure+2)

    else:
        v_ratio_data = TVectorD(len(data), array('d',[data[i]/ggH_powheg[i] for i in range(len(data))]))
        v_ratio_data_hi = TVectorD(len(data), array('d',[data_hi[i]/ggH_powheg[i] for i in range(len(data))]))
        v_ratio_data_lo = TVectorD(len(data), array('d',[data_lo[i]/ggH_powheg[i] for i in range(len(data))]))

        v_ratio_sys_hi = TVectorD(len(data), array('d',[systematics_hi[i]/ggH_powheg[i] for i in range(len(systematics_hi))]))
        v_ratio_sys_lo = TVectorD(len(data), array('d',[systematics_lo[i]/ggH_powheg[i] for i in range(len(systematics_lo))]))

        v_ratio_AC = TVectorD(len(ggH_AC), array('d',[ggH_AC[i]/ggH_powheg[i] for i in range(len(ggH_AC))]))
        v_ratio_AC_hi = TVectorD(len(ggH_AC), array('d',[ggH_AC_NLOunc_hi[i]/ggH_powheg[i] for i in range(len(ggH_AC))]))
        v_ratio_AC_lo = TVectorD(len(ggH_AC), array('d',[ggH_AC_NLOunc_lo[i]/ggH_powheg[i] for i in range(len(ggH_AC))]))
        if acFlagBis:
            v_ratio_ACbis = TVectorD(len(ggH_ACbis), array('d',[ggH_ACbis[i]/ggH_powheg[i] for i in range(len(ggH_ACbis))]))
            v_ratio_ACbis_hi = TVectorD(len(ggH_ACbis), array('d',[ggH_ACbis_NLOunc_hi[i]/ggH_powheg[i] for i in range(len(ggH_ACbis))]))
            v_ratio_ACbis_lo = TVectorD(len(ggH_ACbis), array('d',[ggH_ACbis_NLOunc_lo[i]/ggH_powheg[i] for i in range(len(ggH_ACbis))]))

        v_ratio_powheg = TVectorD(len(ggH_powheg), array('d',[ggH_powheg[i]/ggH_powheg[i] for i in range(len(ggH_powheg))]))
        v_ratio_powheg_hi = TVectorD(len(ggH_powheg), array('d',[ggH_powheg_NLOunc_hi[i]/ggH_powheg[i] for i in range(len(ggH_powheg))]))
        v_ratio_powheg_lo = TVectorD(len(ggH_powheg), array('d',[ggH_powheg_NLOunc_lo[i]/ggH_powheg[i] for i in range(len(ggH_powheg))]))

        g_ratio_data = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_ratio_data_lo,v_ratio_data_hi)
        g_ratio_data.SetMarkerColor(ROOT.kBlack)
        g_ratio_data.SetLineColor(ROOT.kBlack)
        g_ratio_data.SetLineWidth(2)
        g_ratio_data.SetMarkerStyle(20)
        g_ratio_data.SetMarkerSize(1.4)

        g_ratio_datae0 = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_zeros,v_zeros)
        g_ratio_datae0.SetMarkerColor(ROOT.kBlack)
        g_ratio_datae0.SetLineColor(ROOT.kBlack)
        g_ratio_datae0.SetLineWidth(2)
        g_ratio_datae0.SetMarkerStyle(20)
        g_ratio_datae0.SetMarkerSize(1.4)

        g_ratio_sys = TGraphAsymmErrors(v_observable,v_ratio_data,v_zeros,v_zeros,v_ratio_sys_lo,v_ratio_sys_hi)
        g_ratio_sys.SetMarkerColor(ROOT.kRed+2)
        g_ratio_sys.SetLineColor(ROOT.kRed+2)
        g_ratio_sys.SetFillColor(ROOT.kRed+2)
        g_ratio_sys.SetLineWidth(5)

        # g_ratio_powheg = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        # g_ratio_powheg.SetFillStyle(3254);
        # g_ratio_powheg.SetFillColor(ROOT.kAzure+2)
        # g_ratio_powheg.SetLineColor(ROOT.kAzure+2)
        # g_ratio_powheg.SetLineWidth(2)
        # g_ratio_powheg.SetMarkerColor(ROOT.kAzure+2)
        #
        # g_ratio_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        # g_ratio_powhegBorder.SetFillStyle(0);
        # g_ratio_powhegBorder.SetFillColor(ROOT.kAzure+2)
        # g_ratio_powhegBorder.SetLineColor(ROOT.kAzure+2)
        # g_ratio_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        g_ratio_powhege0 = TGraphAsymmErrors(v_observable,v_ratio_powheg,v_dobservable,v_dobservable,v_zeros,v_zeros)
        g_ratio_powhege0.SetFillStyle(3254);
        g_ratio_powhege0.SetFillColor(ROOT.kAzure+2)
        g_ratio_powhege0.SetLineColor(ROOT.kAzure+2)
        g_ratio_powhege0.SetLineWidth(2)
        g_ratio_powhege0.SetMarkerColor(ROOT.kAzure+2)

    if (obsName.startswith("mass4l")):
        g_ratio_powheg = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        g_ratio_powheg.SetFillStyle(3254);
        g_ratio_powheg.SetFillColor(ROOT.kAzure+2)
        g_ratio_powheg.SetLineColor(ROOT.kAzure+2)
        g_ratio_powheg.SetLineWidth(2)
        g_ratio_powheg.SetMarkerColor(ROOT.kAzure+2)

        g_ratio_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        g_ratio_powhegBorder.SetFillStyle(0);
        g_ratio_powhegBorder.SetFillColor(ROOT.kAzure+2)
        g_ratio_powhegBorder.SetLineColor(ROOT.kAzure+2)
        g_ratio_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        h_ratio_powheg = TH1D("h_ratio_powheg","h_ratio_powheg",4, 0, 4)
        for i in range(4):
            h_ratio_powheg.SetBinContent(i+1,v_ratio_powheg[i])
    elif doubleDiff:
        g_ratio_powheg = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        g_ratio_powheg.SetFillStyle(3254);
        g_ratio_powheg.SetFillColor(ROOT.kAzure+2)
        g_ratio_powheg.SetLineColor(ROOT.kAzure+2)
        g_ratio_powheg.SetLineWidth(2)
        g_ratio_powheg.SetMarkerColor(ROOT.kAzure+2)

        g_ratio_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        g_ratio_powhegBorder.SetFillStyle(0);
        g_ratio_powhegBorder.SetFillColor(ROOT.kAzure+2)
        g_ratio_powhegBorder.SetLineColor(ROOT.kAzure+2)
        g_ratio_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        h_ratio_powheg = TH1D("h_ratio_powheg","h_ratio_powheg",4, 0, 4)
        for i in range(4):
            h_ratio_powheg.SetBinContent(i+1,v_ratio_powheg[i])
    else:
        g_ratio_powheg = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        g_ratio_powheg.SetFillStyle(3254);
        g_ratio_powheg.SetFillColor(ROOT.kAzure+2)
        g_ratio_powheg.SetLineColor(ROOT.kAzure+2)
        g_ratio_powheg.SetLineWidth(2)
        g_ratio_powheg.SetMarkerColor(ROOT.kAzure+2)

        g_ratio_powhegBorder = TGraphAsymmErrors(v_observable_1,v_ratio_powheg,v_dobservable_1,v_dobservable_1,v_ratio_powheg_lo,v_ratio_powheg_hi)
        g_ratio_powhegBorder.SetFillStyle(0);
        g_ratio_powhegBorder.SetFillColor(ROOT.kAzure+2)
        g_ratio_powhegBorder.SetLineColor(ROOT.kAzure+2)
        g_ratio_powhegBorder.SetMarkerColor(ROOT.kAzure+2)

        if ("jet" in obsName and (not obsName.startswith("njets"))):
            h_ratio_powheg = TH1D("h_ratio_powheg","h_ratio_powheg",nBins-2, array('d',[float(obs_bins[i]) for i in range(1,len(obs_bins))]) )
            for i in range(1,nBins-1):
                h_ratio_powheg.SetBinContent(i,v_ratio_powheg[i])
        else:
            h_ratio_powheg = TH1D("h_ratio_powheg","h_ratio_powheg",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            for i in range(nBins-1):
                h_ratio_powheg.SetBinContent(i+1,v_ratio_powheg[i])

    h_ratio_powheg.SetLineColor(ROOT.kAzure+2)
    h_ratio_powheg.SetLineWidth(2)

    # g_ratio_minloHJ = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
    # g_ratio_minloHJ.SetFillStyle(3245);
    # g_ratio_minloHJ.SetFillColor(ROOT.kOrange+2)
    # g_ratio_minloHJ.SetLineColor(ROOT.kOrange+2)
    # g_ratio_minloHJ.SetLineWidth(2)
    # g_ratio_minloHJ.SetMarkerColor(ROOT.kOrange+2)
    #
    # g_ratio_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
    # g_ratio_minloHJBorder.SetFillStyle(0);
    # g_ratio_minloHJBorder.SetFillColor(ROOT.kOrange+2)
    # g_ratio_minloHJBorder.SetLineColor(ROOT.kOrange+2)
    # g_ratio_minloHJBorder.SetMarkerColor(ROOT.kOrange+2)

    # g_ratio_minloHJe0 = TGraphAsymmErrors(v_observable,v_ratio_minloHJ,v_dobservable,v_dobservable,v_zeros,v_zeros)
    # g_ratio_minloHJe0.SetFillStyle(3245);
    # g_ratio_minloHJe0.SetFillColor(ROOT.kOrange+2)
    # g_ratio_minloHJe0.SetLineColor(ROOT.kOrange+2)
    # g_ratio_minloHJe0.SetLineWidth(2)
    # g_ratio_minloHJe0.SetMarkerColor(ROOT.kOrange+2)

    if (obsName.startswith("mass4l")):

        g_ratio_minloHJ = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
        g_ratio_minloHJ.SetFillStyle(3245);
        g_ratio_minloHJ.SetFillColor(ROOT.kOrange+2)
        g_ratio_minloHJ.SetLineColor(ROOT.kOrange+2)
        g_ratio_minloHJ.SetLineWidth(2)
        g_ratio_minloHJ.SetMarkerColor(ROOT.kOrange+2)

        g_ratio_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
        g_ratio_minloHJBorder.SetFillStyle(0);
        g_ratio_minloHJBorder.SetFillColor(ROOT.kOrange+2)
        g_ratio_minloHJBorder.SetLineColor(ROOT.kOrange+2)
        g_ratio_minloHJBorder.SetMarkerColor(ROOT.kOrange+2)

        h_ratio_minloHJ = TH1D("h_ratio_minloHJ","h_ratio_minloHJ",4, 0, 4)
        for i in range(4):
            h_ratio_minloHJ.SetBinContent(i+1,v_ratio_minloHJ[i])
    elif doubleDiff:

        g_ratio_minloHJ = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
        g_ratio_minloHJ.SetFillStyle(3245);
        g_ratio_minloHJ.SetFillColor(ROOT.kOrange+2)
        g_ratio_minloHJ.SetLineColor(ROOT.kOrange+2)
        g_ratio_minloHJ.SetLineWidth(2)
        g_ratio_minloHJ.SetMarkerColor(ROOT.kOrange+2)

        g_ratio_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
        g_ratio_minloHJBorder.SetFillStyle(0);
        g_ratio_minloHJBorder.SetFillColor(ROOT.kOrange+2)
        g_ratio_minloHJBorder.SetLineColor(ROOT.kOrange+2)
        g_ratio_minloHJBorder.SetMarkerColor(ROOT.kOrange+2)

        h_ratio_minloHJ = TH1D("h_ratio_minloHJ","h_ratio_minloHJ",4, 0, 4)
        for i in range(4):
            h_ratio_minloHJ.SetBinContent(i+1,v_ratio_minloHJ[i])
    elif acFlag: # We keep the same names as in the last else (just to reduce entropy with if statements)
        g_ratio_minloHJ = TGraphAsymmErrors(v_observable_2,v_ratio_AC,v_dobservable_2,v_dobservable_2,v_ratio_AC_lo,v_ratio_AC_hi)
        g_ratio_minloHJ.SetFillStyle(3245);
        g_ratio_minloHJ.SetFillColor(ROOT.kRed+2)
        g_ratio_minloHJ.SetLineColor(ROOT.kRed+2)
        g_ratio_minloHJ.SetLineWidth(2)
        g_ratio_minloHJ.SetMarkerColor(ROOT.kRed+2)

        g_ratio_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ratio_AC,v_dobservable_2,v_dobservable_2,v_ratio_AC_lo,v_ratio_AC_hi)
        g_ratio_minloHJBorder.SetFillStyle(0);
        g_ratio_minloHJBorder.SetFillColor(ROOT.kRed+2)
        g_ratio_minloHJBorder.SetLineColor(ROOT.kRed+2)
        g_ratio_minloHJBorder.SetMarkerColor(ROOT.kRed+2)

        if acFlagBis:
            g_ratioBis_minloHJ = TGraphAsymmErrors(v_observable_2,v_ratio_ACbis,v_dobservable_2,v_dobservable_2,v_ratio_ACbis_lo,v_ratio_ACbis_hi)
            g_ratioBis_minloHJ.SetFillStyle(3245);
            g_ratioBis_minloHJ.SetFillColor(ROOT.kViolet+2)
            g_ratioBis_minloHJ.SetLineColor(ROOT.kViolet+2)
            g_ratioBis_minloHJ.SetLineWidth(2)
            g_ratioBis_minloHJ.SetMarkerColor(ROOT.kViolet+2)

            g_ratioBis_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ratio_ACbis,v_dobservable_2,v_dobservable_2,v_ratio_ACbis_lo,v_ratio_ACbis_hi)
            g_ratioBis_minloHJBorder.SetFillStyle(0);
            g_ratioBis_minloHJBorder.SetFillColor(ROOT.kViolet+2)
            g_ratioBis_minloHJBorder.SetLineColor(ROOT.kViolet+2)
            g_ratioBis_minloHJBorder.SetMarkerColor(ROOT.kViolet+2)

        if ("jet" in obsName and (not obsName.startswith("njets"))):
            h_ratio_minloHJ = TH1D("h_ratio_minloHJ","h_ratio_minloHJ",nBins-2, array('d',[float(obs_bins[i]) for i in range(1,len(obs_bins))]) )
            for i in range(1,nBins-1):
                h_ratio_minloHJ.SetBinContent(i,v_ratio_AC[i])
        else:
            h_ratio_minloHJ = TH1D("h_ratio_minloHJ","h_ratio_minloHJ",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            for i in range(nBins-1):
                h_ratio_minloHJ.SetBinContent(i+1,v_ratio_AC[i])
            if acFlagBis:
                h_ratioBis_minloHJ = TH1D("h_ratioBis_minloHJ","h_ratioBis_minloHJ",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
                for i in range(nBins-1):
                    h_ratioBis_minloHJ.SetBinContent(i+1,v_ratio_ACbis[i])
    else:
        g_ratio_minloHJ = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
        g_ratio_minloHJ.SetFillStyle(3245);
        g_ratio_minloHJ.SetFillColor(ROOT.kOrange+2)
        g_ratio_minloHJ.SetLineColor(ROOT.kOrange+2)
        g_ratio_minloHJ.SetLineWidth(2)
        g_ratio_minloHJ.SetMarkerColor(ROOT.kOrange+2)

        g_ratio_minloHJBorder = TGraphAsymmErrors(v_observable_2,v_ratio_minloHJ,v_dobservable_2,v_dobservable_2,v_ratio_minloHJ_lo,v_ratio_minloHJ_hi)
        g_ratio_minloHJBorder.SetFillStyle(0);
        g_ratio_minloHJBorder.SetFillColor(ROOT.kOrange+2)
        g_ratio_minloHJBorder.SetLineColor(ROOT.kOrange+2)
        g_ratio_minloHJBorder.SetMarkerColor(ROOT.kOrange+2)

        if ("jet" in obsName and (not obsName.startswith("njets"))):
            h_ratio_minloHJ = TH1D("h_ratio_minloHJ","h_ratio_minloHJ",nBins-2, array('d',[float(obs_bins[i]) for i in range(1,len(obs_bins))]) )
            for i in range(1,nBins-1):
                h_ratio_minloHJ.SetBinContent(i,v_ratio_minloHJ[i])
        else:
            h_ratio_minloHJ = TH1D("h_ratio_minloHJ","h_ratio_minloHJ",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            for i in range(nBins-1):
                h_ratio_minloHJ.SetBinContent(i+1,v_ratio_minloHJ[i])

    if not acFlag: h_ratio_minloHJ.SetLineColor(ROOT.kOrange+2)
    else:
        h_ratio_minloHJ.SetLineColor(ROOT.kRed+2)
        if acFlagBis:
            h_ratioBis_minloHJ.SetLineColor(ROOT.kViolet+2)
            h_ratioBis_minloHJ.SetLineWidth(2)
    h_ratio_minloHJ.SetLineWidth(2)

    g_modeldep = TGraphAsymmErrors(v_observable,v_data,v_twos,v_twos,v_modeldep_lo,v_modeldep_hi)
    #g_modeldep = TGraphAsymmErrors(v_observable,v_data,v_twos,v_twos,v_modeldep_lo,v_modeldep_hi)
    g_systematics = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_systematics_lo,v_systematics_hi)

    #g_modeldep.SetLineWidth(5)
    #g_modeldep.SetMarkerColor(ROOT.kGray)
    g_modeldep.SetFillColor(ROOT.kBlack)
    g_modeldep.SetLineColor(ROOT.kBlack)

    g_systematics.SetLineWidth(5)
    g_systematics.SetMarkerColor(ROOT.kRed)
    g_systematics.SetLineColor(ROOT.kRed)
    g_systematics.SetFillColor(ROOT.kRed)

    g_data_allunc = TGraphAsymmErrors(v_observable,v_data,v_zeros,v_zeros,v_data_lo_allunc,v_data_hi_allunc)
    g_data_allunc.SetMarkerColor(ROOT.kBlack)
    g_data_allunc.SetLineColor(ROOT.kBlack)
    g_data_allunc.SetLineWidth(1)
    g_data_allunc.SetMarkerStyle(20)
    g_data_allunc.SetMarkerSize(1.2)

    if (obsName=="pT4l"):
        label="p_{T}(H)"
        unit="GeV"
    elif (obsName=="massZ2"):
        label = "m(Z_{2})"
        unit = "GeV"
    elif (obsName=="massZ1"):
        label = "m(Z_{1})"
        unit = "GeV"
    elif (obsName=="nJets" or obsName=="njets_pt30_eta4p7"):
        label = "N(jets)"
        unit = ""
    elif (obsName=="njets_pt30_eta2p5"):
        label = "N(jets)"
        unit = ""
    elif (obsName=="pTj1"):
        label = "p_{T}(jet)"
        unit = "GeV"
    elif (obsName=="pTj1jet_pt30_eta2p5"):
        label = "p_{T}(jet)"
        unit = "GeV"
    elif (obsName=="absrapidity_leadingjet_pt30_eta4p7"):
        label = "|y(jet)|"
        unit = ""
    elif (obsName=="absrapidity_leadingjet_pt30_eta2p5"):
        label = "|y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=="absdeltarapidity_hleadingjet_pt30_eta4p7"):
        label = "|y(H)-y(jet)|"
        unit = ""
    elif (obsName=="absdeltarapidity_hleadingjet_pt30_eta2p5"):
        label = "|y(H)-y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=="rapidity4l"):
        label = "|y(H)|"
        unit = ""
    elif (obsName=="costhetastar"):
        label = "|cos#theta*|"
        unit = ""
    elif (obsName=="costhetaZ1"):
        label = "|cos#theta_{1}|"
        unit = ""
    elif (obsName=="costhetaZ2"):
        label = "|cos#theta_{2}|"
        unit = ""
    elif (obsName=="phi"):
        label = "|#Phi|"
        unit = ""
    elif (obsName=="phistar"):
        label = "|#Phi_{1}|"
        unit = ""
    elif (obsName.startswith("mass4l")):
        label = "inclusive"
        unit = ""
    elif obsName=='massZ1_massZ2':
        label = "m(Z_{1})"
        unit = "GeV"
        label_2nd = "m(Z_{2})"
        unit_2nd = "GeV"
    elif obsName=="D0m":
        label = "D_{0-}^{dec}"
        unit = ""
    elif obsName=="D0hp":
        label = "D_{0h+}^{dec}"
        unit = ""
    elif obsName=="Dcp":
        label = "D_{CP}^{dec}"
        unit = ""
    elif obsName=="Dint":
        label = "D_{int}^{dec}"
        unit = ""
    elif obsName=="DL1":
        label = "D_{#Lambda 1}^{dec}"
        unit = ""
    else:
        label = obsName
        unit = ""

    c = TCanvas("c",obsName, 1400, 1400)
    if(opt.SETLOG): c.SetLogy()
    if (not obsName.startswith("mass4l")): c.SetBottomMargin(0.35)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.18)

    if (obsName.startswith("mass4l")):
        dummy = TH1D("dummy","dummy", 4, 0, 4)
        for i in range(1,5):
            dummy.SetBinContent(i,2.5*max(a_ggH_powheg))
        h_XH = TH1D("h_XH","h_XH",4, 0, 4)
        print "mass4l",a_XH
        for i in range(4):
            h_XH.SetBinContent(i+1,a_XH[i])
    elif doubleDiff:
        dummy = TH1D("dummy","dummy", nBins-1, 0, nBins-1)
        for i in range(1,nBins-1):
            dummy.SetBinContent(i,2.5*max(a_ggH_powheg))
        h_XH = TH1D("h_XH","h_XH",nBins-1, 0, nBins-1)
        print "mass4l",a_XH
        for i in range(nBins-1):
            h_XH.SetBinContent(i+1,a_XH[i])
    else:
        if ("jet" in obsName and (not obsName.startswith("njets"))):
            dummy = TH1D("dummy","dummy", int(float(obs_bins[nBins-1])-float(obs_bins[1])), float(obs_bins[1]), float(obs_bins[nBins-1]))
            for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[1]))):
                dummy.SetBinContent(i,2.5*max(a_ggH_powheg))
            #h_XH = TH1D("h_XH","h_XH",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            #for i in range(nBins-1):
            #    h_XH.SetBinContent(i+1,a_XH[i])
            h_XH = TH1D("h_XH","h_XH",nBins-2, array('d',[float(obs_bins[i]) for i in range(1,len(obs_bins))]) )
            for i in range(1,nBins-1):
                h_XH.SetBinContent(i,a_XH[i])
        else:
            dummy = TH1D("dummy","dummy", int(float(obs_bins[nBins-1])-float(obs_bins[0])), float(obs_bins[0]), float(obs_bins[nBins-1]))
            for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[0]))):
                dummy.SetBinContent(i,2.5*max(a_ggH_powheg))
            h_XH = TH1D("h_XH","h_XH",nBins-1, array('d',[float(obs_bins[i]) for i in range(len(obs_bins))]) )
            for i in range(nBins-1):
                h_XH.SetBinContent(i+1,a_XH[i])

    if(opt.SETLOG):
        if (obsName=="massZ2"): dummy.SetMaximum(500.0*max(max(a_data),(max(a_ggH_powheg)))) #AT
        elif (obsName=="massZ1"): dummy.SetMaximum(500.0*max(max(a_data),(max(a_ggH_powheg)))) #AT
        elif doubleDiff: dummy.SetMaximum(300*max(max(a_data),(max(a_ggH_powheg))))
        elif (obsName=="D0m"): dummy.SetMaximum(1100*max(max(a_data),(max(a_ggH_powheg))))
        elif (obsName=="Dcp"): dummy.SetMaximum(1100*max(max(a_data),(max(a_ggH_powheg))))
        elif (obsName=="D0hp"): dummy.SetMaximum(1100*max(max(a_data),(max(a_ggH_powheg))))
        elif (obsName=="DL1"): dummy.SetMaximum(1100*max(max(a_data),(max(a_ggH_powheg))))
        else: dummy.SetMaximum(55.0*max(max(a_data),(max(a_ggH_powheg))))
    else:
        if (obsName.startswith("mass4l")) or doubleDiff: dummy.SetMaximum(1.6*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi)))
        if doubleDiff: dummy.SetMaximum(2.0*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi)))
        elif (obsName=="massZ2"): dummy.SetMaximum(2.0*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi))) #AT
        elif (obsName=="massZ1"): dummy.SetMaximum(2.0*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi))) #AT
        elif (obsName=="D0m"): dummy.SetMaximum(2.0*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi))) #AT
        elif (obsName=="Dcp"): dummy.SetMaximum(2.0*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi))) #AT
        elif (obsName=="D0hp"): dummy.SetMaximum(2.5*(max(max(a_data),(max(a_ggH_powheg)))+max(a_data_hi))) #AT
        elif (obsName=="DL1"): dummy.SetMaximum(2.8*max(max(a_data),(max(a_ggH_powheg))))
	else: dummy.SetMaximum(1.5*(max(max(a_ggH_powheg),(max(a_data)+max(a_data_hi)))))
    if (opt.SETLOG):
        dummy.SetMinimum(0.0601*max(min(a_data),(min(a_ggH_powheg))))
        if (obsName.startswith("pTj1")): dummy.SetMinimum(0.0801*max(min(a_data),(min(a_ggH_powheg))))
    else: dummy.SetMinimum(0.0001)
    dummy.SetLineColor(0)
    dummy.SetMarkerColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.GetXaxis().SetLabelSize(0.0)
    dummy.GetYaxis().SetLabelSize(0.04)
    if (opt.SETLOG and (obsName.startswith('njets') or obsName.startswith('pTj1'))):
        dummy.SetMaximum(200.0*max(max(a_data),(max(a_ggH_powheg))))
        dummy.GetXaxis().SetTitle("")
    else:
        dummy.GetXaxis().SetTitle("")
    if (obsName.startswith('njets')):
        dummy.GetXaxis().SetTitle('')
        dummy.GetXaxis().SetLabelSize(0.0)
        dummy.GetXaxis().SetBinLabel(1,'')
        dummy.GetXaxis().SetBinLabel(2,'')
        dummy.GetXaxis().SetBinLabel(3,'')
        dummy.GetXaxis().SetBinLabel(4,'')
    if (obsName.startswith("mass4l")):
        dummy.GetXaxis().SetLabelSize(0.08)
        dummy.GetXaxis().SetBinLabel(1,'4l')
        dummy.GetXaxis().SetBinLabel(2,'2e2#mu')
        dummy.GetXaxis().SetBinLabel(3,'4#mu')
        dummy.GetXaxis().SetBinLabel(4,'4e')
    if doubleDiff:
        dummy.GetXaxis().SetLabelSize(0.08)
        dummy.GetXaxis().SetTitle('')
        for i in range(1,nBins-1):
            dummy.GetXaxis().SetBinLabel(i,'')
        # dummy.GetXaxis().SetBinLabel(1,'')
        # dummy.GetXaxis().SetBinLabel(2,'')
        # dummy.GetXaxis().SetBinLabel(3,'')
        # dummy.GetXaxis().SetBinLabel(4,'')
        # dummy.GetXaxis().SetBinLabel(5,'')

    dummy.GetXaxis().SetTitleSize(0.0)

    if (label=="inclusive"):
        dummy.GetYaxis().SetTitle("#sigma_{fid} (fb)")
    elif (unit==''):
        dummy.GetYaxis().SetTitle("d#sigma_{fid }/d"+label+" (fb)")
    elif doubleDiff:
        if unit==unit_2nd:
            dummy.GetYaxis().SetTitle("d^{2}#sigma_{fid }/d"+label+"d"+label_2nd+" (fb/"+unit+"^{2})")
        else:
            dummy.GetYaxis().SetTitle("d^{2}#sigma_{fid }/d"+label+"d"+label_2nd+" (fb/"+unit+unit_2nd+")")
        dummy.GetYaxis().SetTitleSize(0.04)
    else:
        dummy.GetYaxis().SetTitle("d#sigma_{fid }/d"+label+" (fb/"+unit+")")
    if (obsName.startswith('njets')):
        dummy.GetYaxis().SetTitle("#sigma_{fid} (fb)")

    dummy.GetYaxis().SetTitleOffset(1.4)
    #dummy.GetYaxis().SetTitleOffset(1.8)
    dummy.Draw("hist")

    h_XH.SetFillColor(kGreen-8)
    h_XH.SetFillStyle(3344)

    legend = TLegend(.28,.65,.85,.90)
    #legend . SetNColumns(2)
    #legend . SetTextSize(0.026)
    #legend . SetTextSize(0.029) # less info for XH
    legend . SetTextSize(0.025) # shrink margins
    if (opt.UNBLIND): legend . AddEntry(g_data , "Data (stat. #oplus sys. unc.)", "ep")
    #else: legend . AddEntry(g_data , "Asimov Data (stat.#oplussys. unc.)", "ep")
    else: legend . AddEntry(g_data , "Toy Data (stat. #oplus sys. unc.)", "ep")
    legend . AddEntry(g_systematics,"Systematic uncertainty","l")
    # legend . AddEntry(g_modeldep,"Model dependence","f")
    if not acFlag:
        legend . AddEntry(g_ggH_minloHJ , "gg#rightarrowH (NNLOPS + JHUGen) + XH", "lf")
        legend . AddEntry(g_ggH_powheg , "gg#rightarrowH (POWHEG + JHUGen) + XH", "lf")
    else:
        legend . AddEntry(g_ggH_powheg , "SM (POWHEG + JHUGen + Pythia)", "lf")
        if (obsName=="D0m"): legend . AddEntry(g_ggH_AC , "AC (POWHEG + JHUGen fa3=1 + Pythia)", "lf")
        if (obsName=="Dcp"): legend . AddEntry(g_ggH_AC , "AC (POWHEG + JHUGen fa3=0.5 + Pythia)", "lf")
        if (obsName=="Dint"): legend . AddEntry(g_ggH_AC , "AC (POWHEG + JHUGen fa2=0.5 + Pythia)", "lf")
        if (obsName=="D0hp"): legend . AddEntry(g_ggH_AC , "AC (POWHEG + JHUGen fa2=1 + Pythia)", "lf")
        if (obsName=="DL1"):
            legend . AddEntry(g_ggH_AC , "AC (POWHEG + JHUGen faL1=1 + Pythia)", "lf")
            legend . AddEntry(g_ggH_ACbis , "AC (POWHEG + JHUGen faL1=0.5 + Pythia)", "lf")
    #legend . AddEntry(g_XH , "XH = VBF + VH + ttH", "l")
    if not acFlag: legend . AddEntry(h_XH , "XH = VBF + VH + ttH (POWHEG + JHUGen)", "f")
    legend . AddEntry(dummy, "(LHCHWG YR4, m_{H}="+opt.THEORYMASS+" GeV)", "")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw()


    h_ggH_powheg.Draw("histsame")
    if acFlag:
        h_ggH_AC.Draw("histsame")
        if acFlagBis:
            h_ggH_ACbis.Draw("histsame")
    else: h_ggH_minloHJ.Draw("histsame")
    g_ggH_powheg.Draw("5same")
    g_ggH_powhegBorder.Draw("5same")
    #g_ggH_powhege0.Draw("epsame")
    if acFlag:
        g_ggH_AC.Draw("5same")
        g_ggH_ACBorder.Draw("5same")
        if acFlagBis:
            g_ggH_ACbis.Draw("5same")
            g_ggH_ACbisBorder.Draw("5same")
    else:
        g_ggH_minloHJ.Draw("5same")
        g_ggH_minloHJBorder.Draw("5same")
    #g_ggH_minloHJe0.Draw("epsame")

    #g_XH.Draw("2same")
    #g_XH.Draw("fsame")
    if not acFlag: h_XH.Draw("histsame")
    # g_modeldep.Draw("2same0")
    # g_modeldep.Draw("psame")
    g_data.Draw("psameZ0")
    g_systematics.Draw("psameZ0")
    g_data_e0.Draw("psame")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    # print opt.LUMISCALE
    # if (not opt.LUMISCALE=="1.0"):
    #     lumi = round(137.1*float(opt.LUMISCALE),1)
    #     latex2.DrawLatex(0.94, 0.94,str(lumi)+" fb^{-1} (13 TeV)")
    # else:
    if(opt.YEAR=='2016'):
        latex2.DrawLatex(0.92, 0.95,"35.9 fb^{-1} (13 TeV)")
    if(opt.YEAR=='2017'):
        latex2.DrawLatex(0.92, 0.95,"41.5 fb^{-1} (13 TeV)")
    if(opt.YEAR=='2018'):
        latex2.DrawLatex(0.92, 0.95,"59.7 fb^{-1} (13 TeV)")
    if(opt.YEAR=='Full'):
        latex2.DrawLatex(0.92, 0.95,"137 fb^{-1} (13 TeV)")
    latex2.SetTextSize(0.7*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.19, 0.94, "CMS #bf{#it{Preliminary}}") #AT Tolta scritta CMS
    latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    #latex2.DrawLatex(0.28, 0.945, "Unpublished")
    #latex2.DrawLatex(0.30, 0.94, "Preliminary")

    latex2.SetTextSize(0.4*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    #latex2.DrawLatex(0.55, 0.67,"N^{3}LO #sigma_{gg#rightarrowH}^{total}")

    if ("2p5" in obsName):
        latex2.SetTextSize(0.4*c.GetTopMargin())
        latex2.SetTextFont(42)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.92, 0.61,"p_{T}(jet) > 30 GeV, |#eta(jet)| < 2.5")
    if ("4p7" in obsName):
        latex2.SetTextSize(0.4*c.GetTopMargin())
        latex2.SetTextFont(42)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.92, 0.61,"p_{T}(jet) > 30 GeV, |#eta(jet)| < 4.7")

    if (obsName=="pT4l"):
        latex2.SetTextSize(0.35*c.GetTopMargin())
        latex2.SetTextFont(42)
        latex2.SetTextAlign(31) # align right
        latex2.SetTextAngle(90)
        latex2.DrawLatex(0.92, 0.77,"#frac{1}{50} #sigma(p_{T}(H) > 200 GeV)")



    dummy.Draw("axissame")

    if (not obsName.startswith("mass4l")):
        pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
        pad.SetTopMargin(0.65)
        pad.SetRightMargin(0.04)
        pad.SetLeftMargin(0.18)
        pad.SetFillColor(0)
        pad.SetGridy(1)
        pad.SetFillStyle(0)
        pad.Draw()
        pad.cd(0)

        if (obsName.startswith("mass4l")):
            dummy2 = TH1D("dummy2","dummy2", 4, 0, 4)
            for i in range(1,5):
                dummy2.SetBinContent(i,1.02)
            dummy2.GetXaxis().SetLabelSize(0.08)
            dummy2.GetXaxis().SetBinLabel(1,'4l')
            dummy2.GetXaxis().SetBinLabel(2,'2e2#mu')
            dummy2.GetXaxis().SetBinLabel(3,'4#mu')
            dummy2.GetXaxis().SetBinLabel(4,'4e')
        elif doubleDiff:
            dummy2 = TH1D("dummy2","dummy2", nBins-1, 0, nBins-1)
            for i in range(1,nBins-1):
                dummy2.SetBinContent(i,1.02)
            dummy2.GetXaxis().SetTitle("")
            dummy2.GetXaxis().SetLabelSize(0.08)
            for i in range(1,nBins):
                dummy2.GetXaxis().SetBinLabel(i,'Bin '+str(i-1))
                # dummy2.GetXaxis().SetBinLabel(1,'0')
                # dummy2.GetXaxis().SetBinLabel(2,'1')
                # dummy2.GetXaxis().SetBinLabel(3,'2')
                # dummy2.GetXaxis().SetBinLabel(4,'3')
                # dummy2.GetXaxis().SetBinLabel(5,'4')
        else:
            if ("jet" in obsName and (not obsName.startswith("njets"))):
                dummy2 = TH1D("dummy2","dummy2", int(float(obs_bins[nBins-1])-float(obs_bins[1])), float(obs_bins[1]), float(obs_bins[nBins-1]))
                for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[1]))):
                    dummy2.SetBinContent(i,1.02)
            else:
                if (obsName=="pT4l"):
                    dummy2 = TH1D("dummy2","dummy2", int(float(obs_bins[nBins-1])-float(obs_bins[0])), float(obs_bins[0]), float(obs_bins[nBins-1]))
                else:
                    dummy2 = TH1D("dummy2","dummy2", int(float(obs_bins[nBins-1])-float(obs_bins[0])), float(obs_bins[0]), float(obs_bins[nBins-1]))
                for i in range(int(float(obs_bins[nBins-1])-float(obs_bins[0]))):
                    dummy2.SetBinContent(i,1.02)
        if doubleDiff:dummy2.GetXaxis().SetLabelSize(0.08)
        else: dummy2.GetXaxis().SetLabelSize(0.04)

        dummy2.SetLineColor(0)
        dummy2.SetMarkerColor(0)
        dummy2.SetLineWidth(0)
        dummy2.SetMarkerSize(0)
        if (obsName.startswith('njets')):
            dummy2.GetXaxis().SetTitle(label)
            dummy2.GetXaxis().SetLabelSize(0.08)
            dummy2.GetXaxis().SetBinLabel(1,'0')
            dummy2.GetXaxis().SetBinLabel(2,'1')
            dummy2.GetXaxis().SetBinLabel(3,'2')
            dummy2.GetXaxis().SetBinLabel(4,' 3')
            dummy2.GetXaxis().SetBinLabel(5,'#geq 4')
        elif (label=="inclusive"):
            dummy2.GetXaxis().SetTitle("")
        elif (unit==""):
            dummy2.GetXaxis().SetTitle(label)
        elif doubleDiff:
            dummy2.GetXaxis().SetTitle("")
        else:
            dummy2.GetXaxis().SetTitle(label+" ("+unit+")")

        dummy2.GetYaxis().SetLabelSize(0.03)
        dummy2.GetYaxis().SetNdivisions(10);
        if (obsName.startswith("pTj1")): dummy2.GetYaxis().SetNdivisions(8);
        dummy2.GetXaxis().SetNdivisions(510)
        ratiomax=1.0
        for i in range(len(data_hi)):
          if ((data[i]+data_hi[i])/ggH_minloHJ[i]>ratiomax): ratiomax=(data[i]+data_hi[i])/ggH_minloHJ[i]
        dummy2.SetMaximum(ratiomax*1.1)
        if (obsName=="pT4l"):
            dummy2.SetMaximum(ratiomax*1.04)
            dummy2.GetYaxis().SetNdivisions(508)
        dummy2.SetMinimum(0.0)
        dummy2.Draw("hist")
        dummy2.GetYaxis().CenterTitle()
        dummy2.GetYaxis().SetTitleSize(0.04)
        dummy2.GetYaxis().SetTitleOffset(1.5)
        dummy2.GetYaxis().SetTitleSize(0.03)
        dummy2.GetYaxis().SetTitleOffset(2.0)
        if acFlag: dummy2.GetYaxis().SetTitle('Ratio to SM')
        else: dummy2.GetYaxis().SetTitle('Ratio to NNLOPS')


        h_ratio_powheg.Draw("histsame")
        h_ratio_minloHJ.Draw("histsame")
        if acFlagBis: h_ratioBis_minloHJ.Draw("histsame")

        g_ratio_minloHJ.Draw("5same")
        g_ratio_minloHJBorder.Draw("5same")
        if acFlagBis:
            g_ratioBis_minloHJ.Draw("5same")
            g_ratioBis_minloHJBorder.Draw("5same")
        #g_ratio_minloHJe0.Draw("epsame")

        g_ratio_powheg.Draw("5same")
        g_ratio_powhegBorder.Draw("5same")
        #g_ratio_powhege0.Draw("epsame")

        g_ratio_data.Draw("psameZ0")
        g_ratio_sys.Draw("psameZ0")
        g_ratio_datae0.Draw("psame")

        dummy2.GetYaxis().SetRangeUser(0,3)
        dummy2.Draw("axissame")

    if (obsName=="pT4l"):
        box = TPaveText(240,-0.4,260,-0.1)
        box.SetLineColor(0)
        box.SetFillColor(0)
        box.AddText("   ")
        box.Draw("same")

    #if (opt.SETLOG): set_log = '_HResTESTlogscale'
    #else: set_log = 'HResTEST'

    if (opt.SETLOG): set_log = '_logscale'
    else: set_log = ''

    if (not opt.UNBLIND): set_log = set_log + '_asimov'

    if (opt.UNBLIND): subdir = 'data'
    elif(not opt.UNBLIND): subdir = 'asimov'

    checkDir('plots/'+obsName)
    checkDir('plots/'+obsName+'/'+subdir)

    c.SaveAs('plots/'+obsName+'/'+subdir+'/'+obsName+'_unfoldwith_'+datamodel+set_log+'.pdf')
    c.SaveAs('plots/'+obsName+'/'+subdir+'/'+obsName+'_unfoldwith_'+datamodel+set_log+'.png')
    c.SaveAs('plots/'+obsName+'/'+subdir+'/'+obsName+'_unfoldwith_'+datamodel+set_log+'.root')
    c.SaveAs('plots/'+obsName+'/'+subdir+'/'+obsName+'_unfoldwith_'+datamodel+set_log+'.C')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.pdf')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.png')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.root')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_unpublished.C')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.pdf')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.png')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.root')
    #c.SaveAs('plots/'+obsName+'_unfoldwith_'+datamodel+set_log+'_preliminary.C')
    # with open('table/table_'+opt.OBSNAME+'_'+subdir+'.py', 'w') as f:
    #      f.write('\\documentclass{article} \n')
    #      f.write('\\begin{document} \n')
    #      f.write('\\begin{table}[!h!tb] \n')
    #      f.write('\\begin{center} \n')
    #      f.write('\\begin{tabular}{|l|')
    #      for columns in range(0,len(obs_bins)-1):
    #          f.write('c|')
    #      f.write('} \\hline \\hline \n')
    #      f.write('Observable & ' )
    #      for obsbin in range(0,len(obs_bins)-1):
    #          f.write('$'+obs_bins[obsbin]+'-'+obs_bins[obsbin+1]+'$')
    #          if (obsbin == len(obs_bins)-2):
    #              f.write(' \\\\ \\hline \n')
    #          else:
    #              f.write (' & ')
    #      f.write(opt.OBSNAME+' & \n')
    #      rbin =0
    #      for obsbin in range(0,len(obs_bins)-1):
    #          f.write('$'+str(round(a_data[obsbin],3))+'^{+'+str(abs(round(sqrt(a_data_hi[obsbin]**2-a_systematics_hi[obsbin]**2),3)))+'+'+str(abs(round(a_systematics_hi[obsbin],3)))+'}'+'_{-'+str(abs(round(sqrt(a_data_lo[obsbin]**2-a_systematics_lo[obsbin]**2),3)))+'-'+str(abs(round(a_systematics_lo[obsbin],3)))+'}$')
    #          if (obsbin == len(obs_bins)-2):
    #              f.write (' \n ')
    #          else:
    #              f.write (' & ')
    #      f.write(' \\\\ \\hline \n')
    #      f.write('\\end{tabular} \n')
    #      f.write('\\end{center} \n')
    #      f.write('\\end{table} \n')
    #      f.write('\\end{document} \n')
doubleDiff = False
obs_name = opt.OBSNAME
if 'vs' in opt.OBSNAME:
    doubleDiff = True
    obs_name_tmp = opt.OBSNAME.split(' vs ')
    obs_name = obs_name_tmp[0]+'_'+obs_name_tmp[1]

_temp = __import__('inputs_sig_'+obs_name+'_'+opt.YEAR, globals(), locals(), ['observableBins'], -1)
obs_bins = _temp.observableBins
print(obs_bins)

if obs_name == 'D0m':
    acFlag = True
    acFlagBis = False
    acSample = '0M'
elif obs_name == 'Dcp':
    acFlag = True
    acFlagBis = False
    acSample = '0Mf05ph0'
elif obs_name == 'D0hp':
    acFlag = True
    acFlagBis = False
    acSample = '0PH'
elif obs_name == 'Dint':
    acFlag = True
    acFlagBis = False
    acSample = '0PHf05ph0'
elif obs_name == 'DL1':
    acFlag = True
    acFlagBis = True
    acSample = '0L1'
    acSampleBis = '0L1f05ph0'
else:
    acFlag = False
    acFlagBis = False

if not doubleDiff:
    if float(obs_bins[len(obs_bins)-1])>300.0:
        obs_bins[len(obs_bins)-1]='250.0'
    if (opt.OBSNAME=="nJets" or opt.OBSNAME.startswith("njets")):
        obs_bins[len(obs_bins)-1]='5'
else:
    obs_bins_boundaries = obs_bins
    obs_bins = [i for i in range(len(obs_bins)+1)]

print('obs_bins', obs_bins)
if doubleDiff: plotXS(obs_name, obs_bins, obs_bins_boundaries)
else: plotXS(obs_name, obs_bins)
sys.path.remove('./inputs')
sys.path.remove('./LHScans')
