from scipy import interpolate
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep
import optparse,sys

# sys.path.append('../inputs/')
from observables import observables
from binning import binning

print ('Welcome in RunInterpolation!')

def parseOptions():

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--year',  dest='YEAR',  type='string', default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--extrapMass',  dest='EXTRAP',  type='float', default=125.38,   help='Mass at which extrapolate')
    parser.add_option('',   '--nnlops', action='store_true', dest='NNLOPS', default=False, help='Flag to extrapolate NNLOPS acceptance at 125.38')
    # The following two options are used together to calculate the acceptance in AC scenario to plot AC predictions on fiducial plot
    parser.add_option('',   '--AC_onlyAcc', action='store_true', dest='AC_ONLYACC', default=False, help='Flag in case we are interested in only the acceptance')
    parser.add_option('',   '--AC_hypothesis', dest='AC_HYP',  type='string',default='',   help='Name of the AC hypothesis, e.g. 0M, 0PM')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
    #     parser.error('Bin boundaries not specified for differential measurement. Exiting...')
    #     sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()


# -----------------------------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------------------------
# -----------------------------------------------------------------------------------------

obsName = opt.OBSNAME

acc = {}
err_acc = {}
eff = {}
err_eff = {}
outinratio = {}
err_outinratio = {}
inc_wrongfrac = {}
binfrac_wrongfrac = {}

mass_points = [124,125,126]

obs_bins, doubleDiff = binning(opt.OBSNAME)
if doubleDiff:
    obs_name = opt.OBSNAME.split(' vs ')[0]
    obs_name_2nd = opt.OBSNAME.split(' vs ')[1]
    obs_name_2d = opt.OBSNAME
else:
    obs_name = opt.OBSNAME

if doubleDiff: obsName = obs_name+'_'+obs_name_2nd
else: obsName = obs_name

if opt.NNLOPS and opt.YEAR!='Full':
    print('NNLOPS only with full run2')
    sys.exit()


sys.path.append('../inputs')
if opt.NNLOPS:
    nnlopsFlag = '_NNLOPS'
    inputs = ['acc']
    pModes = ['ggH']
else:
    nnlopsFlag = ''
    inputs = ['acc', 'err_acc', 'eff', 'err_eff', 'outinratio', 'err_outinratio', 'inc_wrongfrac', 'binfrac_wrongfrac']
    pModes = ['ggH', 'VBFH', 'ZH', 'WH']

for m in mass_points:
    if m == 125:
        _temp = __import__('inputs_sig_'+obsName+nnlopsFlag+'_'+opt.YEAR, globals(), locals(), inputs+['observableBins'])
        observableBins = _temp.observableBins
    else:
        _temp = __import__('inputs_sig_'+str(m)+'_'+obsName+nnlopsFlag+'_'+opt.YEAR, globals(), locals(), inputs)
    acc[m] = _temp.acc
    if opt.NNLOPS:
        continue
    # err_acc[m] = _temp.err_acc
    eff[m] = _temp.eff
    err_eff[m] = _temp.err_eff
    outinratio[m] = _temp.outinratio
    err_outinratio[m] = _temp.err_outinratio
    inc_wrongfrac[m] = _temp.inc_wrongfrac
    binfrac_wrongfrac[m] = _temp.binfrac_wrongfrac
sys.path.remove('../inputs')

spline = {}
extrap_acc = {}
extrap_eff = {}
extrap_outinratio = {}
extrap_inc_wrongfrac = {}
extrap_binfrac_wrongfrac = {}
nBins = len(obs_bins)
if not doubleDiff: nBins = len(obs_bins)-1 #In case of 1D measurement the number of bins is -1 the length of obs_bins(=bin boundaries)
for channel in ['2e2mu', '4e', '4mu']:
    for genBin in range(nBins):
        for recoBin in range(nBins):
            fig,axs = plt.subplots(2, 3, figsize=(30,10), dpi=80)
            axs.ravel()
            plt.style.use(hep.style.CMS)
            # hep.cms.text('Simulation')
            if doubleDiff: hep.cms.lumitext(obs_name+'vs'+obs_name_2d+' bin'+str(genBin)+' '+channel+' Run2')
            else: hep.cms.lumitext(obs_name+' bin'+str(genBin)+' '+channel+' Run2')
            for pMode in pModes:
                #Set the name of the proccesBin to 125 to reduce the changes in the other codes. it should be 125.38 since we are extrapolating
                if doubleDiff:
                    processBin_acc = pMode+'125'+nnlopsFlag+'_'+channel+'_'+obs_name+'_'+obs_name_2nd+'_genbin'+str(genBin)+'_recobin'+str(genBin)
                    processBin = pMode+'125'+nnlopsFlag+'_'+channel+'_'+obs_name+'_'+obs_name_2nd+'_genbin'+str(genBin)+'_recobin'+str(recoBin)
                else:
                    processBin_acc = pMode+'125'+nnlopsFlag+'_'+channel+'_'+obs_name+'_genbin'+str(genBin)+'_recobin'+str(genBin)
                    processBin = pMode+'125'+nnlopsFlag+'_'+channel+'_'+obs_name+'_genbin'+str(genBin)+'_recobin'+str(recoBin)

                #Acceptance
                acc_points = [acc[124][pMode+'124'+nnlopsFlag+'_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)],
                              acc[125][pMode+'125'+nnlopsFlag+'_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)],
                              acc[126][pMode+'126'+nnlopsFlag+'_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)]]
                spline_acc = interpolate.splrep(mass_points, acc_points, k=2)
                extrap_acc[processBin_acc] = float(interpolate.splev(opt.EXTRAP, spline_acc))
                # acc plot interpolation
                x = np.linspace(mass_points[0], mass_points[len(mass_points)-1], 200)
                y = interpolate.splev(x, spline_acc)
                axs[0,0].plot(mass_points, acc_points, 'ro')
                axs[0,0].plot(x, y, label=pMode)
                axs[0,0].set_xlabel('Mass points (GeV)')
                axs[0,0].set_ylabel('acceptance')
                axs[0,0].axvline(opt.EXTRAP, 0, 3000, ls='--', c='black')
                axs[0,0].legend()

                if opt.NNLOPS:
                    continue

                #Efficiency
                eff_points = [eff[124][pMode+'124_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                              eff[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                              eff[126][pMode+'126_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]]
                spline_eff = interpolate.splrep(mass_points, eff_points, k=2)
                extrap_eff[processBin] = float(interpolate.splev(opt.EXTRAP, spline_eff))
                # eff plot interpolation
                x = np.linspace(mass_points[0], mass_points[len(mass_points)-1], 200)
                y = interpolate.splev(x, spline_eff)
                axs[0,1].plot(mass_points, eff_points, 'ro')
                axs[0,1].plot(x, y, label=pMode)
                axs[0,1].set_xlabel('Mass points (GeV)')
                axs[0,1].set_ylabel('eff')
                axs[0,1].axvline(opt.EXTRAP, 0, 3000, ls='--', c='black')
                axs[0,1].legend()

                #Outinratio
                outinratio_points = [outinratio[124][pMode+'124_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                                     outinratio[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                                     outinratio[126][pMode+'126_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]]
                spline_outinratio = interpolate.splrep(mass_points, outinratio_points, k=2)
                extrap_outinratio[processBin] = float(interpolate.splev(opt.EXTRAP, spline_outinratio))
                # outinratio plot interpolation
                x = np.linspace(mass_points[0], mass_points[len(mass_points)-1], 200)
                y = interpolate.splev(x, spline_outinratio)
                axs[0,2].plot(mass_points, outinratio_points, 'ro')
                axs[0,2].plot(x, y, label=pMode)
                axs[0,2].set_xlabel('Mass points (GeV)')
                axs[0,2].set_ylabel('outinratio')
                axs[0,2].axvline(opt.EXTRAP, 0, 3000, ls='--', c='black')
                axs[0,2].legend()

                #inc_wrongfrac
                inc_wrongfrac_points = [inc_wrongfrac[124][pMode+'124_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                                        inc_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                                        inc_wrongfrac[126][pMode+'126_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]]
                spline_inc_wrongfrac = interpolate.splrep(mass_points, inc_wrongfrac_points, k=2)
                extrap_inc_wrongfrac[processBin] = float(interpolate.splev(opt.EXTRAP, spline_inc_wrongfrac))
                # inc_wrongfrac plot interpolation
                x = np.linspace(mass_points[0], mass_points[len(mass_points)-1], 200)
                y = interpolate.splev(x, spline_inc_wrongfrac)
                axs[1,0].plot(mass_points, inc_wrongfrac_points, 'ro')
                axs[1,0].plot(x, y, label=pMode)
                axs[1,0].set_xlabel('Mass points (GeV)')
                axs[1,0].set_ylabel('inc_wrongfrac')
                axs[1,0].axvline(opt.EXTRAP, 0, 3000, ls='--', c='black')
                axs[1,0].legend()

                #binfrac_wrongfrac
                binfrac_wrongfrac_points = [binfrac_wrongfrac[124][pMode+'124_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                                            binfrac_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)],
                                            binfrac_wrongfrac[126][pMode+'126_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]]
                spline_binfrac_wrongfrac = interpolate.splrep(mass_points, binfrac_wrongfrac_points, k=2)
                extrap_binfrac_wrongfrac[processBin] = float(interpolate.splev(opt.EXTRAP, spline_binfrac_wrongfrac))
                # binfrac_wrongfrac plot interpolation
                x = np.linspace(mass_points[0], mass_points[len(mass_points)-1], 200)
                y = interpolate.splev(x, spline_binfrac_wrongfrac)
                axs[1,1].plot(mass_points, binfrac_wrongfrac_points, 'ro')
                axs[1,1].plot(x, y, label=pMode)
                axs[1,1].set_xlabel('Mass points (GeV)')
                axs[1,1].set_ylabel('binfrac_wrongfrac')
                axs[1,1].axvline(opt.EXTRAP, 0, 3000, ls='--', c='black')
                axs[1,1].legend()


            plt.savefig('extrapolation/extrap_%s_%s_%s_Bin%i_%i%s.png' %(opt.YEAR, obsName, channel, genBin, recoBin, nnlopsFlag), bbox_inches='tight')
            # plt.show()
            plt.close()

            if opt.NNLOPS:
                continue
            diff = [(extrap_acc[pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)] - acc[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)]) / acc[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)] for pMode in ['ggH', 'VBFH', 'ZH', 'WH']]
            diff = np.mean(diff)
            extrap_acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)] = (1+diff) * acc[125]['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(genBin)]

            diff = [(extrap_eff[pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] - eff[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]) / eff[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] for pMode in ['ggH', 'VBFH', 'ZH', 'WH']]
            diff = np.mean(diff)
            extrap_eff['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] = (1+diff) * eff[125]['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]

            diff = [(extrap_outinratio[pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] - outinratio[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]) / outinratio[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] for pMode in ['ggH', 'VBFH', 'ZH', 'WH']]
            diff = np.mean(diff)
            extrap_outinratio['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] = (1+diff) * outinratio[125]['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]

            diff = [(extrap_inc_wrongfrac[pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] - inc_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]) / inc_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] for pMode in ['ggH', 'VBFH', 'ZH', 'WH']]
            diff = np.mean(diff)
            extrap_inc_wrongfrac['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] = (1+diff) * inc_wrongfrac[125]['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]

            diff = [(extrap_binfrac_wrongfrac[pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] - binfrac_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]) / binfrac_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] if binfrac_wrongfrac[125][pMode+'125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]>0 else 0 for pMode in ['ggH', 'VBFH', 'ZH', 'WH']]
            diff = np.mean(diff)
            extrap_binfrac_wrongfrac['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)] = (1+diff) * binfrac_wrongfrac[125]['ttH125_'+channel+'_'+obsName+'_genbin'+str(genBin)+'_recobin'+str(recoBin)]

# if doubleDiff: obs_name_dic = obs_name+'_'+obs_name_2nd
# else: obs_name_dic = obs_name
if opt.NNLOPS:
    with open('../inputs/inputs_sig_extrap_'+obsName+'_NNLOPS_'+str(opt.YEAR)+'.py', 'w') as f:
        f.write('observableBins = '+str(observableBins)+' \n')
        f.write('acc = '+str(extrap_acc))
else:
    with open('../inputs/inputs_sig_extrap_'+obsName+'_'+str(opt.YEAR)+'.py', 'w') as f:
        f.write('observableBins = '+str(observableBins)+' \n')
        f.write('acc = '+str(extrap_acc)+' \n')
        f.write('eff = '+str(extrap_eff)+' \n')
        f.write('err_eff = '+str(err_eff[125])+' \n')
        f.write('outinratio = '+str(extrap_outinratio)+' \n')
        f.write('err_outinratio = '+str(err_outinratio[125])+' \n')
        f.write('inc_wrongfrac = '+str(extrap_inc_wrongfrac)+' \n')
        f.write('binfrac_wrongfrac = '+str(extrap_binfrac_wrongfrac))
