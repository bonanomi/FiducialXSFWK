import numpy
import os
import json
import numpy
import optparse, os, sys
from binning import binning
from tabulate import tabulate

print 'Welcome in RunJES!'


def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
    #     parser.error('Bin boundaries not specified for differential measurement. Exiting...')
    #     sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()


obsname = opt.OBSNAME
obsname_out = obsname
doubleDiff = False

if ' vs ' in obsname:
    doubleDiff = True
    obsname_1st = opt.OBSNAME.split(' vs ')[0]
    obsname_2nd = opt.OBSNAME.split(' vs ')[1]
    obsname_out = obsname_1st + '_' + obsname_2nd

year = opt.YEAR

jesNames = ['Total', 'Abs', 'Abs_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavQCD', 'HF', 'HF_year', 'RelBal', 'RelSample_year']

print (obsname, year, obsname_out)
obs_bins, doubleDiff = binning(opt.OBSNAME)

_temp = __import__('JESNP_evts_'+obsname_out, globals(), locals(), ['evts','evts_noWeight'])
evts = _temp.evts
evts_noWeight = _temp.evts_noWeight
_temp = __import__('JESNP_'+obsname_out, globals(), locals(), ['JESNP'])
JESNP = _temp.JESNP

# Check if the folder for tables exist
if not os.path.exists('tables/'+obsname_out):
    os.makedirs('tables/'+obsname_out)

if doubleDiff: nBins = len(obs_bins)
else: nBins = len(obs_bins)-1

# Tables with numerical values
tables = {}
inclusiveJES = {}
for fState in ['2e2mu', '4e', '4mu']:
    for recobin in range(nBins):
        table = []
        for jesName in jesNames:
            # nominal_incl = 0
            # up_incl = 0
            # dn_incl = 0
            # nominal_incl += evts['signal_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            # up_incl += evts['signal_jesup_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            # dn_incl += evts['signal_jesdn_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            table.append(['signal_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin),
                          evts_noWeight['signal_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts_noWeight['signal_jesup_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts_noWeight['signal_jesdn_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['signal_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['signal_jesup_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['signal_jesdn_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          JESNP['signal_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)]])

            # nominal_incl += evts['qqzz_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            # up_incl += evts['qqzz_jesup_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            # dn_incl += evts['qqzz_jesdn_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            table.append(['qqzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin),
                          evts_noWeight['qqzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts_noWeight['qqzz_jesup_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts_noWeight['qqzz_jesdn_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['qqzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['qqzz_jesup_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['qqzz_jesdn_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          JESNP['qqzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)]])

            # nominal_incl += evts['ggzz_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            # up_incl += evts['ggzz_jesup_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            # dn_incl += evts['ggzz_jesdn_'+jesName+'_'+year+'_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
            table.append(['ggzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin),
                          evts_noWeight['ggzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts_noWeight['ggzz_jesup_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts_noWeight['ggzz_jesdn_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['ggzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['ggzz_jesup_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          evts['ggzz_jesdn_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                          JESNP['ggzz_'+jesName+'_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)]])

        # nominal_incl += evts['ZX_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
        # up_incl += evts['ZX_jesup_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
        # dn_incl += evts['ZX_jesdn_'+fState+'_'+obsname_out+'_recobin'+str(recobin)]
        table.append(['ZX_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin),
                      evts_noWeight['ZX_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                      evts_noWeight['ZX_jesup_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                      evts_noWeight['ZX_jesdn_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                      evts['ZX_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                      evts['ZX_jesup_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                      evts['ZX_jesdn_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)],
                      JESNP['ZX_'+fState+'_'+str(year)+'_'+obsname_out+'_recobin'+str(recobin)]])

        # table.append(['INCLUSIVE', nominal_incl, up_incl, dn_incl, str(round(dn_incl/nominal_incl,3))+'/'+str(round(up_incl/nominal_incl,3))])
        table.append([])

        tables[fState, recobin] = tabulate(table, headers=['bin', 'nominal_evts', 'jes_up_evts', 'jes_dn_evts', 'nominal', 'jes_up', 'jes_dn', 'ratio(dn/up)'], numalign="right", tablefmt="github")

        with open('tables/'+obsname_out+'/JESNP_table_'+obsname_out+'_'+fState+'_'+str(year)+'_recobin'+str(recobin)+'.py', 'w') as f:
            f.write(tables[fState,recobin])
