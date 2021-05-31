import numpy
import os
import json
import numpy
import optparse, os, sys
from binning import binning
from tabulate import tabulate

print 'Welcome in RunJES!'

jesNames = ['Total', 'Abs', 'Abs_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavQCD', 'HF', 'HF_year', 'RelBal', 'RelSample_year']

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

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


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

print (obsname, year, obsname_out)
obs_bins, doubleDiff = binning(opt.OBSBINS)

_temp = __import__('JESNP_evts_'+obsname_out+'_'+str(year), globals(), locals(), ['evts'])
evts = _temp.evts
_temp = __import__('JESNP_'+obsname_out+'_'+str(year), globals(), locals(), ['JESNP'])
JESNP = _temp.JESNP

# Check if the folder for tables exist
if not os.path.exists('tables/'+obsname_out):
    os.makedirs('tables/'+obsname_out)

# Tables with numerical values
tables = {}
for fState in ['2e2mu', '4e', '4mu']:
    for jesName in jesNames:
        table = []
        for recobin in range(len(obs_bins)-1):
            nominal_incl = 0
            up_incl = 0
            dn_incl = 0
            for genbin in range(len(obs_bins)-1):
                nominal_incl += evts['fiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)]
                up_incl += evts['fiducial_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)]
                dn_incl += evts['fiducial_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)]
                table.append(['fiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin),
                              evts['fiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)],
                              evts['fiducial_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)],
                              evts['fiducial_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)],
                              JESNP['fiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(genbin)+'_recobin'+str(recobin)]])

            nominal_incl += evts['nonFiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            up_incl += evts['nonFiducial_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            dn_incl += evts['nonFiducial_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            table.append(['nonFiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin),
                          evts['nonFiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['nonFiducial_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['nonFiducial_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          JESNP['nonFiducial_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]])

            nominal_incl += evts['nonResonant_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            up_incl += evts['nonResonant_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            dn_incl += evts['nonResonant_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            table.append(['nonResonant_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin),
                          evts['nonResonant_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['nonResonant_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['nonResonant_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          JESNP['nonResonant_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]])
            table.append(['fid_nonFid_nonRes_INCLUSIVE', nominal_incl, up_incl, dn_incl, str(round(dn_incl/nominal_incl,3))+'/'+str(round(up_incl/nominal_incl,3))])

            nominal_incl += evts['qqzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            up_incl += evts['qqzz_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            dn_incl += evts['qqzz_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            table.append(['qqzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin),
                          evts['qqzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['qqzz_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['qqzz_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          JESNP['qqzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]])

            nominal_incl += evts['ggzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            up_incl += evts['ggzz_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            dn_incl += evts['ggzz_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]
            table.append(['ggzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin),
                          evts['ggzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['ggzz_jesup_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          evts['ggzz_jesdn_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)],
                          JESNP['ggzz_'+jesName+'_'+fState+'_'+obsname_out+'_genbin'+str(recobin)+'_recobin'+str(recobin)]])
            table.append(['INCLUSIVE', nominal_incl, up_incl, dn_incl, str(round(dn_incl/nominal_incl,3))+'/'+str(round(up_incl/nominal_incl,3))])
            table.append([])

        tables[jesName,fState] = tabulate(table, headers=['bin', 'nominal', 'jes_up', 'jes_dn', 'ratio(dn/up)'], numalign="right", tablefmt="github")

        with open('tables/JESNP_table_'+obsname_out+'_'+str(year)+'_'+jesName+'_'+fState+'.py', 'w') as f:
            f.write(tables[jesName,fState])
