import sys, os, string, re, pwd, commands, ast, optparse, shlex, time, copy
import numpy as np
from subprocess import *

sys.path.append('./inputs')
sys.path.append('./LHScans')
from plotUtils import *

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--verbose', action='store_true', dest='VERBOSE', default=False, help='print values')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

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

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

def writeline(obs, bins, sigma, delta, stat, syst, doubleDiff):
	fr_dn = str(round(delta[0]/sigma * 100, 0))
	fr_up = str(round(delta[1]/sigma * 100, 0))

	sigma = str(round(sigma, 2))

	if not doubleDiff:
		dn = str(round(bins[0], 2))
		up = str(round(bins[1], 2))
	else:
		obs_1 = obs[0]; obs_2 = obs[1]
		dn_1 = str(round(bins[0][0], 2))
		up_1 = str(round(bins[0][1], 2))
		dn_2 = str(round(bins[1][0], 2))
		up_2 = str(round(bins[1][1], 2))

	e_dn = str(round(delta[0], 2))
	e_up = str(round(delta[1], 2))

	st_dn = str(round(stat[0], 2))
	st_up = str(round(stat[1], 2))

	sy_dn = str(round(syst[0], 2))
	sy_up = str(round(syst[1], 2))

	if not doubleDiff:
		line = "$%s < %s < %s$ & $%s$ & $_{%s}^{+%s}$ & $_{%s}^{+%s}$ & $_{%s}^{+%s}$ & $_{-%s}^{+%s}$" %(dn, obs, up, sigma, e_dn, e_up, fr_dn, fr_up, st_dn, st_up, sy_dn, sy_dn)
	else:
		# line = "Bin %s & $%s$ & $_{%s}^{+%s}$ & $_{%s}^{+%s}$ & $_{%s}^{+%s}$ & $_{-%s}^{+%s}$" %(dn, sigma, e_dn, e_up, fr_dn, fr_up, st_dn, st_up, sy_dn, sy_dn)
		line = "$%s < %s < %s$, $%s < %s < %s$ & $%s$ & $_{%s}^{+%s}$ & $_{%s}^{+%s}$ & $_{%s}^{+%s}$ & $_{-%s}^{+%s}$" %(dn_1, obs_1, up_1, dn_2, obs_2, up_2, sigma, e_dn, e_up, fr_dn, fr_up, st_dn, st_up, sy_dn, sy_dn)
	return line

def table(results, obsLabel, doubleDiff):
	print("\\begin{table}[hb]")
	print("\\centering")
	print("\\renewcommand{\\arraystretch}{1.5}")
	print("\\begin{tabular}{ |c||c|c|c|c|c| }") #{ |p{3cm}||p{3cm}|p{3cm}|p{3cm}|p{3cm}|  }")
	print("\\hline")
        if not doubleDiff:
		print("\\multicolumn{6}{|c|}{HZZ4L $\sigma_{\\text{fid}}$, $%s$} \\\\[5pt]" %obsLabel)
	else:
		print("\\multicolumn{6}{|c|}{HZZ4L $\sigma_{\\text{fid}}$, $%s$ vs $%s$} \\\\[5pt]" %(obsLabel[0], obsLabel[1]))
	print("\\hline")
	print("Bin range & $d\sigma_{fid}$ (fb) & $\Delta_{tot}$ (fb) & $\delta_{rel}$ (\% fb) & $\Delta_{stat}$ (fb) & $\Delta_{syst}$ (fb)\\\\[5pt]")
	print("\\hline")

	for line in results:
		print("%s \\\\[5pt]" %line)
		print("\\hline")
	print("\\end{tabular}")
	print("\\end{table}")

obsName = opt.OBSNAME
obsLabel = obsName
obsName, doubleDiff = getObsName(obsName)
if doubleDiff:
	obsFirst = getMath(obsLabel.split(' vs ')[0])
	obsSecond = getMath(obsLabel.split(' vs ')[1])

obsBins   = parseBins(opt.OBSBINS)

if doubleDiff: bound_1st, bound_2nd = binBoundaries(obsBins)

label = getMath(obsName)
if doubleDiff: label = [obsFirst, obsSecond]

_temp = __import__('resultsXS_LHScan_expected_'+obsName+'_v3', globals(), locals(), ['resultsXS'], -1)
resultsXS = _temp.resultsXS

lines = []

nBins=len(obsBins)
if not doubleDiff: nBins = nBins -1
obs_reco_high = obsBins[1]

for _bin in range(nBins):
	if not doubleDiff:
		obs_reco_low = obsBins[_bin]
		obs_reco_high = obsBins[_bin+1]
	else:
		down = bound_1st[_bin]
		up = bound_2nd[_bin]
		
	binName = 'SM_125_%s_genbin%i' %(obsName, _bin)
	e_dn = resultsXS[binName]['uncerDn']
	e_up = resultsXS[binName]['uncerUp']
	sigma = resultsXS[binName]['central']

	binName = 'SM_125_%s_genbin%i_statOnly' %(obsName, _bin)
	st_dn = resultsXS[binName]['uncerDn']
	st_up = resultsXS[binName]['uncerUp']
	
	sy_dn = np.sqrt(e_dn**2 - st_dn**2)
	sy_up = np.sqrt(e_up**2 - st_up**2)

	if not doubleDiff:
		lines.append(writeline(label, [obs_reco_low, obs_reco_high], sigma, [e_dn, e_up], [st_dn, st_up], [sy_dn, sy_up], doubleDiff))
	else:
		lines.append(writeline([obsFirst, obsSecond], [down, up], sigma, [e_dn, e_up], [st_dn, st_up], [sy_dn, sy_up], doubleDiff))

table(lines, label, doubleDiff)
