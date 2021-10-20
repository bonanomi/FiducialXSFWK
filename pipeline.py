import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import uproot
import math
import ROOT
import ROOT
import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re
import time
from decimal import *
import json

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--fitOnly', action='store_true', dest='FITONLY', default=False, help='Run only fit')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.OBSBINS=='' and opt.OBSNAME!='inclusive'):
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit()


# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

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

def pipeline():
    # prepare the set of bin boundaries to run over, only 1 bin in case of the inclusive measurement
    observableBins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
    ## Run for the given observable
    obsName = str(opt.OBSNAME)
    obsBins = str(opt.OBSBINS)
    year    = str(opt.YEAR)
    unblind = opt.UNBLIND
    fitOnly = opt.FITONLY

    print('============================================================')
    print('============================================================')
    print('=== RUNNING THE WORKFLOW FOR THE FOLLOWING CONFIGURATION ===')
    print('=== %s  %s  %s  %r ===' %(obsName, year, obsBins, unblind))
    print('============================================================')
    print('============================================================')

    if not fitOnly:
        os.chdir('./coefficients')

        cmd = 'python RunCoefficients.py --obsName "'+obsName+'" --obsBins "'+obsBins+'" --year "'+year+'"'
        print cmd
        output = processCmd(cmd)
        print output

        cmd = 'python pdfUncertainty.py --obsName "'+obsName+'" --obsBins "'+obsBins+'" --year "'+year+'"'
        output = processCmd(cmd)
        print output

        os.chdir('../templates')

        cmd = 'python RunTemplates.py --obsName "'+obsName+'" --obsBins "'+obsBins+'" --year "'+year+'"'
        output = processCmd(cmd)
        print output

        os.chdir('../fit')

    if fitOnly: os.chdir('./fit')
    cmd = 'python RunFiducialXS.py --obsName "'+obsName+'" --obsBins "'+obsBins+'" --year "'+year+'"'
    output = processCmd(cmd)
    print output

    if unblind:
        cmd += ' --unblind True'
        output = processCmd(cmd)
        print output

    os.chdir('../LHScans')

    cmd = 'python plot_LLScan.py --obsName "'+obsName+'" --obsBins "'+obsBins+'" --year "'+year+'"'
    output = processCmd(cmd)
    print output

    if unblind:
        cmd += ' --unblind True'
        output = processCmd(cmd)
        print output

    processCmd('roo')
    cmd = 'python RunPlotCoefficients.py --obsName "'+obsName+'" --obsBins "'+obsBins+'" --year "'+year+'"'
    output = processCmd(cmd)
    print output

    cmd = 'copy_to_www.sh '+obsName
    output = processCmd(cmd)
    print output


pipeline()
