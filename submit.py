import sys, os, pwd, commands
from tqdm import tqdm
from subprocess import *
import optparse, shlex, re

sys.path.append("inputs/")
from binning import binning

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Full',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--submit', action='store_true', dest='SUBMIT', default=False, help='Do only scripts')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--fitOnly', action='store_true', dest='FITONLY', default=False, help='Run only fit')
    parser.add_option('',   '--myID', dest='MYID', type='string', default='', help='Your CERN ID')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()


# Define function for processing of os command
def processCmd(cmd, quiet = 0):
    output = "\n"
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=-1)
    for line in iter(p.stdout.readline, ""):
        output=output+str(line)
        print line,
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError("%r failed, exit status: %d" % (cmd, p.returncode))
    return output

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()

year    = str(opt.YEAR)
unblind = opt.UNBLIND
fitOnly = opt.FITONLY
submit = opt.SUBMIT
my_id   = opt.MYID


def prepareScript(obs_name, obs_bin, name):
    scriptName = './scripts/script_'+name+'.sh'
    with open(scriptName, 'w') as f:
        f.write('export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n')
        f.write('source $VO_CMS_SW_DIR/cmsset_default.sh \n')
        f.write('echo $MYPSW | /opt/exp_soft/cms/t3/eos-login -keytab -init -username ' + my_id + '\n')
        f.write('/opt/exp_soft/cms/t3/eos-login -init -username ' + my_id + '\n')
        f.write('cmsenv \n')
        
        pipeline = 'python pipeline.py --obsName "' + obs_name  + '" --obsBins "' + obs_bin + '" --year "' + year + '"'
        if unblind:
            pipeline += " --unblind True"
        if fitOnly:
            pipeline += " --fitOnly True"     

        f.write(pipeline)

    processCmd('chmod 0744 ' + scriptName)

    return scriptName

for key in binning:
     _bin = binning[key]
     _obs_name = key
     _screen_name = key
 
     if ' vs ' in key:
         first = key.split(' vs ')[0]
         second = key.split(' vs ')[1]
         _screen_name = first + '_' + second
     
     _scriptName = prepareScript(_obs_name, _bin, _screen_name)     
     
     if submit: 
         cmd = "screen -dmS "+ _screen_name + " '" + _scriptName + "'" 
     
         print(cmd)
         processCmd(cmd)
