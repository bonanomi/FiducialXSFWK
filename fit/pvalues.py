import uproot
import binning
import latex_names as tex
import os
from itertools import product
from scipy.stats import chi2

th = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs4l_br','unc_qcd','unc_pdf','unc_acc','unc_br'])

HCOMB_NAMES = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta4p7': 'NJ'}
DECAY_OBS = ['D0m','Dcp','D0hp','Dint','DL1','DL1Zg','costhetaZ1','costhetaZ2','costhetastar','phi','phistar','massZ1','massZ2']
CHANNELS = ['4l', '2e2mu']



FITS_PATH = '../combine_files'
PVAL_PATH = '../pvalues'

def checkDir(folder_path):
    isdir = os.path.isdir(folder_path)
    if not isdir:
        print('Directory {} does not exist. Creating it.' .format(folder_path))
        os.mkdir(folder_path)

class pvalue():
    def __init__(self, obs, ws):
        observable = Observable(obs)
        self.obs_name = obs
        self.ndof = observable.nr_bins-1

        self.ws_name = ws

        self.nll = 0.0
        self.pval = 0.0

        self.chi2pdf = chi2(self.ndof)

        self.get_nll()
        self.set_pvalue()

    def get_nll(self):
        fname = f'{PVAL_PATH}/higgsCombine{self.ws_name}.MultiDimFit.mH125.38.root'#DataSMCompat_{self.obs_name}.MultiDimFit.mH125.38.root'
        if self.obs_name in HCOMB_NAMES: fname = fname.replace(self.obs_name, HCOMB_NAMES[self.obs_name])
        self.nll = uproot.open(fname)['limit'].arrays()
        self.nll = self.nll[b'deltaNLL'][1]
        
    def set_pvalue(self):
        _cdf = self.chi2pdf.cdf(2*self.nll)
        self.pval = 1 - _cdf

class Observable():
    def __init__(self, obs):
        self.bins = []
        self.bins_centers = []
        self.dsigmas = []
        self.acceptance = {}
        self.get_bins(obs)

        self.nr_bins = int(len(self.bins))

        self.get_bins_centers()
        self.set_acceptance(obs)
        
    def get_bins(self, obs):
        self.bins = binning.binning_v2(obs)
        
    def nr_bins(self):
        self.n_bins = len(self.bins)
        
    def get_bins_centers(self):
        bins_centers = []
        dsigmas = []
        for i in range(self.nr_bins-1):
            bin_center = (self.bins[i+1] + self.bins[i])*0.5
            dsigma = self.bins[i+1] - self.bins[i]
            
            if dsigma>100: 
                dsigma = 50
                bin_center = 250
                
            bins_centers.append(bin_center)
            dsigmas.append(dsigma)

        self.bins_centers = bins_centers
        self.dsigmas = dsigmas
        
    def set_acceptance(self, obs):
        acc = __import__(f'inputs_sig_extrap_{obs}_Full', globals(), locals(), ['acc'])
        self.acceptance = acc.acc
        
class XSEC():
    def __init__(self, obs, mh):
        observable = Observable(obs)
        self.obs_name = obs
        self.acceptance = observable.acceptance
        self.nr_bins = observable.nr_bins-1
        
        self.mh = mh
        
        self.xsec_sm = {}
        self.xsec_mh = {}
        self.xsec_4l = {}
        self.xsec_2e2mu = {}

        self.set_xsec()
        self.set_xsec_fs()
        self.xsec_fs = self.xsec_4l
        self.xsec_fs.update(self.xsec_2e2mu)

    def set_xsec(self):
        tmp_xs = {}
        tmp_xs_sm = {}
        h_mass = self.mh
        obs = self.obs_name
        for channel in ['4e','4mu','2e2mu']:
            for obsBin in range(self.nr_bins):
                fidxs_sm = 0
                fidxs_sm += th.higgs_xs['ggH_'+'125.0']*\
                            th.higgs4l_br['125.0'+'_'+channel]*\
                            self.acceptance['ggH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += th.higgs_xs['VBF_'+'125.0']*\
                            th.higgs4l_br['125.0'+'_'+channel]*\
                            self.acceptance['VBFH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += th.higgs_xs['WH_'+'125.0']*\
                            th.higgs4l_br['125.0'+'_'+channel]*\
                            self.acceptance['WH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += th.higgs_xs['ZH_'+'125.0']*\
                            th.higgs4l_br['125.0'+'_'+channel]*\
                            self.acceptance['ZH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs_sm += th.higgs_xs['ttH_'+'125.0']*\
                            th.higgs4l_br['125.0'+'_'+channel]*\
                            self.acceptance['ttH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                fidxs = 0
                fidxs += th.higgs_xs['ggH_'+h_mass]*\
                         th.higgs4l_br[h_mass+'_'+channel]*\
                         self.acceptance['ggH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += th.higgs_xs['VBF_'+h_mass]*\
                         th.higgs4l_br[h_mass+'_'+channel]*\
                         self.acceptance['VBFH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += th.higgs_xs['WH_'+h_mass]*\
                         th.higgs4l_br[h_mass+'_'+channel]*\
                         self.acceptance['WH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += th.higgs_xs['ZH_'+h_mass]*\
                         th.higgs4l_br[h_mass+'_'+channel]*\
                         self.acceptance['ZH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
                fidxs += th.higgs_xs['ttH_'+h_mass]*\
                         th.higgs4l_br[h_mass+'_'+channel]*\
                         self.acceptance['ttH125_'+channel+'_'+obs+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

                tmp_xs_sm[channel+'_genbin'+str(obsBin)] = fidxs_sm
                tmp_xs[channel+'_genbin'+str(obsBin)] = fidxs        
        
        self.xsec_sm = tmp_xs_sm
        self.xsec_mh = tmp_xs

    def set_xsec_fs(self):
        for obsBin in range(self.nr_bins):
            self.xsec_4l['r4lBin'+str(obsBin)] = self.xsec_mh['4e_genbin'+str(obsBin)]+\
                                                 self.xsec_mh['4mu_genbin'+str(obsBin)]
            self.xsec_2e2mu['r2e2muBin'+str(obsBin)] = self.xsec_mh['2e2mu_genbin'+str(obsBin)]
        
class combineCommand():
    def __init__(self, obs, SM: bool, version: str, channel: str):
        observable = Observable(obs)
        self.obs_name = obs
        self.acceptance = observable.acceptance
        self.bins = observable.bins
        self.nr_bins = observable.nr_bins-1
        
        xsec = XSEC(obs, '125.38')
        self.tmp_xs = xsec.xsec_mh
        self.tmp_xs_sm = xsec.xsec_sm
        self.xsec_fs = xsec.xsec_fs

        self.version = version
        self.channel = channel
        
        self.pois  = []
        self.pois_to_float = []
        self.pois_to_freeze = []
        self.kappas = {}
        self.ws = ''
        self.output = ''
        self.command = ''
        self.algo = ''
        
        self.get_ws(obs, SM)
        self.get_pois(obs)
        self.get_kappas(obs)
        self.set_algo(SM)
        self.init_command()
        self.define_pois()
        if SM: self.set_ranges()
        self.set_params(SM)
        if self.version=='v3': self.set_kappas()
        self.freezeParams()
        self.replacement_map()
        
    def get_ws(self, obs, SM):
        if SM:
            self.ws = f'{FITS_PATH}/SM_125_all_13TeV_xs_{self.obs_name}_bin_{self.version}.root'
            self.output = f'bestFit_{self.obs_name}_{self.version}_{self.channel}'
        else:
            self.ws = f'higgsCombinebestFit_{self.obs_name}_{self.version}_{self.channel}'
            self.ws += '.MultiDimFit.mH125.38.root'
            self.output = f'DataSMCompat_{self.obs_name}_{self.version}_{self.channel}'
            
    def set_algo(self, SM):
        if SM:
            self.algo = 'singles'
        else:
            self.algo = 'fixed'
    
    def get_pois(self, obs):
        if self.version=='v3':
            self.pois = [f'r_smH_{obs}_{j}' for j in range(self.nr_bins)]
        elif (self.version=='v4') & (self.channel=='2e2mu'):
            self.pois = [f'r2e2muBin{j}' for j in range(self.nr_bins)]
            self.pois_to_float = self.pois.copy()
            self.pois_to_freeze = [f'r4lBin{j}' for j in range(self.nr_bins)]
            self.pois.extend(self.pois_to_freeze)
        elif (self.version=='v4') & (self.channel=='4l'):
            self.pois = [f'r4lBin{j}' for j in range(self.nr_bins)]
            self.pois_to_float = self.pois.copy()
            self.pois_to_freeze = [f'r2e2muBin{j}' for j in range(self.nr_bins)]
            self.pois.extend(self.pois_to_freeze)
        else:
            raise Exception(f'PhysicsModel {self.version} is not supported! Only v3 and v4 are.')

    def get_kappas(self, obs):
        kappas = [f'K{k}Bin{j}' for k,j in product([1,2],range(self.nr_bins))]
        tmp_xs = self.tmp_xs
        tmp_xs_sm = self.tmp_xs_sm
                
        for obsBin in range(self.nr_bins):
            fidxs4e = tmp_xs['4e_genbin'+str(obsBin)]
            fidxs4mu = tmp_xs['4mu_genbin'+str(obsBin)]
            fidxs2e2mu = tmp_xs['2e2mu_genbin'+str(obsBin)]
            frac4e = fidxs4e/(fidxs4e+fidxs4mu+fidxs2e2mu)
            frac4mu = fidxs4mu/(fidxs4e+fidxs4mu+fidxs2e2mu)
            fidxs4e_sm = tmp_xs_sm['4e_genbin'+str(obsBin)]
            fidxs4mu_sm = tmp_xs_sm['4mu_genbin'+str(obsBin)]
            fidxs2e2mu_sm = tmp_xs_sm['2e2mu_genbin'+str(obsBin)]
            frac4e_sm = fidxs4e_sm/(fidxs4e_sm+fidxs4mu_sm+fidxs2e2mu_sm)
            frac4mu_sm = fidxs4mu_sm/(fidxs4e_sm+fidxs4mu_sm+fidxs2e2mu_sm)
            K1 = frac4e/frac4e_sm
            K2 = frac4mu/frac4mu_sm * (1.0-frac4e_sm)/(1.0-frac4e)
            
            self.kappas[f'K1Bin{obsBin}'] = K1
            self.kappas[f'K2Bin{obsBin}'] = K2
                
    def init_command(self):
        self.command = f'combine -M MultiDimFit {self.ws} --algo={self.algo} -n {self.output} '
        self.command += '-m 125.38 --saveWorkspace'
        
    def define_pois(self):
        self.command += ' '
        self.command += '--redefineSignalPOIs '
        _pois = self.pois if self.version=='v3' else self.pois_to_float

        for poi in _pois:
            self.command += f'{poi},'
        self.command = self.command[:-1]
    
    def set_ranges(self):
        self.command += ' '
        self.command += '--setParameterRanges '
        for poi in self.pois:
            self.command += f'{poi}=0,5:'
        self.command = self.command[:-1]
        
    def set_params(self, SM):
        self.command += ' '
        if SM:
            self.command += '--setParameters '
        else:
            self.command += '--snapshotName "MultiDimFit" --saveNLL --skipInitialFit ' 
            self.command += '--X-rtd MINIMIZER_freezeDisassociatedParams --fixedPointPOIs '
        for poi in self.pois:
            if self.version=='v3':
                self.command += f'{poi}=1,'
            elif self.version=='v4':
                self.command += f'{poi}={self.xsec_fs[poi]},'
        self.command = self.command[:-1]

    def set_kappas(self):
        self.command +=',' #reintroduce comma
        for kappa in self.kappas:
            self.command += f'{kappa}={self.kappas[kappa]},'
        self.command = self.command[:-1]
        
    def freezeParams(self):
        self.command += ' --floatOtherPOIs 0 --freezeParameters MH,'
        if self.version == 'v4':
            for poi in self.pois_to_freeze:
                self.command += f'{poi},'
        self.command = self.command[:-1]

    def replacement_map(self):
        if self.obs_name in HCOMB_NAMES:
            self.command = self.command.replace(self.obs_name, HCOMB_NAMES[self.obs_name])

class Table():
    def __init__(self, _file):
        self.file = _file
    
    def write_header(self):
        self.file.write('\\begin{table*}[hb] \n')
        self.file.write('\\centering \n')
        self.file.write('\\begin{tabular}{cccc} \n')
        self.file.write('\\hline \n')
        self.file.write('\\textbf{Observable} & \\textbf{Model} &  ')
        self.file.write('\\textbf{Final State} & \\textbf{p-value} \\\\ \n')
        self.file.write('\\hline \n')
        self.file.write('\\hline \n')

    
    def write_line(self, _obs, _ver, _fs, _pval):
        _pval = str(_pval)
        line = _obs+' & '+_ver+' & '+_fs+' & '+_pval+'\\\\'
        line = line.replace(_obs, tex.latex_names[_obs])
        line = line.replace(_fs, tex.latex_names[_fs])
        self.file.write(line)
        self.file.write('\n')
        

    def write_footer(self):
        self.file.write('\\hline \n')
        self.file.write('\\end{tabular} \n')
        self.file.write('\\end{table*}')

def fill(_file, _obs, _commands):
    _file.write(_obs+'\n')
    for _cmd in _commands:
      _file.write(_cmd+'\n')
    _file.write('\n\n\n')


if __name__ == '__main__':

  dump_file = open('pval_cmds.txt','w')

  for ver, ch, _obs in product(['v3', 'v4'], CHANNELS, binning.BINS):
    if ((ver=='v3') & (ch=='2e2mu')): continue
    if ((ver=='v4') & (_obs not in DECAY_OBS)): continue
    if 'mass4l' in _obs: continue
    if 'vs' in _obs: continue
    if 'kL' in _obs: continue
    bestFit = combineCommand(_obs, True, ver, ch)
    compatibilitySM = combineCommand(_obs, False, ver, ch)
    commands = [bestFit.command, compatibilitySM.command]
    fill(dump_file, _obs, commands)
    fname = f"higgsCombine{compatibilitySM.output}.MultiDimFit.mH125.38.root"
    if _obs in HCOMB_NAMES: fname = fname.replace(_obs, HCOMB_NAMES[_obs])
    if (os.path.isfile(fname) | os.path.isfile(f"{PVAL_PATH}/{fname}")): 
     print(f"................. Skip {_obs}, fit already done")
     continue

    os.system(bestFit.command)
    os.system(compatibilitySM.command)

  dump_file.close()

  print("\n\n\n")
  checkDir(f"{PVAL_PATH}")
  os.system(f"mv higgsCombine* {PVAL_PATH}")

  latex_table = open('pval_table.txt', 'w')
  table = Table(latex_table)
  table.write_header()

  for ver, ch, _obs in product(['v3', 'v4'], CHANNELS, binning.BINS):
    if ((ver=='v3') & (ch=='2e2mu')): continue
    if ((ver=='v4') & (_obs not in DECAY_OBS)): continue
    if 'mass4l' in _obs: continue
    if 'vs' in _obs: continue
    if 'kL' in _obs: continue
    compatibilitySM = combineCommand(_obs, False, ver, ch)
    combine_file = compatibilitySM.output
    pval = round(pvalue(_obs,combine_file).pval,2)
    print(f"{_obs} ({ver}) SM compatibility with p-val ({ch}) = {pval}")
    table.write_line(_obs,ver,ch,pval)
  table.write_footer()

  latex_table.close()
