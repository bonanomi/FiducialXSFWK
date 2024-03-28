# python3 nano_skimmer.py
# ./nano_skimmer.py
# Requires awkward > 1.4

import uproot
import os
import numpy as np
import awkward as ak
from itertools import product
from collections import defaultdict

branches_map = {
    "event": "EventNumber",
    "ZZCand_mass": "ZZMass",
    "ZZCand_Z1mass": "Z1Mass",
    "ZZCand_Z2mass": "Z2Mass",
    "ZZCand_Z1flav": "Z1Flav",
    "ZZCand_Z2flav": "Z2Flav",
    "puWeight": "PUWeight",
    "overallEventWeight": "overallEventWeight",
    "ZZCand_dataMCWeight": "dataMCWeight",
    "Generator_weight": "genHEPMCweight",
}
b_qqzz_map = {
    # "KFactor_EW_qqZZ_Weight": "KFactor_EW_qqZZ",
    "KFactor_QCD_qqZZ_M_Weight": "KFactor_QCD_qqZZ_M",
}
b_ggzz_map = {
    "KFactor_QCD_ggZZ_Nominal_Weight": "KFactor_QCD_ggZZ_Nominal"
}

def get_genEventSumw(fname):
    with uproot.open(f"{fname}") as f:
        genEventSumw = f['Runs/genEventSumw'].array(library="np")
    return sum(genEventSumw)

class Skimmer:
    '''
       Skimmer class to convert nanoAOD ntuples into miniAOD ones.
       To initialize a `Skimmer` object one needs to specify:
       `process`: name of the process, typically as in the input fname
       `h_mass`: H boson mass, typically as in the input fname
       `year`: data-taking period
       `data_type`: Data or MC
       The input fname expected is therefore:
       `fname = f"{self.process}{self.h_mass}_reducedTree_{self.data_type}_{self.year}.root"`
       The `Skimmer.skim` method creates a file named:
       `fname = f"{self.process}{self.h_mass}_reducedTree_{self.data_type}_{self.year}_skimmed.root"`
    '''
    def __init__(self,
                 process: str,
                 year: str,
                 data_type: str,
                 _dir: str):
        self.process = process
        self.year = year
        self.data_type = data_type
        self._dir = _dir

        self.fname = self._set_fname
        self.out_name = self._set_out_name
        self.counter = 0

        if self._dir != "":
            self.fname = f"{self._dir}/{self.process}/{self.fname}"
            self.out_name = f"{self._dir}/{self.process}/{self.out_name}"
        
        self.b_map = branches_map.copy()
        if "ggTo" in self.process:
            self.b_map.update(b_ggzz_map)
        if "ZZTo" in self.process:
            self.b_map.update(b_qqzz_map)

        if self.year == "2022EE":
            if "ggTo" in self.process:
                cnt_fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/240313_ggZZ/MC_2022EE/{self.process}/ZZ4lAnalysis.root"
            else:
                cnt_fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/240319/MC_2022EE/{self.process}/ZZ4lAnalysis.root"
        if self.year == "2022":
            cnt_fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/240321_NanoMC2022/{self.process}/ZZ4lAnalysis.root"
        self.counter = get_genEventSumw(cnt_fname)

        print(f"... Set up skimmer for {self.fname}")
        print(f"... Will create skimmed file: {self.out_name}\n")
        
        self.b_to_read = list(self.b_map.keys())
        self.b_to_dump = list(self.b_map.values())
        self.b_zz_cands = [b for b in self.b_to_read if 'ZZCand' in b]
        self.b_zz_cands = self.b_zz_cands + [b for b in self.b_to_read if 'index' in b]
        
        b_log = "\n".join("{0:30} {1}".format(r, d) for r, d in zip(self.b_to_read, self.b_to_dump))
        print("{0:30} {1}".format("Reading from:", "Dumping to:"))
        print(f"{b_log}\n")
        
    @property
    def _set_fname(self):
        '''
           Util function that returns the name of the input file,
           i.e. CJLST output from nanoAOD framework.
        '''
        fname = f"{self.process}_reducedTree_{self.data_type}_{self.year}.root"
        return fname

    @property
    def _set_out_name(self):
        '''
           Util function that returns the name of the output file,
           i.e. the file with Run2-like TTrees and branches names.
        '''
        fname = f"{self.process}_reducedTree_{self.data_type}_{self.year}_skimmed.root"
        return fname

    def _get_branches(self, b_in, tree, drop_zzcands=False):
        '''
           Util function that reads from the nanoAOD file
           the branches specified in `branches_map`.
           For failed events the brances with the ZZ candidates
           are skipped.
        '''
        if tree == "AllEvents" and drop_zzcands:
            b_read = list(set(self.b_to_read).difference(self.b_zz_cands))
        else:
            b_read = b_in

        with uproot.open(f"{self.fname}") as f:
            branches = f[tree].arrays(b_read)
            zz_idx = f[tree].arrays("bestCandIdx")

        sel = zz_idx['bestCandIdx']!=-1
        branches = branches[sel]
        return branches

    def _set_branches_and_types(self, tree):
        '''
           Util function that reads ak.Arrays
           and returns dictionaries with the content
           of the arrays (`d_vals`) and their types (`d_types`).
           This output of this function is fundamental to use
           `uproot.mktree` and `.extend` in `self.skim`.
        '''
        d_vals = defaultdict(list)
        d_types = defaultdict(list)
        
        branches = self._get_branches(self.b_to_read, tree, True)

        for b_read, b_write in zip(self.b_to_read, self.b_to_dump):
            if tree == "AllEvents" and b_read in self.b_zz_cands: continue
            d_types[b_write] = branches[b_read].type
            d_vals[b_write]  = branches[b_read]

        return dict(d_vals), dict(d_types)

    def skim(self):
        '''
           Skimming function that reads the nanoAOD input format
           and create an output file consistent with the miniAOD format.
           The output file (`_skimmed.root`) contains two `TTree`s:
           `candTree` and `candTree_failed`, as expected by the analysis framework.
        '''
        print("+++ SKIMMING STARTED! +++\n")

        d_vals_pass, d_types_pass = self._set_branches_and_types("Events")

        d_types_pass["Counter"] = ak.Array(np.ones(len(d_vals_pass["EventNumber"]))).type
        d_vals_pass["Counter"] = ak.Array(np.ones(len(d_vals_pass["EventNumber"])))*self.counter
        
        with uproot.recreate(f"{self.out_name}") as fout:
            fout.mktree("candTree", d_types_pass)
            fout["candTree"].extend(d_vals_pass)

        print("+++ SKIMMING COMPLETED! +++ \n")

if __name__ == "__main__":
    for period in ["2022", "2022EE"]:
        cjlst_dir = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/nanoProd_Run3_{period}/cjlst_trees/"
        for process in ["ZZTo4l", "ggTo2e2mu_Contin_MCFM701", "ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701"]:
            skimmer = Skimmer(process, period, "MC", cjlst_dir)
            if os.path.exists(skimmer.out_name):
                 os.system(f"rm {skimmer.out_name}")
            skimmer.skim()
    