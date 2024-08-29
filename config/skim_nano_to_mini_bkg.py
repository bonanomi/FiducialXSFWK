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
    "ZZCand_pt": "ZZPt",
    "ZZCand_rapidity": "ZZy",
    "rapidity4lAbs": "ZZyAbs",
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
        # TODO: Putback!!!
        if "ggTo" in self.process:
            self.b_map.update(b_ggzz_map)
        if "ZZTo" in self.process:
            self.b_map.update(b_qqzz_map)

        cnt_fname = f"{self._dir}/{self.process}/ZZ4lAnalysis.root"
        '''
        if self.year == "2022EE":
            if "ggTo" in self.process:
                cnt_fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/240313_ggZZ/MC_2022EE/{self.process}/ZZ4lAnalysis.root"
            else:
                cnt_fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/240319/MC_2022EE/{self.process}/ZZ4lAnalysis.root"
        if self.year == "2022":
            cnt_fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/240321_NanoMC2022/{self.process}/ZZ4lAnalysis.root"
        '''
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
    
    @property
    def _get_events(self):
        '''
           Util function that gets the list of events passing
           the ZZ reconstruction (i.e. events in the `Events` TTree).
        '''        
        with uproot.open(f"{self.fname}") as f:
            events = f["Events/event"].array(library = "np")

        return events
    
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

    def _branches_to_array(self, b_in):
        '''
           Util function that takes as input a series of arrays
           and returns a single branch with singletons.
           This is relevant for all those branches that in nanoAOD
           are stored on an event basis while in miniAOD were stored
           at particle level.
           For example, in nanoAOD we have:
           `FidZZ_Z1l1Idx = [[0, 1, 2, 3], [4, 5, 6, 7], ...]`
           `FidZZ_Z1l2Idx = [[0, 1, 2, 3], [4, 5, 6, 7], ...]`
           `FidZZ_Z2l1Idx = [[0, 1, 2, 3], [4, 5, 6, 7], ...]`
           `FidZZ_Z2l2Idx = [[0, 1, 2, 3], [4, 5, 6, 7], ...]`
           while the structure of the `GENlep_Hindex` in miniAOD is:
           `GENlep_Hindex = [[0, 0, 0, 0], [1, 1, 1, 1], ...] `
           The branches parsed in this way are the ones specified in
           the `branches_to_array` dictionary.
        '''
        inputs = [b_in[i] for i in b_in.fields]
        single_branch = ak.concatenate([ak.singletons(arr) for arr in inputs], axis=1)

        return single_branch

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
    
    def _failed_events(self, d_vals, d_types):
        '''
           Util function that removes from `AllEvents`
           the events that are already stored in `Events` since
           they pass the ZZ candidate selection.
        '''
        pass_events = self._get_events
        print(f"... Will remove {len(pass_events)} events (passing ZZ selection) ")
        print(f"... from the {len(d_vals['EventNumber'])} total events")
        sel = ~ak.Array([x in np.array(pass_events) for x in np.array(d_vals['EventNumber'])])
        
        for d in d_vals:
            d_vals[d] = d_vals[d][sel]
            d_types[d] = d_vals[d].type

        print(f"... After filtering, the cand_failed TTree ")
        print(f"... has {len(d_vals['EventNumber'])} events")

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
    for year in ["2022", "2022EE"]:
        cjlst_dir = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/240820/{year}"
        for process in ["ZZTo4l", "ggTo2e2mu_Contin_MCFM701", "ggTo2e2tau_Contin_MCFM701","ggTo2mu2tau_Contin_MCFM701","ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701","ggTo4tau_Contin_MCFM701"]:
            if not ((year == "2022EE") and (process == "ggTo4mu_Contin_MCFM701")): continue
            skimmer = Skimmer(process, year, "MC", cjlst_dir)
            if os.path.exists(skimmer.out_name):
                 os.system(f"rm {skimmer.out_name}")
            skimmer.skim()