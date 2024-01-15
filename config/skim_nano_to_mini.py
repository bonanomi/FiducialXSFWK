import uproot
import numpy as np
import awkward as ak
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
    "FidDressedLeps_pt": "GENlep_pt",
    "FidDressedLeps_eta": "GENlep_eta",
    "FidDressedLeps_phi": "GENlep_phi",
    "FidDressedLeps_mass": "GENlep_mass",
    "FidDressedLeps_id": "GENlep_id",
    "FidDressedLeps_momid": "GENlep_MomId",
    "FidDressedLeps_mommomid": "GENlep_MomMomId",
    "FidDressedLeps_RelIso": "GENlep_RelIso",
    "FidZZ_mass": "GENmass4l",
    "FidZZ_pt": "GENpT4l",
    "FidZZ_eta": "GENeta4l",
    "FidZZ_phi": "GENphi4l",
    "FidZZ_rapidity": "GENrapidity4l",
    "FidZ_DauPdgId": "GENZ_DaughtersId",
    "FidZ_MomPdgId": "GENZ_MomId",
    "ggH_NNLOPS_Weight": "ggH_NNLOPS_weight",
}
branches_to_array = {"GENlep_Hindex": ["FidZZ_Z1l1Idx","FidZZ_Z1l2Idx","FidZZ_Z2l1Idx","FidZZ_Z2l2Idx"]}

class Skimmer:
    def __init__(self,
                 process: str,
                 h_mass: str,
                 year: str,
                 data_type: str):
        self.process = process
        self.h_mass = h_mass
        self.year = year
        self.data_type = data_type
        
        self.fname = self._set_fname
        self.out_name = self._set_out_name
        
        print(f"Set up skimmer for {self.fname}")
        print(f"Will create skimmed file: {self.out_name}")
        
        self.b_to_read = list(branches_map.keys())
        self.b_to_dump = list(branches_map.values())
        
        print(f"Will read the following branches: {self.b_to_read}")
        print(f"Will dump them into: {self.b_to_dump}")

        self.b_array_to_dump = list(branches_to_array.keys())
        self.b_array_to_read = list(branches_to_array.values())

        print(f"Will read the following branches: {self.b_array_to_read}")
        print(f"Will dump them into: {self.b_array_to_dump}")

        
    @property
    def _set_fname(self):
        fname = f"{self.process}{self.h_mass}_reducedTree_{self.data_type}_{self.year}.root"
        return fname

    @property
    def _set_out_name(self):
        fname = f"{self.process}{self.h_mass}_reducedTree_{self.data_type}_{self.year}_skimmed.root"
        return fname
    
    @property
    def _get_events(self):
        events = np.random.randint(low=0, high=13158, size=5000)
        
#         with uproot.open(f"{self.fname}") as f:
#             events = f["Events/event"].array(library = "np")

        return events
    
    def _get_branches(self, b_in, tree):
        with uproot.open(f"{self.fname}") as f:
            branches = f[tree].arrays(b_in)
            
        return branches

    def _branches_to_array(self, b_in):
        inputs = [b_in[i] for i in b_in.fields]
        single_branch = ak.concatenate([ak.singletons(arr) for arr in inputs], axis=1)

        return single_branch

    def _set_branches_and_types(self, tree):
        d_vals = defaultdict(list)
        d_types = defaultdict(list)
        
        branches = self._get_branches(self.b_to_read, tree)

        for b_read, b_write in zip(self.b_to_read, self.b_to_dump):
            d_types[b_write] = branches[b_read].type
            d_vals[b_write]  = branches[b_read]
            
        for b_read, b_write in zip(self.b_array_to_read, self.b_array_to_dump):
            branches_array = self._get_branches(b_read, tree)
            b_array = self._branches_to_array(branches_array)
            d_types[b_write] = b_array.type
            d_vals[b_write]  = b_array

        return dict(d_vals), dict(d_types)
    
    def _failed_events(self, d_vals, d_types):
        pass_events = self._get_events
        sel = ~ak.Array([x in np.array(pass_events) for x in np.array(d_vals['EventNumber'])])
        
        for d in d_vals:
            d_vals[d] = d_vals[d][sel]
            d_types[d] = d_vals[d].type

        return dict(d_vals), dict(d_types)

    def skim(self):
        d_vals_pass, d_types_pass = self._set_branches_and_types("Events")
        d_vals_fail, d_types_fail = self._set_branches_and_types("AllEvents")
        
        d_vals_fail, d_types_fail = self._failed_events(d_vals_fail, d_types_fail)
        
        with uproot.recreate(f"{self.out_name}") as fout:
            fout.mktree("candTree", d_types_pass)
            fout["candTree"].extend(d_vals_pass)

            fout.mktree("candTree_fail", d_types_fail)
            fout["candTree_fail"].extend(d_vals_fail)

if __name__ == "__main__":
    skimmer = Skimmer("ggH", "125", "2018", "MC")
    skimmer.skim()
