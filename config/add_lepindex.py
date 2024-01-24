# python3 add_lepindex.py
# ./add_lepindex.py
# Requires awkward > 1.4
# Requires pyROOT

import ROOT
import uproot
import numpy as np
import awkward as ak
from collections import defaultdict

class Object:
    '''
       Class that allows seeing a set branches plus possibly an index as an Object.
       Class taken from nanoAOD code:
       https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/framework/datamodel.py
    '''

    def __init__(self, event, prefix, index=None):
        self._event = event
        self._prefix = prefix + "_"
        self._index = index

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        if name[:2] == "__" and name[-2:] == "__":
            raise AttributeError
        val = getattr(self._event, self._prefix + name)
        if self._index != None:
            val = val[self._index]
        # convert char to integer number
        val = ord(val) if type(val) == str else val
        self.__dict__[name] = val  # cache
        return val

    def __getitem__(self, attr):
        return self.__getattr__(attr)

    def p4(self, corr_pt=None):
        """Create TLorentzVector for this particle."""
        ret = ROOT.TLorentzVector()
        if corr_pt == None:
            ret.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.mass)
        else:
            ret.SetPtEtaPhiM(corr_pt, self.eta, self.phi, self.mass)
        return ret
    
    def pdgid(self):
        print(self._prefix)
        pdgId = 99
        if self._prefix=='Electron':
            pdgId = 11
        elif self._prefix=='Muon':
            pdgId = 13
        else:
            pdgId = 99
        return pdgId

    def DeltaR(self, other):
        if isinstance(other, ROOT.TLorentzVector):
            deta = abs(other.Eta() - self.eta)
            dphi = abs(other.Phi() - self.phi)
        else:
            deta = abs(other.eta - self.eta)
            dphi = abs(other.phi - self.phi)
        while dphi > math.pi:
            dphi = abs(dphi - 2 * math.pi)
        return math.sqrt(dphi**2 + deta**2)

    def statusflag(self, flag):
        """Find if bit for statusflag is set (for GenPart only)."""
        return (self.statusFlags & statusflags[flag])==statusflags[flag]

    def subObj(self, prefix):
        return Object(self._event, self._prefix + prefix)

    def __repr__(self):
        return ("<%s[%s]>" % (self._prefix[:-1], self._index)) if self._index != None else ("<%s>" % self._prefix[:-1])

    def __str__(self):
        return self.__repr__()


class Collection:
    '''
       Class that creates a Collection of Objects.
       Taken from:
       https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/framework/datamodel.py
    '''
    def __init__(self, event, prefix, lenVar=None):
        self._event = event
        self._prefix = prefix
        if lenVar != None:
            self._len = getattr(event, lenVar)
        else:
            self._len = getattr(event, "n" + prefix)
        self._cache = {}

    def __getitem__(self, index):
        if type(index) == int and index in self._cache:
            return self._cache[index]
        if index >= self._len:
            raise IndexError("Invalid index %r (len is %r) at %s" % (index, self._len, self._prefix))
        elif index < 0:
            raise IndexError("Invalid index %r (negative) at %s" % (index, self._prefix))
        ret = Object(self._event, self._prefix, index=index)
        if type(index) == int:
            self._cache[index] = ret
        return ret

    def __len__(self):
        return self._len

def getLeptons(aCand, event) :
    '''
       Util function that gets the four leptons of a ZZ or ZLL candidate.
       Function taken from:
       https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/tools.py#L33
    '''
    idxs = [aCand.Z1l1Idx, aCand.Z1l2Idx, aCand.Z2l1Idx, aCand.Z2l2Idx]
    electrons = Collection(event, "Electron")
    muons = Collection(event, "Muon")
    leps = list(electrons) + list(muons)
    return [leps[i] for i in idxs]

def get_genEventSumw(input_file, maxEntriesPerSample=None):
    '''
       Util function to get the sum of weights per event.
       Returns the sum of weights, similarly to what we
       stored in Counters->GetBinContent(40) in the miniAODs.
    '''
    f = input_file

    runs  = f.Runs
    event = f.Events
    nRuns = runs.GetEntries()
    nEntries = event.GetEntries()

    iRun = 0
    genEventCount = 0
    genEventSumw = 0.

    while iRun < nRuns and runs.GetEntry(iRun) :
        genEventCount += runs.genEventCount
        genEventSumw += runs.genEventSumw
        iRun +=1
    print ("gen=", genEventCount, "sumw=", genEventSumw)

    if maxEntriesPerSample is not None:
        print(f"Scaling to {maxEntriesPerSample} entries")
        if nEntries>maxEntriesPerSample :
            genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
            nEntries=maxEntriesPerSample
        print("    scaled to:", nEntries, "sumw=", genEventSumw)

    return genEventSumw

def get_p4(input_file):
    '''
       Util function that returns collections of gen-level
       and reco-level leptons (returns their pt, eta, phi, mass)
       that create ZZ candidates.
    '''
    iEntry = 0

    f = input_file

    event = f.Events
    event.SetBranchStatus("*", 0)
    event.SetBranchStatus("run", 1)
    event.SetBranchStatus("luminosityBlock", 1)
    event.SetBranchStatus("*Muon*", 1)
    event.SetBranchStatus("*Electron*", 1)
    event.SetBranchStatus("*ZZCand*", 1)
    event.SetBranchStatus("bestCandIdx", 1)
    event.SetBranchStatus("FidDressedLeps*", 1)
    event.SetBranchStatus("HLT_passZZ4l", 1)
    nEntries = event.GetEntries()

    leps_pt = []
    leps_eta = []
    leps_phi = []
    leps_mass = []

    gen_leps_pt = []
    gen_leps_eta = []
    gen_leps_phi = []
    gen_leps_mass = []

    while iEntry<nEntries and event.GetEntry(iEntry):
        iEntry+=1
        bestCandIdx  = event.bestCandIdx
        if iEntry%10000 == 0: print(f".... Processing {iEntry}th event ....")
        if(bestCandIdx != -1):
            ZZs = Collection(event, 'ZZCand')
            theZZ = ZZs[bestCandIdx]
            
            l1, l2, l3, l4 = getLeptons(theZZ, event)
            gleps = Collection(event,'FidDressedLeps')
            
            leps_pt.append((l1.pt, l2.pt, l3.pt, l4.pt))
            leps_eta.append((l1.eta, l2.eta, l3.eta, l4.eta))
            leps_phi.append((l1.phi, l2.phi, l3.phi, l4.phi))
            leps_mass.append((l1.mass, l2.mass, l3.mass, l4.mass))

            g_leps_pt = []
            g_leps_eta = []
            g_leps_phi = []
            g_leps_mass = []
            
            for gl in gleps:
                g_leps_pt.append(gl.pt)
                g_leps_eta.append(gl.eta)
                g_leps_phi.append(gl.phi)
                g_leps_mass.append(gl.mass)

            gen_leps_pt.append(tuple(g_leps_pt))
            gen_leps_eta.append(tuple(g_leps_eta))
            gen_leps_phi.append(tuple(g_leps_phi))
            gen_leps_mass.append(tuple(g_leps_mass))

    reco_leps = [leps_pt, leps_eta, leps_phi, leps_mass]
    gen_leps  = [gen_leps_pt, gen_leps_eta, gen_leps_phi, gen_leps_mass]

    return reco_leps, gen_leps

def gen_reco_matching(reco_leps, gen_leps):
    '''
       Util function that performs reco-to-gen matching
       for leptons that form ZZ candidates.       
    '''
    leps_pt, leps_eta, leps_phi, leps_mass = reco_leps
    gen_leps_pt, gen_leps_eta, gen_leps_phi, gen_leps_mass = gen_leps
    
    lep_genindex = []
    for i in range(len(leps_pt)):
        _lep_genindex = []
        for lep in leps_pt[i]:
            _lep_genindex.append(-1)

        for idx_l, lep in enumerate(leps_pt[i]):
            minDr = 9999.0
            reco = ROOT.TLorentzVector()
            gen = ROOT.TLorentzVector()
            reco.SetPtEtaPhiM(leps_pt[i][idx_l], leps_eta[i][idx_l], leps_phi[i][idx_l], leps_mass[i][idx_l])
            for idx_g, gl in enumerate(gen_leps_pt[i]):
                gen.SetPtEtaPhiM(gen_leps_pt[i][idx_g], gen_leps_eta[i][idx_g], gen_leps_phi[i][idx_g], gen_leps_mass[i][idx_g])
                dr = reco.DeltaR(gen)
                if((dr < minDr) and (dr<0.5)):
                    _lep_genindex[idx_l] = idx_g
                    minDr = dr
                    
        lep_genindex.append(_lep_genindex)

    return ak.Array(lep_genindex)

def get_h_indices(reco_leps):
    '''
       Util functions that returns the indices of
       4 leptons creating the ZZ candidate.
    '''
    leps_pt, leps_eta, leps_phi, leps_mass = reco_leps

    lep_hindex = []
    for i in range(len(leps_pt)):
        _lep_hindex = []
        for lep in leps_pt[i]:
            _lep_hindex.append(-1)

        lead_Z1 = max(leps_pt[i][0],leps_pt[i][1])
        if (lead_Z1 == leps_pt[i][0]):
            _lep_hindex[0] = 0
            _lep_hindex[1] = 1
        elif(lead_Z1 == leps_pt[i][1]):
            _lep_hindex[0] = 1
            _lep_hindex[1] = 0

        lead_Z2 = max(leps_pt[i][2],leps_pt[i][3])
        if (lead_Z2 == leps_pt[i][2]):
            _lep_hindex[2] = 2
            _lep_hindex[3] = 3
        elif(lead_Z2 == leps_pt[i][3]):
            _lep_hindex[2] = 3
            _lep_hindex[3] = 2

        lep_hindex.append(_lep_hindex)
        
    return ak.Array(lep_hindex)

def add_branches(filename, lep_genindex, lep_hindex):
    '''
       Function that takes a CJLST-processed nanoAOD
       file and adds to the Events and AllEvents TTrees
       the branches for lep_genindex, lep_Hindex and Counter
       created by get_h_indices, gen_reco_matching, and get_genEventSumw.
    '''
    with uproot.open(f"{filename}") as f:
        zz4l_branches = f['Events'].arrays()
        all_branches = f['AllEvents'].arrays()

    d_vals = defaultdict(list)
    d_types = defaultdict(list)

    d_vals_all = defaultdict(list)
    d_types_all = defaultdict(list)

    for b in zz4l_branches.fields:
        d_types[b] = zz4l_branches[b].type
        d_vals[b] = zz4l_branches[b]

    for b in all_branches.fields:
        d_types_all[b] = all_branches[b].type
        d_vals_all[b] = all_branches[b]
        
    d_types['lep_genindex'] = ak.Array(lep_genindex).type
    d_vals['lep_genindex'] = ak.Array(lep_genindex)

    d_types['lep_Hindex'] = ak.Array(lep_hindex).type
    d_vals['lep_Hindex'] = ak.Array(lep_hindex)

    d_vals['Counter'] = ak.ones_like(d_vals['event'])*genEventSumw
    d_types['Counter'] = d_vals['Counter'].type

    d_vals_all['Counter'] = ak.ones_like(d_vals_all['event'])*genEventSumw
    d_types_all['Counter'] = d_vals_all['Counter'].type

    return d_vals, d_types, d_vals_all, d_types_all

if __name__=="main":
    filename = "ggH125_reducedTree_MC_2018.root"
    out_file = "ggH125_reducedTree_MC_2018_lepindex.root"

    f = ROOT.TFile.Open(filename)

    reco_leps, gen_leps = get_p4(f)
    genEventSumw = get_genEventSumw(f)
    lep_genindex = gen_reco_matching(reco_leps, gen_leps)
    lep_hindex = get_h_indices(reco_leps)

    d_vals, d_types, d_vals_all, d_types_all = add_branches(filename, lep_genindex, lep_hindex)

    with uproot.recreate(f"{out_file}") as fout:
        fout.mktree("Events", d_types)
        fout["Events"].extend(d_vals)

        fout.mktree("AllEvents", d_types_all)
        fout["AllEvents"].extend(d_vals_all)
