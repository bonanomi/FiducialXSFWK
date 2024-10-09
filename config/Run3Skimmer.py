## python3 NanoConverter.py
import ROOT
import os,sys
import optparse
from ZZAnalysis.NanoAnalysis.tools import get_genEventSumw

usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
parser = optparse.OptionParser(usage)
parser.add_option('',   '--mc', action='store_true', dest='MC', default=False, help='MC samples')
parser.add_option('',   '--input',dest='INPUT',type='string',default='', help='Name and path to the input file')
parser.add_option('',   '--output',dest='OUTPUT',type='string',default='', help='Name and path to the output file')
(opt, args) = parser.parse_args()

MC = opt.MC

inFileName = opt.INPUT
outFileName = opt.OUTPUT

df = ROOT.RDataFrame('Events', inFileName)
df = df.Range(0, 100)

print(f"df nEntries: {df.Count().GetValue()}")

if "H12" in inFileName:
    df_all = ROOT.RDataFrame('AllEvents', inFileName)
    print(f"df_all nEntries: {df_all.Count().GetValue()}")

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = 'RECREATE'

events = df.Take[df.GetColumnType('event')]('event')
event_values = df.AsNumpy(["event"])["event"]
event_list = list(set(event_values))
 
ROOT.gInterpreter.Declare(f"""
    #include <vector>
    #include <set>

    std::vector<unsigned long long> getEventVector() {{
        return std::vector<unsigned long long>({{ {', '.join(map(str, event_list))} }});
    }}

    auto eventVector = getEventVector();
    auto eventSet = std::set<unsigned long long>(eventVector.begin(), eventVector.end());
    auto isNotInEvents = [](unsigned long long event) {{
        return eventSet.find(event) == eventSet.end();
    }};
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<float> concatenate(ROOT::RVec<float> &A, ROOT::RVec<float> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<short> concatenate(ROOT::RVec<int> &A, ROOT::RVec<int> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<bool> concatenate(ROOT::RVec<bool> &A, ROOT::RVec<bool> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> concatenate(ROOT::RVec<unsigned char> &A, ROOT::RVec<unsigned char> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> addDummyBranch(ROOT::RVec<float> &A){
    int sizeA = A.size();
    ROOT::RVec<unsigned char> B(sizeA);

    for (int i = 0; i < sizeA; ++i) {
        B[i] = 0;
    }

    return B;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<int> getHindex(std::vector<float> LepPt) {
      std::vector<int> _lep_Hindex;
      for(unsigned int i=0;i<LepPt.size();i++){
        _lep_Hindex.push_back(-1);
      }
      float lead_Z1 = max(LepPt.at(0),LepPt.at(1));
      if(lead_Z1 == LepPt.at(0)){
        _lep_Hindex[0] = 0;
        _lep_Hindex[1] = 1;
      }
      else if(lead_Z1 == LepPt.at(1)){
        _lep_Hindex[0] = 1;
        _lep_Hindex[1] = 0;
      }
      float lead_Z2 = max(LepPt.at(2),LepPt.at(3));
      if(lead_Z2 == LepPt.at(2)){
        _lep_Hindex[2] = 2;
        _lep_Hindex[3] = 3;
      }
      else if(lead_Z2 == LepPt.at(3)){
        _lep_Hindex[2] = 3;
        _lep_Hindex[3] = 2;
      }

      return _lep_Hindex;
}

""")

ROOT.gInterpreter.Declare("""
std::vector<int> getGENHindex(int FidZZ_Z1l1Idx, int FidZZ_Z1l2Idx, int FidZZ_Z2l1Idx, int FidZZ_Z2l2Idx) {
    std::vector<int> GENlep_Hindex;
    GENlep_Hindex.push_back(FidZZ_Z1l1Idx);
    GENlep_Hindex.push_back(FidZZ_Z1l2Idx);
    GENlep_Hindex.push_back(FidZZ_Z2l1Idx);
    GENlep_Hindex.push_back(FidZZ_Z2l2Idx);
    return GENlep_Hindex;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<int> getLepGENindex(std::vector<float> LepPt, std::vector<float> LepEta, std::vector<float> LepPhi, ROOT::VecOps::RVec<float> GENlep_id, std::vector<short> LepLepId, ROOT::VecOps::RVec<float> GENlep_pt, ROOT::VecOps::RVec<float> GENlep_eta, ROOT::VecOps::RVec<float> GENlep_phi){
    std::vector<int> lep_genindex;
    for(unsigned int i = 0; i < LepPt.size(); i++){
        lep_genindex.push_back(-1);
    }
    for(unsigned int i = 0; i < LepPt.size(); i++){
        double minDr = 9999.0;
        TLorentzVector reco, gen;
        reco.SetPtEtaPhiM(LepPt.at(i), LepEta.at(i), LepPhi.at(i), 0.0001); 
        for(unsigned int j = 0; j < GENlep_id.size(); j++){
            if (GENlep_id.at(j) != LepLepId.at(i)) continue;
            gen.SetPtEtaPhiM(GENlep_pt.at(j), GENlep_eta.at(j), GENlep_phi.at(j), 0.0001);
            double deta = reco.Eta() - gen.Eta();
            double dphi = abs(reco.Phi() - gen.Phi());
            if (dphi > 3.1415){
                dphi = dphi - 2*3.1415;
            }
            double thisDr = sqrt(deta*deta + dphi*dphi);
            if (thisDr<minDr && thisDr<0.5) {
                lep_genindex[i] = j;
                minDr = thisDr;
            }
        }
    }
    return lep_genindex;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<float> getGENlep_vector(ROOT::VecOps::RVec<float> GENlep_input){
    std::vector<float> GENlep_vector;
    for(unsigned int i = 0; i < GENlep_input.size(); i++){
        GENlep_vector.push_back(GENlep_input.at(i));
    }
    return GENlep_vector;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<int> getGENlep_int(ROOT::VecOps::RVec<int> GENlep_input){
    std::vector<int> GENlep_vector;
    for(unsigned int i = 0; i < GENlep_input.size(); i++){
        GENlep_vector.push_back(GENlep_input.at(i));
    }
    return GENlep_vector;
}
""")

## Variables to store in the output root file
vars = {'RunNumber',
        'EventNumber',
        'LumiNumber',
        'ZZMass',
        'ZZPt',
        'ZZy',
        # 'CRflag',
        'Z1Flav',
        'Z2Flav',
        'Z1Mass',
        'Z2Mass',
        'dataMCWeight',
        # 'PFMET',
        'LepPt',
        'LepEta',
        'LepPhi',
        'Lepdxy',
        'Lepdz',
        'LepLepId',
        'LepSIP',
        'LepCombRelIsoPF',
        # 'LepMissingHit',
        }
if MC:
    vars.add('overallEventWeight')
    # vars.add('xsec')
    # vars.add('L1prefiringWeight')
    # vars.add('LHEPdfWeight')
    # vars.add('LHEScaleWeight')
    vars.add('PUWeight')
    # vars.add('LHEWeight_originalXWGTUP')
    vars.add('lep_Hindex')
    if "H12" in inFileName:
        vars.add('GENlep_pt')
        vars.add('GENlep_eta')
        vars.add('GENlep_phi')
        vars.add('GENlep_mass')
        vars.add('GENlep_id')
        vars.add('GENlep_MomId')
        vars.add('GENlep_MomMomId')
        vars.add('GENlep_RelIso')
        vars.add('GENmass4l')
        vars.add('GENpT4l')
        vars.add('GENeta4l')
        vars.add('GENphi4l')
        vars.add('GENrapidity4l')
        vars.add('GENZ_DaughtersId')
        vars.add('GENZ_MomId')
        vars.add('GENlep_Hindex')
        # vars.add('GENlep_index')
        vars.add('lep_genindex')
        vars.add('passedFiducial')
    vars.add('passedFullSelection')
    vars.add('genHEPMCweight')
    if 'ggH' in inFileName:
        vars.add('ggH_NNLOPS_weight')
    if 'ZZTo' in inFileName:
        # vars.add('KFactor_EW_qqZZ_Weight')
        vars.add('KFactor_QCD_qqZZ_M_Weight')
    if 'ggTo' in inFileName:
        vars.add('KFactor_QCD_ggZZ_Nominal_Weight')

## SR
df_SR = ( df.Filter('bestCandIdx>=0').Define("ZZMass", "ZZCand_mass[bestCandIdx]") ## Dummy
                                     .Define("ZZPt", "ZZCand_pt[bestCandIdx]")
                                     .Define("ZZy", "abs(ZZCand_rapidity[bestCandIdx])")
                                     # .Define("CRflag", "0") ## Dummy
                                     .Define("Z1Flav", "ZZCand_Z1flav[bestCandIdx]")
                                     .Define("Z2Flav", "ZZCand_Z2flav[bestCandIdx]") ## Dummy
                                     .Define("Z1Mass", "ZZCand_Z1mass[bestCandIdx]")
                                     .Define("Z2Mass", "ZZCand_Z2mass[bestCandIdx]") ## Dummy
                                     .Define("dataMCWeight", "ZZCand_dataMCWeight[bestCandIdx]")
                                     .Define('RunNumber', "run")
                                     .Define('EventNumber', "event")
                                     .Define('LumiNumber', "luminosityBlock")
                                     .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                     .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                     .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                     .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                     .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                     .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                     .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                     .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                     ## Need to add the LepMissingHit branch for SS FR method
                                     ## First create a dummy branch for muons filled with zeroes
                                     # .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                     # .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                     ## Variable miniAOD-style
                                     .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_pt[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_pt[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_pt[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepPt")
                                     .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_eta[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_eta[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_eta[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepEta")
                                     .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_phi[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_phi[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_phi[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepPhi")
                                     .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z2l2Idx[bestCandIdx]]}; return Lepdxy")
                                     .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_dz[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_dz[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_dz[ZZCand_Z2l2Idx[bestCandIdx]]}; return Lepdz")
                                     .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_id[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_id[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_id[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepLepId")
                                     .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_sip[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_sip[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_sip[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepSIP")
                                     .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_iso[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_iso[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_iso[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepCombRelIsoPF")
                                     # .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepMissingHit")
                                     # .Define('PFMET', "MET_pt")
                                     .Define('lep_Hindex', "getHindex(LepPt)")
                                     .Define('passedFullSelection', "1")
                                     .Define('genHEPMCweight', "Generator_weight")
                                     .Define('PUWeight', "puWeight")
                                     )

if 'ggH' in inFileName:
    df_SR = df_SR.Define('ggH_NNLOPS_weight', "ggH_NNLOPS_Weight")
            
if "H12" in inFileName:
    df_SR = (df_SR.Define('GENlep_pt', "getGENlep_vector(FidDressedLeps_pt)")
                  .Define('GENlep_eta', "getGENlep_vector(FidDressedLeps_eta)")
                  .Define('GENlep_phi', "getGENlep_vector(FidDressedLeps_phi)")
                  .Define('GENlep_mass', "getGENlep_vector(FidDressedLeps_mass)")
                  .Define('GENlep_id', "getGENlep_vector(FidDressedLeps_id)")
                  .Define('GENlep_MomId', "getGENlep_vector(FidDressedLeps_momid)")
                  .Define('GENlep_MomMomId', "getGENlep_vector(FidDressedLeps_mommomid)")
                  .Define('GENlep_RelIso', "getGENlep_vector(FidDressedLeps_RelIso)")
                  .Define('GENmass4l', "FidZZ_mass")
                  .Define('GENpT4l', "FidZZ_pt")
                  .Define('GENeta4l', "FidZZ_eta")
                  .Define('GENphi4l', "FidZZ_phi")
                  .Define('GENrapidity4l', "FidZZ_rapidity")
                  .Define('GENZ_DaughtersId', "getGENlep_int(FidZ_DauPdgId)")
                  .Define('GENZ_MomId', "getGENlep_int(FidZ_MomPdgId)")
                  .Define('GENlep_Hindex', "getGENHindex(FidZZ_Z1l1Idx, FidZZ_Z1l2Idx, FidZZ_Z2l1Idx, FidZZ_Z2l2Idx)")
                  .Define('lep_genindex', "getLepGENindex(LepPt, LepEta, LepPhi, GENlep_id, LepLepId, GENlep_pt, GENlep_eta, GENlep_phi)")
            )

    df_filter = df_all.Filter("isNotInEvents(event)")
    print(f"df_filter nEntries: {df_filter.Count().GetValue()}")
    vars_fail = {'GENlep_pt',
    'GENlep_eta',
    'GENlep_phi',
    'GENlep_mass',
    'GENlep_id',
    'GENlep_MomId',
    'GENlep_MomMomId',
    'GENlep_RelIso',
    'GENmass4l',
    'GENpT4l',
    'GENeta4l',
    'GENphi4l',
    'GENrapidity4l',
    'GENZ_DaughtersId',
    'GENZ_MomId',
    'GENlep_Hindex',
    'RunNumber',
    'EventNumber',
    'LumiNumber',
    'passedFiducial',
    'passedFullSelection',
    # 'LHEPdfWeight',
    # 'LHEScaleWeight',
    'genHEPMCweight',
    'PUWeight'
    }
    df_fail = (df_filter.Define('GENlep_pt', "getGENlep_vector(FidDressedLeps_pt)")
                        .Define('GENlep_eta', "getGENlep_vector(FidDressedLeps_eta)")
                        .Define('GENlep_phi', "getGENlep_vector(FidDressedLeps_phi)")
                        .Define('GENlep_mass', "getGENlep_vector(FidDressedLeps_mass)")
                        .Define('GENlep_id', "getGENlep_vector(FidDressedLeps_id)")
                        .Define('GENlep_MomId', "getGENlep_vector(FidDressedLeps_momid)")
                        .Define('GENlep_MomMomId', "getGENlep_vector(FidDressedLeps_mommomid)")
                        .Define('GENlep_RelIso', "getGENlep_vector(FidDressedLeps_RelIso)")
                        .Define('GENmass4l', "FidZZ_mass")
                        .Define('GENpT4l', "FidZZ_pt")
                        .Define('GENeta4l', "FidZZ_eta")
                        .Define('GENphi4l', "FidZZ_phi")
                        .Define('GENrapidity4l', "FidZZ_rapidity")
                        .Define('GENZ_DaughtersId', "getGENlep_int(FidZ_DauPdgId)")
                        .Define('GENZ_MomId', "getGENlep_int(FidZ_MomPdgId)")
                        .Define('GENlep_Hindex', "getGENHindex(FidZZ_Z1l1Idx, FidZZ_Z1l2Idx, FidZZ_Z2l1Idx, FidZZ_Z2l2Idx)")
                        .Define('RunNumber', "run")
                        .Define('EventNumber', "event")
                        .Define('LumiNumber', "luminosityBlock")
                        .Define('passedFullSelection', "0")
                        .Define('genHEPMCweight', "Generator_weight")  #LHEWeight_originalXWGTUP")
                        .Define('PUWeight', "puWeight")
                )
    if 'ggH' in inFileName:
        vars_fail.add('ggH_NNLOPS_weight')
        df_fail = df_fail.Define('ggH_NNLOPS_weight', "ggH_NNLOPS_Weight")
    
opts.fMode = 'UPDATE'
print(f"df nEntries: {df_SR.Count().GetValue()}")

df_SR.Snapshot('ZZTree/candTree', outFileName, vars, opts)
if "H12" in inFileName:
    df_fail.Snapshot('ZZTree/candTree_failed', outFileName, vars_fail, opts)

## Add counter only with the 40th entry
counters = ROOT.TH1F("Counters", "Counters", 50, 0, 100)
if opt.MC:
    root = ROOT.TFile.Open(inFileName)
    genEventSumw = get_genEventSumw(root)
    counters.SetBinContent(40, genEventSumw)
    root.Close()

root_file = ROOT.TFile(outFileName, "UPDATE")
counters.Write()
root_file.Close()
