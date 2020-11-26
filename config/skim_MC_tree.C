// root -l skim_MC_tree.C (prod_mod) (year)

#include<iostream>
#include<fstream>
#include<cstring>
#include<math.h>
#include<tuple>
#include<cmath>
#include<TMath.h>
#include<TCanvas.h>
#include<TGraph.h>
#include<TTreeReader.h>
#include<TFile.h>
#include<TTree.h>
#include<TObject.h>
#include<TH2.h>
#include<TF1.h>
#include<TApplication.h>
#include<TGraphErrors.h>
#include<TLorentzVector.h>
#include<TMultiGraph.h>
#include<TLegend.h>
#include<TLatex.h>
#include<TStyle.h>
#include<random>
#include<algorithm>
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

float mass_lep(int flavour){
  if((abs(flavour)) == 11) return 0.0005109989461;
  else if ((abs(flavour)) == 13) return 0.1056583745;
  else if ((abs(flavour)) == 15) return 1.77686;
  // else if ((abs(flavour)) == 0) return 0;
  return 0;
}


pair<vector<TLorentzVector>,vector<Short_t>> sort(vector<TLorentzVector> lep,vector<Short_t>id){
  vector<TLorentzVector> sortedLep;
  vector<Short_t> sortedId;
  for (int k=0;k<4;k++){
    int max = -1;
    float maxPt = 0;
    for (int i=0;i<lep.size();i++){
      if (lep[i].Pt() > maxPt) {
        max = i;
        maxPt = lep[i].Pt();
      }
    }
    sortedLep.push_back(lep[max]);
    sortedId.push_back(id[max]);
    lep.erase(lep.begin()+max);
    id.erase(id.begin()+max);
  }
  return make_pair(sortedLep,sortedId);
}
pair<vector<TLorentzVector>,vector<Short_t>> sort_2(vector<TLorentzVector> lep,vector<Short_t>id){
  vector<TLorentzVector> sortedLep;
  vector<Short_t> sortedId;
  for (int k=0;k<2;k++){
    int max = -1;
    float maxPt = 0;
    for (int i=0;i<lep.size();i++){
      if (lep[i].Pt() > maxPt) {
        max = i;
        maxPt = lep[i].Pt();
      }
    }
    sortedLep.push_back(lep[max]);
    sortedId.push_back(id[max]);
    lep.erase(lep.begin()+max);
  }
  return make_pair(sortedLep,sortedId);
}



int Fiducial(float Z1Flav,float Z2Flav,float Z1Mass,float Z2Mass,float lep1Iso,float lep2Iso,float lep3Iso,float lep4Iso,
              vector<TLorentzVector> lepSorted,vector<Short_t> lepIdSorted,bool iso){
  //Two same-flavour OS lepton pairs
  //If they are greater than zero there two SS leptons
  //If they are equal to zero there are no leptons
  //If they are different from 121 or 169 there are taus
  if(Z1Flav!=-121){
    if(Z1Flav!=-169){
      return false;
    }
  }
  if(Z2Flav!=-121){
    if(Z2Flav!=-169){
      return false;
    }
  }

  //Pseudorapidity + minimal common pt
  for(int i=0;i<lepSorted.size();i++){
    if(abs(lepIdSorted[i]) == 11){ //Electrons
      if((abs(lepSorted[i].Eta())>2.5) || (lepSorted[i].Pt()<7)) return false;
    }
    if(abs(lepIdSorted[i]) == 13){ //Muons
      if((abs(lepSorted[i].Eta())>2.4) || (lepSorted[i].Pt()<5)) return false;
    }
  }

  //Leading lepton
  if(lepSorted[0].Pt()<20) return false;

  //Next-to-leading lepton
  if(lepSorted[1].Pt()<10) return false;

  //Z's invariant mass
  if((Z1Mass<40) || (Z1Mass>120)) return false;
  if((Z2Mass<12) || (Z2Mass>120)) return false;

  //Invariant mass of the selected four leptons
  float m4l = (lepSorted[0]+lepSorted[1]+lepSorted[2]+lepSorted[3]).M();
  if((m4l<105) || (m4l>140)) return false;

  //Distance between the selected four leptons
  if(lepSorted[0].DeltaR(lepSorted[1])<0.02) return false;
  if(lepSorted[0].DeltaR(lepSorted[2])<0.02) return false;
  if(lepSorted[0].DeltaR(lepSorted[3])<0.02) return false;
  if(lepSorted[1].DeltaR(lepSorted[2])<0.02) return false;
  if(lepSorted[1].DeltaR(lepSorted[3])<0.02) return false;
  if(lepSorted[2].DeltaR(lepSorted[3])<0.02) return false;

  //Invariant mass of any opposire sign lepton pair
  if(lepIdSorted[0]!=lepIdSorted[1]){
    float invMass = (lepSorted[0] + lepSorted[1]).M();
    if(invMass<4) return false;
  }
  if(lepIdSorted[0]!=lepIdSorted[2]){
    float invMass = (lepSorted[0] + lepSorted[2]).M();
    if(invMass<4) return false;
  }
  if(lepIdSorted[0]!=lepIdSorted[3]){
    float invMass = (lepSorted[0] + lepSorted[3]).M();
    if(invMass<4) return false;
  }
  if(lepIdSorted[1]!=lepIdSorted[2]){
    float invMass = (lepSorted[1] + lepSorted[2]).M();
    if(invMass<4) return false;
  }
  if(lepIdSorted[1]!=lepIdSorted[3]){
    float invMass = (lepSorted[1] + lepSorted[3]).M();
    if(invMass<4) return false;
  }
  if(lepIdSorted[2]!=lepIdSorted[3]){
    float invMass = (lepSorted[2] + lepSorted[3]).M();
    if(invMass<4) return false;
  }

  //Isolation
  if(iso==true){
    if(lep1Iso>0.35) return false;
    if(lep2Iso>0.35) return false;
    if(lep3Iso>0.35) return false;
    if(lep4Iso>0.35) return false;
  }

  //If all the conditions are satisfied -> events belonging to the fiducial phase space
  return true;
}


  //------------------------------------------------------------------
void add(TString input_dir, TString year, TString prod_mode, bool t_failed=true){
  // Add additional branches
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path = Form("%s/%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  TFile *f = new TFile(new_full_path.Data(),"UPDATE");
  TTree *T;
  if (t_failed) {
    std::cout << "Filling TTree _failed" << std::endl;
    T = (TTree*)f->Get("candTree_failed");
  } else {
    std::cout << "Filling TTree ZZ cands" << std::endl;
    T = (TTree*)f->Get("candTree");
  }
  std::cout << T->GetName() << std::endl;
  float GenLep1Pt,GenLep2Pt,GenLep3Pt,GenLep4Pt,GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta,GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi,GenZ1Flav,GenZ2Flav,
        GenZ1Mass,GenZ2Mass,GenLep1Iso,GenLep2Iso,GenLep3Iso,GenLep4Iso;
  Short_t GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  Short_t _GENnjets_pt30_eta2p5;
  bool _passedFiducialSelection,_passedFiducialSelection_NOISO,_passedFullSelection;
  vector<float> _GenLepPtSorted,_GenLepEtaSorted,_GenLepPhiSorted;
  vector<Short_t> _GenLepIdSorted;
  vector<TLorentzVector> GenLepSorted;
  vector<float> *GENlep_pt = 0;
  vector<float> *GENlep_eta = 0;
  vector<float> *GENlep_phi = 0;
  vector<float> *GENlep_mass = 0;
  vector<float> *GENlep_id = 0;
  vector<float> *GenCleanedJetPt = 0;
  vector<float> *GenCleanedJetEta = 0;
  TBranch *GenLepPtSorted = T->Branch("GenLepPtSorted",&_GenLepPtSorted);
  TBranch *GenLepEtaSorted = T->Branch("GenLepEtaSorted",&_GenLepEtaSorted);
  TBranch *GenLepPhiSorted = T->Branch("GenLepPhiSorted",&_GenLepPhiSorted);
  TBranch *GenLepIdSorted = T->Branch("GenLepIdSorted",&_GenLepIdSorted);
  TBranch *passedFiducialSelection = T->Branch("passedFiducialSelection",&_passedFiducialSelection,"passedFiducialSelection/B");
  TBranch *passedFiducialSelection_NOISO = T->Branch("passedFiducialSelection_NOISO",&_passedFiducialSelection_NOISO,"passedFiducialSelection_NOISO/B");
  TBranch *passedFullSelection = T->Branch("passedFullSelection",&_passedFullSelection,"passedFullSelection/B");
  TBranch *GENnjets_pt30_eta2p5 = T->Branch("GENnjets_pt30_eta2p5",&_GENnjets_pt30_eta2p5,"GENnjets_pt30_eta2p5/S");

  T->SetBranchAddress("GenLep1Pt",&GenLep1Pt);
  T->SetBranchAddress("GenLep2Pt",&GenLep2Pt);
  T->SetBranchAddress("GenLep3Pt",&GenLep3Pt);
  T->SetBranchAddress("GenLep4Pt",&GenLep4Pt);
  T->SetBranchAddress("GenLep1Eta",&GenLep1Eta);
  T->SetBranchAddress("GenLep2Eta",&GenLep2Eta);
  T->SetBranchAddress("GenLep3Eta",&GenLep3Eta);
  T->SetBranchAddress("GenLep4Eta",&GenLep4Eta);
  T->SetBranchAddress("GenLep1Phi",&GenLep1Phi);
  T->SetBranchAddress("GenLep2Phi",&GenLep2Phi);
  T->SetBranchAddress("GenLep3Phi",&GenLep3Phi);
  T->SetBranchAddress("GenLep4Phi",&GenLep4Phi);
  T->SetBranchAddress("GenLep1Id",&GenLep1Id);
  T->SetBranchAddress("GenLep2Id",&GenLep2Id);
  T->SetBranchAddress("GenLep3Id",&GenLep3Id);
  T->SetBranchAddress("GenLep4Id",&GenLep4Id);
  T->SetBranchAddress("GenLep1Iso",&GenLep1Iso);
  T->SetBranchAddress("GenLep2Iso",&GenLep2Iso);
  T->SetBranchAddress("GenLep3Iso",&GenLep3Iso);
  T->SetBranchAddress("GenLep4Iso",&GenLep4Iso);
  T->SetBranchAddress("GenZ1Flav",&GenZ1Flav);
  T->SetBranchAddress("GenZ2Flav",&GenZ2Flav);
  T->SetBranchAddress("GenZ1Mass",&GenZ1Mass);
  T->SetBranchAddress("GenZ2Mass",&GenZ2Mass);
  T->SetBranchAddress("GENlep_pt",&GENlep_pt);
  T->SetBranchAddress("GENlep_eta",&GENlep_eta);
  T->SetBranchAddress("GENlep_phi",&GENlep_phi);
  T->SetBranchAddress("GENlep_mass",&GENlep_mass);
  T->SetBranchAddress("GENlep_id",&GENlep_id);
  T->SetBranchAddress("GenCleanedJetPt",&GenCleanedJetPt);
  T->SetBranchAddress("GenCleanedJetEta",&GenCleanedJetEta);

  float _ZZy,ZZPt,ZZEta;
  Short_t nCleanedJetsPt30,nCleanedJetsPt30_jesUp,nCleanedJetsPt30_jesDn;
  Short_t _njets_pt30_eta2p5,_njets_pt30_eta2p5_jesup,_njets_pt30_eta2p5_jesdn;
  Float_t _pTj1, _GENpTj1;
  vector<float> *LepPt = 0;
  vector<float> *LepPhi = 0;
  vector<float> *LepEta = 0;
  vector<int> *LepLepId = 0;
  vector<float> *ExtraLepPt = 0;
  vector<float> *ExtraLepEta = 0;
  vector<float> *ExtraLepPhi = 0;
  vector<int> *ExtraLepLepId = 0;
  vector<int> _lep_genindex, _lep_Hindex;
  vector<float> *JetPt = 0;
  vector<float> *JetEta = 0;
  vector<float> *JetPt_JESUp = 0;
  vector<float> *JetPt_JESDown = 0;


  TBranch *ZZy = T->Branch("ZZy",&_ZZy,"ZZy/F");
  TBranch *lep_genindex = T->Branch("lep_genindex",&_lep_genindex);
  TBranch *lep_Hindex = T->Branch("lep_Hindex",&_lep_Hindex);
  TBranch *njets_pt30_eta2p5 = T->Branch("njets_pt30_eta2p5",&_njets_pt30_eta2p5,"njets_pt30_eta2p5/S");
  TBranch *njets_pt30_eta2p5_jesup = T->Branch("njets_pt30_eta2p5_jesup",&_njets_pt30_eta2p5_jesup,"_njets_pt30_eta2p5_jesup/S");
  TBranch *njets_pt30_eta2p5_jesdn = T->Branch("njets_pt30_eta2p5_jesdn",&_njets_pt30_eta2p5_jesdn,"_njets_pt30_eta2p5_jesdn/S");

  TBranch *pTj1 = T->Branch("pTj1",&_pTj1,"pTj1/F");
  TBranch *GENpTj1 = T->Branch("GENpTj1",&_GENpTj1,"GENpTj1/F");

  if (!t_failed) {
    T->SetBranchAddress("ZZPt",&ZZPt);
    T->SetBranchAddress("ZZEta",&ZZEta);
    T->SetBranchAddress("LepPt",&LepPt);
    T->SetBranchAddress("LepPhi",&LepPhi);
    T->SetBranchAddress("LepEta",&LepEta);
    T->SetBranchAddress("LepLepId",&LepLepId);
    T->SetBranchAddress("ExtraLepPt",&ExtraLepPt);
    T->SetBranchAddress("ExtraLepEta",&ExtraLepEta);
    T->SetBranchAddress("ExtraLepPhi",&ExtraLepPhi);
    T->SetBranchAddress("ExtraLepLepId",&ExtraLepLepId);
    T->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
    T->SetBranchAddress("nCleanedJetsPt30_jesUp",&nCleanedJetsPt30_jesUp);
    T->SetBranchAddress("nCleanedJetsPt30_jesDn",&nCleanedJetsPt30_jesDn);
    T->SetBranchAddress("JetPt",&JetPt);
    T->SetBranchAddress("JetEta",&JetEta);
    T->SetBranchAddress("JetPt_JESUp",&JetPt_JESUp);
    T->SetBranchAddress("JetPt_JESDown",&JetPt_JESDown);
  }

  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);

    // Sort GenLeptons
    float GenLep1Mass = mass_lep(GenLep1Id);
    float GenLep2Mass = mass_lep(GenLep2Id);
    float GenLep3Mass = mass_lep(GenLep3Id);
    float GenLep4Mass = mass_lep(GenLep4Id);
    vector<Short_t> GenLepId {GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id};
    TLorentzVector t1,t2,t3,t4; // Lorentz vector of the four genleptons
    t1.SetPtEtaPhiM(GenLep1Pt,GenLep1Eta,GenLep1Phi,GenLep1Mass);
    t2.SetPtEtaPhiM(GenLep2Pt,GenLep2Eta,GenLep2Phi,GenLep2Mass);
    t3.SetPtEtaPhiM(GenLep3Pt,GenLep3Eta,GenLep3Phi,GenLep3Mass);
    t4.SetPtEtaPhiM(GenLep4Pt,GenLep4Eta,GenLep4Phi,GenLep4Mass);
    vector<TLorentzVector> GenLep {t1,t2,t3,t4};
    pair<vector<TLorentzVector>,vector<Short_t>> sortedGenLeptons;
    if(t3.Pt()!=0) sortedGenLeptons = sort(GenLep,GenLepId);
    else sortedGenLeptons = sort_2(GenLep,GenLepId); /* Different function for cases in which there are just two GenLeptons, otherwise strange errors with
                                                        the previous function */
    GenLepSorted = sortedGenLeptons.first;
    _GenLepIdSorted = sortedGenLeptons.second;
    for(int i=0;i<GenLepSorted.size();i++){
      _GenLepPtSorted.push_back(GenLepSorted.at(i).Pt());
      _GenLepEtaSorted.push_back(GenLepSorted.at(i).Eta());
      _GenLepPhiSorted.push_back(GenLepSorted.at(i).Phi());
    }

    // Fiducial selections
    _passedFiducialSelection = Fiducial(GenZ1Flav,GenZ2Flav,GenZ1Mass,GenZ2Mass,GenLep1Iso,GenLep2Iso,GenLep3Iso,GenLep4Iso,GenLepSorted,_GenLepIdSorted,true);
    _passedFiducialSelection_NOISO = Fiducial(GenZ1Flav,GenZ2Flav,GenZ1Mass,GenZ2Mass,GenLep1Iso,GenLep2Iso,GenLep3Iso,GenLep4Iso,GenLepSorted,_GenLepIdSorted,false);

    _passedFullSelection = true;
    if (t_failed) {
	     _passedFullSelection = false;
    }

    // GenJet variables (kinematical cuts pT>30 && abs(eta)<2.5)
    _GENnjets_pt30_eta2p5 = 0;
    for(unsigned int i =0;i<GenCleanedJetPt->size();i++){
      if(GenCleanedJetPt->at(i)>30 && abs(GenCleanedJetEta->at(i))<2.5){
        _GENnjets_pt30_eta2p5++;
      }
    }

    // leading GENjet pT
    _GENpTj1 = 0;
    for (unsigned int i = 0; i < GenCleanedJetPt->size(); ++i)
    {
      if(GenCleanedJetPt->at(i) > 30. && abs(GenCleanedJetEta->at(i)) < 2.5 && GenCleanedJetPt->at(i) > _GENpTj1) {
        _GENpTj1 = GenCleanedJetPt->at(i);
      }
    }
    GENpTj1->Fill();
    GenLepPtSorted->Fill();
    GenLepEtaSorted->Fill();
    GenLepPhiSorted->Fill();
    GenLepIdSorted->Fill();
    passedFiducialSelection->Fill();
    passedFiducialSelection_NOISO->Fill();
    passedFullSelection->Fill();
    GENnjets_pt30_eta2p5->Fill();

    GenLepSorted.clear();
    _GenLepPtSorted.clear();
    _GenLepEtaSorted.clear();
    _GenLepPhiSorted.clear();
    _GenLepIdSorted.clear();

    if (t_failed) continue; // From now on reco-only variables

    // leading jet pT
    _pTj1 = 0;
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5 && JetPt->at(i) > _pTj1) {
        _pTj1 = JetPt->at(i);
      }
    }

    // njets
    _njets_pt30_eta2p5 = 0;
    for(unsigned int i=0;i<JetPt->size();i++){
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5){
        _njets_pt30_eta2p5++;
      }
    }
    // JES down
    _njets_pt30_eta2p5_jesdn = 0;
    for(unsigned int i=0;i<JetPt_JESDown->size();i++){
      if(JetPt_JESDown->at(i)>30 && abs(JetEta->at(i))<2.5){ // eta is not affected by jes
        _njets_pt30_eta2p5_jesdn++;
      }
    }
    // JES up
    _njets_pt30_eta2p5_jesup = 0;
    for(unsigned int i=0;i<JetPt_JESUp->size();i++){
      if(JetPt_JESUp->at(i)>30 && abs(JetEta->at(i))<2.5){ // eta is not affected by jes
        _njets_pt30_eta2p5_jesup++;
      }
    }
    pTj1->Fill();
    njets_pt30_eta2p5->Fill();
    njets_pt30_eta2p5_jesdn->Fill();
    njets_pt30_eta2p5_jesup->Fill();

    // Reco-rapidity
    _ZZy = abs(log((sqrt(125*125 + ZZPt*ZZPt*cosh(ZZEta)*cosh(ZZEta))+ZZPt*sinh(ZZEta))/sqrt(125*125+ZZPt*ZZPt)));
    ZZy->Fill();

    // GEN matching
    for(unsigned int i=0;i<LepPt->size();i++){
      _lep_genindex.push_back(-1);
    }
    for(unsigned int i = 0; i < LepPt->size(); i++) {
        double minDr=9999.0;
        TLorentzVector reco, gen;
        reco.SetPtEtaPhiM(LepPt->at(i),LepEta->at(i),LepPhi->at(i),mass_lep(LepLepId->at(i)));
        for (unsigned int j = 0; j < GENlep_id->size(); j++) {
            if (GENlep_id->at(j)!=LepLepId->at(i)) continue;
            gen.SetPtEtaPhiM(GENlep_pt->at(j),GENlep_eta->at(j),GENlep_phi->at(j),GENlep_mass->at(j));
            double thisDr = deltaR(reco.Eta(),reco.Phi(),gen.Eta(),gen.Phi());
            // double thisDr = reco.DeltaR(gen);
            if (thisDr<minDr && thisDr<0.5) {
                _lep_genindex[i]=j;
                minDr=thisDr;
            }
        } // all gen leptons
    } // all reco leptons
    lep_genindex->Fill();
    _lep_genindex.clear();

    for(unsigned int i=0;i<LepPt->size();i++){
      _lep_Hindex.push_back(-1);
    }
    float lead_Z1 = max(LepPt->at(0),LepPt->at(1));
    if(lead_Z1 == LepPt->at(0)){
      _lep_Hindex[0] = 0;
      _lep_Hindex[1] = 1;
    }
    else if(lead_Z1 == LepPt->at(1)){
      _lep_Hindex[0] = 1;
      _lep_Hindex[1] = 0;
    }
    float lead_Z2 = max(LepPt->at(2),LepPt->at(3));
    if(lead_Z2 == LepPt->at(2)){
      _lep_Hindex[2] = 2;
      _lep_Hindex[3] = 3;
    }
    else if(lead_Z2 == LepPt->at(3)){
      _lep_Hindex[2] = 3;
      _lep_Hindex[3] = 2;
    }
    lep_Hindex->Fill();
    _lep_Hindex.clear();

  }
  T->Write("", TObject::kOverwrite);
  delete f;
  return;
}

//---------------------------------------------------------- MAIN ----------------------------------------------------------
void skim_MC_tree (TString prod_mode = "ZH125", TString year = "2017"){

  TString input_dir = "/eos/user/a/atarabin/MC_samples";
  TString full_path = Form("%s/%s/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());

  std::cout << input_dir << " " << year << " " << prod_mode << std::endl;

  auto oldFile = TFile::Open(full_path.Data());
  TTree *oldtree = (TTree*) oldFile->Get("ZZTree/candTree");
  TTree *oldtree_failed = (TTree*) oldFile->Get("ZZTree/candTree_failed");
  TH1F *hCounters = (TH1F*) oldFile->Get("ZZTree/Counters");

  //// candTree
  // Deactivate all branches
  oldtree->SetBranchStatus("*",0);
  // Activate some branches only: our skim
  oldtree->SetBranchStatus("EventNumber",1);
  oldtree->SetBranchStatus("xsec",1);
  oldtree->SetBranchStatus("ZZMass",1);
  oldtree->SetBranchStatus("ZZPt",1);
  oldtree->SetBranchStatus("Z1Mass",1);
  oldtree->SetBranchStatus("Z2Mass",1);
  oldtree->SetBranchStatus("Z1Flav",1);
  oldtree->SetBranchStatus("Z2Flav",1);
  oldtree->SetBranchStatus("ZZPt",1);
  oldtree->SetBranchStatus("ZZEta",1);
  oldtree->SetBranchStatus("nCleanedJetsPt30",1);
  oldtree->SetBranchStatus("nCleanedJetsPt30_jesUp",1);
  oldtree->SetBranchStatus("nCleanedJetsPt30_jesDn",1);
  oldtree->SetBranchStatus("JetPt",1);
  oldtree->SetBranchStatus("JetEta",1);
  oldtree->SetBranchStatus("JetPhi",1);
  oldtree->SetBranchStatus("JetMass",1);
  oldtree->SetBranchStatus("JetPt_JESUp",1);
  oldtree->SetBranchStatus("JetPt_JESDown",1);
  oldtree->SetBranchStatus("GenHMass",1);
  oldtree->SetBranchStatus("GenHPt",1);
  oldtree->SetBranchStatus("GenHRapidity",1);
  oldtree->SetBranchStatus("GenZ1Flav",1);
  oldtree->SetBranchStatus("GenZ2Flav",1);
  oldtree->SetBranchStatus("GenZ1Mass",1);
  oldtree->SetBranchStatus("GenZ2Mass",1);
  oldtree->SetBranchStatus("GenZ1Pt",1);
  oldtree->SetBranchStatus("GenZ2Pt",1);
  oldtree->SetBranchStatus("GenLep1Pt",1);
  oldtree->SetBranchStatus("GenLep2Pt",1);
  oldtree->SetBranchStatus("GenLep3Pt",1);
  oldtree->SetBranchStatus("GenLep4Pt",1);
  oldtree->SetBranchStatus("GenLep1Eta",1);
  oldtree->SetBranchStatus("GenLep2Eta",1);
  oldtree->SetBranchStatus("GenLep3Eta",1);
  oldtree->SetBranchStatus("GenLep4Eta",1);
  oldtree->SetBranchStatus("GenLep1Phi",1);
  oldtree->SetBranchStatus("GenLep2Phi",1);
  oldtree->SetBranchStatus("GenLep3Phi",1);
  oldtree->SetBranchStatus("GenLep4Phi",1);
  oldtree->SetBranchStatus("GenLep1Id",1);
  oldtree->SetBranchStatus("GenLep2Id",1);
  oldtree->SetBranchStatus("GenLep3Id",1);
  oldtree->SetBranchStatus("GenLep4Id",1);
  oldtree->SetBranchStatus("GenLep1Iso",1);
  oldtree->SetBranchStatus("GenLep2Iso",1);
  oldtree->SetBranchStatus("GenLep3Iso",1);
  oldtree->SetBranchStatus("GenLep4Iso",1);
  oldtree->SetBranchStatus("GenAssocLep1Id",1);
  oldtree->SetBranchStatus("GenAssocLep2Id",1);
  oldtree->SetBranchStatus("GenAssocLep1Pt",1);
  oldtree->SetBranchStatus("GenAssocLep2Pt",1);
  oldtree->SetBranchStatus("Gencosthetastar",1);
  oldtree->SetBranchStatus("GenhelcosthetaZ1",1);
  oldtree->SetBranchStatus("GenhelcosthetaZ2",1);
  oldtree->SetBranchStatus("Genhelphi",1);
  oldtree->SetBranchStatus("GenphistarZ1",1);
  oldtree->SetBranchStatus("GenCleanedJetPt",1);
  oldtree->SetBranchStatus("GenCleanedJetMass",1);
  oldtree->SetBranchStatus("GenCleanedJetEta",1);
  oldtree->SetBranchStatus("GenCleanedJetPhi",1);
  oldtree->SetBranchStatus("GenCleanedJetRapidity",1);
  oldtree->SetBranchStatus("nCleanedGenJet",1);
  oldtree->SetBranchStatus("LHEPDFScale",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR1_muF1",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR1_muF2",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR1_muF0p5",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR2_muF1",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR2_muF2",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR2_muF0p5",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR0p5_muF1",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR0p5_muF2",1);
  oldtree->SetBranchStatus("LHEweight_QCDscale_muR0p5_muF0p5",1);
  oldtree->SetBranchStatus("LHEweight_PDFVariation_Up",1);
  oldtree->SetBranchStatus("LHEweight_PDFVariation_Dn",1);
  oldtree->SetBranchStatus("LHEweight_AsMZ_Up",1);
  oldtree->SetBranchStatus("LHEweight_AsMZ_Dn",1);
  oldtree->SetBranchStatus("LepPt",1);
  oldtree->SetBranchStatus("LepEta",1);
  oldtree->SetBranchStatus("LepPhi",1);
  oldtree->SetBranchStatus("LepLepId",1);
  oldtree->SetBranchStatus("ExtraLepPt",1);
  oldtree->SetBranchStatus("ExtraLepPhi",1);
  oldtree->SetBranchStatus("ExtraLepEta",1);
  oldtree->SetBranchStatus("ExtraLepLepId",1);
  oldtree->SetBranchStatus("PUWeight",1);
  oldtree->SetBranchStatus("genHEPMCweight",1);
  oldtree->SetBranchStatus("genHEPMCweight_NNLO",1);
  oldtree->SetBranchStatus("overallEventWeight",1);
  oldtree->SetBranchStatus("L1prefiringWeight",1);
  oldtree->SetBranchStatus("dataMCWeight",1);
  oldtree->SetBranchStatus("trigEffWeight",1);
  oldtree->SetBranchStatus("GENfinalState",1);
  oldtree->SetBranchStatus("passedFiducialSelection_bbf",1);
  oldtree->SetBranchStatus("GENlep_pt",1);
  oldtree->SetBranchStatus("GENlep_eta",1);
  oldtree->SetBranchStatus("GENlep_phi",1);
  oldtree->SetBranchStatus("GENlep_mass",1);
  oldtree->SetBranchStatus("GENlep_id",1);
  oldtree->SetBranchStatus("GENlep_status",1);
  oldtree->SetBranchStatus("GENlep_MomId",1);
  oldtree->SetBranchStatus("GENlep_MomMomId",1);
  oldtree->SetBranchStatus("GENlep_Hindex",1);
  oldtree->SetBranchStatus("GENlep_isoCH",1);
  oldtree->SetBranchStatus("GENlep_isoNH",1);
  oldtree->SetBranchStatus("GENlep_isoPhot",1);
  oldtree->SetBranchStatus("GENlep_RelIso",1);
  oldtree->SetBranchStatus("GENH_pt",1);
  oldtree->SetBranchStatus("GENH_eta",1);
  oldtree->SetBranchStatus("GENH_phi",1);
  oldtree->SetBranchStatus("GENH_mass",1);
  oldtree->SetBranchStatus("GENmass4l",1);
  oldtree->SetBranchStatus("GENmass4mu",1);
  oldtree->SetBranchStatus("GENmass4e",1);
  oldtree->SetBranchStatus("GENmass2e2mu",1);
  oldtree->SetBranchStatus("GENpT4l",1);
  oldtree->SetBranchStatus("GENeta4l",1);
  oldtree->SetBranchStatus("GENrapidity4l",1);
  oldtree->SetBranchStatus("GENcosTheta1",1);
  oldtree->SetBranchStatus("GENcosTheta2",1);
  oldtree->SetBranchStatus("GENcosThetaStar",1);
  oldtree->SetBranchStatus("GENPhi",1);
  oldtree->SetBranchStatus("GENPhi1",1);
  oldtree->SetBranchStatus("GENMH",1);
  oldtree->SetBranchStatus("GENZ_pt",1);
  oldtree->SetBranchStatus("GENZ_eta",1);
  oldtree->SetBranchStatus("GENZ_phi",1);
  oldtree->SetBranchStatus("GENZ_mass",1);
  oldtree->SetBranchStatus("GENZ_DaughtersId",1);
  oldtree->SetBranchStatus("GENZ_MomId",1);
  oldtree->SetBranchStatus("GENmassZ1",1);
  oldtree->SetBranchStatus("GENmassZ2",1);
  oldtree->SetBranchStatus("GENpTZ1",1);
  oldtree->SetBranchStatus("GENpTZ2",1);
  oldtree->SetBranchStatus("GENdPhiZZ",1);
  oldtree->SetBranchStatus("GENmassZZ",1);
  oldtree->SetBranchStatus("GENpTZZ",1);
  if(prod_mode == "ggH125") oldtree->SetBranchStatus("ggH_NNLOPS_weight",1); // Additional entry for the weight in case of ggH

  //// candTree_failed
  // Deactivate all branches
  oldtree_failed->SetBranchStatus("*",0);
  // Activate some branches only: our skim
  oldtree_failed->SetBranchStatus("EventNumber",1);
  oldtree_failed->SetBranchStatus("xsec",1);
  oldtree_failed->SetBranchStatus("GenHMass",1);
  oldtree_failed->SetBranchStatus("GenHPt",1);
  oldtree_failed->SetBranchStatus("GenHRapidity",1);
  oldtree_failed->SetBranchStatus("GenZ1Flav",1);
  oldtree_failed->SetBranchStatus("GenZ2Flav",1);
  oldtree_failed->SetBranchStatus("GenZ1Mass",1);
  oldtree_failed->SetBranchStatus("GenZ2Mass",1);
  oldtree_failed->SetBranchStatus("GenZ1Pt",1);
  oldtree_failed->SetBranchStatus("GenZ2Pt",1);
  oldtree_failed->SetBranchStatus("GenLep1Pt",1);
  oldtree_failed->SetBranchStatus("GenLep2Pt",1);
  oldtree_failed->SetBranchStatus("GenLep3Pt",1);
  oldtree_failed->SetBranchStatus("GenLep4Pt",1);
  oldtree_failed->SetBranchStatus("GenLep1Eta",1);
  oldtree_failed->SetBranchStatus("GenLep2Eta",1);
  oldtree_failed->SetBranchStatus("GenLep3Eta",1);
  oldtree_failed->SetBranchStatus("GenLep4Eta",1);
  oldtree_failed->SetBranchStatus("GenLep1Phi",1);
  oldtree_failed->SetBranchStatus("GenLep2Phi",1);
  oldtree_failed->SetBranchStatus("GenLep3Phi",1);
  oldtree_failed->SetBranchStatus("GenLep4Phi",1);
  oldtree_failed->SetBranchStatus("GenLep1Id",1);
  oldtree_failed->SetBranchStatus("GenLep2Id",1);
  oldtree_failed->SetBranchStatus("GenLep3Id",1);
  oldtree_failed->SetBranchStatus("GenLep4Id",1);
  oldtree_failed->SetBranchStatus("GenLep1Iso",1);
  oldtree_failed->SetBranchStatus("GenLep2Iso",1);
  oldtree_failed->SetBranchStatus("GenLep3Iso",1);
  oldtree_failed->SetBranchStatus("GenLep4Iso",1);
  oldtree_failed->SetBranchStatus("GenAssocLep1Id",1);
  oldtree_failed->SetBranchStatus("GenAssocLep2Id",1);
  oldtree_failed->SetBranchStatus("GenAssocLep1Pt",1);
  oldtree_failed->SetBranchStatus("GenAssocLep2Pt",1);
  oldtree_failed->SetBranchStatus("Gencosthetastar",1);
  oldtree_failed->SetBranchStatus("GenhelcosthetaZ1",1);
  oldtree_failed->SetBranchStatus("GenhelcosthetaZ2",1);
  oldtree_failed->SetBranchStatus("Genhelphi",1);
  oldtree_failed->SetBranchStatus("GenphistarZ1",1);
  oldtree_failed->SetBranchStatus("GenCleanedJetPt",1);
  oldtree_failed->SetBranchStatus("GenCleanedJetMass",1);
  oldtree_failed->SetBranchStatus("GenCleanedJetEta",1);
  oldtree_failed->SetBranchStatus("GenCleanedJetPhi",1);
  oldtree_failed->SetBranchStatus("GenCleanedJetRapidity",1);
  oldtree_failed->SetBranchStatus("nCleanedGenJet",1);
  oldtree_failed->SetBranchStatus("LHEPDFScale",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR1_muF1",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR1_muF2",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR1_muF0p5",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR2_muF1",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR2_muF2",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR2_muF0p5",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR0p5_muF1",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR0p5_muF2",1);
  oldtree_failed->SetBranchStatus("LHEweight_QCDscale_muR0p5_muF0p5",1);
  oldtree_failed->SetBranchStatus("LHEweight_PDFVariation_Up",1);
  oldtree_failed->SetBranchStatus("LHEweight_PDFVariation_Dn",1);
  oldtree_failed->SetBranchStatus("LHEweight_AsMZ_Up",1);
  oldtree_failed->SetBranchStatus("LHEweight_AsMZ_Dn",1);
  oldtree_failed->SetBranchStatus("PUWeight",1);
  oldtree_failed->SetBranchStatus("genHEPMCweight",1);
  oldtree_failed->SetBranchStatus("genHEPMCweight_NNLO",1);
  oldtree_failed->SetBranchStatus("GENfinalState",1);
  oldtree_failed->SetBranchStatus("passedFiducialSelection_bbf",1);
  oldtree_failed->SetBranchStatus("GENlep_pt",1);
  oldtree_failed->SetBranchStatus("GENlep_eta",1);
  oldtree_failed->SetBranchStatus("GENlep_phi",1);
  oldtree_failed->SetBranchStatus("GENlep_mass",1);
  oldtree_failed->SetBranchStatus("GENlep_id",1);
  oldtree_failed->SetBranchStatus("GENlep_status",1);
  oldtree_failed->SetBranchStatus("GENlep_MomId",1);
  oldtree_failed->SetBranchStatus("GENlep_MomMomId",1);
  oldtree_failed->SetBranchStatus("GENlep_Hindex",1);
  oldtree_failed->SetBranchStatus("GENlep_isoCH",1);
  oldtree_failed->SetBranchStatus("GENlep_isoNH",1);
  oldtree_failed->SetBranchStatus("GENlep_isoPhot",1);
  oldtree_failed->SetBranchStatus("GENlep_RelIso",1);
  oldtree_failed->SetBranchStatus("GENH_pt",1);
  oldtree_failed->SetBranchStatus("GENH_eta",1);
  oldtree_failed->SetBranchStatus("GENH_phi",1);
  oldtree_failed->SetBranchStatus("GENH_mass",1);
  oldtree_failed->SetBranchStatus("GENmass4l",1);
  oldtree_failed->SetBranchStatus("GENmass4mu",1);
  oldtree_failed->SetBranchStatus("GENmass4e",1);
  oldtree_failed->SetBranchStatus("GENmass2e2mu",1);
  oldtree_failed->SetBranchStatus("GENpT4l",1);
  oldtree_failed->SetBranchStatus("GENeta4l",1);
  oldtree_failed->SetBranchStatus("GENrapidity4l",1);
  oldtree_failed->SetBranchStatus("GENcosTheta1",1);
  oldtree_failed->SetBranchStatus("GENcosTheta2",1);
  oldtree_failed->SetBranchStatus("GENcosThetaStar",1);
  oldtree_failed->SetBranchStatus("GENPhi",1);
  oldtree_failed->SetBranchStatus("GENPhi1",1);
  oldtree_failed->SetBranchStatus("GENMH",1);
  oldtree_failed->SetBranchStatus("GENZ_pt",1);
  oldtree_failed->SetBranchStatus("GENZ_eta",1);
  oldtree_failed->SetBranchStatus("GENZ_phi",1);
  oldtree_failed->SetBranchStatus("GENZ_mass",1);
  oldtree_failed->SetBranchStatus("GENZ_DaughtersId",1);
  oldtree_failed->SetBranchStatus("GENZ_MomId",1);
  oldtree_failed->SetBranchStatus("GENmassZ1",1);
  oldtree_failed->SetBranchStatus("GENmassZ2",1);
  oldtree_failed->SetBranchStatus("GENpTZ1",1);
  oldtree_failed->SetBranchStatus("GENpTZ2",1);
  oldtree_failed->SetBranchStatus("GENdPhiZZ",1);
  oldtree_failed->SetBranchStatus("GENmassZZ",1);
  oldtree_failed->SetBranchStatus("GENpTZZ",1);

  if(prod_mode == "ggH125") oldtree_failed->SetBranchStatus("ggH_NNLOPS_weight",1); // Additional entry for the weight in case of ggH

  // Copy branches in the new file
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path = Form("%s/%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  TFile *newfile = new TFile(new_full_path.Data(),"RECREATE");
  TTree *newtree = (TTree*) oldtree->CloneTree(0);
  newtree->CopyEntries(oldtree);
  newtree->Write(); // Write candTree
  TTree *newtree_failed = (TTree*) oldtree_failed->CloneTree(0);
  newtree_failed->CopyEntries(oldtree_failed);
  newtree_failed->Write(); // Write candTree_failed
  hCounters->Write(); // Write Counters
  newfile->Close();

  bool t_failed;
  add(input_dir, year, prod_mode);
  add(input_dir, year, prod_mode, t_failed = false);

  // Merge together into a single TTree. Useful for efficiencies calculation.
  TFile* inputfile = TFile::Open(new_full_path.Data(), "READ");
  TTree* tree1 = (TTree*) inputfile->Get("candTree");
  TTree* tree2 = (TTree*) inputfile->Get("candTree_failed");
  TH1F* cnts = (TH1F*) inputfile->Get("Counters");

  TString merged_name = Form("%s_mergedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString merged_path = Form("%s/%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),merged_name.Data());

  TFile* mergedTTree = new TFile(merged_path.Data(), "RECREATE");
  TList* alist = new TList;

  alist->Add(tree1);
  alist->Add(tree2);

  TTree *newtree_single = TTree::MergeTrees(alist);
  newtree_single->SetName("fullTree");
  newtree_single->Write();
  cnts->Write();
  mergedTTree->Close();
  inputfile->Close();

  return 0;
}

