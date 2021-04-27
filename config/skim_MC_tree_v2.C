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

TFile *gConstant_g4 = new TFile("gConstant_HZZ2e2mu_g4.root");
TSpline *spline_g4 = (TSpline*) gConstant_g4->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");

TFile *gConstant_g2 = new TFile("gConstant_HZZ2e2mu_g2.root");
TSpline *spline_g2 = (TSpline*) gConstant_g2->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");

TFile *gConstant_L1 = new TFile("gConstant_HZZ2e2mu_L1.root");
TSpline *spline_L1 = (TSpline*) gConstant_L1->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1");

TFile *gConstant_L1Zgs = new TFile("gConstant_HZZ2e2mu_L1Zgs.root");
TSpline *spline_L1Zgs = (TSpline*) gConstant_L1Zgs->Get("sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs");

vector<TString> jes_name{
  "Total",
  "Abs",
  "Abs_year",
  "BBEC1",
  "BBEC1_year",
  "EC2",
  "EC2_year",
  "FlavQCD",
  "HF",
  "HF_year",
  "RelBal",
  "RelSample_year"
};

float mass_lep(int flavour){
  if((abs(flavour)) == 11) return 0.0005109989461;
  else if ((abs(flavour)) == 13) return 0.1056583745;
  else if ((abs(flavour)) == 15) return 1.77686;
  // else if ((abs(flavour)) == 0) return 0;
  return 0;
}


// pair<vector<TLorentzVector>,vector<Short_t>> sort(vector<TLorentzVector> lep,vector<Short_t>id){
//   vector<TLorentzVector> sortedLep;
//   vector<Short_t> sortedId;
//   for (int k=0;k<4;k++){
//     int max = -1;
//     float maxPt = 0;
//     for (int i=0;i<lep.size();i++){
//       if (lep[i].Pt() > maxPt) {
//         max = i;
//         maxPt = lep[i].Pt();
//       }
//     }
//     sortedLep.push_back(lep[max]);
//     sortedId.push_back(id[max]);
//     lep.erase(lep.begin()+max);
//     id.erase(id.begin()+max);
//   }
//   return make_pair(sortedLep,sortedId);
// }
// pair<vector<TLorentzVector>,vector<Short_t>> sort_2(vector<TLorentzVector> lep,vector<Short_t>id){
//   vector<TLorentzVector> sortedLep;
//   vector<Short_t> sortedId;
//   for (int k=0;k<2;k++){
//     int max = -1;
//     float maxPt = 0;
//     for (int i=0;i<lep.size();i++){
//       if (lep[i].Pt() > maxPt) {
//         max = i;
//         maxPt = lep[i].Pt();
//       }
//     }
//     sortedLep.push_back(lep[max]);
//     sortedId.push_back(id[max]);
//     lep.erase(lep.begin()+max);
//   }
//   return make_pair(sortedLep,sortedId);
// }


//------------------------------------------------------------------
void add(TString input_dir, TString year, TString prod_mode, TString process, bool t_failed=true, bool flag_tmp_2017=false){
  // Add additional branches
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path;
  if(process!="AC") new_full_path = Form("%s/%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  else new_full_path = Form("%s/AC%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  cout << new_full_path << endl;
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
  // Gen-variables
  float GenLep1Pt,GenLep2Pt,GenLep3Pt,GenLep4Pt,GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta,GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi,GenZ1Flav,GenZ2Flav,
        GenZ1Mass,GenZ2Mass,GenLep1Iso,GenLep2Iso,GenLep3Iso,GenLep4Iso;
  Float_t p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz2_1_ghz4_1_JHUGen;
  Float_t p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;
  Float_t p_GEN_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen, p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen;
  Float_t _GEN_Dcp, _GEN_D0m, _GEN_D0hp, _GEN_Dint, _GEN_DL1, _GEN_DL1int, _GEN_DL1Zg, _GEN_DL1Zgint;
  Float_t GENmass4l,GENpT4l,GENeta4l,GENphi4l;
  Short_t GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  Short_t _GENnjets_pt30_eta2p5;
  Short_t _GENnjets_pt30_eta4p7;
  Float_t _GENpTj1, _GENpTj2, _GENpTHj, _GENpTHjj, _GENmHj, _GENmHjj, _GENmjj, _GENdetajj, _GENdphijj, _GENabsdetajj, _GENabsdphijj;
  Float_t _GENTCj1, _GENTBj1;
  Float_t _GENTCj, _GENTBj, _GENTCjmax, _GENTBjmax;
  bool _passedFullSelection, passedFiducialSelection_bbf;
  // vector<float> _GenLepPtSorted,_GenLepEtaSorted,_GenLepPhiSorted;
  // vector<Short_t> _GenLepIdSorted;
  // vector<TLorentzVector> GenLepSorted;
  vector<float> *GENlep_pt = 0;
  vector<float> *GENlep_eta = 0;
  vector<float> *GENlep_phi = 0;
  vector<float> *GENlep_mass = 0;
  vector<float> *GENlep_id = 0;
  vector<float> *GENjetsPt_pt30_eta2p5 = 0;
  vector<float> *GENjetsEta_pt30_eta2p5 = 0;
  vector<float> *GENjetsPhi_pt30_eta2p5 = 0;
  vector<float> *GENjetsMass_pt30_eta2p5 = 0;
  vector<float> *GENjetsPt_pt30_eta4p7 = 0;
  vector<float> *GENjetsEta_pt30_eta4p7 = 0;
  vector<float> *GENjetsPhi_pt30_eta4p7 = 0;
  vector<float> *GENjetsMass_pt30_eta4p7 = 0;
  vector<Short_t> *GENlep_Hindex = 0;
  // TBranch *GenLepPtSorted = T->Branch("GenLepPtSorted",&_GenLepPtSorted);
  // TBranch *GenLepEtaSorted = T->Branch("GenLepEtaSorted",&_GenLepEtaSorted);
  // TBranch *GenLepPhiSorted = T->Branch("GenLepPhiSorted",&_GenLepPhiSorted);
  // TBranch *GenLepIdSorted = T->Branch("GenLepIdSorted",&_GenLepIdSorted);
  TBranch *passedFullSelection = T->Branch("passedFullSelection",&_passedFullSelection,"passedFullSelection/B");
  TBranch *GENnjets_pt30_eta2p5 = T->Branch("GENnjets_pt30_eta2p5",&_GENnjets_pt30_eta2p5,"GENnjets_pt30_eta2p5/S");
  TBranch *GENnjets_pt30_eta4p7 = T->Branch("GENnjets_pt30_eta4p7",&_GENnjets_pt30_eta4p7,"GENnjets_pt30_eta4p7/S");
  TBranch *GENpTj1 = T->Branch("GENpTj1",&_GENpTj1,"GENpTj1/F");
  TBranch *GENpTj2 = T->Branch("GENpTj2",&_GENpTj2,"GENpTj2/F");
  TBranch *GENpTHj = T->Branch("GENpTHj",&_GENpTHj,"GENpTHj/F");
  TBranch *GENpTHjj = T->Branch("GENpTHjj",&_GENpTHjj,"GENpTHjj/F");
  TBranch *GENmHj = T->Branch("GENmHj",&_GENmHj,"GENmHj/F");
  TBranch *GENmHjj = T->Branch("GENmHjj",&_GENmHjj,"GENmHjj/F");
  TBranch *GENmjj = T->Branch("GENmjj",&_GENmjj,"GENmjj/F");
  TBranch *GENdetajj = T->Branch("GENdetajj",&_GENdetajj,"GENdetajj/F");
  TBranch *GENdphijj = T->Branch("GENdphijj",&_GENdphijj,"GENdphijj/F");
  TBranch *GENabsdetajj = T->Branch("GENabsdetajj",&_GENabsdetajj,"GENabsdetajj/F");
  TBranch *GENabsdphijj = T->Branch("GENabsdphijj",&_GENabsdphijj,"GENabsdphijj/F");
  TBranch *GENTCj1 = T->Branch("GENTCj1",&_GENTCj1,"GENTCj1/F");
  TBranch *GENTBj1 = T->Branch("GENTBj1",&_GENTBj1,"GENTBj1/F");
  TBranch *GENTCjmax = T->Branch("GENTCjmax",&_GENTCjmax,"GENTCjmax/F");
  TBranch *GENTBjmax = T->Branch("GENTBjmax",&_GENTBjmax,"GENTBjmax/F");
  TBranch *GEN_Dcp = T->Branch("GEN_Dcp",&_GEN_Dcp,"GEN_Dcp/F");
  TBranch *GEN_D0m = T->Branch("GEN_D0m",&_GEN_D0m,"GEN_D0m/F");
  TBranch *GEN_Dint = T->Branch("GEN_Dint",&_GEN_Dint,"GEN_Dint/F");
  TBranch *GEN_D0hp = T->Branch("GEN_D0hp",&_GEN_D0hp,"GEN_D0hp/F");
  TBranch *GEN_DL1 = T->Branch("GEN_DL1",&_GEN_DL1,"GEN_DL1/F");
  TBranch *GEN_DL1int = T->Branch("GEN_DL1int",&_GEN_DL1int,"GEN_DL1int/F");
  TBranch *GEN_DL1Zg = T->Branch("GEN_DL1Zg",&_GEN_DL1Zg,"GEN_DL1Zg/F");
  TBranch *GEN_DL1Zgint = T->Branch("GEN_DL1Zgint",&_GEN_DL1Zgint,"GEN_DL1Zgint/F");

  if(process=="signal" || process=="AC"){ // Bkgs don't store gen-level information
    T->SetBranchAddress("GENlep_pt",&GENlep_pt);
    T->SetBranchAddress("GENlep_eta",&GENlep_eta);
    T->SetBranchAddress("GENlep_phi",&GENlep_phi);
    T->SetBranchAddress("GENlep_mass",&GENlep_mass);
    T->SetBranchAddress("GENlep_id",&GENlep_id);
    T->SetBranchAddress("GENjetsPt_pt30_eta2p5",&GENjetsPt_pt30_eta2p5);
    T->SetBranchAddress("GENjetsEta_pt30_eta2p5",&GENjetsEta_pt30_eta2p5);
    T->SetBranchAddress("GENjetsPhi_pt30_eta2p5",&GENjetsPhi_pt30_eta2p5);
    T->SetBranchAddress("GENjetsMass_pt30_eta2p5",&GENjetsMass_pt30_eta2p5);
    T->SetBranchAddress("GENjetsPt_pt30_eta4p7",&GENjetsPt_pt30_eta4p7);
    T->SetBranchAddress("GENjetsEta_pt30_eta4p7",&GENjetsEta_pt30_eta4p7);
    T->SetBranchAddress("GENjetsPhi_pt30_eta4p7",&GENjetsPhi_pt30_eta4p7);
    T->SetBranchAddress("GENjetsMass_pt30_eta4p7",&GENjetsMass_pt30_eta4p7);
    T->SetBranchAddress("GENmass4l",&GENmass4l);
    T->SetBranchAddress("GENpT4l",&GENpT4l);
    T->SetBranchAddress("GENeta4l",&GENeta4l);
    T->SetBranchAddress("GENphi4l",&GENphi4l);
    T->SetBranchAddress("GENlep_Hindex",&GENlep_Hindex);
    T->SetBranchAddress("passedFiducialSelection_bbf",&passedFiducialSelection_bbf);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz2_1_ghz4_1_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz2_1_ghz4_1_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",&p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen);
  }

  // Reco-variables and Gen-Reco-matching variables
  float _ZZy,ZZPt,ZZEta,ZZPhi,ZZMass;
  // Short_t nCleanedJetsPt30;
  // Short_t nCleanedJetsPt30_jesUp_Total,nCleanedJetsPt30_jesDn_Total;
  // Short_t nCleanedJetsPt30_jesUp_Abs,nCleanedJetsPt30_jesDn_Abs;
  // Short_t nCleanedJetsPt30_jesUp_Abs_year,nCleanedJetsPt30_jesDn_Abs_year;
  // Short_t nCleanedJetsPt30_jesUp_BBEC1,nCleanedJetsPt30_jesDn_BBEC1;
  // Short_t nCleanedJetsPt30_jesUp_BBEC1_year,nCleanedJetsPt30_jesDn_BBEC1_year;
  // Short_t nCleanedJetsPt30_jesUp_EC2,nCleanedJetsPt30_jesDn_EC2;
  // Short_t nCleanedJetsPt30_jesUp_EC2_year,nCleanedJetsPt30_jesDn_EC2_year;
  // Short_t nCleanedJetsPt30_jesUp_FlavQCD,nCleanedJetsPt30_jesDn_FlavQCD;
  // Short_t nCleanedJetsPt30_jesUp_HF,nCleanedJetsPt30_jesDn_HF;
  // Short_t nCleanedJetsPt30_jesUp_HF_year,nCleanedJetsPt30_jesDn_HF_year;
  // Short_t nCleanedJetsPt30_jesUp_RelBal,nCleanedJetsPt30_jesDn_RelBal;
  // Short_t nCleanedJetsPt30_jesUp_RelSample_year,nCleanedJetsPt30_jesDn_RelSample_year;
  Short_t _njets_pt30_eta2p5, _njets_pt30_eta4p7;
  // Short_t _njets_pt30_eta2p5_jesup_Total, _njets_pt30_eta2p5_jesdn_Total;
  // Short_t _njets_pt30_eta2p5_jesup_Abs, _njets_pt30_eta2p5_jesdn_Abs;
  // Short_t _njets_pt30_eta2p5_jesup_Abs_year, _njets_pt30_eta2p5_jesdn_Abs_year;
  // Short_t _njets_pt30_eta2p5_jesup_BBEC1, _njets_pt30_eta2p5_jesdn_BBEC1;
  // Short_t _njets_pt30_eta2p5_jesup_BBEC1_year, _njets_pt30_eta2p5_jesdn_BBEC1_year;
  // Short_t _njets_pt30_eta2p5_jesup_EC2, _njets_pt30_eta2p5_jesdn_EC2;
  // Short_t _njets_pt30_eta2p5_jesup_EC2_year, _njets_pt30_eta2p5_jesdn_EC2_year;
  // Short_t _njets_pt30_eta2p5_jesup_FlavQCD, _njets_pt30_eta2p5_jesdn_FlavQCD;
  // Short_t _njets_pt30_eta2p5_jesup_HF, _njets_pt30_eta2p5_jesdn_HF;
  // Short_t _njets_pt30_eta2p5_jesup_HF_year, _njets_pt30_eta2p5_jesdn_HF_year;
  // Short_t _njets_pt30_eta2p5_jesup_RelBal, _njets_pt30_eta2p5_jesdn_RelBal;
  // Short_t _njets_pt30_eta2p5_jesup_RelSample_year, _njets_pt30_eta2p5_jesdn_RelSample_year;
  Float_t _pTj1, _pTj2, _pTHj, _pTHjj, _mHj, _mHjj, _mjj, _detajj, _dphijj, _absdetajj, _absdphijj;
  Float_t _pTj1_eta4p7;
  // Float_t _pTj1_jesup, _pTj1_jesdn, _pTj2_jesup, _pTj2_jesdn;
  // Float_t _pTHj_jesup, _pTHj_jesdn, _pTHjj_jesup, _pTHjj_jesdn;
  // Float_t _mHj_jesup, _mHj_jesdn, _mHjj_jesup, _mHjj_jesdn, _mjj_jesup, _mjj_jesdn;
  // Float_t _detajj_jesup, _detajj_jesdn, _dphijj_jesup, _dphijj_jesdn;
  // Float_t _absdetajj_jesup, _absdetajj_jesdn, _absdphijj_jesup, _absdphijj_jesdn;
  Float_t Mj1, ETAj1, PHIj1, Mj2, ETAj2, PHIj2;
  // Float_t Mj1_jesup, ETAj1_jesup, PHIj1_jesup, Mj2_jesup, ETAj2_jesup, PHIj2_jesup;
  // Float_t Mj1_jesdn, ETAj1_jesdn, PHIj1_jesdn, Mj2_jesdn, ETAj2_jesdn, PHIj2_jesdn;
  Float_t _TCj, _TCj1, _TCjmax, _TBj1, _TBj, _TBjmax;
  Float_t tc, tcj, yj;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_GG_SIG_ghg2_1_ghz4_1_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen;
  Float_t p_GG_SIG_ghg2_1_ghz2_1_JHUGen, p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen;
  Float_t _Dcp, _D0m, _D0hp, _Dint, _DL1, _DL1int, _DL1Zg, _DL1Zgint;
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
  vector<float> *JetMass = 0;
  vector<float> *JetPhi = 0;
  // vector<float> *JetPt_JESUp_Total = 0;
  // vector<float> *JetPt_JESUp_Abs = 0;
  // vector<float> *JetPt_JESUp_Abs_year = 0;
  // vector<float> *JetPt_JESUp_BBEC1 = 0;
  // vector<float> *JetPt_JESUp_BBEC1_year = 0;
  // vector<float> *JetPt_JESUp_EC2 = 0;
  // vector<float> *JetPt_JESUp_EC2_year = 0;
  // vector<float> *JetPt_JESUp_FlavQCD = 0;
  // vector<float> *JetPt_JESUp_HF = 0;
  // vector<float> *JetPt_JESUp_HF_year = 0;
  // vector<float> *JetPt_JESUp_RelBal = 0;
  // vector<float> *JetPt_JESUp_RelSample_year = 0;
  // vector<float> *JetPt_JESDown_Total = 0;
  // vector<float> *JetPt_JESDown_Abs = 0;
  // vector<float> *JetPt_JESDown_Abs_year = 0;
  // vector<float> *JetPt_JESDown_BBEC1 = 0;
  // vector<float> *JetPt_JESDown_BBEC1_year = 0;
  // vector<float> *JetPt_JESDown_EC2 = 0;
  // vector<float> *JetPt_JESDown_EC2_year = 0;
  // vector<float> *JetPt_JESDown_FlavQCD = 0;
  // vector<float> *JetPt_JESDown_HF = 0;
  // vector<float> *JetPt_JESDown_HF_year = 0;
  // vector<float> *JetPt_JESDown_RelBal = 0;
  // vector<float> *JetPt_JESDown_RelSample_year = 0;
  TBranch *ZZy = T->Branch("ZZy",&_ZZy,"ZZy/F");
  TBranch *lep_genindex = T->Branch("lep_genindex",&_lep_genindex);
  TBranch *lep_Hindex = T->Branch("lep_Hindex",&_lep_Hindex);
  TBranch *njets_pt30_eta2p5 = T->Branch("njets_pt30_eta2p5",&_njets_pt30_eta2p5,"njets_pt30_eta2p5/S");
  TBranch *njets_pt30_eta4p7 = T->Branch("njets_pt30_eta4p7",&_njets_pt30_eta4p7,"njets_pt30_eta4p7/S");
  // TBranch *njets_pt30_eta2p5_jesup_Total = T->Branch("njets_pt30_eta2p5_jesup_Total",&_njets_pt30_eta2p5_jesup_Total,"_njets_pt30_eta2p5_jesup_Total/S");
  // TBranch *njets_pt30_eta2p5_jesup_Abs = T->Branch("njets_pt30_eta2p5_jesup_Abs",&_njets_pt30_eta2p5_jesup_Abs,"_njets_pt30_eta2p5_jesup_Abs/S");
  // TBranch *njets_pt30_eta2p5_jesup_Abs_year = T->Branch("njets_pt30_eta2p5_jesup_Abs_year",&_njets_pt30_eta2p5_jesup_Abs_year,"_njets_pt30_eta2p5_jesup_Abs_year/S");
  // TBranch *njets_pt30_eta2p5_jesup_BBEC1 = T->Branch("njets_pt30_eta2p5_jesup_BBEC1",&_njets_pt30_eta2p5_jesup_BBEC1,"_njets_pt30_eta2p5_jesup_BBEC1/S");
  // TBranch *njets_pt30_eta2p5_jesup_BBEC1_year = T->Branch("njets_pt30_eta2p5_jesup_BBEC1_year",&_njets_pt30_eta2p5_jesup_BBEC1_year,"_njets_pt30_eta2p5_jesup_BBEC1_year/S");
  // TBranch *njets_pt30_eta2p5_jesup_EC2 = T->Branch("njets_pt30_eta2p5_jesup_EC2",&_njets_pt30_eta2p5_jesup_EC2,"_njets_pt30_eta2p5_jesup_EC2/S");
  // TBranch *njets_pt30_eta2p5_jesup_EC2_year = T->Branch("njets_pt30_eta2p5_jesup_EC2_year",&_njets_pt30_eta2p5_jesup_EC2_year,"_njets_pt30_eta2p5_jesup_EC2_year/S");
  // TBranch *njets_pt30_eta2p5_jesup_FlavQCD = T->Branch("njets_pt30_eta2p5_jesup_FlavQCD",&_njets_pt30_eta2p5_jesup_FlavQCD,"_njets_pt30_eta2p5_jesup_FlavQCD/S");
  // TBranch *njets_pt30_eta2p5_jesup_HF = T->Branch("njets_pt30_eta2p5_jesup_HF",&_njets_pt30_eta2p5_jesup_HF,"_njets_pt30_eta2p5_jesup_HF/S");
  // TBranch *njets_pt30_eta2p5_jesup_HF_year = T->Branch("njets_pt30_eta2p5_jesup_HF_year",&_njets_pt30_eta2p5_jesup_HF_year,"_njets_pt30_eta2p5_jesup_HF_year/S");
  // TBranch *njets_pt30_eta2p5_jesup_RelBal = T->Branch("njets_pt30_eta2p5_jesup_RelBal",&_njets_pt30_eta2p5_jesup_RelBal,"_njets_pt30_eta2p5_jesup_RelBal/S");
  // TBranch *njets_pt30_eta2p5_jesup_RelSample_year = T->Branch("njets_pt30_eta2p5_jesup_RelSample_year",&_njets_pt30_eta2p5_jesup_RelSample_year,"_njets_pt30_eta2p5_jesup_RelSample_year/S");
  // TBranch *njets_pt30_eta2p5_jesdn_Total = T->Branch("njets_pt30_eta2p5_jesdn_Total",&_njets_pt30_eta2p5_jesdn_Total,"_njets_pt30_eta2p5_jesdn_Total/S");
  // TBranch *njets_pt30_eta2p5_jesdn_Abs = T->Branch("njets_pt30_eta2p5_jesdn_Abs",&_njets_pt30_eta2p5_jesdn_Abs,"_njets_pt30_eta2p5_jesdn_Abs/S");
  // TBranch *njets_pt30_eta2p5_jesdn_Abs_year = T->Branch("njets_pt30_eta2p5_jesdn_Abs_year",&_njets_pt30_eta2p5_jesdn_Abs_year,"_njets_pt30_eta2p5_jesdn_Abs_year/S");
  // TBranch *njets_pt30_eta2p5_jesdn_BBEC1 = T->Branch("njets_pt30_eta2p5_jesdn_BBEC1",&_njets_pt30_eta2p5_jesdn_BBEC1,"_njets_pt30_eta2p5_jesdn_BBEC1/S");
  // TBranch *njets_pt30_eta2p5_jesdn_BBEC1_year = T->Branch("njets_pt30_eta2p5_jesdn_BBEC1_year",&_njets_pt30_eta2p5_jesdn_BBEC1_year,"_njets_pt30_eta2p5_jesdn_BBEC1_year/S");
  // TBranch *njets_pt30_eta2p5_jesdn_EC2 = T->Branch("njets_pt30_eta2p5_jesdn_EC2",&_njets_pt30_eta2p5_jesdn_EC2,"_njets_pt30_eta2p5_jesdn_EC2/S");
  // TBranch *njets_pt30_eta2p5_jesdn_EC2_year = T->Branch("njets_pt30_eta2p5_jesdn_EC2_year",&_njets_pt30_eta2p5_jesdn_EC2_year,"_njets_pt30_eta2p5_jesdn_EC2_year/S");
  // TBranch *njets_pt30_eta2p5_jesdn_FlavQCD = T->Branch("njets_pt30_eta2p5_jesdn_FlavQCD",&_njets_pt30_eta2p5_jesdn_FlavQCD,"_njets_pt30_eta2p5_jesdn_FlavQCD/S");
  // TBranch *njets_pt30_eta2p5_jesdn_HF = T->Branch("njets_pt30_eta2p5_jesdn_HF",&_njets_pt30_eta2p5_jesdn_HF,"_njets_pt30_eta2p5_jesdn_HF/S");
  // TBranch *njets_pt30_eta2p5_jesdn_HF_year = T->Branch("njets_pt30_eta2p5_jesdn_HF_year",&_njets_pt30_eta2p5_jesdn_HF_year,"_njets_pt30_eta2p5_jesdn_HF_year/S");
  // TBranch *njets_pt30_eta2p5_jesdn_RelBal = T->Branch("njets_pt30_eta2p5_jesdn_RelBal",&_njets_pt30_eta2p5_jesdn_RelBal,"_njets_pt30_eta2p5_jesdn_RelBal/S");
  // TBranch *njets_pt30_eta2p5_jesdn_RelSample_year = T->Branch("njets_pt30_eta2p5_jesdn_RelSample_year",&_njets_pt30_eta2p5_jesdn_RelSample_year,"_njets_pt30_eta2p5_jesdn_RelSample_year/S");
  // TBranch *njets_pt30_eta2p5_jesup = T->Branch("njets_pt30_eta2p5_jesup",&_njets_pt30_eta2p5_jesup,"_njets_pt30_eta2p5_jesup/S");
  // TBranch *njets_pt30_eta2p5_jesdn = T->Branch("njets_pt30_eta2p5_jesdn",&_njets_pt30_eta2p5_jesdn,"_njets_pt30_eta2p5_jesdn/S");
  TBranch *pTj1 = T->Branch("pTj1",&_pTj1,"pTj1/F");
  TBranch *pTj1_eta4p7 = T->Branch("pTj1_eta4p7",&_pTj1_eta4p7,"pTj1_eta4p7/F");
  // TBranch *pTj1_jesup = T->Branch("pTj1_jesup",&_pTj1_jesup,"pTj1_jesup/F");
  // TBranch *pTj1_jesdn = T->Branch("pTj1_jesdn",&_pTj1_jesdn,"pTj1_jesdn/F");
  TBranch *pTj2 = T->Branch("pTj2",&_pTj2,"pTj2/F");
  // TBranch *pTj2_jesup = T->Branch("pTj2_jesup",&_pTj2_jesup,"pTj2_jesup/F");
  // TBranch *pTj2_jesdn = T->Branch("pTj2_jesdn",&_pTj2_jesdn,"pTj2_jesdn/F");
  TBranch *pTHj = T->Branch("pTHj",&_pTHj,"pTHj/F");
  // TBranch *pTHj_jesup = T->Branch("pTHj_jesup",&_pTHj_jesup,"pTHj_jesup/F");
  // TBranch *pTHj_jesdn = T->Branch("pTHj_jesdn",&_pTHj_jesdn,"pTHj_jesdn/F");
  TBranch *pTHjj = T->Branch("pTHjj",&_pTHjj,"pTHjj/F");
  // TBranch *pTHjj_jesup = T->Branch("pTHjj_jesup",&_pTHjj_jesup,"pTHjj_jesup/F");
  // TBranch *pTHjj_jesdn = T->Branch("pTHjj_jesdn",&_pTHjj_jesdn,"pTHjj_jesdn/F");
  TBranch *mHj = T->Branch("mHj",&_mHj,"mHj/F");
  // TBranch *mHj_jesup = T->Branch("mHj_jesup",&_mHj_jesup,"mHj_jesup/F");
  // TBranch *mHj_jesdn = T->Branch("mHj_jesdn",&_mHj_jesdn,"mHj_jesdn/F");
  TBranch *mHjj = T->Branch("mHjj",&_mHjj,"mHjj/F");
  // TBranch *mHjj_jesup = T->Branch("mHjj_jesup",&_mHjj_jesup,"mHjj_jesup/F");
  // TBranch *mHjj_jesdn = T->Branch("mHjj_jesdn",&_mHjj_jesdn,"mHjj_jesdn/F");
  TBranch *mjj = T->Branch("mjj",&_mjj,"mjj/F");
  // TBranch *mjj_jesup = T->Branch("mjj_jesup",&_mjj_jesup,"mjj_jesup/F");
  // TBranch *mjj_jesdn = T->Branch("mjj_jesdn",&_mjj_jesdn,"mjj_jesdn/F");
  TBranch *detajj = T->Branch("detajj",&_detajj,"detajj/F");
  TBranch *dphijj = T->Branch("dphijj",&_dphijj,"dphijj/F");
  // TBranch *detajj_jesup = T->Branch("detajj_jesup",&_detajj_jesup,"detajj_jesup/F");
  // TBranch *detajj_jesdn = T->Branch("detajj_jesdn",&_detajj_jesdn,"detajj_jesdn/F");
  // TBranch *dphijj_jesup = T->Branch("dphijj_jesup",&_dphijj_jesup,"dphijj_jesup/F");
  // TBranch *dphijj_jesdn = T->Branch("dphijj_jesdn",&_dphijj_jesdn,"dphijj_jesdn/F");
  TBranch *absdetajj = T->Branch("absdetajj",&_absdetajj,"absdetajj/F");
  TBranch *absdphijj = T->Branch("absdphijj",&_absdphijj,"absdphijj/F");
  // TBranch *absdetajj_jesup = T->Branch("absdetajj_jesup",&_absdetajj_jesup,"absdetajj_jesup/F");
  // TBranch *absdetajj_jesdn = T->Branch("absdetajj_jesdn",&_absdetajj_jesdn,"absdetajj_jesdn/F");
  // TBranch *absdphijj_jesup = T->Branch("absdphijj_jesup",&_absdphijj_jesup,"absdphijj_jesup/F");
  // TBranch *absdphijj_jesdn = T->Branch("absdphijj_jesdn",&_absdphijj_jesdn,"absdphijj_jesdn/F");
  TBranch *TCj = T->Branch("TCj",&_TCj,"TCj/F");
  TBranch *TCjmax = T->Branch("TCjmax",&_TCjmax,"TCjmax/F");
  TBranch *TCj1 = T->Branch("TCj1",&_TCj1,"TCj1/F");
  TBranch *TBj1 = T->Branch("TBj1",&_TBj1,"TBj1/F");
  TBranch *TBj = T->Branch("TBj",&_TBj,"TBj/F");
  TBranch *TBjmax = T->Branch("TBjmax", &_TBjmax, "TBjmax/F");
  TBranch *Dcp = T->Branch("Dcp",&_Dcp,"Dcp/F");
  TBranch *D0m = T->Branch("D0m",&_D0m,"D0m/F");
  TBranch *Dint = T->Branch("Dint",&_Dint,"Dint/F");
  TBranch *D0hp = T->Branch("D0hp",&_D0hp,"D0hp/F");
  TBranch *DL1 = T->Branch("DL1",&_DL1,"DL1/F");
  TBranch *DL1int = T->Branch("DL1int",&_DL1int,"DL1int/F");
  TBranch *DL1Zg = T->Branch("DL1Zg",&_DL1Zg,"DL1Zg/F");
  TBranch *DL1Zgint = T->Branch("DL1Zgint",&_DL1Zgint,"DL1Zgint/F");
  if (!t_failed) {
    T->SetBranchAddress("ZZMass",&ZZMass);
    T->SetBranchAddress("ZZPt",&ZZPt);
    T->SetBranchAddress("ZZEta",&ZZEta);
    T->SetBranchAddress("ZZPhi",&ZZPhi);
    T->SetBranchAddress("LepPt",&LepPt);
    T->SetBranchAddress("LepPhi",&LepPhi);
    T->SetBranchAddress("LepEta",&LepEta);
    T->SetBranchAddress("LepLepId",&LepLepId);
    T->SetBranchAddress("ExtraLepPt",&ExtraLepPt);
    T->SetBranchAddress("ExtraLepEta",&ExtraLepEta);
    T->SetBranchAddress("ExtraLepPhi",&ExtraLepPhi);
    T->SetBranchAddress("ExtraLepLepId",&ExtraLepLepId);
    // T->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
    // T->SetBranchAddress("nCleanedJetsPt30_jesUp",&nCleanedJetsPt30_jesUp);
    // T->SetBranchAddress("nCleanedJetsPt30_jesDn",&nCleanedJetsPt30_jesDn);
    T->SetBranchAddress("JetPt",&JetPt);
    T->SetBranchAddress("JetEta",&JetEta);
    T->SetBranchAddress("JetMass",&JetMass);
    T->SetBranchAddress("JetPhi",&JetPhi);
    // T->SetBranchAddress("JetPt_JESUp_Total",&JetPt_JESUp_Total);
    // T->SetBranchAddress("JetPt_JESUp_Abs",&JetPt_JESUp_Abs);
    // T->SetBranchAddress("JetPt_JESUp_Abs_year",&JetPt_JESUp_Abs_year);
    // T->SetBranchAddress("JetPt_JESUp_BBEC1",&JetPt_JESUp_BBEC1);
    // T->SetBranchAddress("JetPt_JESUp_BBEC1_year",&JetPt_JESUp_BBEC1_year);
    // T->SetBranchAddress("JetPt_JESUp_EC2",&JetPt_JESUp_EC2);
    // T->SetBranchAddress("JetPt_JESUp_EC2_year",&JetPt_JESUp_EC2_year);
    // T->SetBranchAddress("JetPt_JESUp_FlavQCD",&JetPt_JESUp_FlavQCD);
    // T->SetBranchAddress("JetPt_JESUp_HF",&JetPt_JESUp_HF);
    // T->SetBranchAddress("JetPt_JESUp_HF_year",&JetPt_JESUp_HF_year);
    // T->SetBranchAddress("JetPt_JESUp_RelBal",&JetPt_JESUp_RelBal);
    // T->SetBranchAddress("JetPt_JESUp_RelSample_year",&JetPt_JESUp_RelSample_year);
    // T->SetBranchAddress("JetPt_JESDown_Total",&JetPt_JESDown_Total);
    // T->SetBranchAddress("JetPt_JESDown_Abs",&JetPt_JESDown_Abs);
    // T->SetBranchAddress("JetPt_JESDown_Abs_year",&JetPt_JESDown_Abs_year);
    // T->SetBranchAddress("JetPt_JESDown_BBEC1",&JetPt_JESDown_BBEC1);
    // T->SetBranchAddress("JetPt_JESDown_BBEC1_year",&JetPt_JESDown_BBEC1_year);
    // T->SetBranchAddress("JetPt_JESDown_EC2",&JetPt_JESDown_EC2);
    // T->SetBranchAddress("JetPt_JESDown_EC2_year",&JetPt_JESDown_EC2_year);
    // T->SetBranchAddress("JetPt_JESDown_FlavQCD",&JetPt_JESDown_FlavQCD);
    // T->SetBranchAddress("JetPt_JESDown_HF",&JetPt_JESDown_HF);
    // T->SetBranchAddress("JetPt_JESDown_HF_year",&JetPt_JESDown_HF_year);
    // T->SetBranchAddress("JetPt_JESDown_RelBal",&JetPt_JESDown_RelBal);
    // T->SetBranchAddress("JetPt_JESDown_RelSample_year",&JetPt_JESDown_RelSample_year);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz4_1_JHUGen",&p_GG_SIG_ghg2_1_ghz4_1_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz2_1_JHUGen",&p_GG_SIG_ghg2_1_ghz2_1_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",&p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",&p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen);
    T->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",&p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen);
  }

  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);

    if(process=="signal" || process=="AC"){
      // // Sort GenLeptons
      // float GenLep1Mass = mass_lep(GenLep1Id);
      // float GenLep2Mass = mass_lep(GenLep2Id);
      // float GenLep3Mass = mass_lep(GenLep3Id);
      // float GenLep4Mass = mass_lep(GenLep4Id);
      // vector<Short_t> GenLepId {GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id};
      // TLorentzVector t1,t2,t3,t4; // Lorentz vector of the four genleptons
      // t1.SetPtEtaPhiM(GenLep1Pt,GenLep1Eta,GenLep1Phi,GenLep1Mass);
      // t2.SetPtEtaPhiM(GenLep2Pt,GenLep2Eta,GenLep2Phi,GenLep2Mass);
      // t3.SetPtEtaPhiM(GenLep3Pt,GenLep3Eta,GenLep3Phi,GenLep3Mass);
      // t4.SetPtEtaPhiM(GenLep4Pt,GenLep4Eta,GenLep4Phi,GenLep4Mass);
      // vector<TLorentzVector> GenLep {t1,t2,t3,t4};
      // pair<vector<TLorentzVector>,vector<Short_t>> sortedGenLeptons;
      // if(t3.Pt()!=0) sortedGenLeptons = sort(GenLep,GenLepId);
      // else sortedGenLeptons = sort_2(GenLep,GenLepId); /* Different function for cases in which there are just two GenLeptons, otherwise strange errors with
      //                                                     the previous function */
      // GenLepSorted = sortedGenLeptons.first;
      // _GenLepIdSorted = sortedGenLeptons.second;
      // for(int i=0;i<GenLepSorted.size();i++){
      //   _GenLepPtSorted.push_back(GenLepSorted.at(i).Pt());
      //   _GenLepEtaSorted.push_back(GenLepSorted.at(i).Eta());
      //   _GenLepPhiSorted.push_back(GenLepSorted.at(i).Phi());
      // }
      _passedFullSelection = true;
      if (t_failed) {
  	     _passedFullSelection = false;
      }

      _GENnjets_pt30_eta2p5 = GENjetsPt_pt30_eta2p5->size();
      _GENnjets_pt30_eta4p7 = GENjetsPt_pt30_eta4p7->size();

      // leading GENjet pT
      _GENpTj1 = 0;
      Mj1 = 0;
      ETAj1 = 0;
      PHIj1 = 0;
      for (unsigned int i = 0; i < GENjetsPt_pt30_eta2p5->size(); ++i){
        if(GENjetsPt_pt30_eta2p5->at(i) > _GENpTj1) {
          _GENpTj1 = GENjetsPt_pt30_eta2p5->at(i);
          Mj1 = GENjetsMass_pt30_eta2p5->at(i);
          ETAj1 = GENjetsEta_pt30_eta2p5->at(i);
          PHIj1 = GENjetsPhi_pt30_eta2p5->at(i);
        }
      }

      // sub-leading GENjet pT
      _GENpTj2 = 0;
      Mj2 = 0;
      ETAj2 = 0;
      PHIj2 = 0;
      for (unsigned int i = 0; i < GENjetsPt_pt30_eta2p5->size(); ++i){
        if(GENjetsPt_pt30_eta2p5->at(i) > _GENpTj2 && _GENpTj1 != GENjetsPt_pt30_eta2p5->at(i)) {
          _GENpTj2 = GENjetsPt_pt30_eta2p5->at(i);
          Mj2 = GENjetsMass_pt30_eta2p5->at(i);
          ETAj2 = GENjetsEta_pt30_eta2p5->at(i);
          PHIj2 = GENjetsPhi_pt30_eta2p5->at(i);
        }
      }

      // H+Njets GENvariables
      // If the event does not pass fiducial selections the variable is equal to one
      _GENpTHj = -1;
      _GENpTHjj = -1;
      _GENmHj = -1;
      _GENmHjj = -1;
      _GENdetajj = -1;
      _GENabsdetajj = -1;
      _GENmjj = -1;
      _GENdphijj = -1;
      _GENabsdphijj = -1;
      _GENTCjmax = 0; _GENTBjmax = 0;
      _GENTCj = -1; _GENTBj = -1;
      _GENTCj1 = 0; _GENTBj1 = 0;
      if(passedFiducialSelection_bbf){

        TLorentzVector GENH;
        TLorentzVector GENj1;
        TLorentzVector GENj2;
        GENH.SetPtEtaPhiM(GENpT4l,GENeta4l,GENphi4l,GENmass4l);
        GENj1.SetPtEtaPhiM(_GENpTj1,ETAj1,PHIj1,Mj1);
        GENj2.SetPtEtaPhiM(_GENpTj2,ETAj2,PHIj2,Mj2);
        _GENTCj1 = sqrt(_GENpTj1*_GENpTj1 + Mj1*Mj1)/(2*cosh(GENj1.Rapidity() - GENH.Rapidity()));
        _GENTBj1 = sqrt(_GENpTj1*_GENpTj1 + Mj1*Mj1)*exp(-1*abs(GENj1.Rapidity() - GENH.Rapidity()));

        for (unsigned int i = 0; i < GENjetsPt_pt30_eta2p5->size(); ++i){
          TLorentzVector theJet;
          theJet.SetPtEtaPhiM(GENjetsPt_pt30_eta2p5->at(i), GENjetsEta_pt30_eta2p5->at(i), GENjetsPhi_pt30_eta2p5->at(i), GENjetsMass_pt30_eta2p5->at(i));
          _GENTCj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*cosh(theJet.Rapidity() - GENH.Rapidity()));
          _GENTBj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*exp(-1*(theJet.Rapidity() - GENH.Rapidity()));
          if (_GENTCj > _GENTCjmax) _GENTCjmax = _GENTCj;
          if (_GENTBj > _GENTBjmax) _GENTBjmax = _GENTBj;
        }

        _GENpTHj = (GENH+GENj1).Pt();
        _GENpTHjj = (GENH+GENj1+GENj2).Pt();
        _GENmHj = (GENH+GENj1).M();
        _GENmHjj = (GENH+GENj1+GENj2).M();
        _GENdetajj = GENj1.Eta()-GENj2.Eta();
        _GENabsdetajj = abs(_GENdetajj);
        _GENmjj = (GENj1+GENj2).M();
        _GENdphijj = PHIj1 - PHIj2;//GENj1.Phi()-GENj2.Phi();
        //if(_GENdphijj>3.14) _GENdphijj -= 3.14;
        //if(_GENdphijj<-3.14) _GENdphijj += 3.14;
        _GENabsdphijj = abs(_GENdphijj);
      }

      // MELA probabilities
      if(passedFiducialSelection_bbf){
        _GEN_D0m = p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + (p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen * pow(spline_g4->Eval(GENmass4l),2)));
        _GEN_Dcp = p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen / (2 * sqrt(p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen * p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen));

        _GEN_D0hp = p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + (p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen * pow(spline_g2->Eval(GENmass4l),2)));
        _GEN_Dint = p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * sqrt(p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen * p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen));

        _GEN_DL1 = p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(spline_L1->Eval(GENmass4l),2)));
        _GEN_DL1int = (p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen/1e4) / (2 * sqrt(p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen * (p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8)));

        _GEN_DL1Zg = p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(spline_L1Zgs->Eval(GENmass4l),2)));
        _GEN_DL1Zgint = (p_GEN_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen/1e4) / (2 * sqrt(p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen * (p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8)));
      }else{
        _GEN_Dcp = -2;
        _GEN_D0m = -2;
        _GEN_Dint = -2;
        _GEN_D0hp = -2;
        _GEN_DL1 = -2;
        _GEN_DL1int = -2;
        _GEN_DL1Zg = -2;
        _GEN_DL1Zgint = -2;
      }

      GENpTj1->Fill();
      // GENpTj1_eta4p7->Fill();
      GENpTj2->Fill();
      GENpTHj->Fill();
      GENpTHjj->Fill();
      GENmHj->Fill();
      GENmHjj->Fill();
      GENdetajj->Fill();
      GENdphijj->Fill();
      GENabsdetajj->Fill();
      GENabsdphijj->Fill();
      GENmjj->Fill();
      GENTCjmax->Fill();
      GENTBjmax->Fill();
      // GenLepPtSorted->Fill();
      // GenLepEtaSorted->Fill();
      // GenLepPhiSorted->Fill();
      // GenLepIdSorted->Fill();
      passedFullSelection->Fill();
      GENnjets_pt30_eta2p5->Fill();
      GENnjets_pt30_eta4p7->Fill();

      GENTCj1->Fill();
      GENTBj1->Fill();

      GEN_Dcp->Fill();
      GEN_D0m->Fill();
      GEN_Dint->Fill();
      GEN_D0hp->Fill();
      GEN_DL1->Fill();
      GEN_DL1int->Fill();
      GEN_DL1Zg->Fill();
      GEN_DL1Zgint->Fill();

      // GenLepSorted.clear();
      // _GenLepPtSorted.clear();
      // _GenLepEtaSorted.clear();
      // _GenLepPhiSorted.clear();
      // _GenLepIdSorted.clear();
    }
    if (t_failed) continue; // From now on reco-only variables

    // leading jet pT
    _pTj1 = 0; Mj1 = 0; ETAj1 = 0; PHIj1 = 0;
    _pTj1_eta4p7 = 0;

    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5 && JetPt->at(i) > _pTj1) {
        _pTj1 = JetPt->at(i);
        Mj1 = JetMass->at(i);
        ETAj1 = JetEta->at(i);
        PHIj1 = JetPhi->at(i);
      }

      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7 && JetPt->at(i) > _pTj1_eta4p7) {
        _pTj1_eta4p7 = JetPt->at(i);
      }

      // if(JetPt_JESUp_Total->at(i)>30 && abs(JetEta->at(i))<2.5 && JetPt_JESUp_Total->at(i) > _pTj1_jesup) {
      //   _pTj1_jesup = JetPt_JESUp_Total->at(i);
      //   Mj1_jesup = JetMass->at(i);
      //   ETAj1_jesup = JetEta->at(i);
      //   PHIj1_jesup = JetPhi->at(i);
      // }
    }
    // sub-leading jet pT
    _pTj2 = 0; Mj2 = 0; ETAj2 = 0; PHIj2 = 0;
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5 && JetPt->at(i) > _pTj2 && _pTj1 != JetPt->at(i)) {
        _pTj2 = JetPt->at(i);
        Mj2 = JetMass->at(i);
        ETAj2 = JetEta->at(i);
        PHIj2 = JetPhi->at(i);
      }
    }

    // H+Njets variables
    TLorentzVector H;
    TLorentzVector j1;
    TLorentzVector j2;
    // TLorentzVector j1_jesup;
    // TLorentzVector j1_jesdn;
    // TLorentzVector j2_jesup;
    // TLorentzVector j2_jesdn;
    H.SetPtEtaPhiM(ZZPt,ZZEta,ZZPhi,ZZMass);
    j1.SetPtEtaPhiM(_pTj1,ETAj1,PHIj1,Mj1);
    // j1_jesup.SetPtEtaPhiM(_pTj1_jesup,ETAj1_jesup,PHIj1_jesup,Mj1_jesup);
    // j1_jesdn.SetPtEtaPhiM(_pTj1_jesdn,ETAj1_jesdn,PHIj1_jesdn,Mj1_jesdn);
    j2.SetPtEtaPhiM(_pTj2,ETAj2,PHIj2,Mj2);
    // j2_jesup.SetPtEtaPhiM(_pTj2_jesup,ETAj2_jesup,PHIj2_jesup,Mj2_jesup);
    // j2_jesdn.SetPtEtaPhiM(_pTj2_jesdn,ETAj2_jesdn,PHIj2_jesdn,Mj2_jesdn);
    _pTHj = (H+j1).Pt();
    // _pTHj_jesup = (H+j1_jesup).Pt();
    // _pTHj_jesdn = (H+j1_jesdn).Pt();
    _pTHjj = (H+j1+j2).Pt();
    // _pTHjj_jesup = (H+j1_jesdn+j2_jesdn).Pt();
    // _pTHjj_jesdn = (H+j1_jesup+j2_jesup).Pt();
    _mHj = (H+j1).M();
    // _mHj_jesup = (H+j1_jesup).M();
    // _mHj_jesdn = (H+j1_jesdn).M();
    _mHjj = (H+j1+j2).M();
    // _mHjj_jesup = (H+j1_jesup+j2_jesup).M();
    // _mHjj_jesdn = (H+j1_jesdn+j2_jesdn).M();
    _detajj = j1.Eta()-j2.Eta();
    // _detajj_jesdn = j1_jesdn.Eta()-j2_jesdn.Eta();
    // _detajj_jesup = j1_jesup.Eta()-j2_jesup.Eta();
    _absdetajj = abs(_detajj);
    // _absdetajj_jesdn = abs(_detajj_jesdn);
    // _absdetajj_jesup = abs(_detajj_jesup);
    _mjj = (j1+j2).M();
    // _mjj_jesup = (j1_jesup+j2_jesup).M();
    // _mjj_jesdn = (j1_jesdn+j2_jesdn).M();
    _dphijj = j1.Phi() - j2.Phi();
    // if(_dphijj_jesdn>3.14) _dphijj_jesdn -= 3.14;
    // if(_dphijj_jesup>3.14) _dphijj_jesup -= 3.14;
    // if(_dphijj_jesup<-3.14) _dphijj_jesup += 3.14;
    // if(_dphijj_jesdn<-3.14) _dphijj_jesdn += 3.14;
    _absdphijj = abs(_dphijj);
    // _absdphijj_jesdn = abs(_dphijj_jesdn);
    // _absdphijj_jesup = abs(_dphijj_jesup);

    // njets
    _njets_pt30_eta2p5 = 0;
    for(unsigned int i=0;i<JetPt->size();i++){
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5){
        _njets_pt30_eta2p5++;
      }
    }
    // njets 4.7 eta
    _njets_pt30_eta4p7 = 0;
    for(unsigned int i=0;i<JetPt->size();i++){
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7){
        _njets_pt30_eta4p7++;
      }
    }
    // JES down
    // _njets_pt30_eta2p5_jesdn = 0;
    // for(unsigned int i=0;i<JetPt_JESDown->size();i++){
    //   if(JetPt_JESDown->at(i)>30 && abs(JetEta->at(i))<2.5){ // eta is not affected by jes
    //     _njets_pt30_eta2p5_jesdn++;
    //   }
    // }
    // // JES up
    // _njets_pt30_eta2p5_jesup = 0;
    // for(unsigned int i=0;i<JetPt_JESUp->size();i++){
    //   if(JetPt_JESUp->at(i)>30 && abs(JetEta->at(i))<2.5){ // eta is not affected by jes
    //     _njets_pt30_eta2p5_jesup++;
    //   }
    // }

    pTj1->Fill();
    pTj1_eta4p7->Fill();
    // pTj1_jesup->Fill();
    // pTj1_jesdn->Fill();
    pTj2->Fill();
    // pTj2_jesup->Fill();
    // pTj2_jesdn->Fill();
    pTHj->Fill();
    // pTHj_jesup->Fill();
    // pTHj_jesdn->Fill();
    pTHjj->Fill();
    // pTHjj_jesup->Fill();
    // pTHjj_jesdn->Fill();
    mHj->Fill();
    // mHj_jesup->Fill();
    // mHj_jesdn->Fill();
    mHjj->Fill();
    // mHjj_jesup->Fill();
    // mHjj_jesdn->Fill();
    detajj->Fill();
    // detajj_jesup->Fill();
    // detajj_jesdn->Fill();
    dphijj->Fill();
    // dphijj_jesup->Fill();
    // dphijj_jesdn->Fill();
    absdetajj->Fill();
    absdphijj->Fill();
    mjj->Fill();
    // mjj_jesup->Fill();
    // mjj_jesdn->Fill();
    njets_pt30_eta2p5->Fill();
    njets_pt30_eta4p7->Fill();
    // njets_pt30_eta2p5_jesdn->Fill();
    // njets_pt30_eta2p5_jesup->Fill();

    // Reco-rapidity
    _ZZy = abs(log((sqrt(125*125 + ZZPt*ZZPt*cosh(ZZEta)*cosh(ZZEta))+ZZPt*sinh(ZZEta))/sqrt(125*125+ZZPt*ZZPt)));
    ZZy->Fill();

    _TCj = 0; _TBj = 0;
    _TCjmax = 0; _TBjmax = 0;
    _TCj1 = 0; _TBj1 = 0;
    _TCj1 = sqrt(_pTj1*_pTj1 + Mj1*Mj1)/(2*cosh(j1.Rapidity() - H.Rapidity()));
    _TBj1 = sqrt(_pTj1*_pTj1 + Mj1*Mj1)*exp(-1*abs(j1.Rapidity() - H.Rapidity()));
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      // Define TCj for all jets with pT > pT_cut
       if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5) {
          TLorentzVector theJet;
          theJet.SetPtEtaPhiM(JetPt->at(i), JetEta->at(i), JetPhi->at(i), JetMass->at(i));
          _TCj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*cosh(theJet.Rapidity() - H.Rapidity())); //theJet.E());
          _TBj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*exp(-1*(theJet.Rapidity() - H.Rapidity()));
          if (_TCj > _TCjmax) _TCjmax = _TCj;
          if (_TBj > _TBjmax) _TBjmax = _TBj;
       }
    }
    TCj->Fill();
    TBj->Fill();
    TCjmax->Fill();
    TBjmax->Fill();
    TCj1->Fill();
    TBj1->Fill();

    // MELA probabilities
    _D0m = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + (p_GG_SIG_ghg2_1_ghz4_1_JHUGen * pow(spline_g4->Eval(ZZMass),2)));
    _Dcp = p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * p_GG_SIG_ghg2_1_ghz4_1_JHUGen));

    _D0hp = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + (p_GG_SIG_ghg2_1_ghz2_1_JHUGen * pow(spline_g2->Eval(ZZMass),2)));
    _Dint = p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * p_GG_SIG_ghg2_1_ghz2_1_JHUGen));

    _DL1 = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(spline_L1->Eval(ZZMass),2)));
    _DL1int = (p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen/1e4) / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * (p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8)));

    _DL1Zg = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(spline_L1Zgs->Eval(ZZMass),2)));
    _DL1Zgint = (p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen/1e4) / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * (p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8)));

    Dcp->Fill();
    D0m->Fill();
    D0hp->Fill();
    Dint->Fill();
    DL1->Fill();
    DL1int->Fill();
    DL1Zg->Fill();
    DL1Zgint->Fill();

    if(process=="signal" || process=="AC"){
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
  }
  T->Write("", TObject::kOverwrite);
  delete f;
  return;
}

//---------------------------------------------------------- MAIN ----------------------------------------------------------
void skim_MC_tree_v2 (TString prod_mode = "VBFH125", TString year = "2018"){

  TString process;
  if(prod_mode=="ZZTo4lext") process = "qqZZ";
  else if(prod_mode=="ggTo2e2mu_Contin_MCFM701" || prod_mode=="ggTo2e2tau_Contin_MCFM701" || prod_mode=="ggTo2mu2tau_Contin_MCFM701" || prod_mode=="ggTo4e_Contin_MCFM701" ||
          prod_mode=="ggTo4mu_Contin_MCFM701" || prod_mode=="ggTo4tau_Contin_MCFM701") process = "ggZZ";
  else if(prod_mode.Contains("H125")) process = "signal"; //If "H125" is in the name of the prod_mode, it is a signal process
  else process = "AC";
  if(prod_mode=="ZZTo4lext" && year=="2018") prod_mode = "ZZTo4lext1"; //Change prod_mode label for qqZZ 2018

  cout << process << endl;

  TString input_dir, full_path;

  if(process=="AC"){
    input_dir = "/eos/user/a/atarabin/MC_samples";
    full_path = Form("%s/AC%s/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());
    cout << full_path << endl;
  }
  else{
    input_dir = "/eos/user/a/atarabin/MC_samples";
    full_path = Form("%s/%s_MELA/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());
    cout << full_path << endl;
  }
  std::cout << input_dir << " " << year << " " << prod_mode << std::endl;

  auto oldFile = TFile::Open(full_path.Data());
  TTree *oldtree = (TTree*) oldFile->Get("ZZTree/candTree");
  TH1F *hCounters = (TH1F*) oldFile->Get("ZZTree/Counters");
  // If it is a bkg process we don't deal with gen-level information. In case of bkgs the oldtree_failed remains empty
  TTree *oldtree_failed;
  if(process=="signal" || process=="AC") oldtree_failed = (TTree*) oldFile->Get("ZZTree/candTree_failed");

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
  oldtree->SetBranchStatus("ZZPhi",1);
  oldtree->SetBranchStatus("costhetastar",1);
  oldtree->SetBranchStatus("helcosthetaZ1",1);
  oldtree->SetBranchStatus("helcosthetaZ2",1);
  oldtree->SetBranchStatus("helphi",1);
  oldtree->SetBranchStatus("phistarZ1",1);
  // oldtree->SetBranchStatus("nCleanedJetsPt30",1);
  // oldtree->SetBranchStatus("nCleanedJetsPt30_jesUp",1);
  // oldtree->SetBranchStatus("nCleanedJetsPt30_jesDn",1);
  oldtree->SetBranchStatus("JetPt",1);
  oldtree->SetBranchStatus("JetEta",1);
  oldtree->SetBranchStatus("JetPhi",1);
  oldtree->SetBranchStatus("JetMass",1);
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
  // JESUp
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree->SetBranchStatus(Form("JetPt_JESUp_%s", jes_name.at(i).Data()),1);
    // oldtree->SetBranchStatus(Form("nCleanedJetsPt30_jesUp_%s", jes_name.at(i).Data()),1);
  }
  // JESDown
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree->SetBranchStatus(Form("JetPt_JESDown_%s", jes_name.at(i).Data()),1);
    // oldtree->SetBranchStatus(Form("nCleanedJetsPt30_jesDn_%s", jes_name.at(i).Data()),1);
  }
  if(process=="qqZZ") {
    oldtree->SetBranchStatus("KFactor_EW_qqZZ",1);
    oldtree->SetBranchStatus("KFactor_QCD_qqZZ_M",1);
  }
  if(process=="ggZZ") {
    oldtree->SetBranchStatus("KFactor_QCD_ggZZ_Nominal",1);
  }
  oldtree->SetBranchStatus("PUWeight",1);
  oldtree->SetBranchStatus("genHEPMCweight",1);
  if(year!="2016") oldtree->SetBranchStatus("genHEPMCweight_NNLO",1);
  oldtree->SetBranchStatus("overallEventWeight",1);
  oldtree->SetBranchStatus("L1prefiringWeight",1);
  oldtree->SetBranchStatus("dataMCWeight",1);
  oldtree->SetBranchStatus("trigEffWeight",1);
  //MELA probabilities
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz4_1_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz2_1_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",1);
  oldtree->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",1);
  if(process=="signal"){
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
    oldtree->SetBranchStatus("GENlep_RelIso",1);
    oldtree->SetBranchStatus("GENH_pt",1);
    oldtree->SetBranchStatus("GENH_eta",1);
    oldtree->SetBranchStatus("GENH_phi",1);
    oldtree->SetBranchStatus("GENH_mass",1);
    oldtree->SetBranchStatus("GENmass4l",1);
    oldtree->SetBranchStatus("GENpT4l",1);
    oldtree->SetBranchStatus("GENeta4l",1);
    oldtree->SetBranchStatus("GENphi4l",1);
    oldtree->SetBranchStatus("GENrapidity4l",1);
    oldtree->SetBranchStatus("GENcosTheta1",1);
    oldtree->SetBranchStatus("GENcosTheta2",1);
    oldtree->SetBranchStatus("GENcosThetaStar",1);
    oldtree->SetBranchStatus("GENPhi",1);
    oldtree->SetBranchStatus("GENPhi1",1);
    oldtree->SetBranchStatus("GENmassZ1",1);
    oldtree->SetBranchStatus("GENmassZ2",1);
    oldtree->SetBranchStatus("GENZ_pt",1);
    oldtree->SetBranchStatus("GENZ_eta",1);
    oldtree->SetBranchStatus("GENZ_phi",1);
    oldtree->SetBranchStatus("GENZ_mass",1);
    oldtree->SetBranchStatus("GENZ_DaughtersId",1);
    oldtree->SetBranchStatus("GENZ_MomId",1);
    //GENjet
    oldtree->SetBranchStatus("GENjetsPt_pt30_eta4p7",1);
    oldtree->SetBranchStatus("GENjetsEta_pt30_eta4p7",1);
    oldtree->SetBranchStatus("GENjetsPhi_pt30_eta4p7",1);
    oldtree->SetBranchStatus("GENjetsMass_pt30_eta4p7",1);
    oldtree->SetBranchStatus("GENjetsPt_pt30_eta2p5",1);
    oldtree->SetBranchStatus("GENjetsEta_pt30_eta2p5",1);
    oldtree->SetBranchStatus("GENjetsPhi_pt30_eta2p5",1);
    oldtree->SetBranchStatus("GENjetsMass_pt30_eta2p5",1);
    //MELA probabilities
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz2_1_ghz4_1_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",1);
    oldtree->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",1);
  }
  if(prod_mode == "ggH125") oldtree->SetBranchStatus("ggH_NNLOPS_weight",1); // Additional entry for the weight in case of ggH

  //skim oldtree_failed for signal only
  if(process=="signal" || process=="AC"){
    //// candTree_failed
    // Deactivate all branches
    oldtree_failed->SetBranchStatus("*",0);
    // Activate some branches only: our skim
    oldtree_failed->SetBranchStatus("EventNumber",1);
    oldtree_failed->SetBranchStatus("xsec",1);
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
    if(year!="2016") oldtree_failed->SetBranchStatus("genHEPMCweight_NNLO",1);
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
    oldtree_failed->SetBranchStatus("GENlep_RelIso",1);
    oldtree_failed->SetBranchStatus("GENH_pt",1);
    oldtree_failed->SetBranchStatus("GENH_eta",1);
    oldtree_failed->SetBranchStatus("GENH_phi",1);
    oldtree_failed->SetBranchStatus("GENH_mass",1);
    oldtree_failed->SetBranchStatus("GENmass4l",1);
    oldtree_failed->SetBranchStatus("GENpT4l",1);
    oldtree_failed->SetBranchStatus("GENeta4l",1);
    oldtree_failed->SetBranchStatus("GENphi4l",1);
    oldtree_failed->SetBranchStatus("GENrapidity4l",1);
    oldtree_failed->SetBranchStatus("GENcosTheta1",1);
    oldtree_failed->SetBranchStatus("GENcosTheta2",1);
    oldtree_failed->SetBranchStatus("GENcosThetaStar",1);
    oldtree_failed->SetBranchStatus("GENPhi",1);
    oldtree_failed->SetBranchStatus("GENPhi1",1);
    oldtree_failed->SetBranchStatus("GENZ_pt",1);
    oldtree_failed->SetBranchStatus("GENZ_eta",1);
    oldtree_failed->SetBranchStatus("GENZ_phi",1);
    oldtree_failed->SetBranchStatus("GENZ_mass",1);
    oldtree_failed->SetBranchStatus("GENZ_DaughtersId",1);
    oldtree_failed->SetBranchStatus("GENZ_MomId",1);
    oldtree_failed->SetBranchStatus("GENmassZ1",1);
    oldtree_failed->SetBranchStatus("GENmassZ2",1);
    //GENjet
    oldtree_failed->SetBranchStatus("GENjetsPt_pt30_eta4p7",1);
    oldtree_failed->SetBranchStatus("GENjetsEta_pt30_eta4p7",1);
    oldtree_failed->SetBranchStatus("GENjetsPhi_pt30_eta4p7",1);
    oldtree_failed->SetBranchStatus("GENjetsMass_pt30_eta4p7",1);
    oldtree_failed->SetBranchStatus("GENjetsPt_pt30_eta2p5",1);
    oldtree_failed->SetBranchStatus("GENjetsEta_pt30_eta2p5",1);
    oldtree_failed->SetBranchStatus("GENjetsPhi_pt30_eta2p5",1);
    oldtree_failed->SetBranchStatus("GENjetsMass_pt30_eta2p5",1);
    //MELA probabilities
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz4_1_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz2_1_ghz4_1_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz2_1_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",1);
    oldtree_failed->SetBranchStatus("p_GEN_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",1);
    if(prod_mode == "ggH125") oldtree_failed->SetBranchStatus("ggH_NNLOPS_weight",1); // Additional entry for the weight in case of ggH
  }

  // Copy branches in the new file
  if(process!="signal" && process!="AC") input_dir = "/eos/user/a/atarabin/MC_samples"; //Bkg only (At the moment bkgs original root file are not stored in our folder but in CJLST's)
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path;
  if(process!="AC") new_full_path = Form("%s/%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  else new_full_path = Form("%s/AC%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  cout << new_full_path << endl;
  TFile *newfile = new TFile(new_full_path.Data(),"RECREATE");
  TTree *newtree = (TTree*) oldtree->CloneTree(0);
  newtree->CopyEntries(oldtree);
  newtree->Write(); // Write candTree
  if(process=="signal" || process=="AC"){
    TTree *newtree_failed = (TTree*) oldtree_failed->CloneTree(0);
    newtree_failed->CopyEntries(oldtree_failed);
    newtree_failed->Write(); // Write candTree_failed
  }
  hCounters->Write(); // Write Counters
  newfile->Close();

  bool t_failed;
  bool year_tmp = false;
  if(year=="2017" && process=="signal") year_tmp = true;
  if(process=="signal" || process=="AC") add(input_dir, year, prod_mode, process, t_failed = true, year_tmp);
  add(input_dir, year, prod_mode, process, t_failed = false, year_tmp);

  if(process=="signal" || process=="AC"){
    // Merge together into a single TTree. Useful for efficiencies calculation.
    TFile* inputfile = TFile::Open(new_full_path.Data(), "READ");
    TTree* tree1 = (TTree*) inputfile->Get("candTree");
    TTree* tree2 = (TTree*) inputfile->Get("candTree_failed");
    TH1F* cnts = (TH1F*) inputfile->Get("Counters");

    TString merged_name = Form("%s_mergedTree_MC_%s.root", prod_mode.Data(), year.Data());
    TString merged_path;
    if(process!="AC") merged_path = Form("%s/%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),merged_name.Data());
    else merged_path = Form("%s/AC%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),merged_name.Data());
    cout << merged_path << endl;

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
  }

  return 0;
}
