// root -q -b skim_nano.C"(\"ggH125\",\"2018\")"

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
// #include "DataFormats/Math/interface/deltaR.h"

using namespace std;

// float mass_lep(int flavour){
//   if((abs(flavour)) == 11) return 0.0005109989461;
//   else if ((abs(flavour)) == 13) return 0.1056583745;
//   else if ((abs(flavour)) == 15) return 1.77686;
//   return 0;
// }

// //------------------------------------------------------------------
void add(TString input_dir, TString year, TString prod_mode, TString process, bool t_failed=true, bool flag_tmp_2017=false){
//   // Add additional branches
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path;
  if(process!="AC") new_full_path = Form("%s", new_name.Data());
  else new_full_path = Form("%s/AC%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  cout << new_full_path << endl;
  TFile *f = new TFile(new_full_path.Data(),"UPDATE");
  TTree *T;
  if (t_failed) {
    std::cout << "Filling TTree _failed" << std::endl;
    T = (TTree*)f->Get("AllEvents");
  } else {
    std::cout << "Filling TTree ZZ cands" << std::endl;
    T = (TTree*)f->Get("Events");
  }
  std::cout << T->GetName() << std::endl;
//   // Gen-variables
//   float GenLep1Pt,GenLep2Pt,GenLep3Pt,GenLep4Pt,GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta,GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi,GenZ1Flav,GenZ2Flav,
//         GenZ1Mass,GenZ2Mass,GenLep1Iso,GenLep2Iso,GenLep3Iso,GenLep4Iso;
//   Float_t GENmass4l,GENpT4l,GENeta4l,GENphi4l;
//   Short_t GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  // bool _passedFullSelection, passedFiducialSelection_bbf;
  bool _passedFullSelection;

  // Bool_t passedFiducial;

//   vector<float> *GENlep_pt = 0;
//   vector<float> *GENlep_eta = 0;
//   vector<float> *GENlep_phi = 0;
//   vector<float> *GENlep_mass = 0;
//   vector<float> *GENlep_id = 0;
//   vector<Short_t> *GENlep_Hindex = 0;
  TBranch *passedFullSelection = T->Branch("passedFullSelection",&_passedFullSelection,"passedFullSelection/B");

  // if(process=="signal" || process=="AC"){ // Bkgs don't store gen-level information
//     T->SetBranchAddress("GENlep_pt",&GENlep_pt);
//     T->SetBranchAddress("GENlep_eta",&GENlep_eta);
//     T->SetBranchAddress("GENlep_phi",&GENlep_phi);
//     T->SetBranchAddress("GENlep_mass",&GENlep_mass);
//     T->SetBranchAddress("GENlep_id",&GENlep_id);
//     T->SetBranchAddress("GENmass4l",&GENmass4l);
//     T->SetBranchAddress("GENpT4l",&GENpT4l);
//     T->SetBranchAddress("GENeta4l",&GENeta4l);
//     T->SetBranchAddress("GENphi4l",&GENphi4l);
//     T->SetBranchAddress("GENlep_Hindex",&GENlep_Hindex);
    // T->SetBranchAddress("passedFiducial",&passedFiducial);
  // }

//   // Reco-variables and Gen-Reco-matching variables
//   float _ZZy,ZZPt,ZZEta,ZZPhi,ZZMass;
//   vector<float> *LepPt = 0;
//   vector<float> *LepPhi = 0;
//   vector<float> *LepEta = 0;
//   vector<int> *LepLepId = 0;
//   vector<float> *ExtraLepPt = 0;
//   vector<float> *ExtraLepEta = 0;
//   vector<float> *ExtraLepPhi = 0;
//   vector<int> *ExtraLepLepId = 0;
//   vector<int> _lep_genindex, _lep_Hindex;

//   TBranch *ZZy = T->Branch("ZZy",&_ZZy,"ZZy/F");
//   TBranch *lep_genindex = T->Branch("lep_genindex",&_lep_genindex);
//   TBranch *lep_Hindex = T->Branch("lep_Hindex",&_lep_Hindex);

//   if (!t_failed) {
//     T->SetBranchAddress("ZZMass",&ZZMass);
//     T->SetBranchAddress("ZZPt",&ZZPt);
//     T->SetBranchAddress("ZZEta",&ZZEta);
//     T->SetBranchAddress("ZZPhi",&ZZPhi);
//     T->SetBranchAddress("LepPt",&LepPt);
//     T->SetBranchAddress("LepPhi",&LepPhi);
//     T->SetBranchAddress("LepEta",&LepEta);
//     T->SetBranchAddress("LepLepId",&LepLepId);
//     T->SetBranchAddress("ExtraLepPt",&ExtraLepPt);
//     T->SetBranchAddress("ExtraLepEta",&ExtraLepEta);
//     T->SetBranchAddress("ExtraLepPhi",&ExtraLepPhi);
//     T->SetBranchAddress("ExtraLepLepId",&ExtraLepLepId);
//   }

  _event = 0;

  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);

    if(process=="signal" || process=="AC"){

      _passedFullSelection = true;
      if (t_failed) {
  	     _passedFullSelection = false;
      }

      passedFullSelection->Fill();
    }

    if (t_failed) continue; // From now on reco-only variables

//     // Reco-rapidity
//     _ZZy = abs(log((sqrt(125*125 + ZZPt*ZZPt*cosh(ZZEta)*cosh(ZZEta))+ZZPt*sinh(ZZEta))/sqrt(125*125+ZZPt*ZZPt)));
//     ZZy->Fill();

//     if(process=="signal" || process=="AC"){
//       // GEN matching
//       for(unsigned int i=0;i<LepPt->size();i++){
//         _lep_genindex.push_back(-1);
//       }
//       for(unsigned int i = 0; i < LepPt->size(); i++) {
//           double minDr=9999.0;
//           TLorentzVector reco, gen;
//           reco.SetPtEtaPhiM(LepPt->at(i),LepEta->at(i),LepPhi->at(i),mass_lep(LepLepId->at(i)));
//           for (unsigned int j = 0; j < GENlep_id->size(); j++) {
//               if (GENlep_id->at(j)!=LepLepId->at(i)) continue;
//               gen.SetPtEtaPhiM(GENlep_pt->at(j),GENlep_eta->at(j),GENlep_phi->at(j),GENlep_mass->at(j));
//               double thisDr = deltaR(reco.Eta(),reco.Phi(),gen.Eta(),gen.Phi());
//               // double thisDr = reco.DeltaR(gen);
//               if (thisDr<minDr && thisDr<0.5) {
//                   _lep_genindex[i]=j;
//                   minDr=thisDr;
//               }
//           } // all gen leptons
//       } // all reco leptons
//       lep_genindex->Fill();
//       _lep_genindex.clear();


//       for(unsigned int i=0;i<LepPt->size();i++){
//         _lep_Hindex.push_back(-1);
//       }
//       float lead_Z1 = max(LepPt->at(0),LepPt->at(1));
//       if(lead_Z1 == LepPt->at(0)){
//         _lep_Hindex[0] = 0;
//         _lep_Hindex[1] = 1;
//       }
//       else if(lead_Z1 == LepPt->at(1)){
//         _lep_Hindex[0] = 1;
//         _lep_Hindex[1] = 0;
//       }
//       float lead_Z2 = max(LepPt->at(2),LepPt->at(3));
//       if(lead_Z2 == LepPt->at(2)){
//         _lep_Hindex[2] = 2;
//         _lep_Hindex[3] = 3;
//       }
//       else if(lead_Z2 == LepPt->at(3)){
//         _lep_Hindex[2] = 3;
//         _lep_Hindex[3] = 2;
//       }
//       lep_Hindex->Fill();
//       _lep_Hindex.clear();
//     }
  } // End loop on all events
  T->Write("", TObject::kOverwrite);
  delete f;
  return;
}

//---------------------------------------------------------- MAIN ----------------------------------------------------------
void skim_nano (TString prod_mode = "VBFH125", TString year = "2018"){

  TString process;
  if(prod_mode=="ZZTo4lext") process = "qqZZ";
  else if(prod_mode=="ggTo2e2mu_Contin_MCFM701" || prod_mode=="ggTo2e2tau_Contin_MCFM701" || prod_mode=="ggTo2mu2tau_Contin_MCFM701" || prod_mode=="ggTo4e_Contin_MCFM701" ||
          prod_mode=="ggTo4mu_Contin_MCFM701" || prod_mode=="ggTo4tau_Contin_MCFM701") process = "ggZZ";
  else if(prod_mode.Contains("H12")) process = "signal"; //If "H125" is in the name of the prod_mode, it is a signal process
  else process = "AC";
  if(prod_mode=="ZZTo4lext" && year=="2018") prod_mode = "ZZTo4lext1"; //Change prod_mode label for qqZZ 2018

  cout << process << endl;

  TString input_dir, full_path;

  input_dir = ".";
  // full_path = Form("%s/%s_MELA/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());
  full_path = Form("ZZ4lAnalysis.root");
  cout << full_path << endl;

  std::cout << input_dir << " " << year << " " << prod_mode << std::endl;

  auto oldFile = TFile::Open(full_path.Data());
  TTree *oldtree = (TTree*) oldFile->Get("Events");

  // If it is a bkg process we don't deal with gen-level information. In case of bkgs the oldtree_failed remains empty
  TTree *oldtree_failed;
  if(process=="signal") oldtree_failed = (TTree*) oldFile->Get("AllEvents");

  // candTree
  // Deactivate all branches
  oldtree->SetBranchStatus("*",0);
  // Activate some branches only: our skim
  oldtree->SetBranchStatus("event",1);
  // TODO: Do we have xsec in the nanoAOD processing?
  // TODO: Either save it or add it in the new uproot skimmer
  // TODO: Doesn't make sense not to save it IMHO
  // oldtree->SetBranchStatus("xsec",1);
  oldtree->SetBranchStatus("ZZCand_mass",1);
  oldtree->SetBranchStatus("ZZCand_Z1mass",1);
  oldtree->SetBranchStatus("ZZCand_Z2mass",1);
  oldtree->SetBranchStatus("ZZCand_Z1flav",1);
  oldtree->SetBranchStatus("ZZCand_Z2flav",1);
  // oldtree->SetBranchStatus("ZZCand_pt",1);
  // oldtree->SetBranchStatus("ZZCand_eta",1);
  // oldtree->SetBranchStatus("ZZCand_phi",1);
  // TODO: Leptons are at ZZCand_Z1l1Idx, how do I access them?
  // TODO: Was LepPt from the sorted leptons? If so, only in pT?
  // oldtree->SetBranchStatus("LepPt",1);
  // oldtree->SetBranchStatus("LepEta",1);
  // oldtree->SetBranchStatus("LepPhi",1);
  // oldtree->SetBranchStatus("LepLepId",1);
  // oldtree->SetBranchStatus("ExtraLepPt",1);
  // oldtree->SetBranchStatus("ExtraLepPhi",1);
  // oldtree->SetBranchStatus("ExtraLepEta",1);
  // oldtree->SetBranchStatus("ExtraLepLepId",1);

  if(process=="qqZZ") {
    oldtree->SetBranchStatus("KFactor_EW_qqZZ_Weight",1);
    oldtree->SetBranchStatus("KFactor_QCD_qqZZ_M_Weight",1);
  }
  if(process=="ggZZ") {
    oldtree->SetBranchStatus("KFactor_QCD_ggZZ_Nominal_Weight",1);
  }

  oldtree->SetBranchStatus("puWeight",1);
  oldtree->SetBranchStatus("genWeight",1); // TODO: Check if this is genHEPMCweight
  // if(year!="2016") oldtree->SetBranchStatus("genHEPMCweight_NNLO",1);
  oldtree->SetBranchStatus("overallEventWeight",1);
  // TODO: Not available in Nano02Apr2020
  // oldtree->SetBranchStatus("L1prefiringWeight",1);
  oldtree->SetBranchStatus("ZZCand_dataMCWeight",1);
  // oldtree->SetBranchStatus("trigEffWeight",1);

  if(process=="signal"){
    oldtree->SetBranchStatus("passedFiducial",1);
    oldtree->SetBranchStatus("FidDressedLeps_pt",1);
    oldtree->SetBranchStatus("FidDressedLeps_eta",1);
    oldtree->SetBranchStatus("FidDressedLeps_phi",1);
    oldtree->SetBranchStatus("FidDressedLeps_mass",1);
    oldtree->SetBranchStatus("FidDressedLeps_id",1);
    oldtree->SetBranchStatus("FidDressedLeps_momid",1);
    oldtree->SetBranchStatus("FidDressedLeps_mommomid",1);
    oldtree->SetBranchStatus("FidDressedLeps_RelIso",1);
    oldtree->SetBranchStatus("FidZZ_Z1l1Idx",1);
    oldtree->SetBranchStatus("FidZZ_Z1l2Idx",1);
    oldtree->SetBranchStatus("FidZZ_Z2l1Idx",1);
    oldtree->SetBranchStatus("FidZZ_Z2l2Idx",1);
    oldtree->SetBranchStatus("FidZZ_mass",1);
    oldtree->SetBranchStatus("FidZZ_pt",1);
    oldtree->SetBranchStatus("FidZZ_eta",1);
    oldtree->SetBranchStatus("FidZZ_phi",1);
    oldtree->SetBranchStatus("FidZZ_rapidity",1);
    oldtree->SetBranchStatus("FidZ_DauPdgId",1);
    oldtree->SetBranchStatus("FidZ_MomPdgId",1);
  }
  if(prod_mode == "ggH125") oldtree->SetBranchStatus("ggH_NNLOPS_Weight",1); // Additional entry for the weight in case of ggH

  // //skim oldtree_failed for signal only
  if(process=="signal" || process=="AC"){
  //   //// candTree_failed
  //   // Deactivate all branches
    oldtree_failed->SetBranchStatus("*",0);
  //   // Activate some branches only: our skim
  //   oldtree_failed->SetBranchStatus("event",1);
  //   oldtree_failed->SetBranchStatus("xsec",1);
  //   oldtree_failed->SetBranchStatus("puWeight",1);
    oldtree_failed->SetBranchStatus("genWeight",1);
    oldtree_failed->SetBranchStatus("passedFiducial",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_pt",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_eta",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_phi",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_mass",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_id",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_momid",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_mommomid",1);
    oldtree_failed->SetBranchStatus("FidDressedLeps_RelIso",1);
    oldtree_failed->SetBranchStatus("FidZZ_Z1l1Idx",1);
    oldtree_failed->SetBranchStatus("FidZZ_Z1l2Idx",1);
    oldtree_failed->SetBranchStatus("FidZZ_Z2l1Idx",1);
    oldtree_failed->SetBranchStatus("FidZZ_Z2l2Idx",1);
    oldtree_failed->SetBranchStatus("FidZZ_mass",1);
    oldtree_failed->SetBranchStatus("FidZZ_pt",1);
    oldtree_failed->SetBranchStatus("FidZZ_eta",1);
    oldtree_failed->SetBranchStatus("FidZZ_phi",1);
    oldtree_failed->SetBranchStatus("FidZZ_rapidity",1);
    oldtree_failed->SetBranchStatus("FidZ_DauPdgId",1);
    oldtree_failed->SetBranchStatus("FidZ_MomPdgId",1);
    if(prod_mode == "ggH125") oldtree_failed->SetBranchStatus("ggH_NNLOPS_Weight",1); // Additional entry for the weight in case of ggH
  }

  // Copy branches in the new file
  if(process!="signal" && process!="AC") input_dir = ".";
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path;

  // if(process!="AC") new_full_path = Form("%s/%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  // else new_full_path = Form("%s/AC%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());

  new_full_path = Form("%s", new_name.Data());
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

  newfile->Close();

  bool t_failed;
  bool year_tmp = false;
  if(year=="2017" && process=="signal") year_tmp = true;
  if(process=="signal" || process=="AC") add(input_dir, year, prod_mode, process, t_failed = true, year_tmp);
  add(input_dir, year, prod_mode, process, t_failed = false, year_tmp);

  return 0;
}
