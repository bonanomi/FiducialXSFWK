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

//------------------------------------------------------------------
void add(TString input_dir, TString year, TString prod_mode, TString process, bool t_failed=true, bool flag_tmp_2017=false){
  // Add additional branches
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

  bool _passedFullSelection;

  TBranch *passedFullSelection = T->Branch("passedFullSelection",&_passedFullSelection,"passedFullSelection/B");

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

    if (t_failed) continue;
  }

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
  oldtree->SetBranchStatus("bestCandIdx",1);
  oldtree->SetBranchStatus("Counter",1);
  oldtree->SetBranchStatus("ZZCand_mass",1);
  oldtree->SetBranchStatus("ZZCand_Z1mass",1);
  oldtree->SetBranchStatus("ZZCand_Z2mass",1);
  oldtree->SetBranchStatus("ZZCand_Z1flav",1);
  oldtree->SetBranchStatus("ZZCand_Z2flav",1);
  oldtree->SetBranchStatus("ZZCand_Z1l1Idx",1);
  oldtree->SetBranchStatus("ZZCand_Z1l2Idx",1);
  oldtree->SetBranchStatus("ZZCand_Z2l1Idx",1);
  oldtree->SetBranchStatus("ZZCand_Z2l2Idx",1);
  oldtree->SetBranchStatus("lep_genindex",1);
  oldtree->SetBranchStatus("lep_Hindex",1);

  if(process=="qqZZ") {
    oldtree->SetBranchStatus("KFactor_EW_qqZZ_Weight",1);
    oldtree->SetBranchStatus("KFactor_QCD_qqZZ_M_Weight",1);
  }
  if(process=="ggZZ") {
    oldtree->SetBranchStatus("KFactor_QCD_ggZZ_Nominal_Weight",1);
  }

  oldtree->SetBranchStatus("puWeight",1);
  oldtree->SetBranchStatus("Generator_weight",1);
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
    // Counters of Z and GENLeps in the event
    oldtree->SetBranchStatus("nFidZ",1);
    oldtree->SetBranchStatus("nFidDressedLeps",1);
  }
  if(prod_mode == "ggH125") oldtree->SetBranchStatus("ggH_NNLOPS_Weight",1); // Additional entry for the weight in case of ggH

  //skim oldtree_failed for signal only
  if(process=="signal" || process=="AC"){
    // candTree_failed
    // Deactivate all branches
    oldtree_failed->SetBranchStatus("*",0);
    // Activate some branches only: our skim
    oldtree_failed->SetBranchStatus("event",1);
    oldtree_failed->SetBranchStatus("Counter",1);
    oldtree_failed->SetBranchStatus("puWeight",1);
    oldtree_failed->SetBranchStatus("Generator_weight",1);
    oldtree_failed->SetBranchStatus("overallEventWeight", 1);
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
    // Counters of Z and GENLeps in the event
    oldtree_failed->SetBranchStatus("nFidZ",1);
    oldtree_failed->SetBranchStatus("nFidDressedLeps",1);
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
