// root -l skim_data_tree\(year\)

#include<iostream>
#include<fstream>
#include<cstring>
#include<math.h>
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
#include<TMultiGraph.h>
#include<TLegend.h>
#include<TLatex.h>
#include<TStyle.h>
#include<random>
#include<algorithm>

using namespace std;

void add(TString newpath, int year, TString tree_dir){

  TFile *f = new TFile(Form("%s/reducedTree_AllData_%d.root", newpath.Data(), year),"UPDATE");
  cout << Form("%s/reducedTree_AllData_%d.root", newpath.Data(), year) << endl;
  TTree *T = (TTree*)f->Get(tree_dir);
  Short_t Z1Flav,Z2Flav;
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float _chan, _CMS_zz4l_mass, _ZZy;
  float _njets_pt30_eta2p5, _pTj1, _pTj2, _pTHj, _pTHjj, _mHj, _mHjj, _mjj, _detajj, _dphijj;
  Float_t Mj1, ETAj1, PHIj1, Mj2, ETAj2, PHIj2;
  vector<float> *JetPt = 0;
  vector<float> *JetEta = 0;
  vector<float> *JetPhi = 0;
  vector<float> *JetMass = 0;
  TBranch *chan = T->Branch("chan",&_chan,"chan/F");
  TBranch *ZZy = T->Branch("ZZy",&_ZZy,"ZZy/F");
  TBranch *CMS_zz4l_mass = T->Branch("CMS_zz4l_mass",&_CMS_zz4l_mass,"CMS_zz4l_mass/F");
  TBranch *njets_pt30_eta2p5 = T->Branch("njets_pt30_eta2p5",&_njets_pt30_eta2p5,"njets_pt30_eta2p5/F");
  TBranch *pTj1 = T->Branch("pTj1",&_pTj1,"pTj1/F");
  TBranch *pTj2 = T->Branch("pTj2",&_pTj2,"pTj2/F");
  TBranch *pTHj = T->Branch("pTHj",&_pTHj,"pTHj/F");
  TBranch *pTHjj = T->Branch("pTHjj",&_pTHjj,"pTHjj/F");
  TBranch *mHj = T->Branch("mHj",&_mHj,"mHj/F");
  TBranch *mHjj = T->Branch("mHjj",&_mHjj,"mHjj/F");
  TBranch *mjj = T->Branch("mjj",&_mjj,"mjj/F");
  TBranch *detajj = T->Branch("detajj",&_detajj,"detajj/F");
  TBranch *dphijj = T->Branch("dphijj",&_dphijj,"dphijj/F");
  T->SetBranchAddress("Z1Flav",&Z1Flav);
  T->SetBranchAddress("Z2Flav",&Z2Flav);
  T->SetBranchAddress("ZZMass",&ZZMass);
  T->SetBranchAddress("ZZPt",&ZZPt);
  T->SetBranchAddress("ZZEta",&ZZEta);
  T->SetBranchAddress("ZZPhi",&ZZPhi);
  T->SetBranchAddress("JetPt",&JetPt);
  T->SetBranchAddress("JetEta",&JetEta);
  T->SetBranchAddress("JetMass",&JetMass);
  T->SetBranchAddress("JetPhi",&JetPhi);

  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);
    if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121){
      _chan=2;
    }
    else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==169){
      _chan=1;
    }
    else{
      _chan=3;
    }
    _CMS_zz4l_mass = ZZMass;
    _ZZy = abs(log((sqrt(125*125 + ZZPt*ZZPt*cosh(ZZEta)*cosh(ZZEta))+ZZPt*sinh(ZZEta))/sqrt(125*125+ZZPt*ZZPt)));

    // njets
    _njets_pt30_eta2p5 = 0;
    for(unsigned int i=0;i<JetPt->size();i++){
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5){
        _njets_pt30_eta2p5 = _njets_pt30_eta2p5 + 1;
      }
    }

    // leading jet pT
    _pTj1 = 0;
    Mj1 = 0;
    ETAj1 = 0;
    PHIj1 = 0;
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5 && JetPt->at(i) > _pTj1) {
        _pTj1 = JetPt->at(i);
        Mj1 = JetMass->at(i);
        ETAj1 = JetEta->at(i);
        PHIj1 = JetPhi->at(i);
      }
    }
    // sub-leading jet pT
    _pTj2 = 0;
    Mj2 = 0;
    ETAj2 = 0;
    PHIj2 = 0;
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
    H.SetPtEtaPhiM(ZZPt,ZZEta,ZZPhi,ZZMass);
    j1.SetPtEtaPhiM(_pTj1,ETAj1,PHIj1,Mj1);
    j2.SetPtEtaPhiM(_pTj2,ETAj2,PHIj2,Mj2);
    _pTHj = (H+j1).Pt();
    _pTHjj = (H+j1+j2).Pt();
    _mHj = (H+j1).M();
    _mHjj = (H+j1+j2).M();
    _detajj = j1.Eta()-j2.Eta();
    _dphijj = j1.Phi()-j2.Phi();

    njets_pt30_eta2p5->Fill();
    pTj1->Fill();
    pTj2->Fill();
    pTHj->Fill();
    pTHjj->Fill();
    mHj->Fill();
    mHjj->Fill();
    detajj->Fill();
    dphijj->Fill();
    mjj->Fill();
    ZZy->Fill();
    chan->Fill();
    CMS_zz4l_mass->Fill();
  }
  T->Print();
  T->Write("", TObject::kOverwrite);
  delete f;
}


void skim_data_tree (int year = 2016){

  TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased";
  auto oldFile = TFile::Open(Form("%s/Data_%d/AllData/ZZ4lAnalysis.root", path.Data(), year));
  TTree *oldtree = (TTree*) oldFile->Get("ZZTree/candTree");
  TTree *oldtree_CR = (TTree*) oldFile->Get("CRZLLTree/candTree");

  // Deactivate all branches
  oldtree->SetBranchStatus("*",0);
  // Activate some branches only: our skim
  oldtree->SetBranchStatus("ZZMass",1);
  oldtree->SetBranchStatus("ZZPhi",1);
  oldtree->SetBranchStatus("Z1Flav",1);
  oldtree->SetBranchStatus("Z2Flav",1);
  oldtree->SetBranchStatus("ZZPt",1);
  oldtree->SetBranchStatus("Z1Mass",1);
  oldtree->SetBranchStatus("Z2Mass",1);
  oldtree->SetBranchStatus("ZZEta",1);
  oldtree->SetBranchStatus("JetPt",1);
  oldtree->SetBranchStatus("JetEta",1);
  oldtree->SetBranchStatus("JetPhi",1);
  oldtree->SetBranchStatus("JetMass",1);
  oldtree->SetBranchStatus("costhetastar",1);
  oldtree->SetBranchStatus("helcosthetaZ1",1);
  oldtree->SetBranchStatus("helcosthetaZ2",1);
  oldtree->SetBranchStatus("helphi",1);
  oldtree->SetBranchStatus("phistarZ1",1);
  oldtree->SetBranchStatus("LepLepId",1);
  oldtree->SetBranchStatus("LepEta",1);
  oldtree->SetBranchStatus("LepPt",1);


  // Deactivate all branches
  oldtree_CR->SetBranchStatus("*",0);
  // Activate some branches only: our skim
  oldtree_CR->SetBranchStatus("ZZMass",1);
  oldtree_CR->SetBranchStatus("ZZPhi",1);
  oldtree_CR->SetBranchStatus("Z1Flav",1);
  oldtree_CR->SetBranchStatus("Z2Flav",1);
  oldtree_CR->SetBranchStatus("ZZPt",1);
  oldtree_CR->SetBranchStatus("Z1Mass",1);
  oldtree_CR->SetBranchStatus("Z2Mass",1);
  oldtree_CR->SetBranchStatus("ZZEta",1);
  oldtree_CR->SetBranchStatus("JetPt",1);
  oldtree_CR->SetBranchStatus("JetEta",1);
  oldtree_CR->SetBranchStatus("JetPhi",1);
  oldtree_CR->SetBranchStatus("JetMass",1);
  oldtree_CR->SetBranchStatus("costhetastar",1);
  oldtree_CR->SetBranchStatus("helcosthetaZ1",1);
  oldtree_CR->SetBranchStatus("helcosthetaZ2",1);
  oldtree_CR->SetBranchStatus("helphi",1);
  oldtree_CR->SetBranchStatus("phistarZ1",1);
  oldtree_CR->SetBranchStatus("LepLepId",1);
  oldtree_CR->SetBranchStatus("LepEta",1);
  oldtree_CR->SetBranchStatus("LepPt",1);




  TString newpath = "/eos/user/a/atarabin/Data";
  TFile *newfile = new TFile(Form("%s/reducedTree_AllData_%d.root", newpath.Data(), year),"RECREATE");
  auto *newtree = oldtree->CloneTree(0);
  auto *newtree_CR = oldtree_CR->CloneTree(0);
  newtree->SetName("SR");
  newtree->CopyEntries(oldtree);
  newtree->Write();
  newtree_CR->SetName("CRZLL");
  newtree_CR->CopyEntries(oldtree_CR);
  newtree_CR->Write();
  newfile->Close();

  add(newpath, year, "SR");
  add(newpath, year, "CRZLL");

  return 0;
}
