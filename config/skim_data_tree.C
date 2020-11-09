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

void skim_data_tree (int year = 2018){

  TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased";
  auto oldFile = TFile::Open(Form("%s/Data_%d/AllData/ZZ4lAnalysis.root", path.Data(), year));
  TTree *oldtree = (TTree*) oldFile->Get("ZZTree/candTree");

  // Deactivate all branches
  oldtree->SetBranchStatus("*",0);
  // Activate 6 branches only: our skim
  oldtree->SetBranchStatus("ZZMass",1);
  oldtree->SetBranchStatus("Z1Flav",1);
  oldtree->SetBranchStatus("Z2Flav",1);
  oldtree->SetBranchStatus("ZZPt",1);
  oldtree->SetBranchStatus("Z1Mass",1);
  oldtree->SetBranchStatus("Z2Mass",1);
  oldtree->SetBranchStatus("ZZEta",1);

  TString newpath = "/eos/user/a/atarabin/Data";
  TFile *newfile = new TFile(Form("%s/reducedTree_AllData_%d.root", newpath.Data(), year),"RECREATE");
  auto *newtree = oldtree->CloneTree(0);
  newtree->CopyEntries(oldtree);
  newtree->Write();
  newfile->Close();

  TFile *f = new TFile(Form("%s/reducedTree_AllData_%d.root", newpath.Data(), year),"UPDATE");
  TTree *T = (TTree*)f->Get("candTree");
  Short_t Z1Flav,Z2Flav;
  float ZZMass, ZZPt, ZZEta;
  float _chan, _CMS_zz4l_mass, _ZZy;
  TBranch *chan = T->Branch("chan",&_chan,"chan/F");
  TBranch *ZZy = T->Branch("ZZy",&_ZZy,"ZZy/F");
  TBranch *CMS_zz4l_mass = T->Branch("CMS_zz4l_mass",&_CMS_zz4l_mass,"CMS_zz4l_mass/F");
  T->SetBranchAddress("Z1Flav",&Z1Flav);
  T->SetBranchAddress("Z2Flav",&Z2Flav);
  T->SetBranchAddress("ZZMass",&ZZMass);
  T->SetBranchAddress("ZZPt",&ZZPt);
  T->SetBranchAddress("ZZEta",&ZZEta);
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

    ZZy->Fill();
    chan->Fill();
    CMS_zz4l_mass->Fill();
  }
  T->Print();
  T->Write("", TObject::kOverwrite);
  delete f;

  return 0;
}
