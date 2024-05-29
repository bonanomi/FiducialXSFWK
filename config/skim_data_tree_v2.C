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

TFile *gConstant_g4 = new TFile("gConstant_HZZ2e2mu_g4.root");
TSpline *spline_g4 = (TSpline*) gConstant_g4->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");

TFile *gConstant_g2 = new TFile("gConstant_HZZ2e2mu_g2.root");
TSpline *spline_g2 = (TSpline*) gConstant_g2->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");

TFile *gConstant_L1 = new TFile("gConstant_HZZ2e2mu_L1.root");
TSpline *spline_L1 = (TSpline*) gConstant_L1->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1");

TFile *gConstant_L1Zgs = new TFile("gConstant_HZZ2e2mu_L1Zgs.root");
TSpline *spline_L1Zgs = (TSpline*) gConstant_L1Zgs->Get("sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs");

TFile *gConstant_bkg_2e2mu = new TFile("SmoothKDConstant_m4l_Dbkgkin_2e2mu13TeV.root");
TSpline *DbkgkinSpline2e2mu = (TSpline*) gConstant_bkg_2e2mu->Get("sp_gr_varReco_Constant_Smooth");

TFile *gConstant_bkg_4mu = new TFile("SmoothKDConstant_m4l_Dbkgkin_4mu13TeV.root");
TSpline *DbkgkinSpline4mu = (TSpline*) gConstant_bkg_4mu->Get("sp_gr_varReco_Constant_Smooth");

TFile *gConstant_bkg_4e = new TFile("SmoothKDConstant_m4l_Dbkgkin_4e13TeV.root");
TSpline *DbkgkinSpline4e = (TSpline*) gConstant_bkg_4e->Get("sp_gr_varReco_Constant_Smooth");

// -------------------------------------- Dkin discriminats -------------------------------------- //
float getDbkgkinConstant(int ZZflav, float ZZMass){ // ZZflav==id1*id2*id3*id4
  if (abs(ZZflav)==11*11*11*11 || abs(ZZflav)==2*11*11*11*11 || abs(ZZflav)==2*11*11*2*11*11) return DbkgkinSpline4e->Eval(ZZMass);
  if (abs(ZZflav)==11*11*13*13 || abs(ZZflav)==2*11*11*13*13 || abs(ZZflav)==2*11*11*2*13*13) return DbkgkinSpline2e2mu->Eval(ZZMass);
  if (abs(ZZflav)==13*13*13*13 || abs(ZZflav)==2*13*13*13*13 || abs(ZZflav)==2*13*13*2*13*13) return DbkgkinSpline4mu->Eval(ZZMass);
  std::cout << "Invalid ZZflav " << ZZflav << std::endl; assert(0); return 0;
}

float deltaphi (TLorentzVector tetra1, TLorentzVector tetra2){
  //Direction of the two jets - vectors in the lab frame
  TVector3 j1dir(tetra1.X(), tetra1.Y(), tetra1.Z());
  TVector3 j2dir(tetra2.X(), tetra2.Y(), tetra2.Z());

  //Transverse component in the xy plane
  TVector3 jt1(tetra1.X(), tetra1.Y(), 0);
  TVector3 jt2(tetra2.X(), tetra2.Y(), 0);

  //Unit vectors of the transverse components
  TVector3 jt1_norm   = jt1 * (1/jt1.Mag());
  TVector3 jt2_norm   = jt2 * (1/jt2.Mag());

  //Unit vector of the z axis
  TVector3 z(0,0,1);

  //Cross product between transverse components
  Double_t cross      = jt1_norm.Cross(jt2_norm) * z;
  Double_t cross_norm = cross * (1 / abs(cross));

  //Dot product between transverse components
  Double_t dot         = jt1_norm * jt2_norm;

  //Difference between the direction of the two jets
  Double_t diff       = (j1dir - j2dir) * z;
  Double_t diff_norm  = diff * (1 / abs(diff));

  return acos(dot) * diff_norm * cross_norm;
}

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


void add(TString newpath, TString year, TString tree_dir){

  TFile *f = new TFile(Form("%s/reducedTree_AllData_%s.root", newpath.Data(), year.Data()),"UPDATE");
  cout << Form("%s/reducedTree_AllData_%s.root", newpath.Data(), year.Data()) << endl;

  // TFile *f = new TFile(Form("%s/%s/reducedTree.root", newpath.Data(), year.Data()),"UPDATE");
  // cout << Form("%s/%s/reducedTree.root", newpath.Data(), year.Data()) << endl;

  TTree *T = (TTree*)f->Get(tree_dir);
  Short_t Z1Flav,Z2Flav;
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float _chan, _CMS_zz4l_mass, _ZZy;
  float _njets_pt30_eta2p5, _pTj1, _pTj2, _pTHj, _pTHjj, _mHj, _mHjj, _mjj, _detajj, _dphijj, _absdphijj, _absdetajj;
  Float_t Mj1, ETAj1, PHIj1, Mj2, ETAj2, PHIj2;
  Float_t _TCjmax, _TBjmax, _TCj, _TBj, _TCj1, _TBj1;
  Float_t _Dcp, _D0m, _D0hp, _Dint, _DL1, _DL1int, _DL1Zg, _DL1Zgint, _Dbkg, _Dbkg_kin;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_GG_SIG_ghg2_1_ghz4_1_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen;
  Float_t p_GG_SIG_ghg2_1_ghz2_1_JHUGen, p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen;
  Float_t p_m4l_SIG, p_m4l_BKG, p_QQB_BKG_MCFM;
  // Float_t _pTj1_eta4p7;
  Int_t _njets_pt30_eta4p7;
  vector<float> *JetPt = 0;
  vector<float> *JetEta = 0;
  vector<float> *JetPhi = 0;
  vector<float> *JetMass = 0;
  TBranch *chan = T->Branch("chan",&_chan,"chan/F");
  TBranch *CMS_zz4l_mass = T->Branch("CMS_zz4l_mass",&_CMS_zz4l_mass,"CMS_zz4l_mass/F");

  T->SetBranchAddress("Z1Flav",&Z1Flav);
  T->SetBranchAddress("Z2Flav",&Z2Flav);
  T->SetBranchAddress("ZZMass",&ZZMass);
  T->SetBranchAddress("ZZPt",&ZZPt);
  T->SetBranchAddress("ZZEta",&ZZEta);
  T->SetBranchAddress("ZZPhi",&ZZPhi);

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
    chan->Fill();
    CMS_zz4l_mass->Fill();
  }
  T->Print();
  T->Write("", TObject::kOverwrite);
  delete f;
}


void skim_data_tree_v2 (TString year = "2022"){

  TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/"+year+"/Data";
  auto oldFile = TFile::Open(Form("%s/AllData_%s.root", path.Data(), year.Data()));

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

  TString newpath = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/"+year+"/Data";
  TFile *newfile = new TFile(Form("%s/reducedTree_AllData_%s.root", newpath.Data(), year.Data()),"RECREATE");

  auto *newtree = oldtree->CloneTree(0);
  newtree->SetName("SR");
  newtree->CopyEntries(oldtree);
  newtree->Write();
  newfile->Close();

  add(newpath, year, "SR");

  return 0;
}
