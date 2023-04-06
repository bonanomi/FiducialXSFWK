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

// #include "LeptonSFHelper.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

// -------------------------------------- constants for discriminats -------------------------------------- //
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


// -------------------------------------- Muons SFs (from ReReco to UL) -------------------------------------- //
// 2016 Muons
TString fipMu_2016 = Form("final_HZZ_SF_2016UL_mupogsysts_newLoose.root");
TFile *root_file_16 = TFile::Open(fipMu_2016.Data(),"READ");
TH2D *h_Mu_SF_2016  = (TH2D*)root_file_16->Get("FINAL")->Clone();
TH2D *h_Mu_Unc_2016 = (TH2D*)root_file_16->Get("ERROR")->Clone();

// 2017 Muons
TString fipMu_2017 = Form("final_HZZ_SF_2017UL_mupogsysts_newLoose.root");
TFile *root_file_17 = TFile::Open(fipMu_2017.Data(),"READ");
TH2D *h_Mu_SF_2017  = (TH2D*)root_file_17->Get("FINAL")->Clone();
TH2D *h_Mu_Unc_2017 = (TH2D*)root_file_17->Get("ERROR")->Clone();

// 2018 Muons
TString fipMu_2018 = Form("final_HZZ_SF_2018UL_mupogsysts_newLoose.root");
TFile *root_file_18 = TFile::Open(fipMu_2018.Data(),"READ");
TH2D *h_Mu_SF_2018  = (TH2D*)root_file_18->Get("FINAL")->Clone();
TH2D *h_Mu_Unc_2018 = (TH2D*)root_file_18->Get("ERROR")->Clone();

float getSF(int year, int flav, float pt, float eta, float SCeta, bool isCrack)
{
   float SelSF = 1.0;

    //Muon SF
    if(year == 2016)
    {
       SelSF = h_Mu_SF_2016->GetBinContent(h_Mu_SF_2016->GetXaxis()->FindBin(eta),h_Mu_SF_2016->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
    }
    else if(year == 2017)
    {
       SelSF = h_Mu_SF_2017->GetBinContent(h_Mu_SF_2017->GetXaxis()->FindBin(eta),h_Mu_SF_2017->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
    }
    else if(year == 2018)
    {
       SelSF = h_Mu_SF_2018->GetBinContent(h_Mu_SF_2018->GetXaxis()->FindBin(eta),h_Mu_SF_2018->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
    }
    else {
       std::cout << "Muon SFs for " << year << " is not supported!" << std::endl;
       abort();
    }

    float SF = SelSF;

    return SF;
}

float getSFError(int year, int flav, float pt, float eta, float SCeta, bool isCrack)
{
   float RecoSF = 1.0;
   float SelSF = 1.0;

   float RecoSF_Unc = 0.0;
   float SelSF_Unc = 0.0;
   float SFError = 0.0;

   //Muon SF
   if(abs(flav) == 13 )
   {
      if(year == 2016)
      {
         SelSF = h_Mu_SF_2016->GetBinContent(h_Mu_SF_2016->GetXaxis()->FindBin(eta),h_Mu_SF_2016->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
         SelSF_Unc = h_Mu_Unc_2016->GetBinContent(h_Mu_Unc_2016->GetXaxis()->FindBin(eta),h_Mu_Unc_2016->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else if(year == 2017)
      {
         SelSF = h_Mu_SF_2017->GetBinContent(h_Mu_SF_2017->GetXaxis()->FindBin(eta),h_Mu_SF_2017->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
         SelSF_Unc = h_Mu_Unc_2017->GetBinContent(h_Mu_Unc_2017->GetXaxis()->FindBin(eta),h_Mu_Unc_2017->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else if(year == 2018)
      {
         SelSF = h_Mu_SF_2018->GetBinContent(h_Mu_SF_2018->GetXaxis()->FindBin(eta),h_Mu_SF_2018->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
         SelSF_Unc = h_Mu_Unc_2018->GetBinContent(h_Mu_Unc_2018->GetXaxis()->FindBin(eta),h_Mu_Unc_2018->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else {
         std::cout << "Muon SFs for " << year << " is not supported!" << std::endl;
         abort();
      }

      SFError = SelSF_Unc/SelSF; // assume full correlation between different muons (and uncorrelated reco and sel uncertainties)
   }

   return SFError;
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

float mass_lep(int flavour){
  if((abs(flavour)) == 11) return 0.0005109989461;
  else if ((abs(flavour)) == 13) return 0.1056583745;
  else if ((abs(flavour)) == 15) return 1.77686;
  // else if ((abs(flavour)) == 0) return 0;
  return 0;
}


//------------------------------------------------------------------
void add(TString input_dir, TString year, TString prod_mode, TString process, bool t_failed=true){
  // Add additional branches
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path;
  if(process!="AC") {
    new_full_path = Form("%s/%sUL/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  }
  else {
    new_full_path = Form("%s/AC%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  }
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
  Float_t _GEN_Dcp, _GEN_D0m, _GEN_D0hp, _GEN_Dint, _GEN_DL1, _GEN_DL1int, _GENrapidity4lAbs, _GEN_DL1Zg, _GEN_DL1Zgint;
  Float_t GENmass4l,GENpT4l,GENeta4l,GENphi4l;
  Short_t GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  Float_t GENrapidity4l;
  Short_t _GENnjets_pt30_eta2p5;
  Short_t _GENnjets_pt30_eta4p7;
  Float_t _GENpTj1, _GENpTj2, _GENpTHj, _GENpTHjj, _GENmHj, _GENmHjj, _GENmjj, _GENdetajj, _GENdphijj, _GENabsdetajj, _GENabsdphijj;
  Float_t GENMj1, GENETAj1, GENPHIj1, GENMj2, GENETAj2, GENPHIj2;
  Float_t _GENTCj1, _GENTBj1;
  Float_t _GENTCj, _GENTBj, _GENTCjmax, _GENTBjmax;
  bool _passedFullSelection, passedFiducialSelection_bbf;
  Float_t _weight, _SFcorr;
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
  TBranch *passedFullSelection = T->Branch("passedFullSelection",&_passedFullSelection,"passedFullSelection/B");
  TBranch *weight = T->Branch("weight",&_weight,"weight/F");
  TBranch *SFcorr = T->Branch("SFcorr",&_SFcorr,"SFcorr/F");
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
  TBranch *GENrapidity4lAbs = T->Branch("GENrapidity4lAbs",&_GENrapidity4lAbs,"GENrapidity4lAbs/F");

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
    T->SetBranchAddress("GENrapidity4l",&GENrapidity4l);
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
  float _ZZy,ZZPt,ZZEta,ZZPhi,ZZMass, L1prefiringWeight, overallEventWeight, dataMCWeight;
  Short_t Z1Flav, Z2Flav, ZZFlav;
  Short_t _njets_pt30_eta2p5, _njets_pt30_eta4p7;
  Float_t _pTj1, _pTj2, _pTHj, _pTHjj, _mHj, _mHjj, _mjj, _detajj, _dphijj, _absdetajj, _absdphijj;
  Float_t _Mj1, _ETAj1, _PHIj1, _Mj2, _ETAj2, _PHIj2;
  Float_t _TCj, _TCj1, _TCjmax, _TBj1, _TBj, _TBjmax;
  Float_t tc, tcj, yj;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_GG_SIG_ghg2_1_ghz4_1_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen;
  Float_t p_GG_SIG_ghg2_1_ghz2_1_JHUGen, p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen;
  Float_t p_QQB_BKG_MCFM, p_m4l_SIG, p_m4l_BKG;
  Float_t _Dcp, _D0m, _D0hp, _Dint, _DL1, _DL1int, _DL1Zg, _DL1Zgint, _Dbkg, _Dbkg_kin;
  vector<float> _LepSF_new, _LepSF_Unc_new;
  vector<float> *LepSF = 0;
  vector<float> *LepPt = 0;
  vector<float> *LepPhi = 0;
  vector<float> *LepEta = 0;
  vector<int> *LepLepId = 0;
  vector<bool> *LepisCrack = 0;
  vector<float> *ExtraLepPt = 0;
  vector<float> *ExtraLepEta = 0;
  vector<float> *ExtraLepPhi = 0;
  vector<int> *ExtraLepLepId = 0;
  vector<int> _lep_genindex, _lep_Hindex;
  vector<float> *JetPt = 0;
  vector<float> *JetEta = 0;
  vector<float> *JetMass = 0;
  vector<float> *JetPhi = 0;
  // // JetSigma
  // vector<float> *JetSigma_Total = 0;
  // vector<float> *JetSigma_Abs = 0;
  // vector<float> *JetSigma_Abs_year = 0;
  // vector<float> *JetSigma_BBEC1 = 0;
  // vector<float> *JetSigma_BBEC1_year = 0;
  // vector<float> *JetSigma_EC2 = 0;
  // vector<float> *JetSigma_EC2_year = 0;
  // vector<float> *JetSigma_FlavQCD = 0;
  // vector<float> *JetSigma_HF = 0;
  // vector<float> *JetSigma_HF_year = 0;
  // vector<float> *JetSigma_RelBal = 0;
  // vector<float> *JetSigma_RelSample_year = 0;
  // JetJESUp
  // vector<float> *_JetPt_JESUp_Total = 0;
  // vector<float> *_JetPt_JESUp_Abs = 0;
  // vector<float> *_JetPt_JESUp_Abs_year = 0;
  // vector<float> *_JetPt_JESUp_BBEC1 = 0;
  // vector<float> *_JetPt_JESUp_BBEC1_year = 0;
  // vector<float> *_JetPt_JESUp_EC2 = 0;
  // vector<float> *_JetPt_JESUp_EC2_year = 0;
  // vector<float> *_JetPt_JESUp_FlavQCD = 0;
  // vector<float> *_JetPt_JESUp_HF = 0;
  // vector<float> *_JetPt_JESUp_HF_year = 0;
  // vector<float> *_JetPt_JESUp_RelBal = 0;
  // vector<float> *_JetPt_JESUp_RelSample_year = 0;
  // JetJESDn
  // vector<float> *_JetPt_JESDown_Total = 0;
  // vector<float> *_JetPt_JESDown_Abs = 0;
  // vector<float> *_JetPt_JESDown_Abs_year = 0;
  // vector<float> *_JetPt_JESDown_BBEC1 = 0;
  // vector<float> *_JetPt_JESDown_BBEC1_year = 0;
  // vector<float> *_JetPt_JESDown_EC2 = 0;
  // vector<float> *_JetPt_JESDown_EC2_year = 0;
  // vector<float> *_JetPt_JESDown_FlavQCD = 0;
  // vector<float> *_JetPt_JESDown_HF = 0;
  // vector<float> *_JetPt_JESDown_HF_year = 0;
  // vector<float> *_JetPt_JESDown_RelBal = 0;
  // vector<float> *_JetPt_JESDown_RelSample_year = 0;
  TBranch *ZZy = T->Branch("ZZy",&_ZZy,"ZZy/F");
  TBranch *lep_genindex = T->Branch("lep_genindex",&_lep_genindex);
  TBranch *lep_Hindex = T->Branch("lep_Hindex",&_lep_Hindex);
  TBranch *LepSF_new = T->Branch("LepSF_new",&_LepSF_new);
  TBranch *LepSF_Unc_new = T->Branch("LepSF_Unc_new",&_LepSF_Unc_new);
  TBranch *njets_pt30_eta2p5 = T->Branch("njets_pt30_eta2p5",&_njets_pt30_eta2p5,"njets_pt30_eta2p5/S");
  TBranch *njets_pt30_eta4p7 = T->Branch("njets_pt30_eta4p7",&_njets_pt30_eta4p7,"njets_pt30_eta4p7/S");
  TBranch *pTj1 = T->Branch("pTj1",&_pTj1,"pTj1/F");
  TBranch *Mj1 = T->Branch("Mj1",&_Mj1,"Mj1/F");
  TBranch *ETAj1 = T->Branch("ETAj1",&_ETAj1,"ETAj1/F");
  TBranch *PHIj1 = T->Branch("PHIj1",&_PHIj1,"PHIj1/F");
  TBranch *pTj2 = T->Branch("pTj2",&_pTj2,"pTj2/F");
  TBranch *Mj2 = T->Branch("Mj2",&_Mj2,"Mj2/F");
  TBranch *ETAj2 = T->Branch("ETAj2",&_ETAj2,"ETAj2/F");
  TBranch *PHIj2 = T->Branch("PHIj2",&_PHIj2,"PHIj2/F");
  TBranch *pTHj = T->Branch("pTHj",&_pTHj,"pTHj/F");
  TBranch *pTHjj = T->Branch("pTHjj",&_pTHjj,"pTHjj/F");
  TBranch *mHj = T->Branch("mHj",&_mHj,"mHj/F");
  TBranch *mHjj = T->Branch("mHjj",&_mHjj,"mHjj/F");
  TBranch *mjj = T->Branch("mjj",&_mjj,"mjj/F");
  TBranch *detajj = T->Branch("detajj",&_detajj,"detajj/F");
  TBranch *dphijj = T->Branch("dphijj",&_dphijj,"dphijj/F");
  TBranch *absdetajj = T->Branch("absdetajj",&_absdetajj,"absdetajj/F");
  TBranch *absdphijj = T->Branch("absdphijj",&_absdphijj,"absdphijj/F");
  TBranch *TCj = T->Branch("TCj",&_TCj,"TCj/F");
  TBranch *TCjmax = T->Branch("TCjmax",&_TCjmax,"TCjmax/F");
  TBranch *TCj1 = T->Branch("TCj1",&_TCj1,"TCj1/F");
  TBranch *TBj1 = T->Branch("TBj1",&_TBj1,"TBj1/F");
  TBranch *TBj = T->Branch("TBj",&_TBj,"TBj/F");
  TBranch *TBjmax = T->Branch("TBjmax", &_TBjmax, "TBjmax/F");
  TBranch *Dbkg = T->Branch("Dbkg",&_Dbkg,"Dbkg/F");
  TBranch *Dbkg_kin = T->Branch("Dbkg_kin",&_Dbkg_kin,"Dbkg_kin/F");
  TBranch *Dcp = T->Branch("Dcp",&_Dcp,"Dcp/F");
  TBranch *D0m = T->Branch("D0m",&_D0m,"D0m/F");
  TBranch *Dint = T->Branch("Dint",&_Dint,"Dint/F");
  TBranch *D0hp = T->Branch("D0hp",&_D0hp,"D0hp/F");
  TBranch *DL1 = T->Branch("DL1",&_DL1,"DL1/F");
  TBranch *DL1int = T->Branch("DL1int",&_DL1int,"DL1int/F");
  TBranch *DL1Zg = T->Branch("DL1Zg",&_DL1Zg,"DL1Zg/F");
  TBranch *DL1Zgint = T->Branch("DL1Zgint",&_DL1Zgint,"DL1Zgint/F");
  // JetJESUp
  // TBranch *JetPt_JESUp_Total = T->Branch("JetPt_JESUp_Total",&_JetPt_JESUp_Total);
  // TBranch *JetPt_JESUp_Abs = T->Branch("JetPt_JESUp_Abs",&_JetPt_JESUp_Abs);
  // TBranch *JetPt_JESUp_Abs_year = T->Branch("JetPt_JESUp_Abs_year",&_JetPt_JESUp_Abs_year);
  // TBranch *JetPt_JESUp_BBEC1 = T->Branch("JetPt_JESUp_BBEC1",&_JetPt_JESUp_BBEC1);
  // TBranch *JetPt_JESUp_BBEC1_year = T->Branch("JetPt_JESUp_BBEC1_year",&_JetPt_JESUp_BBEC1_year);
  // TBranch *JetPt_JESUp_EC2 = T->Branch("JetPt_JESUp_EC2",&_JetPt_JESUp_EC2);
  // TBranch *JetPt_JESUp_EC2_year = T->Branch("JetPt_JESUp_EC2_year",&_JetPt_JESUp_EC2_year);
  // TBranch *JetPt_JESUp_FlavQCD = T->Branch("JetPt_JESUp_FlavQCD",&_JetPt_JESUp_FlavQCD);
  // TBranch *JetPt_JESUp_HF = T->Branch("JetPt_JESUp_HF",&_JetPt_JESUp_HF);
  // TBranch *JetPt_JESUp_HF_year = T->Branch("JetPt_JESUp_HF_year",&_JetPt_JESUp_HF_year);
  // TBranch *JetPt_JESUp_RelBal = T->Branch("JetPt_JESUp_RelBal",&_JetPt_JESUp_RelBal);
  // TBranch *JetPt_JESUp_RelSample_year = T->Branch("JetPt_JESUp_RelSample_year",&_JetPt_JESUp_RelSample_year);
  // JetJESDn
  // TBranch *JetPt_JESDown_Total = T->Branch("JetPt_JESDown_Total",&_JetPt_JESDown_Total);
  // TBranch *JetPt_JESDown_Abs = T->Branch("JetPt_JESDown_Abs",&_JetPt_JESDown_Abs);
  // TBranch *JetPt_JESDown_Abs_year = T->Branch("JetPt_JESDown_Abs_year",&_JetPt_JESDown_Abs_year);
  // TBranch *JetPt_JESDown_BBEC1 = T->Branch("JetPt_JESDown_BBEC1",&_JetPt_JESDown_BBEC1);
  // TBranch *JetPt_JESDown_BBEC1_year = T->Branch("JetPt_JESDown_BBEC1_year",&_JetPt_JESDown_BBEC1_year);
  // TBranch *JetPt_JESDown_EC2 = T->Branch("JetPt_JESDown_EC2",&_JetPt_JESDown_EC2);
  // TBranch *JetPt_JESDown_EC2_year = T->Branch("JetPt_JESDown_EC2_year",&_JetPt_JESDown_EC2_year);
  // TBranch *JetPt_JESDown_FlavQCD = T->Branch("JetPt_JESDown_FlavQCD",&_JetPt_JESDown_FlavQCD);
  // TBranch *JetPt_JESDown_HF = T->Branch("JetPt_JESDown_HF",&_JetPt_JESDown_HF);
  // TBranch *JetPt_JESDown_HF_year = T->Branch("JetPt_JESDown_HF_year",&_JetPt_JESDown_HF_year);
  // TBranch *JetPt_JESDown_RelBal = T->Branch("JetPt_JESDown_RelBal",&_JetPt_JESDown_RelBal);
  // TBranch *JetPt_JESDown_RelSample_year = T->Branch("JetPt_JESDown_RelSample_year",&_JetPt_JESDown_RelSample_year);
  if (!t_failed) {
    T->SetBranchAddress("ZZMass",&ZZMass);
    T->SetBranchAddress("Z1Flav",&Z1Flav);
    T->SetBranchAddress("Z2Flav",&Z2Flav);
    T->SetBranchAddress("ZZPt",&ZZPt);
    T->SetBranchAddress("ZZEta",&ZZEta);
    T->SetBranchAddress("ZZPhi",&ZZPhi);
    T->SetBranchAddress("LepSF",&LepSF);
    T->SetBranchAddress("LepPt",&LepPt);
    T->SetBranchAddress("LepPhi",&LepPhi);
    T->SetBranchAddress("LepEta",&LepEta);
    T->SetBranchAddress("LepLepId",&LepLepId);
    T->SetBranchAddress("LepisCrack",&LepisCrack);
    T->SetBranchAddress("ExtraLepPt",&ExtraLepPt);
    T->SetBranchAddress("ExtraLepEta",&ExtraLepEta);
    T->SetBranchAddress("ExtraLepPhi",&ExtraLepPhi);
    T->SetBranchAddress("ExtraLepLepId",&ExtraLepLepId);
    T->SetBranchAddress("JetPt",&JetPt);
    T->SetBranchAddress("JetEta",&JetEta);
    T->SetBranchAddress("JetMass",&JetMass);
    T->SetBranchAddress("JetPhi",&JetPhi);
    T->SetBranchAddress("overallEventWeight",&overallEventWeight);
    T->SetBranchAddress("dataMCWeight",&dataMCWeight);
    T->SetBranchAddress("L1prefiringWeight",&L1prefiringWeight);
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
    T->SetBranchAddress("p_QQB_BKG_MCFM",&p_QQB_BKG_MCFM);
    T->SetBranchAddress("p_m4l_SIG",&p_m4l_SIG);
    T->SetBranchAddress("p_m4l_BKG",&p_m4l_BKG);

    // // JetSigma
    // T->SetBranchAddress("JetSigma_Total",&JetSigma_Total);
    // T->SetBranchAddress("JetSigma_Abs",&JetSigma_Abs);
    // T->SetBranchAddress("JetSigma_Abs_year",&JetSigma_Abs_year);
    // T->SetBranchAddress("JetSigma_BBEC1",&JetSigma_BBEC1);
    // T->SetBranchAddress("JetSigma_BBEC1_year",&JetSigma_BBEC1_year);
    // T->SetBranchAddress("JetSigma_EC2",&JetSigma_EC2);
    // T->SetBranchAddress("JetSigma_EC2_year",&JetSigma_EC2_year);
    // T->SetBranchAddress("JetSigma_FlavQCD",&JetSigma_FlavQCD);
    // T->SetBranchAddress("JetSigma_HF",&JetSigma_HF);
    // T->SetBranchAddress("JetSigma_HF_year",&JetSigma_HF_year);
    // T->SetBranchAddress("JetSigma_RelBal",&JetSigma_RelBal);
    // T->SetBranchAddress("JetSigma_RelSample_year",&JetSigma_RelSample_year);
  }

  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);

    if(process=="signal" || process=="AC"){
      _passedFullSelection = true;
      if (t_failed) {
  	     _passedFullSelection = false;
      }

      // H+Njets GENvariables
      // If the event does not pass fiducial selections the variable is equal to -99
      _GENpTHj = -99;
      _GENpTHjj = -99;
      _GENmHj = -99;
      _GENmHjj = -99;
      _GENdetajj = -99;
      _GENabsdetajj = -99;
      _GENmjj = -99;
      _GENdphijj = -99;
      _GENabsdphijj = -99;
      _GENTCjmax = -99; _GENTBjmax = -99;
      _GENTCj = -99; _GENTBj = -99;
      _GENTCj1 = -99; _GENTBj1 = -99;
      _GENpTj1 = -99; _GENpTj2 = -99;
      _GENrapidity4lAbs = -99;
      if(passedFiducialSelection_bbf){

        _GENrapidity4lAbs = abs(GENrapidity4l);

        _GENnjets_pt30_eta2p5 = GENjetsPt_pt30_eta2p5->size();
        _GENnjets_pt30_eta4p7 = GENjetsPt_pt30_eta4p7->size();

        // leading GENjet pT
        _GENpTj1 = -1;
        GENMj1 = 0;
        GENETAj1 = 0;
        GENPHIj1 = 0;
        for (unsigned int i = 0; i < GENjetsPt_pt30_eta4p7->size(); ++i){
          if(GENjetsPt_pt30_eta4p7->at(i) > _GENpTj1) {
            _GENpTj1 = GENjetsPt_pt30_eta4p7->at(i);
            GENMj1 = GENjetsMass_pt30_eta4p7->at(i);
            GENETAj1 = GENjetsEta_pt30_eta4p7->at(i);
            GENPHIj1 = GENjetsPhi_pt30_eta4p7->at(i);
          }
        }

        // sub-leading GENjet pT
        _GENpTj2 = -1;
        GENMj2 = 0;
        GENETAj2 = 0;
        GENPHIj2 = 0;
        for (unsigned int i = 0; i < GENjetsPt_pt30_eta4p7->size(); ++i){
          if(GENjetsPt_pt30_eta4p7->at(i) > _GENpTj2 && _GENpTj1 != GENjetsPt_pt30_eta4p7->at(i)) {
            _GENpTj2 = GENjetsPt_pt30_eta4p7->at(i);
            GENMj2 = GENjetsMass_pt30_eta4p7->at(i);
            GENETAj2 = GENjetsEta_pt30_eta4p7->at(i);
            GENPHIj2 = GENjetsPhi_pt30_eta4p7->at(i);
          }
        }

        TLorentzVector GENH;
        TLorentzVector GENj1;
        TLorentzVector GENj2;
        GENH.SetPtEtaPhiM(GENpT4l,GENeta4l,GENphi4l,GENmass4l);
        GENj1.SetPtEtaPhiM(_GENpTj1,GENETAj1,GENPHIj1,GENMj1);
        GENj2.SetPtEtaPhiM(_GENpTj2,GENETAj2,GENPHIj2,GENMj2);
        _GENTCj1 = sqrt(_GENpTj1*_GENpTj1 + GENMj1*GENMj1)/(2*cosh(GENj1.Rapidity() - GENH.Rapidity()));
        _GENTBj1 = sqrt(_GENpTj1*_GENpTj1 + GENMj1*GENMj1)*exp(-1*abs(GENj1.Rapidity() - GENH.Rapidity()));

        _GENTCjmax = -1; _GENTBjmax = -1;
        for (unsigned int i = 0; i < GENjetsPt_pt30_eta4p7->size(); ++i){
          TLorentzVector theJet;
          theJet.SetPtEtaPhiM(GENjetsPt_pt30_eta4p7->at(i), GENjetsEta_pt30_eta4p7->at(i), GENjetsPhi_pt30_eta4p7->at(i), GENjetsMass_pt30_eta4p7->at(i));
          _GENTCj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*cosh(theJet.Rapidity() - GENH.Rapidity()));
          _GENTBj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*exp(-1*abs(theJet.Rapidity() - GENH.Rapidity()));
          if (_GENTCj > _GENTCjmax) _GENTCjmax = _GENTCj;
          if (_GENTBj > _GENTBjmax) _GENTBjmax = _GENTBj;
        }

        if(_GENpTj1>0){
          _GENpTHj = (GENH+GENj1).Pt();
          _GENmHj = (GENH+GENj1).M();
        }else{
          _GENpTHj = -1;
          _GENmHj = -1;
        }

        if(_GENpTj2>0){
          _GENpTHjj = (GENH+GENj1+GENj2).Pt();
          _GENmHjj = (GENH+GENj1+GENj2).M();
          _GENdetajj = GENj1.Eta()-GENj2.Eta();
          _GENabsdetajj = abs(_GENdetajj);
          _GENmjj = (GENj1+GENj2).M();
          _GENdphijj = deltaphi(GENj1, GENj2);
          _GENabsdphijj = abs(_GENdphijj);
        }else{
          _GENpTHjj = -1;
          _GENmHjj = -1;
          _GENdetajj = -1;
          _GENabsdetajj = -1;
          _GENmjj = -1;
          _GENdphijj = -10;
          _GENabsdphijj = -1;
        }
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
        _GEN_Dcp = -99;
        _GEN_D0m = -99;
        _GEN_Dint = -99;
        _GEN_D0hp = -99;
        _GEN_DL1 = -99;
        _GEN_DL1int = -99;
        _GEN_DL1Zg = -99;
        _GEN_DL1Zgint = -99;
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

      GENrapidity4lAbs->Fill();


    }
    if (t_failed) continue; // From now on reco-only variables

    // leading jet pT
    _pTj1 = -1; _Mj1 = 0; _ETAj1 = 0; _PHIj1 = 0;
    // _pTj1_eta4p7 = -1;

    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7 && JetPt->at(i) > _pTj1) {
        _pTj1 = JetPt->at(i);
        _Mj1 = JetMass->at(i);
        _ETAj1 = JetEta->at(i);
        _PHIj1 = JetPhi->at(i);
      }

    }
    // sub-leading jet pT
    _pTj2 = -1; _Mj2 = 0; _ETAj2 = 0; _PHIj2 = 0;
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7 && JetPt->at(i) > _pTj2 && _pTj1 != JetPt->at(i)) {
        _pTj2 = JetPt->at(i);
        _Mj2 = JetMass->at(i);
        _ETAj2 = JetEta->at(i);
        _PHIj2 = JetPhi->at(i);
      }
    }

    // H+Njets variables
    TLorentzVector H;
    TLorentzVector j1;
    TLorentzVector j2;
    H.SetPtEtaPhiM(ZZPt,ZZEta,ZZPhi,ZZMass);
    j1.SetPtEtaPhiM(_pTj1,_ETAj1,_PHIj1,_Mj1);
    j2.SetPtEtaPhiM(_pTj2,_ETAj2,_PHIj2,_Mj2);
    if (_pTj1 > 0) { // Variables that can be built only if there is a leading jet
      _pTHj = (H+j1).Pt();
      _mHj = (H+j1).M();
    }
    else {
      _pTHj = -1;
      _mHj = -1;
    }
    if (_pTj2 > 0) { // Variables that can be built only if there is a leading jet and a subleading jet (the latter exists only if the former does)
      _pTHjj = (H+j1+j2).Pt();
      _mHjj = (H+j1+j2).M();
      _detajj = j1.Eta()-j2.Eta();
      _absdetajj = abs(_detajj);
      _mjj = (j1+j2).M();
      _dphijj = deltaphi(j1, j2);
      _absdphijj = abs(_dphijj);
    }
    else {
      _pTHjj = -1;
      _mHjj = -1;
      _detajj = -99;
      _absdetajj = -99;
      _mjj = -1;
      _dphijj = -99;
      _absdphijj = -99;
    }

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

    _weight = overallEventWeight * L1prefiringWeight;

    float LepSF_new_tmp = 1.0;
    float updatedSF = 1.0;
    for(int lepl = 0; lepl < 4; ++lepl){
      float lid = LepLepId->at(lepl);
      if(abs(lid) == 11) {
        _LepSF_new.push_back(LepSF->at(lepl));
        updatedSF *= LepSF->at(lepl);
        continue;
      } //The correction is just for muons
      float lpt = LepPt->at(lepl);
      float leta = LepEta->at(lepl);
      bool isCrack = LepisCrack->at(lepl);
      int year_int = year.Atoi();
      LepSF_new_tmp = getSF(year_int, lid, lpt, leta, leta, isCrack);
      _LepSF_new.push_back(LepSF_new_tmp);
      // _LepSF_Unc_new = getSFError(year_int, lid, lpt, leta, leta, isCrack);
      updatedSF *= LepSF_new_tmp;
    }
    if (updatedSF == 1.0) {
      updatedSF = dataMCWeight;
    }
    _SFcorr = updatedSF/dataMCWeight;

    Float_t mj1 = _Mj1; //Maybe tmp, I do not know if Fill make the variable empty

    pTj1->Fill();
    ETAj1->Fill();
    PHIj1->Fill();
    Mj1->Fill();
    pTj2->Fill();
    ETAj2->Fill();
    PHIj2->Fill();
    Mj2->Fill();
    pTHj->Fill();
    pTHjj->Fill();
    mHj->Fill();
    mHjj->Fill();
    detajj->Fill();
    dphijj->Fill();
    absdetajj->Fill();
    absdphijj->Fill();
    mjj->Fill();
    njets_pt30_eta2p5->Fill();
    njets_pt30_eta4p7->Fill();
    weight->Fill();
    SFcorr->Fill();
    LepSF_new->Fill();
    LepSF_Unc_new->Fill();

    _LepSF_new.clear();
    _LepSF_Unc_new.clear();


    // Reco-rapidity
    // _ZZy = abs(log((sqrt(125*125 + ZZPt*ZZPt*cosh(ZZEta)*cosh(ZZEta))+ZZPt*sinh(ZZEta))/sqrt(125*125+ZZPt*ZZPt)));
    _ZZy = abs(H.Rapidity());
    ZZy->Fill();

    _TCj = -1; _TBj = -1;
    _TCjmax = -1; _TBjmax = -1;
    _TCj1 = -1; _TBj1 = -1;
    _TCj1 = sqrt(_pTj1*_pTj1 + mj1*mj1)/(2*cosh(j1.Rapidity() - H.Rapidity()));
    _TBj1 = sqrt(_pTj1*_pTj1 + mj1*mj1)*exp(-1*abs(j1.Rapidity() - H.Rapidity()));
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      // Define TCj for all jets with pT > pT_cut
       if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7) {
          TLorentzVector theJet;
          theJet.SetPtEtaPhiM(JetPt->at(i), JetEta->at(i), JetPhi->at(i), JetMass->at(i));
          _TCj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*cosh(theJet.Rapidity() - H.Rapidity())); //theJet.E());
          _TBj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*exp(-1*abs(theJet.Rapidity() - H.Rapidity()));
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
    ZZFlav = Z1Flav * Z2Flav;
    _Dbkg_kin = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav, ZZMass));
    _Dbkg     = p_GG_SIG_ghg2_1_ghz1_1_JHUGen*p_m4l_SIG/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen*p_m4l_SIG + p_m4l_BKG*p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav, ZZMass));

    _D0m = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + (p_GG_SIG_ghg2_1_ghz4_1_JHUGen * pow(spline_g4->Eval(ZZMass),2)));
    _Dcp = p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * p_GG_SIG_ghg2_1_ghz4_1_JHUGen));

    _D0hp = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + (p_GG_SIG_ghg2_1_ghz2_1_JHUGen * pow(spline_g2->Eval(ZZMass),2)));
    _Dint = p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * p_GG_SIG_ghg2_1_ghz2_1_JHUGen));

    _DL1 = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(spline_L1->Eval(ZZMass),2)));
    _DL1int = (p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen/1e4) / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * (p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8)));

    _DL1Zg = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (p_GG_SIG_ghg2_1_ghz1_1_JHUGen + ((p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(spline_L1Zgs->Eval(ZZMass),2)));
    _DL1Zgint = (p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen/1e4) / (2 * sqrt(p_GG_SIG_ghg2_1_ghz1_1_JHUGen * (p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8)));

    Dbkg->Fill();
    Dbkg_kin->Fill();
    Dcp->Fill();
    D0m->Fill();
    D0hp->Fill();
    Dint->Fill();
    DL1->Fill();
    DL1int->Fill();
    DL1Zg->Fill();
    DL1Zgint->Fill();

    // // JES variations
    // for(unsigned int i=0; i<JetPt->size(); i++){
    //
    //   _JetPt_JESUp_Total->push_back(JetPt->at(i) * (1.0 + JetSigma_Total->at(i)));
    //   _JetPt_JESDown_Total->push_back(JetPt->at(i) * (1.0 - JetSigma_Total->at(i)));
    //
    //   _JetPt_JESUp_Abs->push_back(JetPt->at(i) * (1.0 + JetSigma_Abs->at(i)));
    //   _JetPt_JESDown_Abs->push_back(JetPt->at(i) * (1.0 - JetSigma_Abs->at(i)));
    //
    //   _JetPt_JESUp_Abs_year->push_back(JetPt->at(i) * (1.0 + JetSigma_Abs_year->at(i)));
    //   _JetPt_JESDown_Abs_year->push_back(JetPt->at(i) * (1.0 - JetSigma_Abs_year->at(i)));
    //
    //   _JetPt_JESUp_BBEC1->push_back(JetPt->at(i) * (1.0 + JetSigma_BBEC1->at(i)));
    //   _JetPt_JESDown_BBEC1->push_back(JetPt->at(i) * (1.0 - JetSigma_BBEC1->at(i)));
    //
    //   _JetPt_JESUp_BBEC1_year->push_back(JetPt->at(i) * (1.0 + JetSigma_BBEC1_year->at(i)));
    //   _JetPt_JESDown_BBEC1_year->push_back(JetPt->at(i) * (1.0 - JetSigma_BBEC1_year->at(i)));
    //
    //   _JetPt_JESUp_EC2->push_back(JetPt->at(i) * (1.0 + JetSigma_EC2->at(i)));
    //   _JetPt_JESDown_EC2->push_back(JetPt->at(i) * (1.0 - JetSigma_EC2->at(i)));
    //
    //   _JetPt_JESUp_EC2_year->push_back(JetPt->at(i) * (1.0 + JetSigma_EC2_year->at(i)));
    //   _JetPt_JESDown_EC2_year->push_back(JetPt->at(i) * (1.0 - JetSigma_EC2_year->at(i)));
    //
    //   _JetPt_JESUp_FlavQCD->push_back(JetPt->at(i) * (1.0 + JetSigma_FlavQCD->at(i)));
    //   _JetPt_JESDown_FlavQCD->push_back(JetPt->at(i) * (1.0 - JetSigma_FlavQCD->at(i)));
    //
    //   _JetPt_JESUp_HF->push_back(JetPt->at(i) * (1.0 + JetSigma_HF->at(i)));
    //   _JetPt_JESDown_HF->push_back(JetPt->at(i) * (1.0 - JetSigma_HF->at(i)));
    //
    //   _JetPt_JESUp_HF_year->push_back(JetPt->at(i) * (1.0 + JetSigma_HF_year->at(i)));
    //   _JetPt_JESDown_HF_year->push_back(JetPt->at(i) * (1.0 - JetSigma_HF_year->at(i)));
    //
    //   _JetPt_JESUp_RelBal->push_back(JetPt->at(i) * (1.0 + JetSigma_RelBal->at(i)));
    //   _JetPt_JESDown_RelBal->push_back(JetPt->at(i) * (1.0 - JetSigma_RelBal->at(i)));
    //
    //   _JetPt_JESUp_RelSample_year->push_back(JetPt->at(i) * (1.0 + JetSigma_RelSample_year->at(i)));
    //   _JetPt_JESDown_RelSample_year->push_back(JetPt->at(i) * (1.0 - JetSigma_RelSample_year->at(i)));
    // }
    // JetPt_JESUp_Total->Fill();
    // JetPt_JESDown_Total->Fill();
    // JetPt_JESUp_Abs->Fill();
    // JetPt_JESDown_Abs->Fill();
    // JetPt_JESUp_Abs_year->Fill();
    // JetPt_JESDown_Abs_year->Fill();
    // JetPt_JESUp_BBEC1->Fill();
    // JetPt_JESDown_BBEC1->Fill();
    // JetPt_JESUp_BBEC1_year->Fill();
    // JetPt_JESDown_BBEC1_year->Fill();
    // JetPt_JESUp_EC2->Fill();
    // JetPt_JESDown_EC2->Fill();
    // JetPt_JESUp_EC2_year->Fill();
    // JetPt_JESDown_EC2_year->Fill();
    // JetPt_JESUp_FlavQCD->Fill();
    // JetPt_JESDown_FlavQCD->Fill();
    // JetPt_JESUp_HF->Fill();
    // JetPt_JESDown_HF->Fill();
    // JetPt_JESUp_HF_year->Fill();
    // JetPt_JESDown_HF_year->Fill();
    // JetPt_JESUp_RelBal->Fill();
    // JetPt_JESDown_RelBal->Fill();
    // JetPt_JESUp_RelSample_year->Fill();
    // JetPt_JESDown_RelSample_year->Fill();
    //
    // _JetPt_JESUp_Total->clear();
    // _JetPt_JESDown_Total->clear();
    // _JetPt_JESUp_Abs->clear();
    // _JetPt_JESDown_Abs->clear();
    // _JetPt_JESUp_Abs_year->clear();
    // _JetPt_JESDown_Abs_year->clear();
    // _JetPt_JESUp_BBEC1->clear();
    // _JetPt_JESDown_BBEC1->clear();
    // _JetPt_JESUp_BBEC1_year->clear();
    // _JetPt_JESDown_BBEC1_year->clear();
    // _JetPt_JESUp_EC2->clear();
    // _JetPt_JESDown_EC2->clear();
    // _JetPt_JESUp_EC2_year->clear();
    // _JetPt_JESDown_EC2_year->clear();
    // _JetPt_JESUp_FlavQCD->clear();
    // _JetPt_JESDown_FlavQCD->clear();
    // _JetPt_JESUp_HF->clear();
    // _JetPt_JESDown_HF->clear();
    // _JetPt_JESUp_HF_year->clear();
    // _JetPt_JESDown_HF_year->clear();
    // _JetPt_JESUp_RelBal->clear();
    // _JetPt_JESDown_RelBal->clear();
    // _JetPt_JESUp_RelSample_year->clear();
    // _JetPt_JESDown_RelSample_year->clear();


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
  T->Write(0, TObject::kOverwrite);
  delete f;
  return;
}

//---------------------------------------------------------- MAIN ----------------------------------------------------------
void skim_MC_tree_v2 (TString prod_mode = "VBFH125", TString year = "2018"){

  //process flag
  TString process;
  if(prod_mode=="ZZTo4l") {
    process = "qqZZ";
  }
  else if(prod_mode=="ggTo2e2mu_Contin_MCFM701" || prod_mode=="ggTo2e2tau_Contin_MCFM701" || prod_mode=="ggTo2mu2tau_Contin_MCFM701" || prod_mode=="ggTo4e_Contin_MCFM701" ||
          prod_mode=="ggTo4mu_Contin_MCFM701" || prod_mode=="ggTo4tau_Contin_MCFM701") {
    process = "ggZZ";
  }
  else if (prod_mode.Contains("H1")) {
    process = "signal"; //If "H125" is in the name of the prod_mode, it is a signal process
  }
  else {
    process = "AC";
  }
  // //Change prod_mode label for qqZZ 2018
  // if(prod_mode=="ZZTo4l" && year=="2018") {
  //   prod_mode = "ZZTo4lext";
  // }

  cout << process << endl;

  TString input_dir, full_path;
  if(process=="AC"){
    input_dir = "/eos/user/a/atarabin/MC_samples";
    full_path = Form("%s/AC%s_MELA/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());
    cout << full_path << endl;
  }
  else if((prod_mode.Contains("H125")) || (process.Contains("ZZ"))){ //We use the ones with the correct JES implementation for signal 125 and bkgs
    input_dir = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIIUL";
    full_path = Form("%s/MC%s_JES/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());
    cout << full_path << endl;
  }
  else{
    // input_dir = "/eos/user/a/atarabin";
    input_dir = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIIUL";
    full_path = Form("%s/MC%s/%s/ZZ4lAnalysis.root", input_dir.Data(), year.Data(), prod_mode.Data());
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
  oldtree->SetBranchStatus("p_m4l_SIG",1);
  oldtree->SetBranchStatus("p_m4l_BKG",1);
  oldtree->SetBranchStatus("costhetastar",1);
  oldtree->SetBranchStatus("helcosthetaZ1",1);
  oldtree->SetBranchStatus("helcosthetaZ2",1);
  oldtree->SetBranchStatus("helphi",1);
  oldtree->SetBranchStatus("phistarZ1",1);
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
  oldtree->SetBranchStatus("LepisCrack",1);
  oldtree->SetBranchStatus("LepSF",1);
  oldtree->SetBranchStatus("LepSF_Unc",1);
  oldtree->SetBranchStatus("ExtraLepPt",1);
  oldtree->SetBranchStatus("ExtraLepPhi",1);
  oldtree->SetBranchStatus("ExtraLepEta",1);
  oldtree->SetBranchStatus("ExtraLepLepId",1);
  // JESUp
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree->SetBranchStatus(Form("JetPt_JESUp_%s", jes_name.at(i).Data()),1);
  }
  // JESDown
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree->SetBranchStatus(Form("JetPt_JESDown_%s", jes_name.at(i).Data()),1);
  }
  // // JetSigma
  // for(unsigned int i=0;i<jes_name.size();i++){
  //   oldtree->SetBranchStatus(Form("JetSigma_%s", jes_name.at(i).Data()),1);
  //   // oldtree->SetBranchStatus(Form("nCleanedJetsPt30_jesUp_%s", jes_name.at(i).Data()),1);
  // }
  if(process=="qqZZ") {
    oldtree->SetBranchStatus("KFactor_EW_qqZZ",1);
    oldtree->SetBranchStatus("KFactor_QCD_qqZZ_M",1);
  }
  if(process=="ggZZ") {
    oldtree->SetBranchStatus("KFactor_QCD_ggZZ_Nominal",1);
  }
  oldtree->SetBranchStatus("PUWeight",1);
  oldtree->SetBranchStatus("genHEPMCweight",1);
  // if(year!="2016") oldtree->SetBranchStatus("genHEPMCweight_NNLO",1);
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
  oldtree->SetBranchStatus("p_QQB_BKG_MCFM",1);
  if(process=="signal" || process=="AC"){
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
  if(prod_mode.Contains("ggH")) oldtree->SetBranchStatus("ggH_NNLOPS_weight",1); // Additional entry for the weight in case of ggH

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
    // oldtree_failed->SetBranchStatus("ggH_NNLOPS_weight",1);
    oldtree_failed->SetBranchStatus("PUWeight",1);
    oldtree_failed->SetBranchStatus("genHEPMCweight",1);
    // if(year!="2016") oldtree_failed->SetBranchStatus("genHEPMCweight_NNLO",1);
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
    if(prod_mode.Contains("ggH")) oldtree_failed->SetBranchStatus("ggH_NNLOPS_weight",1); // Additional entry for the weight in case of ggH
  }

  // Copy branches in the new file
  if(process!="AC") {
    input_dir = "/eos/user/a/atarabin/MC_samples";
  }
  TString new_name = Form("%s_reducedTree_MC_%s.root", prod_mode.Data(), year.Data());
  TString new_full_path;
  if(process!="AC") {
    new_full_path = Form("%s/%sUL/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  }
  else {
    new_full_path = Form("%s/AC%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),new_name.Data());
  }
  cout << new_full_path << endl;
  TFile *newfile = new TFile(new_full_path.Data(),"RECREATE");
  TTree *newtree = (TTree*) oldtree->CloneTree(0);
  newtree->CopyEntries(oldtree);
  newtree->Write(0, TObject::kOverwrite); // Write candTree
  if(process=="signal" || process=="AC"){
    TTree *newtree_failed = (TTree*) oldtree_failed->CloneTree(0);
    newtree_failed->CopyEntries(oldtree_failed);
    newtree_failed->Write(0, TObject::kOverwrite); // Write candTree_failed
  }
  hCounters->Write(); // Write Counters
  newfile->Close();

  bool t_failed;
  if(process=="signal" || process=="AC") add(input_dir, year, prod_mode, process, t_failed = true);
  add(input_dir, year, prod_mode, process, t_failed = false);

  // if(process=="signal" || process=="AC"){
  //   // Merge together into a single TTree. Useful for efficiencies calculation.
  //   TFile* inputfile = TFile::Open(new_full_path.Data(), "READ");
  //   TTree* tree1 = (TTree*) inputfile->Get("candTree");
  //   TTree* tree2 = (TTree*) inputfile->Get("candTree_failed");
  //   TH1F* cnts = (TH1F*) inputfile->Get("Counters");
  //
  //   TString merged_name = Form("%s_mergedTree_MC_%s.root", prod_mode.Data(), year.Data());
  //   TString merged_path;
  //   if(process!="AC") merged_path = Form("%s/MC%s/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),merged_name.Data());
  //   else merged_path = Form("%s/AC%s_MELA/%s/%s", input_dir.Data(),year.Data(),prod_mode.Data(),merged_name.Data());
  //   cout << merged_path << endl;
  //
  //   TFile* mergedTTree = new TFile(merged_path.Data(), "RECREATE");
  //   TList* alist = new TList;
  //
  //   alist->Add(tree1);
  //   alist->Add(tree2);
  //
  //   TTree *newtree_single = TTree::MergeTrees(alist);
  //   newtree_single->SetName("fullTree");
  //   newtree_single->Write();
  //   cnts->Write();
  //   mergedTTree->Close();
  //   inputfile->Close();
  // }

  return 0;
}
