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
  TBranch *absdphijj = T->Branch("absdphijj",&_absdphijj,"absdphijj/F");
  TBranch *absdetajj = T->Branch("absdetajj",&_absdetajj,"absdetajj/F");
  TBranch *TCjmax = T->Branch("TCjmax",&_TCjmax,"TCjmax/F");
  TBranch *TBjmax = T->Branch("TBjmax",&_TBjmax,"TBjmax/F");
  TBranch *TCj1 = T->Branch("TCj1",&_TCj1,"TCj1/F");
  TBranch *TBj1 = T->Branch("TBj1",&_TBj1,"TBj1/F");
  // TBranch *pTj1_eta4p7 = T->Branch("pTj1_eta4p7",&_pTj1_eta4p7,"pTj1_eta4p7/F");
  TBranch *njets_pt30_eta4p7 = T->Branch("njets_pt30_eta4p7",&_njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
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

    // njets
    _njets_pt30_eta2p5 = 0;
    _njets_pt30_eta4p7 = 0;
    for(unsigned int i=0;i<JetPt->size();i++){
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<2.5){
        _njets_pt30_eta2p5 = _njets_pt30_eta2p5 + 1;
      }
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7){
        _njets_pt30_eta4p7++;
      }
    }

    // leading jet pT
    _pTj1 = 0;
    // _pTj1_eta4p7 = 0;
    Mj1 = 0;
    ETAj1 = 0;
    PHIj1 = 0;
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7 && JetPt->at(i) > _pTj1) {
        _pTj1 = JetPt->at(i);
        Mj1 = JetMass->at(i);
        ETAj1 = JetEta->at(i);
        PHIj1 = JetPhi->at(i);
      }
      // if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7 && JetPt->at(i) > _pTj1_eta4p7) {
      //   _pTj1_eta4p7 = JetPt->at(i);
      // }
    }
    // sub-leading jet pT
    _pTj2 = 0;
    Mj2 = 0;
    ETAj2 = 0;
    PHIj2 = 0;
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
      if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7 && JetPt->at(i) > _pTj2 && _pTj1 != JetPt->at(i)) {
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
    _TCj = -1; _TBj = -1;
    _TCj1 = 0; _TBj1 = 0;
    _TCjmax = -1; _TBjmax = -1;
    _TCj1 = sqrt(_pTj1*_pTj1 + Mj1*Mj1)/(2*cosh(j1.Rapidity() - H.Rapidity()));
    _TBj1 = sqrt(_pTj1*_pTj1 + Mj1*Mj1)*exp(-1*abs(j1.Rapidity() - H.Rapidity()));
    for (unsigned int i = 0; i < JetPt->size(); ++i)
    {
       if(JetPt->at(i)>30 && abs(JetEta->at(i))<4.7) {
          TLorentzVector theJet;
          theJet.SetPtEtaPhiM(JetPt->at(i), JetEta->at(i), JetPhi->at(i), JetMass->at(i));
          _TCj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))/(2*cosh(theJet.Rapidity() - H.Rapidity())); //theJet.E());
          _TBj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2))*exp(-1*abs(theJet.Rapidity() - H.Rapidity()));
          if (_TCj > _TCjmax) _TCjmax = _TCj;
          if (_TBj > _TBjmax) _TBjmax = _TBj;
       }
    }

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
      // _dphijj = j1.Phi() - j2.Phi();
      _dphijj = deltaphi(j1,j2);
      _absdphijj = abs(_dphijj);
    }
    else {
      _pTHjj = -1;
      _mHjj = -1;
      _detajj = -99;
      _absdetajj = -99;
      _mjj = -1;
      _dphijj = -99;
      _absdphijj = -1;
    }

    // _ZZy = abs(log((sqrt(125*125 + ZZPt*ZZPt*cosh(ZZEta)*cosh(ZZEta))+ZZPt*sinh(ZZEta))/sqrt(125*125+ZZPt*ZZPt)));
    _ZZy = abs(H.Rapidity());

    // MELA probabilities
    int ZZFlav = Z1Flav * Z2Flav;
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
    njets_pt30_eta2p5->Fill();
    pTj1->Fill();
    pTj2->Fill();
    pTHj->Fill();
    pTHjj->Fill();
    mHj->Fill();
    mHjj->Fill();
    detajj->Fill();
    absdetajj->Fill();
    dphijj->Fill();
    mjj->Fill();
    ZZy->Fill();
    chan->Fill();
    CMS_zz4l_mass->Fill();
    TCjmax->Fill();
    TBjmax->Fill();
    TCj1->Fill();
    TBj1->Fill();
    njets_pt30_eta4p7->Fill();
    // pTj1_eta4p7->Fill();
  }
  T->Print();
  T->Write("", TObject::kOverwrite);
  delete f;
}


void skim_data_tree_v2 (TString year = "2018"){

  TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIIUL";
  // auto oldFile = TFile::Open(Form("%s/Data_%s/AllData/ZZ4lAnalysis.root", path.Data(), year.Data())); //2016
  // auto oldFile = TFile::Open(Form("%s/data%s_JES/AllData/ZZ4lAnalysis.root", path.Data(), year.Data())); //2017
  auto oldFile = TFile::Open(Form("%s/data%s_JES/AllData_noZTree/ZZ4lAnalysis.root", path.Data(), year.Data())); //2018
  // TString path = "/eos/user/a/atarabin/2018_JES";
  // auto oldFile = TFile::Open(Form("%s/%s/ZZ4lAnalysis.root", path.Data(), year.Data()));
  // TString path = "/eos/user/a/atarabin/data_newProd/2016";
  // auto oldFile = TFile::Open(Form("%s/AllData/ZZ4lAnalysis.root", path.Data()));

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
  // oldtree->SetBranchStatus("JetPt_JESUp",1);
  // oldtree->SetBranchStatus("JetPt_JESDown",1);
  // JESUp
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree->SetBranchStatus(Form("JetPt_JESUp_%s", jes_name.at(i).Data()),1);
  }
  // JESDown
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree->SetBranchStatus(Form("JetPt_JESDown_%s", jes_name.at(i).Data()),1);
  }
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
  oldtree->SetBranchStatus("p_m4l_SIG",1);
  oldtree->SetBranchStatus("p_m4l_BKG",1);

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
  // oldtree_CR->SetBranchStatus("JetPt_JESUp",1);
  // oldtree_CR->SetBranchStatus("JetPt_JESDown",1);
  // JESUp
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree_CR->SetBranchStatus(Form("JetPt_JESUp_%s", jes_name.at(i).Data()),1);
  }
  // JESDown
  for(unsigned int i=0;i<jes_name.size();i++){
    oldtree_CR->SetBranchStatus(Form("JetPt_JESDown_%s", jes_name.at(i).Data()),1);
  }
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
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz4_1_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz2_1_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",1);
  oldtree_CR->SetBranchStatus("p_QQB_BKG_MCFM",1);
  oldtree_CR->SetBranchStatus("p_m4l_SIG",1);
  oldtree_CR->SetBranchStatus("p_m4l_BKG",1);


  TString newpath = "/eos/user/a/atarabin/Data";
  TFile *newfile = new TFile(Form("%s/reducedTree_AllData_%s.root", newpath.Data(), year.Data()),"RECREATE");

  // TString newpath = "/eos/user/a/atarabin/2018_JES";
  // TFile *newfile = new TFile(Form("%s/%s/reducedTree.root", newpath.Data(), year.Data()),"RECREATE");

  // TString newpath = "/eos/user/a/atarabin/data_newProd/2016";
  // TFile *newfile = new TFile(Form("%s/%s/reducedTree.root", newpath.Data(), year.Data()),"RECREATE");

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
