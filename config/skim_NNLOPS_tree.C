// root -l skim_NNLOPS_tree.C

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


float mass_lep(int flavour){
  if((abs(flavour)) == 11) return 0.0005109989461;
  else if ((abs(flavour)) == 13) return 0.1056583745;
  else if ((abs(flavour)) == 15) return 1.77686;
  // else if ((abs(flavour)) == 0) return 0;
  return 0;
}

//------------------------------------------------------------------
void add(){

  TFile *f = new TFile("ggH_NNLOPS.root","UPDATE");
  TTree *T;
  T = (TTree*)f->Get("passedEvents");

  bool passedFiducialSelection;
  Float_t GENrapidity4l, GENpT4l, GENeta4l, GENmass4l;

  vector<float> *GENlep_pt = 0;
  vector<float> *GENlep_eta = 0;
  vector<float> *GENlep_phi = 0;
  vector<float> *GENlep_id = 0;

  vector<float> *GENjet_pt = 0;
  vector<float> *GENjet_eta = 0;
  vector<float> *GENjet_phi = 0;
  vector<float> *GENjet_mass = 0;

  T->SetBranchAddress("GENlep_pt",&GENlep_pt);
  T->SetBranchAddress("GENlep_eta",&GENlep_eta);
  T->SetBranchAddress("GENlep_phi",&GENlep_phi);
  T->SetBranchAddress("GENlep_id",&GENlep_id);
  T->SetBranchAddress("GENjet_pt",&GENjet_pt);
  T->SetBranchAddress("GENjet_eta",&GENjet_eta);
  T->SetBranchAddress("GENjet_phi",&GENjet_phi);
  T->SetBranchAddress("GENjet_mass",&GENjet_mass);
  T->SetBranchAddress("passedFiducialSelection",&passedFiducialSelection);
  T->SetBranchAddress("GENrapidity4l",&GENrapidity4l);
  T->SetBranchAddress("GENpT4l",&GENpT4l);
  T->SetBranchAddress("GENeta4l",&GENeta4l);
  T->SetBranchAddress("GENmass4l",&GENmass4l);


  Float_t _GENpTHj,_GENpTHjj,_GENmHj,_GENmHjj,_GENdetajj,_GENabsdetajj,_GENmjj,_GENdphijj,_GENabsdphijj,_GENTCjmax,_GENTBjmax,
          _GENTCj,_GENTBj,_GENTCj1,_GENTBj1,_GENpTj1,_GENpTj2,_GENrapidity4lAbs;
  Float_t GENMj1, GENETAj1, GENPHIj1, GENMj2, GENETAj2, GENPHIj2, GENphi4l;

  TBranch *GENpTj1 = T->Branch("GENpTj1",&_GENpTj1,"GENpTj1/F");
  TBranch *GENpTj2 = T->Branch("GENpTj2",&_GENpTj2,"GENpTj2/F");
  TBranch *GENpTHj = T->Branch("GENpTHj",&_GENpTHj,"GENpTHj/F");
  TBranch *GENpTHjj = T->Branch("GENpTHjj",&_GENpTHjj,"GENpTHjj/F");
  TBranch *GENmHj = T->Branch("GENmHj",&_GENmHj,"GENmHj/F");
  TBranch *GENmHjj = T->Branch("GENmHjj",&_GENmHjj,"GENmHjj/F");
  TBranch *GENdetajj = T->Branch("GENdetajj",&_GENdetajj,"GENdetajj/F");
  TBranch *GENabsdetajj = T->Branch("GENabsdetajj",&_GENabsdetajj,"GENabsdetajj/F");
  TBranch *GENmjj = T->Branch("GENmjj",&_GENmjj,"GENmjj/F");
  TBranch *GENdphijj = T->Branch("GENdphijj",&_GENdphijj,"GENdphijj/F");
  TBranch *GENabsdphijj = T->Branch("GENabsdphijj",&_GENabsdphijj,"GENabsdphijj/F");
  TBranch *GENTCjmax = T->Branch("GENTCjmax",&_GENTCjmax,"GENTCjmax/F");
  TBranch *GENTBjmax = T->Branch("GENTBjmax",&_GENTBjmax,"GENTBjmax/F");
  TBranch *GENrapidity4lAbs = T->Branch("GENrapidity4lAbs",&_GENrapidity4lAbs,"GENrapidity4lAbs/F");


  Long64_t nentries = T->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);

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

    if(!passedFiducialSelection) continue;

    _GENrapidity4lAbs = abs(GENrapidity4l);


    // leading GENjet pT
    _GENpTj1 = -1;
    GENMj1 = 0;
    GENETAj1 = 0;
    GENPHIj1 = 0;
    for (unsigned int i = 0; i < GENjet_pt->size(); ++i){
      if(GENjet_pt->at(i) > _GENpTj1) {
        _GENpTj1 = GENjet_pt->at(i);
        GENMj1 = GENjet_mass->at(i);
        GENETAj1 = GENjet_eta->at(i);
        GENPHIj1 = GENjet_phi->at(i);
      }
    }


    // sub-leading GENjet pT
    _GENpTj2 = -1;
    GENMj2 = 0;
    GENETAj2 = 0;
    GENPHIj2 = 0;
    for (unsigned int i = 0; i < GENjet_pt->size(); ++i){
      if(GENjet_pt->at(i) > _GENpTj2 && _GENpTj1 != GENjet_pt->at(i)) {
        _GENpTj2 = GENjet_pt->at(i);
        GENMj2 = GENjet_mass->at(i);
        GENETAj2 = GENjet_eta->at(i);
        GENPHIj2 = GENjet_phi->at(i);
      }
    }

    // Build GENphi4l (not saved in original ROOT files)
    TLorentzVector lep_1;
    lep_1.SetPtEtaPhiM(GENlep_pt->at(0), GENlep_eta->at(0), GENlep_phi->at(0), mass_lep(GENlep_id->at(0)));
    TLorentzVector lep_2;
    lep_2.SetPtEtaPhiM(GENlep_pt->at(1), GENlep_eta->at(1), GENlep_phi->at(1), mass_lep(GENlep_id->at(1)));
    TLorentzVector lep_3;
    lep_3.SetPtEtaPhiM(GENlep_pt->at(2), GENlep_eta->at(2), GENlep_phi->at(2), mass_lep(GENlep_id->at(2)));
    TLorentzVector lep_4;
    lep_4.SetPtEtaPhiM(GENlep_pt->at(3), GENlep_eta->at(3), GENlep_phi->at(3), mass_lep(GENlep_id->at(3)));
    GENphi4l = (lep_1+lep_2+lep_3+lep_4).Phi();

    TLorentzVector GENH;
    TLorentzVector GENj1;
    TLorentzVector GENj2;
    GENH.SetPtEtaPhiM(GENpT4l,GENeta4l,GENphi4l,GENmass4l);
    GENj1.SetPtEtaPhiM(_GENpTj1,GENETAj1,GENPHIj1,GENMj1);
    GENj2.SetPtEtaPhiM(_GENpTj2,GENETAj2,GENPHIj2,GENMj2);
    _GENTCj1 = sqrt(_GENpTj1*_GENpTj1 + GENMj1*GENMj1)/(2*cosh(GENj1.Rapidity() - GENH.Rapidity()));
    _GENTBj1 = sqrt(_GENpTj1*_GENpTj1 + GENMj1*GENMj1)*exp(-1*abs(GENj1.Rapidity() - GENH.Rapidity()));

    _GENTCjmax = -1; _GENTBjmax = -1;
    for (unsigned int i = 0; i < GENjet_pt->size(); ++i){
      TLorentzVector theJet;
      theJet.SetPtEtaPhiM(GENjet_pt->at(i), GENjet_eta->at(i), GENjet_phi->at(i), GENjet_mass->at(i));
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
      _GENdphijj = deltaphi(GENj1,GENj2);
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

    GENpTj1->Fill();
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
    GENrapidity4lAbs->Fill();
  }

  T->Write(0, TObject::kOverwrite);
  delete f;
  return;

}

//---------------------------------------------------------- MAIN ----------------------------------------------------------
void skim_NNLOPS_tree (){

  auto oldFile = TFile::Open("testGGH_nnlops_GENonly.root");
  TTree *oldtree = (TTree*) oldFile->Get("Ana/passedEvents");

  // Deactivate all branches
  oldtree->SetBranchStatus("*",0);
  // Activate some branches only: our skim
  oldtree->SetBranchStatus("genWeight",1);
  oldtree->SetBranchStatus("qcdWeights",1);
  oldtree->SetBranchStatus("nnloWeights",1);
  oldtree->SetBranchStatus("pdfENVup",1);
  oldtree->SetBranchStatus("pdfENVdown",1);
  oldtree->SetBranchStatus("passedFiducialSelection",1);
  oldtree->SetBranchStatus("GENlep_pt",1);
  oldtree->SetBranchStatus("GENlep_eta",1);
  oldtree->SetBranchStatus("GENlep_phi",1);
  oldtree->SetBranchStatus("GENlep_id",1);
  oldtree->SetBranchStatus("GENlep_Hindex",1);
  oldtree->SetBranchStatus("GENpT4l",1);
  oldtree->SetBranchStatus("GENeta4l",1);
  oldtree->SetBranchStatus("GENrapidity4l",1);
  oldtree->SetBranchStatus("GENcosTheta1",1);
  oldtree->SetBranchStatus("GENcosTheta2",1);
  oldtree->SetBranchStatus("GENcosThetaStar",1);
  oldtree->SetBranchStatus("GENPhi",1);
  oldtree->SetBranchStatus("GENPhi1",1);
  oldtree->SetBranchStatus("GENZ_DaughtersId",1);
  oldtree->SetBranchStatus("GENmass4l",1);
  oldtree->SetBranchStatus("GENmassZ1",1);
  oldtree->SetBranchStatus("GENmassZ2",1);
  oldtree->SetBranchStatus("GENjet_pt",1);
  oldtree->SetBranchStatus("GENjet_eta",1);
  oldtree->SetBranchStatus("GENjet_phi",1);
  oldtree->SetBranchStatus("GENjet_mass",1);
  oldtree->SetBranchStatus("GENnjets_pt30_eta4p7",1);
  oldtree->SetBranchStatus("GENpt_leadingjet_pt30_eta4p7",1);

  // Copy branches in the new file
  TFile *newfile = new TFile("ggH_NNLOPS.root","RECREATE");
  TTree *newtree = (TTree*) oldtree->CloneTree(0);
  newtree->CopyEntries(oldtree);
  newtree->Write(0, TObject::kOverwrite); // Write candTree
  newfile->Close();

  add();

  return 0;
}
