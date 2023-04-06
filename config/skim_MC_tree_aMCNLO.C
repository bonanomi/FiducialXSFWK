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


void skim_MC_tree_aMCNLO (TString prod_mode = "VBFH125", TString year = "2018"){

  TString nfile = Form("/eos/user/m/mbonanom/HIG-21-009/MC_samples/2016_MELA/ggH125_aMCNLO/ggH_amcatnloFXFX_slimmed.root");
  TFile *f = new TFile(nfile.Data(),"UPDATE");
  TTree *T;
  T = (TTree*)f->Get("Ana/passedEvents");
  std::cout << T->GetName() << std::endl;

  return 0;
}
