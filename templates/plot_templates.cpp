// c++ -o  plot_templates plot_templates.cpp `root-config --cflags --glibs`

// ROOT include
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include "TH2.h"
#include "TChain.h"
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TCanvas.h>
#include "TRandom.h"

// C include
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include "TLatex.h"

#include "setTDRStyle.C"
using namespace std;

const TString sPlotsStore = "plots";


void analysisInit() {
    gErrorIgnoreLevel = kWarning;
    gErrorIgnoreLevel = kError;
}

// ---------------------------------------------------------------------------------------------------------------------------------

void setCavasAndStyles(TString canvasName, TCanvas* &c, TString stat = "", double leftMaring = 0.15, double rightMaring = 0.05, double bottomMaring = 0.15, double topMaring = 0.05){
    // setup environment
    //gROOT->ProcessLine(".L setTDRStyle.C");
    setTDRStyle();
    // setup canvas
    c = new TCanvas(canvasName,"myPlots",0,0,800,800);
    c->cd(1); c->SetLogy(0);
    gStyle->SetOptStat(stat);
    gStyle->SetPalette(1);

    c->GetPad(0)->SetRightMargin(rightMaring);
    c->GetPad(0)->SetLeftMargin(leftMaring);
    c->GetPad(0)->SetTopMargin(topMaring);
    c->GetPad(0)->SetBottomMargin(bottomMaring);
}

// ---------------------------------------------------------------------------------------------------------------------------------

void cmsPreliminary(TCanvas* &c, TString top){
    c->cd();

    TLatex *CMSPrelim = new TLatex();
    CMSPrelim->SetNDC(kTRUE);

    CMSPrelim->SetTextSize(0.8*c->GetTopMargin());
    CMSPrelim->SetTextFont(42);
    CMSPrelim->SetTextAlign(31); // align right
    CMSPrelim->DrawLatex(0.93, 0.96,"#bf{"+top+"}");
}

// ---------------------------------------------------------------------------------------------------------------------------------

int setHistProperties(TH1D* &hist, Width_t lineWidth, Style_t lineStyle, Color_t lineColor, Style_t fillStyle=0, Color_t fillColor=0, TString xAxisTitle = "skip", TString yAxisTitle = "skip"){
    if (!hist) return -1;
    // line
    hist->SetLineWidth(lineWidth);
    hist->SetLineStyle(lineStyle);
    hist->SetLineColor(lineColor);
    // fill
    hist->SetFillStyle(fillStyle);
    hist->SetFillColor(fillColor);
    // divisions, offsets, sizes
    hist->GetXaxis()->SetNdivisions(510);
    hist->GetYaxis()->SetNdivisions(510);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleOffset(1.24);
    // titles
    if (xAxisTitle!="skip") hist->GetXaxis()->SetTitle(xAxisTitle);
    if (yAxisTitle!="skip") hist->GetYaxis()->SetTitle(yAxisTitle);
    // return
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------------------------

int setLegendProperties(TLegend* &leg, TString sHeader = "skip", Style_t fillStyle=0, Color_t fillColor=0){
    // sanity-check
    if (!leg) return -1;
    // titles
    if (sHeader!="skip") leg->SetHeader(sHeader);;
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    // return
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------------------------

int main (int argc, char ** argv){

  analysisInit();

  // setup environment & canvas
  TCanvas *c1;
  setCavasAndStyles("c1",c1,"");

  TString obsTag = argv[1];
  const int N_BINS = argc-3;
  TString binRange[N_BINS] = {};
  TString binRangeLow[N_BINS] = {};
  TString binRangeHigh[N_BINS] = {};
  TString binRangeLeg[N_BINS] = {};
  for (int i=0; i<N_BINS; i++){
    cout << i << endl;
    TString tmp = argv[i+2];
    TString tmpBis = argv[i+3];
    cout << tmp << endl;
    cout << tmpBis << endl;
    binRange[i] = tmp + TString("_") + tmpBis;
    binRangeLow[i] = tmp;
    binRangeHigh[i] = tmpBis;
    binRangeLeg[i] = tmp + TString(" < ") + obsTag + TString(" < ") + tmpBis;
  }

  // TString binRange[N_BINS]     = {"0_10", "10_20", "20_30", "30_45", "45_80", "80_120", "120_200", "200_1300"};
  // TString binRangeLow[N_BINS]  = {"0", "10", "20", "30", "45", "80", "120", "200"};
  // TString binRangeHigh[N_BINS] = {"10", "20", "30", "45", "80", "120", "200", "1300"};
  // TString binRangeLeg[N_BINS]  = {"0 < pT < 10 GeV", "10 < pT < 15 GeV", "20 < pT < 30 GeV", "30 < pT < 45 GeV", "45 < pT < 80 GeV", "80 < pT < 120 GeV", "120 < pT < 200 GeV", "200 < pT < 1300 GeV"};


  // TString obsTag = "pT4l";
  // const int N_BINS = 8;
  // TString binRange[N_BINS]     = {"0_10", "10_20", "20_30", "30_45", "45_80", "80_120", "120_200", "200_1300"};
  // TString binRangeLow[N_BINS]  = {"0", "10", "20", "30", "45", "80", "120", "200"};
  // TString binRangeHigh[N_BINS] = {"10", "20", "30", "45", "80", "120", "200", "1300"};
  // TString binRangeLeg[N_BINS]  = {"0 < pT < 10 GeV", "10 < pT < 15 GeV", "20 < pT < 30 GeV", "30 < pT < 45 GeV", "45 < pT < 80 GeV", "80 < pT < 120 GeV", "120 < pT < 200 GeV", "200 < pT < 1300 GeV"};

  // TString obsTag = "D0m";
  // const int N_BINS = 5;
  // TString binRange[N_BINS]     = {"0.0_0.4", "0.4_0.55", "0.55_0.7", "0.7_0.85", "0.85_1.0"};
  // TString binRangeLow[N_BINS]  = {"0.0", "0.4", "0.55", "0.7", "0.85"};
  // TString binRangeHigh[N_BINS] = {"0.4", "0.55", "0.7", "0.85", "1.0"};
  // TString binRangeLeg[N_BINS]  = {"0 < D0m < 0.4 ", "0.4 < D0m < 0.55", "0.55 < pT < 0.7", "0.7 < pT < 0.85", "0.85 < pT < 1"};


  // TString obsTag = "rapidity4l";
  // const int N_BINS = 6;
  // TString binRange[N_BINS]     = {"0.0_0.15", "0.15_0.3", "0.3_0.6", "0.6_0.9", "0.9_1.2", "1.2_2.5"};
  // TString binRangeLow[N_BINS]  = {"0.0", "0.15", "0.3", "0.6", "0.9", "1.2"};
  // TString binRangeHigh[N_BINS] = {"0.15", "0.3", "0.6", "0.9", "1.2", "2.5"};
  // TString binRangeLeg[N_BINS]  = {"0 <  |y(H)| < 0.15", "0.15 < |y(H)| < 0.3", "0.3 < |y(H)| < 0.6 GeV", "0.6 < |y(H)| < 0.9 GeV", "0.9 < |y(H)| < 1.2 GeV", "1.2 < |y(H)| < 2.5 GeV"};

  // TString obsTag = "massZ1";
  // const int N_BINS = 6;
  // TString binRange[N_BINS]     = {"40_65", "65_73", "73_80", "80_85", "85_90", "90_120"};
  // TString binRangeLow[N_BINS]  = {"40", "65", "73", "80", "85", "90"};
  // TString binRangeHigh[N_BINS] = {"65", "73", "80", "85", "90", "120"};
  // TString binRangeLeg[N_BINS]  = {"40 <  m(Z1) < 65 GeV", "65 < m(Z1) < 73 GeV", "73 < m(Z1) < 80 GeV", "80 < m(Z1) < 85 GeV", "85 < m(Z1) < 90 GeV", "90 < m(Z1) < 120 GeV"};
  //
  // TString obsTag = "massZ2";
  // const int N_BINS = 6;
  // TString binRange[N_BINS]     = {"12_24", "24_28", "28_32", "32_40", "40_48", "48_65"};
  // TString binRangeLow[N_BINS]  = {"12", "24", "28", "32", "40", "48"};
  // TString binRangeHigh[N_BINS] = {"24", "28", "32", "40", "48", "65"};
  // TString binRangeLeg[N_BINS]  = {"12 <  m(Z2) < 24 GeV", "24 < m(Z2) < 28 GeV", "28 < m(Z2) < 32 GeV", "32 < m(Z2) < 40 GeV", "40 < m(Z2) < 48 GeV", "48 < m(Z2) < 65 GeV"};

  // TString obsTag = "njets_pt30_eta2p5";
  // const int N_BINS = 5;
  // TString binRange[N_BINS]     = {"0_1", "1_2", "2_3", "3_4", "4_20"};
  // TString binRangeLow[N_BINS]  = {"0", "1", "2", "3", "4"};
  // TString binRangeHigh[N_BINS] = {"1", "2", "3", "4", "20"};
  // TString binRangeLeg[N_BINS]  = {"nJet=0", "nJet=1", "nJet=2", "nJet=3", "nJet>=4"};

  // TString obsTag = "pTj1";
  // const int N_BINS = 4;
  // TString binRange[N_BINS]     = {"30_55", "55_95", "95_200", "200_13000"};
  // TString binRangeLow[N_BINS]  = {"30", "55", "95", "200"};
  // TString binRangeHigh[N_BINS] = {"55", "95", "200", "13000"};
  // TString binRangeLeg[N_BINS]  = {"30 < pTj1 < 55 GeV", "55 < pTj1 < 95 GeV", "95 < pTj1 < 200 GeV", "200 < pTj1 < 1300 GeV"};

//   TString obsTag = "costhetastar";
//   const int N_BINS = 5;
//   TString binRange[N_BINS]     = {"0.0_0.2", "0.2_0.4", "0.4_0.6", "0.6_0.8", "0.8_1.0"};
//   TString binRangeLow[N_BINS]  = {"0.0", "0.2", "0.4", "0.6", "0.8"};
//   TString binRangeHigh[N_BINS] = {"0.2", "0.4", "0.6", "0.8", "1.0"};
//   TString binRangeLeg[N_BINS]  = {"0.0 < cos(#theta) < 0.2 ", "55 < cos(#theta) < 0.4 ", "0.4 < cos(#theta) < 0.6 ", "0.6 < cos(#theta) < 0.8 ", "0.8 < cos(#theta) < 1.0 "};

  // TString obsTag = "massZ1_massZ2";
  // const int N_BINS = 5;
  // TString binRange[N_BINS]     = {"40_82_12_32", "40_74_32_65", "74_120_32_65", "82_120_24_32", "82_120_12_24"};
  // TString binRangeLow[N_BINS]  = {"40_12", "40_32", "74_32", "82_24", "82_12"};
  // TString binRangeHigh[N_BINS] = {"82_32", "74_65", "120_65", "120_32", "120_24"};
  // TString binRangeLeg[N_BINS]  = {"40 < m(Z1) < 82 / 12 < m(Z2) < 32 GeV", "40 < m(Z1) < 74 / 32 < m(Z2) 65 GeV", "74 < m(Z1) < 120 / 32 < m(Z2) 65 GeV", "82 < m(Z1) < 120 / 24 < m(Z2) 32 GeV", "82 < m(Z1) < 120 / 12 < m(Z2) 24 GeV"};


  const int N_BKGS = 3;
  TString bkgName[N_BKGS]   = {"qqzz", "ggzz", "ZJetsCR"};

  const int N_YEAR = 3;
  TString year[N_YEAR]      = {"2016", "2018", "2017"}; //, "2017", "2018"};

  for (int iYear = 0; iYear<N_YEAR; iYear++){
    const TString sTemplateDirName = year[iYear]+"/"+obsTag;

    TFile* fTemplateFile_2e2mu[N_BKGS][N_BINS];
    TFile* fTemplateFile_4mu[N_BKGS][N_BINS];
    TFile* fTemplateFile_4e[N_BKGS][N_BINS];
    TH1D* h1D_2e2mu[N_BKGS][N_BINS];
    TH1D* h1D_4mu[N_BKGS][N_BINS];
    TH1D* h1D_4e[N_BKGS][N_BINS];

    cout << "obsTag: " << obsTag << endl;
    for (int iBin = 0; iBin<N_BINS; iBin++ ) {
        for (int iBkg = 0; iBkg<N_BKGS; iBkg++ ) {
          // int nSmooth = 2;

          TString sTemplateFileName = "XSBackground_"+bkgName[iBkg]+"_2e2mu_"+obsTag+"_"+binRange[iBin]+".root";
          fTemplateFile_2e2mu[iBkg][iBin] = new TFile(sTemplateDirName+"/"+sTemplateFileName, "READ");
          h1D_2e2mu[iBkg][iBin] = (TH1D*) fTemplateFile_2e2mu[iBkg][iBin]->Get("m4l_"+obsTag+"_"+binRange[iBin]);
          // cout << "sTemplateDirName/sTemplateFileName: " << sTemplateDirName+"/"+sTemplateFileName << endl;
          // cout << "h1D_2e2mu["<<bkgName[iBkg]<<"]["<<binRange[iBin]<<"]->GetEntries(): " << h1D_2e2mu[iBkg][iBin]->GetEntries() << endl;
          // for (int k = 0; k < nSmooth; k++) smoothAndNormaliseTemplate1D(h1D_2e2mu[iBkg][iBin]);

          sTemplateFileName = "XSBackground_"+bkgName[iBkg]+"_4mu_"+obsTag+"_"+binRange[iBin]+".root";
          fTemplateFile_4mu[iBkg][iBin] = new TFile(sTemplateDirName+"/"+sTemplateFileName, "READ");
          h1D_4mu[iBkg][iBin] = (TH1D*) fTemplateFile_4mu[iBkg][iBin]->Get("m4l_"+obsTag+"_"+binRange[iBin]);
          // for (int k = 0; k < nSmooth; k++) smoothAndNormaliseTemplate1D(h1D_4mu[iBkg][iBin]);

          sTemplateFileName = "XSBackground_"+bkgName[iBkg]+"_4e_"+obsTag+"_"+binRange[iBin]+".root";
          fTemplateFile_4e[iBkg][iBin] = new TFile(sTemplateDirName+"/"+sTemplateFileName, "READ");
          h1D_4e[iBkg][iBin] = (TH1D*) fTemplateFile_4e[iBkg][iBin]->Get("m4l_"+obsTag+"_"+binRange[iBin]);
          // for (int k = 0; k < nSmooth; k++) smoothAndNormaliseTemplate1D(h1D_4e[iBkg][iBin]);
        } // iBkg<N_BKGS
      } // iBin<N_BINS

      // prepare dummy
      double var_plotHigh = 140, var_plotLow = 105;
      int var_nBins = 20;
      TString varAxLabel = "m_{4l} (GeV)";
      double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
      TString sUnit = (varAxLabel.Contains(" (GeV)"))?"(GeV)":" ";
      TString sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
      TH1D* h1D_dummy = new TH1D("dummy", "dummy", var_nBins, var_plotLow, var_plotHigh);
      setHistProperties(h1D_dummy,1,1,kBlue-7,0,0,varAxLabel,"Events/"+sBinWidth);

      // common proeprties
      Width_t lineWidth = 2.5;
      double leg_xl = 0.52, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;

      // plot hists
      int kBkg_qqZZ = 0, kBkg_ggZZ = 1, kBkg_ZJets = 2;
      c1->cd();
      for (int iBin = 0; iBin<N_BINS; iBin++ ) {
          /////// 2e2mu /////
          // qqZZZ + ggZZ +ZX //
          h1D_dummy->SetMaximum(2.0*h1D_2e2mu[kBkg_qqZZ][iBin]->GetMaximum());
          h1D_dummy->Draw(); cmsPreliminary(c1, binRangeLeg[iBin]+"      2e2#mu      "+year[iYear]); TLegend* leg1 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg1);
          setHistProperties(h1D_2e2mu[kBkg_qqZZ][iBin],lineWidth,1,kBlack); h1D_2e2mu[kBkg_qqZZ][iBin]->Draw("histsame"); leg1->AddEntry(h1D_2e2mu[kBkg_qqZZ][iBin], "q#bar{q} #rightarrow ZZ","L");
          setHistProperties(h1D_2e2mu[kBkg_ggZZ][iBin],lineWidth,1,kBlue-7); h1D_2e2mu[kBkg_ggZZ][iBin]->Draw("histsame"); leg1->AddEntry(h1D_2e2mu[kBkg_ggZZ][iBin], "gg #rightarrow ZZ","L");
          setHistProperties(h1D_2e2mu[kBkg_ZJets][iBin],lineWidth,1,kRed-7); h1D_2e2mu[kBkg_ZJets][iBin]->Draw("histsame"); leg1->AddEntry(h1D_2e2mu[kBkg_ZJets][iBin], "Z + X","L");
          leg1->Draw();
          c1->SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_2e2mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".pdf");
          c1->SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_2e2mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".png");


          /////// 4mu /////
          // qqZZZ + ggZZ +ZX //
          h1D_dummy->SetMaximum(2.0*h1D_4mu[kBkg_qqZZ][iBin]->GetMaximum());
          h1D_dummy->Draw(); cmsPreliminary(c1, binRangeLeg[iBin]+"      4#mu      "+year[iYear]); TLegend* leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2);
          setHistProperties(h1D_4mu[kBkg_qqZZ][iBin],lineWidth,1,kBlack); h1D_4mu[kBkg_qqZZ][iBin]->Draw("histsame"); leg2->AddEntry(h1D_4mu[kBkg_qqZZ][iBin], "q#bar{q} #rightarrow ZZ","L");
          setHistProperties(h1D_4mu[kBkg_ggZZ][iBin],lineWidth,1,kBlue-7); h1D_4mu[kBkg_ggZZ][iBin]->Draw("histsame"); leg2->AddEntry(h1D_4mu[kBkg_ggZZ][iBin], "gg #rightarrow ZZ","L");
          setHistProperties(h1D_4mu[kBkg_ZJets][iBin],lineWidth,1,kRed-7); h1D_4mu[kBkg_ZJets][iBin]->Draw("histsame"); leg2->AddEntry(h1D_4mu[kBkg_ZJets][iBin], "Z + X","L");
          leg2->Draw();
          c1->SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".pdf");
          c1->SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4mu_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".png");

          /////// 4e /////
          // qqZZZ + ggZZ +ZX //
          h1D_dummy->SetMaximum(2.0*h1D_4e[kBkg_qqZZ][iBin]->GetMaximum());
          h1D_dummy->Draw(); cmsPreliminary(c1, binRangeLeg[iBin]+"      4e      "+year[iYear]); TLegend* leg3 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg3);
          setHistProperties(h1D_4e[kBkg_qqZZ][iBin],lineWidth,1,kBlack); h1D_4e[kBkg_qqZZ][iBin]->Draw("histsame"); leg3->AddEntry(h1D_4e[kBkg_qqZZ][iBin], "q#bar{q} #rightarrow ZZ","L");
          setHistProperties(h1D_4e[kBkg_ggZZ][iBin],lineWidth,1,kBlue-7); h1D_4e[kBkg_ggZZ][iBin]->Draw("histsame"); leg3->AddEntry(h1D_4e[kBkg_ggZZ][iBin], "gg #rightarrow ZZ","L");
          setHistProperties(h1D_4e[kBkg_ZJets][iBin],lineWidth,1,kRed-7); h1D_4e[kBkg_ZJets][iBin]->Draw("histsame"); leg3->AddEntry(h1D_4e[kBkg_ZJets][iBin], "Z + X","L");
          leg3->Draw();
          c1->SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4e_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".pdf");
          c1->SaveAs(sPlotsStore+"/"+year[iYear]+"/"+obsTag+"/XSTemplates_4e_"+obsTag+"_"+year[iYear]+"_"+binRange[iBin]+"_"+bkgName[kBkg_qqZZ]+"_"+bkgName[kBkg_ggZZ]+"_"+bkgName[kBkg_ZJets]+".png");
      } // iBin<N_BINS

    }
  return 0;
}
