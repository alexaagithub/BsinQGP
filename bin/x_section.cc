/////////////////////////////////////////////////////////////////////////
//
// X-section, ratio and systematics
//
//
// August 2019
//
/////////////////////////////////////////////////////////////////////////


#include <sstream>
#include <vector>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBifurGauss.h>
#include "TMath.h"
#include <RooGenericPdf.h>
//#include "TRatioPlot.h"
#include <RooBifurGauss.h>
#include <RooProduct.h>
#include <RooHist.h>
#include "RooStats/SPlot.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include <iostream>
#include <TF1.h>
#include <RooPolynomial.h>
#include <fstream>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include "TMultiGraph.h"
#include <TAxis.h>

using namespace RooStats;
using namespace RooFit;
using namespace std;

void plot_xsection(int bin_n,double* pt_m, double* pt_l, double* pt_h, double* x_sec, double* stat, double* syst, TString meson);
double* error_syst_final(int n_bins,int n_sources_error,double (*source_error_bin)[4]);
double get_error(int bin, TGraphAsymmErrors * graph);
double get_pt_low(int bin, TGraphAsymmErrors * graph);
double get_pt_high(int bin, TGraphAsymmErrors * graph);


int main(){

  //Raw yield files
  //TFile* f_raw_yield_Bs = new TFile("./results/Bs/Bpt/pT.root");
  TFile* f_raw_yield_Bs = new TFile("/lstore/cms/ev19u032/pT_Bs_meson.root"); //file updated
  TFile* f_raw_yield_Bu = new TFile("/lstore/cms/ev19u032/pT_Bu_meson.root");
  
  //Efficiency files
  TFile* f_efficiency_Bs = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/for_students/MCstudiesPbPbPtBin.root");
  TFile* f_efficiency_Bu = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/for_students/MCstudiesPbPb_Bsbin.root");

  //Efficiency systematic error files
   TFile* f_eff_syst_Bs = new TFile("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/results/Bs/efficiency/root_files_Bpt/efficiency_systematic_errors.root");
  TFile* f_eff_syst_Bu = new TFile("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/results/Bu/efficiency/root_files_Bpt/efficiency_systematic_errors.root");

  const double branching_fraction_Bs = 0.0000313;
  const double branching_fraction_Bu = 0.0000599;

  const double branching_fraction_error_Bs = 0.0000030;
  const double branching_fraction_error_Bu = 0.0000023;

  const double luminosity = 1.5e-9;
  const double luminosity_error = 0.02;  

  double pt_bins[] = {5,10, 15, 20, 50};
  const int n_pt_bins = 4;
  const int n_s_errors = 4;
  //number of sources of systematic error
  //signal_yield, branching fraction, luminosity, efficiency+acceptance
  double s_errors_bin_Bu [n_pt_bins][n_s_errors];
  double s_errors_bin_Bs [n_pt_bins][n_s_errors];
  //arrays with the values of systematic errors for each bin

  //raw_yield

  TGraphAsymmErrors* raw_yield_Bs = (TGraphAsymmErrors*)f_raw_yield_Bs->Get("Graph");
  TGraphAsymmErrors* raw_yield_Bu = (TGraphAsymmErrors*)f_raw_yield_Bu->Get("Graph");

  //systematic errors 
  TGraphAsymmErrors* raw_yield_Bs_syst = (TGraphAsymmErrors*)f_raw_yield_Bs->Get("Graph;2");
  TGraphAsymmErrors* raw_yield_Bu_syst = (TGraphAsymmErrors*)f_raw_yield_Bu->Get("Graph;2");
  //statistical errors
  TGraphAsymmErrors* raw_yield_Bs_stat = (TGraphAsymmErrors*)f_raw_yield_Bs->Get("Graph;1");
  TGraphAsymmErrors* raw_yield_Bu_stat = (TGraphAsymmErrors*)f_raw_yield_Bu->Get("Graph;1");

  //efficiency
  TH1D* efficiency_Bs = new TH1D("efficiency_Bs", "efficiency_Bs", n_pt_bins, pt_bins);
  efficiency_Bs = (TH1D*)f_efficiency_Bs->Get("hEff");
  TH1D* efficiency_Bu = new TH1D("efficiency_Bu", "efficiency_Bu", n_pt_bins, pt_bins);
  efficiency_Bu = (TH1D*)f_efficiency_Bu->Get("hEff");

  TGraphErrors* eff_syst_Bs = (TGraphErrors*)f_eff_syst_Bs->Get("Graph");
  TGraphErrors* eff_syst_Bu = (TGraphErrors*)f_eff_syst_Bu->Get("Graph;1");


  //yield values
  double n;
  double* raw_Bs_y = raw_yield_Bs->GetY();
  double* raw_Bu_y = raw_yield_Bu->GetY();
  //pT values
  double* pT_Bs = raw_yield_Bs->GetX();
  double* pT_Bu = raw_yield_Bu->GetX();

  //yield and pT errors
  double syst_Bs[n_pt_bins];
  double syst_Bu[n_pt_bins];
  double stat_Bs[n_pt_bins]; 
  double stat_Bu[n_pt_bins];
  
  double pT_min_Bs[n_pt_bins];
  double pT_max_Bs[n_pt_bins];
  double pT_min_Bu[n_pt_bins];
  double pT_max_Bu[n_pt_bins];

  for(int i = 0; i < n_pt_bins; i++){
    syst_Bs[i] = get_error(i,raw_yield_Bs_syst);
    syst_Bu[i] = get_error(i,raw_yield_Bu_syst);
    stat_Bs[i] = get_error(i,raw_yield_Bs_stat);
    stat_Bu[i] = get_error(i,raw_yield_Bu_stat);
    
    pT_min_Bs[i] = get_pt_low(i, raw_yield_Bs_stat);
    pT_max_Bs[i] = get_pt_high(i, raw_yield_Bs_stat);
    pT_min_Bu[i] = get_pt_low(i, raw_yield_Bu_stat);
    pT_max_Bu[i] = get_pt_high(i, raw_yield_Bu_stat);

  } 


  //efficiency value
  double eff;
  double* eff_s_Bs = eff_syst_Bs->GetY();
  cout<<"eff_Bs: "<<eff_s_Bs<<endl; 
  double* eff_s_Bu = eff_syst_Bu->GetY();

  //ARRAYS//


  for(int i = 0; i < n_pt_bins; i++){
    s_errors_bin_Bu[i][0] = luminosity_error;
    s_errors_bin_Bs[i][0] = luminosity_error;

    s_errors_bin_Bu[i][1] = (branching_fraction_error_Bu/branching_fraction_Bu);
    s_errors_bin_Bs[i][1] = (branching_fraction_error_Bs/branching_fraction_Bs);

    s_errors_bin_Bu[i][2] = syst_Bu[i]/raw_Bu_y[i];
    s_errors_bin_Bs[i][2] = syst_Bs[i]/raw_Bs_y[i];

    cout<<"yieldbu: "<< raw_Bu_y[i]<<endl;
    cout<<"yieldbs: "<< raw_Bs_y[i]<<endl;

    s_errors_bin_Bu[i][3] = eff_s_Bu[i]/efficiency_Bu->GetBinContent(i+1);
    s_errors_bin_Bs[i][3] = eff_s_Bs[i]/efficiency_Bs->GetBinContent(i+1);

    //s_errors_bin_Bu[i][3] = eff_s_Bu[i]/efficiency_Bu->GetBinContent(i+1);
    // s_errors_bin_Bs[i][3] = 0.00005;

  }



  ////////

  //x-section values
  double x_sec_Bu[n_pt_bins];
  double x_sec_Bs[n_pt_bins];
  double x_sec0;

  for(int i = 0; i < n_pt_bins; i++)
    {
      //Bu
      n = raw_Bu_y[i];
      eff = efficiency_Bu->GetBinContent(i+1);
      x_sec0 = n/(eff*branching_fraction_Bu*luminosity);
      x_sec_Bu[i] = x_sec0;
      cout << x_sec_Bu[i] << endl;
      cout << endl;

      //Bs
      n = raw_Bs_y[i];
      eff = efficiency_Bs->GetBinContent(i+1);
      x_sec0 = n/(eff*branching_fraction_Bs*luminosity);
      x_sec_Bs[i] = x_sec0;
      cout << x_sec_Bs[i] << endl;
      cout << endl;
    }

  //x-section systematics
  //double x_sec_Bu_syst[n_pt_bins] = error_syst_final(n_pt_bins,n_s_errors,s_errors_bin_Bu);
  //double x_sec_Bs_syst[n_pt_bins] = error_syst_final(n_pt_bins,n_s_errors,s_errors_bin_Bs);
  double *x_sec_Bu_syst = error_syst_final(n_pt_bins,n_s_errors,s_errors_bin_Bu);
  double *x_sec_Bs_syst = error_syst_final(n_pt_bins,n_s_errors,s_errors_bin_Bs);
  
  double syst_Bu_final[n_pt_bins];
  double syst_Bs_final[n_pt_bins];
  
  for (int i=0;i<n_pt_bins;i++){
    syst_Bu_final[i] = x_sec_Bu_syst[i]*x_sec_Bu[i];
    syst_Bs_final[i] = x_sec_Bs_syst[i]*x_sec_Bs[i];
    }
  
  double stat_Bu_final[n_pt_bins];
  double stat_Bs_final[n_pt_bins];
  
  for (int i=0;i<n_pt_bins;i++){
    stat_Bu_final[i] = (stat_Bu[i]/raw_Bu_y[i])* x_sec_Bu[i];
    stat_Bs_final[i] = (stat_Bs[i]/raw_Bs_y[i])* x_sec_Bs[i];
    }
  

  //PLOT B+ X-section
  plot_xsection(n_pt_bins,pT_Bu,pT_min_Bu,pT_max_Bu,x_sec_Bu,stat_Bu_final,syst_Bu_final,"Bu");
  //PLOT Bs X-section
  plot_xsection(n_pt_bins,pT_Bs,pT_min_Bs,pT_max_Bs,x_sec_Bs,stat_Bs_final,syst_Bs_final,"Bs");

 //PRINT VALUES
  for (int i=0;i<n_pt_bins;i++){
    cout<<"bin: "<<i<<" x-sec: "<<x_sec_Bu[i]<<" syst error: "<<syst_Bu_final[i]<<" stat error: "<<stat_Bu_final[i]<<endl;
    cout<<"bin: "<<i<<" x-sec: "<<x_sec_Bs[i]<<" syst error: "<<syst_Bs_final[i]<<" stat error: "<<stat_Bs_final[i]<<endl;

  }

}

//main ends

//returns the histogram's errors per bin
double get_error(int bin, TGraphAsymmErrors * graph_name){
  double error = graph_name->GetErrorY(bin);
  return error;
}
//get_error ends


double get_pt_low(int bin, TGraphAsymmErrors * graph_name){
  double low = graph_name->GetErrorXlow(bin);
  return low;
}

//get_pt_low ends

double get_pt_high(int bin, TGraphAsymmErrors * graph_name){
  double high = graph_name->GetErrorXhigh(bin);
  return high;
}

//get_pt_high ends
 

void plot_xsection(int bin_n,double* pt_m, double* pt_l, double* pt_h, double* x_sec, double* stat, double* syst, TString meson){
  TCanvas c;
  // TMultiGraph* mg = new TMultiGraph();
  
  TGraphAsymmErrors* g_stat = new TGraphAsymmErrors(bin_n,pt_m,x_sec,pt_l,pt_h,stat,stat);
  g_stat->SetTitle("");
  g_stat->SetMarkerColor(4);
  g_stat->SetMarkerStyle(1);
  g_stat->SetLineColor(1);
  //g_stat->GetXaxis()->SetLimits(10,52);
  //g_stat->GetXaxis()->SetRangeUser(10,52);
  g_stat->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  g_stat->GetYaxis()->SetTitle("X-section [nb/GeV]");
  g_stat->SetMinimum(-5e15);
  g_stat->SetMaximum(15e17);
  g_stat->Draw("AP");
  c.Modified();
 

  /*
  double pt_zero[bin_n];
  for (int i=0;i<bin_n;i++) pt_zero[i]= 0.;
  
  TGraphAsymmErrors* g_syst= new TGraphAsymmErrors(bin_n,pt_m,x_sec,pt_zero,pt_zero,syst,syst);
  g_syst->SetTitle("");
  g_syst->SetFillColor(2);
  g_syst->SetFillStyle(3001);
  g_syst->SetMarkerColor(4);
  g_syst->SetMarkerStyle(1);
  g_syst->SetLineColor(2);
  g_stat->GetXaxis()->SetLimits(10,52);
  //g_syst->Draw("2 same");

  
  mg->Add(g_stat);
  mg->Add(g_syst);
  //mg->GetXaxis()->SetLimits(10,52);
  //mg->GetXaxis()->SetRange(10,52);
  mg->Draw("AP");
  mg->GetXaxis()->SetLimits(10,52);
  mg->SetMinimum(-5e15);
  mg->SetMaximum(15e15);
  c.Modified();
  mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  mg->GetYaxis()->SetTitle("X-section [nb/GeV]");
  */
  c.SaveAs(Form("./xresults/x_section_" + meson + ".pdf"));
}
//plot_xsection ends



//function that evaluates the final systematic error
double* error_syst_final(int n_bins,int n_sources_error,double (*source_error_bin)[4]){
  double *final_syst = new double[n_bins];
  //value of the systematic error per bin
  double sum_pow_syst[n_bins];
  //sum of the squares ofthe syst errors
   for(int s = 0; s<n_bins; s++){
    sum_pow_syst[s]=0;
  } 
 
  //loop through the number of bins
  for(int i=0;i<n_bins;i++){
    //loop through the array of sources of error
    for(int k=0;k<n_sources_error;k++){
      sum_pow_syst[i] += pow(source_error_bin[i][k],2); 
    }
    final_syst[i] = sqrt(sum_pow_syst[i]);
  }
  return final_syst;
}
//error_syst_final ends

//function that plots the ratio between the Bu plot and the Bs

