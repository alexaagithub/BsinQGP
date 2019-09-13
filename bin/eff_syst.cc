#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <iostream>

using namespace std;

#define particle 1

//0 = B+;
//1 = Bs;

int main(){
  
  TFile* f_eff0 = new TFile(particle ?  "./results_eff/Bs_efficiency0.root" :  "./results_eff/Bu_efficiency0.root");
  TFile* f_eff1 = new TFile(particle ?  "./results_eff/Bs_efficiency1.root" : "./results_eff/Bu_efficiency1.root");
 
  //const int n_pt_bins = 8; //run B+ again now for 8 bins
  const int n_pt_bins = 4;
  //double pt_bins[] = {5, 10, 15, 20, 50};

  TEfficiency* efficiency0; //= new TEfficiency("efficiency0", "efficiency0", n_pt_bins, pt_bins);
  TEfficiency* efficiency1; //= new TEfficiency("efficiency1", "efficiency1", n_pt_bins, pt_bins);
  efficiency0 = (TEfficiency*)f_eff0->Get("hist_tot_noweights_clone");
  efficiency1 = (TEfficiency*)f_eff1->Get("hist_tot_weights_clone");

  double eff0;
  double eff1;
  double syst;

  double y_values[n_pt_bins];
   
  double x_values[] = {7.5, 12.5, 17.5, 35};

  //double x_values[] = {6, 8.5, 12.5, 17.5, 25, 35, 45, 55};
 
  double x_errors[] = {2.5, 2.5, 2.5, 15};

  //double x_errors[] = {1, 1.5, 2.5, 2.5, 5, 5, 5, 5};

  double y_errors[] = {0, 0, 0, 0};
 
  //double y_errors[] = {0, 0, 0, 0, 0, 0, 0, 0};
  

  for(int i = 0; i < n_pt_bins; i++) {
    eff0 = efficiency0->GetEfficiency(i + 1);
    eff1 = efficiency1->GetEfficiency(i + 1);
    syst = (eff1 - eff0) / eff0;
    y_values[i] = syst;
    
    cout << "rel err:" << syst <<  endl;
  }

  TGraphErrors* systematic_errors = new TGraphErrors(n_pt_bins, x_values, y_values, x_errors, y_errors);
  TCanvas c;
  systematic_errors->SetMarkerColor(4);
  systematic_errors->SetMarkerStyle(5);
  systematic_errors->Draw("AP");
  systematic_errors->SetTitle("");
  systematic_errors->GetYaxis()->SetTitle("systematic uncertainty");
  systematic_errors->GetXaxis()->SetTitle("p_{T} [GeV]");

  c.SaveAs(particle ? "./results_eff/Bs_systematic_error.gif" : "./results_eff/Bu_systematic_error.gif");
  c.SaveAs(particle ? "./results_eff/Bs_systematic_error.pdf" : "./results_eff/Bu_systematic_error.pdf");
  

  TFile* f1 = new TFile(particle ? "./results_eff/Bs_efficiency_systematic_errors.root" : "./results_eff/Bu_efficiency_systematic_errors.root", "recreate");
  f1->cd();
  systematic_errors->Write();
  f1->Write();
  f1->ls();
  f1->Close();
  
}
