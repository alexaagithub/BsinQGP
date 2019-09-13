#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TString.h>
#include <iostream>

using namespace std;

//double read_weights(TString var, double var_value, int meson);
//double getWeight(double var_value, TH1D* h_weight);

double read_weights(TH1D* histo, double var_value);
double getWeight(double var_value, TH1D* histo);
TH1D* read_histos_weights(TString var, int meson);

#define particle 1
//0 = B+;   
//1 = Bs;

/*
int main(int argc, char *argv[]){

  int particle; 
  std::string argument = argv[1];
  std::stringstream convert;
  
  convert << argv[1];
  convert >> particle;

  cout << "---" << particle << endl;
*/

int main(){

  
  //MC cuts
 
  TFile* f_mc_cuts = new TFile(particle ? "/lstore/cms/ev19u032/prefiltered_trees_final/selected_mc_ntphi_PbPb_2018_corrected_nocuts_BDT.root" : "/lstore/cms/ev19u032/prefiltered_trees_final/selected_mc_ntKp_PbPb_2018_corrected_BDT.root");
  TTree* t_cuts = (TTree*)f_mc_cuts->Get(particle ? "ntphi" : "ntKp");
  //TH1D* histo_BDT = (TH1D*)t_cuts->Get("BDT_total");

  

  //acceptance cut only
  TFile* f_mc_nocuts = new TFile(particle ? "/lstore/cms/ev19u032/prefiltered_trees_final/acceptance_only_selected_mc_ntphi_PbPb_2018_corrected_nocuts_BDT.root" : "/lstore/cms/ev19u032/prefiltered_trees_final/acceptance_only_selected_mc_ntKp_PbPb_2018_corrected_BDT.root");
  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get(particle ? "ntphi" : "ntKp");
  //TH1D* histo_BDT_nocuts = (TH1D*)t_nocuts->Get("BDT_total");
  
  //variable's histogram
  
  TH1D* histos_BDT = read_histos_weights("BDT_total",particle);
  //return 0;
  // TCanvas a;
  // histos_BDT->Draw();
  //a.SaveAs("./teste2.pdf");
  //  return 0;
  
  int n_pt_bins = 4;
  //int n_pt_bins = 8;
  double pt_bins[] = {5, 10, 15, 20, 50};
  //double pt_bins[] = {5 ,7, 10, 15, 20, 30, 40, 50, 60};
 
  
  TH1D* hist_tot_noweights = new TH1D("hist_tot_noweights", "hist_tot_noweights", n_pt_bins, pt_bins);
  TH1D* hist_passed_noweights = new TH1D("hist_passed_noweights", "hist_passed_noweights", n_pt_bins, pt_bins);
 

  TH1D* hist_tot_weights = new TH1D("hist_tot_weights", "hist_tot_weights", n_pt_bins, pt_bins);
  TH1D* hist_passed_weights = new TH1D("hist_passed_weights", "hist_passed_weights", n_pt_bins, pt_bins);
 
 
 
  //BDT - acceptance cut
  double bdt_1;
  float bpt1;

  t_nocuts->SetBranchAddress("BDT_total", &bdt_1);
  t_nocuts->SetBranchAddress("Bpt", &bpt1);

  double weight = 1;

  cout << endl;  
  int nevt1 = t_nocuts->GetEntries();
  for(int evt = 0; evt < nevt1; evt++)    {
      t_nocuts->GetEntry(evt);
      hist_tot_noweights->Fill(bpt1);
      //weight = read_weights("BDT_total", bdt_1, particle);
      weight = read_weights(histos_BDT, bdt_1);
      //cout<<"weight: "<<weight<<endl;
      hist_tot_weights->Fill(bpt1, weight);
      //if (weight!=1){
      //cout<<"o 1 é diferente!"<<endl;
      //}
      if(evt%1000==0) cout << "." << std::flush;
    }
  cout << endl;  

 
  
  TCanvas tot_noweights;
  hist_tot_noweights->Draw();
  tot_noweights.SaveAs(particle ? "./results_eff/Bs_tot_noweights.pdf" : "./results_eff/Bu_tot_noweights.pdf");
 

  TCanvas tot_weights;
  hist_tot_weights->Draw();
  tot_weights.SaveAs(particle ? "./results_eff/Bs_totweights.pdf" : "./results_eff/Bu_totweights.pdf");
 
  //BDT - cuts
  double bdt_2;
  float bpt2;

  t_cuts->SetBranchAddress("BDT_total", &bdt_2);
  t_cuts->SetBranchAddress("Bpt", &bpt2);

  double weight2 = 1;

  cout << endl;  

  int nevt2 = t_cuts->GetEntries();
  for(int evt = 0; evt < nevt2; evt++)    {
      t_cuts->GetEntry(evt);
      hist_passed_noweights->Fill(bpt2);
      //cout<<"bdt: "<<bdt_2<<endl;
      // weight2 = read_weights("BDT_total", bdt_2, particle);
      weight2 = read_weights(histos_BDT, bdt_2);
      // cout<<"weight2: "<<weight2<<endl;
      hist_passed_weights->Fill(bpt2, weight2);
      //if (weight2 !=1){
      //cout<<"o 2 é diferente!"<<endl;
      // }

      if(evt%1000==0) cout << "." << std::flush;

    }

  cout << endl;  


  
  TCanvas passed_noweights;
  hist_passed_noweights->Draw();
  passed_noweights.SaveAs(particle ? "./results_eff/Bs_passed_noweights.pdf" : "./results_eff/Bu_passed_noweights.pdf");
 
  
  TCanvas passed_weights;
  hist_passed_weights->Draw();
  passed_weights.SaveAs(particle ? "./results_eff/Bs_passed_weights.pdf" : "./results_eff/Bu_passed_weights.pdf");
  
  
  TEfficiency* efficiency0 = new TEfficiency(*hist_passed_noweights, *hist_tot_noweights);
  TCanvas c0;
  efficiency0->Draw("AP");
  c0.SaveAs(particle ? "./results_eff/Bs_efficiency0.pdf" : "./results_eff/Bu_efficiency0.pdf");
 

  TEfficiency* efficiency1 = new TEfficiency(*hist_passed_weights, *hist_tot_weights);
  TCanvas c1;
  efficiency1->Draw("AP");
  c1.SaveAs(particle ? "./results_eff/Bs_efficiency1.pdf" : "./results_eff/Bu_efficiency1.pdf");

  //create root files
 
  TFile* f0 = new TFile(particle ? "./results_eff/Bs_efficiency0.root" : "./results_eff/Bu_efficiency0.root" , "recreate");
  
  f0->cd();
  efficiency0->Write();
  f0->Write();

  TFile* f1 = new TFile(particle ? "./results_eff/Bs_efficiency1.root" : "./results_eff/Bu_efficiency1.root" , "recreate");
 
  f1->cd();
  efficiency1->Write();
  f1->Write();

  for(int i = 1; i < n_pt_bins + 1; i++)    {
    cout << "eff0 " << i << " " << efficiency0->GetEfficiency(i) << endl;
    cout << "eff1 " << i << " " << efficiency1->GetEfficiency(i) << endl;
  }

  delete hist_tot_noweights;
  delete hist_passed_noweights;
  delete hist_tot_weights;
  delete hist_passed_weights;

  f_mc_cuts->Close();
  delete f_mc_cuts;
  f_mc_nocuts->Close();
  delete f_mc_nocuts;

  f0->ls();
  f0->Close();
  f1->ls();
  f1->Close();
 
  
  return 0;
  
}

  
//receives the variable, returns its histogram.
  TH1D* read_histos_weights(TString variable, int meson){
 
  TString input_file = meson ? "/lstore/cms/ev19u032/weights/weights_Bs.root" :"/lstore/cms/ev19u032/weights/weights_Bu.root";
 
  TFile* f_wei = new TFile(input_file);

  TH1D* histo_variable = (TH1D*)f_wei->Get(Form("weights_"+variable));


  TCanvas c;
  histo_variable->Draw();
  c.SaveAs("histo_variable.png");
  //  f_wei->Close();
  //delete f_wei;

  return histo_variable;

  }


//receives the variable's histogram and the event;  returns the weights
double read_weights(TH1D* histo, double var_value){

  double weight;
  double variable_min;
  double variable_max;

  variable_min = histo->GetXaxis()->GetXmin();
  variable_max = histo->GetXaxis()->GetXmax();  
  
  //if the event is not in the range its weight is 1.
  if(var_value>=variable_min && var_value<=variable_max){  
    weight = getWeight(var_value,histo);
  }
  else {
    weight = 1;
  }


  return weight;
}

//auxiliar function
double getWeight(double var_value, TH1D* histo){
  int bin = histo->FindBin(var_value);
  return histo->GetBinContent(bin);
}


/*
double read_weights(TString variable, double var_value, int meson){
 
  TString input_file = meson ? "/lstore/cms/ev19u032/weights/weights_Bs.root" :"/lstore/cms/ev19u032/weights/weights_Bu.root";
 
  TFile* f_wei = new TFile(input_file, "read");

  TH1D* histo_variable = (TH1D*)f_wei->Get(Form("weights_"+variable));

  double weight;
  double variable_min;
  double variable_max;

  variable_min = histo_variable->GetXaxis()->GetXmin();
  variable_max = histo_variable->GetXaxis()->GetXmax();  
  
  //if the event is not in the range its weight is 1.
  if(var_value>=variable_min && var_value<=variable_max){  
    weight = getWeight(var_value,histo_variable);
  }
  else {
    weight = 1;
  }

  f_wei->Close();
  delete f_wei;

  return weight;
}

//auxiliar function
double getWeight(double var_value, TH1D* h_weight){
  int bin = h_weight->FindBin(var_value);
  return h_weight->GetBinContent(bin);
}

 */

