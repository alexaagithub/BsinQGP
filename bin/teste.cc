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

double read_weights(TString var, double var_value, int meson);
double getWeight(double var_value, TH1D* h_weight);

#define particle 0
//0 = B+;   
//1 = Bs;

int main(){

  //MC cuts
 
  TFile* f_mc_cuts = new TFile(particle ? "/lstore/cms/ev19u032/prefiltered_trees_final/selected_mc_ntphi_PbPb_2018_corrected_nocuts_BDT.root" : "/lstore/cms/ev19u032/prefiltered_trees_final/selected_mc_ntKp_PbPb_2018_corrected_BDT.root");
  TTree* t_cuts = (TTree*)f_mc_cuts->Get(particle ? "ntphi" : "ntKp");
  //TH1D* histo_BDT = (TH1D*)t_cuts->Get("BDT_total");

  

  //acceptance cut only
  TFile* f_mc_nocuts = new TFile(particle ? "/lstore/cms/ev19u032/prefiltered_trees_final/acceptance_only_selected_mc_ntphi_PbPb_2018_corrected_BDT.root" : "/lstore/cms/ev19u032/prefiltered_trees_final/acceptance_only_selected_mc_ntKp_PbPb_2018_corrected_BDT.root");
  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get(particle ? "ntphi" : "ntKp");
  //TH1D* histo_BDT_nocuts = (TH1D*)t_nocuts->Get("BDT_total");

 
  double pt_bins[] = {5, 10, 15, 20, 50};
 
  int n_pt_bins = 4;
  

  //TH1D* hist_tot_noweights = new TH1D("hist_tot_noweights", "hist_tot_noweights", n_pt_bins, pt_bins);
  // TH1D* hist_passed_noweights = new TH1D("hist_passed_noweights", "hist_passed_noweights", n_pt_bins, pt_bins);
  TH1D* histo_teste = new TH1D("histo_teste", "histo_teste", n_pt_bins, pt_bins);
  TH1D* histo_teste3 = new TH1D("histo_teste3", "histo_teste3", n_pt_bins, pt_bins);

  // TH1D* hist_tot_weights = new TH1D("hist_tot_weights", "hist_tot_weights", n_pt_bins, pt_bins);
  //TH1D* hist_passed_weights = new TH1D("hist_passed_weights", "hist_passed_weights", n_pt_bins, pt_bins);
  TH1D* histo_teste2 = new TH1D("histo_teste2", "histo_teste2", n_pt_bins, pt_bins);
  TH1D* histo_teste4 = new TH1D("histo_teste4", "histo_teste4", n_pt_bins, pt_bins);

  //BDT - acceptance cut
  double bdt_1;
  //float bdt_1;
  float bpt1;

  t_nocuts->SetBranchAddress("BDT_total", &bdt_1);
  t_nocuts->SetBranchAddress("Bpt", &bpt1);

  double weight = 1;

for(int evt = 0; evt < 5; evt++)
  //for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);
      // hist_tot_noweights->Fill(bdt_1);
      weight = read_weights("BDT_total", bdt_1, particle);
      // hist_tot_weights->Fill(bdt_1, weight);
      histo_teste->Fill(weight);
      histo_teste3->Fill(6.0);
    }

  /*
  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);
      hist_tot_noweights->Fill(bpt1);
      cout<<"a partícula: "<<particle<<endl;
      weight = read_weights("Bpt", bpt1, particle);
      cout<<"weight: "<<weight<<endl;
      hist_tot_weights->Fill(bpt1, weight);
    }
  */
  //descomentar dps
  /*
  TCanvas tot_noweights;
  hist_tot_noweights->Draw();
  tot_noweights.SaveAs(particle ? "./results_eff/Bs_tot_noweights.pdf" : "./results_eff/Bu_tot_noweights.pdf");
 

  TCanvas tot_weights;
  hist_tot_weights->Draw();
  tot_weights.SaveAs(particle ? "./results_eff/Bs_totweights.pdf" : "./results_eff/Bu_totweights.pdf");
  */
  //BDT - cuts
  double bdt_2;
  //float bdt_2;
  float bpt2;

  t_cuts->SetBranchAddress("BDT_total", &bdt_2);
  t_cuts->SetBranchAddress("Bpt", &bpt2);

  double weight2 = 1;

for(int evt = 0; evt < 50; evt++)
  //for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);
      //hist_passed_noweights->Fill(bdt_2);
      //cout<<"bdt: "<<bdt_2<<endl;
      weight2 = read_weights("BDT_total", bdt_2, particle);
      // cout<<"weight2: "<<weight2<<endl;
      // hist_passed_weights->Fill(bdt_2, weight2);
      histo_teste2->Fill(bpt2,weight2);
      histo_teste4->Fill(30);

      cout << "test2: " << histo_teste2->GetBinContent(3) << endl;
      cout << "test4: " << histo_teste4->GetBinContent(4) << endl;
    }



  
  /*
  for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);
      hist_passed_noweights->Fill(bpt2);
      cout<<"particle: "<<particle<<endl;
      weight2 = read_weights("Bpt", bpt2, particle);
      cout<<"weight2"<<weight2<<endl;
      hist_passed_weights->Fill(bpt2, weight2);
    }
  */

  //descomentar dps
  /*
  TCanvas passed_noweights;
  hist_passed_noweights->Draw();
  passed_noweights.SaveAs(particle ? "./results_testinho/Bs_passed_noweights.pdf" : "./results_testinho/Bu_passed_noweights.pdf");
 
  
  TCanvas passed_weights;
  hist_passed_weights->Draw();
  passed_weights.SaveAs(particle ? "./results_testinho/Bs_passed_weights.pdf" : "./results_testinho/Bu_passed_weights.pdf");
  
  
  TEfficiency* efficiency0 = new TEfficiency(*hist_passed_noweights, *hist_tot_noweights);
  TCanvas c0;
  efficiency0->Draw("AP");
  c0.SaveAs(particle ? "./results_testinho/Bs_efficiency0.pdf" : "./results_testinho/Bu_efficiency0.pdf");
 

  TEfficiency* efficiency1 = new TEfficiency(*hist_passed_weights, *hist_tot_weights);
  TCanvas c1;
  efficiency1->Draw("AP");
  c1.SaveAs(particle ? "./results_testinho/Bs_efficiency1.pdf" : "./results_testinho/Bu_efficiency1.pdf");

  //create root files
 
  TFile* f0 = new TFile(particle ? "./results_testinho/Bs_efficiency0.root" : "./results_testinho/Bu_efficiency0.root" , "recreate");
  
  f0->cd();
  efficiency0->Write();
  f0->Write();

  TFile* f1 = new TFile(particle ? "./results_testinho/Bs_efficiency1.root" : "./results_testinho/Bu_efficiency1.root" , "recreate");
 
  f1->cd();
  efficiency1->Write();
  f1->Write();

  for(int i = 1; i < n_pt_bins + 1; i++)
    {
      cout << efficiency0->GetEfficiency(i) << endl;
    }

  cout << endl;

  for(int i = 1; i < n_pt_bins + 1; i++)
    {
      cout << efficiency1->GetEfficiency(i) << endl;
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
  */

  TCanvas histos_teste;
  histo_teste->Draw();
  histos_teste.SaveAs(particle ? "./results_testinho/Bs_histoteste.pdf" : "./results_testinho/histoteste.pdf");
 
  
  TCanvas histos_teste2;
  histo_teste2->Draw();
  histos_teste2.SaveAs(particle ? "./results_testinho/Bs_histoteste2.pdf" : "./results_testinho/histoteste2.pdf");
 
  TCanvas histos_teste3;
  histo_teste3->Draw();
  histos_teste3.SaveAs(particle ? "./results_testinho/Bs_histoteste3.pdf" : "./results_testinho/histoteste3.pdf");
  
  TCanvas histos_teste4;
  histo_teste4->Draw();
  histos_teste4.SaveAs(particle ? "./results_testinho/Bs_histoteste4.pdf" : "./results_testinho/histoteste4.pdf");

  delete histo_teste;
  delete histo_teste2;
  delete histo_teste3;
  delete histo_teste4;
 
  return 0;
  
}

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
