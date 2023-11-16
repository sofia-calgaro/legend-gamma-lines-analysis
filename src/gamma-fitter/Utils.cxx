#include "Utils.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "../settings/json.hpp"
using namespace nlohmann;
#include <sys/stat.h>

bool fileExists(const std::string& filePath) {
    struct stat buffer;
    return (stat(filePath.c_str(), &buffer) == 0);
}


void write_to_json(std::string filePath, std::string name_fit, ordered_json foutput) {
/*  ---> this part is having problems, adding too '}}' at the eof
  // check if the json file aready exists
  if (fileExists(filePath)) {
    std::fstream json_file;
    json_file.open(filePath, std::ios::in | std::ios::out);
    // parse the existing JSON from the file
    ordered_json foutput_tot = json::parse(json_file);
    // check if the 'name_fit' key already exists; if so, remove it
    //if ( foutput_tot.find(name_fit) != foutput_tot.end() ) {
    //  foutput_tot.erase(name_fit);
    //}
    // add the new key
    //foutput_tot.update(foutput); //, true);
    foutput_tot[name_fit] = foutput;
    // the array is empty again
    json_file.clear();
    // let's write from the top of the file
    json_file.seekp(0);
    // set indent to 2
    json_file << foutput_tot.dump(2);
    json_file.close();
  }
  // json file does not exist, let's write in it for the first time!
  else{
    std::fstream json_file;
    json_file.open(filePath, std::ios::out);
    json foutput_tot;
    foutput_tot[name_fit] = foutput;
    json_file << foutput_tot.dump(2);
    json_file.close();
  }
*/

  std::fstream json_file;
  json_file.open(filePath, std::ios::out);
  json foutput_tot;
  foutput_tot[name_fit] = foutput;
  json_file << foutput_tot.dump(2);
  json_file.close();

  std::ifstream jsonFile(filePath);
  if (jsonFile.is_open()) {
    std::string fileContent((std::istreambuf_iterator<char>(jsonFile)), std::istreambuf_iterator<char>());
    std::cout << "File Content:\n" << fileContent << std::endl;
  } else {
    std::cerr << "Error: Could not open file for reading." << std::endl;
  }
}


std::vector<int> CheckPosterior(std::string fit_name, std::string path_to_output, int bkg) {
  // bkg = 0,1,2 - the posterior check for the step background is still not implemented
  // The idea is to perform the check for the peak search analysis only, where we don't expect to fit with a step function
  
  std::string name = path_to_output + "/histo_marginalized." + fit_name + ".root";
  TFile *file = new TFile(name.c_str()); 

  // E0 marginalzied posterior histogram
  std::string histo_E0 = "h1_histo_fitter_model_parameter_intensity0";
  TH1D *hE0 = (TH1D*) file->Get(histo_E0.c_str());

  // p0 marginalzied posterior histogram
  std::string histo_p0 = "h1_histo_fitter_model_parameter_par00";
  TH1D *hp0 = (TH1D*) file->Get(histo_p0.c_str());

  // p1 marginalzied posterior histogram (if exists)
  std::string histo_p1 = "h1_histo_fitter_model_parameter_par1";
  TH1D *hp1 = new TH1D();
  if ( bkg>=1 ) { hp1 = (TH1D*) file->Get(histo_p1.c_str()); }

  // p2 marginalzied posterior histogram (if exists)
  std::string histo_p2 = "h1_histo_fitter_model_parameter_par2";
  TH1D *hp2 = new TH1D();
  if ( bkg==2 ) { hp2 = (TH1D*) file->Get(histo_p2.c_str()); }
  
  // study of the 1st and last bin content
  int p0BinMax=0, p1BinMax=0, p2BinMax=0;
  p0BinMax = hp0->GetMaximumBin();
  if ( bkg>=1 ) {
    p1BinMax = hp1->GetMaximumBin();
    if ( bkg==2 ) { p2BinMax = hp2->GetMaximumBin(); }
  }

  // Get bins that correspond to 5% of all bins
  int tot_bins = hp0->GetNbinsX();
  int bin_5_perc = tot_bins*0.05;

  // Mean of first 5% bins 
  double p0FirstSum=0, p1FirstSum=0, p2FirstSum=0;
  double p0FirstMean=0, p1FirstMean=0, p2FirstMean=0;
  for ( int i=1; i<=bin_5_perc; i++ ) {
    p0FirstSum += hp0->GetBinContent(i);
    if ( bkg>=1 ) {
      p1FirstSum += hp1->GetBinContent(i);
      if ( bkg==2 ) { p2FirstSum += hp2->GetBinContent(i); }
    }
  }
  p0FirstMean = p0FirstSum/bin_5_perc;
  if ( bkg>=1 ) {
    p1FirstMean = p1FirstSum/bin_5_perc;
    if ( bkg==2 ) { p2FirstMean = p2FirstSum/bin_5_perc; }
  }

  // Mean of last 5% bins 
  double p0LastSum=0, p1LastSum=0, p2LastSum=0;
  double p0LastMean=0, p1LastMean=0, p2LastMean=0;
  for ( int i=(hp0->GetNbinsX())-(bin_5_perc-1); i<=hp0->GetNbinsX(); i++ ) { p0LastSum += hp0->GetBinContent(i); }
  if ( bkg>=1 ) { 
    for ( int i=(hp1->GetNbinsX())-(bin_5_perc-1); i<=hp1->GetNbinsX(); i++ ) { p1LastSum += hp1->GetBinContent(i); } }
    if ( bkg==2 ) { for ( int i=(hp2->GetNbinsX())-(bin_5_perc-1); i<=hp2->GetNbinsX(); i++ ) { p2LastSum += hp2->GetBinContent(i); }
  }
  p0LastMean = p0LastSum/bin_5_perc;
  if ( bkg>=1 ) {
    p1LastMean = p1LastSum/bin_5_perc;
    if ( bkg==2 ) { p2LastMean = p2LastSum/bin_5_perc; }
  }

  double p0FirstBin=0, p1FirstBin=0, p2FirstBin=0;
  double p0LastBin=0, p1LastBin=0, p2LastBin=0;
  p0FirstBin = hp0->GetBinContent(1);
  p0LastBin = hp0->GetBinContent(hp0->GetNbinsX());
  p1FirstBin = hp1->GetBinContent(1);
  p1LastBin = hp1->GetBinContent(hp1->GetNbinsX());
  p2FirstBin = hp2->GetBinContent(0); p2LastBin = hp2->GetBinContent(hp2->GetNbinsX()); 
  double E0LastBin = hE0->GetBinContent(hE0->GetNbinsX());

  int ck0=0, ck1=0, ck2=0, ckE0=0;
  if ( p0FirstMean/(hp0->GetBinContent(p0BinMax))>0.1 || p0LastMean/(hp0->GetBinContent(p0BinMax))>0.1 ) { ck0 = 1; }
  if ( p1FirstMean/(hp1->GetBinContent(p1BinMax))>0.1 || p1LastMean/(hp1->GetBinContent(p1BinMax))>0.1 ) { ck1 = 1; }
  if ( p2FirstMean/(hp2->GetBinContent(p2BinMax))>0.1 || p2LastMean/(hp2->GetBinContent(p2BinMax))>0.1 ) { ck2 = 1; }
  if ( E0LastBin!=0 ) { ckE0 = 1; }
  std::vector<int> results = {ck0, ck1, ck2, ckE0};
      
  return results;
}
