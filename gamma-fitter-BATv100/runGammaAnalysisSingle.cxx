// MODELS
#include "GammaFitFunctions.h"
#include "GammaLineFit.h"

// ROOT
#include <TMath.h>
#include <TString.h>
#include <TF1.h>
#include <TParameter.h>

//C++
#include <string>


void GammaAnalysisSingle( TString name, TH1D* histo, TF1* res, TString outputDir )
{
  std::GammaLineFit fitter(name,histo,res, outputDir);

  // K lines
  fitter.Fit("K42_1525", {1524.7}, {1505.,1545.}, std::kLinear, 0.5, 0.2, BCEngineMCMC::kMedium); // 18.1%, kStep?, pos priot 0.2 -> 0.5
  fitter.Fit("K40_1461",    {1460.8}, {1441.,1481.}, std::kStep,   0.2, 0.2, BCEngineMCMC::kMedium); // 10.7%, kStep?
}



int main (int argc, char ** argv)
{

  TString histo_file = (("./tmp/"+std::string(argv[1]))+("-"+std::string(argv[2])))+".root";
  TFile *file = TFile::Open(histo_file, "READ");
  TNamed *output = (TNamed*)file->Get("outputDir");
  TString outputDir = (TString)output->GetTitle();
  TH1D* histo = (TH1D*)file->Get("Spectrum");
  TParameter<double>* a_res = (TParameter<double>*)file->Get("a_res");
  TParameter<double>* b_res = (TParameter<double>*)file->Get("b_res");
  double a_res_double = (double) a_res->GetVal();
  double b_res_double = (double) b_res->GetVal();
  TF1* f_res = new TF1("f_res","sqrt([0]+[1]*x)");
  f_res->SetParameter(0,a_res_double);
  f_res->SetParameter(1,b_res_double);
  GammaAnalysisSingle("histo" , histo, f_res, outputDir);
  
}


