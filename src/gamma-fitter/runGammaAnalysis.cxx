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
#include <iomanip> 
#include <iostream>


void GammaAnalysis( TString name, TH1D* histo, TF1* res, TString outputDir )
{
  std::GammaLineFit fitter(name, histo, res, outputDir);

  //K lines
  fitter.Fit("K42_1525", {1524.7}, {1505.,1545.}, std::kLinear, 0.5, 0.2, BCEngineMCMC::kMedium); // 18.1%, kStep?, pos priot 0.2 -> 0.5
  fitter.Fit("K40_1461",    {1460.8}, {1441.,1481.}, std::kStep,   0.2, 0.2, BCEngineMCMC::kMedium); // 10.7%, kStep?
  //Co lines
  fitter.Fit("Co60_1332",   {1332.5}, {1313.,1353.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); //100.0%
  fitter.Fit("Co60_1173",   {1173.2}, {1153.,1193.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 99.9%
  //Th chain
  fitter.Fit("Ac228_911" ,  { 911.2}, { 891., 931.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 25.8%
  fitter.Fit("Ac228_1588" ,  { 1588.2}, { 1568., 1608.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 3.2%
  //fitter.Fit("Ac228_969" ,  { 969.0}, { 949., 989.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 15.8%, -> in overlap!
  //fitter.Fit("Ac228_338" ,  { 338.3}, { 318., 348.}, std::kQuadratic); // 11.3% check coax!!                   -> in overlap!
  fitter.Fit("Bi212_727",   { 727.3}, { 707., 747.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); //  6.7%
  fitter.Fit("Tl208_2614",  {2614.5}, {2595.,2635.}, std::kFlat,   0.2, 0.2, BCEngineMCMC::kMedium); // 99.8%
  fitter.Fit("Tl208_583",   { 583.2}, { 563., 603.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 85.0% 
  fitter.Fit("Tl208_861",   { 860.6}, { 841., 881.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 12.5% 
  //U chain
  fitter.Fit("Pa234m_1001", {1001.0}, { 981.,1021.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); //  0.8% 
  fitter.Fit("Pb214_352",   { 351.9}, { 332., 372.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 35.6%
  fitter.Fit("Pb214_295",   { 295.2}, { 275., 315.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 35.6%
  fitter.Fit("Bi214_609",   { 609.3}, { 589., 629.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 45.5%
  fitter.Fit("Bi214_1378",   { 1377.7}, { 1358., 1398.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); //  4%
  fitter.Fit("Bi214_1730" ,  { 1729.6}, { 1710., 1750.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 2.9%
  fitter.Fit("Bi214_1764",  {1764.5}, {1744.,1784.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 15.3%
  //fitter.Fit("Bi214_1120",  {1120.3}, {1100.,1140.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 14.9% -> in overlap!
  fitter.Fit("Bi214_1238",  {1238.1}, {1218.,1258.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); //  5.8%
  fitter.Fit("Bi214_2204",  {2204.1}, {2184.,2224.}, std::kFlat,      0.2, 0.2, BCEngineMCMC::kMedium); //  4.9%  
  fitter.Fit("Bi214_2448",   { 2447.9}, { 2428., 2468.}, std::kFlat, 0.2, 0.2, BCEngineMCMC::kMedium); //  1.5%
  //overlap
  fitter.Fit("e+e-_Kr85_514",      { 511.0, 514.0}, { 491., 534.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 0.4%
  fitter.Fit("Pb212_239_Pb214_242",{ 238.6, 242.0}, { 218., 262.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 43.6%,  7.3%
  fitter.Fit("Ac228_338_Pb214_352",{ 338.3, 351.9}, { 318., 372.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 11.3%, 35.6%
  fitter.Fit("Ac228_965_Ac228_969",{ 964.8, 969.0}, { 945., 989.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 5%, 25.8%
  fitter.Fit("Bi214_1120_Zn65_1125",  { 1120.3, 1125.0}, { 1100.,1145.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 14.9%, ?%

}



void GammaAnalysisSingle( TString name, TH1D* histo, TF1* res, TString outputDir )
{
  std::GammaLineFit fitter(name, histo, res, outputDir);

  fitter.Fit("K42_1525", {1524.7}, {1505.,1545.}, std::kLinear, 0.5, 0.2, BCEngineMCMC::kMedium); // 18.1%, kStep?, pos priot 0.2 -> 0.5
  fitter.Fit("K40_1461",    {1460.8}, {1441.,1481.}, std::kStep,   0.2, 0.2, BCEngineMCMC::kMedium); // 10.7%, kStep?
}



void PeakSearchAnalysis( TString name, TH1D* histo, TF1* res, TString outputDir, double min_E, double max_E, double step, double fit_window )
{
  std::GammaLineFit fitter(name, histo, res, outputDir);

  for ( double E0=min_E; E0<=max_E; E0+=step ) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << E0;
    std::string fit_name = stream.str();
    double low_edge = E0-fit_window/2;
    double upp_edge = E0+fit_window/2;
    if (E0 <= 500)            { fitter.Fit(fit_name, {E0}, {low_edge, upp_edge}, std::kQuadratic, 0, 0, BCEngineMCMC::kMedium); }
    if (E0 > 500 && E0<=2000) { fitter.Fit(fit_name, {E0}, {low_edge, upp_edge}, std::kLinear, 0, 0, BCEngineMCMC::kMedium); }
    if (E0 >2000)             { fitter.Fit(fit_name, {E0}, {low_edge, upp_edge}, std::kFlat, 0, 0, BCEngineMCMC::kMedium); }
  }
}



int main (int argc, char ** argv)
{

  // get spectrum
  TString histo_file = (("./tmp/"+std::string(argv[1]))+("-"+std::string(argv[2])))+".root";
  TFile *file = TFile::Open(histo_file, "READ");
  TNamed *output = (TNamed*)file->Get("outputDir");
  TString outputDir = (TString)output->GetTitle();
  TH1D* histo = (TH1D*)file->Get("Spectrum");

  // get resolution curve
  TParameter<double>* a_res = (TParameter<double>*)file->Get("a_res");
  TParameter<double>* b_res = (TParameter<double>*)file->Get("b_res");
  double a_res_double = (double) a_res->GetVal();
  double b_res_double = (double) b_res->GetVal();
  TF1* f_res = new TF1("f_res","sqrt([0]+[1]*x)");
  f_res->SetParameter(0,a_res_double);
  f_res->SetParameter(1,b_res_double);

  std::string single = argv[3];
  std::string analysis = argv[4];
  // run analysis
  if ( single == "true" ) { GammaAnalysisSingle("histo" , histo, f_res, outputDir); }
  else { 
    if (analysis == "gamma-lines") { GammaAnalysis("histo" , histo, f_res, outputDir); }
    if (analysis == "peak-search") { 
        double min_E = std::stod(argv[5]);
        double max_E = std::stod(argv[6]);
        double step  = std::stod(argv[7]);
        double fit_window = std::stod(argv[8]);
        PeakSearchAnalysis("histo" , histo, f_res, outputDir, min_E, max_E, step, fit_window); 
    }
  }
  
}

  
