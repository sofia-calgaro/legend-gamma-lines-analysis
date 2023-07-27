// MODELS
#include "GammaFitFunctions.h"
#include "GammaLineFit.h"

// ROOT
#include <TMath.h>
#include <TString.h>
#include <TF1.h>
#include <TParameter.h>

//C++
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>



void GammaAnalysis( TString name, TH1D* histo, TF1* res, TString outputDir )
{
  std::GammaLineFit fitter(name,histo,res, outputDir);
  //misc
  fitter.Fit("K42_1525", {1524.7}, {1505.,1545.}, std::kLinear, 0.5, 0.2, BCEngineMCMC::kMedium); // 18.1%, kStep?, pos priot 0.2 -> 0.5
  /*
  fitter.Fit("K40_1461",    {1460.8}, {1441.,1481.}, std::kStep,   0.2, 0.2, BCEngineMCMC::kMedium); // 10.7%, kStep?
  fitter.Fit("Co60_1332",   {1332.5}, {1313.,1353.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); //100.0%
  fitter.Fit("Co60_1173",   {1173.2}, {1153.,1193.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 99.9%
  //th chain
  fitter.Fit("Ac228_911" ,  { 911.2}, { 891., 931.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 25.8%
  fitter.Fit("Ac228_969" ,  { 969.0}, { 949., 989.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 15.8%, 948!!!!!!
  //fitter.Fit("Ac228_338" ,  { 338.3}, { 318., 348.}, std::kQuadratic); // 11.3% check coax!!
  fitter.Fit("Bi212_727",   { 727.3}, { 707., 747.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); //  6.7%
  fitter.Fit("Tl208_2614",  {2614.5}, {2595.,2635.}, std::kFlat,   0.2, 0.2, BCEngineMCMC::kMedium); // 99.8%
  fitter.Fit("Tl208_583",   { 583.2}, { 563., 603.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 85.0% 
  fitter.Fit("Tl208_861",   { 860.6}, { 841., 881.}, std::kLinear, 0.2, 0.2, BCEngineMCMC::kMedium); // 12.5% 
  //u chain
  fitter.Fit("Pa234m_1001", {1001.0}, { 981.,1021.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); //  0.8% check coax!
  fitter.Fit("Pb214_352",   { 351.9}, { 332., 372.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 35.6%
  fitter.Fit("Pb214_295",   { 295.2}, { 275., 315.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 35.6%
  fitter.Fit("Bi214_609",   { 609.3}, { 589., 629.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 45.5%
  fitter.Fit("Bi214_1764",  {1764.5}, {1744.,1784.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 15.3%
  fitter.Fit("Bi214_1120",  {1120.3}, {1100.,1140.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 14.9%
  fitter.Fit("Bi214_1238",  {1238.1}, {1218.,1258.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); //  5.8%
  fitter.Fit("Bi214_2204",  {2204.1}, {2184.,2224.}, std::kFlat,      0.2, 0.2, BCEngineMCMC::kMedium); //  4.9%  
  //overlap
  fitter.Fit("e+e-_Kr85_514",      { 511.0, 514.0}, { 491., 534.}, std::kLinear,    0.2, 0.2, BCEngineMCMC::kMedium); // 0.4%
  fitter.Fit("Pb212_239_Pb214_242",{ 238.6, 242.0}, { 218., 262.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 43.6%,  7.3%
  fitter.Fit("Ac228_338_Pb214_352",{ 338.3, 351.9}, { 318., 372.}, std::kQuadratic, 0.2, 0.2, BCEngineMCMC::kMedium); // 11.3%, 35.6%
  */
}



int main() 
{

  TString histo_file = "./tmp/tmp-spectra.root";
  TFile *file = TFile::Open(histo_file, "READ");
  TNamed *output = (TNamed*)file->Get("outputDir");
  TString outputDir = (TString)output->GetTitle();
  TH1D* hEnrBEGeRaw = (TH1D*)file->Get("Spectrum");
  TParameter<double>* a_res = (TParameter<double>*)file->Get("a_res");
  TParameter<double>* b_res = (TParameter<double>*)file->Get("b_res");
  double a_res_double = (double) a_res->GetVal();
  double b_res_double = (double) b_res->GetVal();
  TF1* f_res = new TF1("f_res","sqrt([0]+[1]*x)");
  f_res->SetParameter(0,a_res_double);
  f_res->SetParameter(1,b_res_double);
  GammaAnalysis("EnrBEGeRaw" , hEnrBEGeRaw, f_res, outputDir);
  
}


