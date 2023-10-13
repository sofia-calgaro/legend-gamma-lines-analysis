#include "GammaLineFit.h"
#include <iomanip>

#include "../settings/json.hpp"
using namespace nlohmann;

namespace std {

int GammaLineFit::Fit( vector<double> lines, pair<double,double> range,
  GammaBackground backgroundType, double linePosPrior, double fwhmPrior,
  BCEngineMCMC::Precision precision )
{
  if(fFitPerformed) Reset();

  //! create proto fit function
  GammaProtoPolynom* protoPolynom;
  switch( backgroundType) {
    case kStep:      protoPolynom = new GammaProtoPolynom(0,0,lines);     break;
    case kLinear:    protoPolynom = new GammaProtoPolynom(1,range.first); break;
    case kQuadratic: protoPolynom = new GammaProtoPolynom(2,range.first); break;
    default:         protoPolynom = new GammaProtoPolynom(0);             break;
  }
  GammaProtoFunctionSum* protoFitFunction = new GammaProtoFunctionSum({protoPolynom});
  for(unsigned int i=0;i<lines.size();i++) protoFitFunction->AddFunction(new GammaProtoGauss(i));

  //! create proto background function
  vector<pair<double,double>> skip;
  for(auto kv : lines) skip.push_back({kv-fRes->Eval(kv)/2.35*3,kv+fRes->Eval(kv)/2.35*3}); //+-3 sigma
  GammaProtoInterruptedFunction* protoBkgFunction = new GammaProtoInterruptedFunction(protoPolynom,skip);

  //! collect proto functions
  fProtoGarbage = protoFitFunction->GetFunctions();
  fProtoGarbage.push_back(protoFitFunction);
  fProtoGarbage.push_back(protoBkgFunction);

  //! create functions
  fFitFunction = new GammaFitFunction("fitfunc",protoFitFunction,range);
  fBkgFunction = new GammaFitFunction("bkgfunc",protoBkgFunction,range);
  int nBkgPars = fBkgFunction->GetNpar();

  //! create fit histo
  pair<int,int> rangeBins = { fHist->GetXaxis()->FindBin(range.first),
                              fHist->GetXaxis()->FindBin(range.second) };
  int nBins = rangeBins.second - rangeBins.first;
  TH1D* fFitHist = new TH1D("fitHist","",nBins,range.first,range.second);
  for(int i=1;i<=nBins;i++) fFitHist->SetBinContent(i,fHist->GetBinContent(rangeBins.first+i-1));
  double binning = fFitHist->GetBinWidth(1); //BAT treats as cts / 1 keV

  //! evaluate preliminary background parameters
  fFitHist->Fit(fBkgFunction,"RL");
 
  //! evaluate preliminary line intensities
  for(int i=0;i<nBkgPars;i++) fFitFunction->FixParameter(i,fBkgFunction->GetParameter(i));
  for(unsigned int i=0;i<lines.size();i++) {
    fFitFunction->FixParameter(nBkgPars+3*i+0,lines.at(i));
    fFitFunction->FixParameter(nBkgPars+3*i+1,fRes->Eval(lines.at(i)));
    fFitFunction->SetParLimits(nBkgPars+3*i+2,0.,max(5.,1.5*( //counts in +-3sigma w/o bkg +50%
      fFitHist->Integral( fFitHist->FindBin(lines.at(i)-fRes->Eval(lines.at(i))/2.35*3) ,
                          fFitHist->FindBin(lines.at(i)+fRes->Eval(lines.at(i))/2.35*3) ) -
      fBkgFunction->Integral( lines.at(i)-fRes->Eval(lines.at(i))/2.35*3 ,
                              lines.at(i)+fRes->Eval(lines.at(i))/2.35*3 ))));
  }
  fFitHist->Fit(fFitFunction,"RL");

  //! set parameter ranges using preliminaries
  for(int i=0;i<nBkgPars;i++) {
    fFitFunction->SetParLimits(i, //+-4*uncertainty
      (fBkgFunction->GetParameter(i)-4*fBkgFunction->GetParError(i))/binning, 
      (fBkgFunction->GetParameter(i)+4*fBkgFunction->GetParError(i))/binning);
  }
  for(unsigned int i=0;i<lines.size();i++) {
    fFitFunction->SetParLimits(nBkgPars+3*i+0, //+-4*prior width
     lines.at(i)-linePosPrior*3 ,
     lines.at(i)+linePosPrior*3 );
    fFitFunction->SetParLimits(nBkgPars+3*i+1, //+-4*prior width
     fRes->Eval(lines.at(i))-fwhmPrior*3, 
     fRes->Eval(lines.at(i))+fwhmPrior*3 );
    fFitFunction->SetParLimits(nBkgPars+3*i+2, //preliminary parameters +-4*uncertainty
     max(0.,(fFitFunction->GetParameter(nBkgPars+3*i+2)-
     4*fFitFunction->GetParError(nBkgPars+3*i+2))/binning),
     (fFitFunction->GetParameter(nBkgPars+3*i+2)+
     4*fFitFunction->GetParError(nBkgPars+3*i+2))/binning); 
  }

  //! create hist fitter and set priors
  // Create a copy of the TH1D histogram (needed for v1.0.0)
  TH1D histCopy = *fFitHist;
  // Create a copy of the GammaFitFunction (needed for v1.0.0)
  GammaFitFunction funcCopy = *fFitFunction;
  fHistFitter = new BCHistogramFitter(histCopy, funcCopy, "histo_fitter_model");

  for(unsigned int i=0;i<lines.size();i++) {
    fHistFitter->SetPriorGauss(nBkgPars+3*i+0,lines.at(i),linePosPrior);
    fHistFitter->SetPriorGauss(nBkgPars+3*i+1,fRes->Eval(lines.at(i)),fwhmPrior);
  };

  //! set mc options
  fHistFitter->SetOptimizationMethod(BCModel::kOptMinuit);
  fHistFitter->SetPrecision(precision);
  for(unsigned int i=0;i<lines.size();i++) {  
    fHistFitter->GetParameter(nBkgPars+3*i+0).SetNbins(100);   
    fHistFitter->GetParameter(nBkgPars+3*i+1).SetNbins(100);   
    fHistFitter->GetParameter(nBkgPars+3*i+2).SetNbins(100);
  }
  for(int i=0;i<nBkgPars;i++) fHistFitter->GetParameter(i).SetNbins(100);

  //! perform fit
  fHistFitter->Fit();

  fFitPerformed = true;
  return 0;
}

int GammaLineFit::Fit( TString name, vector<double> lines, pair<double,double> range,
  GammaBackground backgroundType, double linePosPrior, double fwhmPrior,
  BCEngineMCMC::Precision precision ) 
{
  Fit(lines,range,backgroundType,linePosPrior,fwhmPrior,precision);

  //! provide log output
  vector<double> bestFitPars      = fHistFitter->GetBestFitParameters();
  vector<double> bestFitParErrors = fHistFitter->GetBestFitParameterErrors();
  int nBkgPars = fBkgFunction->GetNpar();


  ordered_json foutput;

  string name_fit = name.Data();
  foutput[name_fit]["range_in_keV"] = {range.first,range.second};
  foutput[name_fit]["bin_width_keV"] = fHist->GetBinWidth(1);

  string units = "";
  for(int i=0;i<fBkgFunction->GetNpar();i++) { 
    if (i==0) units = "_in_cts";
    else if (i>0&&backgroundType!=kStep) units = "_in_cts/keV";
    else if (i>1) units = Form("_in_cts/keV^%d",i);
    foutput[name_fit]["fit_parameters"]["background"][fBkgFunction->GetParName(i) + units]["value"] = bestFitPars[i];
    foutput[name_fit]["fit_parameters"]["background"][fBkgFunction->GetParName(i)+ units]["err"] = bestFitParErrors[i];
  }

  for(unsigned int i=0;i<lines.size();i++) {    
    foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+0) + string("_in_keV")]["value"] = bestFitPars[nBkgPars+3*i+0];
    foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+0) + string("_in_keV")]["err"] = bestFitParErrors[nBkgPars+3*i+0];
    foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+0) + string("_in_keV")]["line"] = lines.at(i);

    foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+1) + string("_in_keV")]["value"] = bestFitPars[nBkgPars+3*i+1];
    foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+1) + string("_in_keV")]["err"] = bestFitParErrors[nBkgPars+3*i+1];
    foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+1) + string("_in_keV")]["resolution"] = fRes->Eval(lines.at(i));
    
    // Perform the fit and get the marginalized posterior
    fHistFitter->MarginalizeAll();

    // Get the marginalized posterior for the intensity parameter
    BCH1D intensity = fHistFitter->GetMarginalized(nBkgPars + 3 * i + 2);

    // fix the bands for each posterior at 68.3%, 95.4%, 99.7% [not there in the original src code]
    // fHistFitter->GetBCH1DdrawingOptions().SetBandType(BCH1D::kCentralInterval);
    
    // Get the mode for the intensity of the peak
    double mode = intensity.GetBestFitParameters();

    /*
    // Compute quantiles through BAT (sometimes it fails, ie high>mode when it's not the case)...
    double low=0, high=0;
    BCH1D::BCH1DSmallestInterval SI = intensity.GetSmallestIntervals(0.682689);
    if (SI.intervals.empty()) {
      log << "SI intervals are empty!" << std::endl;
    }
    low = SI.intervals.front().xmin;
    high = SI.intervals.front().xmax;
    */

    // ...another way to compute quantiles (not through BAT) 
    double low=0, high=0;
    low = intensity.GetQuantile(0.16);
    high = intensity.GetQuantile(0.84);

  //************************************************************************
    if(low>0 and mode>0) {  
      foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+2) + string("_in_cts")]["mode"] = mode;
      foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+2) + string("_in_cts")]["range_min"] = low;
      foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+2) + string("_in_cts")]["range_max"] = high;
      foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+2) + string("_in_cts")]["comments"] = "at 68%";
    } 
    else {
      foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+2) + string("_in_cts")]["upper_limit"] = intensity.GetQuantile(0.90);
      foutput[name_fit]["fit_parameters"]["line"][fFitFunction->GetParName(nBkgPars+3*i+2) + string("_in_cts")]["comments"] = "at 90%";
    }
  }
  foutput[name_fit]["fit_parameters"]["p-value"] = fHistFitter->GetPValue();
  
   
  fstream json_file; 
  json_file.open(Form("%s/%s.gamma.json",fOutputDir.Data(),GetName()), ios::in | ios::out);
  
  if (json_file){
    ordered_json foutput_tot = json::parse(json_file);
    foutput_tot.update(foutput); //, true);
    json_file.clear();
    json_file.seekp(0);
    json_file << foutput_tot.dump(2);
  }
  else{
    json_file.close();
    json_file.open(Form("%s/%s.gamma.json",fOutputDir.Data(),GetName()), ios::out);
    json_file << foutput.dump(2);
  }
  
  //! save fit to files
  TCanvas canvas;
   fHistFitter->DrawFit("HIST",false);
   canvas.Print(Form("%s/%s.%s.png",fOutputDir.Data(),GetName(),name.Data()));
  TFile* file = new TFile(Form("%s/%s.%s.root",fOutputDir.Data(),GetName(), name.Data()),"RECREATE");
   fHistFitter->GetHistogram().Write(Form("%s_hist",name.Data()));
   //fHistFitter->GetFitFunction()->Write(Form("%s_bestfit",name.Data()));
   fHistFitter->GetFitFunctionGraph(bestFitPars)->Write(Form("%s_bestfit",name.Data()));
   //fHistFitter->GetErrorBand(0.682689)/*Graph(0.16,0.84)*/.Write(Form("%s_band",name.Data()));  // not able to convert from 0.9.4 to 1.0.0
   
  file->Close();

  //! provide pdf file with posterior distributions
  fHistFitter->PrintAllMarginalized(Form("%s/%s.%s.pdf",fOutputDir.Data(),GetName(),name.Data()), 2,2);
  //! provide pdf file with prior+posterior distributions
  fHistFitter->PrintKnowledgeUpdatePlots(Form("%s/%s.priorsANDposteriors.%s.pdf",fOutputDir.Data(),GetName(),name.Data()), 2,2);
  fHistFitter->PrintSummary();
  fHistFitter->WriteMarginalizedDistributions(Form("%s/%s_marginalized.%s.root",fOutputDir.Data(),GetName(), name.Data()), "RECREATE");

  return 0;
}

} // namespace
